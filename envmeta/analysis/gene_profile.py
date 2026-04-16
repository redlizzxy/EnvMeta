"""MAG 元素循环基因谱热图。

对标 scripts/python/06_MAG_gene_profile.py。核心：MAG × KO 拷贝数热图，按
元素分块（As / N / S / Fe），色深 = log1p(copies)。

S6-fix2 统一化：与 mag_heatmap / pathway / mag_quality 共享 4 层参数；行标签
改用 Genus species / Genus sp. Mx_XX / MAG_id；左侧门彩条 + 右侧门图例；
配色默认 **viridis**（感知均匀、色域宽）替换原单色 YlOrBr。

输入：
    ko_annotation_df: MAG × KO 长表（MAG + KEGG_ko）
    taxonomy_df:      可选 MAG + classification
    keystone_df:      可选 MAG
    abundance_df:     可选 MAG × sample（排序 / Top-N 筛选）

输出：
    figure — 主热图（MAG 行 × KO 列，按元素排序，log1p 拷贝数）
    stats  — MAG × KO 拷贝数长表 + 元数据 + gene_count
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from envmeta.analysis import _mag_common as _mc
from envmeta.analysis._mag_common import PHYLUM_COLORS  # re-export
from envmeta.analysis.base import AnalysisResult
from envmeta.geocycle.knowledge_base import (
    element_colors, element_pathway_ko_order, flat_ko_map,
)

# BC helpers
_parse_ko_annotation = None  # 不再导出；pathway 里仍有同名函数可用


DEFAULTS = {
    # Layer 1 — MAG 子集过滤
    "filter_mode": "top_plus_keystone",
    "top_n_by": "mean",
    "top_n_count": 30,
    "max_mags": 60,
    # Layer 2 — 视觉
    "highlight_keystones": True,
    "show_phylum_bar": True,
    "show_phylum_legend": True,
    "width_mm": 340,
    "height_mm": 220,
    # Layer 3 — 行排序
    "row_order": "phylum_cluster",
    "linkage_method": "average",
    # Layer 4 — 基因谱特有
    "element_filter": None,
    "cmap_name": "viridis",               # ← 改默认（原 YlOrBr）
    "show_gene_names": True,
    "drop_zero_kos": True,
    "show_element_bar": True,
    # 0 值视觉处理（解决"深色占大部分视觉"问题）
    "blank_zeros": True,                  # True → 0 值显白色（缺失 vs 低表达一眼分）
    "zero_threshold": 0.0,                # 原始拷贝数 ≤ 此值视为 0；设 >0 可把"低"也当缺失
    "sort_ko_by_coverage": False,         # True → KO 列按覆盖率降序（左密右稀）
    # 向后兼容旧 key
    "sort_by": None,
    "annotate_keystone": None,
    "top_abundance_n": None,              # 旧 S6-fix name
}


def _normalize_deprecated_params(p: dict) -> dict:
    if p.get("sort_by") is not None:
        mapping = {
            "phylum_then_count": "phylum_cluster",
            "phylum_cluster": "phylum_cluster",
            "count": "metric_desc",
            "metric_desc": "metric_desc",
            "abundance": "abundance",
        }
        p["row_order"] = mapping.get(p["sort_by"],
                                     p.get("row_order", "phylum_cluster"))
    if p.get("annotate_keystone") is not None:
        p["highlight_keystones"] = bool(p["annotate_keystone"])
    if p.get("top_abundance_n") is not None:
        p["top_n_count"] = int(p["top_abundance_n"])
    return p


def _parse_ko_copies(df: pd.DataFrame) -> dict[str, dict[str, int]]:
    """把 MAG+KEGG_ko 长表解析为 {MAG: {KO: copies}}。"""
    mag_c = _mc.mag_col(df)
    ko_c = next((c for c in df.columns
                 if c.lower() in ("kegg_ko", "kegg ko", "ko", "ko_id")),
                df.columns[1] if len(df.columns) > 1 else None)
    if ko_c is None:
        raise ValueError("无法定位 KEGG_ko / KO 列")
    out: dict[str, dict[str, int]] = {}
    for _, r in df.iterrows():
        mag = str(r[mag_c])
        raw = str(r[ko_c])
        if raw in ("", "nan", "None"):
            continue
        for ko in raw.split(","):
            k = ko.strip().replace("ko:", "")
            if k and k.startswith("K"):
                out.setdefault(mag, {})[k] = out.get(mag, {}).get(k, 0) + 1
    return out


def analyze(
    ko_annotation_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    abundance_df: pd.DataFrame | None = None,
    params: dict | None = None,
) -> AnalysisResult:
    p = _normalize_deprecated_params({**DEFAULTS, **(params or {})})

    # KO × 元素 × 通路 顺序
    mag_copies = _parse_ko_copies(ko_annotation_df)
    ko_order_triples = element_pathway_ko_order()
    ko_map = flat_ko_map()
    elem_color = element_colors()
    if p["element_filter"]:
        keep = set(p["element_filter"])
        ko_order_triples = [(e, pw, ko) for (e, pw, ko) in ko_order_triples
                            if e in keep]
    kb_ko_list = [ko for (_, _, ko) in ko_order_triples]
    kb_ko_element = {ko: e for (e, _, ko) in ko_order_triples}

    # MAG 全集
    if taxonomy_df is not None and not taxonomy_df.empty:
        tax = taxonomy_df.copy()
        tax = tax.rename(columns={_mc.mag_col(tax): "MAG"})
        all_mags = sorted(tax["MAG"].astype(str).unique().tolist())
    else:
        all_mags = sorted(mag_copies.keys())
    if not all_mags:
        raise ValueError("无 MAG 可分析")

    mat = np.zeros((len(all_mags), len(kb_ko_list)), dtype=float)
    for i, mag in enumerate(all_mags):
        owned = mag_copies.get(mag, {})
        for j, ko in enumerate(kb_ko_list):
            mat[i, j] = owned.get(ko, 0)

    # 过滤全零 KO 列
    ko_indices_keep = list(range(len(kb_ko_list)))
    if p["drop_zero_kos"]:
        col_sum = mat.sum(axis=0)
        ko_indices_keep = [j for j, s in enumerate(col_sum) if s > 0]
        mat = mat[:, ko_indices_keep]
    active_kos = [kb_ko_list[j] for j in ko_indices_keep]
    active_elements = [kb_ko_element[k] for k in active_kos]

    df = pd.DataFrame(mat, index=all_mags, columns=active_kos).reset_index()
    df = df.rename(columns={"index": "MAG"})

    # Phylum / Genus / Species / label
    df = _mc.annotate_taxonomy(df, taxonomy_df)
    df = _mc.annotate_keystone(df, keystone_df)

    # 丰度均值
    if abundance_df is not None and not abundance_df.empty:
        ab = abundance_df.copy()
        ab = ab.rename(columns={_mc.mag_col(ab): "MAG"})
        scols = [c for c in ab.columns if c != "MAG"]
        ab["abundance_mean"] = ab[scols].apply(
            pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)
        df = df.merge(ab[["MAG", "abundance_mean"]], on="MAG", how="left")
        df["abundance_mean"] = df["abundance_mean"].fillna(0.0)
    else:
        df["abundance_mean"] = 0.0

    df["gene_count"] = df[active_kos].sum(axis=1)

    # Layer 1 — filter_mode
    df = _mc.apply_filter_mode(
        df,
        mode=p["filter_mode"],
        top_n_count=int(p["top_n_count"]),
        top_n_by=p["top_n_by"],
        score_col="abundance_mean" if p["top_n_by"] != "variance" else "gene_count",
    )

    # Layer 3 — 行排序
    df = _mc.order_rows(
        df,
        mode=p["row_order"],
        metric_col="gene_count",
        abundance_col="abundance_mean",
        cluster_matrix=df[active_kos].to_numpy(dtype=float),
        linkage_method=p["linkage_method"],
        log_transform_cluster=True,
    )

    # max_mags 硬截断
    if p.get("max_mags"):
        df = df.head(int(p["max_mags"])).reset_index(drop=True)

    # 可选：KO 列按当前 MAG 子集里的覆盖率降序（左密右稀）
    if p.get("sort_ko_by_coverage"):
        sub = df[active_kos].to_numpy(dtype=float)
        coverage = (sub > 0).sum(axis=0)
        order = np.argsort(-coverage, kind="stable")
        active_kos = [active_kos[i] for i in order]
        active_elements = [active_elements[i] for i in order]

    mat_sorted = df[active_kos].to_numpy()
    mat_log = np.log1p(mat_sorted)

    fig = _draw(df, active_kos, active_elements, mat_sorted, mat_log,
                elem_color, ko_map, p)

    stats_df = df[["MAG", "label", "Phylum", "Genus", "Species",
                   "is_keystone", "abundance_mean",
                   "gene_count"] + active_kos].copy()
    return AnalysisResult(figure=fig, stats=stats_df, params=p)


def _draw(df, active_kos, active_elements, mat_raw, mat_log,
          elem_color, ko_map, p) -> plt.Figure:
    n_mag, n_ko = mat_log.shape
    # gridspec 三列布局（同 pathway）
    fig = plt.figure(
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4),
        constrained_layout=True,
    )
    show_phy = p["show_phylum_bar"]
    show_leg = p["show_phylum_legend"]
    wr_phy = 0.04 if show_phy else 0.0001
    wr_leg = 0.25 if show_leg else 0.0001
    gs = fig.add_gridspec(1, 3, width_ratios=[wr_phy, 1.0, wr_leg], wspace=0.02)
    ax_phy = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    ax_leg = fig.add_subplot(gs[0, 2])
    ax_leg.axis("off")

    cmap = plt.get_cmap(p["cmap_name"]).copy()
    vmax = mat_log.max() if mat_log.max() > 0 else 1.0
    # 0 值（或 ≤ zero_threshold）留白 —— 解决"深色占大部分视觉"问题
    if p.get("blank_zeros", True):
        cmap.set_bad("white")
        thresh = float(p.get("zero_threshold", 0.0))
        display_mat = np.ma.masked_where(mat_raw <= thresh, mat_log)
        # vmin 设为下一档色度避免最低彩色块与白色混淆
        vmin = 1e-3
    else:
        display_mat = mat_log
        vmin = 0
    im = ax.imshow(display_mat, aspect="auto", cmap=cmap,
                   vmin=vmin, vmax=vmax, interpolation="nearest")

    # 元素色带（顶部）
    if p["show_element_bar"]:
        for j, el in enumerate(active_elements):
            ax.add_patch(Rectangle(
                (j - 0.5, -0.9), 1.0, 0.35,
                color=elem_color.get(el, "#888"), clip_on=False, zorder=3))

    # 左侧门彩条（独立 Axes）
    if show_phy:
        _mc.draw_phylum_bar(ax_phy, df["Phylum"].tolist())
    else:
        ax_phy.axis("off")

    # 列标签
    if p["show_gene_names"]:
        col_labels = [ko_map.get(ko, (ko,))[0] for ko in active_kos]
    else:
        col_labels = active_kos
    ax.set_xticks(range(n_ko))
    ax.set_xticklabels(col_labels, rotation=75, ha="right",
                       fontsize=max(4, min(8, 400 // max(n_ko, 1))))

    # 行标签（Genus / Genus sp. / MAG_id）
    y_labels = [
        ("★ " if (p["highlight_keystones"] and r["is_keystone"]) else "")
        + str(r.get("label", r["MAG"]))
        for _, r in df.iterrows()
    ]
    ax.set_yticks(range(n_mag))
    ax.set_yticklabels(y_labels,
                       fontsize=max(3, min(8, 380 // max(n_mag, 1))),
                       fontstyle="italic")

    # 元素分组竖线
    if p["show_element_bar"] and n_ko > 1:
        prev = active_elements[0]
        for j, el in enumerate(active_elements):
            if el != prev:
                ax.axvline(j - 0.5, color="white", lw=1.0, zorder=4)
                prev = el

    ax.set_title(
        f"MAG Element-Cycle Gene Profile  ({n_mag} MAGs × {n_ko} KOs,"
        f" color = log1p(copies))",
        fontsize=11, fontweight="bold", pad=12,
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # 右侧图例区：Phylum + keystone 说明 + colorbar
    if show_leg:
        _mc.draw_phylum_legend(ax_leg, df["Phylum"].tolist())
        if p["highlight_keystones"] and df["is_keystone"].any():
            _mc.draw_keystone_note(ax_leg)
    cax = ax_leg.inset_axes([0.12, 0.0, 0.25, 0.18])
    cbar = fig.colorbar(im, cax=cax, orientation="vertical")
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label("log1p(copies)", fontsize=6.5)
    return fig
