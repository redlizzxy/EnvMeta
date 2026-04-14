"""MAG 元素循环基因谱热图。

对标 `scripts/python/06_MAG_gene_profile.py`。核心：MAG × KO 拷贝数热图，按
元素分块（As / N / S / Fe），色深 = log1p(copies)。

输入：
    ko_annotation_df: MAG × KO 长表（MAG + KEGG_ko），可含重复 → 拷贝数
    taxonomy_df:      可选 MAG + classification
    keystone_df:      可选 MAG + Genus
    abundance_df:     可选 MAG × sample（排序 / keystone 丰度筛选）

输出：
    figure — 主热图（MAG 行 × KO 列，按元素排序，log1p 拷贝数）
    stats  — MAG × KO 拷贝数长表 + 元数据 + 基因总数
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from envmeta.analysis.base import AnalysisResult
from envmeta.geocycle.knowledge_base import (
    element_colors, element_pathway_ko_order, flat_ko_map,
)
from envmeta.analysis.pathway import (
    PHYLUM_COLORS, _extract_phylum, _mag_col, _parse_ko_annotation,
)

DEFAULTS = {
    "width_mm": 340,
    "height_mm": 220,
    "max_mags": None,            # None → 全部
    "sort_by": "phylum_then_count",  # | "count" | "abundance"
    "element_filter": None,
    "annotate_keystone": True,
    "cmap_name": "YlOrBr",
    "drop_zero_kos": True,       # True → 过滤全 0 KO 列
    "show_gene_names": True,     # KO 列标签是否显示基因名
    "show_element_bar": True,
    "show_phylum_bar": True,
}


def _parse_ko_copies(df: pd.DataFrame) -> dict[str, dict[str, int]]:
    """把 MAG+KEGG_ko 长表解析为 {MAG: {KO: copies}}（同一 MAG 多次出现 = 拷贝 > 1）。"""
    mag_c = _mag_col(df)
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
    p = {**DEFAULTS, **(params or {})}

    # 拷贝数
    mag_copies = _parse_ko_copies(ko_annotation_df)
    # 知识库：按 element → pathway → KO 的自然顺序
    ko_order_triples = element_pathway_ko_order()
    ko_map = flat_ko_map()  # {ko: (gene_name, pathway, element)}
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
        tax = tax.rename(columns={_mag_col(tax): "MAG"})
        cls_col = next((c for c in tax.columns
                        if "classif" in c.lower() or "taxonomy" in c.lower()),
                       tax.columns[1] if len(tax.columns) > 1 else None)
        tax["Phylum"] = tax[cls_col].apply(_extract_phylum) if cls_col else "Unknown"
        all_mags = sorted(tax["MAG"].astype(str).unique().tolist())
    else:
        tax = None
        all_mags = sorted(mag_copies.keys())
    if not all_mags:
        raise ValueError("无 MAG 可分析")

    # 拷贝数矩阵
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

    # DataFrame
    df = pd.DataFrame(mat, index=all_mags, columns=active_kos).reset_index()
    df = df.rename(columns={"index": "MAG"})
    if tax is not None:
        df = df.merge(tax[["MAG", "Phylum"]], on="MAG", how="left")
    df["Phylum"] = df.get("Phylum", pd.Series(["Unknown"] * len(df))).fillna("Unknown")

    if keystone_df is not None and not keystone_df.empty:
        ks = keystone_df.copy()
        ks = ks.rename(columns={_mag_col(ks): "MAG"})
        df["is_keystone"] = df["MAG"].isin(set(ks["MAG"].astype(str)))
    else:
        df["is_keystone"] = False

    if abundance_df is not None and not abundance_df.empty:
        ab = abundance_df.copy()
        ab = ab.rename(columns={_mag_col(ab): "MAG"})
        scols = [c for c in ab.columns if c != "MAG"]
        ab["abundance_mean"] = ab[scols].apply(
            pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)
        df = df.merge(ab[["MAG", "abundance_mean"]], on="MAG", how="left")
        df["abundance_mean"] = df["abundance_mean"].fillna(0.0)
    else:
        df["abundance_mean"] = 0.0

    df["gene_count"] = df[active_kos].sum(axis=1)

    # 排序
    if p["sort_by"] == "count":
        df = df.sort_values("gene_count", ascending=False)
    elif p["sort_by"] == "abundance":
        df = df.sort_values("abundance_mean", ascending=False)
    else:  # phylum_then_count
        phy_rank = {ph: i for i, ph in
                    enumerate(df["Phylum"].value_counts().index.tolist())}
        df["_phy_rank"] = df["Phylum"].map(phy_rank)
        df = df.sort_values(["_phy_rank", "gene_count"],
                            ascending=[True, False]).drop(columns=["_phy_rank"])
    df = df.reset_index(drop=True)

    if p["max_mags"] and len(df) > p["max_mags"]:
        df = df.head(p["max_mags"]).reset_index(drop=True)

    # 绘图用矩阵（按排序后的 MAG 顺序）
    mat_sorted = df[active_kos].to_numpy()
    mat_log = np.log1p(mat_sorted)

    fig = _draw(df, active_kos, active_elements, mat_log, elem_color,
                ko_map, p)

    stats_df = df[["MAG", "Phylum", "is_keystone", "abundance_mean",
                   "gene_count"] + active_kos].copy()
    return AnalysisResult(figure=fig, stats=stats_df, params=p)


def _draw(df, active_kos, active_elements, mat_log, elem_color, ko_map, p) -> plt.Figure:
    n_mag, n_ko = mat_log.shape
    fig, ax = plt.subplots(
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4),
        constrained_layout=True,
    )
    cmap = plt.get_cmap(p["cmap_name"])
    vmax = mat_log.max() if mat_log.max() > 0 else 1.0
    im = ax.imshow(mat_log, aspect="auto", cmap=cmap, vmin=0, vmax=vmax,
                   interpolation="nearest")

    # 元素色带（顶部）
    if p["show_element_bar"]:
        for j, el in enumerate(active_elements):
            ax.add_patch(Rectangle(
                (j - 0.5, -0.9), 1.0, 0.35,
                color=elem_color.get(el, "#888"), clip_on=False, zorder=3))

    # 门色带（左侧）
    if p["show_phylum_bar"]:
        for i, (_, r) in enumerate(df.iterrows()):
            ax.add_patch(Rectangle(
                (-0.6, i - 0.5), 0.4, 1.0,
                color=PHYLUM_COLORS.get(r["Phylum"], "#888"),
                clip_on=False, zorder=3))

    # 列标签：gene_name (ko) 或 ko
    if p["show_gene_names"]:
        col_labels = []
        for ko in active_kos:
            gene = ko_map.get(ko, (ko,))[0]
            col_labels.append(f"{gene}")
    else:
        col_labels = active_kos
    ax.set_xticks(range(n_ko))
    ax.set_xticklabels(col_labels, rotation=75, ha="right",
                       fontsize=max(4, min(8, 400 // max(n_ko, 1))))

    # 行标签（MAG + ★）
    y_labels = [
        ("★ " if (p["annotate_keystone"] and r["is_keystone"]) else "") + r["MAG"]
        for _, r in df.iterrows()
    ]
    ax.set_yticks(range(n_mag))
    ax.set_yticklabels(y_labels, fontsize=max(3, min(8, 380 // max(n_mag, 1))))

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
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, label="log1p(copies)")
    cbar.ax.tick_params(labelsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    return fig
