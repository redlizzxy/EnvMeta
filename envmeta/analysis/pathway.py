"""MAG 通路完整度（heatmap / bubble）。

对标 scripts/python/08_pathway_completeness.py。核心：以知识库（elements.json）
定义的 18 个元素循环通路为列，每个 MAG 按"拥有 KO / 通路 KO 总数"计算完整度。

S6-fix2 统一化：与 mag_heatmap / gene_profile / mag_quality 共享 4 层参数
（filter_mode / top_n_by / top_n_count / max_mags + Layer 2 视觉 + Layer 3
行排序）+ 共享 Phylum 彩条 + 门图例 + Genus 标签规则。

输入：
    ko_annotation_df: MAG × KO 注释（长表：MAG + KEGG_ko）
    taxonomy_df:      可选，MAG + GTDB classification
    keystone_df:      可选，MAG + Genus
    abundance_df:     可选，MAG × sample（用于 Top-N 过滤与气泡图大小）

输出：
    figure — heatmap 或 bubble（按 style 参数）
    stats  — MAG × pathway 完整度长表 + Phylum / Genus / Species / label /
             is_keystone + abundance_mean + total_completeness
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

from envmeta.analysis import _mag_common as _mc
from envmeta.analysis._mag_common import PHYLUM_COLORS  # re-export for BC
from envmeta.analysis.base import AnalysisResult
from envmeta.geocycle.knowledge_base import (
    element_colors, pathway_display, pathway_element_map, pathway_ko_sets,
)

# 向后兼容的辅助：旧代码可能直接导入这两个
_mag_col = _mc.mag_col
_extract_phylum = _mc.extract_phylum


DEFAULTS = {
    # Layer 1 — MAG 子集过滤（共享）
    "filter_mode": "top_plus_keystone",   # all|top_n|keystone_only|top_plus_keystone
    "top_n_by": "mean",                   # mean|sum|variance
    "top_n_count": 30,
    "max_mags": 50,                       # 硬截断；0/None → 不截
    # Layer 2 — 视觉（共享）
    "highlight_keystones": True,
    "show_phylum_bar": True,
    "show_phylum_legend": True,
    "width_mm": 320,
    "height_mm": 220,
    # Layer 3 — 行排序（共享）
    "row_order": "phylum_cluster",        # phylum_cluster|metric_desc|abundance
    "linkage_method": "average",
    # Layer 4 — 通路特有
    "style": "heatmap",                   # heatmap|bubble
    "element_filter": None,
    "min_completeness": 0.0,
    "bubble_scale": 4.0,
    "cmap_name": "YlGnBu",
    # 向后兼容旧 key（会被规范化）
    "sort_by": None,                      # deprecated
    "show_phylum_stripe": None,           # deprecated
    "annotate_keystone": None,            # deprecated
    "annotate_top_n": None,               # deprecated
}


def _normalize_deprecated_params(p: dict) -> dict:
    """把旧 param key 迁移到新名字（保证老 bundle / 测试不破）。"""
    if p.get("sort_by") is not None:
        mapping = {
            "phylum_then_total": "phylum_cluster",
            "phylum_cluster": "phylum_cluster",
            "total": "metric_desc",
            "metric_desc": "metric_desc",
            "abundance": "abundance",
        }
        p["row_order"] = mapping.get(p["sort_by"], p.get("row_order", "phylum_cluster"))
    if p.get("show_phylum_stripe") is not None:
        p["show_phylum_bar"] = bool(p["show_phylum_stripe"])
    if p.get("annotate_keystone") is not None:
        p["highlight_keystones"] = bool(p["annotate_keystone"])
    return p


def _parse_ko_annotation(df: pd.DataFrame) -> dict[str, set[str]]:
    """MAG+KEGG_ko 长表 → {MAG: {KO}}；KEGG_ko 可含 'ko:' 前缀与 ',' 分隔。"""
    mag_c = _mc.mag_col(df)
    ko_c = next((c for c in df.columns
                 if c.lower() in ("kegg_ko", "kegg ko", "ko", "ko_id")),
                df.columns[1] if len(df.columns) > 1 else None)
    if ko_c is None:
        raise ValueError("无法定位 KEGG_ko / KO 列")
    mag_kos: dict[str, set[str]] = {}
    for _, r in df.iterrows():
        mag = str(r[mag_c])
        raw = str(r[ko_c])
        if raw in ("", "nan", "None"):
            continue
        for ko in raw.split(","):
            k = ko.strip().replace("ko:", "")
            if k and k.startswith("K"):
                mag_kos.setdefault(mag, set()).add(k)
    return mag_kos


def _compute_completeness(
    mag_kos: dict[str, set[str]],
    pw_kos: dict[str, list[str]],
    all_mags: list[str],
) -> pd.DataFrame:
    rows = []
    for m in all_mags:
        owned = mag_kos.get(m, set())
        row = {"MAG": m}
        for pw_id, kos in pw_kos.items():
            detected = owned & set(kos)
            row[pw_id] = len(detected) / len(kos) * 100 if kos else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def analyze(
    ko_annotation_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    abundance_df: pd.DataFrame | None = None,
    params: dict | None = None,
) -> AnalysisResult:
    p = _normalize_deprecated_params({**DEFAULTS, **(params or {})})

    # === KO → 通路 ===
    mag_kos = _parse_ko_annotation(ko_annotation_df)
    pw_kos = pathway_ko_sets()
    pw_elem = pathway_element_map()
    pw_display = pathway_display(lang="en")
    elem_color = element_colors()
    if p["element_filter"]:
        keep = set(p["element_filter"])
        pw_kos = {k: v for k, v in pw_kos.items() if pw_elem[k] in keep}
    pw_order = list(pw_kos.keys())

    # === MAG 全集 ===
    if taxonomy_df is not None and not taxonomy_df.empty:
        tax = taxonomy_df.copy()
        tax = tax.rename(columns={_mc.mag_col(tax): "MAG"})
        all_mags = sorted(tax["MAG"].astype(str).unique().tolist())
    else:
        all_mags = sorted(mag_kos.keys())
    if not all_mags:
        raise ValueError("无 MAG 可分析（检查 KO 注释表非空）")

    df = _compute_completeness(mag_kos, pw_kos, all_mags)

    # 注释 Phylum / Genus / Species / label（统一规则）
    df = _mc.annotate_taxonomy(df, taxonomy_df)
    df = _mc.annotate_keystone(df, keystone_df)

    # 丰度均值
    if abundance_df is not None and not abundance_df.empty:
        ab = abundance_df.copy()
        ab = ab.rename(columns={_mc.mag_col(ab): "MAG"})
        sample_cols = [c for c in ab.columns if c != "MAG"]
        ab["abundance_mean"] = ab[sample_cols].apply(
            pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)
        df = df.merge(ab[["MAG", "abundance_mean"]], on="MAG", how="left")
        df["abundance_mean"] = df["abundance_mean"].fillna(0.0)
    else:
        df["abundance_mean"] = 0.0

    # 总完整度（行主 metric）
    df["total_completeness"] = df[pw_order].sum(axis=1)

    # 最低总完整度过滤
    df = df[df["total_completeness"] >= p["min_completeness"]].copy()
    if df.empty:
        raise ValueError(f"无 MAG 总完整度 ≥ {p['min_completeness']}")
    df = df.reset_index(drop=True)

    # Layer 1 — filter_mode 子集过滤
    df = _mc.apply_filter_mode(
        df,
        mode=p["filter_mode"],
        top_n_count=int(p["top_n_count"]),
        top_n_by=p["top_n_by"],
        score_col="abundance_mean" if p["top_n_by"] != "variance" else "total_completeness",
    )

    # Layer 3 — 行排序（phylum_cluster 用完整度向量做聚类）
    df = _mc.order_rows(
        df,
        mode=p["row_order"],
        metric_col="total_completeness",
        abundance_col="abundance_mean",
        cluster_matrix=df[pw_order].to_numpy(dtype=float),
        linkage_method=p["linkage_method"],
        log_transform_cluster=False,  # 完整度本身 0-100 无需 log1p
    )

    # max_mags 硬截断
    if p.get("max_mags"):
        df = df.head(int(p["max_mags"])).reset_index(drop=True)

    # === 绘图 ===
    if p["style"] == "bubble":
        fig = _draw_bubble(df, pw_order, pw_display, pw_elem, elem_color, p)
    else:
        fig = _draw_heatmap(df, pw_order, pw_display, pw_elem, elem_color, p)

    stats_df = df[["MAG", "label", "Phylum", "Genus", "Species",
                   "is_keystone", "abundance_mean",
                   "total_completeness"] + pw_order].copy()
    return AnalysisResult(figure=fig, stats=stats_df, params=p)


# ── 共享布局辅助 ────────────────────────────────────────────

def _make_layout(p: dict):
    """3 列 gridspec：[门彩条 | 主图 | 右侧图例]。"""
    fig = plt.figure(
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4),
        constrained_layout=True,
    )
    show_phy = p["show_phylum_bar"]
    show_leg = p["show_phylum_legend"]
    wr_phy = 0.05 if show_phy else 0.0001
    wr_leg = 0.30 if show_leg else 0.0001
    gs = fig.add_gridspec(1, 3, width_ratios=[wr_phy, 1.0, wr_leg], wspace=0.02)
    ax_phy = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    ax_leg = fig.add_subplot(gs[0, 2])
    ax_leg.axis("off")
    return fig, ax_phy, ax, ax_leg


def _row_labels(df, highlight_keystones: bool) -> list[str]:
    labels = []
    for _, r in df.iterrows():
        prefix = "★ " if (highlight_keystones and r.get("is_keystone", False)) else ""
        labels.append(prefix + str(r.get("label", r["MAG"])))
    return labels


# ── heatmap ─────────────────────────────────────────────────

def _draw_heatmap(df, pw_order, pw_display, pw_elem, elem_color, p) -> plt.Figure:
    n_mag = len(df)
    n_pw = len(pw_order)
    fig, ax_phy, ax, ax_leg = _make_layout(p)

    mat = df[pw_order].to_numpy()
    cmap = plt.get_cmap(p["cmap_name"])
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=100,
                   interpolation="nearest")

    # 左侧 Phylum 彩条（独立 Axes）
    if p["show_phylum_bar"]:
        _mc.draw_phylum_bar(ax_phy, df["Phylum"].tolist())
    else:
        ax_phy.axis("off")

    # 顶部元素彩条（画在主图 Axes 内）
    elem_of_col = [pw_elem[pw] for pw in pw_order]
    for j, el in enumerate(elem_of_col):
        ax.add_patch(Rectangle((j - 0.5, -0.9), 1.0, 0.4,
                               color=elem_color.get(el, "#888"),
                               clip_on=False, zorder=3))

    # 行标签（Genus species / Genus sp. Mx_XX / MAG_id）
    ax.set_yticks(range(n_mag))
    ax.set_yticklabels(_row_labels(df, p["highlight_keystones"]),
                       fontsize=max(3, min(8, 360 // max(n_mag, 1))),
                       fontstyle="italic")

    ax.set_xticks(range(n_pw))
    ax.set_xticklabels([pw_display.get(pw, pw) for pw in pw_order],
                       rotation=45, ha="right", fontsize=8)
    ax.set_title(
        f"Pathway Completeness (MAG × pathway)  n={n_mag} MAGs × {n_pw} pathways",
        fontsize=11, fontweight="bold", pad=12,
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # 右侧图例 + colorbar
    if p["show_phylum_legend"]:
        _mc.draw_phylum_legend(ax_leg, df["Phylum"].tolist())
        if p["highlight_keystones"] and df["is_keystone"].any():
            _mc.draw_keystone_note(ax_leg)
    cax = ax_leg.inset_axes([0.1, 0.0, 0.22, 0.18])
    cbar = fig.colorbar(im, cax=cax, orientation="vertical")
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label("Completeness (%)", fontsize=6.5)
    return fig


# ── bubble ──────────────────────────────────────────────────

def _draw_bubble(df, pw_order, pw_display, pw_elem, elem_color, p) -> plt.Figure:
    n_mag = len(df)
    n_pw = len(pw_order)
    fig, ax_phy, ax, ax_leg = _make_layout(p)

    mat = df[pw_order].to_numpy()
    abund = df["abundance_mean"].to_numpy()
    ab_scale = abund / (abund.max() + 1e-9) * p["bubble_scale"] * 80 + 8

    cmap = plt.get_cmap(p["cmap_name"])
    norm = Normalize(vmin=0, vmax=100)
    for i in range(n_mag):
        for j in range(n_pw):
            val = mat[i, j]
            if val <= 0:
                continue
            ax.scatter(j, i, s=ab_scale[i] * (val / 100) + 10,
                       c=[cmap(val / 100)], edgecolors="white", linewidths=0.4,
                       alpha=0.85, zorder=3)

    # 左侧 Phylum 彩条
    if p["show_phylum_bar"]:
        _mc.draw_phylum_bar(ax_phy, df["Phylum"].tolist())
    else:
        ax_phy.axis("off")

    # 顶部元素彩条
    elem_of_col = [pw_elem[pw] for pw in pw_order]
    for j, el in enumerate(elem_of_col):
        ax.add_patch(Rectangle((j - 0.5, -0.9), 1.0, 0.4,
                               color=elem_color.get(el, "#888"),
                               clip_on=False, zorder=2))

    ax.set_xticks(range(n_pw))
    ax.set_xticklabels([pw_display.get(pw, pw) for pw in pw_order],
                       rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n_mag))
    ax.set_yticklabels(_row_labels(df, p["highlight_keystones"]),
                       fontsize=max(4, min(9, 360 // max(n_mag, 1))),
                       fontstyle="italic")
    ax.set_xlim(-1, n_pw)
    ax.set_ylim(-1.5, n_mag)
    ax.invert_yaxis()
    ax.set_title(
        "Pathway Completeness Bubble  "
        "(color=completeness %, bubble size ≈ abundance × completeness)",
        fontsize=10, fontweight="bold", pad=12,
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # 右侧图例 + colorbar（completeness）+ size legend
    if p["show_phylum_legend"]:
        _mc.draw_phylum_legend(ax_leg, df["Phylum"].tolist())
        if p["highlight_keystones"] and df["is_keystone"].any():
            _mc.draw_keystone_note(ax_leg)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cax = ax_leg.inset_axes([0.1, 0.0, 0.22, 0.18])
    cbar = fig.colorbar(sm, cax=cax, orientation="vertical")
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label("Completeness (%)", fontsize=6.5)

    from matplotlib.lines import Line2D
    legend_vals = [20, 50, 100]
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#555",
               markersize=np.sqrt(v * 0.4 * p["bubble_scale"]),
               markeredgecolor="white", markeredgewidth=0.3,
               label=f"{v}% × top-abund.")
        for v in legend_vals
    ]
    leg2 = ax_leg.legend(
        handles=handles,
        loc="lower left", bbox_to_anchor=(0.0, 0.22),
        fontsize=6, frameon=False,
        title="bubble size", title_fontsize=6.5,
    )
    ax_leg.add_artist(leg2)

    return fig
