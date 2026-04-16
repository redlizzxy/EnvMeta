"""Top-N MAG 丰度热图（对标 scripts/python/07_MAG_abundance_heatmap.py）。

S6-fix2 统一化：4 张 MAG 图共享 Layer 1-3 参数（filter_mode / top_n_by /
top_n_count / max_mags / row_order）+ 共享 Phylum 彩条 + 门图例 + Genus 标签。

保留 mag_heatmap 特有视觉：
- 三段非线性配色（Blues→YlGn→YlOrRd，适配少数极高 + 多数低值的长尾分布）
- 顶部组彩条（SampleID × Group）+ 组间竖线
- 可选行/列聚类 + 聚类前 log1p 变换
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.patches import Rectangle

from envmeta.analysis import _mag_common as _mc
from envmeta.analysis._mag_common import (GROUP_COLORS, PHYLUM_COLORS)
from envmeta.analysis.base import AnalysisResult


DEFAULTS = {
    # Layer 1 — MAG 子集过滤（共享）
    "filter_mode": "top_plus_keystone",
    "top_n_by": "mean",                   # mean|sum|variance
    "top_n_count": 30,
    "max_mags": 0,                        # 0 → 不截（Top-N 已经是子集）
    # Layer 2 — 视觉（共享）
    "highlight_keystones": True,
    "show_phylum_bar": True,
    "show_phylum_legend": True,
    "width_mm": 240,
    "height_mm": 220,
    # Layer 3 — 行排序（共享）
    "row_order": "phylum_cluster",        # phylum_cluster|metric_desc|abundance
    "linkage_method": "average",
    # Layer 4 — mag_heatmap 特有
    "color_breakpoints": (0.2, 0.5),
    "log_transform": True,
    "cluster_cols": False,
    "show_group_bar": True,
    # 向后兼容旧 key
    "top_n": None,                        # 老 API: top_n_count 别名
    "selection_by": None,                 # 老 API: top_n_by 别名
    "cluster_rows": None,                 # 已由 row_order 覆盖
    "cluster_within_phylum": None,        # 已由 row_order=phylum_cluster 覆盖
}


def _normalize_deprecated_params(p: dict) -> dict:
    if p.get("top_n") is not None:
        p["top_n_count"] = int(p["top_n"])
    if p.get("selection_by") is not None:
        p["top_n_by"] = p["selection_by"]
    # cluster_rows=False 表示不想聚类 → row_order=metric_desc（按 selection_score）
    if p.get("cluster_rows") is False and p.get("row_order") == "phylum_cluster":
        p["row_order"] = "metric_desc"
    # cluster_within_phylum=False 且 cluster_rows=True 表示用全局聚类（不再支持；
    # 映射为 phylum_cluster 仍可，本次不单独支持）
    return p


def _build_tricolor_cmap(mat: np.ndarray, breakpoints: tuple[float, float]):
    """三段非线性配色：Blues（低频）→ YlGn（中）→ YlOrRd（高）。"""
    vmax = np.ceil(mat.max() * 10) / 10 if mat.size else 1.0
    lo_bp, hi_bp = breakpoints
    if vmax < hi_bp + 0.1:
        vmax = hi_bp + 0.1
    bounds_lo = list(np.arange(0, lo_bp + 1e-9, 0.025))
    bounds_mid = list(np.arange(lo_bp, hi_bp + 1e-9, 0.05))
    bounds_hi = list(np.arange(hi_bp, min(vmax, 1.0) + 1e-9, 0.1))
    if vmax > 1.0:
        bounds_hi += list(np.arange(1.0, vmax + 0.01, 0.5))
    bounds = sorted({round(b, 3) for b in bounds_lo + bounds_mid + bounds_hi})
    if bounds[-1] < vmax:
        bounds.append(round(float(vmax), 2))
    n_lo = sum(1 for b in bounds if b <= lo_bp) - 1
    n_mid = sum(1 for b in bounds if lo_bp < b <= hi_bp)
    n_hi = sum(1 for b in bounds if b > hi_bp)
    seg_lo = plt.cm.Blues(np.linspace(0.08, 0.55, max(n_lo + 1, 2)))
    seg_mid = plt.cm.YlGn(np.linspace(0.18, 0.65, max(n_mid + 1, 2)))
    seg_hi = plt.cm.YlOrRd(np.linspace(0.40, 0.95, max(n_hi + 1, 2)))
    all_colors = (list(seg_lo[:n_lo + 1])
                  + list(seg_mid[1:n_mid + 1])
                  + list(seg_hi[1:n_hi + 1]))
    cmap = LinearSegmentedColormap.from_list("tricolor", all_colors, N=256)
    norm = BoundaryNorm(bounds, ncolors=256)
    return cmap, norm, bounds, vmax


def analyze(
    abundance_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    metadata_df: pd.DataFrame | None = None,
    params: dict | None = None,
) -> AnalysisResult:
    p = _normalize_deprecated_params({**DEFAULTS, **(params or {})})

    # ── 丰度矩阵 ────────────────────────────────────────────
    ab = abundance_df.copy()
    ab = ab.rename(columns={_mc.mag_col(ab): "MAG"})
    ab = ab[ab["MAG"].astype(str) != "unmapped"]
    sample_cols = [c for c in ab.columns if c != "MAG"]
    if not sample_cols:
        raise ValueError("abundance 表缺少样本列")
    for c in sample_cols:
        ab[c] = pd.to_numeric(ab[c], errors="coerce")
    ab = ab.dropna(subset=sample_cols, how="all").fillna(0.0)
    ab["MAG"] = ab["MAG"].astype(str)
    ab = ab.drop_duplicates(subset="MAG").reset_index(drop=True)
    if ab.empty:
        raise ValueError("abundance 表无有效 MAG")

    mat_full = ab[sample_cols].to_numpy(dtype=float)
    df = ab.copy()
    df["selection_score"] = _mc._select_score(mat_full, p["top_n_by"])
    df["abundance_mean"] = mat_full.mean(axis=1)

    # 注释 Phylum / Genus / Species / label
    df = _mc.annotate_taxonomy(df, taxonomy_df)
    df = _mc.annotate_keystone(df, keystone_df)

    # ── metadata → 样本按组排序 ──────────────────────────────
    groups_per_sample: dict[str, str | None] = {s: None for s in sample_cols}
    if metadata_df is not None and not metadata_df.empty:
        md = metadata_df.copy()
        sid_col = next((c for c in md.columns
                        if c.lower() in ("sampleid", "sample", "sample_id")),
                       md.columns[0])
        grp_col = next((c for c in md.columns if c.lower() == "group"), None)
        if grp_col is not None:
            groups_per_sample = dict(zip(
                md[sid_col].astype(str), md[grp_col].astype(str)))
            by_grp: dict[str, list[str]] = {}
            for s in sample_cols:
                g = groups_per_sample.get(s)
                by_grp.setdefault(g if g is not None else "", []).append(s)
            new_order: list[str] = []
            known = [g for g in ("CK", "A", "B") if g in by_grp]
            extra = [g for g in by_grp if g not in known]
            for g in known + extra:
                new_order.extend(by_grp[g])
            if new_order and new_order != sample_cols:
                sample_cols = new_order

    # ── Layer 1 — filter_mode ────────────────────────────────
    df = _mc.apply_filter_mode(
        df,
        mode=p["filter_mode"],
        top_n_count=int(p["top_n_count"]),
        top_n_by=p["top_n_by"],
        score_matrix=df[sample_cols].to_numpy(dtype=float),
    )

    # ── Layer 3 — 行排序 ────────────────────────────────────
    df = _mc.order_rows(
        df,
        mode=p["row_order"],
        metric_col="selection_score",
        abundance_col="abundance_mean",
        cluster_matrix=df[sample_cols].to_numpy(dtype=float),
        linkage_method=p["linkage_method"],
        log_transform_cluster=p["log_transform"],
    )
    df["row_order"] = range(len(df))

    # 列聚类（可选）
    mat = df[sample_cols].to_numpy(dtype=float)
    if p["cluster_cols"]:
        mat_T = np.log1p(mat.T) if p["log_transform"] else mat.T
        col_leaf = _mc._cluster_order(mat_T, p["linkage_method"])
        sample_cols = [sample_cols[k] for k in col_leaf]
        mat = df[sample_cols].to_numpy(dtype=float)

    # max_mags 硬截断
    if p.get("max_mags"):
        df = df.head(int(p["max_mags"])).reset_index(drop=True)
        mat = df[sample_cols].to_numpy(dtype=float)

    fig = _draw(df, sample_cols, mat, groups_per_sample, p)

    stats_df = df[["MAG", "label", "Phylum", "Genus", "Species",
                   "is_keystone", "row_order", "selection_score",
                   "abundance_mean"] + sample_cols].copy()
    return AnalysisResult(figure=fig, stats=stats_df, params=p)


def _draw(df, sample_cols, mat, groups_per_sample, p) -> plt.Figure:
    n_mag, n_sample = mat.shape
    cmap, norm, bounds, vmax = _build_tricolor_cmap(
        mat, tuple(p["color_breakpoints"]))

    fig = plt.figure(
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4),
        constrained_layout=True,
    )
    show_phy = p["show_phylum_bar"]
    show_leg = p["show_phylum_legend"]
    wr_phy = 0.05 if show_phy else 0.0001
    wr_leg = 0.35 if show_leg else 0.0001
    gs = fig.add_gridspec(1, 3, width_ratios=[wr_phy, 1.0, wr_leg], wspace=0.02)
    ax_phy = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    ax_leg = fig.add_subplot(gs[0, 2])
    ax_leg.axis("off")

    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm,
                   interpolation="nearest")

    # 顶部组彩条
    has_group_info = any(groups_per_sample.values())
    if p["show_group_bar"] and has_group_info:
        for j, s in enumerate(sample_cols):
            g = groups_per_sample.get(s)
            color = GROUP_COLORS.get(g, "#CCCCCC")
            ax.add_patch(Rectangle(
                (j - 0.5, -1.6), 1.0, 0.9,
                facecolor=color, edgecolor="none",
                clip_on=False, alpha=0.9, zorder=3,
            ))
        prev = None
        for j, s in enumerate(sample_cols):
            g = groups_per_sample.get(s)
            if prev is not None and g != prev:
                ax.axvline(j - 0.5, color="black", lw=1.2, zorder=4)
            prev = g

    # 左侧门彩条
    if show_phy:
        _mc.draw_phylum_bar(ax_phy, df["Phylum"].tolist())
    else:
        ax_phy.axis("off")

    # 行标签
    labels = []
    for _, r in df.iterrows():
        prefix = "★ " if (p["highlight_keystones"] and r["is_keystone"]) else ""
        labels.append(prefix + str(r.get("label", r["MAG"])))
    ax.set_yticks(range(n_mag))
    ax.set_yticklabels(labels,
                       fontsize=max(4, min(8, 260 // max(n_mag, 1))),
                       fontstyle="italic")

    ax.set_xticks(range(n_sample))
    ax.set_xticklabels(sample_cols, rotation=45, ha="right",
                       fontsize=max(5, min(9, 280 // max(n_sample, 1))))

    # 门边界横线
    prev_phy = None
    for i, phy in enumerate(df["Phylum"].tolist()):
        if prev_phy is not None and phy != prev_phy:
            ax.axhline(i - 0.5, color="#666", lw=0.5, zorder=3)
        prev_phy = phy

    ax.set_title(
        f"Top {n_mag} MAG Abundance  ({n_sample} samples,"
        f" color = relative abundance %)",
        fontsize=10, fontweight="bold", pad=12,
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # 右侧图例：Phylum + Group + ★ + colorbar
    if show_leg:
        _mc.draw_phylum_legend(ax_leg, df["Phylum"].tolist())
        # Group 图例
        if p["show_group_bar"] and has_group_info:
            import matplotlib.patches as mpatches
            grp_seen = [g for g in ("CK", "A", "B")
                        if g in set(groups_per_sample.values())]
            extra_g = sorted({g for g in groups_per_sample.values()
                              if g is not None and g not in ("CK", "A", "B")})
            grp_handles = [
                mpatches.Patch(facecolor=GROUP_COLORS.get(g, "#CCC"),
                               edgecolor="none", label=g)
                for g in (grp_seen + extra_g)
            ]
            if grp_handles:
                leg2 = ax_leg.legend(
                    handles=grp_handles,
                    loc="center left", bbox_to_anchor=(0.0, 0.38),
                    fontsize=6, title="Group", title_fontsize=7,
                    frameon=False, handlelength=1.2, handleheight=0.8,
                    labelspacing=0.4,
                )
                ax_leg.add_artist(leg2)
        if p["highlight_keystones"] and df["is_keystone"].any():
            _mc.draw_keystone_note(ax_leg)

    cax = ax_leg.inset_axes([0.1, 0.0, 0.22, 0.18])
    cbar = fig.colorbar(im, cax=cax, orientation="vertical")
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label("Rel. abundance (%)", fontsize=6.5)
    return fig
