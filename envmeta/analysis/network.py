"""共现网络辅助 — Degree vs Betweenness 散点图 + Gephi 预处理。

EnvMeta 不画共现网络图本体（Gephi 更专业），而是提供：
1. **Degree vs Betweenness 散点图**（keystone 筛选标准可视化）
2. **Gephi 预处理**（标签清理 / 格式校验）→ 交给 gephi_prep.py 处理
3. **Gephi 推荐参数指南**（写进 app expander）

输入：
    nodes_df: 节点表（Id + Degree + Betweenness + Phylum 等，Gephi nodes CSV 格式）
    edges_df: 边表（Source + Target + Weight + Spearman_r + p_value 等）
    taxonomy_df: 可选 GTDB 分类（补 Genus 标签）
    keystone_df: 可选 keystone 列表

输出：
    figure — Degree vs Betweenness 散点图
    stats  — 节点表（Id / label / Degree / Betweenness / Phylum / is_keystone / Module）
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from envmeta.analysis import _mag_common as _mc
from envmeta.analysis.base import AnalysisResult


DEFAULTS = {
    # Layer 1 — MAG 子集过滤（共享；网络默认 all 保持拓扑完整）
    "filter_mode": "all",
    "top_n_by": "mean",
    "top_n_count": 30,
    "max_mags": 0,
    # Layer 2 — 视觉
    "highlight_keystones": True,
    "show_phylum_legend": True,
    "width_mm": 260,
    "height_mm": 200,
    # Layer 4 — 网络特有
    "degree_threshold": 10,
    "betweenness_threshold": 200,
    "node_color_by": "keystone",       # "keystone" | "phylum"
    "color_keystone": "#1B3A5C",
    "color_non_keystone": "#B0D4E8",
}


def analyze(
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    params: dict | None = None,
) -> AnalysisResult:
    p = {**DEFAULTS, **(params or {})}

    # ── 节点表规范化 ────────────────────────────────────────
    df = nodes_df.copy()
    # Id 列嗅探
    id_col = next((c for c in df.columns if c.lower() in ("id", "mag", "name", "genome")),
                  df.columns[0])
    df = df.rename(columns={id_col: "MAG"})
    df["MAG"] = df["MAG"].astype(str)

    # Degree / Betweenness 列嗅探
    deg_col = next((c for c in df.columns if c.lower() == "degree"), None)
    bet_col = next((c for c in df.columns if c.lower() in ("betweenness", "betweenness centrality")), None)
    if deg_col is None or bet_col is None:
        raise ValueError("nodes 表需要 Degree + Betweenness 列（可在 Gephi 统计面板计算）")
    df["Degree"] = pd.to_numeric(df[deg_col], errors="coerce").fillna(0).astype(int)
    df["Betweenness"] = pd.to_numeric(df[bet_col], errors="coerce").fillna(0.0)

    # Module（可选）
    mod_col = next((c for c in df.columns if c.lower() in ("module", "modularity_class")), None)
    if mod_col:
        df["Module"] = df[mod_col]
    else:
        df["Module"] = 0

    # Phylum / Genus / Species / label（统一规则）
    df = _mc.annotate_taxonomy(df, taxonomy_df)
    df = _mc.annotate_keystone(df, keystone_df)

    # 若 nodes CSV 自带 is_keystone 列且没传 keystone_df，用 CSV 里的
    if keystone_df is None and "is_keystone" in nodes_df.columns:
        ks_col = nodes_df["is_keystone"]
        if ks_col.dtype == object:
            ks_col = ks_col.str.lower().isin(("true", "1", "yes"))
        df["is_keystone"] = ks_col.values[:len(df)]

    # Layer 1 — filter_mode
    df = _mc.apply_filter_mode(
        df, mode=p["filter_mode"],
        top_n_count=int(p["top_n_count"]),
        top_n_by=p["top_n_by"],
        score_col="Degree",
    )

    # max_mags
    if p.get("max_mags"):
        df = df.head(int(p["max_mags"])).reset_index(drop=True)

    # ── 绘图 ────────────────────────────────────────────────
    fig = _draw_scatter(df, p)

    # ── stats ────────────────────────────────────────────────
    stats_cols = ["MAG", "label", "Degree", "Betweenness", "Phylum",
                  "Family", "Genus", "is_keystone", "Module"]
    stats_cols = [c for c in stats_cols if c in df.columns]
    stats_df = df[stats_cols].copy()

    return AnalysisResult(figure=fig, stats=stats_df, params=p)


def _draw_scatter(df, p) -> plt.Figure:
    """Degree vs Betweenness 散点图（keystone 筛选标准可视化）。"""
    fig = plt.figure(
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4),
        constrained_layout=True,
    )
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.25], wspace=0.02)
    ax = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])
    ax_leg.axis("off")

    color_by = p["node_color_by"]

    # 先画非 keystone（zorder 低）→ 再画 keystone（zorder 高）
    for is_ks in [False, True]:
        sub = df[df["is_keystone"] == is_ks]
        if sub.empty:
            continue
        if color_by == "phylum":
            colors = [_mc.PHYLUM_COLORS.get(phy, "#888") for phy in sub["Phylum"]]
        else:
            colors = [p["color_keystone"] if is_ks else p["color_non_keystone"]]
            colors = colors * len(sub)
        sizes = np.where(is_ks, 80, 30)
        alpha = 0.9 if is_ks else 0.5
        zorder = 4 if is_ks else 2
        ax.scatter(sub["Degree"], sub["Betweenness"],
                   s=sizes, c=colors, alpha=alpha,
                   edgecolors="white", linewidths=0.5, zorder=zorder)

        # keystone 标注 genus label
        if is_ks and p["highlight_keystones"]:
            for _, r in sub.iterrows():
                label = str(r.get("label", r["MAG"]))
                if len(label) > 18:
                    label = label[:16] + ".."
                ax.annotate(
                    label, (r["Degree"], r["Betweenness"]),
                    fontsize=6, fontstyle="italic", fontweight="bold",
                    xytext=(5, 5), textcoords="offset points",
                    color=p["color_keystone"], zorder=5,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                              edgecolor="none", alpha=0.7),
                )

    # 阈值线
    deg_thr = p["degree_threshold"]
    bet_thr = p["betweenness_threshold"]
    ax.axvline(x=deg_thr, color="#E74C3C", ls="--", alpha=0.4, lw=1)
    ax.axhline(y=bet_thr, color="#E74C3C", ls="--", alpha=0.4, lw=1)
    ax.text(0.97, 0.97,
            f"Degree >= {deg_thr} or\nBetweenness >= {bet_thr}",
            transform=ax.transAxes, fontsize=7, va="top", ha="right",
            color="#E74C3C", alpha=0.7)

    ax.set_xlabel("Degree", fontsize=10)
    ax.set_ylabel("Betweenness Centrality", fontsize=10)
    ax.tick_params(labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("#ccc")

    n_nodes = len(df)
    n_ks = int(df["is_keystone"].sum())
    ax.set_title(
        f"Keystone Species Selection  ({n_nodes} nodes, {n_ks} keystone)",
        fontsize=11, fontweight="bold", pad=10,
    )

    # 右侧图例
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=p["color_keystone"], markersize=9,
               markeredgecolor="white", markeredgewidth=0.5,
               label="Keystone"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=p["color_non_keystone"], markersize=6,
               markeredgecolor="white", markeredgewidth=0.5,
               label="Non-keystone"),
    ]
    leg = ax_leg.legend(handles=handles, loc="upper left",
                        bbox_to_anchor=(0.0, 1.0),
                        fontsize=7, frameon=False,
                        title="Node type", title_fontsize=8)
    ax_leg.add_artist(leg)

    # Phylum 图例（color_by=phylum 时才画）
    if color_by == "phylum" and p["show_phylum_legend"]:
        _mc.draw_phylum_legend(ax_leg, df["Phylum"].tolist(),
                               loc="center left",
                               bbox_to_anchor=(0.0, 0.5))

    # 网络统计
    ax_leg.text(0.0, 0.15,
                f"n = {n_nodes}\nkeystone = {n_ks}\n"
                f"avg degree = {df['Degree'].mean():.1f}",
                transform=ax_leg.transAxes,
                fontsize=6.5, va="top", ha="left",
                color="#666")

    return fig
