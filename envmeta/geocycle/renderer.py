"""CycleData → matplotlib 静态渲染（v1）。

布局：2×2 元素象限（As/N/S/Fe），每象限列出活跃通路（横向柱）+ 前几名 MAG；
若有环境相关性，右侧追加环境-通路连线摘要条。
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

from envmeta.geocycle.model import CycleData, ElementCycle

DEFAULTS = {
    "width_mm": 360,
    "height_mm": 260,
    "show_env_panel": True,
    "max_pathways_per_element": 8,
    "max_contributors_shown": 3,
    "title": "Biogeochemical Cycle Diagram (v1)",
}


def _draw_element_quadrant(ax, ec: ElementCycle, cfg: dict) -> None:
    """在给定 ax 上画一个元素的通路+贡献 MAG 布局。"""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("#FAFAFA")
    for spine in ax.spines.values():
        spine.set_edgecolor(ec.color)
        spine.set_linewidth(1.5)

    # 标题
    ax.text(
        0.02, 0.98,
        f"{ec.display_name}  ({ec.n_active_pathways}/{len(ec.pathways)} active)",
        transform=ax.transAxes, fontsize=11, fontweight="bold",
        color=ec.color, va="top", ha="left",
    )

    # 活跃通路
    active = [p for p in ec.pathways if p.n_active_mags > 0]
    active = active[: cfg["max_pathways_per_element"]]
    if not active:
        ax.text(0.5, 0.5, "No active pathway",
                transform=ax.transAxes, fontsize=10, color="#999",
                ha="center", va="center")
        return

    max_contrib = max(p.total_contribution for p in active) or 1.0
    n = len(active)
    row_h = 0.85 / max(n, 1)
    y_top = 0.90
    bar_x_start = 0.03
    bar_x_end = 0.60
    for i, pw in enumerate(active):
        y = y_top - (i + 1) * row_h
        width = (pw.total_contribution / max_contrib) * (bar_x_end - bar_x_start)
        # 柱
        ax.add_patch(Rectangle(
            (bar_x_start, y + 0.02), width, row_h * 0.55,
            transform=ax.transAxes, color=ec.color, alpha=0.85, zorder=2,
        ))
        # 通路名
        ax.text(bar_x_start, y + 0.02 + row_h * 0.58,
                f"{pw.display_name}  (n={pw.n_active_mags}, c̄={pw.mean_completeness:.0f}%)",
                transform=ax.transAxes, fontsize=8, color="#222", va="bottom")
        # 反应式
        if pw.reaction:
            ax.text(bar_x_start, y + 0.005,
                    pw.reaction, transform=ax.transAxes,
                    fontsize=7, color="#555", va="bottom", fontstyle="italic")
        # Top 贡献 MAG
        labels = [c.label for c in pw.contributors[: cfg["max_contributors_shown"]]]
        if labels:
            ax.text(bar_x_end + 0.02, y + 0.02 + row_h * 0.25,
                    "→ " + ", ".join(labels),
                    transform=ax.transAxes, fontsize=7.5, color="#333",
                    va="center", ha="left")


def _draw_env_panel(fig, data: CycleData, bbox) -> None:
    """底部水平条：显著的 env-pathway 相关性。"""
    ax = fig.add_axes(bbox)
    ax.axis("off")
    if not data.env_correlations:
        ax.text(0.5, 0.5, "No env-pathway correlation above threshold",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=9, color="#888")
        return
    # 按 |rho| 降序
    sorted_corrs = sorted(data.env_correlations,
                          key=lambda c: abs(c.rho), reverse=True)[:12]
    n = len(sorted_corrs)
    row_h = 0.8 / max(n, 1)
    ax.text(0.01, 0.96, f"Environmental coupling (|ρ| ≥ {data.params.get('env_rho_min', 0.5)}, p < {data.params.get('env_p_max', 0.05)})",
            transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
    # 可信度标签颜色
    conf_colors = {
        "strong":     "#27AE60",
        "suggestive": "#F39C12",
        "weak":       "#95A5A6",
        "spurious?":  "#E74C3C",
        "none":       "#CCCCCC",
        "unknown":    "#CCCCCC",
    }
    for i, ec in enumerate(sorted_corrs):
        y = 0.88 - (i + 1) * row_h
        color = "#C0392B" if ec.rho > 0 else "#2980B9"
        bar_w = abs(ec.rho) * 0.3
        ax.add_patch(Rectangle((0.35, y), bar_w, row_h * 0.7,
                               transform=ax.transAxes, color=color, alpha=0.8))
        pp = ec.perm_p if ec.perm_p is not None else ec.p_value
        sig = "***" if pp < 0.001 else "**" if pp < 0.01 else "*" if pp < 0.05 else ""
        ax.text(0.01, y + row_h * 0.3, ec.pathway_id,
                transform=ax.transAxes, fontsize=8, va="center")
        ax.text(0.22, y + row_h * 0.3, ec.env_factor,
                transform=ax.transAxes, fontsize=8, va="center",
                fontweight="bold")
        ax.text(0.35 + bar_w + 0.01, y + row_h * 0.3,
                f"ρ={ec.rho:.2f} permP={pp:.3f} {sig}",
                transform=ax.transAxes, fontsize=7.5, va="center",
                color=color)
        # 可信度徽章
        conf = ec.confidence or "unknown"
        ax.text(0.91, y + row_h * 0.3, conf,
                transform=ax.transAxes, fontsize=7,
                fontweight="bold", color="white",
                ha="center", va="center",
                bbox=dict(facecolor=conf_colors.get(conf, "#888"),
                          alpha=0.95, pad=2,
                          boxstyle="round,pad=0.25", edgecolor="none"))


def render(data: CycleData, params: dict | None = None) -> plt.Figure:
    cfg = {**DEFAULTS, **(params or {})}

    fig = plt.figure(
        figsize=(cfg["width_mm"] / 25.4, cfg["height_mm"] / 25.4),
    )
    # 主标题
    fig.suptitle(
        cfg["title"],
        fontsize=13, fontweight="bold", y=0.985,
    )
    subtitle = (
        f"{data.meta.get('n_mags', 0)} MAGs  |  "
        f"{data.meta.get('n_pathways_active', 0)} / "
        f"{data.meta.get('n_pathways_total', 0)} pathways active  |  "
        f"{data.meta.get('n_env_correlations', 0)} env couplings"
    )
    fig.text(0.5, 0.955, subtitle, ha="center", fontsize=9, color="#555")

    has_env = cfg["show_env_panel"] and data.env_correlations is not None
    # 主体：2×2 象限；若 has_env，底部留 env 面板
    left, right = 0.04, 0.98
    top, bottom = 0.93, 0.27 if has_env else 0.05
    qw = (right - left) / 2 - 0.02
    qh = (top - bottom) / 2 - 0.02

    elements_ordered = sorted(
        data.elements,
        key=lambda e: {"arsenic": 0, "nitrogen": 1, "sulfur": 2, "iron": 3}.get(
            e.element_id, 99),
    )
    positions = [
        (left,         top - qh),       # 左上 As
        (left + qw + 0.04, top - qh),   # 右上 N
        (left,         bottom),         # 左下 S
        (left + qw + 0.04, bottom),     # 右下 Fe
    ]
    for i, ec in enumerate(elements_ordered[:4]):
        x, y = positions[i]
        ax = fig.add_axes([x, y, qw, qh])
        _draw_element_quadrant(ax, ec, cfg)

    if has_env:
        _draw_env_panel(fig, data, bbox=[left, 0.03, right - left, 0.20])

    return fig
