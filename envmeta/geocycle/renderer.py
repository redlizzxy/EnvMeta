"""CycleData → matplotlib 渲染（v2 — Mockup 10 级联细胞）。

布局：2×2 元素象限（As/N/S/Fe）；每象限画最多 ``max_cells_per_element`` 个
"合并细胞"，每个细胞 = 一个 MAG 在一条通路上的完整基因级联（底物→酶→
中间产物→…→最终产物）。底部保留 env-pathway 相关性摘要。

相较 v1（2×2 横向柱图）：
  - 细胞替代抽象柱
  - 基因级联显示化学物中间产物（KB substrate/product 字段驱动）
  - 细胞外左=起始底物、右=最终产物，箭头穿膜

传入 ``cell_mode="bars"`` 可回退 v1 样式。
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from envmeta.geocycle.cell_renderer import draw_cascade_cell, genes_to_steps
from envmeta.geocycle.model import CycleData, ElementCycle, PathwayActivity

DEFAULTS = {
    "width_mm": 460,
    "height_mm": 320,
    "show_env_panel": True,
    "max_cells_per_element": 3,           # 每元素象限最多画几个细胞
    "max_pathways_per_element": 8,        # v1 回退模式用
    "max_contributors_shown": 3,          # v1 回退模式用
    "cell_mode": "cascade",               # "cascade"（v2，默认）/ "bars"（v1 回退）
    "title": "Biogeochemical Cycle Diagram (v2)",
    "cell_height_ratio": 0.22,            # 单细胞占象限高的比例
}


# =============================================================================
# v2 — 级联细胞象限
# =============================================================================

def _select_top_cells(
    ec: ElementCycle, max_cells: int,
) -> list[tuple[PathwayActivity, object]]:
    """为一个元素挑选至多 max_cells 个 (pathway, contributor) 组合。

    优先覆盖不同通路 + 优先高贡献度。返回 [(pathway, contributor), ...]。
    """
    # 按 pathway 排序（inference 已按 total_contribution 降序）
    picked: list[tuple[PathwayActivity, object]] = []
    for pw in ec.pathways:
        if len(picked) >= max_cells:
            break
        if pw.n_active_mags == 0 or not pw.contributors:
            continue
        top_mag = pw.contributors[0]
        if not top_mag.genes:
            continue
        picked.append((pw, top_mag))

    # 若活跃通路 < max_cells，从同一通路补第 2/3 contributor
    if len(picked) < max_cells:
        for pw in ec.pathways:
            if pw.n_active_mags == 0:
                continue
            for c in pw.contributors[1:]:
                if len(picked) >= max_cells:
                    break
                if not c.genes:
                    continue
                picked.append((pw, c))
            if len(picked) >= max_cells:
                break
    return picked


def _draw_element_quadrant_cascade(
    ax, ec: ElementCycle, cfg: dict,
) -> list[dict]:
    """v2 级联细胞布局。返回每个细胞绘制的锚点字典列表。"""
    # 数据坐标系：统一 18×10
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 18, 10, fill=False,
                           edgecolor=ec.color, linewidth=1.6, alpha=0.6))

    # 标题（元素名 + 活跃通路计数）
    ax.text(0.3, 9.55,
            f"{ec.display_name}  ({ec.n_active_pathways}/{len(ec.pathways)} active)",
            fontsize=13, fontweight="bold", color=ec.color, zorder=10)

    picked = _select_top_cells(ec, cfg["max_cells_per_element"])
    if not picked:
        ax.text(9, 5, "No active pathway", ha="center", va="center",
                fontsize=11, color="#999")
        return []

    n = len(picked)
    # 垂直均匀分布；顶部留 1.2 给标题
    top, bot = 8.4, 0.6
    cell_h_max = 1.4
    row_span = (top - bot) / n
    cell_h = min(cell_h_max, row_span * 0.75)
    ys = [top - (i + 0.5) * row_span for i in range(n)]

    anchors: list[dict] = []
    for (pw, contrib), cy in zip(picked, ys):
        n_genes = len(contrib.genes)
        # 细胞宽随基因数伸缩（2 基因≈5 宽；4 基因≈8.5 宽；上限 10）
        cell_w = float(np.clip(3.5 + n_genes * 1.4, 4.5, 10.0))
        cell_x0 = (18 - cell_w) / 2

        steps = genes_to_steps(
            contrib.genes,
            default_color=ec.color,
        )
        title = pw.display_name
        mag_label = contrib.label or contrib.mag
        r = draw_cascade_cell(
            ax,
            cy=cy, cell_x0=cell_x0, cell_w=cell_w, cell_h=cell_h,
            title=title, mag_label=mag_label,
            steps=steps,
            element_color=ec.color,
            show_heatmap=False,   # 上方 mini 热图留待 S2.5-3（需 KO-env 相关性数据）
            show_outside_chems=True,
        )
        # 小信息条：通路完整度 + 贡献度（右下）
        ax.text(
            17.7, cy - cell_h / 2 - 0.08,
            f"c={contrib.completeness:.0f}%  ab̄={contrib.abundance_mean:.2f}",
            fontsize=6.5, color="#555", ha="right", va="top",
        )
        r.update({
            "pathway_id": pw.pathway_id,
            "pathway_name": pw.display_name,
            "mag": contrib.mag,
            "element": ec.element_id,
            "ax": ax,
        })
        anchors.append(r)
    return anchors


# =============================================================================
# v1 回退 — 抽象柱图（保留以支持老 validation / 回退）
# =============================================================================

def _draw_element_quadrant_bars(ax, ec: ElementCycle, cfg: dict) -> None:
    """v1 布局：元素框 + 横向柱列出活跃通路 + Top 贡献 MAG。"""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("#FAFAFA")
    for spine in ax.spines.values():
        spine.set_edgecolor(ec.color)
        spine.set_linewidth(1.5)

    ax.text(
        0.02, 0.98,
        f"{ec.display_name}  ({ec.n_active_pathways}/{len(ec.pathways)} active)",
        transform=ax.transAxes, fontsize=11, fontweight="bold",
        color=ec.color, va="top", ha="left",
    )
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
        ax.add_patch(Rectangle(
            (bar_x_start, y + 0.02), width, row_h * 0.55,
            transform=ax.transAxes, color=ec.color, alpha=0.85, zorder=2,
        ))
        ax.text(bar_x_start, y + 0.02 + row_h * 0.58,
                f"{pw.display_name}  (n={pw.n_active_mags}, c̄={pw.mean_completeness:.0f}%)",
                transform=ax.transAxes, fontsize=8, color="#222", va="bottom")
        if pw.reaction:
            ax.text(bar_x_start, y + 0.005,
                    pw.reaction, transform=ax.transAxes,
                    fontsize=7, color="#555", va="bottom", fontstyle="italic")
        labels = [c.label for c in pw.contributors[: cfg["max_contributors_shown"]]]
        if labels:
            ax.text(bar_x_end + 0.02, y + 0.02 + row_h * 0.25,
                    "→ " + ", ".join(labels),
                    transform=ax.transAxes, fontsize=7.5, color="#333",
                    va="center", ha="left")


# =============================================================================
# env 面板（v1/v2 共用）
# =============================================================================

def _draw_env_panel(fig, data: CycleData, bbox) -> None:
    ax = fig.add_axes(bbox)
    ax.axis("off")
    if not data.env_correlations:
        ax.text(0.5, 0.5, "No env-pathway correlation above threshold",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=9, color="#888")
        return
    sorted_corrs = sorted(data.env_correlations,
                          key=lambda c: abs(c.rho), reverse=True)[:12]
    n = len(sorted_corrs)
    row_h = 0.8 / max(n, 1)
    ax.text(0.01, 0.96,
            f"Environmental coupling (|ρ| ≥ {data.params.get('env_rho_min', 0.5)}, "
            f"p < {data.params.get('env_p_max', 0.05)})",
            transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
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
        conf = ec.confidence or "unknown"
        ax.text(0.91, y + row_h * 0.3, conf,
                transform=ax.transAxes, fontsize=7,
                fontweight="bold", color="white",
                ha="center", va="center",
                bbox=dict(facecolor=conf_colors.get(conf, "#888"),
                          alpha=0.95, pad=2,
                          boxstyle="round,pad=0.25", edgecolor="none"))


# =============================================================================
# 主入口
# =============================================================================

def render(data: CycleData, params: dict | None = None) -> plt.Figure:
    cfg = {**DEFAULTS, **(params or {})}

    fig = plt.figure(
        figsize=(cfg["width_mm"] / 25.4, cfg["height_mm"] / 25.4),
    )
    fig.suptitle(cfg["title"],
                 fontsize=14, fontweight="bold", y=0.985)
    subtitle = (
        f"{data.meta.get('n_mags', 0)} MAGs  |  "
        f"{data.meta.get('n_pathways_active', 0)} / "
        f"{data.meta.get('n_pathways_total', 0)} pathways active  |  "
        f"{data.meta.get('n_env_correlations', 0)} env couplings"
    )
    fig.text(0.5, 0.958, subtitle, ha="center", fontsize=9, color="#555")

    has_env = cfg["show_env_panel"] and data.env_correlations is not None
    left, right = 0.04, 0.98
    top, bottom = 0.93, 0.27 if has_env else 0.05
    qw = (right - left) / 2 - 0.02
    qh = (top - bottom) / 2 - 0.02

    elements_ordered = sorted(
        data.elements,
        key=lambda e: {"arsenic": 0, "nitrogen": 1, "sulfur": 2, "iron": 3}
        .get(e.element_id, 99),
    )
    positions = [
        (left,             top - qh),      # 左上 As
        (left + qw + 0.04, top - qh),      # 右上 N
        (left,             bottom),        # 左下 S
        (left + qw + 0.04, bottom),        # 右下 Fe
    ]

    mode = cfg.get("cell_mode", "cascade")
    all_anchors: list[dict] = []
    for i, ec in enumerate(elements_ordered[:4]):
        x, y = positions[i]
        ax = fig.add_axes([x, y, qw, qh])
        if mode == "cascade":
            anchors = _draw_element_quadrant_cascade(ax, ec, cfg)
            all_anchors.extend(anchors)
        else:
            _draw_element_quadrant_bars(ax, ec, cfg)

    if has_env:
        _draw_env_panel(fig, data, bbox=[left, 0.03, right - left, 0.20])

    # 附加到 figure 便于上层（S2.5-3 化学物耦合）检索细胞锚点
    fig._envmeta_cycle_anchors = all_anchors   # type: ignore[attr-defined]
    return fig
