"""Mockup 06 & 07 — 数据驱动布局（无预设 FeOOH）。

设计原则：
- 布局由 KB 元素数量决定（1 元素=1 格，4 元素=2×2，6 元素=2×3）
- 每个元素象限内：化学物流程链（节点 + 箭头）+ 关键步骤嵌入式细胞
- 中心不强制特定物质，不同研究可复用同一布局
- 元素间耦合可通过"共享化学物节点"或"外部虚线耦合"展示（可选）

输出：
  06_single_group_neutral.png — 单组详细版（4 元素 2×2，化学链 + 细胞）
  07_three_groups_horizontal.png — 3 组左右并排（简化细胞）
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import (
    Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle, Polygon
)
import numpy as np
from pathlib import Path

OUT = Path(__file__).parent
plt.rcParams["font.family"] = "DejaVu Sans"

ELEM = {"As": "#E74C3C", "As2": "#C0392B",
        "N": "#3498DB", "N2": "#2874A6",
        "S": "#229954", "S2": "#1E8449",
        "Fe": "#CA8A1E", "Fe2": "#A77215"}
GROUP_COLOR = {"CK": "#1c9cbd", "A": "#e3943d", "B": "#92181e"}
ENV_LABELS = ["pH", "Eh", "TOC", "As"]


def mini_heatmap(ax, x, y, vals, size=0.09):
    for i, r in enumerate(vals):
        ax.add_patch(Rectangle(
            (x + i * size, y), size, size,
            facecolor=plt.cm.RdBu_r((r + 1) / 2),
            edgecolor="#333", linewidth=0.3, zorder=8,
        ))


def gene_oval(ax, x, y, name, color, corr=None, w=0.6, h=0.32, show_hm=True):
    ax.add_patch(Ellipse((x, y), w, h, facecolor=color,
                         edgecolor="#222", linewidth=0.9, zorder=6))
    ax.text(x, y, name, ha="center", va="center",
            fontsize=6.5, fontweight="bold", color="white", zorder=7)
    if show_hm and corr:
        hm_size = 0.1
        mini_heatmap(ax, x - 2 * hm_size, y + h / 2 + 0.05, corr, hm_size)


def draw_cell(ax, cx, cy, w, h, name, mag, genes, show_hm=True):
    """细胞代谢小图。"""
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.08",
        facecolor="#FDF6E3", edgecolor="#6B4F2A",
        linewidth=1.3, zorder=4,
    ))
    ax.text(cx, cy + h / 2 - 0.18, name, ha="center", va="top",
            fontsize=8, fontweight="bold", color="#222", zorder=7)
    ax.text(cx, cy + h / 2 - 0.42, mag, ha="center", va="top",
            fontsize=6, fontstyle="italic", color="#555", zorder=7)
    n = len(genes)
    gy = cy - 0.08
    if n:
        xs = np.linspace(cx - (w - 0.9) / 2, cx + (w - 0.9) / 2, n)
        for xg, (g, c, corr) in zip(xs, genes):
            gene_oval(ax, xg, gy, g, c, corr, show_hm=show_hm)


def chem_node(ax, x, y, txt, color, size=0.55):
    """化学物节点（圆形）。"""
    ax.add_patch(Ellipse((x, y), size * 1.5, size, facecolor="white",
                         edgecolor=color, linewidth=1.5, zorder=5))
    ax.text(x, y, txt, ha="center", va="center", fontsize=8.5,
            fontweight="bold", color="#222", zorder=6)


def arrow(ax, x0, y0, x1, y1, color="#555", lw=1.2, style="-|>",
          alpha=0.85, rad=0.0):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle=style, mutation_scale=14,
        linewidth=lw, color=color, alpha=alpha, zorder=3,
        connectionstyle=f"arc3,rad={rad}",
    ))


# ======================================================================
# 元素象限绘制器（复用）
# ======================================================================
def draw_as_quadrant(ax, detail=True, group_widths=None):
    """As 象限：化学链 As(V) ↔ As(III) → MAs(III) + 嵌入式细胞。"""
    # 化学物节点（布局：横向排布）
    chem_node(ax, 1.2, 3.5, "As(V)", ELEM["As"])
    chem_node(ax, 4.0, 3.5, "As(III)", ELEM["As"])
    chem_node(ax, 6.8, 3.5, "MAs(III)", ELEM["As2"])

    w = group_widths or {"ox": 2.0, "red": 2.0, "me": 1.5, "resp": 1.5}

    # 箭头（化学转化）
    arrow(ax, 1.9, 3.7, 3.3, 3.7, ELEM["As"], lw=w["ox"], rad=0.0)    # ox
    arrow(ax, 3.3, 3.3, 1.9, 3.3, ELEM["As"], lw=w["red"], rad=0.0)   # red
    arrow(ax, 4.7, 3.5, 6.1, 3.5, ELEM["As2"], lw=w["me"], rad=0.0)   # methylation

    # 箭头旁标签：通路名 + Top MAG
    ax.text(2.6, 4.0, "Arsenite ox.", fontsize=6.5, ha="center",
            color=ELEM["As"], fontweight="bold")
    ax.text(2.6, 4.25, "→ Gallionella", fontsize=6, ha="center",
            fontstyle="italic", color="#333")
    ax.text(2.6, 3.0, "Arsenate red.", fontsize=6.5, ha="center",
            color=ELEM["As"], fontweight="bold")
    ax.text(2.6, 2.75, "→ Gallionella", fontsize=6, ha="center",
            fontstyle="italic", color="#333")
    ax.text(5.4, 3.8, "As methylation", fontsize=6.5, ha="center",
            color=ELEM["As2"], fontweight="bold")
    ax.text(5.4, 4.05, "→ Fen-1038", fontsize=6, ha="center",
            fontstyle="italic", color="#333")

    if detail:
        # 嵌入细胞：Arsenite oxidizer + Arsenate reducer（2 个小细胞）
        draw_cell(ax, cx=2.0, cy=1.5, w=1.7, h=0.95,
                  name="NAOB", mag="Gallionella",
                  genes=[("aioA", ELEM["As"], [0.6, 0.2, -0.1, 0.8]),
                         ("aioB", ELEM["As"], [0.55, 0.15, -0.15, 0.75])])
        draw_cell(ax, cx=4.0, cy=1.5, w=1.7, h=0.95,
                  name="DARPs", mag="Desulfo.",
                  genes=[("arrA", ELEM["As2"], [0.8, 0.1, -0.2, 0.9]),
                         ("arsC", ELEM["As"], [0.5, 0.2, -0.1, 0.7])])
        draw_cell(ax, cx=6.0, cy=1.5, w=1.7, h=0.95,
                  name="As-methylator", mag="Fen-1038",
                  genes=[("arsM", ELEM["As2"], [0.4, 0.3, 0.1, 0.65])])

    # 标题
    ax.text(0.2, 4.6, "As cycle", fontsize=11, fontweight="bold",
            color=ELEM["As"])


def draw_n_quadrant(ax, detail=True, group_widths=None):
    """N 象限：NH4 ← N2 ↓ NO3 → NO2 → NO → N2O → N2"""
    chem_node(ax, 0.8, 3.5, "NH₄⁺", ELEM["N"])
    chem_node(ax, 1.6, 2.0, "N₂", ELEM["N"])
    chem_node(ax, 2.5, 3.5, "NO₃⁻", ELEM["N"])
    chem_node(ax, 4.0, 3.5, "NO₂⁻", ELEM["N"])
    chem_node(ax, 5.5, 3.5, "NO", ELEM["N"])
    chem_node(ax, 6.5, 2.5, "N₂O", ELEM["N"])
    chem_node(ax, 7.4, 1.5, "N₂", ELEM["N"])

    w = group_widths or {"nar": 2.0, "nir": 2.0, "nor": 1.8,
                          "nos": 1.5, "amo": 1.0, "nif": 1.5}

    arrow(ax, 3.1, 3.5, 3.4, 3.5, ELEM["N"], lw=w["nar"])    # NO3→NO2
    arrow(ax, 4.6, 3.5, 5.0, 3.5, ELEM["N"], lw=w["nir"])    # NO2→NO
    arrow(ax, 5.8, 3.3, 6.2, 2.8, ELEM["N"], lw=w["nor"])    # NO→N2O
    arrow(ax, 6.8, 2.2, 7.1, 1.8, ELEM["N"], lw=w["nos"])    # N2O→N2
    arrow(ax, 1.4, 3.3, 1.0, 3.3, ELEM["N"], lw=w["amo"], rad=-0.2)  # NH4→NO2 (AMO)
    arrow(ax, 2.0, 2.3, 0.9, 3.2, ELEM["N"], lw=w["nif"], rad=0.2)   # N2→NH4 (Nif)

    # 通路名简化（省地方，只标主要）
    ax.text(3.2, 3.8, "nar", fontsize=6.5, color=ELEM["N"])
    ax.text(4.7, 3.8, "nir", fontsize=6.5, color=ELEM["N"])
    ax.text(5.8, 3.0, "nor", fontsize=6.5, color=ELEM["N"])
    ax.text(7.0, 2.0, "nos", fontsize=6.5, color=ELEM["N"])

    if detail:
        draw_cell(ax, cx=3.3, cy=1.2, w=1.7, h=0.95,
                  name="Denitrifier", mag="Gallionella",
                  genes=[("narG", ELEM["N"], [-0.3, 0.4, 0.6, 0.5]),
                         ("nirK", ELEM["N2"], [0.1, 0.7, -0.1, 0.5])])
        draw_cell(ax, cx=5.3, cy=1.2, w=1.7, h=0.95,
                  name="NO/N2O-red.", mag="Gallionella",
                  genes=[("norB", ELEM["N"], [0.0, 0.5, 0.1, 0.4]),
                         ("nosZ", ELEM["N"], [0.3, 0.4, 0.2, 0.6])])
        draw_cell(ax, cx=1.3, cy=1.0, w=1.7, h=0.95,
                  name="N-fixer", mag="JACRMN01",
                  genes=[("nifD", ELEM["N2"], [0.5, 0.76, -0.3, 0.2]),
                         ("nifK", ELEM["N2"], [0.4, 0.72, -0.25, 0.2])])

    ax.text(0.2, 4.6, "N cycle", fontsize=11, fontweight="bold",
            color=ELEM["N"])


def draw_s_quadrant(ax, detail=True, group_widths=None):
    """S 象限：SO4 ↔ H2S + S2O3 + Sulfide ox."""
    chem_node(ax, 1.5, 3.5, "SO₄²⁻", ELEM["S"])
    chem_node(ax, 4.5, 3.5, "S⁰", ELEM["S"])
    chem_node(ax, 7.0, 3.5, "H₂S", ELEM["S2"])
    chem_node(ax, 4.5, 1.8, "S₂O₃²⁻", ELEM["S"])

    w = group_widths or {"dsr": 2.0, "sox": 2.0, "sqr": 1.5, "tsd": 1.0}

    arrow(ax, 2.2, 3.3, 6.3, 3.3, ELEM["S"], lw=w["dsr"], rad=0.2)   # SO4→H2S
    arrow(ax, 6.3, 3.7, 2.2, 3.7, ELEM["S"], lw=w["sox"], rad=-0.2)  # H2S→SO4 (Sox)
    arrow(ax, 4.5, 3.0, 4.5, 2.3, ELEM["S"], lw=w["tsd"], rad=0.0)   # S0↔S2O3

    ax.text(4.2, 4.1, "Sulfide oxidation (sox)", fontsize=6.5, ha="center",
            color=ELEM["S"], fontweight="bold")
    ax.text(4.2, 4.35, "→ Sulfuricaulis", fontsize=6, ha="center",
            fontstyle="italic", color="#333")
    ax.text(4.2, 2.8, "Dissim. sulfate red. (dsr)", fontsize=6.5, ha="center",
            color=ELEM["S"], fontweight="bold")
    ax.text(4.2, 2.55, "→ Sulfuricaulis", fontsize=6, ha="center",
            fontstyle="italic", color="#333")

    if detail:
        draw_cell(ax, cx=2.5, cy=0.9, w=1.7, h=0.95,
                  name="SRB", mag="Sulfuricaulis",
                  genes=[("dsrA", ELEM["S"], [0.3, -0.2, 0.1, 0.5]),
                         ("aprA", ELEM["S2"], [0.2, -0.1, 0.2, 0.4])])
        draw_cell(ax, cx=5.5, cy=0.9, w=1.7, h=0.95,
                  name="SOB", mag="Thiobacillus",
                  genes=[("soxA", ELEM["S"], [0.1, 0.3, -0.2, 0.6]),
                         ("soxB", ELEM["S"], [0.0, 0.4, -0.2, 0.6])])

    ax.text(0.2, 4.6, "S cycle", fontsize=11, fontweight="bold",
            color=ELEM["S"])


def draw_fe_quadrant(ax, detail=True, group_widths=None):
    """Fe 象限：Fe(II) ↔ Fe(III) + 铁氧还原/转运。"""
    chem_node(ax, 2.0, 3.8, "Fe(III)", ELEM["Fe"])
    chem_node(ax, 6.0, 3.8, "Fe(II)", ELEM["Fe"])
    chem_node(ax, 4.0, 2.0, "Fe-organic", ELEM["Fe2"])

    w = group_widths or {"dirb": 2.0, "fox": 1.5, "uptake": 2.0}

    arrow(ax, 2.8, 3.6, 5.2, 3.6, ELEM["Fe"], lw=w["dirb"])    # Fe(III)→Fe(II) (DIRB)
    arrow(ax, 5.2, 4.0, 2.8, 4.0, ELEM["Fe"], lw=w["fox"])     # Fe(II)→Fe(III) (FeOB)
    arrow(ax, 5.8, 3.4, 4.5, 2.3, ELEM["Fe2"], lw=w["uptake"]) # uptake

    ax.text(4.0, 4.4, "Fe(III) reduction", fontsize=6.5, ha="center",
            color=ELEM["Fe"], fontweight="bold")
    ax.text(4.0, 4.65, "→ Geobacter", fontsize=6, ha="center",
            fontstyle="italic", color="#333")
    ax.text(4.0, 3.2, "Fe(II) oxidation", fontsize=6.5, ha="center",
            color=ELEM["Fe"], fontweight="bold")
    ax.text(4.0, 3.4, "→ Sideroxydans", fontsize=6, ha="center",
            fontstyle="italic", color="#333")

    if detail:
        draw_cell(ax, cx=2.2, cy=0.7, w=1.7, h=0.95,
                  name="DIRB", mag="Geobacter",
                  genes=[("mtrA", ELEM["Fe"], [-0.5, 0.6, 0.3, 0.4]),
                         ("mtrC", ELEM["Fe"], [-0.45, 0.55, 0.28, 0.38])])
        draw_cell(ax, cx=5.8, cy=0.7, w=1.7, h=0.95,
                  name="FeOB", mag="Sideroxydans",
                  genes=[("foxC", ELEM["Fe"], [0.3, 0.5, -0.1, 0.45])])

    ax.text(0.2, 4.6, "Fe cycle", fontsize=11, fontweight="bold",
            color=ELEM["Fe"])


# ======================================================================
# Mockup 06: 单组详细版（无 FeOOH 预设）
# ======================================================================
fig = plt.figure(figsize=(16, 11))
fig.suptitle(
    "Mockup 06 — Data-driven 2×2 element grid (Group CK, no forced central structure)",
    fontsize=13, fontweight="bold", y=0.98,
)

for i, (draw_fn, name) in enumerate([
    (draw_as_quadrant, "As"),
    (draw_n_quadrant,  "N"),
    (draw_s_quadrant,  "S"),
    (draw_fe_quadrant, "Fe"),
]):
    row, col = divmod(i, 2)
    ax = fig.add_axes([0.05 + col * 0.48, 0.52 - row * 0.45, 0.44, 0.40])
    ax.set_xlim(0, 8); ax.set_ylim(0, 5)
    ax.set_aspect("auto"); ax.axis("off")
    ax.set_facecolor("#FAFAFA")
    for spine in ax.spines.values():
        spine.set_visible(False)
    draw_fn(ax, detail=True)
    # 象限边框
    rect = Rectangle((0, 0), 8, 5, fill=False, edgecolor="#CCC",
                     linewidth=1.0, transform=ax.transData)
    ax.add_patch(rect)

# 元素间耦合线（KB 里的 coupling 字段）
fig.add_artist(plt.Line2D([0.53, 0.53], [0.55, 0.88],
                          color="#999", linestyle="--", alpha=0.5,
                          transform=fig.transFigure))
fig.text(0.53, 0.925, "coupling", rotation=0, ha="center", fontsize=8,
         color="#888", style="italic")

# colorbar + 图例
cax = fig.add_axes([0.35, 0.03, 0.3, 0.015])
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-1, 1))
plt.colorbar(sm, cax=cax, orientation="horizontal")
cax.set_title("Spearman ρ (gene × env)", fontsize=8, pad=2)
cax.tick_params(labelsize=7)

fig.text(0.5, 0.005,
         "Layout: KB elements determine grid (2×2 here).  "
         "Central area is neutral — no forced mineral/cell assumption.  "
         "Couplings drawn from KB coupling field (optional dashed lines).",
         ha="center", fontsize=8, color="#555", fontstyle="italic")

plt.savefig(OUT / "06_single_group_neutral.png", dpi=140,
            bbox_inches="tight", facecolor="white")
plt.close(fig)
print("[OK] 06_single_group_neutral.png")


# ======================================================================
# Mockup 07: 3 组横向对比（简化版）
# ======================================================================
fig = plt.figure(figsize=(24, 8))
fig.suptitle(
    "Mockup 07 — 3 groups side-by-side (CK / A / B, simplified)",
    fontsize=14, fontweight="bold", y=0.98,
)

# 每组不同宽度模拟组间差异
group_w_sets = {
    "CK": {
        "As": {"ox": 1.0, "red": 1.0, "me": 0.8, "resp": 0.5},
        "N":  {"nar": 1.0, "nir": 1.0, "nor": 0.8, "nos": 0.6, "amo": 0.3, "nif": 0.4},
        "S":  {"dsr": 1.0, "sox": 1.0, "sqr": 0.8, "tsd": 0.5},
        "Fe": {"dirb": 1.0, "fox": 0.8, "uptake": 1.0},
    },
    "A": {
        "As": {"ox": 1.5, "red": 2.0, "me": 1.2, "resp": 0.8},
        "N":  {"nar": 1.5, "nir": 1.5, "nor": 1.2, "nos": 1.0, "amo": 0.5, "nif": 0.8},
        "S":  {"dsr": 1.5, "sox": 1.5, "sqr": 1.2, "tsd": 0.8},
        "Fe": {"dirb": 1.5, "fox": 1.2, "uptake": 1.5},
    },
    "B": {
        "As": {"ox": 2.5, "red": 3.0, "me": 2.0, "resp": 1.5},  # B 组 As 循环最活跃
        "N":  {"nar": 2.0, "nir": 2.5, "nor": 2.2, "nos": 1.8, "amo": 0.2, "nif": 1.2},
        "S":  {"dsr": 3.0, "sox": 2.0, "sqr": 1.8, "tsd": 1.2},  # B 组 SRB 强
        "Fe": {"dirb": 2.0, "fox": 1.8, "uptake": 2.5},
    },
}

for i, g in enumerate(["CK", "A", "B"]):
    # 组标题条
    left = 0.03 + i * 0.325
    fig.add_axes([left, 0.92, 0.32, 0.04]).axis("off")
    fig.text(left + 0.16, 0.94, f"Group {g}", ha="center", va="center",
             fontsize=13, fontweight="bold", color="white",
             bbox=dict(facecolor=GROUP_COLOR[g], pad=5,
                       boxstyle="round,pad=0.3", edgecolor="none"))

    gw = group_w_sets[g]
    for j, (fn, name) in enumerate([
        (draw_as_quadrant, "As"),
        (draw_n_quadrant,  "N"),
        (draw_s_quadrant,  "S"),
        (draw_fe_quadrant, "Fe"),
    ]):
        row, col = divmod(j, 2)
        ax = fig.add_axes([left + col * 0.16, 0.48 - row * 0.43,
                           0.155, 0.40])
        ax.set_xlim(0, 8); ax.set_ylim(0, 5)
        ax.axis("off")
        ax.set_facecolor("#FAFAFA")
        fn(ax, detail=False, group_widths=gw[name])  # 简化：不画细胞，只化学链
        rect = Rectangle((0, 0), 8, 5, fill=False, edgecolor="#DDD",
                         linewidth=0.8, transform=ax.transData)
        ax.add_patch(rect)

fig.text(0.5, 0.01,
         "Arrow thickness = group-specific total contribution  •  "
         "Drill-down: click a group to see cell-level detail (Mockup 06 style)",
         ha="center", fontsize=9, color="#555", fontstyle="italic")

plt.savefig(OUT / "07_three_groups_horizontal.png", dpi=130,
            bbox_inches="tight", facecolor="white")
plt.close(fig)
print("[OK] 07_three_groups_horizontal.png")
print("\n完成。")
