"""Mockup 08 & 09 — 化学物流经细胞（Fig 1-d 风格）+ 元素耦合线。

设计：
- 每个 MAG/功能菌群 = 一个细胞矩形
- 底物（化学物节点）在细胞外左侧 → 箭头穿入细胞
- 基因椭圆在细胞内
- 产物在细胞外右侧 ← 箭头穿出细胞
- 元素象限之间用虚线连接（KB coupling 字段驱动）

输出：
  08_through_cells_single.png — 单组详细版（4 元素 + 耦合线）
  09_through_cells_3groups.png — 3 组横向对比（简化）
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import (
    Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle
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


def mini_heatmap(ax, x, y, vals, size=0.07):
    for i, r in enumerate(vals):
        ax.add_patch(Rectangle(
            (x + i * size, y), size, size,
            facecolor=plt.cm.RdBu_r((r + 1) / 2),
            edgecolor="#333", linewidth=0.3, zorder=9,
        ))


def gene_oval(ax, x, y, name, color, corr=None, w=0.7, h=0.32, show_hm=True):
    ax.add_patch(Ellipse((x, y), w, h, facecolor=color,
                         edgecolor="#222", linewidth=0.9, zorder=7))
    ax.text(x, y, name, ha="center", va="center",
            fontsize=6.5, fontweight="bold", color="white", zorder=8)
    if show_hm and corr is not None:
        hm_size = 0.08
        mini_heatmap(ax, x - 2 * hm_size, y + h / 2 + 0.04, corr, hm_size)


def draw_cell_through(ax, cy, cell_x0, cell_x1, cell_h,
                       name, mag, genes, reactions, show_hm=True):
    """
    绘制"化学流经细胞"的单行：
    外部底物 → [穿膜箭头] → 细胞内基因 → [穿膜箭头] → 外部产物

    reactions: list of (substrate, product, arrow_y_offset)
      示例：[("As(V)", "As(III)", 0), ("As(III)", "As-org", 0.2)]
    """
    # 细胞膜
    ax.add_patch(FancyBboxPatch(
        (cell_x0, cy - cell_h / 2), cell_x1 - cell_x0, cell_h,
        boxstyle="round,pad=0.1",
        facecolor="#FDF6E3", edgecolor="#6B4F2A", linewidth=1.5,
        zorder=5,
    ))
    # 菌群/MAG 名称
    ax.text(cell_x0 + 0.15, cy + cell_h / 2 - 0.3, name,
            fontsize=8, fontweight="bold", color="#333", zorder=8)
    ax.text(cell_x1 - 0.15, cy + cell_h / 2 - 0.3, f"[{mag}]",
            fontsize=7, fontstyle="italic", color="#555",
            ha="right", zorder=8)

    # 基因横向均匀布局
    n = len(genes)
    if n:
        gxs = np.linspace(cell_x0 + 0.6, cell_x1 - 0.6, n)
        gy = cy - 0.2
        for gx, (g, c, corr) in zip(gxs, genes):
            gene_oval(ax, gx, gy, g, c, corr, show_hm=show_hm)

    # 底物箭头（左侧穿入）+ 产物箭头（右侧穿出）
    for sub, prod, yoff in reactions:
        y = cy + yoff
        # 底物节点 + 箭头
        ax.add_patch(Ellipse((cell_x0 - 1.6, y), 0.95, 0.32,
                              facecolor="white", edgecolor="#222",
                              linewidth=1.2, zorder=6))
        ax.text(cell_x0 - 1.6, y, sub, ha="center", va="center",
                fontsize=8, fontweight="bold", zorder=7)
        ax.add_patch(FancyArrowPatch(
            (cell_x0 - 1.1, y), (cell_x0, y),
            arrowstyle="-|>", mutation_scale=14,
            linewidth=1.2, color="#555", zorder=4,
        ))
        # 产物节点 + 箭头
        ax.add_patch(Ellipse((cell_x1 + 1.6, y), 0.95, 0.32,
                              facecolor="white", edgecolor="#222",
                              linewidth=1.2, zorder=6))
        ax.text(cell_x1 + 1.6, y, prod, ha="center", va="center",
                fontsize=8, fontweight="bold", zorder=7)
        ax.add_patch(FancyArrowPatch(
            (cell_x1, y), (cell_x1 + 1.1, y),
            arrowstyle="-|>", mutation_scale=14,
            linewidth=1.2, color="#555", zorder=4,
        ))


# ======================================================================
# 4 元素象限定义
# ======================================================================

def draw_as_quadrant(ax, show_hm=True):
    """As 象限：3 行 = 3 个功能菌群。"""
    draw_cell_through(
        ax, cy=7.5, cell_x0=4, cell_x1=7, cell_h=1.2,
        name="Arsenite oxidizer", mag="Gallionella",
        genes=[("aioA", ELEM["As"], [0.6, 0.2, -0.1, 0.8]),
               ("aioB", ELEM["As"], [0.55, 0.15, -0.15, 0.75])],
        reactions=[("As(III)", "As(V)", 0)],
        show_hm=show_hm,
    )
    draw_cell_through(
        ax, cy=5.0, cell_x0=4, cell_x1=7, cell_h=1.2,
        name="DARPs", mag="Desulfobacterota",
        genes=[("arrA", ELEM["As2"], [0.8, 0.1, -0.2, 0.9]),
               ("arrB", ELEM["As2"], [0.75, 0.08, -0.18, 0.85]),
               ("arsC", ELEM["As"], [0.5, 0.2, -0.1, 0.7])],
        reactions=[("As(V)", "As(III)", 0)],
        show_hm=show_hm,
    )
    draw_cell_through(
        ax, cy=2.5, cell_x0=4, cell_x1=7, cell_h=1.2,
        name="As-methylator", mag="Fen-1038",
        genes=[("arsM", ELEM["As2"], [0.4, 0.3, 0.1, 0.65]),
               ("arsH", ELEM["As2"], [0.35, 0.25, 0.05, 0.6])],
        reactions=[("As(III)", "MAs(III)", 0)],
        show_hm=show_hm,
    )

    # 标题
    ax.text(0.3, 8.6, "As cycle", fontsize=13, fontweight="bold",
            color=ELEM["As"])


def draw_n_quadrant(ax, show_hm=True):
    draw_cell_through(
        ax, cy=7.5, cell_x0=4, cell_x1=7.5, cell_h=1.2,
        name="Denitrifier (nar/nir)", mag="Gallionella",
        genes=[("narG", ELEM["N"], [-0.3, 0.4, 0.6, 0.5]),
               ("narH", ELEM["N"], [-0.25, 0.35, 0.55, 0.48]),
               ("nirK", ELEM["N2"], [0.1, 0.7, -0.1, 0.5]),
               ("nirS", ELEM["N2"], [0.15, 0.75, -0.15, 0.55])],
        reactions=[("NO₃⁻", "NO₂⁻", 0.25),
                   ("NO₂⁻", "NO", -0.25)],
        show_hm=show_hm,
    )
    draw_cell_through(
        ax, cy=5.0, cell_x0=4, cell_x1=7.5, cell_h=1.2,
        name="Denitrifier (nor/nos)", mag="Gallionella",
        genes=[("norB", ELEM["N"], [0.0, 0.5, 0.1, 0.4]),
               ("norC", ELEM["N"], [-0.05, 0.45, 0.08, 0.38]),
               ("nosZ", ELEM["N"], [0.3, 0.4, 0.2, 0.6])],
        reactions=[("NO", "N₂O", 0.25),
                   ("N₂O", "N₂", -0.25)],
        show_hm=show_hm,
    )
    draw_cell_through(
        ax, cy=2.5, cell_x0=4, cell_x1=7.5, cell_h=1.2,
        name="N-fixer", mag="JACRMN01",
        genes=[("nifD", ELEM["N2"], [0.5, 0.76, -0.3, 0.2]),
               ("nifK", ELEM["N2"], [0.4, 0.72, -0.25, 0.2])],
        reactions=[("N₂", "NH₄⁺", 0)],
        show_hm=show_hm,
    )
    ax.text(0.3, 8.6, "N cycle", fontsize=13, fontweight="bold",
            color=ELEM["N"])


def draw_s_quadrant(ax, show_hm=True):
    draw_cell_through(
        ax, cy=7.5, cell_x0=4, cell_x1=7, cell_h=1.2,
        name="SRB (dissim.)", mag="Sulfuricaulis",
        genes=[("dsrA", ELEM["S"], [0.3, -0.2, 0.1, 0.5]),
               ("dsrB", ELEM["S"], [0.28, -0.22, 0.08, 0.48]),
               ("aprA", ELEM["S2"], [0.2, -0.1, 0.2, 0.4])],
        reactions=[("SO₄²⁻", "H₂S", 0)],
        show_hm=show_hm,
    )
    draw_cell_through(
        ax, cy=5.0, cell_x0=4, cell_x1=7, cell_h=1.2,
        name="SOB (Sox)", mag="Thiobacillus",
        genes=[("soxA", ELEM["S"], [0.1, 0.3, -0.2, 0.6]),
               ("soxB", ELEM["S"], [0.08, 0.32, -0.22, 0.58]),
               ("soxY", ELEM["S2"], [0.12, 0.28, -0.2, 0.55])],
        reactions=[("H₂S", "SO₄²⁻", 0.25),
                   ("S⁰", "SO₄²⁻", -0.25)],
        show_hm=show_hm,
    )
    draw_cell_through(
        ax, cy=2.5, cell_x0=4, cell_x1=7, cell_h=1.2,
        name="S-assim.", mag="Gp6-AA40",
        genes=[("cysJ", ELEM["S2"], [0.2, 0.1, 0.3, 0.4]),
               ("cysI", ELEM["S2"], [0.18, 0.12, 0.28, 0.38])],
        reactions=[("SO₄²⁻", "S-cys", 0)],
        show_hm=show_hm,
    )
    ax.text(0.3, 8.6, "S cycle", fontsize=13, fontweight="bold",
            color=ELEM["S"])


def draw_fe_quadrant(ax, show_hm=True):
    draw_cell_through(
        ax, cy=7.5, cell_x0=4, cell_x1=7, cell_h=1.2,
        name="DIRB (Fe reducer)", mag="Geobacteraceae",
        genes=[("mtrA", ELEM["Fe"], [-0.5, 0.6, 0.3, 0.4]),
               ("mtrC", ELEM["Fe"], [-0.45, 0.55, 0.28, 0.38]),
               ("omcE", ELEM["Fe2"], [-0.3, 0.5, 0.2, 0.3])],
        reactions=[("Fe(III)", "Fe(II)", 0)],
        show_hm=show_hm,
    )
    draw_cell_through(
        ax, cy=5.0, cell_x0=4, cell_x1=7, cell_h=1.2,
        name="FeOB (Fe oxidizer)", mag="Sideroxydans",
        genes=[("foxC", ELEM["Fe"], [0.3, 0.5, -0.1, 0.45]),
               ("cyc2", ELEM["Fe2"], [0.28, 0.48, -0.08, 0.42])],
        reactions=[("Fe(II)", "Fe(III)", 0)],
        show_hm=show_hm,
    )
    draw_cell_through(
        ax, cy=2.5, cell_x0=4, cell_x1=7, cell_h=1.2,
        name="Fe uptake", mag="DASXPG01",
        genes=[("fbpA", ELEM["Fe2"], [0.1, 0.2, 0.0, 0.3]),
               ("tonB", ELEM["Fe"], [0.15, 0.25, 0.0, 0.3])],
        reactions=[("Fe-ext", "Fe-int", 0)],
        show_hm=show_hm,
    )
    ax.text(0.3, 8.6, "Fe cycle", fontsize=13, fontweight="bold",
            color=ELEM["Fe"])


# ======================================================================
# Mockup 08: 单组 2×2 + 元素耦合线
# ======================================================================
fig = plt.figure(figsize=(22, 14))
fig.suptitle(
    "Mockup 08 — Substrate flows THROUGH cells (Fig 1-d style) + element couplings",
    fontsize=14, fontweight="bold", y=0.99,
)

quad_positions = {
    "As": (0.04, 0.52, 0.46, 0.40),   # 左上
    "N":  (0.52, 0.52, 0.46, 0.40),   # 右上
    "S":  (0.04, 0.06, 0.46, 0.40),   # 左下
    "Fe": (0.52, 0.06, 0.46, 0.40),   # 右下
}
quad_draws = {"As": draw_as_quadrant, "N": draw_n_quadrant,
              "S": draw_s_quadrant, "Fe": draw_fe_quadrant}

axes_by_elem = {}
for e, pos in quad_positions.items():
    ax = fig.add_axes(pos)
    ax.set_xlim(0, 11); ax.set_ylim(0, 9)
    ax.set_aspect("auto"); ax.axis("off")
    ax.set_facecolor("#FAFAFA")
    # 象限边框（淡色以示分区）
    rect = Rectangle((0, 0), 11, 9, fill=False, edgecolor=ELEM[e],
                     linewidth=1.5, alpha=0.5)
    ax.add_patch(rect)
    quad_draws[e](ax, show_hm=True)
    axes_by_elem[e] = ax

# ── 元素耦合线（KB coupling 字段驱动）─────────────────
# 在 fig 坐标系画连线 + 标签（穿过象限间空白区）
coupling_lines = [
    ("As", "S",  "As-S coprecipitation\n(As₂S₃, AsS)", "#B03A2E"),
    ("As", "Fe", "As adsorption\non FeOOH",            "#B9770E"),
    ("Fe", "S",  "FeS / FeS₂\nformation",              "#7D6608"),
    ("S",  "N",  "S oxidation\ncoupled w/ NO₃",        "#1D8348"),
]

def center_of(e):
    x, y, w, h = quad_positions[e]
    return (x + w / 2, y + h / 2)

for e1, e2, label, color in coupling_lines:
    x1, y1 = center_of(e1)
    x2, y2 = center_of(e2)
    line = plt.Line2D([x1, x2], [y1, y2], linestyle="--", linewidth=1.8,
                      color=color, alpha=0.55, transform=fig.transFigure,
                      zorder=1)
    fig.lines.append(line)
    # 标签放在线中点
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    fig.text(mx, my, label, ha="center", va="center",
             fontsize=8, color=color, fontweight="bold",
             bbox=dict(facecolor="white", alpha=0.85, edgecolor=color,
                       boxstyle="round,pad=0.3", linewidth=1))

# Colorbar (Spearman ρ)
cax = fig.add_axes([0.82, 0.015, 0.15, 0.018])
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-1, 1))
plt.colorbar(sm, cax=cax, orientation="horizontal")
cax.set_title("Spearman ρ (gene × env)", fontsize=8, pad=2)
cax.tick_params(labelsize=7)

# Footer
fig.text(0.5, 0.005,
         "Each row = one MAG.  Substrate (left node) → [cross membrane] → "
         "genes catalyze → [cross membrane] → product (right node).  "
         "Dashed lines = element coupling from KB.",
         ha="center", fontsize=9, color="#555", fontstyle="italic")

plt.savefig(OUT / "08_through_cells_single.png", dpi=130,
            bbox_inches="tight", facecolor="white")
plt.close(fig)
print("[OK] 08_through_cells_single.png")


# ======================================================================
# Mockup 09: 3 组横向对比（简化，每组 2×2，无 mini heatmap）
# ======================================================================
fig = plt.figure(figsize=(32, 12))
fig.suptitle(
    "Mockup 09 — 3 groups horizontal comparison (substrate-through-cell style)",
    fontsize=15, fontweight="bold", y=0.99,
)

for gi, g in enumerate(["CK", "A", "B"]):
    # 组色条
    left = 0.02 + gi * 0.33
    fig.text(left + 0.16, 0.965, f"Group {g}",
             ha="center", va="center", fontsize=14, fontweight="bold",
             color="white",
             bbox=dict(facecolor=GROUP_COLOR[g], pad=6,
                       boxstyle="round,pad=0.4", edgecolor="none"))

    # 2x2 元素网格（每组内）
    for ei, e in enumerate(["As", "N", "S", "Fe"]):
        row, col = divmod(ei, 2)
        ax = fig.add_axes([
            left + col * 0.16,
            0.50 - row * 0.44,
            0.155, 0.40
        ])
        ax.set_xlim(0, 11); ax.set_ylim(0, 9)
        ax.axis("off")
        ax.set_facecolor("#FAFAFA")
        rect = Rectangle((0, 0), 11, 9, fill=False, edgecolor=ELEM[e],
                         linewidth=1.0, alpha=0.4)
        ax.add_patch(rect)
        quad_draws[e](ax, show_hm=False)  # 简化：不画 mini heatmap

fig.text(0.5, 0.003,
         "Horizontal group comparison.  Drill-down any group → Mockup 08 "
         "(full detail with gene × env correlation grids).",
         ha="center", fontsize=10, color="#555", fontstyle="italic")

plt.savefig(OUT / "09_through_cells_3groups.png", dpi=110,
            bbox_inches="tight", facecolor="white")
plt.close(fig)
print("[OK] 09_through_cells_3groups.png")
print("\n完成。")
