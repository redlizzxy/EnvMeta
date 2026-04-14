"""Mockup 4 — 生态系统骨架 + 嵌入式细胞代谢图（融合用户手绘 Fig1+Fig2 思路）。

Figure 2 骨架：中心 FeOOH + 化学物 + 功能菌群位置
Figure 1 填充：每个菌群位置替换为其 Top MAG 的单细胞代谢图
  - 细胞内基因按元素配色（As=红/橙，N=蓝，S=绿，Fe=棕）
  - 每个基因上方加 4 格 mini 相关性热图（vs pH/Eh/TOC/As）
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import (
    Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle, Circle, Polygon
)
import numpy as np
from pathlib import Path

OUT = Path(__file__).parent
plt.rcParams["font.family"] = "DejaVu Sans"

# 元素配色（一级：元素。二级 alpha 细分通路，可按需做到 shades）
ELEM = {
    "As":  "#E74C3C",   # 红
    "As2": "#C0392B",   # 深红（As 另一通路示意）
    "N":   "#3498DB",   # 蓝
    "N2":  "#2874A6",
    "S":   "#229954",   # 绿
    "S2":  "#1E8449",
    "Fe":  "#CA8A1E",   # 棕
    "Fe2": "#A77215",
}
ENV_LABELS = ["pH", "Eh", "TOC", "As"]


def mini_heatmap(ax, x, y, vals, size=0.1):
    """4 格相关性热图。"""
    for i, r in enumerate(vals):
        ax.add_patch(Rectangle(
            (x + i * size, y), size, size,
            facecolor=plt.cm.RdBu_r((r + 1) / 2),
            edgecolor="#333", linewidth=0.3, zorder=8,
        ))


def gene_oval(ax, x, y, name, color, corr, w=0.55, h=0.34):
    """基因椭圆（=催化酶）+ 上方 4 格 mini 热图。"""
    ax.add_patch(Ellipse((x, y), w, h, facecolor=color,
                         edgecolor="#222", linewidth=0.9, zorder=6))
    ax.text(x, y, name, ha="center", va="center",
            fontsize=6.5, fontweight="bold", color="white", zorder=7)
    # 上方 mini heatmap
    hm_size = 0.11
    hm_x = x - 2 * hm_size
    hm_y = y + h / 2 + 0.03
    mini_heatmap(ax, hm_x, hm_y, corr, size=hm_size)


def draw_cell(ax, cx, cy, w, h, name, mag_label, genes, ec="#6B4F2A"):
    """绘制一个功能菌群 = 细胞代谢图。

    genes: list of (gene_name, color, corr)
    """
    # 细胞外框（柔和黄底 = Fig1 风格）
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.08",
        facecolor="#FDF6E3", edgecolor=ec, linewidth=1.5, zorder=4,
    ))
    # 顶部名字 + MAG 物种
    ax.text(cx, cy + h / 2 - 0.22, name, ha="center", va="top",
            fontsize=9, fontweight="bold", color="#222", zorder=7)
    ax.text(cx, cy + h / 2 - 0.55, mag_label, ha="center", va="top",
            fontsize=7, fontstyle="italic", color="#555", zorder=7)
    # 基因横向布局
    n = len(genes)
    gy = cy - 0.15
    if n > 0:
        gx_step = (w - 0.6) / max(n, 1)
        for i, (g, c, corr) in enumerate(genes):
            gx = cx - (w - 0.6) / 2 + i * gx_step + gx_step / 2 \
                 if n > 1 else cx
            gene_oval(ax, gx, gy, g, c, corr)


# ======================================================================
# 主图
# ======================================================================
fig, ax = plt.subplots(figsize=(17, 11))
ax.set_xlim(0, 22)
ax.set_ylim(0, 14)
ax.set_aspect("equal")
ax.axis("off")
ax.set_facecolor("#FCFCFC")

fig.suptitle(
    "Mockup 4 — Ecosystem-scaffolded cell diagrams (Group CK, Illustrative)",
    fontsize=14, fontweight="bold", y=0.98,
)

# ── 中心 FeOOH 矿物团 ───────────────────────────────────
# 用一个不规则多边形模拟 FeOOH "云团" 效果
feooh_points = np.array([
    [10.5, 7.8], [11.2, 8.3], [12.0, 8.2], [12.7, 7.9], [13.1, 7.3],
    [13.0, 6.6], [12.5, 6.0], [11.8, 5.7], [11.0, 5.8], [10.3, 6.1],
    [9.9, 6.8], [10.1, 7.5],
])
ax.add_patch(Polygon(feooh_points, facecolor="#D68910", alpha=0.55,
                     edgecolor="#7D6608", linewidth=1.2, zorder=2))
ax.text(11.5, 7.0, "FeOOH", ha="center", va="center", fontsize=14,
        fontweight="bold", color="#4D3A06", zorder=3)

# FeOOH 表面吸附的化学物种（模拟手绘红色斑块）
for x, y, label in [(10.6, 7.6, "As(V)"), (12.4, 7.6, "As(III)"),
                    (12.6, 6.4, "Fe(III)"), (10.7, 6.2, "Fe(II)")]:
    ax.add_patch(Ellipse((x, y), 0.55, 0.25, facecolor="#C0392B",
                         alpha=0.65, edgecolor="none", zorder=3))
    ax.text(x, y, label, ha="center", va="center", fontsize=7,
            fontweight="bold", color="white", zorder=4)

# ── 外部化学物池 ────────────────────────────────────────
def chem_label(x, y, txt, elem_color):
    ax.text(x, y, txt, ha="center", va="center", fontsize=9,
            fontweight="bold", color="#222",
            bbox=dict(facecolor="white", edgecolor=elem_color,
                      linewidth=1.2, boxstyle="round,pad=0.2"),
            zorder=5)

chem_label(2.5, 13.2, "NO₃⁻", ELEM["N"])
chem_label(1.5, 11.5, "N₂/NO₂⁻", ELEM["N"])
chem_label(1.5, 10.0, "N₂/NH₄⁺", ELEM["N"])
chem_label(1.5, 3.8, "SO₄²⁻", ELEM["S"])
chem_label(9.0, 1.5, "H₂S", ELEM["S"])
chem_label(13.5, 1.5, "S₂O₃²⁻", ELEM["S"])
chem_label(20.5, 3.0, "As-S", ELEM["As"])
chem_label(20.5, 10.0, "Fe(II)", ELEM["Fe"])

# ── 功能菌群（细胞代谢图）─────────────────────────────
# 左上：NAOB
draw_cell(ax, cx=4.5, cy=11.8, w=4.5, h=2.0,
          name="NAOB  (As(III) oxidizer + denitrifier)",
          mag_label="Top MAG: Gallionella sp. [Mx_All_63]",
          genes=[("aioA", ELEM["As"],  [0.60, 0.20, -0.10, 0.80]),
                 ("aioB", ELEM["As"],  [0.55, 0.15, -0.15, 0.75]),
                 ("narG", ELEM["N"],   [-0.30, 0.40, 0.60, 0.50]),
                 ("narH", ELEM["N"],   [-0.25, 0.35, 0.55, 0.48])])

# 右上：DARPs
draw_cell(ax, cx=18.0, cy=11.8, w=4.5, h=2.0,
          name="DARPs  (As(V) respiring reducer)",
          mag_label="Top MAG: Desulfobacterota [Mx_All_158]",
          genes=[("arrA", ELEM["As2"], [0.80, 0.10, -0.20, 0.90]),
                 ("arrB", ELEM["As2"], [0.75, 0.08, -0.18, 0.85]),
                 ("arsC", ELEM["As"],  [0.50, 0.20, -0.10, 0.70]),
                 ("arsB", ELEM["As"],  [0.40, 0.10, 0.00, 0.60])])

# 中左：DIRB
draw_cell(ax, cx=4.5, cy=7.5, w=4.5, h=2.0,
          name="DIRB  (Fe(III) reducer)",
          mag_label="Top MAG: Geobacteraceae [Mx_All_48]",
          genes=[("mtrA", ELEM["Fe"],  [-0.50, 0.60, 0.30, 0.40]),
                 ("mtrC", ELEM["Fe"],  [-0.45, 0.55, 0.28, 0.38]),
                 ("omcE", ELEM["Fe2"], [-0.30, 0.50, 0.20, 0.30])])

# 中右：NIOB / 含 Fe + N 耦合
draw_cell(ax, cx=18.0, cy=7.5, w=4.5, h=2.0,
          name="NIOB  (N-coupled Fe(II) oxidizer)",
          mag_label="Top MAG: Sideroxydans [Mx_All_76]",
          genes=[("foxC", ELEM["Fe"],  [0.30, 0.50, -0.10, 0.45]),
                 ("nirK", ELEM["N"],   [0.10, 0.80, -0.20, 0.60]),
                 ("nirS", ELEM["N2"],  [0.15, 0.75, -0.15, 0.55])])

# 左下：SRB
draw_cell(ax, cx=4.5, cy=3.2, w=4.5, h=2.0,
          name="SRB  (Dissimilatory sulfate reducer)",
          mag_label="Top MAG: Sulfuricaulis [Mx_All_158]",
          genes=[("dsrA", ELEM["S"],   [0.30, -0.20, 0.10, 0.50]),
                 ("dsrB", ELEM["S"],   [0.28, -0.22, 0.08, 0.48]),
                 ("aprA", ELEM["S2"],  [0.20, -0.10, 0.20, 0.40]),
                 ("aprB", ELEM["S2"],  [0.18, -0.12, 0.18, 0.38])])

# 右下：SOB + SOAsR（合并，标注双功能）
draw_cell(ax, cx=18.0, cy=3.2, w=4.8, h=2.2,
          name="SOB / SOAsR  (S oxidizer ± As coupling)",
          mag_label="Top MAG: Thiobacillus [Mx_All_48]",
          genes=[("soxA", ELEM["S"],   [0.10, 0.30, -0.20, 0.60]),
                 ("soxB", ELEM["S"],   [0.08, 0.32, -0.22, 0.58]),
                 ("soxY", ELEM["S2"],  [0.12, 0.28, -0.20, 0.55]),
                 ("arsC", ELEM["As"],  [0.40, 0.25, -0.15, 0.70])])

# ── 化学流向箭头（生态耦合）──────────────────────────
def flow(ax, x0, y0, x1, y1, color="#555", lw=1.2, style="-|>", alpha=0.85,
         rad=0.0):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle=style, mutation_scale=15,
        linewidth=lw, color=color, alpha=alpha, zorder=2,
        connectionstyle=f"arc3,rad={rad}",
    ))


# NAOB: NO3- → NAOB → N2/NO2- (denitrification) + As(III) → As(V) (oxidation)
flow(ax, 2.5, 12.9, 3.0, 12.2, ELEM["N"])       # NO3 in
flow(ax, 3.0, 11.4, 2.0, 11.5, ELEM["N"])       # out N2/NO2
flow(ax, 6.7, 11.6, 10.2, 7.6, ELEM["As"], rad=-0.15)  # As(III) oxidation onto surface

# DARPs: As(V) → As(III)
flow(ax, 15.8, 11.4, 12.6, 7.8, ELEM["As2"], rad=-0.15)  # As(V) into cell
flow(ax, 15.8, 11.7, 11.7, 7.5, ELEM["As"], lw=1.5, rad=0.15)  # As(III) out

# DIRB: Fe(III) → Fe(II)
flow(ax, 6.7, 7.4, 10.0, 6.5, ELEM["Fe"])
flow(ax, 10.4, 6.1, 6.7, 7.1, ELEM["Fe2"], style="-|>")  # back line with Fe(II)
# Fe(II) reservoir
flow(ax, 6.7, 7.1, 20.5, 10.0, ELEM["Fe2"], rad=-0.25, alpha=0.4)

# NIOB: Fe(II) + NO3/NO2 coupling
flow(ax, 15.8, 7.3, 13.0, 7.0, ELEM["Fe"], rad=0.1)
flow(ax, 16.2, 8.5, 2.5, 13.2, ELEM["N"], rad=-0.4, alpha=0.3)  # NO3- source

# SRB: SO4 → H2S
flow(ax, 2.5, 4.0, 3.0, 3.5, ELEM["S"])
flow(ax, 6.7, 3.0, 9.0, 1.8, ELEM["S2"])   # H2S out

# SOB: H2S → SO4/S2O3 (via Sox)
flow(ax, 9.5, 1.8, 15.8, 2.8, ELEM["S"])   # H2S → SOB
flow(ax, 18.0, 2.1, 13.5, 1.8, ELEM["S2"]) # S2O3 → ...

# SOB-AsS precipitation
flow(ax, 20.2, 2.8, 20.5, 2.9, ELEM["As2"])  # short As-S out

# FeOOH -> AsS formation (abiotic)
flow(ax, 13.1, 6.5, 20.3, 3.2, "#999", rad=-0.3, alpha=0.6)

# ── 图例 ────────────────────────────────────────────
leg_ax = fig.add_axes([0.05, 0.01, 0.9, 0.04])
leg_ax.axis("off")
leg_items = [
    ("As pathways",   ELEM["As"]),
    ("N pathways",    ELEM["N"]),
    ("S pathways",    ELEM["S"]),
    ("Fe pathways",   ELEM["Fe"]),
]
for i, (lbl, c) in enumerate(leg_items):
    x = 0.1 + i * 0.15
    leg_ax.add_patch(Ellipse((x, 0.5), 0.03, 0.4,
                             facecolor=c, edgecolor="#222",
                             transform=leg_ax.transAxes))
    leg_ax.text(x + 0.02, 0.5, lbl, transform=leg_ax.transAxes,
                va="center", fontsize=9)

# 相关性 colorbar
cax = fig.add_axes([0.75, 0.015, 0.18, 0.02])
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-1, 1))
plt.colorbar(sm, cax=cax, orientation="horizontal")
cax.set_title("Spearman ρ (gene × env)", fontsize=8, pad=2)
cax.tick_params(labelsize=7)

# 说明
ax.text(11, 0.3,
        "Each cell = a functional group's top-completeness MAG.  "
        "Gene ovals = catalyzing enzymes (colored by element/pathway).  "
        "Mini 4-grid above each gene = Spearman ρ vs pH / Eh / TOC / As",
        ha="center", fontsize=8, color="#555", fontstyle="italic")

plt.savefig(OUT / "04_ecosystem_embedded_cells.png",
            dpi=140, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("[OK] 04_ecosystem_embedded_cells.png")


# ======================================================================
# Mockup 5 — 3 组纵向对比（缩略版）
# ======================================================================
fig = plt.figure(figsize=(20, 15))
fig.suptitle(
    "Mockup 5 — Group stratification (CK / A / B, simplified)",
    fontsize=14, fontweight="bold", y=0.99,
)

def draw_mini_ecosystem(ax, group, color_seed=0):
    """简化版生态系统图（去掉 mini heatmap，保留结构）。"""
    ax.set_xlim(0, 22); ax.set_ylim(0, 14)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_facecolor("#FCFCFC")

    # 组标签条
    ax.add_patch(Rectangle((0, 13.3), 22, 0.6,
                           facecolor={"CK": "#1c9cbd", "A": "#e3943d",
                                      "B": "#92181e"}[group],
                           alpha=0.9))
    ax.text(11, 13.6, f"Group {group}", ha="center", va="center",
            fontsize=12, fontweight="bold", color="white")

    # FeOOH
    pts = np.array([
        [10.5, 8.3], [11.2, 8.8], [12.0, 8.7], [12.7, 8.4], [13.1, 7.8],
        [13.0, 7.1], [12.5, 6.5], [11.8, 6.2], [11.0, 6.3], [10.3, 6.6],
        [9.9, 7.3], [10.1, 8.0],
    ])
    ax.add_patch(Polygon(pts, facecolor="#D68910", alpha=0.55,
                         edgecolor="#7D6608"))
    ax.text(11.5, 7.5, "FeOOH", ha="center", va="center",
            fontsize=12, fontweight="bold", color="#4D3A06")

    # 组间差异：箭头粗细 & MAG 标签按 group 变化
    widths = {"CK": 1.0, "A": 2.0, "B": 3.0}
    top_mag_variants = {
        "NAOB": {"CK": "Gallionella", "A": "Gallionella", "B": "Defluviilinea"},
        "DARPs": {"CK": "Desulfo.sp1", "A": "Desulfo.sp2", "B": "Desulfo.sp3"},
        "SRB":   {"CK": "Sulfuricaulis", "A": "Sulfuricaulis", "B": "Sulfuricaulis"},
        "SOB":   {"CK": "Thiobacillus", "A": "Thiobacillus", "B": "Thiobacillus"},
        "DIRB":  {"CK": "Geobacter", "A": "Geobacter", "B": "DYIB01"},
        "NIOB":  {"CK": "Sideroxydans", "A": "Sideroxydans", "B": "Sideroxydans"},
    }

    positions = {
        "NAOB":  (4.5, 11.5),
        "DARPs": (18.0, 11.5),
        "DIRB":  (4.5, 8.0),
        "NIOB":  (18.0, 8.0),
        "SRB":   (4.5, 3.8),
        "SOB":   (18.0, 3.8),
    }
    elem_map = {"NAOB": ELEM["As"], "DARPs": ELEM["As2"],
                "DIRB": ELEM["Fe"], "NIOB": ELEM["Fe2"],
                "SRB": ELEM["S"], "SOB": ELEM["S2"]}

    for fg, (x, y) in positions.items():
        ax.add_patch(FancyBboxPatch(
            (x - 1.8, y - 0.7), 3.6, 1.4,
            boxstyle="round,pad=0.05",
            facecolor="#FDF6E3", edgecolor=elem_map[fg], linewidth=1.3))
        ax.text(x, y + 0.3, fg, ha="center", va="center",
                fontsize=9, fontweight="bold")
        ax.text(x, y - 0.2, top_mag_variants[fg][group],
                ha="center", va="center",
                fontsize=7, fontstyle="italic", color="#444")
        # 内部 3-4 个小点示意基因
        for i in range(3):
            gx = x - 1.3 + i * 1.3
            gy = y - 0.5
            ax.add_patch(Circle((gx, gy), 0.12,
                                facecolor=elem_map[fg], edgecolor="#222",
                                linewidth=0.5))

    # 化学流向箭头（按 group 粗细变化）
    w = widths[group]
    for (x0, y0), (x1, y1), c in [
        ((6.3, 11.5), (10.2, 8.0), ELEM["As"]),
        ((15.7, 11.5), (12.5, 8.0), ELEM["As2"]),
        ((6.3, 8.0), (10.0, 7.3), ELEM["Fe"]),
        ((15.7, 8.0), (13.0, 7.5), ELEM["Fe2"]),
        ((6.3, 3.8), (9.0, 2.5), ELEM["S"]),
        ((15.7, 3.8), (13.0, 2.5), ELEM["S2"]),
    ]:
        ax.add_patch(FancyArrowPatch(
            (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=16,
            linewidth=w, color=c, alpha=0.8, zorder=2,
        ))

axes = [fig.add_subplot(3, 1, i + 1) for i in range(3)]
for ax, g in zip(axes, ["CK", "A", "B"]):
    draw_mini_ecosystem(ax, g)

fig.text(0.5, 0.01,
         "Arrow thickness = group-specific total contribution  •  "
         "MAG labels in each cell change per group when data differs  •  "
         "(Gene interiors shown as simple dots in overview; detail view in Mockup 4)",
         ha="center", fontsize=9, color="#555", fontstyle="italic")

plt.savefig(OUT / "05_groups_stacked_ecosystem.png",
            dpi=130, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("[OK] 05_groups_stacked_ecosystem.png")


print("\n完成：2 张新 mockup 保存到 paper/figures/mockups/")
