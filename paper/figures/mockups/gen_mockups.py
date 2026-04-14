"""生成循环图 3 层可视化的 mockup 示意图（非真实数据，仅示意）。

输出 3 张：
  01_macro_3groups.png — 宏观层 × 3 组对比（As 循环）
  02_micro_cell_N.png  — 微观层（亚细胞 N 循环，Fig 2A 风格）
  03_hierarchy_demo.png — 3 层整合概览（宏观 + 微观 + 组对比）
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle, Circle
from pathlib import Path

plt.rcParams["font.family"] = "DejaVu Sans"
OUT = Path(__file__).parent

ELEM_COLOR = {"As": "#E74C3C", "N": "#2C3E87", "S": "#1E8449", "Fe": "#B7600B"}
GROUP_COLOR = {"CK": "#1c9cbd", "A": "#e3943d", "B": "#92181e"}


# ======================================================================
# Mockup 1: 宏观层 — 3 组对比 As 循环
# ======================================================================
def draw_as_cycle(ax, group, top_mags, arrow_widths):
    """在 ax 上画 As 循环示意：As(V) ↔ As(III) → As-org。
    top_mags: dict pathway → top MAG
    arrow_widths: dict pathway → normalized 0.5–3.0
    """
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_facecolor("#FAFAFA")

    # 标题带组色
    ax.add_patch(Rectangle((0.3, 9.1), 9.4, 0.7,
                           facecolor=GROUP_COLOR[group], alpha=0.85,
                           transform=ax.transData, clip_on=False))
    ax.text(5, 9.45, f"Group {group} — Arsenic cycle",
            ha="center", va="center", fontsize=11, fontweight="bold",
            color="white")

    # 化学物节点
    nodes = {
        "As(V)":   (2, 6),
        "As(III)": (2, 2.5),
        "As-org":  (7.5, 2.5),
    }
    for name, (x, y) in nodes.items():
        circ = Circle((x, y), 0.8, facecolor="#FDEBD0",
                      edgecolor=ELEM_COLOR["As"], linewidth=2, zorder=3)
        ax.add_patch(circ)
        ax.text(x, y, name, ha="center", va="center",
                fontsize=10, fontweight="bold", zorder=4)

    # 箭头：化学转化
    arrows = [
        {"pw": "Arsenate red.",      "from": "As(V)",   "to": "As(III)",
         "xy_offset": (-0.5, 0)},
        {"pw": "Arsenite ox.",       "from": "As(III)", "to": "As(V)",
         "xy_offset": (0.5, 0)},
        {"pw": "As methylation",     "from": "As(III)", "to": "As-org",
         "xy_offset": (0, -0.5)},
    ]
    for arr in arrows:
        x0, y0 = nodes[arr["from"]]
        x1, y1 = nodes[arr["to"]]
        dx, dy = arr["xy_offset"]
        lw = arrow_widths.get(arr["pw"], 1.5)
        fa = FancyArrowPatch(
            (x0 + dx, y0 + dy * 0.5), (x1 + dx, y1 + dy * 0.5),
            arrowstyle="-|>", mutation_scale=18,
            linewidth=lw, color=ELEM_COLOR["As"],
            connectionstyle="arc3,rad=0.2" if dx != 0 else "arc3,rad=0",
            zorder=2, alpha=0.85,
        )
        ax.add_patch(fa)
        # 通路标签
        mx = (x0 + x1) / 2 + dx * 1.5
        my = (y0 + y1) / 2 + dy * 1.2
        ax.text(mx, my, arr["pw"], fontsize=7.5, fontweight="bold",
                color=ELEM_COLOR["As"], ha="center", zorder=5,
                bbox=dict(facecolor="white", edgecolor="none",
                          alpha=0.9, pad=1.5))
        # Top MAG 标注
        top = top_mags.get(arr["pw"])
        if top:
            ax.text(mx, my - 0.4, f"→ {top}", fontsize=7,
                    fontstyle="italic", color="#333",
                    ha="center", zorder=5,
                    bbox=dict(facecolor="white", edgecolor="#AAA",
                              alpha=0.9, pad=1.5, boxstyle="round,pad=0.2"))


fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), constrained_layout=True)
fig.suptitle(
    "Mockup 1 — Macro layer: per-group element cycle comparison (As)",
    fontsize=13, fontweight="bold", y=1.02,
)

# 模拟三组数据（不真实，示意）
data_per_group = {
    "CK": {
        "tops": {"Arsenate red.": "Gallionella",
                 "Arsenite ox.":  "Gallionella",
                 "As methylation": "Fen-1038"},
        "widths": {"Arsenate red.": 1.5, "Arsenite ox.": 2.5,
                   "As methylation": 1.0},
    },
    "A": {
        "tops": {"Arsenate red.": "Gallionella",
                 "Arsenite ox.":  "SPC001",
                 "As methylation": "Fen-1038"},
        "widths": {"Arsenate red.": 2.2, "Arsenite ox.": 2.0,
                   "As methylation": 1.8},
    },
    "B": {
        "tops": {"Arsenate red.": "Defluviilinea",  # 注意阈值敏感
                 "Arsenite ox.":  "Gallionella",
                 "As methylation": "SPC001"},
        "widths": {"Arsenate red.": 3.0, "Arsenite ox.": 1.2,
                   "As methylation": 2.5},
    },
}
for ax, g in zip(axes, ["CK", "A", "B"]):
    draw_as_cycle(ax, g, data_per_group[g]["tops"],
                  data_per_group[g]["widths"])

# Footer 注释
fig.text(0.5, -0.02,
         "Arrow thickness = group-specific total contribution  •  "
         "Labels = top-completeness MAG (descriptive, not causal)",
         ha="center", fontsize=9, color="#555", fontstyle="italic")

plt.savefig(OUT / "01_macro_3groups.png", dpi=140, bbox_inches="tight",
            facecolor="white")
plt.close(fig)
print("[OK]01_macro_3groups.png")


# ======================================================================
# Mockup 2: 微观层 — N 循环 Fig 2A 风格
# ======================================================================
fig, ax = plt.subplots(figsize=(13, 7))
ax.set_xlim(0, 13); ax.set_ylim(0, 7.5)
ax.set_aspect("equal"); ax.axis("off")
ax.set_facecolor("#FAFAFA")
ax.text(6.5, 7.2, "Mockup 2 — Micro layer: subcellular N-cycle diagram (Fig 2A style)",
        ha="center", fontsize=12, fontweight="bold")

# 细胞结构 — 横向分层
# 0-1.5: periplasm (上)
# 1.5-2: inner membrane (上)
# 2-5: cytoplasm
# 5-5.5: inner membrane (下)
# 5.5-7: periplasm (下)
def band(y0, y1, color, label, label_x=0.3):
    ax.add_patch(Rectangle((0.5, y0), 12, y1 - y0,
                           facecolor=color, alpha=0.35, edgecolor="none"))
    ax.text(label_x, (y0 + y1) / 2, label, fontsize=8.5,
            ha="left", va="center", color="#555", fontstyle="italic",
            fontweight="bold")

band(5.8, 6.8, "#D4E5EF", "Periplasm")
band(5.5, 5.8, "#999999", "Inner membrane")
band(2.0, 5.5, "#FFFFFF", "Cytoplasm", label_x=0.3)
band(1.7, 2.0, "#999999", "Inner membrane", label_x=0.3)
band(0.7, 1.7, "#D4E5EF", "Periplasm", label_x=0.3)

# 化学物（底物/产物）
chem_nodes = {
    "NO3⁻": (1.7, 6.3, "#D4E5EF"),
    "NO2⁻": (5.0, 6.3, "#D4E5EF"),
    "NO":   (8.0, 6.3, "#D4E5EF"),
    "N2O":  (10.5, 6.3, "#D4E5EF"),
    "N2":   (12.0, 6.3, "#D4E5EF"),
    "NO3⁻(c)": (1.7, 3.5, "#FFFFFF"),
    "NO2⁻(c)": (5.0, 3.5, "#FFFFFF"),
}
for name, (x, y, bg) in chem_nodes.items():
    ax.add_patch(FancyBboxPatch((x - 0.55, y - 0.25), 1.1, 0.5,
                                 boxstyle="round,pad=0.05",
                                 facecolor=bg, edgecolor="#333", linewidth=1.2))
    ax.text(x, y, name, ha="center", va="center", fontsize=8.5,
            fontweight="bold")

# 酶（KO）的亚细胞定位
#   periplasm-facing: NapA, NirK/S, NosZ
#   inner membrane cytoplasmic: NarG, NarH, NorB
enzymes = [
    # name, (x, y), location, color, correlation_vec (fake)
    ("NapA", (3.0, 6.0), "periplasm",   "#3498DB",  [ 0.6, 0.2, -0.1,  0.8]),
    ("NapB", (3.0, 5.65),"inner_mem",   "#5FAEDC", [ 0.5, 0.1, -0.2,  0.7]),
    ("NarG", (3.0, 3.0), "cytoplasm",   "#8E44AD", [-0.3, 0.4,  0.6,  0.5]),
    ("NarH", (3.6, 3.0), "cytoplasm",   "#A569BD", [-0.2, 0.3,  0.5,  0.5]),
    ("NirS", (6.5, 6.0), "periplasm",   "#2ECC71", [ 0.1, 0.8, -0.2,  0.6]),
    ("NirK", (6.5, 6.4), "periplasm",   "#58D68D", [ 0.2, 0.7, -0.1,  0.5]),
    ("NorB", (9.0, 5.65),"inner_mem",   "#F39C12", [ 0.0, 0.5,  0.1,  0.4]),
    ("NosZ", (11.0, 6.0),"periplasm",   "#E67E22", [ 0.3, 0.4,  0.2,  0.6]),
]

env_labels = ["pH", "Eh", "TOC", "As"]

def draw_enzyme(ax, name, x, y, color, corr):
    # 酶本体
    ax.add_patch(FancyBboxPatch((x - 0.38, y - 0.17), 0.76, 0.34,
                                 boxstyle="round,pad=0.03",
                                 facecolor=color, edgecolor="#222",
                                 linewidth=1.2, zorder=5))
    ax.text(x, y, name, ha="center", va="center", fontsize=7,
            fontweight="bold", color="white", zorder=6)
    # 相关性 mini 热图（4 格）
    hm_x = x - 0.38
    hm_y = y + 0.22
    for i, r in enumerate(corr):
        cx = hm_x + i * 0.2
        val_color = plt.cm.RdBu_r((r + 1) / 2)
        ax.add_patch(Rectangle((cx, hm_y), 0.19, 0.19,
                               facecolor=val_color, edgecolor="#333",
                               linewidth=0.5, zorder=6))
        ax.text(cx + 0.095, hm_y + 0.095, env_labels[i],
                ha="center", va="center", fontsize=4.5,
                color="white" if abs(r) > 0.3 else "#333", zorder=7)

for name, pos, loc, color, corr in enzymes:
    draw_enzyme(ax, name, pos[0], pos[1], color, corr)

# 化学流向箭头
flow_arrows = [
    ((1.7, 6.05), (2.6, 6.0), "#555"),    # NO3 → NapA
    ((3.4, 6.0), (4.4, 6.15), "#555"),    # NapA → NO2
    ((5.0, 6.05), (6.1, 6.2), "#555"),    # NO2 → NirS/K
    ((6.9, 6.2), (7.6, 6.3), "#555"),     # NirS/K → NO
    ((8.0, 6.05), (8.8, 5.85), "#555"),   # NO → NorB
    ((9.2, 5.7), (10.0, 6.1), "#555"),    # NorB → N2O
    ((10.5, 6.05), (10.8, 6.05), "#555"), # N2O → NosZ (short)
    ((11.2, 6.0), (11.7, 6.15), "#555"),  # NosZ → N2
    # cytoplasmic leg
    ((1.7, 3.8), (2.6, 3.0), "#777"),     # NO3(c) → NarG
    ((3.4, 3.0), (4.4, 3.25), "#777"),    # NarG → NO2(c)
    # transport across membrane
    ((1.7, 6.0), (1.7, 3.8), "#AAA"),     # NO3 periplasm → cytoplasm
    ((5.0, 3.6), (5.0, 6.0), "#AAA"),     # NO2(c) → periplasm
]
for (x0, y0), (x1, y1), color in flow_arrows:
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=13,
        linewidth=1.2, color=color, alpha=0.8, zorder=3,
    ))

# 图例：关联性色阶 + enzymes
import matplotlib.cm as cm
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-1, 1))
cax = fig.add_axes([0.35, 0.07, 0.3, 0.015])
plt.colorbar(sm, cax=cax, orientation="horizontal", label="Spearman ρ (enzyme vs env)")

# 反硝化 / DNRA 说明
ax.text(6.5, 0.3, "Solid = denitrification  ·  Dashed = DNRA  ·  "
                  "Enzyme color = MAG abundance × copy  ·  "
                  "Mini 4-grid = Spearman vs pH/Eh/TOC/As",
        ha="center", fontsize=8, color="#555", fontstyle="italic")

plt.savefig(OUT / "02_micro_cell_N.png", dpi=140, bbox_inches="tight",
            facecolor="white")
plt.close(fig)
print("[OK]02_micro_cell_N.png")


# ======================================================================
# Mockup 3: 3 层整合概览（缩略版）
# ======================================================================
fig = plt.figure(figsize=(16, 10))
fig.suptitle("Mockup 3 — 3-layer cycle diagram architecture (hierarchy overview)",
             fontsize=14, fontweight="bold", y=0.99)

# --- Layer 1: 宏观 × 3 组（顶部）-------------------------
for i, g in enumerate(["CK", "A", "B"]):
    ax = fig.add_axes([0.04 + i * 0.32, 0.62, 0.29, 0.32])
    draw_as_cycle(ax, g, data_per_group[g]["tops"],
                  data_per_group[g]["widths"])

fig.text(0.5, 0.96,
         "Layer 1 — Macro: element-microbe interactions per group",
         ha="center", fontsize=11, fontweight="bold", color="#333")

# --- Layer 2: 微观（中部）-------------------------------
ax_m = fig.add_axes([0.05, 0.21, 0.9, 0.35])
ax_m.set_xlim(0, 13); ax_m.set_ylim(0, 7.5); ax_m.set_aspect("auto")
ax_m.axis("off"); ax_m.set_facecolor("#FAFAFA")

for y0, y1, c, lbl in [(5.8, 6.8, "#D4E5EF", "Periplasm"),
                        (5.5, 5.8, "#999",    "IM"),
                        (2.0, 5.5, "#FFF",    "Cytoplasm"),
                        (1.7, 2.0, "#999",    "IM"),
                        (0.7, 1.7, "#D4E5EF", "Periplasm")]:
    ax_m.add_patch(Rectangle((0.5, y0), 12, y1 - y0,
                             facecolor=c, alpha=0.35))
    ax_m.text(0.3, (y0 + y1) / 2, lbl, fontsize=7,
              ha="left", va="center", color="#555", fontstyle="italic")

for name, pos, loc, color, corr in enzymes:
    draw_enzyme(ax_m, name, pos[0], pos[1], color, corr)

# chem nodes
for name, (x, y, bg) in chem_nodes.items():
    ax_m.add_patch(FancyBboxPatch((x - 0.55, y - 0.25), 1.1, 0.5,
                                   boxstyle="round,pad=0.05",
                                   facecolor=bg, edgecolor="#333",
                                   linewidth=1.0))
    ax_m.text(x, y, name, ha="center", va="center", fontsize=7.5,
              fontweight="bold")

for (x0, y0), (x1, y1), color in flow_arrows:
    ax_m.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11,
        linewidth=1.0, color=color, alpha=0.8,
    ))

fig.text(0.5, 0.585,
         "Layer 2 — Micro: subcellular enzyme activity inside key MAGs",
         ha="center", fontsize=11, fontweight="bold", color="#333")

# --- Layer 3 (连接示意图 + 底部 env 耦合条) -----------
ax_e = fig.add_axes([0.1, 0.04, 0.8, 0.12])
ax_e.axis("off")
fig.text(0.5, 0.18,
         "Layer 3 — Group comparison + env coupling "
         "(shown via arrow widths / colors + correlation mini-heatmaps)",
         ha="center", fontsize=10, fontweight="bold", color="#333")

env_corr_demo = [
    ("Arsenate reduction", "Total_As", 0.85, 0.002),
    ("Ammonia oxidation",  "Eh",       0.84, 0.003),
    ("As methylation",     "Total_As", 0.76, 0.010),
    ("N fixation",         "Eh",       0.76, 0.010),
]
for i, (pw, env, r, p) in enumerate(env_corr_demo):
    y = 0.7 - i * 0.18
    ax_e.add_patch(Rectangle((0.30, y - 0.05), abs(r) * 0.4,
                             0.10, facecolor="#C0392B" if r > 0 else "#2980B9",
                             alpha=0.75, transform=ax_e.transAxes))
    ax_e.text(0.01, y, pw, transform=ax_e.transAxes,
              fontsize=8, va="center")
    ax_e.text(0.24, y, env, transform=ax_e.transAxes,
              fontsize=8, va="center", fontweight="bold")
    ax_e.text(0.72, y, f"ρ={r:.2f}  p={p:.3f}",
              transform=ax_e.transAxes, fontsize=8, va="center",
              color="#555")

plt.savefig(OUT / "03_hierarchy_demo.png", dpi=140, bbox_inches="tight",
            facecolor="white")
plt.close(fig)
print("[OK]03_hierarchy_demo.png")

print("\n完成：3 张 mockup 保存到 paper/figures/mockups/")
