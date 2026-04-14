"""Mockup 10 — 细胞内基因级联 + 化学物-化学物耦合线。

关键改进：
1. 每个 MAG = 1 个大细胞，内部多个基因**按代谢顺序排列**，基因之间标注**中间产物**
2. 元素耦合线连**具体化学物节点**（如 As(III) ─ H₂S → As-S），不是笼统元素标签
3. 入出细胞：首个底物在细胞左外、最终产物在细胞右外
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


def mini_heatmap(ax, x, y, vals, size=0.09):
    for i, r in enumerate(vals):
        ax.add_patch(Rectangle(
            (x + i * size, y), size, size,
            facecolor=plt.cm.RdBu_r((r + 1) / 2),
            edgecolor="#333", linewidth=0.3, zorder=9,
        ))


def chem_outside(ax, x, y, txt, size_w=1.0, size_h=0.38):
    """外部化学物节点（白底椭圆）。"""
    ax.add_patch(Ellipse((x, y), size_w, size_h,
                         facecolor="white", edgecolor="#222",
                         linewidth=1.3, zorder=6))
    ax.text(x, y, txt, ha="center", va="center",
            fontsize=9, fontweight="bold", zorder=7)
    return (x, y)


def intermediate(ax, x, y, txt):
    """细胞内部中间产物（浅黄小标签）。"""
    ax.text(x, y, txt, ha="center", va="center",
            fontsize=7.5, fontweight="bold", color="#7D5F1A",
            bbox=dict(facecolor="#FFF9C4", edgecolor="#7D5F1A",
                      linewidth=0.6, boxstyle="round,pad=0.18"),
            zorder=7)
    return (x, y)


def gene(ax, x, y, name, color, corr=None, w=0.75, h=0.36, show_hm=True):
    """基因椭圆（酶）+ 上方 mini heatmap。"""
    ax.add_patch(Ellipse((x, y), w, h, facecolor=color,
                         edgecolor="#222", linewidth=0.9, zorder=7))
    ax.text(x, y, name, ha="center", va="center",
            fontsize=7, fontweight="bold", color="white", zorder=8)
    if show_hm and corr is not None:
        hm_size = 0.09
        mini_heatmap(ax, x - 2 * hm_size, y + h / 2 + 0.06,
                     corr, hm_size)
    return (x, y)


def draw_merged_cell(ax, cy, x_left_chem, cell_x0, cell_w, cell_h,
                     name, mag, cascade, show_hm=True):
    """
    细胞内基因级联。
    cascade: [(gene_name, color, corr_vec, intermediate_after), ...]
      第 1 个 gene 的底物 = 第 1 行的 substrate（外部）
      第 k 个 gene 的产物 = 下一 gene 的底物 = intermediate（内部）
      最后 gene 的产物 = 外部产物

    返回 (input_node_pos, output_node_pos) 以便连耦合线。
    """
    cell_x1 = cell_x0 + cell_w

    # 细胞膜框
    ax.add_patch(FancyBboxPatch(
        (cell_x0, cy - cell_h / 2), cell_w, cell_h,
        boxstyle="round,pad=0.1",
        facecolor="#FDF6E3", edgecolor="#6B4F2A", linewidth=1.5,
        zorder=4,
    ))
    # 标题
    ax.text(cell_x0 + 0.2, cy + cell_h / 2 - 0.28, name,
            fontsize=8.5, fontweight="bold", color="#333", zorder=8)
    ax.text(cell_x1 - 0.2, cy + cell_h / 2 - 0.28, f"[{mag}]",
            fontsize=7, fontstyle="italic", color="#555",
            ha="right", zorder=8)

    # 基因 + 中间产物均匀布局
    n = len(cascade)
    inner_x0 = cell_x0 + 0.6
    inner_x1 = cell_x1 - 0.6
    # 位置：gene, intermediate, gene, intermediate, ..., gene
    # 若 n 个基因，有 n-1 个中间产物
    # 点数：2n - 1（但入/出端不在细胞内）
    positions = np.linspace(inner_x0, inner_x1, 2 * n - 1) if n > 1 else [
        (inner_x0 + inner_x1) / 2]
    gene_y = cy - 0.2

    gene_positions = []
    intermed_positions = []
    for i, pos in enumerate(positions):
        if i % 2 == 0:  # gene
            gi = i // 2
            gname, gcolor, corr, _inter = cascade[gi]
            gene(ax, pos, gene_y, gname, gcolor, corr, show_hm=show_hm)
            gene_positions.append(pos)
        else:  # intermediate between gene_{(i-1)/2} and gene_{(i+1)/2}
            gi = (i - 1) // 2
            _, _, _, inter_label = cascade[gi]
            intermediate(ax, pos, gene_y, inter_label)
            intermed_positions.append(pos)

    # 内部流向箭头（gene → inter → gene）
    # 短箭头连接相邻元素
    if len(positions) > 1:
        for i in range(len(positions) - 1):
            ax.add_patch(FancyArrowPatch(
                (positions[i] + 0.38, gene_y),
                (positions[i + 1] - 0.38, gene_y),
                arrowstyle="-|>", mutation_scale=10,
                linewidth=0.9, color="#555", zorder=3,
            ))

    # 外部底物（左）→ 穿膜 → 第 1 基因
    sub_txt = cascade[0][3]  # 这里复用 intermediate 字段作 substrate 名（见下 API 约定）
    # 修正：重构 cascade 数据，第 0 元素应含 substrate。改参数结构：
    # 我们在调用端传一个单独的 substrate + products 列表。
    # 这里使用 substrate 另传：
    pass  # 入/出箭头由调用端补充

    return {
        "cell_x0": cell_x0,
        "cell_x1": cell_x1,
        "cy": cy,
        "gene_positions": gene_positions,
        "intermed_positions": intermed_positions,
        "gene_y": gene_y,
    }


def full_cascade_row(ax, cy, x_left_chem, cell_x0, cell_w, cell_h,
                     name, mag, substrate, cascade_steps,
                     show_hm=True):
    """
    cascade_steps: [(gene_name, color, corr_vec), ...] 按顺序
                   最后一步的产物 = 细胞外最终产物（从 products_chain 尾元素）
    substrate: 外部底物（细胞左侧）
    """
    # cascade_steps 提供：gene 信息
    # 构造 intermediates: 除最后，每步后一个中间产物
    # 由调用者在 cascade_steps 后添加 final_product
    # 我们改一下参数：cascade = list of (gene, color, corr, produces)
    pass


# 重构 — 更简洁的 API
def draw_cascade_cell(ax, cy, cell_x0, cell_w, cell_h,
                      name, mag, substrate, steps,
                      show_hm=True):
    """
    steps: list of dicts {"gene": str, "color": hex, "corr": [4 floats], "product": str}
      steps[0]'s gene 把 substrate → steps[0]['product']
      steps[1]'s gene 把 steps[0]['product'] → steps[1]['product']
      ...
    最终：最后一步的 product 离开细胞

    返回 {"substrate_node": (x,y), "final_product_node": (x,y)}
    """
    cell_x1 = cell_x0 + cell_w

    # 细胞膜
    ax.add_patch(FancyBboxPatch(
        (cell_x0, cy - cell_h / 2), cell_w, cell_h,
        boxstyle="round,pad=0.1",
        facecolor="#FDF6E3", edgecolor="#6B4F2A", linewidth=1.5,
        zorder=4,
    ))
    # 标题
    ax.text(cell_x0 + 0.2, cy + cell_h / 2 - 0.28, name,
            fontsize=8.5, fontweight="bold", color="#333", zorder=8)
    ax.text(cell_x1 - 0.2, cy + cell_h / 2 - 0.28, f"[{mag}]",
            fontsize=7, fontstyle="italic", color="#555",
            ha="right", zorder=8)

    n_steps = len(steps)
    # 内部有 n_steps 个基因 + (n_steps-1) 个中间产物
    # 总共 2n-1 个元素均匀排列
    inner_x0 = cell_x0 + 0.5
    inner_x1 = cell_x1 - 0.5
    slots = 2 * n_steps - 1
    if slots == 1:
        xs = [(inner_x0 + inner_x1) / 2]
    else:
        xs = np.linspace(inner_x0, inner_x1, slots)
    gy = cy - 0.2

    for idx, x in enumerate(xs):
        if idx % 2 == 0:
            # gene
            si = idx // 2
            s = steps[si]
            gene(ax, x, gy, s["gene"], s["color"], s["corr"],
                 show_hm=show_hm)
        else:
            # intermediate = 前一个 gene 的 product（且是下一个 gene 的底物）
            si = (idx - 1) // 2
            inter = steps[si]["product"]
            intermediate(ax, x, gy, inter)

    # 内部 step 箭头
    if slots > 1:
        for i in range(slots - 1):
            ax.add_patch(FancyArrowPatch(
                (xs[i] + 0.38, gy), (xs[i + 1] - 0.38, gy),
                arrowstyle="-|>", mutation_scale=10,
                linewidth=0.9, color="#555", zorder=3,
            ))

    # 外部底物节点（左）
    sub_x = cell_x0 - 1.5
    sub_pos = chem_outside(ax, sub_x, cy, substrate)
    ax.add_patch(FancyArrowPatch(
        (sub_x + 0.5, cy), (cell_x0 - 0.05, cy),
        arrowstyle="-|>", mutation_scale=14,
        linewidth=1.4, color="#444", zorder=4,
    ))

    # 外部最终产物（右）
    final_product = steps[-1]["product"]
    prod_x = cell_x1 + 1.5
    prod_pos = chem_outside(ax, prod_x, cy, final_product)
    ax.add_patch(FancyArrowPatch(
        (cell_x1 + 0.05, cy), (prod_x - 0.5, cy),
        arrowstyle="-|>", mutation_scale=14,
        linewidth=1.4, color="#444", zorder=4,
    ))

    return {"substrate": sub_pos, "product": prod_pos,
            "cell": (cell_x0, cell_x1, cy)}


# ======================================================================
# 绘图主函数
# ======================================================================
fig = plt.figure(figsize=(26, 16))
fig.suptitle(
    "Mockup 10 — Merged cells + internal cascade + chemical-to-chemical couplings",
    fontsize=14, fontweight="bold", y=0.99,
)

quad_positions = {
    "As": (0.03, 0.52, 0.46, 0.42),
    "N":  (0.51, 0.52, 0.46, 0.42),
    "S":  (0.03, 0.06, 0.46, 0.42),
    "Fe": (0.51, 0.06, 0.46, 0.42),
}

# 存储关键化学物节点坐标（figure 坐标）→ 后续耦合用
chem_anchors = {}


def quadrant_ax(elem):
    pos = quad_positions[elem]
    ax = fig.add_axes(pos)
    ax.set_xlim(0, 18); ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_facecolor("#FAFAFA")
    rect = Rectangle((0, 0), 18, 10, fill=False, edgecolor=ELEM[elem],
                     linewidth=1.5, alpha=0.55)
    ax.add_patch(rect)
    return ax, pos


def ax_to_fig(ax, x, y):
    """把 ax 内数据坐标转 figure 坐标。"""
    return fig.transFigure.inverted().transform(
        ax.transData.transform((x, y)))


# ── As 象限 ────────────────────────────────────────────
ax, pos = quadrant_ax("As")
ax.text(0.3, 9.4, "As cycle", fontsize=14, fontweight="bold",
        color=ELEM["As"])

r1 = draw_cascade_cell(
    ax, cy=7.8, cell_x0=3.5, cell_w=3.5, cell_h=1.3,
    name="Arsenite oxidizer", mag="Gallionella",
    substrate="As(III)",
    steps=[{"gene": "aioA", "color": ELEM["As"],
            "corr": [0.6, 0.2, -0.1, 0.8], "product": "As(V)"}],
)
r2 = draw_cascade_cell(
    ax, cy=5.0, cell_x0=3.5, cell_w=5.5, cell_h=1.3,
    name="DARPs (respiratory + detox)", mag="Desulfobacterota",
    substrate="As(V)",
    steps=[{"gene": "arrA", "color": ELEM["As2"],
            "corr": [0.8, 0.1, -0.2, 0.9], "product": "As(III)_int"},
           {"gene": "arsB", "color": ELEM["As"],
            "corr": [0.4, 0.1, 0.0, 0.6], "product": "As(III)"}],
)
r3 = draw_cascade_cell(
    ax, cy=2.2, cell_x0=3.5, cell_w=5.5, cell_h=1.3,
    name="As-methylator", mag="Fen-1038",
    substrate="As(III)",
    steps=[{"gene": "arsM", "color": ELEM["As2"],
            "corr": [0.4, 0.3, 0.1, 0.65], "product": "MAs(III)"},
           {"gene": "arsH", "color": ELEM["As2"],
            "corr": [0.35, 0.25, 0.05, 0.6], "product": "DMAs(III)"}],
)

# 记录锚点：As(III) 的产物节点位置（来自 DARPs 或 Gallionella 输出）
chem_anchors["As(III)_out"] = ax_to_fig(
    ax, r2["product"][0], r2["product"][1])
chem_anchors["As(V)_out"] = ax_to_fig(
    ax, r1["product"][0], r1["product"][1])


# ── N 象限 ────────────────────────────────────────────
ax, pos = quadrant_ax("N")
ax.text(0.3, 9.4, "N cycle", fontsize=14, fontweight="bold",
        color=ELEM["N"])

# 一个大细胞：完整反硝化 + DNRA 由同一 MAG
rN1 = draw_cascade_cell(
    ax, cy=7.5, cell_x0=3.5, cell_w=10.5, cell_h=1.3,
    name="Complete denitrifier", mag="Gallionella",
    substrate="NO₃⁻",
    steps=[{"gene": "narG", "color": ELEM["N"],
            "corr": [-0.3, 0.4, 0.6, 0.5], "product": "NO₂⁻"},
           {"gene": "nirK", "color": ELEM["N2"],
            "corr": [0.1, 0.7, -0.1, 0.5], "product": "NO"},
           {"gene": "norB", "color": ELEM["N"],
            "corr": [0.0, 0.5, 0.1, 0.4], "product": "N₂O"},
           {"gene": "nosZ", "color": ELEM["N"],
            "corr": [0.3, 0.4, 0.2, 0.6], "product": "N₂"}],
)
rN2 = draw_cascade_cell(
    ax, cy=4.5, cell_x0=3.5, cell_w=4.5, cell_h=1.3,
    name="N-fixer", mag="JACRMN01",
    substrate="N₂",
    steps=[{"gene": "nifD", "color": ELEM["N2"],
            "corr": [0.5, 0.76, -0.3, 0.2], "product": "NH₄⁺"}],
)
rN3 = draw_cascade_cell(
    ax, cy=2.0, cell_x0=3.5, cell_w=4.5, cell_h=1.3,
    name="Ammonia oxidizer", mag="Nitrosotenuis",
    substrate="NH₄⁺",
    steps=[{"gene": "amoA", "color": ELEM["N"],
            "corr": [0.61, 0.84, -0.55, 0.22], "product": "NO₂⁻"}],
)

# ── S 象限 ────────────────────────────────────────────
ax, pos = quadrant_ax("S")
ax.text(0.3, 9.4, "S cycle", fontsize=14, fontweight="bold",
        color=ELEM["S"])

rS1 = draw_cascade_cell(
    ax, cy=7.5, cell_x0=3.5, cell_w=6.5, cell_h=1.3,
    name="SRB (dissim.)", mag="Sulfuricaulis",
    substrate="SO₄²⁻",
    steps=[{"gene": "aprA", "color": ELEM["S"],
            "corr": [0.3, -0.2, 0.1, 0.5], "product": "APS"},
           {"gene": "dsrA", "color": ELEM["S2"],
            "corr": [0.28, -0.22, 0.08, 0.48], "product": "SO₃²⁻"},
           {"gene": "dsrB", "color": ELEM["S2"],
            "corr": [0.25, -0.2, 0.08, 0.45], "product": "H₂S"}],
)
rS2 = draw_cascade_cell(
    ax, cy=4.8, cell_x0=3.5, cell_w=6.5, cell_h=1.3,
    name="SOB (Sox complex)", mag="Thiobacillus",
    substrate="H₂S",
    steps=[{"gene": "soxY", "color": ELEM["S"],
            "corr": [0.12, 0.28, -0.2, 0.55], "product": "S⁰"},
           {"gene": "soxA", "color": ELEM["S"],
            "corr": [0.1, 0.3, -0.2, 0.6], "product": "S₂O₃²⁻"},
           {"gene": "soxB", "color": ELEM["S2"],
            "corr": [0.08, 0.32, -0.22, 0.58], "product": "SO₄²⁻"}],
)
rS3 = draw_cascade_cell(
    ax, cy=2.0, cell_x0=3.5, cell_w=6.5, cell_h=1.3,
    name="S-assimilator", mag="Gp6-AA40",
    substrate="SO₄²⁻",
    steps=[{"gene": "cysJ", "color": ELEM["S"],
            "corr": [0.2, 0.1, 0.3, 0.4], "product": "SO₃²⁻"},
           {"gene": "cysI", "color": ELEM["S2"],
            "corr": [0.18, 0.12, 0.28, 0.38], "product": "S-cys"}],
)

chem_anchors["H2S_out"] = ax_to_fig(
    ax, rS1["product"][0], rS1["product"][1])
chem_anchors["SO4_out"] = ax_to_fig(
    ax, rS2["product"][0], rS2["product"][1])


# ── Fe 象限 ───────────────────────────────────────────
ax, pos = quadrant_ax("Fe")
ax.text(0.3, 9.4, "Fe cycle", fontsize=14, fontweight="bold",
        color=ELEM["Fe"])

rF1 = draw_cascade_cell(
    ax, cy=7.5, cell_x0=3.5, cell_w=5.5, cell_h=1.3,
    name="DIRB (Fe reducer)", mag="Geobacteraceae",
    substrate="Fe(III)",
    steps=[{"gene": "mtrA", "color": ELEM["Fe"],
            "corr": [-0.5, 0.6, 0.3, 0.4], "product": "Fe(III)-int"},
           {"gene": "mtrC", "color": ELEM["Fe"],
            "corr": [-0.45, 0.55, 0.28, 0.38], "product": "Fe(II)"}],
)
rF2 = draw_cascade_cell(
    ax, cy=4.5, cell_x0=3.5, cell_w=5.5, cell_h=1.3,
    name="FeOB (Fe oxidizer)", mag="Sideroxydans",
    substrate="Fe(II)",
    steps=[{"gene": "cyc2", "color": ELEM["Fe2"],
            "corr": [0.28, 0.48, -0.08, 0.42], "product": "e⁻→O₂"},
           {"gene": "foxC", "color": ELEM["Fe"],
            "corr": [0.3, 0.5, -0.1, 0.45], "product": "Fe(III)"}],
)

chem_anchors["Fe(II)_out"] = ax_to_fig(
    ax, rF1["product"][0], rF1["product"][1])
chem_anchors["Fe(III)_out"] = ax_to_fig(
    ax, rF2["product"][0], rF2["product"][1])


# ======================================================================
# 元素耦合线：连具体化学物 + 生成产物标签
# ======================================================================
def coupling_curve(p1, p2, product_label, color):
    """从 p1 (化学物节点 fig 坐标) 画曲线到 p2，中间放 product 节点。"""
    mx = (p1[0] + p2[0]) / 2
    my = (p1[1] + p2[1]) / 2
    # 曲线（用 PathPatch 不好跨 axes，直接用两段 Line2D）
    for (a, b) in [(p1, (mx, my)), ((mx, my), p2)]:
        line = plt.Line2D([a[0], b[0]], [a[1], b[1]],
                          linestyle="--", linewidth=1.8,
                          color=color, alpha=0.75,
                          transform=fig.transFigure, zorder=1)
        fig.lines.append(line)
    # 产物节点
    fig.text(mx, my, product_label, ha="center", va="center",
             fontsize=8.5, fontweight="bold", color="white",
             bbox=dict(facecolor=color, alpha=0.9, pad=6,
                       boxstyle="round,pad=0.3", edgecolor="black",
                       linewidth=1), zorder=10)


# 耦合 1: As(III) + H₂S → As₂S₃ / AsS
coupling_curve(chem_anchors["As(III)_out"], chem_anchors["H2S_out"],
               "As₂S₃ /\nAsS ↓", "#8E44AD")

# 耦合 2: Fe(III) + H₂S → FeS / FeS₂ + S⁰
coupling_curve(chem_anchors["Fe(III)_out"], chem_anchors["H2S_out"],
               "FeS /\nFeS₂ ↓", "#7D6608")

# 耦合 3: H₂S + NO₃⁻ (via N cycle) → SO₄²⁻ + N₂（S-ox 驱动反硝化）
# 这需要 N 的 NO₃ 锚点。从 rN1 substrate 位置拿。
chem_anchors["NO3_in"] = ax_to_fig(
    fig.axes[1], rN1["substrate"][0], rN1["substrate"][1])  # axes[1] = N quadrant
coupling_curve(chem_anchors["H2S_out"], chem_anchors["NO3_in"],
               "SO₄²⁻ +\nN₂ ↑", "#117A65")

# 耦合 4: As(V) + FeOOH 吸附（作为 surface 标识，连 As(V) 到 Fe(III)）
coupling_curve(chem_anchors["As(V)_out"], chem_anchors["Fe(III)_out"],
               "As-Fe\nadsorption", "#B9770E")


# Colorbar
cax = fig.add_axes([0.83, 0.01, 0.15, 0.015])
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-1, 1))
plt.colorbar(sm, cax=cax, orientation="horizontal")
cax.set_title("Spearman ρ (gene × env)", fontsize=8, pad=2)
cax.tick_params(labelsize=7)

# Footer
fig.text(0.5, 0.005,
         "Yellow labels = intermediate products inside cell  •  "
         "Dashed lines connect specific chemical species, labeled with coupling product  •  "
         "Each cell = 1 MAG with full enzyme cascade",
         ha="center", fontsize=9, color="#555", fontstyle="italic")

plt.savefig(OUT / "10_cascade_with_coupling.png", dpi=130,
            bbox_inches="tight", facecolor="white")
plt.close(fig)
print("[OK] 10_cascade_with_coupling.png")
