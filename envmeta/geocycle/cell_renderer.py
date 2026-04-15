"""级联细胞渲染 primitive（Mockup 10 抽取）。

本模块提供单个"MAG 合并细胞"的 matplotlib 绘制函数，供主 renderer 反复调用：
底物（细胞外）→ 穿膜 → 基因 → 中间产物 → … → 最终产物（细胞外）。

这是**纯视觉原语**，不关心整体布局。主 renderer 负责：
  - 在每个元素象限里决定画几个细胞 / 位置
  - 调用 ``draw_cascade_cell()`` 画单个细胞
  - 收集返回的 ``substrate_pos`` / ``product_pos`` 用于后续化学物耦合连线

设计要点（S2.5-6）：
  - 化学物（外部 + 细胞内中间产物）**只有文字**，不画椭圆 / 边框 —— 与用户
    手绘 Fig 1-d 的视觉一致
  - 酶/基因仍以彩色椭圆呈现（醒目区分）
  - 化学式通过 ``_pretty_formula`` 转为 matplotlib mathtext，正确显示下标上标
    （SO4-2 → SO₄²⁻；As(III)/Fe(II) 保持罗马数字；As2S3 → As₂S₃ 等）
  - 连续相同的中间产物（常见：转运链里所有基因产物都是 Fe_internal）
    会被**折叠**，只画酶不画重复中间产物
"""
from __future__ import annotations

import re
from typing import Any

import numpy as np
from matplotlib.axes import Axes
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle


CELL_FILL = "#FDF6E3"
CELL_EDGE = "#6B4F2A"
INTERMEDIATE_COLOR = "#7D5F1A"
EXTERNAL_COLOR = "#222"
ARROW_COLOR = "#555"
OUTSIDE_ARROW_COLOR = "#444"


# --- 化学式 → matplotlib mathtext ------------------------------------------

_FORMULA_OVERRIDES: dict[str, str] = {
    "SO4-2":         r"$\mathrm{SO_4^{2-}}$",
    "SO4^2-":        r"$\mathrm{SO_4^{2-}}$",
    "SO4":           r"$\mathrm{SO_4}$",
    "SO3-2":         r"$\mathrm{SO_3^{2-}}$",
    "S-2":           r"$\mathrm{S^{2-}}$",
    "S2O3-2":        r"$\mathrm{S_2O_3^{2-}}$",
    "NO3-":          r"$\mathrm{NO_3^{-}}$",
    "NO2-":          r"$\mathrm{NO_2^{-}}$",
    "N2O":           r"$\mathrm{N_2O}$",
    "N2":            r"$\mathrm{N_2}$",
    "NH3":           r"$\mathrm{NH_3}$",
    "NH4+":          r"$\mathrm{NH_4^{+}}$",
    "NH2OH":         r"$\mathrm{NH_2OH}$",
    "H2S":           r"$\mathrm{H_2S}$",
    "As(V)":         r"$\mathrm{As(V)}$",
    "As(III)":       r"$\mathrm{As(III)}$",
    "As(V)-adduct":  r"$\mathrm{As(V)}$-adduct",
    "As(III)-organic": r"$\mathrm{As(III)}$-organic",
    "As(III)_out":   r"$\mathrm{As(III)_{out}}$",
    "As-glutathione": r"As-glutathione",
    "As-glutathione_out": r"As-glutathione$_{\mathrm{out}}$",
    "MMA/DMA":       r"MMA/DMA",
    "As2S3":         r"$\mathrm{As_2S_3}$",
    "FeS":           r"$\mathrm{FeS}$",
    "Fe(II)":        r"$\mathrm{Fe(II)}$",
    "Fe(III)":       r"$\mathrm{Fe(III)}$",
    "Fe_internal":   r"$\mathrm{Fe_{int}}$",
    "Fe-As_surface": r"Fe–As$_{\mathrm{surf}}$",
    "APS":           r"APS",
    "cysteine":      "cysteine",
    "S-cys":         "S-cys",
    "S0":            r"$\mathrm{S^{0}}$",
}


_RE_CHARGE = re.compile(r"^(?P<base>[A-Za-z][A-Za-z0-9()/.-]*?)(?P<chg>[-+]\d?|\^\{?[-+0-9]+\}?)$")
_RE_SUBSCRIPT = re.compile(r"(?P<elem>[A-Z][a-z]?)(?P<sub>\d+)")


def _pretty_formula(raw: str | None) -> str:
    """把原始化学字符串转成 matplotlib mathtext（自动识别下标/上标）。"""
    if raw is None:
        return ""
    s = str(raw).strip()
    if not s:
        return ""
    if s in _FORMULA_OVERRIDES:
        return _FORMULA_OVERRIDES[s]
    # 尝试 charge 后缀：NO3-, SO4-2, NH4+
    m = _RE_CHARGE.match(s)
    if m:
        base = m.group("base")
        chg = m.group("chg").lstrip("^").strip("{}")
        # 归一化 "-2" → "2-"；"+" → "+"
        if len(chg) > 1 and chg[0] in "-+":
            chg = chg[1:] + chg[0]
        # 数字转下标
        base_mm = _RE_SUBSCRIPT.sub(lambda mm: mm.group("elem") + "_{" + mm.group("sub") + "}", base)
        return rf"$\mathrm{{{base_mm}^{{{chg}}}}}$"
    # 只有下标（H2O, N2）
    if _RE_SUBSCRIPT.search(s):
        mm = _RE_SUBSCRIPT.sub(lambda x: x.group("elem") + "_{" + x.group("sub") + "}", s)
        return rf"$\mathrm{{{mm}}}$"
    return s


def mini_heatmap(ax: Axes, x: float, y: float, vals: list[float],
                 size: float = 0.09, cmap_name: str = "RdBu_r") -> None:
    """绘制 n×1 的 mini 相关热图。vals ∈ [-1, 1]。"""
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap(cmap_name)
    for i, r in enumerate(vals):
        if r is None or (isinstance(r, float) and np.isnan(r)):
            color = "#EEEEEE"
        else:
            color = cmap(max(-1.0, min(1.0, r)) / 2.0 + 0.5)
        ax.add_patch(Rectangle(
            (x + i * size, y), size, size,
            facecolor=color, edgecolor="#333", linewidth=0.3, zorder=9,
        ))


def _chem_outside(ax: Axes, x: float, y: float, txt: str) -> tuple[float, float]:
    """外部化学物节点（纯文字，无椭圆）。返回中心坐标。"""
    ax.text(x, y, _pretty_formula(txt), ha="center", va="center",
            fontsize=10, fontweight="bold", color=EXTERNAL_COLOR, zorder=7)
    return (x, y)


def _intermediate(ax: Axes, x: float, y: float, txt: str) -> tuple[float, float]:
    """细胞内部中间产物（纯文字，无黄底）。"""
    ax.text(x, y, _pretty_formula(txt), ha="center", va="center",
            fontsize=8, fontweight="bold", color=INTERMEDIATE_COLOR, zorder=7)
    return (x, y)


def _gene(ax: Axes, x: float, y: float, name: str, color: str,
          corr: list[float] | None = None,
          w: float = 0.75, h: float = 0.36, show_hm: bool = True) -> None:
    """基因椭圆（酶）+ 可选上方 mini heatmap。"""
    ax.add_patch(Ellipse((x, y), w, h, facecolor=color,
                         edgecolor="#222", linewidth=0.9, zorder=7))
    ax.text(x, y, name, ha="center", va="center",
            fontsize=7, fontweight="bold", color="white", zorder=8)
    if show_hm and corr is not None and len(corr) > 0:
        hm_size = 0.09
        mini_heatmap(ax, x - len(corr) * hm_size / 2.0,
                     y + h / 2 + 0.06, corr, hm_size)


def draw_cascade_cell(
    ax: Axes,
    *,
    cy: float,
    cell_x0: float,
    cell_w: float,
    cell_h: float,
    title: str,
    mag_label: str,
    steps: list[dict[str, Any]],
    element_color: str = "#888888",
    show_heatmap: bool = True,
    show_outside_chems: bool = True,
) -> dict[str, Any]:
    """画一个"合并细胞"带基因级联 (supports Mockup 10)。

    参数
    ------
    cy           细胞纵向中心
    cell_x0      细胞左边 x
    cell_w       细胞宽
    cell_h       细胞高
    title        细胞标题（常用 "通路名" 或 "功能角色"）
    mag_label    MAG 展示名（右上角，斜体）
    steps        [{"gene": str, "substrate": str|None, "product": str|None,
                   "corr": list[float]|None, "color": str|None}, ...]
                 steps[0]["substrate"] 作为细胞外起始底物
                 steps[i]["product"] == steps[i+1]["substrate"]（调用端保证）
    element_color 元素主色，用作基因默认填充
    show_heatmap  是否画每基因上方的 mini 相关性热图
    show_outside_chems 是否画细胞外的起始底物 + 最终产物节点 + 穿膜箭头

    返回
    ------
    {
      "substrate_pos": (x, y) | None  外部底物节点
      "product_pos":   (x, y) | None  外部最终产物节点
      "cell_box":      (x0, x1, cy)
      "gene_positions": [(x, y), ...]
    }
    """
    cell_x1 = cell_x0 + cell_w

    ax.add_patch(FancyBboxPatch(
        (cell_x0, cy - cell_h / 2), cell_w, cell_h,
        boxstyle="round,pad=0.1",
        facecolor=CELL_FILL, edgecolor=CELL_EDGE, linewidth=1.5, zorder=4,
    ))

    pad = min(0.2, cell_h * 0.2)
    ax.text(cell_x0 + pad, cy + cell_h / 2 - pad * 1.2, title,
            fontsize=8.5, fontweight="bold", color="#333", zorder=8)
    if mag_label:
        ax.text(cell_x1 - pad, cy + cell_h / 2 - pad * 1.2, f"[{mag_label}]",
                fontsize=7, fontstyle="italic", color="#555",
                ha="right", zorder=8)

    if not steps:
        ax.text((cell_x0 + cell_x1) / 2, cy - 0.1,
                "No active enzyme",
                ha="center", va="center",
                fontsize=9, color="#999", zorder=7)
        return {
            "substrate_pos": None,
            "product_pos": None,
            "cell_box": (cell_x0, cell_x1, cy),
            "gene_positions": [],
        }

    n_steps = len(steps)
    # 折叠连续相同的中间产物：若 steps[i].product == steps[i+1].substrate（常态）
    # 并且该产物已作为上一 intermediate 呈现，则不再重复画。
    # 检查最终产物链的 distinct 数（不含最后一步的最终产物）。
    inner_products = [str(s.get("product") or "") for s in steps[:-1]]
    all_same_intermediate = (
        len(inner_products) >= 1 and len(set(inner_products)) == 1
        and inner_products[0] != ""
    )

    inner_pad = 0.5
    inner_x0 = cell_x0 + inner_pad
    inner_x1 = cell_x1 - inner_pad
    gy = cy - (cell_h * 0.15)

    gene_positions: list[tuple[float, float]] = []
    if all_same_intermediate and n_steps >= 2:
        # 酶紧挨排成一行；中间产物仅显示一次（在酶群中央上方稍小字）
        xs_genes = list(np.linspace(inner_x0, inner_x1, n_steps))
        for x, s in zip(xs_genes, steps):
            _gene(ax, x, gy, name=s.get("gene", "?"),
                  color=s.get("color") or element_color,
                  corr=s.get("corr"), show_hm=show_heatmap)
            gene_positions.append((x, gy))
        # 酶群上方提示中间产物（小字一次）
        mid_x = (inner_x0 + inner_x1) / 2
        ax.text(mid_x, gy + cell_h * 0.28,
                "→ " + _pretty_formula(inner_products[0]) + " →",
                ha="center", va="center",
                fontsize=7, fontstyle="italic",
                color=INTERMEDIATE_COLOR, zorder=7)
        # 酶之间连接箭头（短）
        for i in range(len(xs_genes) - 1):
            ax.add_patch(FancyArrowPatch(
                (xs_genes[i] + 0.38, gy), (xs_genes[i + 1] - 0.38, gy),
                arrowstyle="-|>", mutation_scale=8,
                linewidth=0.8, color=ARROW_COLOR, zorder=3,
            ))
    else:
        slots = 2 * n_steps - 1
        if slots == 1:
            xs = [(inner_x0 + inner_x1) / 2]
        else:
            xs = list(np.linspace(inner_x0, inner_x1, slots))
        for idx, x in enumerate(xs):
            if idx % 2 == 0:
                si = idx // 2
                s = steps[si]
                _gene(ax, x, gy, name=s.get("gene", "?"),
                      color=s.get("color") or element_color,
                      corr=s.get("corr"), show_hm=show_heatmap)
                gene_positions.append((x, gy))
            else:
                si = (idx - 1) // 2
                inter = steps[si].get("product") or "?"
                _intermediate(ax, x, gy, str(inter))
        if slots > 1:
            for i in range(slots - 1):
                ax.add_patch(FancyArrowPatch(
                    (xs[i] + 0.38, gy), (xs[i + 1] - 0.38, gy),
                    arrowstyle="-|>", mutation_scale=10,
                    linewidth=0.9, color=ARROW_COLOR, zorder=3,
                ))

    substrate_pos = None
    product_pos = None
    if show_outside_chems:
        sub_txt = steps[0].get("substrate")
        if sub_txt:
            sub_x = cell_x0 - 1.1
            substrate_pos = _chem_outside(ax, sub_x, cy, str(sub_txt))
            ax.add_patch(FancyArrowPatch(
                (sub_x + 0.5, cy), (cell_x0 - 0.05, cy),
                arrowstyle="-|>", mutation_scale=14,
                linewidth=1.4, color=OUTSIDE_ARROW_COLOR, zorder=4,
            ))
        final = steps[-1].get("product")
        if final:
            prod_x = cell_x1 + 1.1
            product_pos = _chem_outside(ax, prod_x, cy, str(final))
            ax.add_patch(FancyArrowPatch(
                (cell_x1 + 0.05, cy), (prod_x - 0.5, cy),
                arrowstyle="-|>", mutation_scale=14,
                linewidth=1.4, color=OUTSIDE_ARROW_COLOR, zorder=4,
            ))

    return {
        "substrate_pos": substrate_pos,
        "product_pos": product_pos,
        "cell_box": (cell_x0, cell_x1, cy),
        "gene_positions": gene_positions,
    }


def genes_to_steps(
    gene_list: list[dict[str, Any]],
    *,
    default_color: str = "#888",
    env_corr: dict[str, list[float]] | None = None,
) -> list[dict[str, Any]]:
    """把 MAGContribution.genes → draw_cascade_cell 所需的 steps 格式。

    gene_list: [{"ko", "name", "substrate", "product"}, ...]
    env_corr:  可选 {ko: [rho_pH, rho_Eh, rho_TOC, rho_As]}，每项提供 mini 热图
    """
    out: list[dict[str, Any]] = []
    for g in gene_list:
        step = {
            "gene": g.get("name") or g.get("ko") or "?",
            "substrate": g.get("substrate"),
            "product": g.get("product"),
            "color": default_color,
            "corr": env_corr.get(g["ko"]) if env_corr else None,
        }
        out.append(step)
    return out
