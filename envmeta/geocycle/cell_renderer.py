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
    "SO4-2":         r"$\mathbf{SO_4^{2-}}$",
    "SO4^2-":        r"$\mathbf{SO_4^{2-}}$",
    "SO4":           r"$\mathbf{SO_4}$",
    "SO3-2":         r"$\mathbf{SO_3^{2-}}$",
    "S-2":           r"$\mathbf{S^{2-}}$",
    "S2O3-2":        r"$\mathbf{S_2O_3^{2-}}$",
    "NO3-":          r"$\mathbf{NO_3^{-}}$",
    "NO2-":          r"$\mathbf{NO_2^{-}}$",
    "N2O":           r"$\mathbf{N_2O}$",
    "N2":            r"$\mathbf{N_2}$",
    "NH3":           r"$\mathbf{NH_3}$",
    "NH4+":          r"$\mathbf{NH_4^{+}}$",
    "NH2OH":         r"$\mathbf{NH_2OH}$",
    "H2S":           r"$\mathbf{H_2S}$",
    "As(V)":         r"$\mathbf{As(V)}$",
    "As(III)":       r"$\mathbf{As(III)}$",
    "As(V)-adduct":  r"$\mathbf{As(V)}$-adduct",
    "As(III)-organic": r"$\mathbf{As(III)}$-organic",
    "As(III)_out":   r"$\mathbf{As(III)_{out}}$",
    "As-glutathione": r"As-glutathione",
    "As-glutathione_out": r"As-glutathione$_{\mathbf{out}}$",
    "MMA/DMA":       r"MMA/DMA",
    "As2S3":         r"$\mathbf{As_2S_3}$",
    "FeS":           r"$\mathbf{FeS}$",
    "Fe(II)":        r"$\mathbf{Fe(II)}$",
    "Fe(III)":       r"$\mathbf{Fe(III)}$",
    "Fe_internal":   r"$\mathbf{Fe_{int}}$",
    "Fe-As_surface": r"Fe–As$_{\mathbf{surf}}$",
    "APS":           r"APS",
    "cysteine":      "cysteine",
    "S-cys":         "S-cys",
    "S0":            r"$\mathbf{S^{0}}$",
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
        return rf"$\mathbf{{{base_mm}^{{{chg}}}}}$"
    # 只有下标（H2O, N2）
    if _RE_SUBSCRIPT.search(s):
        mm = _RE_SUBSCRIPT.sub(lambda x: x.group("elem") + "_{" + x.group("sub") + "}", s)
        return rf"$\mathbf{{{mm}}}$"
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
    ax.text(x, y, _pretty_formula(txt), ha="center", va="center_baseline",
            fontsize=10, fontweight="bold", color=EXTERNAL_COLOR, zorder=7)
    return (x, y)


def _intermediate(ax: Axes, x: float, y: float, txt: str) -> tuple[float, float]:
    """细胞内部中间产物（纯文字，无黄底）。"""
    ax.text(x, y, _pretty_formula(txt), ha="center", va="center_baseline",
            fontsize=9, fontweight="bold", color=INTERMEDIATE_COLOR, zorder=7)
    return (x, y)


def _gene_size(name: str) -> tuple[float, float, int]:
    """按基因名长度算椭圆宽 / 高 / 字号。"""
    n = len(name or "")
    # 椭圆宽度动态：短名 0.75，长名最多拓到 1.5 数据单位
    w = max(0.75, min(1.55, 0.14 * n + 0.05))
    h = 0.36
    if n <= 6:
        fs = 7
    elif n <= 10:
        fs = 6
    else:
        fs = 5
    return w, h, fs


def _gene(ax: Axes, x: float, y: float, name: str, color: str,
          corr: list[float] | None = None,
          w: float | None = None, h: float | None = None,
          show_hm: bool = True) -> None:
    """基因椭圆（酶）+ 可选上方 mini heatmap。椭圆宽度按名长自适应。"""
    auto_w, auto_h, fs = _gene_size(name)
    w = w if w is not None else auto_w
    h = h if h is not None else auto_h
    ax.add_patch(Ellipse((x, y), w, h, facecolor=color,
                         edgecolor="#222", linewidth=0.9, zorder=7))
    ax.text(x, y, name, ha="center", va="center_baseline",
            fontsize=fs, fontweight="bold", color="white", zorder=8)
    if show_hm and corr is not None and len(corr) > 0:
        hm_size = 0.09
        mini_heatmap(ax, x - len(corr) * hm_size / 2.0,
                     y + h / 2 + 0.06, corr, hm_size)


def _bundle_label(steps: list[dict[str, Any]]) -> str:
    """对多基因 bundle 生成短标签，不截断。

    启发式：
    - 所有名字共享同一 3+ 字符前缀 → 前缀保留 + 各自后缀拼 '/'
      narG/narH/narI → narG/H/I
    - 否则全拼 'aprA/aprB/dsrA/dsrB'（无截断，靠 _draw_bundle_ellipse 自适应宽度）
    """
    names = [s.get("gene") or s.get("ko") or "?" for s in steps]
    if len(names) == 1:
        return names[0]
    pref = _common_prefix(names)
    if len(pref) >= 3 and all(len(n) > len(pref) for n in names):
        tails = [n[len(pref):] for n in names]
        return pref + "/".join(tails)
    return "/".join(names)


def _common_prefix(names: list[str]) -> str:
    if not names:
        return ""
    pref = names[0]
    for n in names[1:]:
        while not n.startswith(pref):
            pref = pref[:-1]
            if not pref:
                return ""
    return pref


def _draw_bundle_ellipse(
    ax: Axes, x: float, y: float, label: str, color: str,
    n_genes: int, w: float | None = None, h: float | None = None,
) -> None:
    """合并 bundle 椭圆：代替 N 个小椭圆的一个大椭圆。不截断文字，宽度自适应。

    超过一行负担时自动 wrap 成 2 行（椭圆高度加大）。
    """
    txt_len = len(label)
    # 超长标签自动换行（>18 字符 或 > 3 个 / 分隔符时换 2 行）
    if txt_len > 18 or label.count("/") > 3:
        display_label = _wrap_bundle_label(label)
        lines = display_label.count("\n") + 1
    else:
        display_label = label
        lines = 1
    # 椭圆宽度按每行实际显示长度
    max_line_len = max((len(s) for s in display_label.split("\n")), default=txt_len)
    auto_w = max(1.3, min(3.6, 0.13 * max_line_len + 0.25))
    auto_h = 0.46 if lines == 1 else 0.72
    w = w if w is not None else auto_w
    h = h if h is not None else auto_h
    fs_table = {1: (7, 6, 5), 2: (6, 5, 5)}
    fs_opts = fs_table.get(lines, (5, 5, 5))
    fs = fs_opts[0] if max_line_len <= 10 else fs_opts[1] if max_line_len <= 16 else fs_opts[2]
    ax.add_patch(Ellipse((x, y), w, h, facecolor=color,
                         edgecolor="#222", linewidth=1.1, zorder=7))
    ax.text(x, y, display_label, ha="center", va="center_baseline",
            fontsize=fs, fontweight="bold", color="white", zorder=8,
            linespacing=1.1)


def _wrap_bundle_label(label: str) -> str:
    """把长 bundle 标签沿 '/' 分隔平分成 2 行。"""
    parts = label.split("/")
    if len(parts) < 2:
        return label
    mid = len(parts) // 2
    return "/".join(parts[:mid]) + "/\n" + "/".join(parts[mid:])


def _draw_chain_row(
    ax: Axes,
    chain: list[list[int]],
    steps: list[dict[str, Any]],
    *,
    inner_x0: float, inner_x1: float, row_y: float,
    cell_h: float, element_color: str,
    show_heatmap: bool = False,
    external_chems: bool = False,
    cell_x0: float | None = None,
    cell_x1: float | None = None,
    gene_positions: list[tuple[float, float]] | None = None,
) -> tuple[tuple[float, float] | None, tuple[float, float] | None]:
    """画一条"链"在指定 row_y 上。

    chain: list of segments，段间化学物自然衔接
    inline_chems=True 时：链的起始 substrate 画在 inline_left_x 处（细胞内），
                          最终 product 画在 inline_right_x 处。
    """
    if gene_positions is None:
        gene_positions = []

    # 每段 1 slot；段间有 intermediate slot
    n_seg = len(chain)
    slots = 2 * n_seg - 1 if n_seg > 0 else 0
    if slots == 0:
        return (None, None)
    if slots == 1:
        xs = [(inner_x0 + inner_x1) / 2]
    else:
        xs = list(np.linspace(inner_x0, inner_x1, slots))

    for idx, x in enumerate(xs):
        if idx % 2 == 0:
            seg = chain[idx // 2]
            if len(seg) == 1:
                s = steps[seg[0]]
                _gene(ax, x, row_y, name=s.get("gene", "?"),
                      color=s.get("color") or element_color,
                      corr=s.get("corr"), show_hm=show_heatmap)
                gene_positions.append((x, row_y))
            else:
                seg_steps = [steps[i] for i in seg]
                label = _bundle_label(seg_steps)
                bundle_color = seg_steps[0].get("color") or element_color
                _draw_bundle_ellipse(ax, x, row_y, label, bundle_color,
                                      n_genes=len(seg))
                gene_positions.append((x, row_y))
                complex_id = seg_steps[0].get("complex") or ""
                sub = (f"({complex_id} complex, {len(seg)} subunits)"
                       if complex_id
                       else f"({len(seg)} isozymes)")
                ax.text(x, row_y - cell_h * 0.30, sub,
                        ha="center", va="center",
                        fontsize=6.5, fontstyle="italic",
                        color="#666", zorder=7)
        else:
            # 段间中间产物（链内一定衔接）
            left_seg = chain[(idx - 1) // 2]
            left_prod = str(steps[left_seg[-1]].get("product") or "?")
            _intermediate(ax, x, row_y, left_prod)

    if slots > 1:
        for i in range(slots - 1):
            ax.add_patch(FancyArrowPatch(
                (xs[i] + 0.38, row_y), (xs[i + 1] - 0.38, row_y),
                arrowstyle="-|>", mutation_scale=10,
                linewidth=0.9, color=ARROW_COLOR, zorder=3,
            ))

    # 外部 substrate / product（和单链一致的外部格式）
    substrate_pos = None
    product_pos = None
    if external_chems and cell_x0 is not None and cell_x1 is not None:
        chain_sub = steps[chain[0][0]].get("substrate")
        chain_prod = steps[chain[-1][-1]].get("product")
        if chain_sub:
            sub_x = cell_x0 - 1.1
            substrate_pos = _chem_outside(ax, sub_x, row_y, str(chain_sub))
            ax.add_patch(FancyArrowPatch(
                (sub_x + 0.5, row_y), (cell_x0 - 0.05, row_y),
                arrowstyle="-|>", mutation_scale=14,
                linewidth=1.4, color=OUTSIDE_ARROW_COLOR, zorder=4,
            ))
        if chain_prod:
            prod_x = cell_x1 + 1.1
            product_pos = _chem_outside(ax, prod_x, row_y, str(chain_prod))
            ax.add_patch(FancyArrowPatch(
                (cell_x1 + 0.05, row_y), (prod_x - 0.5, row_y),
                arrowstyle="-|>", mutation_scale=14,
                linewidth=1.4, color=OUTSIDE_ARROW_COLOR, zorder=4,
            ))
    return substrate_pos, product_pos


def _split_into_chains(segments: list[list[int]],
                       steps: list[dict[str, Any]]) -> list[list[list[int]]]:
    """把段列表分成链（chain）：连续段之间若化学物衔接则同链，否则新链。

    例：segments=[[0:sqr], [1..5:Sox]]，steps[0].product=S-0 vs
        steps[1].substrate=S2O3-2 → 两段不连通 → 两条独立链。
    """
    if not segments:
        return []
    chains: list[list[list[int]]] = [[segments[0]]]
    for i in range(1, len(segments)):
        prev_prod = str(steps[segments[i - 1][-1]].get("product") or "")
        cur_sub = str(steps[segments[i][0]].get("substrate") or "")
        if prev_prod and cur_sub and prev_prod == cur_sub:
            chains[-1].append(segments[i])
        else:
            chains.append([segments[i]])
    return chains


def _segment_by_complex(steps: list[dict[str, Any]]) -> list[list[int]]:
    """把 steps 按连续相同非空 complex 分段。

    返回 [[step_idx, ...], ...]。每个子列表是一个段：
    - 单元素段：该 step 的 complex 为 None 或与前后不同
    - 多元素段：连续 step 共享相同 complex（多亚基酶复合体）

    例：steps = [(sqr,None), (soxA,M595), (soxB,M595), (dummy,None), (foo,M999)]
        → [[0], [1, 2], [3], [4]]
    """
    out: list[list[int]] = []
    cur: list[int] = []
    cur_cx: str | None = None
    for i, s in enumerate(steps):
        cx = s.get("complex")
        if cx is not None and cx == cur_cx:
            cur.append(i)
        else:
            if cur:
                out.append(cur)
            cur = [i]
            cur_cx = cx
    if cur:
        out.append(cur)
    return out


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
    title_y = cy + cell_h / 2 - pad * 1.2
    # S2.5-8：若标题带 "★ " 前缀，红色渲染该星，其余黑色
    if title.startswith("★ "):
        ax.text(cell_x0 + pad, title_y, "★",
                fontsize=11, fontweight="bold", color="#C0392B", zorder=8)
        ax.text(cell_x0 + pad + 0.35, title_y, title[2:],
                fontsize=8.5, fontweight="bold", color="#333", zorder=8)
    else:
        ax.text(cell_x0 + pad, title_y, title,
                fontsize=8.5, fontweight="bold", color="#333", zorder=8)
    if mag_label:
        # keystone 标记：sentinel " ✦" → 金色四角星（DejaVu Sans 支持）
        if mag_label.endswith(" ✦"):
            plain = mag_label[:-2]
            ax.text(cell_x1 - pad, title_y, f"[{plain}]",
                    fontsize=7, fontstyle="italic", color="#555",
                    ha="right", zorder=8)
            ax.text(cell_x1 - pad - len(plain) * 0.08 - 0.32,
                    title_y, "✦",
                    fontsize=9, fontweight="bold", color="#D4A017",
                    ha="right", zorder=8)
        else:
            ax.text(cell_x1 - pad, title_y, f"[{mag_label}]",
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
    # —— 布局判据（三分枝） ————————————————————————————————
    # parallel_complex: 所有 substrate 相同 AND 所有 product 相同
    #   → 多亚基 / 同工酶做同一反应；不画内部箭头
    # all_same_intermediate: step[i].product 全一致但 substrate 不一致
    #   → 退化的折叠显示
    # else: 典型级联（denitrification 等）
    subs = [str(s.get("substrate") or "") for s in steps]
    prods = [str(s.get("product") or "") for s in steps]
    parallel_complex = (
        n_steps >= 2
        and len(set(subs)) == 1 and subs[0] != ""
        and len(set(prods)) == 1 and prods[0] != ""
    )
    inner_products = prods[:-1]   # 除最后一步外的中间产物
    all_same_intermediate = (
        not parallel_complex
        and len(inner_products) >= 1
        and len(set(inner_products)) == 1
        and inner_products[0] != ""
    )

    inner_pad = 0.5
    inner_x0 = cell_x0 + inner_pad
    inner_x1 = cell_x1 - inner_pad
    gy = cy - (cell_h * 0.15)

    gene_positions: list[tuple[float, float]] = []
    if parallel_complex:
        # 整条通路是一个并联簇（多亚基 or 同工酶）→ 合并成单大椭圆
        mid_x = (inner_x0 + inner_x1) / 2
        label = _bundle_label(steps)
        bundle_color = steps[0].get("color") or element_color
        _draw_bundle_ellipse(ax, mid_x, gy, label, bundle_color, n_genes=n_steps)
        gene_positions.append((mid_x, gy))
        # subtitle：有共同 complex → "(M00530 complex)"；无 → "(N isozymes)"
        complexes = {s.get("complex") for s in steps if s.get("complex")}
        if len(complexes) == 1:
            subtitle = f"({next(iter(complexes))} complex, {n_steps} subunits)"
        else:
            subtitle = f"({n_steps} isozymes)"
        ax.text(mid_x, gy - cell_h * 0.30,
                subtitle,
                ha="center", va="center",
                fontsize=6.5, fontstyle="italic",
                color="#666", zorder=7)
    elif all_same_intermediate and n_steps >= 2:
        # 酶紧挨排成一行；中间产物仅显示一次（在酶群中央上方稍小字）
        xs_genes = list(np.linspace(inner_x0, inner_x1, n_steps))
        for x, s in zip(xs_genes, steps):
            _gene(ax, x, gy, name=s.get("gene", "?"),
                  color=s.get("color") or element_color,
                  corr=s.get("corr"), show_hm=show_heatmap)
            gene_positions.append((x, gy))
        mid_x = (inner_x0 + inner_x1) / 2
        ax.text(mid_x, gy + cell_h * 0.28,
                "→ " + _pretty_formula(inner_products[0]) + " →",
                ha="center", va="center",
                fontsize=7, fontstyle="italic",
                color=INTERMEDIATE_COLOR, zorder=7)
        for i in range(len(xs_genes) - 1):
            ax.add_patch(FancyArrowPatch(
                (xs_genes[i] + 0.38, gy), (xs_genes[i + 1] - 0.38, gy),
                arrowstyle="-|>", mutation_scale=8,
                linewidth=0.8, color=ARROW_COLOR, zorder=3,
            ))
    else:
        # S2.5-10d: 按 complex 分段；再按化学物衔接分"链"
        segments = _segment_by_complex(steps)
        chains = _split_into_chains(segments, steps)
        if len(chains) == 1:
            _draw_chain_row(
                ax, chains[0], steps,
                inner_x0=inner_x0, inner_x1=inner_x1, row_y=gy,
                cell_h=cell_h, element_color=element_color,
                show_heatmap=show_heatmap,
                external_chems=False,
                gene_positions=gene_positions,
            )
        else:
            # 多链：每链一行；底物/产物画在 cell 外部（和单链一致格式）
            n_chains = len(chains)
            row_h = (cell_h * 0.92) / n_chains
            ys = [cy + (n_chains - 1) / 2 * row_h - i * row_h
                   for i in range(n_chains)]
            multichain_sub = None
            multichain_prod = None
            multichain_sub_map: dict[str, tuple[float, float]] = {}
            multichain_prod_map: dict[str, tuple[float, float]] = {}
            for ci, (chain, row_y) in enumerate(zip(chains, ys)):
                sub_pos, prod_pos = _draw_chain_row(
                    ax, chain, steps,
                    inner_x0=inner_x0, inner_x1=inner_x1, row_y=row_y,
                    cell_h=row_h * 0.9, element_color=element_color,
                    show_heatmap=False,
                    external_chems=True,
                    cell_x0=cell_x0, cell_x1=cell_x1,
                    gene_positions=gene_positions,
                )
                chain_sub_sp = steps[chain[0][0]].get("substrate")
                chain_prod_sp = steps[chain[-1][-1]].get("product")
                if sub_pos and chain_sub_sp:
                    multichain_sub_map.setdefault(str(chain_sub_sp), sub_pos)
                if prod_pos and chain_prod_sp:
                    multichain_prod_map.setdefault(str(chain_prod_sp), prod_pos)
                # fallback：优先记录首个带底物 / 末个带产物 chain 的位置
                if sub_pos and multichain_sub is None:
                    multichain_sub = sub_pos
                if prod_pos:
                    multichain_prod = prod_pos

    substrate_pos = None
    product_pos = None
    substrate_pos_map: dict[str, tuple[float, float]] = {}
    product_pos_map: dict[str, tuple[float, float]] = {}
    # 多链：已在 _draw_chain_row(external_chems=True) 中画了外部 substrate/product
    is_multi_chain = False
    if not (parallel_complex or (all_same_intermediate and n_steps >= 2)):
        _segs = _segment_by_complex(steps)
        _chains = _split_into_chains(_segs, steps)
        is_multi_chain = len(_chains) > 1
        if is_multi_chain:
            substrate_pos = locals().get("multichain_sub")
            product_pos = locals().get("multichain_prod")
            substrate_pos_map = dict(locals().get("multichain_sub_map") or {})
            product_pos_map = dict(locals().get("multichain_prod_map") or {})

    if show_outside_chems and not is_multi_chain:
        # 取首个非空 substrate（跳过纯调控基因 fur/tonB 的 None substrate）
        sub_txt = next((s.get("substrate") for s in steps if s.get("substrate")), None)
        if sub_txt:
            sub_x = cell_x0 - 1.1
            substrate_pos = _chem_outside(ax, sub_x, cy, str(sub_txt))
            substrate_pos_map[str(sub_txt)] = substrate_pos
            ax.add_patch(FancyArrowPatch(
                (sub_x + 0.5, cy), (cell_x0 - 0.05, cy),
                arrowstyle="-|>", mutation_scale=14,
                linewidth=1.4, color=OUTSIDE_ARROW_COLOR, zorder=4,
            ))
        # 取末个非空 product
        final = next((s.get("product") for s in reversed(steps) if s.get("product")), None)
        if final:
            prod_x = cell_x1 + 1.1
            product_pos = _chem_outside(ax, prod_x, cy, str(final))
            product_pos_map[str(final)] = product_pos
            ax.add_patch(FancyArrowPatch(
                (cell_x1 + 0.05, cy), (prod_x - 0.5, cy),
                arrowstyle="-|>", mutation_scale=14,
                linewidth=1.4, color=OUTSIDE_ARROW_COLOR, zorder=4,
            ))

    return {
        "substrate_pos": substrate_pos,
        "product_pos": product_pos,
        "substrate_pos_map": substrate_pos_map,
        "product_pos_map": product_pos_map,
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
            "complex": g.get("complex"),   # S2.5-10 KEGG MODULE
        }
        out.append(step)
    return out
