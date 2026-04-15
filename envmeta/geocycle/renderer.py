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

from envmeta.geocycle.cell_renderer import (
    _pretty_formula, _segment_by_complex, _split_into_chains,
    draw_cascade_cell, genes_to_steps,
)
from envmeta.geocycle.knowledge_base import couplings as kb_couplings
from envmeta.geocycle.model import CycleData, ElementCycle, PathwayActivity

DEFAULTS = {
    "width_mm": 460,
    "height_mm": 320,
    "show_env_panel": True,
    "max_cells_per_element": 3,           # 每元素象限最多画几个细胞
    "max_pathways_per_element": 8,        # v1 回退模式用
    "max_contributors_shown": 3,          # v1 回退模式用
    "cell_mode": "cascade",               # "cascade"（v2，默认）/ "bars"（v1 回退）
    "show_couplings": True,               # 画跨元素化学物耦合线（S2.5-3）
    "most_active_pathways": set(),        # S2.5-8：当前组里跨组对比最活的 pathway_id
    "show_inactive_pathways": True,       # S2.5-10d：象限底部列 KB 有但数据里 n=0 的通路名
    "title": "Biogeochemical Cycle Diagram (v2)",
    "cell_height_ratio": 0.22,            # 单细胞占象限高的比例
}


# 化学物名标准化（KB 里写法 vs cell 实际绘制值）—— 统一到同一键做匹配
_SPECIES_ALIASES = {
    "sulfide":     "S-2",
    "H2S":         "S-2",
    "H₂S":         "S-2",
    "SO4-2":       "SO4-2",
    "SO4^2-":      "SO4-2",
    "As(iii)":     "As(III)",
    "As(v)":       "As(V)",
    "Fe(iii)":     "Fe(III)",
    "Fe(ii)":      "Fe(II)",
}


def _norm_species(s: str | None) -> str | None:
    if not s:
        return None
    s = str(s).strip()
    return _SPECIES_ALIASES.get(s, s)


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
    top, bot = 8.4, 0.6
    available = top - bot
    # 预计算每 cell 期望高度：多链 cell（如 sqr + Sox）给 1.6× 基础高
    base_cell_h = 1.4
    weights: list[float] = []
    for pw, contrib in picked:
        steps_preview = genes_to_steps(contrib.genes, default_color=ec.color)
        _segs = _segment_by_complex(steps_preview)
        _chains = _split_into_chains(_segs, steps_preview)
        weights.append(1.6 if len(_chains) >= 2 else 1.0)
    total_weight = sum(weights) or 1.0
    # 缩放到可用高度
    cell_hs = [base_cell_h * w for w in weights]
    total_h = sum(cell_hs)
    if total_h > available * 0.95:
        scale = (available * 0.95) / total_h
        cell_hs = [h * scale for h in cell_hs]
    # 每 cell 纵向中点 cy
    gap = max(0.05, (available - sum(cell_hs)) / max(n + 1, 1))
    ys: list[float] = []
    y_cur = top - gap
    for h in cell_hs:
        ys.append(y_cur - h / 2)
        y_cur -= h + gap

    most_active = set(cfg.get("most_active_pathways") or ())
    anchors: list[dict] = []
    for (pw, contrib), cy, this_cell_h in zip(picked, ys, cell_hs):
        n_genes = len(contrib.genes)
        # 细胞宽随基因数伸缩（2 基因≈5 宽；4 基因≈8.5 宽；上限 10）
        cell_w = float(np.clip(3.5 + n_genes * 1.4, 4.5, 10.0))
        cell_x0 = (18 - cell_w) / 2

        steps = genes_to_steps(
            contrib.genes,
            default_color=ec.color,
        )
        # S2.5-8 通路 ★（跨组对比最活）：前缀红星
        prefix = "★ " if pw.pathway_id in most_active else ""
        title = prefix + pw.display_name
        mag_label = contrib.label or contrib.mag
        # S2.5-8 MAG 关键物种标记：label 尾部加 " ✦" sentinel（cell_renderer 识别后着金色）
        if getattr(contrib, "is_keystone", False):
            mag_label = mag_label + " ✦"
        r = draw_cascade_cell(
            ax,
            cy=cy, cell_x0=cell_x0, cell_w=cell_w, cell_h=this_cell_h,
            title=title, mag_label=mag_label,
            steps=steps,
            element_color=ec.color,
            show_heatmap=False,   # 上方 mini 热图留待 S2.5-3（需 KO-env 相关性数据）
            show_outside_chems=True,
        )
        # 小信息条：通路完整度 + 贡献度（右下）
        ax.text(
            17.7, cy - this_cell_h / 2 - 0.08,
            f"c={contrib.completeness:.0f}%  ab̄={contrib.abundance_mean:.2f}",
            fontsize=6.5, color="#555", ha="right", va="top",
        )
        # 化学物种名（用于跨元素耦合匹配）
        sub_species = _norm_species(contrib.genes[0].get("substrate")) if contrib.genes else None
        prod_species = _norm_species(contrib.genes[-1].get("product")) if contrib.genes else None
        r.update({
            "pathway_id": pw.pathway_id,
            "pathway_name": pw.display_name,
            "mag": contrib.mag,
            "element": ec.element_id,
            "ax": ax,
            "substrate_species": sub_species,
            "product_species": prod_species,
        })
        anchors.append(r)

    # S2.5-10d.3：象限底部列出 KB 里有但数据里完全没承载者的通路
    if cfg.get("show_inactive_pathways", True):
        inactive = [pw for pw in ec.pathways if pw.n_active_mags == 0]
        if inactive:
            names = ", ".join(p.display_name for p in inactive)
            # 长文本自动换行（大约每 50 字符一行）
            if len(names) > 45:
                import textwrap
                names = "\n  ".join(textwrap.wrap(names, width=45))
            ax.text(
                0.3, 0.35,
                f"[inactive: {names}]",
                fontsize=6.5, color="#999", fontstyle="italic",
                va="bottom", ha="left",
            )
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
# 化学物耦合连线（S2.5-3）
# =============================================================================

def _anchor_fig_point(fig, anchor: dict, which: str) -> tuple[float, float] | None:
    """把 anchor 的 substrate/product 数据坐标转 fig 坐标。

    which ∈ {"substrate_pos", "product_pos"}
    """
    key = which
    pos = anchor.get(key)
    ax = anchor.get("ax")
    if pos is None or ax is None:
        return None
    x, y = pos
    try:
        fx, fy = fig.transFigure.inverted().transform(
            ax.transData.transform((x, y)))
        return (float(fx), float(fy))
    except Exception:
        return None


def _find_best_pair(
    anchors: list[dict], species_a: str, species_b: str,
) -> tuple[dict, str, dict, str] | None:
    """在 anchors 中找一对"产物匹配 species_a / species_b"的细胞。

    优先：两端都是产物（product_species）；若找不到，退化为 substrate。
    元素不同时优先（跨元素耦合更有意义）。
    返回 (anchor_a, pos_key_a, anchor_b, pos_key_b) 或 None。
    """
    a_norm = _norm_species(species_a)
    b_norm = _norm_species(species_b)

    def _candidates(species_norm: str) -> list[tuple[dict, str]]:
        out: list[tuple[dict, str]] = []
        for a in anchors:
            if a.get("product_species") == species_norm and a.get("product_pos"):
                out.append((a, "product_pos"))
        for a in anchors:
            if a.get("substrate_species") == species_norm and a.get("substrate_pos"):
                # 产物优先（已加在前），这里做 substrate 兜底，但避免重复
                if not any(x[0] is a and x[1] == "product_pos" for x in out):
                    out.append((a, "substrate_pos"))
        return out

    cand_a = _candidates(a_norm)
    cand_b = _candidates(b_norm)
    if not cand_a or not cand_b:
        return None

    # 优先挑"不同元素"的一对
    for (aa, ka) in cand_a:
        for (ab, kb) in cand_b:
            if aa is ab:
                continue
            if aa.get("element") != ab.get("element"):
                return (aa, ka, ab, kb)
    # 无跨元素时退化
    for (aa, ka) in cand_a:
        for (ab, kb) in cand_b:
            if aa is ab:
                continue
            return (aa, ka, ab, kb)
    return None


def _redox_label(cp: dict, species_a: str, species_b: str) -> str:
    """redox 型耦合用"substrate→product"呈现；其他类型用 product 本身。"""
    product = cp.get("product", "") or ""
    ctype = (cp.get("type") or "").lower()
    if ctype == "redox":
        # 判断哪个是反应物（被氧化/还原）；约定 species_b 通常是被作用的
        # 但 KB 里 NO3-+As(III)→As(V) 写法是 (NO3-, As(III))→As(V)，As(III) 被氧化
        # 用启发式：被变化的物种（species_a/b 之一）包含在产物名里视作"终态"，
        # 另一个 = 起态
        reactant = species_b if species_b in product else species_a
        return (_pretty_formula(reactant) + r"$\rightarrow$"
                + _pretty_formula(product))
    return _pretty_formula(product)


def _deconflict_midpoints(raw_items: list[dict], threshold: float = 0.08,
                          offset: float = 0.05) -> list[dict]:
    """对 midpoint 距离 <threshold 的耦合做垂直错开。

    raw_items: [{"p1", "p2", "mid", ...}, ...]
    原地更新 mid；对重叠组按索引奇偶上/下偏移 offset。
    """
    n = len(raw_items)
    # 简化的 N^2 扫描（耦合条数 ≤10 足够）
    for i in range(n):
        overlap_idxs: list[int] = []
        mxi, myi = raw_items[i]["mid"]
        for j in range(n):
            if i == j:
                continue
            mxj, myj = raw_items[j]["mid"]
            d = ((mxi - mxj) ** 2 + (myi - myj) ** 2) ** 0.5
            if d < threshold:
                overlap_idxs.append(j)
        if overlap_idxs:
            # 按 i 在全部重叠组里的位置决定偏移方向 / 幅度
            group = sorted(set([i] + overlap_idxs))
            rank = group.index(i)
            # rank 0: -offset；1: +offset；2: -2*offset；...
            sign = 1 if rank % 2 else -1
            magnitude = ((rank + 1) // 2) * offset + offset * 0.5
            mx, my = raw_items[i]["mid"]
            raw_items[i]["mid"] = (mx, my + sign * magnitude)
    return raw_items


def _draw_couplings(fig, anchors: list[dict],
                    couplings_list: list[dict]) -> list[dict]:
    """跨象限画耦合线：化学物 A ─┐ (产物节点) ┌─ 化学物 B。

    S2.5-7c 改造：两阶段流程 —— 先计算所有 midpoints 并做去叠，再统一绘制；
    redox 类型产物节点显示 "reactant→product"。
    """
    # 阶段 1：收集 p1/p2/mid 不绘制
    raw: list[dict] = []
    for cp in couplings_list:
        a = cp.get("species_a")
        b = cp.get("species_b")
        if not a or not b:
            continue
        pair = _find_best_pair(anchors, a, b)
        if pair is None:
            continue
        anchor_a, key_a, anchor_b, key_b = pair
        p1 = _anchor_fig_point(fig, anchor_a, key_a)
        p2 = _anchor_fig_point(fig, anchor_b, key_b)
        if p1 is None or p2 is None:
            continue
        mx = (p1[0] + p2[0]) / 2
        my = (p1[1] + p2[1]) / 2
        raw.append({
            "coupling": cp,
            "species_a": a, "species_b": b,
            "p1": p1, "p2": p2, "mid": (mx, my),
        })

    _deconflict_midpoints(raw)

    # 阶段 2：绘制
    drawn: list[dict] = []
    for item in raw:
        cp = item["coupling"]
        color = cp.get("color", "#555")
        p1 = item["p1"]; p2 = item["p2"]; mid = item["mid"]
        label = _redox_label(cp, item["species_a"], item["species_b"])
        ctype = (cp.get("type") or "").lower()

        # 两段折线经过偏移后的 midpoint
        for (pa, pb) in [(p1, mid), (mid, p2)]:
            line = plt.Line2D([pa[0], pb[0]], [pa[1], pb[1]],
                              linestyle="--", linewidth=1.6,
                              color=color, alpha=0.75,
                              transform=fig.transFigure, zorder=1)
            fig.lines.append(line)
        # 产物节点
        fig.text(mid[0], mid[1], label, ha="center", va="center",
                 fontsize=7.5, fontweight="bold", color="white",
                 bbox=dict(facecolor=color, alpha=0.9, pad=4,
                           boxstyle="round,pad=0.3", edgecolor="black",
                           linewidth=0.8), zorder=12)
        # redox 在下方加 (ox) 提示
        if ctype == "redox":
            fig.text(mid[0], mid[1] - 0.018, "(ox)",
                     ha="center", va="top", fontsize=6.5,
                     fontstyle="italic", color=color, zorder=12)

        drawn.append({
            "species_a": item["species_a"], "species_b": item["species_b"],
            "product": cp.get("product", ""), "type": cp.get("type"),
            "color": color, "from": p1, "to": p2, "mid": mid,
            "label": label,
        })
    return drawn


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
    gf = data.meta.get("group_filter")
    group_tag = f"group={gf}  |  " if gf and str(gf).lower() not in ("all", "none") else ""
    rk = data.meta.get("contributor_ranking", "abundance")
    rank_tag = f"ranking={rk}  |  " if rk and rk != "abundance" else ""
    subtitle = (
        f"{group_tag}{rank_tag}"
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

    # 化学物耦合（跨象限连线）— 在 env 面板之前画，以免虚线压到条形
    drawn_couplings: list[dict] = []
    if mode == "cascade" and cfg.get("show_couplings", True) and all_anchors:
        drawn_couplings = _draw_couplings(fig, all_anchors, kb_couplings())

    if has_env:
        _draw_env_panel(fig, data, bbox=[left, 0.03, right - left, 0.20])

    # 附加到 figure 便于上层（S2.5-3 化学物耦合）检索细胞锚点
    fig._envmeta_cycle_anchors = all_anchors   # type: ignore[attr-defined]
    fig._envmeta_cycle_couplings = drawn_couplings   # type: ignore[attr-defined]
    return fig
