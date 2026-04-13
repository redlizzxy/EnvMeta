"""元素循环基因热图（4 元素合一版）。

输入：
    ko_abundance_df: 列含 `KEGG_ko` + 若干样本列（TPM 或丰度）
    metadata_df:     含 SampleID + Group

输出：
    figure — 一张热图（KO × Group），左侧有元素色块 + 通路色块
    stats  — Z-score 矩阵（KO × Group）
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from envmeta.analysis.base import AnalysisResult
from envmeta.geocycle.knowledge_base import (
    element_colors,
    element_display,
    element_pathway_ko_order,
    flat_ko_map,
)

DEFAULTS = {
    "width_mm": 180,
    "height_mm": 240,
    "cmap": "RdBu_r",
    "zscore": True,              # 按 KO（行）Z-score
    "group_col": "Group",
    "element_filter": None,      # None → 全部 4 元素；list → 子集
    "show_gene_names": True,
    "group_order": None,         # None → 用 metadata 中出现顺序
    "title_size": 11,
    "label_size": 9,
    "tick_size": 7,
}

# 通路色块用柔和配色，同一元素内通路颜色相近
_PATHWAY_COLORS_BY_ELEMENT = {
    "arsenic": ["#fdeaea", "#fcd0d0", "#f8a0a0", "#f47070", "#e74c3c", "#c0392b"],
    "nitrogen": ["#ebf5fb", "#d6eaf8", "#aed6f1", "#7fb3d5", "#5499c7", "#2e86c1"],
    "sulfur": ["#fef9e7", "#fcf3cf", "#f9e79f", "#f7dc6f", "#f1c40f", "#d4ac0d"],
    "iron": ["#fdf2e9", "#fae5d3", "#f5cba7", "#edbb99", "#e67e22", "#b9770e"],
}


# ==============================================================================
# 数据准备
# ==============================================================================

def _extract_ko_samples(df: pd.DataFrame) -> pd.DataFrame:
    """把 wide 格式（带 KEGG_ko 列 + 注释列 + 样本列）归一成 KO × 样本 的 DataFrame。"""
    ko_col = next((c for c in df.columns if c.lower() in ("kegg_ko", "kegg ko")), None)
    if ko_col is None:
        raise ValueError("输入缺少 KEGG_ko 列")

    # 数值列视为样本列
    numeric_cols = []
    for c in df.columns:
        if c == ko_col:
            continue
        nan_ratio = pd.to_numeric(df[c], errors="coerce").isna().mean()
        if nan_ratio < 0.05:
            numeric_cols.append(c)
    if len(numeric_cols) < 2:
        raise ValueError("找不到 ≥2 个数值样本列")

    ko_to_row = df[[ko_col] + numeric_cols].copy()
    for c in numeric_cols:
        ko_to_row[c] = pd.to_numeric(ko_to_row[c], errors="coerce").fillna(0.0)
    # 同一 KO 多行（如不同 EC 注释）取均值
    agg = ko_to_row.groupby(ko_col, as_index=True)[numeric_cols].mean()
    return agg


def _group_mean(df: pd.DataFrame, metadata: pd.DataFrame, group_col: str,
                group_order: list[str] | None) -> pd.DataFrame:
    mapping = dict(zip(metadata["SampleID"].astype(str), metadata[group_col].astype(str)))
    sample_to_group = {s: mapping.get(s) for s in df.columns if s in mapping}
    known = [s for s, g in sample_to_group.items() if g is not None]
    if not known:
        raise ValueError("metadata 中没有匹配到任何样本 ID")
    groups = pd.Series({s: sample_to_group[s] for s in known})
    means = df[known].T.groupby(groups).mean().T
    if group_order is None:
        seen = []
        for g in metadata[group_col].astype(str):
            if g not in seen and g in means.columns:
                seen.append(g)
        group_order = seen
    return means[[g for g in group_order if g in means.columns]]


def _zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    def z(row):
        mu, sigma = row.mean(), row.std()
        if sigma == 0 or np.isnan(sigma):
            return row * 0.0
        return (row - mu) / sigma
    return df.apply(z, axis=1)


# ==============================================================================
# 绘图
# ==============================================================================

def _assign_pathway_colors(present_rows: list[tuple[str, str, str]]) -> dict[tuple[str, str], str]:
    """为 (element, pathway) 分配颜色。"""
    colors: dict[tuple[str, str], str] = {}
    by_element: dict[str, list[str]] = {}
    for el, pw, _ in present_rows:
        by_element.setdefault(el, [])
        if pw not in by_element[el]:
            by_element[el].append(pw)
    for el, pws in by_element.items():
        palette = _PATHWAY_COLORS_BY_ELEMENT.get(el, ["#cccccc"] * len(pws))
        for i, pw in enumerate(pws):
            colors[(el, pw)] = palette[min(i, len(palette) - 1)]
    return colors


def _build_figure(z_df: pd.DataFrame, row_meta: list[tuple[str, str, str]],
                  *, width_mm: int, height_mm: int, cmap: str,
                  show_gene_names: bool,
                  title_size: int, label_size: int, tick_size: int) -> plt.Figure:
    """绘热图 + 左侧元素/通路色块 + 右侧基因名。"""
    kb_flat = flat_ko_map()
    el_colors = element_colors()
    el_display = element_display(lang="en")
    pw_colors = _assign_pathway_colors(row_meta)

    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), constrained_layout=False)
    # 左：两列窄 strip（通路 + 元素），中：热图，右：基因名
    gs = fig.add_gridspec(
        nrows=1, ncols=4,
        width_ratios=[0.3, 0.3, 5.0, 1.2],
        wspace=0.05, left=0.03, right=0.96, top=0.92, bottom=0.1,
    )
    ax_el = fig.add_subplot(gs[0, 0])
    ax_pw = fig.add_subplot(gs[0, 1], sharey=ax_el)
    ax_heat = fig.add_subplot(gs[0, 2], sharey=ax_el)
    ax_gene = fig.add_subplot(gs[0, 3], sharey=ax_el)

    data = z_df.to_numpy()
    n_rows, n_cols = data.shape
    y_pos = np.arange(n_rows)

    # 热图
    vmax = np.nanmax(np.abs(data))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    im = ax_heat.imshow(data, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax,
                        extent=(-0.5, n_cols - 0.5, n_rows - 0.5, -0.5))
    ax_heat.set_xticks(np.arange(n_cols))
    ax_heat.set_xticklabels(z_df.columns, fontsize=label_size)
    ax_heat.set_yticks([])
    ax_heat.tick_params(top=False, bottom=True, left=False, right=False,
                        labeltop=False, labelbottom=True)
    for spine in ax_heat.spines.values():
        spine.set_linewidth(0.4)

    # 元素色块
    for i, (el, _, _) in enumerate(row_meta):
        ax_el.add_patch(mpatches.Rectangle((0, i - 0.5), 1, 1,
                                           facecolor=el_colors.get(el, "#ccc"),
                                           edgecolor="white", linewidth=0.3))
    ax_el.set_xlim(0, 1)
    ax_el.set_ylim(n_rows - 0.5, -0.5)
    ax_el.set_xticks([])
    ax_el.set_yticks([])
    for spine in ax_el.spines.values():
        spine.set_visible(False)

    # 通路色块
    for i, (el, pw, _) in enumerate(row_meta):
        ax_pw.add_patch(mpatches.Rectangle((0, i - 0.5), 1, 1,
                                           facecolor=pw_colors.get((el, pw), "#ddd"),
                                           edgecolor="white", linewidth=0.3))
    ax_pw.set_xlim(0, 1)
    ax_pw.set_ylim(n_rows - 0.5, -0.5)
    ax_pw.set_xticks([])
    ax_pw.set_yticks([])
    for spine in ax_pw.spines.values():
        spine.set_visible(False)

    # 基因名右栏
    ax_gene.set_xlim(0, 1)
    ax_gene.set_ylim(n_rows - 0.5, -0.5)
    ax_gene.set_xticks([])
    ax_gene.set_yticks([])
    for spine in ax_gene.spines.values():
        spine.set_visible(False)
    if show_gene_names:
        for i, (_, _, ko) in enumerate(row_meta):
            gene, _, _ = kb_flat.get(ko, (ko, "?", "?"))
            ax_gene.text(0.02, i, f"{ko}  {gene}", fontsize=tick_size,
                         va="center", ha="left", style="italic", color="#333")

    # 元素分割线（在通路色块右侧跨到热图，分隔元素块）
    prev_el = None
    for i, (el, _, _) in enumerate(row_meta):
        if prev_el is not None and el != prev_el:
            ax_heat.axhline(i - 0.5, color="white", linewidth=1.2)
        prev_el = el

    # 标题 + 色块图例（元素）
    ax_heat.set_title("Element-cycling gene abundance (row-wise Z-score)",
                      fontsize=title_size, pad=8)
    legend_handles = [
        mpatches.Patch(facecolor=el_colors[eid], edgecolor="white",
                       label=el_display.get(eid, eid))
        for eid in el_colors if any(rm[0] == eid for rm in row_meta)
    ]
    ax_heat.legend(handles=legend_handles, loc="upper left",
                   bbox_to_anchor=(1.15, -0.02), frameon=False,
                   fontsize=tick_size, ncol=min(4, len(legend_handles)))

    # 颜色 bar
    cbar = fig.colorbar(im, ax=ax_heat, shrink=0.5, pad=0.12, location="right")
    cbar.set_label("Z-score", fontsize=label_size)
    cbar.ax.tick_params(labelsize=tick_size)

    return fig


# ==============================================================================
# 公共入口
# ==============================================================================

def analyze(ko_abundance_df: pd.DataFrame, metadata_df: pd.DataFrame,
            params: dict | None = None) -> AnalysisResult:
    p = {**DEFAULTS, **(params or {})}

    ko_samples = _extract_ko_samples(ko_abundance_df)
    kb_order = element_pathway_ko_order()          # [(el, pw, ko), ...]
    if p["element_filter"]:
        keep_els = set(p["element_filter"])
        kb_order = [x for x in kb_order if x[0] in keep_els]

    # 保留 abundance 中存在的 KO，按 kb_order 排序
    present_order = [x for x in kb_order if x[2] in ko_samples.index]
    if not present_order:
        raise ValueError("输入中没有任何知识库已知 KO")
    ko_ordered = ko_samples.loc[[x[2] for x in present_order]]

    grouped = _group_mean(ko_ordered, metadata_df, p["group_col"], p["group_order"])
    z = _zscore_rows(grouped) if p["zscore"] else grouped

    fig = _build_figure(
        z, present_order,
        width_mm=p["width_mm"], height_mm=p["height_mm"], cmap=p["cmap"],
        show_gene_names=p["show_gene_names"],
        title_size=p["title_size"], label_size=p["label_size"], tick_size=p["tick_size"],
    )

    return AnalysisResult(figure=fig, stats=z,
                          params={**p, "_raw_means": grouped,
                                  "_row_meta": present_order})
