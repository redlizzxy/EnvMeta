"""物种组成堆叠图（Taxonomy stacked bar）。

对标 scripts/R/01_tax_stackplot.R，用 matplotlib 重写。

输入：
    abundance_df: 丰度表 (Taxonomy × 样本)，首列 "Taxonomy"，其余列为样本名
    metadata_df:  样本分组表，必含 "SampleID" + "Group"

输出（AnalysisResult）：
    figure — matplotlib Figure
    stats  — 百分比矩阵 (Top-N 物种 × 样本或分组)
    params — 实际使用的参数（用于后续代码生成器）
"""
from __future__ import annotations

import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from envmeta.analysis.base import AnalysisResult

# ==============================================================================
# 常量（迭代 2 迁移到 envmeta/params/presets.py）
# ==============================================================================

STACK_12 = [
    "#bbc4e4", "#d7e7af", "#bf9d6d", "#a4d6c1", "#bce2e8", "#b5b5b6",
    "#0084c2", "#3ab483", "#00978c", "#21825e", "#6aa3d2", "#3d62ad",
]
OTHERS_COLOR = "#cccccc"
UNCLASSIFIED_REGEX = re.compile(r"(?i)(unclassif|unassign|unknown)")

DEFAULTS = {
    "top_n": 10,
    "style": "sample",           # "sample" | "group"
    "group_col": "Group",
    "width_mm": 160,
    "height_mm": 100,
    "drop_unclassified": True,
    "palette": None,             # None → 用 STACK_12
    "others_label": "Others",
    "bar_width": 0.8,
    "show_legend": True,
    "sort_by": "mean",           # "mean" | "max" | "median"
    "reverse_stack": False,      # True → 高丰度在柱顶（默认在柱底）
}


_AGG_FUNCS = {"mean": "mean", "max": "max", "median": "median"}


# ==============================================================================
# 核心
# ==============================================================================

def _to_numeric_matrix(abundance_df: pd.DataFrame) -> pd.DataFrame:
    """首列做行索引，其余列转 float。"""
    df = abundance_df.copy()
    df = df.set_index(df.columns[0])
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return df


def _drop_unclassified(df: pd.DataFrame) -> pd.DataFrame:
    mask = df.index.to_series().apply(lambda s: bool(UNCLASSIFIED_REGEX.search(str(s))))
    return df.loc[~mask]


def _normalize_to_percent(df: pd.DataFrame) -> pd.DataFrame:
    """列和归一到 100%。"""
    col_sums = df.sum(axis=0).replace(0, np.nan)
    return df.div(col_sums, axis=1).fillna(0.0) * 100.0


def _topn_with_others(df: pd.DataFrame, top_n: int, others_label: str, sort_by: str) -> pd.DataFrame:
    """按给定聚合函数取 Top-N，其余聚合为 Others；仅当有剩余时才追加 Others 行。"""
    agg = _AGG_FUNCS.get(sort_by, "mean")
    if df.shape[0] <= top_n:
        ranked = getattr(df, agg)(axis=1).sort_values(ascending=False)
        return df.loc[ranked.index]     # 按排序返回（无 Others）
    ranked = getattr(df, agg)(axis=1).sort_values(ascending=False)
    top = df.loc[ranked.head(top_n).index]
    rest = df.loc[ranked.iloc[top_n:].index]
    top.loc[others_label] = rest.sum(axis=0)
    return top


def _aggregate_by_group(df: pd.DataFrame, metadata: pd.DataFrame, group_col: str) -> pd.DataFrame:
    """把样本列按 Group 取均值，返回 (Taxonomy × Group) 矩阵。"""
    mapping = dict(zip(metadata["SampleID"].astype(str), metadata[group_col].astype(str)))
    sample_to_group = {s: mapping.get(s) for s in df.columns if s in mapping}
    known = [s for s, g in sample_to_group.items() if g is not None]
    if not known:
        raise ValueError("metadata 中无法匹配到任何样本 ID")
    sub = df[known]
    groups = pd.Series(sample_to_group).reindex(known)
    agg = sub.T.groupby(groups).mean().T
    # 保留 metadata 里的组别顺序
    ordered = [g for g in metadata[group_col].drop_duplicates().tolist() if g in agg.columns]
    return agg[ordered]


def _build_figure(pct_df: pd.DataFrame, *, width_mm: int, height_mm: int,
                  palette: list[str], others_label: str,
                  bar_width: float, show_legend: bool,
                  reverse_stack: bool = False) -> plt.Figure:
    """绘制百分比堆叠柱状图。pct_df: Taxonomy × X（样本或分组）。
    reverse_stack=True 时高丰度项堆在柱顶。"""
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4), constrained_layout=True)

    # 绘图顺序：默认（reverse_stack=False）高丰度在底，所以从 df 顺序开始画（最先画的在底）
    taxa = pct_df.index.tolist()
    if reverse_stack:
        taxa = list(reversed(taxa))
    x_labels = pct_df.columns.tolist()
    x_pos = np.arange(len(x_labels))

    # 颜色：给非-Others 的 taxa 循环分配调色板，Others 用灰
    colors = {}
    palette_cycle = list(palette) if palette else STACK_12
    non_others = [t for t in taxa if t != others_label]
    for i, t in enumerate(non_others):
        colors[t] = palette_cycle[i % len(palette_cycle)]
    if others_label in taxa:
        colors[others_label] = OTHERS_COLOR

    bottom = np.zeros(len(x_labels))
    for t in taxa:
        vals = pct_df.loc[t].to_numpy()
        ax.bar(x_pos, vals, bar_width, bottom=bottom, color=colors[t], label=t,
               edgecolor="white", linewidth=0.3)
        bottom += vals

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, rotation=45, ha="right")
    ax.set_ylabel("Relative abundance (%)")
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.5, len(x_labels) - 0.5)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    if show_legend:
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
                  frameon=False, fontsize=8, title="Taxonomy")

    return fig


# ==============================================================================
# 公共入口
# ==============================================================================

def analyze(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    params: dict | None = None,
) -> AnalysisResult:
    p = {**DEFAULTS, **(params or {})}

    df = _to_numeric_matrix(abundance_df)
    if p["drop_unclassified"]:
        df = _drop_unclassified(df)
    df = _normalize_to_percent(df)

    if p["style"] == "group":
        df = _aggregate_by_group(df, metadata_df, p["group_col"])
    elif p["style"] == "sample":
        # 只保留 metadata 中出现的样本，并按其顺序
        known_samples = [s for s in metadata_df["SampleID"].astype(str).tolist() if s in df.columns]
        if known_samples:
            df = df[known_samples]
    else:
        raise ValueError(f"style 必须是 'sample' 或 'group'，收到 {p['style']!r}")

    pct = _topn_with_others(df, p["top_n"], p["others_label"], p["sort_by"])
    pct = _normalize_to_percent(pct)   # Top-N + Others 再归一一次防精度漂移

    fig = _build_figure(
        pct,
        width_mm=p["width_mm"], height_mm=p["height_mm"],
        palette=p["palette"], others_label=p["others_label"],
        bar_width=p["bar_width"], show_legend=p["show_legend"],
        reverse_stack=p["reverse_stack"],
    )

    return AnalysisResult(figure=fig, stats=pct, params=p)
