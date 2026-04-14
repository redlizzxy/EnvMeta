"""元素循环基因 log2FC 差异柱图。

对标 scripts/python/05_gene_heatmap_log2fc.py 后半段 (Fig2-9)。

输入：
    ko_abundance_df: KO × sample 丰度（含 KEGG_ko 列），e.g. eggnog TPM spf
    metadata_df:     SampleID + Group 表

输出：
    figure — 2×2（4 元素）子图，水平条形 log2FC + 显著性
    stats  — 扁平表：ko, gene, pathway, element, mean_a, mean_b, log2fc, t, p, padj, significance
"""
from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from statsmodels.stats.multitest import multipletests

from envmeta.analysis.base import AnalysisResult
from envmeta.geocycle.knowledge_base import (
    element_colors, element_display, flat_ko_map, load_kb,
)

DEFAULTS = {
    "width_mm": 220,
    "height_mm": 200,
    "group_col": "Group",
    "group_a": None,
    "group_b": None,
    "element_filter": None,
    "alpha": 0.05,
    "log2fc_threshold": 1.0,
    "pseudocount": 0.5,
    "color_up_sig": "#C0392B",
    "color_down_sig": "#2980B9",
    "color_up_nonsig": "#E8A0A0",
    "color_down_nonsig": "#A0B8D8",
    "show_gene_names": True,
    "n_cols": 2,
}


def _sig_label(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""


def _find_ko_col(df: pd.DataFrame) -> str:
    for c in df.columns:
        if c.lower() in ("kegg_ko", "kegg ko", "ko"):
            return c
    # 兜底：第二列
    return df.columns[1]


def analyze(ko_abundance_df: pd.DataFrame, metadata_df: pd.DataFrame,
            params: dict | None = None) -> AnalysisResult:
    p = {**DEFAULTS, **(params or {})}
    if not p["group_a"] or not p["group_b"]:
        raise ValueError("必须指定 group_a 和 group_b")
    if p["group_a"] == p["group_b"]:
        raise ValueError("group_a 和 group_b 不能相同")

    ko_map = flat_ko_map()
    kb = load_kb()
    colors_el = element_colors(kb)
    name_el = element_display(kb, "en")

    # 样本 → 组
    md = metadata_df.copy()
    sid_col = "SampleID" if "SampleID" in md.columns else "Sample_ID"
    md[sid_col] = md[sid_col].astype(str).str.strip()
    group_col = p["group_col"]
    groups = md.groupby(group_col)[sid_col].apply(list).to_dict()
    if p["group_a"] not in groups or p["group_b"] not in groups:
        raise ValueError(f"metadata 里找不到组 {p['group_a']} 或 {p['group_b']}")

    # 丰度表转数值矩阵（行 = KO，列 = 样本）
    ab = ko_abundance_df.copy()
    ko_col = _find_ko_col(ab)
    ab[ko_col] = ab[ko_col].astype(str).str.strip()
    # 保留已知 KO
    ab = ab[ab[ko_col].isin(ko_map)]
    sample_cols = [c for c in ab.columns if c in set(md[sid_col])]
    if not sample_cols:
        raise ValueError("丰度表的列名和 metadata 的 SampleID 无交集")
    ab_num = ab[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    ab_num.index = ab[ko_col].values
    # 聚合重复 KO 取均值
    ab_num = ab_num.groupby(level=0).mean()

    samples_a = [s for s in groups[p["group_a"]] if s in ab_num.columns]
    samples_b = [s for s in groups[p["group_b"]] if s in ab_num.columns]
    if len(samples_a) < 2 or len(samples_b) < 2:
        raise ValueError(f"每组至少 2 个样本（现 A={len(samples_a)}, B={len(samples_b)}）")

    # 统计
    pc = p["pseudocount"]
    rows = []
    for ko in ab_num.index:
        va = ab_num.loc[ko, samples_a].astype(float).values
        vb = ab_num.loc[ko, samples_b].astype(float).values
        mean_a = float(np.mean(va))
        mean_b = float(np.mean(vb))
        log2fc = math.log2((mean_a + pc) / (mean_b + pc))
        try:
            t, pv = sp_stats.ttest_ind(va, vb, equal_var=False)
            if np.isnan(pv):
                pv = 1.0
                t = 0.0
        except Exception:
            t, pv = 0.0, 1.0
        gene, pathway, element = ko_map[ko]
        rows.append({
            "ko": ko, "gene": gene, "pathway": pathway, "element": element,
            "mean_a": mean_a, "mean_b": mean_b, "log2fc": log2fc,
            "t": float(t), "p": float(pv),
        })

    stats_df = pd.DataFrame(rows)
    if stats_df.empty:
        raise ValueError("丰度表和知识库 KO 无交集")

    _, padj, _, _ = multipletests(stats_df["p"].values, method="fdr_bh")
    stats_df["padj"] = padj
    stats_df["significant"] = (stats_df["padj"] < p["alpha"]) & (
        stats_df["log2fc"].abs() >= p["log2fc_threshold"])
    stats_df["significance"] = stats_df["padj"].apply(_sig_label)

    # 过滤元素（同时影响 stats 和绘图）
    elements = p["element_filter"] or list(kb["elements"].keys())
    elements = [e for e in elements if e in set(stats_df["element"])]
    if not elements:
        raise ValueError("过滤后无元素可画")
    stats_df = stats_df[stats_df["element"].isin(elements)].reset_index(drop=True)

    # === 绘图 ===
    n = len(elements)
    n_cols = min(p["n_cols"], n)
    n_rows = (n + n_cols - 1) // n_cols
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4 * n_rows / 2),
        constrained_layout=True, squeeze=False,
    )

    for idx, el in enumerate(elements):
        ax = axes[idx // n_cols][idx % n_cols]
        sub = stats_df[stats_df["element"] == el].sort_values("log2fc")
        if sub.empty:
            ax.set_visible(False)
            continue

        y = np.arange(len(sub))
        colors = [
            (p["color_up_sig"] if lfc > 0 else p["color_down_sig"]) if sig else
            (p["color_up_nonsig"] if lfc > 0 else p["color_down_nonsig"])
            for lfc, sig in zip(sub["log2fc"], sub["significant"])
        ]
        ax.barh(y, sub["log2fc"].values, color=colors, edgecolor="white", linewidth=0.3)
        ax.axvline(0, color="black", lw=0.6)
        ax.axvline(p["log2fc_threshold"], color="gray", lw=0.5, ls="--")
        ax.axvline(-p["log2fc_threshold"], color="gray", lw=0.5, ls="--")

        ax.set_yticks(y)
        if p["show_gene_names"]:
            ax.set_yticklabels(sub["gene"].tolist(), fontsize=7, style="italic")
        else:
            ax.set_yticklabels(sub["ko"].tolist(), fontsize=6)
        ax.set_xlabel(f"log2({p['group_a']}/{p['group_b']})")

        # 元素色块标题
        ax.set_title(name_el.get(el, el), color=colors_el.get(el, "#333"),
                     fontsize=10, fontweight="bold")

        # 星号标注
        for yi, (_, r) in zip(y, sub.iterrows()):
            if r["significance"]:
                x = r["log2fc"]
                ha = "left" if x >= 0 else "right"
                offset = 0.1 if x >= 0 else -0.1
                ax.text(x + offset, yi, r["significance"],
                        ha=ha, va="center", fontsize=7, color="black")

        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    for idx in range(n, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)

    return AnalysisResult(figure=fig, stats=stats_df, params=p)
