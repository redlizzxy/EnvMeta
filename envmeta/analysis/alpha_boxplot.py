"""α 多样性箱线图。

对标 scripts/R/02_alpha_diversity.R（单源简化版）。

输入：
    alpha_df:    α 多样性指数表，含 SampleID/Sample_ID 列 + 多个数值指数列
    metadata_df: 样本分组表，含 SampleID + Group

输出（AnalysisResult）：
    figure — 多子图箱线图（每指数一个子图）
    stats  — 扁平表：metric, group_a, group_b, n_a, n_b, stat, p, padj, kw_p
    params — 实际使用的参数
"""
from __future__ import annotations

from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from statsmodels.stats.multitest import multipletests

from envmeta.analysis.base import AnalysisResult

DEFAULT_PALETTE = {"CK": "#1c9cbd", "A": "#e3943d", "B": "#92181e"}

DEFAULTS = {
    "width_mm": 180,
    "height_mm": 120,
    "metrics": None,             # None → 自动选 alpha 表里所有数值列
    "group_col": "Group",
    "sample_col": None,          # None → 自动从 alpha 表推断（SampleID / Sample_ID）
    "group_order": None,
    "palette": None,             # None → 用 DEFAULT_PALETTE，组不在其中回退为自动色
    "show_points": True,
    "show_pvalues": True,
    "alpha": 0.05,
    "n_cols": 3,
}


def _find_sample_col(df: pd.DataFrame) -> str:
    for c in df.columns:
        if c.lower().replace("_", "").replace(" ", "") in ("sampleid",):
            return c
    raise ValueError(f"alpha 表未找到 SampleID 列（候选：{list(df.columns)}）")


def _auto_metrics(merged: pd.DataFrame, group_col: str, sample_col: str) -> list[str]:
    skip = {group_col.lower(), sample_col.lower(), "replicate"}
    out = []
    for c in merged.columns:
        if c.lower() in skip:
            continue
        if pd.to_numeric(merged[c], errors="coerce").notna().mean() > 0.95:
            out.append(c)
    return out


def _sig_label(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""


def analyze(alpha_df: pd.DataFrame, metadata_df: pd.DataFrame,
            params: dict | None = None) -> AnalysisResult:
    p = {**DEFAULTS, **(params or {})}

    sample_col = p["sample_col"] or _find_sample_col(alpha_df)
    group_col = p["group_col"]

    alpha = alpha_df.copy()
    alpha[sample_col] = alpha[sample_col].astype(str).str.strip()
    meta = metadata_df.copy()
    meta_sample_col = "SampleID" if "SampleID" in meta.columns else "Sample_ID"
    meta[meta_sample_col] = meta[meta_sample_col].astype(str).str.strip()

    merged = alpha.merge(meta[[meta_sample_col, group_col]],
                         left_on=sample_col, right_on=meta_sample_col, how="inner")
    if merged.empty:
        raise ValueError("alpha 表与 metadata 样本 ID 无交集")

    metrics = p["metrics"] or _auto_metrics(merged, group_col, sample_col)
    if not metrics:
        raise ValueError("没有可画的数值指数列")

    group_order = p["group_order"] or list(dict.fromkeys(meta[group_col].tolist()))
    groups_in_data = [g for g in group_order if g in set(merged[group_col])]
    if len(groups_in_data) < 2:
        raise ValueError(f"至少需要 2 个组才能做统计，现有 {groups_in_data}")

    palette = {**DEFAULT_PALETTE, **(p["palette"] or {})}
    default_colors = plt.cm.tab10.colors
    for i, g in enumerate(groups_in_data):
        palette.setdefault(g, default_colors[i % len(default_colors)])

    # === 统计 ===
    rows = []
    for metric in metrics:
        vals = merged[metric] = pd.to_numeric(merged[metric], errors="coerce")
        per_group = {g: vals[merged[group_col] == g].dropna().values
                     for g in groups_in_data}
        if all(len(v) == 0 for v in per_group.values()):
            continue
        try:
            kw_stat, kw_p = sp_stats.kruskal(*[v for v in per_group.values() if len(v) > 0])
        except Exception:
            kw_p = np.nan

        pair_rows = []
        for a, b in combinations(groups_in_data, 2):
            va, vb = per_group[a], per_group[b]
            if len(va) < 1 or len(vb) < 1:
                continue
            try:
                u, pv = sp_stats.mannwhitneyu(va, vb, alternative="two-sided")
            except ValueError:
                u, pv = np.nan, 1.0
            pair_rows.append({
                "metric": metric, "group_a": a, "group_b": b,
                "n_a": len(va), "n_b": len(vb),
                "stat": u, "p": pv, "kw_p": kw_p,
            })
        if pair_rows:
            ps = np.array([r["p"] for r in pair_rows], dtype=float)
            _, padj, _, _ = multipletests(ps, method="fdr_bh")
            for r, pa in zip(pair_rows, padj):
                r["padj"] = pa
                r["significance"] = _sig_label(pa)
            rows.extend(pair_rows)

    stats_df = pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["metric", "group_a", "group_b", "n_a", "n_b",
                 "stat", "p", "padj", "kw_p", "significance"])

    # === 绘图 ===
    n = len(metrics)
    n_cols = min(p["n_cols"], n)
    n_rows = (n + n_cols - 1) // n_cols
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4 * n_rows),
        constrained_layout=True, squeeze=False,
    )

    for idx, metric in enumerate(metrics):
        ax = axes[idx // n_cols][idx % n_cols]
        data = [merged.loc[merged[group_col] == g, metric].dropna().values
                for g in groups_in_data]
        box = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False)
        for patch, g in zip(box["boxes"], groups_in_data):
            patch.set_facecolor(palette[g])
            patch.set_alpha(0.4)
            patch.set_edgecolor(palette[g])
        for median in box["medians"]:
            median.set_color("black")

        if p["show_points"]:
            for i, g in enumerate(groups_in_data):
                vals = data[i]
                jitter = np.random.RandomState(42 + i).normal(0, 0.06, size=len(vals))
                ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                           color=palette[g], s=18, edgecolor="white",
                           linewidth=0.5, zorder=3)

        ax.set_xticks(range(1, len(groups_in_data) + 1))
        ax.set_xticklabels(groups_in_data)
        ax.set_ylabel(metric)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

        if p["show_pvalues"] and rows:
            sub = stats_df[stats_df["metric"] == metric].sort_values("padj")
            y_max = np.max([v.max() if len(v) else 0 for v in data])
            y_min = np.min([v.min() if len(v) else 0 for v in data])
            gap = (y_max - y_min) * 0.08 if y_max > y_min else 1.0
            step = 0
            for _, r in sub.iterrows():
                if r["padj"] >= p["alpha"]:
                    continue
                i_a = groups_in_data.index(r["group_a"]) + 1
                i_b = groups_in_data.index(r["group_b"]) + 1
                y = y_max + gap * (1 + step)
                ax.plot([i_a, i_a, i_b, i_b], [y, y + gap * 0.3, y + gap * 0.3, y],
                        color="black", lw=0.7)
                ax.text((i_a + i_b) / 2, y + gap * 0.4, r["significance"],
                        ha="center", va="bottom", fontsize=10)
                step += 1

    for idx in range(n, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)

    return AnalysisResult(figure=fig, stats=stats_df, params=p)
