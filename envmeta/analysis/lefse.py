"""LEfSe LDA 差异特征条形图。

对标 scripts/R/04_LEfSe.R（R 脚本只绘制 Galaxy LEfSe 预跑的 input.res；本模块
直接从丰度表做完整 LEfSe 流程：Kruskal-Wallis 筛选 + LDA 效应量 + 柱图）。

简化算法（Segata et al. 2011 思路）：
1. 每个 taxon 做 Kruskal-Wallis 检验（所有组），保留 p < alpha_kw 的特征
2. 对每个显著特征：
   - Enriched 组 = 组均值最大的组
   - LDA 效应量 = log10(1 + 1e6 × |max_mean − min_mean|)
3. 过滤 |LDA| > lda_threshold
4. 水平柱图，按组分段 + 组内按 LDA 降序排列

输入：
    abundance_df: 首列 Taxonomy + 样本列（可为 Species/Genus/Phylum 任一层级）
    metadata_df:  SampleID + Group

输出：
    figure — 水平 LDA 条形图（按组上色）
    stats  — feature, short_name, tax_level, group, lda, kw_p, mean_by_group...
"""
from __future__ import annotations

import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from envmeta.analysis.base import AnalysisResult

UNCLASSIFIED_REGEX = re.compile(r"(?i)(unclassif|unassign|unknown)")
TAX_PREFIX = {"k": "Kingdom", "p": "Phylum", "c": "Class", "o": "Order",
              "f": "Family", "g": "Genus", "s": "Species", "t": "Strain"}

DEFAULT_PALETTE = {"CK": "#1c9cbd", "A": "#e3943d", "B": "#92181e"}

DEFAULTS = {
    "width_mm": 180,
    "height_mm": 200,
    "alpha_kw": 0.05,
    "lda_threshold": 2.0,
    "group_col": "Group",
    "group_order": None,
    "palette": None,
    "tax_levels": None,       # None → 全部；e.g. ["Genus", "Species"]
    "max_features": None,     # None → 不截断
    "drop_unclassified": True,
    "lda_scale": 1e6,         # 把 abs(diff) 放大到 log10 ~2-6 范围
}


def _tax_level_from_name(name: str) -> str:
    """从 `s__Foo_bar` 或 `k__X.p__Y.g__Z` 提取最后一级的分类层级。"""
    last = str(name).split(".")[-1].split(";")[-1].strip()
    m = re.match(r"^([kpcofgst])__", last)
    if m:
        return TAX_PREFIX[m.group(1)]
    return "Unknown"


def _short_name(name: str) -> str:
    last = str(name).split(".")[-1].split(";")[-1].strip()
    return re.sub(r"^[kpcofgst]__", "", last)


def _find_sample_col(df: pd.DataFrame) -> str:
    for c in df.columns:
        if c.lower().replace("_", "") == "sampleid":
            return c
    raise ValueError(f"metadata 缺少 SampleID 列（候选：{list(df.columns)}）")


def analyze(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    params: dict | None = None,
) -> AnalysisResult:
    p = {**DEFAULTS, **(params or {})}
    sample_col = _find_sample_col(metadata_df)
    meta = metadata_df[[sample_col, p["group_col"]]].copy()
    meta[sample_col] = meta[sample_col].astype(str)

    # 丰度：行=taxon, 列=样本
    ab = abundance_df.copy()
    tax_col = ab.columns[0]
    ab = ab.set_index(tax_col)
    if p["drop_unclassified"]:
        mask = ab.index.to_series().apply(lambda s: not UNCLASSIFIED_REGEX.search(str(s)))
        ab = ab.loc[mask]
    ab = ab.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # 对齐样本
    shared = [c for c in ab.columns if c in set(meta[sample_col])]
    if len(shared) < 4:
        raise ValueError(f"丰度表与 metadata 共享样本过少（{len(shared)}），至少 4 个")
    ab = ab[shared]
    group_map = dict(zip(meta[sample_col], meta[p["group_col"]].astype(str)))
    groups_per_sample = np.array([group_map[s] for s in shared])

    # 组顺序
    if p["group_order"]:
        group_order = [g for g in p["group_order"] if g in set(groups_per_sample)]
    else:
        group_order = list(dict.fromkeys(groups_per_sample.tolist()))
    if len(group_order) < 2:
        raise ValueError("LEfSe 至少需要 2 个组")

    # === 每个 taxon：KW + 组均值 + LDA 效应量 ===
    rows = []
    X = ab.to_numpy()  # (n_taxa, n_samples)
    for i, feature in enumerate(ab.index):
        values = X[i]
        if values.sum() == 0:
            continue
        by_group = [values[groups_per_sample == g] for g in group_order]
        if any(len(v) < 2 for v in by_group):
            continue
        try:
            kw_stat, kw_p = sp_stats.kruskal(*by_group)
        except ValueError:
            continue
        means = np.array([v.mean() for v in by_group])
        enriched_idx = int(np.argmax(means))
        enriched = group_order[enriched_idx]
        diff = float(means.max() - means.min())
        lda = float(np.log10(1.0 + p["lda_scale"] * abs(diff))) if diff > 0 else 0.0
        rows.append({
            "feature": str(feature),
            "short_name": _short_name(feature),
            "tax_level": _tax_level_from_name(feature),
            "group": enriched,
            "lda": lda,
            "kw_stat": float(kw_stat),
            "kw_p": float(kw_p),
            **{f"mean_{g}": float(means[j]) for j, g in enumerate(group_order)},
        })

    stats_df = pd.DataFrame(rows)
    if stats_df.empty:
        raise ValueError("未得到任何可检验特征（检查数据非空 + 每组 ≥ 2 个样本）")

    # KW 筛选
    sig = stats_df[stats_df["kw_p"] < p["alpha_kw"]].copy()
    # LDA 阈值
    sig = sig[sig["lda"] >= p["lda_threshold"]].copy()
    # tax 层级过滤
    if p["tax_levels"]:
        sig = sig[sig["tax_level"].isin(p["tax_levels"])].copy()
    if sig.empty:
        raise ValueError(
            f"无特征通过筛选（KW p<{p['alpha_kw']} + LDA≥{p['lda_threshold']}"
            + (f" + levels={p['tax_levels']}" if p["tax_levels"] else "") + "）"
        )

    # 组内降序 + 组间按 group_order
    sig["group"] = pd.Categorical(sig["group"], categories=group_order, ordered=True)
    sig = sig.sort_values(["group", "lda"], ascending=[True, False]).reset_index(drop=True)
    if p["max_features"] and len(sig) > p["max_features"]:
        # 每组按比例截断
        per = max(1, p["max_features"] // len(group_order))
        sig = sig.groupby("group", observed=True).head(per).reset_index(drop=True)

    # === 绘图 ===
    fig, ax = plt.subplots(
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4),
        constrained_layout=True,
    )
    # 由下至上：组内从小到大（横向条图 invert）
    plot_df = sig.iloc[::-1].reset_index(drop=True)
    palette = {**DEFAULT_PALETTE, **(p["palette"] or {})}
    cmap_fallback = plt.cm.tab10.colors
    for i, g in enumerate(group_order):
        palette.setdefault(g, cmap_fallback[i % len(cmap_fallback)])
    colors = [palette[str(g)] for g in plot_df["group"]]
    y_pos = np.arange(len(plot_df))
    ax.barh(y_pos, plot_df["lda"], color=colors, edgecolor="white",
            height=0.7, linewidth=0.4)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df["short_name"], fontsize=7)
    ax.set_xlabel("LDA Score (log 10)", fontsize=9)
    ax.set_ylabel("")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.axvline(p["lda_threshold"], color="#999", lw=0.5, ls="--")

    # 图例
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=palette[g], label=g) for g in group_order]
    ax.legend(handles=handles, title="Enriched in", loc="lower right",
              frameon=False, fontsize=8, title_fontsize=9)

    return AnalysisResult(figure=fig, stats=sig, params=p)
