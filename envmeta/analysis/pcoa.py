"""PCoA (Principal Coordinates Analysis) + PERMANOVA。

对标 scripts/R/02_beta_PCoA.R。输入预计算距离矩阵（如 Bray-Curtis），
不自己算 distance。
"""
from __future__ import annotations

from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa

from envmeta.analysis.base import AnalysisResult

# ==============================================================================
# 常量
# ==============================================================================

DEFAULT_PALETTE = {
    "CK": "#1c9cbd", "A": "#e3943d", "B": "#92181e",
    # 兜底：后续组别从下面 extra 里循环取
}
_EXTRA_COLORS = ["#7a3b69", "#2f5d50", "#8a6d3b", "#5b3a29", "#36454F"]

DEFAULTS = {
    "width_mm": 120,
    "height_mm": 100,
    "marker_size": 90,
    "axes": (0, 1),              # PC1 vs PC2
    "show_labels": True,
    "show_permanova": True,
    "n_permutations": 999,
    "group_col": "Group",
    "palette": None,             # None → 用 DEFAULT_PALETTE + EXTRA 兜底
    "random_seed": 123,
    "title_size": 11,
    "label_size": 9,
    "tick_size": 8,
}


# ==============================================================================
# 核心
# ==============================================================================

def _align_samples(dist_df: pd.DataFrame, meta_df: pd.DataFrame, group_col: str) -> tuple[list[str], pd.Series]:
    """返回 (共同样本列表, 样本→组别 Series)。"""
    meta = meta_df.copy()
    meta[meta.columns[0]] = meta[meta.columns[0]].astype(str)
    meta = meta.set_index(meta.columns[0])       # SampleID / Sample_ID

    # 距离矩阵的行/列应一致。我们从"除首列外"的样本名 × 首列值取交集
    dist_samples = [str(s) for s in dist_df.columns[1:]]
    meta_samples = meta.index.astype(str).tolist()
    common = [s for s in meta_samples if s in dist_samples]
    if not common:
        raise ValueError("距离矩阵和 metadata 没有共同样本")

    groups = meta.loc[common, group_col].astype(str)
    return common, groups


def _subset_distance(dist_df: pd.DataFrame, samples: list[str]) -> np.ndarray:
    """从 pandas DataFrame（首列是样本名）取出方阵；确保 C-contiguous（skbio 要求）。"""
    idx = dist_df.set_index(dist_df.columns[0])
    mat = idx.loc[samples, samples].apply(pd.to_numeric, errors="coerce").to_numpy()
    return np.ascontiguousarray(mat, dtype=np.float64)


def _resolve_palette(user_palette: dict | None, groups_order: list[str]) -> dict[str, str]:
    pal = dict(DEFAULT_PALETTE)
    if user_palette:
        pal.update(user_palette)
    extras = iter(_EXTRA_COLORS)
    for g in groups_order:
        if g not in pal:
            pal[g] = next(extras, "#888888")
    return pal


def _pairwise_permanova(mat: np.ndarray, samples: list[str], groups: pd.Series,
                        n_perm: int, seed: int) -> pd.DataFrame:
    rows = []
    unique_groups = list(dict.fromkeys(groups.tolist()))
    for g1, g2 in combinations(unique_groups, 2):
        mask = groups.isin([g1, g2]).to_numpy()
        if mask.sum() < 3:
            continue
        sub_samples = [s for s, m in zip(samples, mask) if m]
        sub_groups = groups[mask].tolist()
        sub_mat = np.ascontiguousarray(mat[np.ix_(mask, mask)], dtype=np.float64)
        np.random.seed(seed)
        dm = DistanceMatrix(sub_mat, ids=sub_samples)
        res = permanova(dm, grouping=sub_groups, permutations=n_perm)
        rows.append({"group1": g1, "group2": g2,
                     "pseudo_F": float(res["test statistic"]),
                     "p_value": float(res["p-value"])})
    return pd.DataFrame(rows)


def _build_figure(coords: pd.DataFrame, groups: pd.Series, explained: pd.Series,
                  axes: tuple[int, int], palette: dict[str, str], *,
                  marker_size: int, show_labels: bool,
                  permanova_text: str | None,
                  width_mm: int, height_mm: int,
                  title_size: int, label_size: int, tick_size: int) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4), constrained_layout=True)

    i, j = axes
    x = coords.iloc[:, i]
    y = coords.iloc[:, j]
    for g in dict.fromkeys(groups.tolist()):
        mask = (groups == g).to_numpy()
        ax.scatter(x[mask], y[mask], s=marker_size, c=palette[g],
                   alpha=0.9, edgecolors="white", linewidths=0.8, label=g, zorder=3)

    # 轴标签：PC1 (variance%)
    ax.set_xlabel(f"PC{i+1} ({explained.iloc[i]*100:.1f}%)", fontsize=label_size)
    ax.set_ylabel(f"PC{j+1} ({explained.iloc[j]*100:.1f}%)", fontsize=label_size)
    ax.tick_params(labelsize=tick_size)
    ax.axhline(0, color="#aaaaaa", linewidth=0.5, linestyle="--", zorder=1)
    ax.axvline(0, color="#aaaaaa", linewidth=0.5, linestyle="--", zorder=1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    if show_labels:
        try:
            from adjustText import adjust_text
            texts = [ax.text(xi, yi, lbl, fontsize=tick_size, color="#444")
                     for xi, yi, lbl in zip(x, y, coords.index)]
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#888", lw=0.3))
        except ImportError:
            for xi, yi, lbl in zip(x, y, coords.index):
                ax.annotate(lbl, (xi, yi), fontsize=tick_size, color="#444",
                            xytext=(3, 3), textcoords="offset points")

    if permanova_text:
        ax.text(0.98, 0.98, permanova_text, transform=ax.transAxes,
                ha="right", va="top", fontsize=tick_size,
                style="italic", color="#555",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8))

    ax.legend(loc="best", frameon=False, fontsize=label_size, title_fontsize=label_size)

    return fig


# ==============================================================================
# 公共入口
# ==============================================================================

def analyze(distance_df: pd.DataFrame, metadata_df: pd.DataFrame,
            params: dict | None = None) -> AnalysisResult:
    p = {**DEFAULTS, **(params or {})}

    common, groups = _align_samples(distance_df, metadata_df, p["group_col"])
    mat = _subset_distance(distance_df, common)

    np.random.seed(p["random_seed"])
    dm = DistanceMatrix(mat, ids=common)
    ord_result = pcoa(dm)

    coords = ord_result.samples                # DataFrame (n_samples × n_axes)
    explained = ord_result.proportion_explained  # Series (n_axes,)

    # 全局 PERMANOVA
    perm_res = permanova(dm, grouping=groups.tolist(), permutations=p["n_permutations"])
    pseudo_F = float(perm_res["test statistic"])
    p_val = float(perm_res["p-value"])
    # R² 近似：R² ≈ (pseudo_F × (N-k)) / (pseudo_F × (N-k) + k-1)
    # 其中 k=组数, N=样本数（此处给出一个"总变差解释"的近似值）
    n = len(common)
    k = len(set(groups.tolist()))
    r2 = (pseudo_F * (k - 1)) / (pseudo_F * (k - 1) + (n - k)) if (n > k) else float("nan")

    pairwise = _pairwise_permanova(mat, common, groups, p["n_permutations"], p["random_seed"])

    palette = _resolve_palette(p["palette"], list(dict.fromkeys(groups.tolist())))
    perm_text = f"PERMANOVA  R²={r2:.3f}, P={p_val:.3f}" if p["show_permanova"] else None

    fig = _build_figure(
        coords, groups, explained, axes=p["axes"], palette=palette,
        marker_size=p["marker_size"], show_labels=p["show_labels"],
        permanova_text=perm_text,
        width_mm=p["width_mm"], height_mm=p["height_mm"],
        title_size=p["title_size"], label_size=p["label_size"], tick_size=p["tick_size"],
    )

    # 汇总 stats：explained + PERMANOVA + pairwise 三张表拼在一个 DataFrame 里
    summary_rows = [
        {"metric": "PERMANOVA_pseudo_F", "value": pseudo_F},
        {"metric": "PERMANOVA_R2_approx", "value": r2},
        {"metric": "PERMANOVA_p_value", "value": p_val},
        {"metric": "n_permutations", "value": p["n_permutations"]},
        {"metric": "n_samples", "value": n},
        {"metric": "n_groups", "value": k},
    ]
    for i, (ax_name, val) in enumerate(explained.head(5).items()):
        summary_rows.append({"metric": f"PC{i+1}_explained", "value": float(val)})
    summary = pd.DataFrame(summary_rows)
    stats_dict = {"summary": summary, "pairwise_permanova": pairwise, "coords": coords}

    return AnalysisResult(figure=fig, stats=summary, params={**p, "_stats_tables": stats_dict})
