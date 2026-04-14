"""RDA 约束排序图（环境因子约束物种）。

对标 scripts/R/03_RDA.R：Hellinger 变换 + skbio RDA + Mantel 检验。

输入：
    abundance_df: 丰度表（首列 Taxonomy，其余列为样本）
    env_df:       环境因子表（含 SampleID + Group + 数值因子列）
    metadata_df:  SampleID + Group 分组（用于颜色与样本 ID 对齐）

输出（AnalysisResult）：
    figure — 样本散点（按组色）+ 环境因子箭头（带 Mantel 显著性）
    stats  — 扁平表：variable, type, value（含 explained_variance + mantel_r/p + 箭头坐标）
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skbio.stats.distance import DistanceMatrix, mantel
from skbio.stats.ordination import rda as skbio_rda
from sklearn.preprocessing import StandardScaler

from envmeta.analysis.base import AnalysisResult

UNCLASSIFIED_REGEX = re.compile(r"(?i)(unclassif|unassign|unknown)")

DEFAULTS = {
    "width_mm": 180,
    "height_mm": 140,
    "n_permutations": 999,
    "show_env_arrows": True,
    "show_sample_labels": True,
    "arrow_scale": 1.0,
    "group_col": "Group",
    "palette": None,
    "drop_unclassified": True,
    "random_seed": 123,
}

DEFAULT_PALETTE = {"CK": "#1c9cbd", "A": "#e3943d", "B": "#92181e"}


def _sig_label(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""


def _hellinger(mat: np.ndarray) -> np.ndarray:
    """Hellinger 变换：sqrt(x_ij / row_sum_i)。mat 行=样本，列=物种。"""
    row_sums = mat.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    return np.sqrt(mat / row_sums)


def _align_env_to_abundance(
    ab_cols: list[str], env_df: pd.DataFrame, metadata_df: pd.DataFrame,
    group_col: str,
) -> tuple[list[str], pd.DataFrame]:
    """返回 (abundance 列的子集, 对应顺序的 env DataFrame)。

    优先按 SampleID 直接交集；否则按 (Group, 出现顺序) 匹配（兼容 env 用
    别名 SampleID 的情况，如 abundance=2_1 而 env=CK_1）。
    """
    env_sample_col = "SampleID" if "SampleID" in env_df.columns else "Sample_ID"
    meta_sample_col = "SampleID" if "SampleID" in metadata_df.columns else "Sample_ID"
    env_ids = env_df[env_sample_col].astype(str).tolist()

    # Case 1：直接交集
    common = [c for c in ab_cols if c in set(env_ids)]
    if len(common) >= 3:
        ab_cols_kept = common
        env_aligned = env_df.set_index(env_sample_col).loc[common].reset_index()
        return ab_cols_kept, env_aligned

    # Case 2：按 Group + 出现顺序匹配
    meta = metadata_df.copy()
    meta[meta_sample_col] = meta[meta_sample_col].astype(str)
    # 丰度表列 → 组（按 metadata）
    ab_to_group = dict(zip(meta[meta_sample_col], meta[group_col].astype(str)))
    # env 按 group 出现顺序
    env_per_group: dict[str, list[int]] = {}
    for i, g in enumerate(env_df[group_col].astype(str)):
        env_per_group.setdefault(g, []).append(i)

    # 丰度表按 group 出现顺序过滤
    ab_cols_kept = []
    env_rows_idx = []
    ab_per_group_seen: dict[str, int] = {}
    for c in ab_cols:
        g = ab_to_group.get(c)
        if g is None:
            continue
        pos = ab_per_group_seen.get(g, 0)
        if g in env_per_group and pos < len(env_per_group[g]):
            ab_cols_kept.append(c)
            env_rows_idx.append(env_per_group[g][pos])
            ab_per_group_seen[g] = pos + 1
    if len(ab_cols_kept) < 3:
        raise ValueError(
            f"丰度表列和 env 样本无法对齐（直接交集={len(common)}, "
            f"按组对齐={len(ab_cols_kept)}），至少需要 3 个样本"
        )
    env_aligned = env_df.iloc[env_rows_idx].reset_index(drop=True)
    return ab_cols_kept, env_aligned


def analyze(
    abundance_df: pd.DataFrame,
    env_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    params: dict | None = None,
) -> AnalysisResult:
    p = {**DEFAULTS, **(params or {})}
    np.random.seed(p["random_seed"])

    # 丰度：行=样本，列=物种（Hellinger 之前转置）
    ab = abundance_df.copy()
    ab = ab.set_index(ab.columns[0])
    if p["drop_unclassified"]:
        mask = ab.index.to_series().apply(lambda s: not UNCLASSIFIED_REGEX.search(str(s)))
        ab = ab.loc[mask]
    ab = ab.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # env 对齐到丰度表样本
    kept_cols, env_aligned = _align_env_to_abundance(
        ab.columns.tolist(), env_df, metadata_df, p["group_col"]
    )
    ab = ab[kept_cols]                          # 物种 × 样本
    species_by_sample = ab.T.to_numpy()         # 样本 × 物种

    env_num_cols = [
        c for c in env_aligned.columns
        if c not in (p["group_col"], "SampleID", "Sample_ID", "Replicate")
        and pd.to_numeric(env_aligned[c], errors="coerce").notna().mean() > 0.95
    ]
    if len(env_num_cols) < 2:
        raise ValueError(f"env 表至少需要 2 个数值因子列（当前 {env_num_cols}）")
    env_mat = env_aligned[env_num_cols].astype(float).to_numpy()

    # 变换 + 标准化
    Y = _hellinger(species_by_sample)
    X = StandardScaler().fit_transform(env_mat)

    # skbio RDA
    ord_res = skbio_rda(Y, X, sample_ids=kept_cols, feature_ids=ab.index.tolist(),
                        constraint_ids=env_num_cols)

    # 坐标
    site = ord_res.samples.iloc[:, :2].to_numpy()  # 样本得分（RDA1, RDA2）
    biplot = ord_res.biplot_scores.iloc[:, :2].to_numpy()  # 环境因子箭头
    explained = ord_res.proportion_explained.iloc[:2].to_numpy()

    # Bray-Curtis 距离 + Mantel（按环境因子）
    from scipy.spatial.distance import pdist, squareform
    bc = squareform(pdist(species_by_sample, metric="braycurtis"))
    dm_species = DistanceMatrix(bc, ids=kept_cols)

    mantel_rows = []
    for i, factor in enumerate(env_num_cols):
        env_vec = env_mat[:, i]
        env_dist = squareform(pdist(env_vec.reshape(-1, 1), metric="euclidean"))
        dm_env = DistanceMatrix(env_dist, ids=kept_cols)
        try:
            r, pv, _ = mantel(dm_species, dm_env, method="pearson",
                              permutations=p["n_permutations"],
                              seed=p["random_seed"])
        except Exception:
            r, pv = np.nan, 1.0
        mantel_rows.append({
            "variable": factor, "type": "mantel",
            "mantel_r": float(r), "mantel_p": float(pv),
            "significance": _sig_label(pv),
            "arrow_x": float(biplot[i, 0]),
            "arrow_y": float(biplot[i, 1]),
        })

    # stats：解释度 + Mantel
    stats_rows = [
        {"variable": "RDA1", "type": "explained_variance",
         "value": float(explained[0])},
        {"variable": "RDA2", "type": "explained_variance",
         "value": float(explained[1])},
    ]
    stats_rows.extend(mantel_rows)
    stats_df = pd.DataFrame(stats_rows)

    # === 绘图 ===
    fig, ax = plt.subplots(
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4),
        constrained_layout=True,
    )

    # 样本点（按组颜色）
    meta_idx = metadata_df.set_index(
        "SampleID" if "SampleID" in metadata_df.columns else "Sample_ID"
    )
    groups = [str(meta_idx.loc[s, p["group_col"]]) if s in meta_idx.index else "?"
              for s in kept_cols]
    palette = {**DEFAULT_PALETTE, **(p["palette"] or {})}
    group_order = list(dict.fromkeys(groups))
    cmap_fallback = plt.cm.tab10.colors
    for i, g in enumerate(group_order):
        palette.setdefault(g, cmap_fallback[i % len(cmap_fallback)])

    for g in group_order:
        idx = [i for i, gg in enumerate(groups) if gg == g]
        ax.scatter(site[idx, 0], site[idx, 1], s=80,
                   color=palette[g], label=g, edgecolor="white",
                   linewidth=0.8, zorder=3)

    # 环境因子箭头
    if p["show_env_arrows"]:
        scale = p["arrow_scale"] * max(np.abs(site).max() / max(np.abs(biplot).max(), 1e-9), 1.0)
        for i, factor in enumerate(env_num_cols):
            x, y = biplot[i, 0] * scale, biplot[i, 1] * scale
            ax.annotate(
                "", xy=(x, y), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color="#555", lw=1.2),
            )
            sig = mantel_rows[i]["significance"]
            label = f"{factor}{sig}" if sig else factor
            ax.text(x * 1.08, y * 1.08, label, fontsize=8, color="#333",
                    ha="left" if x >= 0 else "right",
                    va="bottom" if y >= 0 else "top")

    if p["show_sample_labels"]:
        for i, s in enumerate(kept_cols):
            ax.text(site[i, 0], site[i, 1], f" {s}", fontsize=7,
                    color="#333", va="center")

    ax.axhline(0, color="#999", lw=0.5, ls="--")
    ax.axvline(0, color="#999", lw=0.5, ls="--")
    ax.set_xlabel(f"RDA1 ({explained[0]*100:.1f}%)")
    ax.set_ylabel(f"RDA2 ({explained[1]*100:.1f}%)")
    ax.legend(loc="best", frameon=False, fontsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    return AnalysisResult(figure=fig, stats=stats_df, params=p)
