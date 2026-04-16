"""Top-N MAG 丰度热图（对标 scripts/python/07_MAG_abundance_heatmap.py）。

核心：在长尾分布的 MAG × sample 丰度矩阵中突出 Top-N，提供门分组 + 关键物种
标注 + 组顶彩条 + 三段非线性配色（适配少数极高 + 多数低值的场景）。

输入：
    abundance_df: MAG × sample 宽表（首列 MAG / Genome / Name；其余为样本丰度，%）
    taxonomy_df:  可选 MAG + GTDB classification（用于门彩条、门内聚类）
    keystone_df:  可选 MAG 列表（MAG + 可选 Genus），打 ★ 高亮
    metadata_df:  可选 SampleID + Group（用于组彩条 + 样本按组排序）

输出：
    figure — 主热图 + 门色带 + 组色带 + 关键物种 ★ + 色标 + 图例
    stats  — 合并表：MAG × sample 丰度百分比 + Phylum + is_keystone +
             row_order（聚类后行序）+ selection_score（mean/sum/var 指标）
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from envmeta.analysis.base import AnalysisResult
from envmeta.analysis.pathway import PHYLUM_COLORS, _extract_phylum, _mag_col

GROUP_COLORS = {"CK": "#4DAF4A", "A": "#377EB8", "B": "#E41A1C"}

DEFAULTS = {
    "top_n": 30,
    "selection_by": "mean",               # "mean" | "sum" | "variance"
    "log_transform": True,                # 仅用于聚类距离，不改变色板映射
    "cluster_rows": True,
    "cluster_cols": False,
    "cluster_within_phylum": True,        # True → 先按门分组，门内聚类
    "linkage_method": "average",          # "average" | "ward" | "complete"
    "color_breakpoints": (0.2, 0.5),      # 三段非线性配色边界（%）
    "show_phylum_bar": True,
    "show_group_bar": True,
    "highlight_keystones": True,
    "width_mm": 180,
    "height_mm": 220,
}


def _select_top_n(mat: np.ndarray, top_n: int, by: str) -> np.ndarray:
    """返回按 mean/sum/variance 排序的行索引（降序，截取 Top-N）。"""
    if by == "sum":
        score = mat.sum(axis=1)
    elif by == "variance":
        score = mat.var(axis=1)
    else:
        score = mat.mean(axis=1)
    n = min(top_n, len(score))
    order = np.argsort(-score, kind="stable")[:n]
    return order, score


def _build_tricolor_cmap(mat: np.ndarray, breakpoints: tuple[float, float]):
    """三段非线性配色：低段 Blues（频繁低值） → 中段 YlGn → 高段 YlOrRd。

    目的：把视觉分辨率集中在"少数极高值"与"低值背景之上的中等值"两段上，
    而非被一两个极值压缩到全图蓝灰。
    """
    vmax = np.ceil(mat.max() * 10) / 10 if mat.size else 1.0
    lo_bp, hi_bp = breakpoints
    if vmax < hi_bp + 0.1:
        vmax = hi_bp + 0.1

    # 低段细分，高段粗分
    bounds_lo = list(np.arange(0, lo_bp + 1e-9, 0.025))
    bounds_mid = list(np.arange(lo_bp, hi_bp + 1e-9, 0.05))
    bounds_hi = list(np.arange(hi_bp, min(vmax, 1.0) + 1e-9, 0.1))
    if vmax > 1.0:
        bounds_hi += list(np.arange(1.0, vmax + 0.01, 0.5))
    bounds = sorted({round(b, 3) for b in bounds_lo + bounds_mid + bounds_hi})
    if bounds[-1] < vmax:
        bounds.append(round(float(vmax), 2))

    n_lo = sum(1 for b in bounds if b <= lo_bp) - 1
    n_mid = sum(1 for b in bounds if lo_bp < b <= hi_bp)
    n_hi = sum(1 for b in bounds if b > hi_bp)

    seg_lo = plt.cm.Blues(np.linspace(0.08, 0.55, max(n_lo + 1, 2)))
    seg_mid = plt.cm.YlGn(np.linspace(0.18, 0.65, max(n_mid + 1, 2)))
    seg_hi = plt.cm.YlOrRd(np.linspace(0.40, 0.95, max(n_hi + 1, 2)))
    all_colors = (list(seg_lo[:n_lo + 1])
                  + list(seg_mid[1:n_mid + 1])
                  + list(seg_hi[1:n_hi + 1]))
    cmap = LinearSegmentedColormap.from_list("tricolor", all_colors, N=256)
    norm = BoundaryNorm(bounds, ncolors=256)
    return cmap, norm, bounds, vmax


def _cluster_order(mat: np.ndarray, method: str) -> list[int]:
    """返回层次聚类的叶子顺序；≤2 行时返回原序。"""
    if len(mat) <= 2:
        return list(range(len(mat)))
    try:
        dist = pdist(mat, metric="euclidean")
        link = linkage(dist, method=method)
        return list(leaves_list(link))
    except Exception:
        return list(range(len(mat)))


def analyze(
    abundance_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    metadata_df: pd.DataFrame | None = None,
    params: dict | None = None,
) -> AnalysisResult:
    p = {**DEFAULTS, **(params or {})}

    # ── 丰度矩阵 ─────────────────────────────────────────────
    ab = abundance_df.copy()
    mag_c = _mag_col(ab)
    ab = ab.rename(columns={mag_c: "MAG"})
    ab = ab[ab["MAG"].astype(str) != "unmapped"]
    sample_cols = [c for c in ab.columns if c != "MAG"]
    if not sample_cols:
        raise ValueError("abundance 表缺少样本列")
    for c in sample_cols:
        ab[c] = pd.to_numeric(ab[c], errors="coerce")
    ab = ab.dropna(subset=sample_cols, how="all").fillna(0.0)
    ab["MAG"] = ab["MAG"].astype(str)
    ab = ab.drop_duplicates(subset="MAG").reset_index(drop=True)
    if ab.empty:
        raise ValueError("abundance 表无有效 MAG")

    mat_full = ab[sample_cols].to_numpy(dtype=float)
    row_order, scores = _select_top_n(mat_full, p["top_n"], p["selection_by"])
    df = ab.iloc[row_order].copy().reset_index(drop=True)
    df["selection_score"] = scores[row_order]
    mat = df[sample_cols].to_numpy(dtype=float)

    # ── 分类 ────────────────────────────────────────────────
    if taxonomy_df is not None and not taxonomy_df.empty:
        tax = taxonomy_df.copy()
        tax = tax.rename(columns={_mag_col(tax): "MAG"})
        tax["MAG"] = tax["MAG"].astype(str)
        cls_col = next(
            (c for c in tax.columns
             if "classif" in c.lower() or "taxonomy" in c.lower()),
            tax.columns[1] if len(tax.columns) > 1 else None,
        )
        if cls_col is not None:
            tax["Phylum"] = tax[cls_col].apply(_extract_phylum)
            df = df.merge(tax[["MAG", "Phylum"]], on="MAG", how="left")
    df["Phylum"] = df.get("Phylum", pd.Series(["Unknown"] * len(df))).fillna("Unknown")

    # ── keystone ────────────────────────────────────────────
    if keystone_df is not None and not keystone_df.empty:
        ks = keystone_df.copy()
        ks = ks.rename(columns={_mag_col(ks): "MAG"})
        df["is_keystone"] = df["MAG"].isin(set(ks["MAG"].astype(str)))
    else:
        df["is_keystone"] = False

    # ── metadata → group / sample 排序 ───────────────────────
    groups_per_sample: dict[str, str | None] = {s: None for s in sample_cols}
    if metadata_df is not None and not metadata_df.empty:
        md = metadata_df.copy()
        sid_col = next((c for c in md.columns if c.lower() in ("sampleid", "sample", "sample_id")), md.columns[0])
        grp_col = next((c for c in md.columns if c.lower() == "group"), None)
        if grp_col is not None:
            groups_per_sample = dict(zip(md[sid_col].astype(str), md[grp_col].astype(str)))
            # 按组分组后保持原样本在组内的相对顺序
            by_grp: dict[str, list[str]] = {}
            for s in sample_cols:
                g = groups_per_sample.get(s)
                by_grp.setdefault(g if g is not None else "", []).append(s)
            new_order: list[str] = []
            known = [g for g in ("CK", "A", "B") if g in by_grp]
            extra = [g for g in by_grp if g not in known]
            for g in known + extra:
                new_order.extend(by_grp[g])
            if new_order and new_order != sample_cols:
                sample_cols = new_order
                df = df[["MAG"] + sample_cols + [c for c in df.columns
                                                   if c not in {"MAG", *sample_cols}]]
                mat = df[sample_cols].to_numpy(dtype=float)

    # ── 行聚类（门内 or 全局）──────────────────────────────
    if p["cluster_rows"]:
        mat_for_clust = np.log1p(mat) if p["log_transform"] else mat
        if p["cluster_within_phylum"] and "Phylum" in df.columns:
            phy_order = df["Phylum"].value_counts().index.tolist()
            reordered: list[int] = []
            for phy in phy_order:
                idx = df.index[df["Phylum"] == phy].tolist()
                if not idx:
                    continue
                sub = mat_for_clust[idx]
                leaf = _cluster_order(sub, p["linkage_method"])
                reordered.extend([idx[k] for k in leaf])
            if reordered:
                df = df.iloc[reordered].reset_index(drop=True)
                mat = df[sample_cols].to_numpy(dtype=float)
        else:
            leaf = _cluster_order(mat_for_clust, p["linkage_method"])
            df = df.iloc[leaf].reset_index(drop=True)
            mat = df[sample_cols].to_numpy(dtype=float)
    df["row_order"] = range(len(df))

    # ── 列聚类（可选）──────────────────────────────────────
    if p["cluster_cols"]:
        mat_T = np.log1p(mat.T) if p["log_transform"] else mat.T
        col_leaf = _cluster_order(mat_T, p["linkage_method"])
        sample_cols = [sample_cols[k] for k in col_leaf]
        df = df[["MAG"] + sample_cols + [c for c in df.columns
                                           if c not in {"MAG", *sample_cols}]]
        mat = df[sample_cols].to_numpy(dtype=float)

    fig = _draw(df, sample_cols, mat, groups_per_sample, p)

    stats_df = df[["MAG", "Phylum", "is_keystone", "row_order",
                   "selection_score"] + sample_cols].copy()
    return AnalysisResult(figure=fig, stats=stats_df, params=p)


def _draw(df, sample_cols, mat, groups_per_sample, p) -> plt.Figure:
    n_mag, n_sample = mat.shape
    cmap, norm, bounds, vmax = _build_tricolor_cmap(mat, tuple(p["color_breakpoints"]))

    fig, ax = plt.subplots(
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4),
        constrained_layout=True,
    )
    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm,
                   interpolation="nearest")

    # 顶部组色带
    if p["show_group_bar"] and any(groups_per_sample.values()):
        for j, s in enumerate(sample_cols):
            g = groups_per_sample.get(s)
            color = GROUP_COLORS.get(g, "#CCCCCC")
            ax.add_patch(Rectangle(
                (j - 0.5, -1.6), 1.0, 0.9,
                facecolor=color, edgecolor="none",
                clip_on=False, alpha=0.9, zorder=3,
            ))
        # 组间竖线分隔
        prev = None
        for j, s in enumerate(sample_cols):
            g = groups_per_sample.get(s)
            if prev is not None and g != prev:
                ax.axvline(j - 0.5, color="black", lw=1.2, zorder=4)
            prev = g

    # 左侧门色带
    if p["show_phylum_bar"]:
        for i, phy in enumerate(df["Phylum"].tolist()):
            ax.add_patch(Rectangle(
                (-1.5, i - 0.5), 0.8, 1.0,
                facecolor=PHYLUM_COLORS.get(phy, "#888"),
                edgecolor="none", clip_on=False, zorder=3,
            ))

    # 行标签（MAG + ★）
    labels = []
    for _, r in df.iterrows():
        prefix = "★ " if (p["highlight_keystones"] and r["is_keystone"]) else ""
        labels.append(prefix + r["MAG"])
    ax.set_yticks(range(n_mag))
    ax.set_yticklabels(labels, fontsize=max(4, min(8, 260 // max(n_mag, 1))))

    # 列标签
    ax.set_xticks(range(n_sample))
    ax.set_xticklabels(sample_cols, rotation=45, ha="right",
                       fontsize=max(5, min(9, 280 // max(n_sample, 1))))

    # 门边界横线
    prev_phy = None
    for i, phy in enumerate(df["Phylum"].tolist()):
        if prev_phy is not None and phy != prev_phy:
            ax.axhline(i - 0.5, color="#666", lw=0.5, zorder=3)
        prev_phy = phy

    ax.set_title(
        f"Top {n_mag} MAG Abundance  ({n_sample} samples,"
        f" color = relative abundance %)",
        fontsize=10, fontweight="bold", pad=12,
    )
    cbar = fig.colorbar(im, ax=ax, shrink=0.5,
                        label="Relative abundance (%)")
    cbar.ax.tick_params(labelsize=7)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    return fig
