"""4 张 MAG-based 分析图共享的辅助函数与常量。

目的：让 mag_quality / mag_heatmap / pathway / gene_profile 四张图的
- MAG 标签规则（Genus species / Genus sp. Mx_XX / MAG_id）
- Phylum 颜色映射 + 彩条 + 右侧图例
- MAG 子集过滤（filter_mode：all / top_n / keystone_only / top_plus_keystone）
- 行排序（phylum_cluster / metric_desc / abundance）
完全一致，跨图对照时视觉元素（物种名格式、门色块、keystone 标记）稳定。
"""
from __future__ import annotations

from typing import Literal

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

# ── 常量 ─────────────────────────────────────────────────────

PHYLUM_COLORS: dict[str, str] = {
    "Pseudomonadota":      "#E74C3C",
    "Acidobacteriota":     "#3498DB",
    "Chloroflexota":       "#2ECC71",
    "Bacteroidota_A":      "#F39C12",
    "Patescibacteriota":   "#9B59B6",
    "Desulfobacterota":    "#1ABC9C",
    "Actinomycetota":      "#E67E22",
    "Gemmatimonadota":     "#34495E",
    "Planctomycetota":     "#D35400",
    "Myxococcota":         "#8E44AD",
    "Nitrospirota":        "#27AE60",
    "Zixibacteria":        "#C0392B",
    "Desulfobacterota_B":  "#E78AC3",
    "Desulfobacterota_E":  "#FC8D59",
    "Verrucomicrobiota":   "#5AAE61",
    "Methylomirabilota":   "#878787",
    "Bacteroidota":        "#FEE0B6",
    "Thermoproteota":      "#B35806",
    "Halobacteriota":      "#8073AC",
    "Other":               "#BDC3C7",
    "Unknown":             "#888888",
}

GROUP_COLORS: dict[str, str] = {
    "CK": "#4DAF4A",
    "A":  "#377EB8",
    "B":  "#E41A1C",
}

FilterMode = Literal["all", "top_n", "keystone_only", "top_plus_keystone"]
TopNBy = Literal["mean", "sum", "variance"]
RowOrder = Literal["phylum_cluster", "metric_desc", "abundance"]


# ── 列名嗅探 + GTDB 解析 ────────────────────────────────────

def mag_col(df: pd.DataFrame) -> str:
    """在 MAG 相关表里嗅探 MAG 列名（宽松匹配）。"""
    for c in ("MAG", "Name", "user_genome", "Genome", "Bin"):
        if c in df.columns:
            return c
    return df.columns[0]


def extract_rank(cls: str, prefix: str) -> str:
    """GTDB classification 字符串里取某级（例 'g__' → Genus）。"""
    if not isinstance(cls, str):
        return ""
    for part in cls.split(";"):
        part = part.strip()
        if part.startswith(prefix):
            return part[len(prefix):].strip()
    return ""


def extract_phylum(cls: str) -> str:
    phy = extract_rank(cls, "p__")
    return phy or "Unknown"


def mag_display_label(mag: str, genus: str = "", species: str = "",
                      family: str = "") -> str:
    """四档 fallback（S2.5-13 规则扩展）：
    - Genus species         （两者都有）
    - Genus sp. Mx_XX       （只有 Genus）
    - Family sp. Mx_XX      （只有 Family；比如 GTDB `g__;` 但 `f__` 有值）
    - Mx_XX                 （Genus / Family 都空）
    """
    g = (genus or "").strip()
    s = (species or "").strip()
    f = (family or "").strip()
    tail = str(mag).split("_")[-1] if "_" in str(mag) else str(mag)
    if g and s:
        if s.startswith(g + " "):
            s = s[len(g) + 1:]
        return f"{g} {s}".strip()
    if g:
        return f"{g} sp. Mx_{tail}"
    if f:
        return f"{f} sp. Mx_{tail}"
    return str(mag)


def annotate_taxonomy(
    df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None,
) -> pd.DataFrame:
    """把 taxonomy classification 解析成 Phylum/Genus/Species/label 四列加进 df。

    df 必须有 "MAG" 列；若传入 taxonomy_df=None，Phylum="Unknown"，
    label=MAG id。
    """
    out = df.copy()
    out["MAG"] = out["MAG"].astype(str)
    if taxonomy_df is not None and not taxonomy_df.empty:
        tax = taxonomy_df.copy()
        tax = tax.rename(columns={mag_col(tax): "MAG"})
        tax["MAG"] = tax["MAG"].astype(str)
        cls_col = next(
            (c for c in tax.columns
             if "classif" in c.lower() or "taxonomy" in c.lower()),
            tax.columns[1] if len(tax.columns) > 1 else None,
        )
        if cls_col is not None:
            tax["Phylum"] = tax[cls_col].apply(extract_phylum)
            tax["Family"] = tax[cls_col].apply(lambda x: extract_rank(x, "f__"))
            tax["Genus"] = tax[cls_col].apply(lambda x: extract_rank(x, "g__"))
            tax["Species"] = tax[cls_col].apply(lambda x: extract_rank(x, "s__"))
            out = out.merge(tax[["MAG", "Phylum", "Family", "Genus", "Species"]],
                            on="MAG", how="left")
    for col, default in (("Phylum", "Unknown"), ("Family", ""),
                          ("Genus", ""), ("Species", "")):
        if col not in out.columns:
            out[col] = default
        out[col] = out[col].fillna(default)
    out["label"] = out.apply(
        lambda r: mag_display_label(r["MAG"], r["Genus"], r["Species"],
                                     r["Family"]),
        axis=1)
    return out


def annotate_keystone(
    df: pd.DataFrame,
    keystone_df: pd.DataFrame | None,
) -> pd.DataFrame:
    """添加 is_keystone 列（True/False），基于 keystone_df 里的 MAG 集合。"""
    out = df.copy()
    if keystone_df is not None and not keystone_df.empty:
        ks = keystone_df.copy()
        ks = ks.rename(columns={mag_col(ks): "MAG"})
        ks_set = set(ks["MAG"].astype(str))
        out["is_keystone"] = out["MAG"].astype(str).isin(ks_set)
    else:
        out["is_keystone"] = False
    return out


# ── MAG 子集过滤 ────────────────────────────────────────────

def _select_score(mat: np.ndarray, by: TopNBy) -> np.ndarray:
    if mat.size == 0:
        return np.zeros(0)
    if by == "sum":
        return mat.sum(axis=1)
    if by == "variance":
        return mat.var(axis=1)
    return mat.mean(axis=1)


def apply_filter_mode(
    df: pd.DataFrame,
    *,
    mode: FilterMode = "top_plus_keystone",
    top_n_count: int = 30,
    top_n_by: TopNBy = "mean",
    score_matrix: np.ndarray | None = None,
    score_col: str = "abundance_mean",
) -> pd.DataFrame:
    """按 filter_mode 从 df 里挑子集 MAG。

    - all: 全部 MAG
    - top_n: 按 score 选 Top-N（score 来源：score_matrix 优先，按 top_n_by 打分；
      否则用 df[score_col]）
    - keystone_only: 仅 is_keystone=True 的 MAG
    - top_plus_keystone: Top-N ∪ 全部 keystone

    df 需要有 is_keystone 列（调用前先 annotate_keystone）。
    """
    if mode == "all":
        return df.reset_index(drop=True)

    # 计算 score：优先用 matrix（mag_heatmap 传丰度矩阵）
    if score_matrix is not None and len(score_matrix) == len(df):
        scores = _select_score(score_matrix, top_n_by)
    elif score_col in df.columns:
        scores = df[score_col].to_numpy(dtype=float)
    else:
        scores = np.zeros(len(df))

    # Top-N 索引
    n = min(top_n_count, len(df))
    top_idx = set(np.argsort(-scores, kind="stable")[:n].tolist())
    ks_idx = set(df.index[df.get("is_keystone", False).fillna(False)].tolist()) \
        if "is_keystone" in df.columns else set()

    if mode == "top_n":
        keep = top_idx
    elif mode == "keystone_only":
        keep = ks_idx
    elif mode == "top_plus_keystone":
        keep = top_idx | ks_idx
    else:
        raise ValueError(f"未知 filter_mode={mode!r}")

    out = df.loc[sorted(keep)].reset_index(drop=True)
    if out.empty:
        raise ValueError(
            f"filter_mode={mode!r} 过滤后无 MAG 可显示"
            f"（keystone 数 / top_n_count 过小）")
    return out


# ── 行排序 ──────────────────────────────────────────────────

def _cluster_order(mat: np.ndarray, method: str) -> list[int]:
    """层次聚类叶子顺序；≤2 行时返回原序。"""
    if len(mat) <= 2:
        return list(range(len(mat)))
    try:
        dist = pdist(mat, metric="euclidean")
        link = linkage(dist, method=method)
        return list(leaves_list(link))
    except Exception:
        return list(range(len(mat)))


def order_rows(
    df: pd.DataFrame,
    *,
    mode: RowOrder = "phylum_cluster",
    metric_col: str = "total_completeness",
    abundance_col: str = "abundance_mean",
    cluster_matrix: np.ndarray | None = None,
    linkage_method: str = "average",
    log_transform_cluster: bool = True,
) -> pd.DataFrame:
    """统一行排序：phylum_cluster / metric_desc / abundance。

    - phylum_cluster: 按 Phylum 分组（value_counts 降序）→ 门内用
      cluster_matrix（若提供）做层次聚类
    - metric_desc: 按 df[metric_col] 降序（各图主 metric 不同）
    - abundance: 按 df[abundance_col] 降序
    """
    if mode == "metric_desc" and metric_col in df.columns:
        return df.sort_values(metric_col, ascending=False).reset_index(drop=True)

    if mode == "abundance" and abundance_col in df.columns:
        return df.sort_values(abundance_col, ascending=False).reset_index(drop=True)

    # phylum_cluster（默认兜底）
    if "Phylum" not in df.columns:
        return df.reset_index(drop=True)
    phy_order = df["Phylum"].value_counts().index.tolist()
    reordered: list[int] = []
    have_matrix = (cluster_matrix is not None and len(cluster_matrix) == len(df))
    mat_for_clust = None
    if have_matrix:
        mat_for_clust = (np.log1p(cluster_matrix) if log_transform_cluster
                         else cluster_matrix)
    for phy in phy_order:
        idx = df.index[df["Phylum"] == phy].tolist()
        if not idx:
            continue
        if have_matrix and len(idx) > 2:
            sub = mat_for_clust[idx]
            leaf = _cluster_order(sub, linkage_method)
            reordered.extend([idx[k] for k in leaf])
        else:
            reordered.extend(idx)
    if not reordered:
        return df.reset_index(drop=True)
    return df.iloc[reordered].reset_index(drop=True)


# ── Phylum 彩条 + 图例绘制 ──────────────────────────────────

def draw_phylum_bar(ax_phy, phylum_list: list[str]) -> None:
    """在独立子 Axes 上绘制左侧门彩条。ax_phy 会被关闭坐标轴。"""
    from matplotlib.patches import Rectangle
    for i, phy in enumerate(phylum_list):
        ax_phy.add_patch(Rectangle(
            (0, i - 0.5), 1.0, 1.0,
            facecolor=PHYLUM_COLORS.get(phy, "#888"),
            edgecolor="none", zorder=3,
        ))
    ax_phy.set_xlim(0, 1)
    ax_phy.set_ylim(len(phylum_list) - 0.5, -0.5)
    ax_phy.axis("off")


def draw_phylum_legend(
    ax_leg,
    phylum_list: list[str],
    *,
    title: str = "Phylum",
    loc: str = "upper left",
    bbox_to_anchor: tuple[float, float] = (0.0, 1.0),
    fontsize: float = 6.5,
    title_fontsize: float = 7.5,
) -> plt.matplotlib.legend.Legend:
    """在 ax_leg 上画 Phylum 图例（按出现顺序去重 + 计数）。"""
    seen: list[str] = []
    for phy in phylum_list:
        if phy not in seen:
            seen.append(phy)
    counts = {phy: phylum_list.count(phy) for phy in seen}
    handles = [
        mpatches.Patch(facecolor=PHYLUM_COLORS.get(phy, "#888"),
                       edgecolor="none",
                       label=f"{phy} ({counts[phy]})")
        for phy in seen
    ]
    leg = ax_leg.legend(
        handles=handles, loc=loc, bbox_to_anchor=bbox_to_anchor,
        fontsize=fontsize, title=title, title_fontsize=title_fontsize,
        frameon=False, handlelength=1.2, handleheight=0.8,
        labelspacing=0.4,
    )
    ax_leg.add_artist(leg)
    return leg


def draw_keystone_note(ax_leg, y: float = 0.22,
                       *, fontsize: float = 6.5) -> None:
    """右侧图例区底部加一条 '★ = keystone species' 说明。"""
    ax_leg.text(
        0.0, y, "★ = keystone species",
        transform=ax_leg.transAxes,
        fontsize=fontsize, fontweight="bold",
        va="center", ha="left",
    )
