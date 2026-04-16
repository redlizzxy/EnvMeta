"""MAG 通路完整度（heatmap / bubble）。

对标 scripts/python/08_pathway_completeness.py。核心：以知识库（elements.json）
定义的 18 个元素循环通路为列，每个 MAG 按"拥有 KO / 通路 KO 总数"计算完整度。

输入：
    ko_annotation_df: MAG × KO 注释（长表：MAG + KEGG_ko 两列，KEGG_ko 可含
                       `ko:` 前缀与 `,` 分隔）
    taxonomy_df:      可选，MAG + classification（GTDB）
    keystone_df:      可选，MAG + Genus
    abundance_df:     可选，MAG × sample（用于 Top-N 过滤与气泡图大小权重）

输出：
    figure — heatmap 或 bubble（按 style 参数）
    stats  — MAG × pathway 完整度长表 + 元素标签 + keystone 标记
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

from envmeta.analysis.base import AnalysisResult
from envmeta.geocycle.knowledge_base import (
    element_colors, pathway_display, pathway_element_map, pathway_ko_sets,
)

DEFAULTS = {
    "width_mm": 280,
    "height_mm": 200,
    "style": "heatmap",          # "heatmap" | "bubble"
    "max_mags": None,            # None → 全部；int → 按丰度/完整度取 Top-N
    "min_completeness": 0.0,     # 过滤：MAG 的总完整度低于此值不显示
    "sort_by": "phylum_then_total",  # "phylum_then_total" | "total" | "abundance"
    "element_filter": None,      # None → 全部；["arsenic", "sulfur", ...]
    "annotate_keystone": True,
    "annotate_top_n": 10,        # Top-N 丰度 MAG 在热图前缀标星
    "bubble_scale": 4.0,
    "cmap_name": "YlGnBu",
    "show_phylum_stripe": True,
}

# Phylum → 颜色（和 mag_quality 共享）
PHYLUM_COLORS = {
    "Pseudomonadota":    "#E74C3C",
    "Acidobacteriota":   "#3498DB",
    "Chloroflexota":     "#2ECC71",
    "Bacteroidota_A":    "#F39C12",
    "Patescibacteriota": "#9B59B6",
    "Desulfobacterota":  "#1ABC9C",
    "Actinomycetota":    "#E67E22",
    "Gemmatimonadota":   "#34495E",
    "Planctomycetota":   "#D35400",
    "Myxococcota":       "#8E44AD",
    "Nitrospirota":      "#27AE60",
    "Zixibacteria":      "#C0392B",
    "Other":             "#BDC3C7",
    "Unknown":           "#888888",
}


def _mag_col(df: pd.DataFrame) -> str:
    for c in ("MAG", "Name", "user_genome", "Genome", "Bin"):
        if c in df.columns:
            return c
    return df.columns[0]


def _extract_phylum(cls: str) -> str:
    if not isinstance(cls, str):
        return "Unknown"
    for part in cls.split(";"):
        part = part.strip()
        if part.startswith("p__"):
            return part[3:] or "Unknown"
    return "Unknown"


def _parse_ko_annotation(df: pd.DataFrame) -> dict[str, set[str]]:
    """把长表 MAG+KEGG_ko 解析为 {MAG: {ko, ...}}。
    KEGG_ko 可能含 `ko:K00001` 前缀与 `,` 分隔多个。"""
    mag_c = _mag_col(df)
    ko_c = next((c for c in df.columns
                 if c.lower() in ("kegg_ko", "kegg ko", "ko", "ko_id")),
                df.columns[1] if len(df.columns) > 1 else None)
    if ko_c is None:
        raise ValueError("无法定位 KEGG_ko / KO 列")
    mag_kos: dict[str, set[str]] = {}
    for _, r in df.iterrows():
        mag = str(r[mag_c])
        raw = str(r[ko_c])
        if raw in ("", "nan", "None"):
            continue
        for ko in raw.split(","):
            k = ko.strip().replace("ko:", "")
            if k and k.startswith("K"):
                mag_kos.setdefault(mag, set()).add(k)
    return mag_kos


def _compute_completeness(
    mag_kos: dict[str, set[str]],
    pw_kos: dict[str, list[str]],
    all_mags: list[str],
) -> pd.DataFrame:
    rows = []
    for m in all_mags:
        owned = mag_kos.get(m, set())
        row = {"MAG": m}
        for pw_id, kos in pw_kos.items():
            detected = owned & set(kos)
            row[pw_id] = len(detected) / len(kos) * 100 if kos else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def analyze(
    ko_annotation_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    abundance_df: pd.DataFrame | None = None,
    params: dict | None = None,
) -> AnalysisResult:
    p = {**DEFAULTS, **(params or {})}

    # === 解析 KO 注释 → {MAG: {KO}} ===
    mag_kos = _parse_ko_annotation(ko_annotation_df)

    # === 通路定义（从 KB）===
    pw_kos = pathway_ko_sets()
    pw_elem = pathway_element_map()
    pw_display = pathway_display(lang="en")
    elem_color = element_colors()
    if p["element_filter"]:
        keep = set(p["element_filter"])
        pw_kos = {k: v for k, v in pw_kos.items() if pw_elem[k] in keep}
    pw_order = list(pw_kos.keys())  # KB 定义顺序（按元素聚类）

    # === 所有 MAG 列表：优先用 taxonomy；兜底用 KO 注释里出现的 MAG ===
    if taxonomy_df is not None and not taxonomy_df.empty:
        tax = taxonomy_df.copy()
        tmag = _mag_col(tax)
        tax = tax.rename(columns={tmag: "MAG"})
        class_col = next((c for c in tax.columns
                          if "classif" in c.lower() or "taxonomy" in c.lower()),
                         tax.columns[1] if len(tax.columns) > 1 else None)
        if class_col:
            tax["Phylum"] = tax[class_col].apply(_extract_phylum)
        else:
            tax["Phylum"] = "Unknown"
        all_mags = sorted(tax["MAG"].astype(str).unique().tolist())
    else:
        tax = None
        all_mags = sorted(mag_kos.keys())

    if not all_mags:
        raise ValueError("无 MAG 可分析（检查 KO 注释表非空）")

    df = _compute_completeness(mag_kos, pw_kos, all_mags)

    # 合并 phylum
    if tax is not None:
        df = df.merge(tax[["MAG", "Phylum"]], on="MAG", how="left")
    df["Phylum"] = df.get("Phylum", pd.Series(["Unknown"] * len(df))).fillna("Unknown")

    # keystone
    if keystone_df is not None and not keystone_df.empty:
        ks = keystone_df.copy()
        ks = ks.rename(columns={_mag_col(ks): "MAG"})
        df["is_keystone"] = df["MAG"].isin(set(ks["MAG"].astype(str)))
    else:
        df["is_keystone"] = False

    # 丰度均值（用于 Top-N 筛选与气泡图）
    if abundance_df is not None and not abundance_df.empty:
        ab = abundance_df.copy()
        ab_mag = _mag_col(ab)
        ab = ab.rename(columns={ab_mag: "MAG"})
        sample_cols = [c for c in ab.columns if c != "MAG"]
        ab["abundance_mean"] = ab[sample_cols].apply(
            pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)
        df = df.merge(ab[["MAG", "abundance_mean"]], on="MAG", how="left")
        df["abundance_mean"] = df["abundance_mean"].fillna(0.0)
    else:
        df["abundance_mean"] = 0.0

    # 总完整度
    df["total_completeness"] = df[pw_order].sum(axis=1)

    # 过滤
    df = df[df["total_completeness"] >= p["min_completeness"]].copy()
    if df.empty:
        raise ValueError(f"无 MAG 总完整度 ≥ {p['min_completeness']}")

    # 排序
    if p["sort_by"] == "abundance":
        df = df.sort_values("abundance_mean", ascending=False)
    elif p["sort_by"] == "total":
        df = df.sort_values("total_completeness", ascending=False)
    else:  # phylum_then_total
        phy_rank = {ph: i for i, ph in
                    enumerate(df["Phylum"].value_counts().index.tolist())}
        df["_phy_rank"] = df["Phylum"].map(phy_rank)
        df = df.sort_values(["_phy_rank", "total_completeness"],
                            ascending=[True, False]).drop(columns=["_phy_rank"])
    df = df.reset_index(drop=True)

    # Top-N 过滤
    if p["max_mags"] and len(df) > p["max_mags"]:
        df = df.head(p["max_mags"]).reset_index(drop=True)

    # === 绘图 ===
    if p["style"] == "bubble":
        fig = _draw_bubble(df, pw_order, pw_display, pw_elem, elem_color, p)
    else:
        fig = _draw_heatmap(df, pw_order, pw_display, pw_elem, elem_color, p)

    # === stats ===
    stats_df = df[["MAG", "Phylum", "is_keystone", "abundance_mean",
                   "total_completeness"] + pw_order].copy()
    return AnalysisResult(figure=fig, stats=stats_df, params=p)


def _draw_heatmap(df, pw_order, pw_display, pw_elem, elem_color, p) -> plt.Figure:
    n_mag = len(df)
    n_pw = len(pw_order)
    fig, ax = plt.subplots(
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4),
        constrained_layout=True,
    )
    mat = df[pw_order].to_numpy()
    cmap = plt.get_cmap(p["cmap_name"])
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=100,
                   interpolation="nearest")

    # Phylum 彩条
    if p["show_phylum_stripe"]:
        for i, (_, r) in enumerate(df.iterrows()):
            color = PHYLUM_COLORS.get(r["Phylum"], "#888")
            ax.add_patch(Rectangle((-0.6, i - 0.5), 0.4, 1.0,
                                   color=color, clip_on=False, zorder=3))

    # 元素彩条（列顶部）
    elem_of_col = [pw_elem[pw] for pw in pw_order]
    for j, el in enumerate(elem_of_col):
        ax.add_patch(Rectangle((j - 0.5, -0.9), 1.0, 0.4,
                               color=elem_color.get(el, "#888"),
                               clip_on=False, zorder=3))

    # keystone 标记（在 MAG label 前加 ★）
    y_labels = []
    for _, r in df.iterrows():
        tag = "★ " if (p["annotate_keystone"] and r["is_keystone"]) else ""
        y_labels.append(f"{tag}{r['MAG']}")
    ax.set_yticks(range(n_mag))
    ax.set_yticklabels(y_labels, fontsize=max(3, min(8, 400 // max(n_mag, 1))))

    ax.set_xticks(range(n_pw))
    ax.set_xticklabels([pw_display.get(pw, pw) for pw in pw_order],
                       rotation=45, ha="right", fontsize=8)
    ax.set_title(
        f"Pathway Completeness (MAG × pathway)  n={n_mag} MAGs × {n_pw} pathways",
        fontsize=11, fontweight="bold", pad=12,
    )
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, label="Completeness (%)")
    cbar.ax.tick_params(labelsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    return fig


def _draw_bubble(df, pw_order, pw_display, pw_elem, elem_color, p) -> plt.Figure:
    n_mag = len(df)
    n_pw = len(pw_order)
    fig, ax = plt.subplots(
        figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4),
        constrained_layout=True,
    )

    mat = df[pw_order].to_numpy()          # 完整度 → 颜色
    abund = df["abundance_mean"].to_numpy()  # 丰度 → 大小
    ab_scale = abund / (abund.max() + 1e-9) * p["bubble_scale"] * 80 + 8

    cmap = plt.get_cmap(p["cmap_name"])
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=0, vmax=100)
    for i in range(n_mag):
        for j in range(n_pw):
            val = mat[i, j]
            if val <= 0:
                continue
            ax.scatter(j, i, s=ab_scale[i] * (val / 100) + 10,
                       c=[cmap(val / 100)], edgecolors="white", linewidths=0.4,
                       alpha=0.85, zorder=3)

    # 元素彩条
    for j, pw in enumerate(pw_order):
        el = pw_elem[pw]
        ax.add_patch(Rectangle((j - 0.5, -0.9), 1.0, 0.4,
                               color=elem_color.get(el, "#888"),
                               clip_on=False, zorder=2))

    ax.set_xticks(range(n_pw))
    ax.set_xticklabels([pw_display.get(pw, pw) for pw in pw_order],
                       rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n_mag))
    y_labels = [("★ " if r["is_keystone"] else "") + r["MAG"] for _, r in df.iterrows()]
    ax.set_yticklabels(y_labels, fontsize=max(4, min(9, 380 // max(n_mag, 1))))
    ax.set_xlim(-1, n_pw)
    ax.set_ylim(-1.5, n_mag)
    ax.invert_yaxis()
    ax.set_title(
        "Pathway Completeness Bubble  "
        "(color=completeness %, bubble size ≈ abundance × completeness)",
        fontsize=10, fontweight="bold", pad=12,
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # 颜色 colorbar（completeness %）
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.02,
                        label="Completeness (%)")
    cbar.ax.tick_params(labelsize=8)

    # 图例：大小（与 abundance × completeness 的代表点对应）
    from matplotlib.lines import Line2D
    legend_vals = [20, 50, 100]
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#555",
               markersize=np.sqrt(v * 0.4 * p["bubble_scale"]),
               markeredgecolor="white", markeredgewidth=0.3,
               label=f"{v}% × top-abund.")
        for v in legend_vals
    ]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.18, 1.0),
              fontsize=7, frameon=False,
              title="bubble size", title_fontsize=7)
    return fig
