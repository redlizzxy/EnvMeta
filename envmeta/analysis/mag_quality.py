"""MAG 质量评估散点图（对标 scripts/python/06_MAG_quality.py）。

X=Completeness, Y=Contamination，点按门着色，大小按基因组大小，Keystone
物种用菱形高亮 + 标签。

输入：
    quality_df:   CheckM2 质量表（含 Name/MAG + Completeness + Contamination + Genome_Size）
    taxonomy_df:  可选，MAG + GTDB classification 字符串（`d__...;p__...;...`）
    keystone_df:  可选，MAG + Genus + Phylum（keystone 物种列表）

输出：
    figure — 散点图 + 质量分区背景 + 右侧图例（门 + 基因组大小 + keystone 标记）
    stats  — 按 High/Medium/Low 汇总 + 全 MAG 质量分类明细
"""
from __future__ import annotations

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from envmeta.analysis.base import AnalysisResult

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
}

DEFAULTS = {
    "width_mm": 260,
    "height_mm": 140,
    "high_completeness": 90.0,
    "high_contamination": 5.0,
    "med_completeness": 50.0,
    "med_contamination": 10.0,
    "size_scale": 12.0,         # 基因组大小（Mb）→ 点大小的倍率
    "min_phylum_count": 3,      # 少于此数的门归为 Other
    "show_keystone_labels": True,
    "xlim": (35.0, 102.0),
    "ylim": None,               # None → 自动
}


def _mag_col(df: pd.DataFrame) -> str:
    for c in ("MAG", "Name", "Bin", "Genome"):
        if c in df.columns:
            return c
    return df.columns[0]


def _extract_phylum(classification: str) -> str:
    if not isinstance(classification, str):
        return "Unknown"
    for part in classification.split(";"):
        part = part.strip()
        if part.startswith("p__"):
            return part[3:] or "Unknown"
    return "Unknown"


def _quality_class(row, hc, hcon, mc, mcon) -> str:
    if row["Completeness"] >= hc and row["Contamination"] < hcon:
        return "High"
    if row["Completeness"] >= mc and row["Contamination"] < mcon:
        return "Medium"
    return "Low"


def analyze(
    quality_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    params: dict | None = None,
) -> AnalysisResult:
    p = {**DEFAULTS, **(params or {})}

    df = quality_df.copy()
    mag_c = _mag_col(df)
    df = df.rename(columns={mag_c: "MAG"})
    required = {"Completeness", "Contamination", "Genome_Size"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"quality 表缺少必需列：{missing}")
    for c in required:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=list(required)).reset_index(drop=True)
    if df.empty:
        raise ValueError("quality 表无有效数据")

    # Phylum
    if taxonomy_df is not None and not taxonomy_df.empty:
        tax = taxonomy_df.copy()
        tax_mag = _mag_col(tax)
        tax = tax.rename(columns={tax_mag: "MAG"})
        class_col = next((c for c in tax.columns
                          if "classif" in c.lower() or "taxonomy" in c.lower()),
                         tax.columns[1] if len(tax.columns) > 1 else None)
        if class_col is None:
            raise ValueError("taxonomy 表需要 classification / Taxonomy 列")
        tax["Phylum"] = tax[class_col].apply(_extract_phylum)
        df = df.merge(tax[["MAG", "Phylum"]], on="MAG", how="left")
    df["Phylum"] = df.get("Phylum", pd.Series(["Unknown"] * len(df))).fillna("Unknown")

    # 小门合并为 Other
    phy_counts = df["Phylum"].value_counts()
    minor = set(phy_counts[phy_counts < p["min_phylum_count"]].index)
    df["Phylum_group"] = df["Phylum"].apply(
        lambda x: "Other" if (x in minor or x not in PHYLUM_COLORS) else x
    )

    # Keystone
    if keystone_df is not None and not keystone_df.empty:
        ks = keystone_df.copy()
        ks_mag = _mag_col(ks)
        ks = ks.rename(columns={ks_mag: "MAG"})
        df["is_keystone"] = df["MAG"].isin(set(ks["MAG"]))
        if "Genus" in ks.columns:
            df["ks_genus"] = df["MAG"].map(dict(zip(ks["MAG"], ks["Genus"])))
        else:
            df["ks_genus"] = np.nan
    else:
        df["is_keystone"] = False
        df["ks_genus"] = np.nan

    # Quality
    df["Quality"] = df.apply(
        lambda r: _quality_class(r, p["high_completeness"], p["high_contamination"],
                                 p["med_completeness"], p["med_contamination"]),
        axis=1,
    )

    # === 绘图 ===
    fig = plt.figure(figsize=(p["width_mm"] / 25.4, p["height_mm"] / 25.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.3], wspace=0.05)
    ax = fig.add_subplot(gs[0])
    ax_leg = fig.add_subplot(gs[1])
    ax_leg.axis("off")

    # 质量分区背景 + 分界线
    ax.axhspan(0, p["high_contamination"], alpha=0.06, color="#27AE60", zorder=0)
    ax.axhspan(p["high_contamination"], p["med_contamination"],
               alpha=0.04, color="#F39C12", zorder=0)
    ax.axhline(p["high_contamination"], color="#27AE60", ls="--", lw=0.8, alpha=0.6)
    ax.axhline(p["med_contamination"], color="#E74C3C", ls="--", lw=0.8, alpha=0.6)
    ax.axvline(p["high_completeness"], color="#27AE60", ls="--", lw=0.8, alpha=0.6)
    ax.axvline(p["med_completeness"], color="#E74C3C", ls="--", lw=0.8, alpha=0.6)
    ax.text(p["high_completeness"] + 2, 0.3, "High Quality", fontsize=7,
            color="#27AE60", fontstyle="italic", alpha=0.7)
    ax.text((p["med_completeness"] + p["high_completeness"]) / 2, p["high_contamination"] + 1,
            "Medium Quality", fontsize=7, color="#F39C12", alpha=0.7)

    # 散点（非 keystone）
    for phy, grp in df[~df["is_keystone"]].groupby("Phylum_group"):
        color = PHYLUM_COLORS.get(phy, "#BDC3C7")
        sizes = grp["Genome_Size"] / 1e6 * p["size_scale"]
        ax.scatter(grp["Completeness"], grp["Contamination"],
                   s=sizes, c=color, alpha=0.55, edgecolors="white",
                   linewidths=0.3, zorder=2)

    # 散点（keystone，菱形 + 黑边 + 标签）
    for _, row in df[df["is_keystone"]].iterrows():
        color = PHYLUM_COLORS.get(row["Phylum_group"], "#BDC3C7")
        size = row["Genome_Size"] / 1e6 * p["size_scale"] * 1.8
        ax.scatter(row["Completeness"], row["Contamination"],
                   s=size, c=color, alpha=0.9, edgecolors="black",
                   linewidths=1.2, zorder=4, marker="D")
        if p["show_keystone_labels"]:
            label = row["ks_genus"] if pd.notna(row.get("ks_genus")) else row["MAG"]
            if label == "Unknown" or pd.isna(label):
                label = row["MAG"]
            ax.annotate(
                label, (row["Completeness"], row["Contamination"]),
                xytext=(row["Completeness"] + 1.5, row["Contamination"] + 0.4),
                fontsize=6, fontstyle="italic", fontweight="bold", color="#222",
                arrowprops=dict(arrowstyle="-", color="#888", lw=0.4),
                zorder=5,
            )

    ax.set_xlim(*p["xlim"])
    if p["ylim"]:
        ax.set_ylim(*p["ylim"])
    else:
        ax.set_ylim(-0.5, max(df["Contamination"].max() + 1, 12))
    ax.set_xlabel("Completeness (%)", fontsize=10)
    ax.set_ylabel("Contamination (%)", fontsize=10)
    hq = (df["Quality"] == "High").sum()
    mq = (df["Quality"] == "Medium").sum()
    lq = (df["Quality"] == "Low").sum()
    ax.set_title(
        f"MAG Quality Assessment  (n={len(df)}; High {hq} / Medium {mq} / Low {lq})",
        fontsize=11, fontweight="bold", pad=10,
    )
    ax.tick_params(labelsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # 图例 ── 门
    phy_order = df["Phylum_group"].value_counts().index.tolist()
    phy_handles = [
        mpatches.Patch(facecolor=PHYLUM_COLORS.get(phy, "#BDC3C7"), edgecolor="none",
                       label=f"{phy} ({(df['Phylum_group'] == phy).sum()})")
        for phy in phy_order
    ]
    leg1 = ax_leg.legend(handles=phy_handles, loc="upper left",
                         bbox_to_anchor=(0.0, 1.0), fontsize=7, title="Phylum",
                         title_fontsize=8, frameon=False,
                         handlelength=1.2, handleheight=0.8)
    ax_leg.add_artist(leg1)

    # 图例 ── 基因组大小
    size_vals = [1, 3, 6, 9]
    size_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#999",
               markersize=np.sqrt(s * p["size_scale"] * 0.8),
               label=f"{s} Mb", markeredgecolor="white", markeredgewidth=0.3)
        for s in size_vals
    ]
    leg2 = ax_leg.legend(handles=size_handles, loc="center left",
                         bbox_to_anchor=(0.0, 0.35), fontsize=7,
                         title="Genome size", title_fontsize=8, frameon=False,
                         handlelength=1.5)
    ax_leg.add_artist(leg2)

    # 图例 ── Keystone
    if df["is_keystone"].any():
        ks_h = Line2D([0], [0], marker="D", color="none", markerfacecolor="#999",
                      markersize=8, markeredgecolor="black", markeredgewidth=1.0,
                      label="Keystone species")
        ax_leg.legend(handles=[ks_h], loc="lower left",
                      bbox_to_anchor=(0.0, 0.05), fontsize=7, frameon=False)


    # === stats：summary + 明细 ===
    summary_rows = []
    for q in ("High", "Medium", "Low"):
        sub = df[df["Quality"] == q]
        summary_rows.append({
            "Quality": q,
            "Count": int(len(sub)),
            "Completeness_mean": float(sub["Completeness"].mean()) if len(sub) else np.nan,
            "Contamination_mean": float(sub["Contamination"].mean()) if len(sub) else np.nan,
            "Genome_Size_mean_Mb": float(sub["Genome_Size"].mean() / 1e6) if len(sub) else np.nan,
        })
    detail_cols = [c for c in
                   ["MAG", "Completeness", "Contamination", "Genome_Size",
                    "Phylum", "Quality", "is_keystone"]
                   if c in df.columns]
    stats_df = pd.concat([
        pd.DataFrame(summary_rows).assign(type="summary"),
        df[detail_cols].assign(type="detail").sort_values(
            ["Quality", "Completeness"], ascending=[True, False]),
    ], ignore_index=True, sort=False)

    return AnalysisResult(figure=fig, stats=stats_df, params=p)
