"""
Figure X — Calibration vs Stress overall_score gap, three KEGG-curated datasets.

Paper 3 Results §X.2 supporting figure. Horizontal bar plot showing the
calibration `STRONG` (1.000) → stress label gap for the three KEGG-curated
datasets (Liu 2023 cold seep, Grettenberger 2021 AMD stream, Ayala 2020 pit
lake), with the cross-topic n=0 rejection annotated for the two non-arsenic
datasets.

Outputs (publication-grade, matplotlib defaults aligned with EnvMeta house style):
- figure_x_calibration_vs_stress.pdf  (vector, primary publication format)
- figure_x_calibration_vs_stress.png  (raster, 600 dpi)
- figure_x_calibration_vs_stress.svg  (editable vector)

Run:
    conda activate envmeta
    python paper/figures/paper3_hypothesis_scoring/figure_x_generate.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path(__file__).parent

# ── 数据（来自 Table 1 + Table 2）──────────────────────────────
DATASETS = [
    {
        "label": "Liu 2023\ncold seep\n(same-topic)",
        "calibration": 1.000,
        "stress_v1": 0.625,
        "stress_v2": 0.250,
        "cross_topic_n0": False,
        "annotation": "v1: B-tier (binary limit)\nv2: A-tier ★",
    },
    {
        "label": "Grettenberger 2021\nAMD stream\n(cross-topic) ★",
        "calibration": 1.000,
        "stress_v1": 0.250,
        "stress_v2": 0.250,  # already A-tier; same as v1
        "cross_topic_n0": True,
        "annotation": "A-tier (already)\ncross-topic n=0 ★",
    },
    {
        "label": "Ayala 2020\npit lake\n(cross-topic) ★",
        "calibration": 1.000,
        "stress_v1": 0.455,
        "stress_v2": 0.182,
        "cross_topic_n0": True,
        "annotation": "v1: B-tier (binary limit)\nv2: A-tier ★",
    },
]

# ── 出版级 matplotlib 配置 ─────────────────────────────────────
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,  # TrueType for editable PDF
    "ps.fonttype": 42,
})

CALIB_COLOR = "#2E86AB"      # 深蓝 (calibration)
STRESS_V1_COLOR = "#F4A261"  # 橙 (stress v1, default thresholds)
STRESS_V2_COLOR = "#E63946"  # 红 (stress v2, dominance-aware)
THRESHOLD_GREY = "#888888"


def make_figure():
    n = len(DATASETS)
    fig, ax = plt.subplots(figsize=(8.5, 4.5))

    bar_height = 0.25
    y_calib = np.arange(n) + bar_height
    y_v1 = np.arange(n)
    y_v2 = np.arange(n) - bar_height

    calib_vals = [d["calibration"] for d in DATASETS]
    v1_vals = [d["stress_v1"] for d in DATASETS]
    v2_vals = [d["stress_v2"] for d in DATASETS]

    # bars
    ax.barh(y_calib, calib_vals, height=bar_height,
            color=CALIB_COLOR, label="Calibration",
            edgecolor="white", linewidth=0.5, zorder=3)
    ax.barh(y_v1, v1_vals, height=bar_height,
            color=STRESS_V1_COLOR, label="Stress v1 (default thresholds)",
            edgecolor="white", linewidth=0.5, zorder=3)
    ax.barh(y_v2, v2_vals, height=bar_height,
            color=STRESS_V2_COLOR,
            label="Stress v2 (dominance-aware, v0.9.x)",
            edgecolor="white", linewidth=0.5, zorder=3)

    # value labels on bars
    for i, (cv, v1, v2) in enumerate(zip(calib_vals, v1_vals, v2_vals)):
        ax.text(cv + 0.015, y_calib[i], f"{cv:.3f}",
                va="center", fontsize=7.5, color=CALIB_COLOR, fontweight="bold")
        ax.text(v1 + 0.015, y_v1[i], f"{v1:.3f}",
                va="center", fontsize=7.5, color=STRESS_V1_COLOR, fontweight="bold")
        ax.text(v2 + 0.015, y_v2[i], f"{v2:.3f}",
                va="center", fontsize=7.5, color=STRESS_V2_COLOR, fontweight="bold")

    # threshold lines
    ax.axvline(0.75, color=THRESHOLD_GREY, linestyle="--", linewidth=0.7,
               alpha=0.6, zorder=1)
    ax.axvline(0.40, color=THRESHOLD_GREY, linestyle=":", linewidth=0.7,
               alpha=0.6, zorder=1)
    ax.text(0.75, n - 0.3, "strong\n0.75", ha="center", fontsize=7,
            color=THRESHOLD_GREY)
    ax.text(0.40, n - 0.3, "suggestive\n0.40", ha="center", fontsize=7,
            color=THRESHOLD_GREY)

    # discrimination grade annotations (right side)
    for i, d in enumerate(DATASETS):
        ax.text(1.08, i, d["annotation"], va="center", fontsize=7.5,
                color="#444444",
                fontweight="bold" if d["cross_topic_n0"] else "normal")

    # axes
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels([d["label"] for d in DATASETS], fontsize=8.5)
    ax.set_xlabel("overall_score", fontsize=9)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(-0.7, n - 0.3 + 0.3)
    ax.invert_yaxis()  # top dataset first

    ax.legend(loc="lower right", frameon=False, fontsize=7.5)

    ax.set_title(
        "Figure X. Calibration vs stress overall_score across three KEGG-curated datasets.\n"
        "v1 = default thresholds (binary mean_completeness ≥ 50%); v2 = dominance-aware (v0.9.x min_dominance_fraction = 0.20).\n"
        "Cross-topic arsenate_reduction rejected with n=0 active MAGs in 2/2 non-arsenic datasets ★.",
        fontsize=9, loc="left", pad=10,
    )

    plt.subplots_adjust(left=0.16, right=0.74, top=0.80, bottom=0.13)
    return fig


def main():
    fig = make_figure()
    base = OUT_DIR / "figure_x_calibration_vs_stress"
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    fig.savefig(f"{base}.png", bbox_inches="tight", dpi=600)
    fig.savefig(f"{base}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved:")
    for ext in ("pdf", "png", "svg"):
        p = base.with_suffix(f".{ext}")
        size = p.stat().st_size / 1024
        print(f"  {p.name}: {size:.1f} KB")


if __name__ == "__main__":
    main()
