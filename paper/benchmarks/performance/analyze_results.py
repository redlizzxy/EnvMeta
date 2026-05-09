"""Combine all benchmark results into final tables + scaling figure.

Inputs (in paper/benchmarks/performance/results/):
- sample_baseline_runtime.tsv  (sample_data, 14 figures x 3 repeats)
- liu_full_runtime.tsv         (Liu 1084x87, 11 figures x 3 repeats; 3 fail)
- liu_sweep_m{N}s{S}.tsv       (Liu subsample sweep, 9 cells x 3 figures x 2 repeats)
- (optional) liu_synth_dense_*.tsv (synthetic dense annotation tests)

Outputs (in paper/benchmarks/performance/):
- table_per_figure.tsv         (all figure x dataset combinations, status + runtime)
- scaling_curve.{pdf,png,svg}  (log-log: N_MAG x N_sample x runtime, per figure)
- hardware_sizing.md           (text table: data scale -> RAM/CPU recommendation)
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "paper" / "benchmarks" / "performance" / "results"
OUT = ROOT / "paper" / "benchmarks" / "performance"


def load_all() -> pd.DataFrame:
    files = sorted(RES.glob("*_runtime.tsv")) + sorted(RES.glob("liu_sweep_*.tsv"))
    files = list(dict.fromkeys(files))  # dedupe
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        df["source"] = f.stem
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def per_figure_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Wide table: figure x dataset combos."""
    keep = df[["dataset", "figure", "n_mags", "n_samples", "n_groups",
               "wall_s_median", "peak_mem_mb_max", "status", "error"]].copy()
    keep = keep.sort_values(["figure", "n_mags", "n_samples"])
    return keep


def _classify_dataset(ds: str) -> str:
    """Classify a dataset string into: 'sample-dense', 'liu-sparse', 'liu-dense-synth'."""
    if "synth_dense" in ds:
        return "Liu (synthetic dense, 25 KO/MAG, 4 env)"
    if ds == "sample":
        return "sample (real dense, 30 KO/MAG, 4 env)"
    return "Liu (real sparse, 1.5 KO/MAG, 1 env)"


def make_scaling_figure(df: pd.DataFrame, out_basename: str) -> None:
    """Plot wall_s vs (n_mags x n_samples) per figure for cycle/mag_heatmap/pathway."""
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "DejaVu Sans"

    sub = df[(df["figure"].isin(["cycle_diagram", "mag_heatmap", "pathway"]))
             & (df["status"] == "ok")].copy()
    sub["scale"] = sub["n_mags"] * sub["n_samples"]
    sub = sub.dropna(subset=["wall_s_median", "scale"])
    sub["category"] = sub["dataset"].apply(_classify_dataset)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.7), constrained_layout=True)
    figures = ["cycle_diagram", "mag_heatmap", "pathway"]
    color_map = {
        "sample (real dense, 30 KO/MAG, 4 env)":         "#2E86AB",
        "Liu (real sparse, 1.5 KO/MAG, 1 env)":          "#9E9E9E",
        "Liu (synthetic dense, 25 KO/MAG, 4 env)":       "#E63946",
    }
    marker_map = {
        "sample (real dense, 30 KO/MAG, 4 env)":         "o",
        "Liu (real sparse, 1.5 KO/MAG, 1 env)":          "^",
        "Liu (synthetic dense, 25 KO/MAG, 4 env)":       "s",
    }

    for ax, fig_name in zip(axes, figures):
        sf = sub[sub["figure"] == fig_name].copy()
        if sf.empty:
            ax.set_title(f"{fig_name} (no data)")
            continue
        for cat, sg in sf.groupby("category"):
            sg = sg.sort_values("scale")
            ax.scatter(sg["scale"], sg["wall_s_median"],
                       label=cat, c=color_map.get(cat, "#000"),
                       marker=marker_map.get(cat, "o"),
                       s=90, alpha=0.85, edgecolor="white", linewidth=1.3)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("N_MAGs x N_samples", fontsize=11)
        ax.set_ylabel("Wall time (s, median of 2-3 reps)", fontsize=11)
        ax.set_title(fig_name, fontsize=12, fontweight="bold")
        ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.5)

    # Legend at bottom (shared)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9,
               bbox_to_anchor=(0.5, -0.04), frameon=False)

    fig.suptitle("EnvMeta runtime scaling: dense annotation x N_env dominates cycle_diagram cost",
                 fontsize=12.5, fontweight="bold")
    for ext in ("pdf", "png", "svg"):
        path = OUT / f"{out_basename}.{ext}"
        plt.savefig(path, dpi=300 if ext == "png" else None,
                    bbox_inches="tight")
    print(f"[fig] -> {out_basename}.{{pdf,png,svg}}")


def hardware_sizing_md(df: pd.DataFrame) -> str:
    """Build the hardware sizing markdown table from observed data."""
    ok = df[df["status"] == "ok"].copy()
    agg = ok.groupby(["dataset", "n_mags", "n_samples"]).agg(
        total_wall_s=("wall_s_median", "sum"),
        max_mem_mb=("peak_mem_mb_max", "max"),
        n_figures_ok=("figure", "count"),
    ).reset_index().sort_values(["n_mags", "n_samples"])

    # Manual markdown rendering (avoid tabulate dependency)
    cols = ["dataset", "n_mags", "n_samples", "total_wall_s", "max_mem_mb", "n_figures_ok"]
    header = "| " + " | ".join(cols) + " |"
    sep = "|" + "|".join(["---"] * len(cols)) + "|"
    rows = []
    for _, r in agg.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (int, np.integer)):
                cells.append(str(int(v)))
            elif isinstance(v, (float, np.floating)):
                cells.append(f"{v:.2f}")
            else:
                cells.append(str(v))
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join([header, sep] + rows)


def main() -> int:
    if not RES.exists():
        print("[err] no results dir found")
        return 1
    df = load_all()
    if df.empty:
        print("[err] no benchmark results loaded")
        return 1
    print(f"[loaded] {len(df)} rows from {df['source'].nunique()} files")

    # Per-figure summary table
    summary = per_figure_summary(df)
    summary.to_csv(OUT / "table_per_figure.tsv", sep="\t", index=False)
    print(f"[csv]  -> table_per_figure.tsv  ({len(summary)} rows)")

    # Scaling figure
    make_scaling_figure(df, out_basename="scaling_curve")

    # Hardware sizing
    md = hardware_sizing_md(df)
    print("\n=== Per-dataset aggregate ===")
    print(md)
    (OUT / "hardware_sizing_observed.md").write_text(
        "# Per-dataset observed totals\n\n" + md + "\n", encoding="utf-8"
    )
    print(f"\n[md] -> hardware_sizing_observed.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
