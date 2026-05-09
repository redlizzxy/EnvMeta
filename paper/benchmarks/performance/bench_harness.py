"""EnvMeta performance benchmark harness.

Measures wall-clock time + peak RSS delta for each of EnvMeta's 14 figures
on a configurable dataset. Background sampling thread captures peak RSS
during each call.

Usage:
    python paper/benchmarks/performance/bench_harness.py --dataset sample
    python paper/benchmarks/performance/bench_harness.py --dataset liu
    python paper/benchmarks/performance/bench_harness.py --dataset liu --subsample 200,30
    python paper/benchmarks/performance/bench_harness.py --dataset liu --figures cycle_diagram

Outputs:
    paper/benchmarks/performance/results/<dataset>[_<subsample>]_runtime.tsv
"""
from __future__ import annotations

import argparse
import gc
import json
import sys
import threading
import time
import traceback
from contextlib import contextmanager
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Headless rendering
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from envmeta.analysis import (
    alpha_boxplot, cycle_diagram, gene_heatmap, gene_profile, lefse, log2fc,
    mag_heatmap, mag_quality, network, pathway, pcoa, rda, stackplot,
)
from envmeta.geocycle import hypothesis as hyp_mod


# ──────────────────────────────────────────────────────────────────
# Memory profiler (background sampling thread)
# ──────────────────────────────────────────────────────────────────

class MemSampler:
    def __init__(self, interval: float = 0.05):
        self.interval = interval
        self.proc = psutil.Process()
        self.peak = 0
        self.baseline = 0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def _sample(self):
        while not self._stop.is_set():
            try:
                rss = self.proc.memory_info().rss
                if rss > self.peak:
                    self.peak = rss
            except Exception:
                pass
            self._stop.wait(self.interval)

    def __enter__(self):
        gc.collect()
        self.baseline = self.proc.memory_info().rss
        self.peak = self.baseline
        self._stop.clear()
        self._thread = threading.Thread(target=self._sample, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=2)
        plt.close("all")
        gc.collect()

    @property
    def delta_mb(self) -> float:
        return (self.peak - self.baseline) / 1024 / 1024


@contextmanager
def timed():
    t0 = time.perf_counter()
    yield (lambda: time.perf_counter() - t0)


# ──────────────────────────────────────────────────────────────────
# Dataset loaders
# ──────────────────────────────────────────────────────────────────

SAMPLE = ROOT / "tests" / "sample_data"
LIU_BASE = ROOT / "paper" / "benchmarks" / "external" / "liu_2023_coldseep" / "input_data_local"
LIU_RESHAPE = ROOT / "paper" / "benchmarks" / "performance" / "liu_reshaped"


def load_sample():
    """Load sample_data EnvMeta inputs into a dict."""
    return {
        "abundance_taxonomy": pd.read_csv(SAMPLE / "Phylum.txt", sep="\t"),
        "metadata": pd.read_csv(SAMPLE / "metadata.txt", sep="\t"),
        "alpha": pd.read_csv(SAMPLE / "alpha.txt", sep="\t"),
        "distance": pd.read_csv(SAMPLE / "beta_bray.txt", sep="\t"),
        "env": pd.read_csv(SAMPLE / "env_factors.txt", sep="\t"),
        "ko_wide": pd.read_csv(SAMPLE / "ko_tpm.spf", sep="\t"),
        "abundance_mag": pd.read_csv(SAMPLE / "abundance.tsv", sep="\t"),
        "quality": pd.read_csv(SAMPLE / "quality_report.tsv", sep="\t"),
        "taxonomy": pd.read_csv(SAMPLE / "mag_taxonomy_labels.tsv", sep="\t", header=None,
                                 names=["Genome", "Taxonomy"]),
        "keystone": pd.read_csv(SAMPLE / "keystone_species.txt", sep="\t"),
        "ko_long": pd.read_csv(SAMPLE / "kegg_target_only.tsv", sep="\t"),
        "gephi_nodes": pd.read_csv(SAMPLE / "gephi_nodes.csv"),
        "gephi_edges": pd.read_csv(SAMPLE / "gephi_edges.csv"),
        "hypothesis_yaml": SAMPLE / "sample_hypothesis.yaml",
    }


def load_liu(subsample: tuple[int, int] | None = None):
    """Load Liu 2023 EnvMeta inputs, optionally subsample to (n_mags, n_samples).

    LEfSe / log2FC / network / RDA-with-multiple-env are not loadable from Liu
    (single Group + missing files). Those are skipped at runtime.
    """
    abundance_mag = pd.read_csv(LIU_BASE / "abundance.tsv", sep="\t")
    metadata = pd.read_csv(LIU_BASE / "metadata.tsv", sep="\t")
    env = pd.read_csv(LIU_BASE / "env_factors.tsv", sep="\t")
    ko_long = pd.read_csv(LIU_BASE / "kegg_target_only.tsv", sep="\t")
    taxonomy = pd.read_csv(LIU_BASE / "mag_taxonomy_labels.tsv", sep="\t", header=None,
                            names=["Genome", "Taxonomy"])
    quality = pd.read_csv(LIU_BASE / "quality_report.tsv", sep="\t")
    alpha = pd.read_csv(LIU_RESHAPE / "alpha.tsv", sep="\t")
    distance = pd.read_csv(LIU_RESHAPE / "distance_bray.tsv", sep="\t")
    ko_wide = pd.read_csv(LIU_RESHAPE / "ko_wide.tsv", sep="\t")
    phylum = pd.read_csv(LIU_RESHAPE / "phylum_abundance.tsv", sep="\t")

    if subsample:
        n_mag_target, n_sample_target = subsample
        # Subsample MAGs (take first N by abundance sum desc to keep top contributors)
        sample_cols = [c for c in abundance_mag.columns if c != "Genome"]
        mag_sums = abundance_mag[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
        top_mag_idx = mag_sums.nlargest(min(n_mag_target, len(abundance_mag))).index
        abundance_mag = abundance_mag.loc[top_mag_idx].reset_index(drop=True)
        keep_genomes = set(abundance_mag["Genome"])
        ko_long = ko_long[ko_long["MAG"].isin(keep_genomes)].reset_index(drop=True)
        taxonomy = taxonomy[taxonomy["Genome"].isin(keep_genomes)].reset_index(drop=True)
        quality = quality[quality["Name"].isin(keep_genomes)].reset_index(drop=True)

        # Subsample samples
        keep_samples = sample_cols[:min(n_sample_target, len(sample_cols))]
        abundance_mag = abundance_mag[["Genome"] + keep_samples]
        env_keep = env[env["SampleID"].isin(keep_samples)].reset_index(drop=True)
        env = env_keep
        metadata = metadata[metadata["SampleID"].isin(keep_samples)].reset_index(drop=True)
        alpha = alpha[alpha["SampleID"].isin(keep_samples)].reset_index(drop=True)
        # distance: filter rows + cols
        distance = distance[distance["SampleID"].isin(keep_samples)].reset_index(drop=True)
        keep_cols = ["SampleID"] + [c for c in distance.columns[1:] if c in keep_samples]
        distance = distance[keep_cols]
        # ko_wide: keep KO + sample columns
        ko_keep = ["KO"] + [c for c in ko_wide.columns[1:] if c in keep_samples]
        ko_wide = ko_wide[ko_keep]
        phy_keep = ["Taxonomy"] + [c for c in phylum.columns[1:] if c in keep_samples]
        phylum = phylum[phy_keep]

    return {
        "abundance_taxonomy": phylum,
        "metadata": metadata,
        "alpha": alpha,
        "distance": distance,
        "env": env,
        "ko_wide": ko_wide,
        "abundance_mag": abundance_mag,
        "quality": quality,
        "taxonomy": taxonomy,
        "keystone": None,
        "ko_long": ko_long,
        "gephi_nodes": None,
        "gephi_edges": None,
        "hypothesis_yaml": ROOT / "paper" / "benchmarks" / "external" / "liu_2023_coldseep"
                            / "liu2023_hypothesis.yaml",
    }


# ──────────────────────────────────────────────────────────────────
# Figure runners (each is `def run(D)` returning AnalysisResult or score)
# ──────────────────────────────────────────────────────────────────

FIGURES = {}


def figure(name):
    def _wrap(fn):
        FIGURES[name] = fn
        return fn
    return _wrap


@figure("stackplot")
def _f_stackplot(D):
    return stackplot.analyze(D["abundance_taxonomy"], D["metadata"])


@figure("alpha_boxplot")
def _f_alpha(D):
    return alpha_boxplot.analyze(D["alpha"], D["metadata"])


@figure("pcoa")
def _f_pcoa(D):
    return pcoa.analyze(D["distance"], D["metadata"])


@figure("rda")
def _f_rda(D):
    return rda.analyze(D["abundance_taxonomy"], D["env"], D["metadata"])


@figure("lefse")
def _f_lefse(D):
    return lefse.analyze(D["abundance_taxonomy"], D["metadata"])


@figure("gene_heatmap")
def _f_gene_heatmap(D):
    return gene_heatmap.analyze(D["ko_wide"], D["metadata"])


@figure("log2fc")
def _f_log2fc(D):
    groups = sorted(D["metadata"]["Group"].unique())
    if len(groups) < 2:
        raise RuntimeError(f"log2fc needs >= 2 groups (got {groups})")
    return log2fc.analyze(D["ko_wide"], D["metadata"],
                          params={"group_a": groups[0], "group_b": groups[1]})


@figure("mag_quality")
def _f_mag_quality(D):
    return mag_quality.analyze(D["quality"], taxonomy_df=D["taxonomy"],
                                keystone_df=D["keystone"])


@figure("mag_heatmap")
def _f_mag_heatmap(D):
    return mag_heatmap.analyze(D["abundance_mag"], taxonomy_df=D["taxonomy"],
                                keystone_df=D["keystone"], metadata_df=D["metadata"])


@figure("pathway")
def _f_pathway(D):
    return pathway.analyze(D["ko_long"], taxonomy_df=D["taxonomy"],
                            keystone_df=D["keystone"], abundance_df=D["abundance_mag"])


@figure("gene_profile")
def _f_gene_profile(D):
    return gene_profile.analyze(D["ko_long"], taxonomy_df=D["taxonomy"],
                                 keystone_df=D["keystone"], abundance_df=D["abundance_mag"])


@figure("network")
def _f_network(D):
    if D["gephi_nodes"] is None or D["gephi_edges"] is None:
        raise RuntimeError("network requires gephi_nodes + gephi_edges (Liu lacks)")
    return network.analyze(D["gephi_nodes"], D["gephi_edges"],
                            taxonomy_df=D["taxonomy"], keystone_df=D["keystone"])


@figure("cycle_diagram")
def _f_cycle(D):
    return cycle_diagram.analyze(
        D["ko_long"], taxonomy_df=D["taxonomy"], keystone_df=D["keystone"],
        abundance_df=D["abundance_mag"], env_df=D["env"], metadata_df=D["metadata"],
    )


@figure("hypothesis")
def _f_hypothesis(D):
    res = cycle_diagram.analyze(
        D["ko_long"], taxonomy_df=D["taxonomy"], keystone_df=D["keystone"],
        abundance_df=D["abundance_mag"], env_df=D["env"], metadata_df=D["metadata"],
    )
    h = hyp_mod.load_hypothesis(D["hypothesis_yaml"])
    return hyp_mod.score(h, res.data, run_null=True, null_n=999, run_sensitivity=True)


# ──────────────────────────────────────────────────────────────────
# Driver
# ──────────────────────────────────────────────────────────────────

def run_one(name, D, repeats: int = 3) -> dict:
    """Run figure `name` `repeats` times, return median wall-time + peak mem."""
    fn = FIGURES[name]
    times = []
    mems = []
    last_err = None
    for i in range(repeats):
        try:
            with MemSampler() as ms, timed() as elapsed:
                fn(D)
                t = elapsed()
            times.append(t)
            mems.append(ms.delta_mb)
        except Exception as e:
            last_err = e
            tb = traceback.format_exc(limit=3)
            return {"figure": name, "status": "fail",
                    "error": f"{type(e).__name__}: {e}",
                    "traceback": tb,
                    "wall_s_median": float("nan"),
                    "peak_mem_mb_median": float("nan"),
                    "n_success": 0, "n_repeats": repeats}
    return {"figure": name, "status": "ok",
            "wall_s_median": float(np.median(times)),
            "wall_s_min": float(np.min(times)),
            "wall_s_max": float(np.max(times)),
            "peak_mem_mb_median": float(np.median(mems)),
            "peak_mem_mb_max": float(np.max(mems)),
            "n_success": len(times), "n_repeats": repeats}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True, choices=["sample", "liu"])
    p.add_argument("--subsample", type=str, default=None,
                    help="N_mags,N_samples (for liu only)")
    p.add_argument("--figures", type=str, default="all",
                    help="comma-separated figure names or 'all'")
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--out", type=str, default=None)
    args = p.parse_args(argv)

    sub = None
    sub_label = ""
    if args.subsample:
        a, b = args.subsample.split(",")
        sub = (int(a), int(b))
        sub_label = f"_sub{a}m{b}s"

    print(f"[load] dataset={args.dataset} subsample={sub}")
    if args.dataset == "sample":
        D = load_sample()
    else:
        D = load_liu(subsample=sub)

    n_mags = len(D["abundance_mag"]) if D["abundance_mag"] is not None else 0
    sample_cols = [c for c in D["abundance_mag"].columns if c != "Genome"]
    n_samples = len(sample_cols)
    n_groups = D["metadata"]["Group"].nunique() if D["metadata"] is not None else 0
    print(f"  n_mags={n_mags} n_samples={n_samples} n_groups={n_groups}")

    if args.figures == "all":
        target_figs = list(FIGURES.keys())
    else:
        target_figs = [f.strip() for f in args.figures.split(",")]

    results = []
    for fname in target_figs:
        if fname not in FIGURES:
            print(f"  [skip] {fname} (not a figure name)")
            continue
        print(f"  [run]  {fname:<16}", end="", flush=True)
        r = run_one(fname, D, repeats=args.repeats)
        r["dataset"] = args.dataset + sub_label
        r["n_mags"] = n_mags
        r["n_samples"] = n_samples
        r["n_groups"] = n_groups
        if r["status"] == "ok":
            print(f" {r['wall_s_median']:7.2f} s  peak_dRSS={r['peak_mem_mb_median']:6.1f} MB")
        else:
            print(f" FAIL ({r['error'][:60]})")
        results.append(r)

    out_dir = ROOT / "paper" / "benchmarks" / "performance" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = Path(args.out) if args.out else out_dir / f"{args.dataset}{sub_label}_runtime.tsv"
    df = pd.DataFrame(results)
    cols = ["dataset", "figure", "n_mags", "n_samples", "n_groups",
            "wall_s_median", "wall_s_min", "wall_s_max",
            "peak_mem_mb_median", "peak_mem_mb_max",
            "n_success", "n_repeats", "status", "error"]
    df = df.reindex(columns=[c for c in cols if c in df.columns])
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\n[done] -> {out_path}  ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
