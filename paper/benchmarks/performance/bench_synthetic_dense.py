"""Synthetic dense-annotation scaling test.

Liu's published kegg_target_only.tsv is pre-filtered to 8 As-cycle KOs (sparse).
This script augments Liu's MAG-KO annotation by randomly assigning each MAG a
realistic-density set of KOs from EnvMeta's 57-KO knowledge base, simulating
what a full KofamScan/DRAM/METABOLIC-annotated 1000+ MAG dataset would look
like. Also synthesizes 4 env factors (vs Liu's 1) to mirror sample_data env
breadth (4 factors). Then runs cycle_diagram across the 9-cell sweep.

Output: paper/benchmarks/performance/results/liu_synth_dense_*.tsv
"""
from __future__ import annotations

import gc
import sys
import threading
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from envmeta.analysis import cycle_diagram, mag_heatmap, pathway
from envmeta.geocycle import knowledge_base as kb_mod

LIU_BASE = ROOT / "paper" / "benchmarks" / "external" / "liu_2023_coldseep" / "input_data_local"
RESHAPE = ROOT / "paper" / "benchmarks" / "performance" / "liu_reshaped"
OUT = ROOT / "paper" / "benchmarks" / "performance" / "results"
OUT.mkdir(parents=True, exist_ok=True)

KOS_PER_MAG = 25  # realistic mean KO/MAG hits for 4-element KB at typical density
RNG_SEED = 42


def get_full_kb_kos() -> list[str]:
    """Get the full set of 57 KOs from EnvMeta's KB via flat_ko_map()."""
    flat = kb_mod.flat_ko_map()  # ko -> {element, pathway, ...}
    return sorted(flat.keys())


def build_synthetic_ko_long(abundance: pd.DataFrame, seed: int = RNG_SEED) -> pd.DataFrame:
    """Generate dense MAG-KO annotation: each MAG gets ~KOS_PER_MAG random KOs from KB."""
    rng = np.random.default_rng(seed)
    all_kos = get_full_kb_kos()
    rows = []
    for genome in abundance["Genome"]:
        n_kos = rng.integers(KOS_PER_MAG - 5, KOS_PER_MAG + 5)
        chosen = rng.choice(all_kos, size=min(n_kos, len(all_kos)), replace=False)
        for ko in chosen:
            rows.append({
                "MAG": genome,
                "Gene_ID": f"{genome}_synth_{ko}",
                "KEGG_ko": f"ko:{ko}",
                "Description": "synthetic dense annotation",
            })
    return pd.DataFrame(rows)


def build_synthetic_env(metadata: pd.DataFrame, seed: int = RNG_SEED) -> pd.DataFrame:
    """Generate 4 numeric env factors (correlated/uncorrelated mix) to match sample_data breadth."""
    rng = np.random.default_rng(seed)
    n = len(metadata)
    return pd.DataFrame({
        "SampleID": metadata["SampleID"].values,
        "Group": metadata["Group"].values,
        "pH":       rng.uniform(5.0, 8.5, n),
        "Eh":       rng.uniform(-200, 300, n),
        "TOC":      rng.uniform(0.5, 8.0, n),
        "Total_S":  rng.uniform(50, 800, n),
    })


# --- runtime + memory measurement ----------------------------------

class MemSampler:
    def __init__(self, interval: float = 0.05):
        self.interval = interval
        self.proc = psutil.Process()
        self.peak = 0
        self.baseline = 0
        self._stop = threading.Event()
        self._thread = None

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


def run_cell(n_mag: int, n_sample: int, abundance, taxonomy, quality,
             metadata, repeats: int = 2) -> list[dict]:
    """Subsample, build synthetic dense KO + 4 env, run 3 figures."""
    sample_cols = [c for c in abundance.columns if c != "Genome"]
    keep_samples = sample_cols[:n_sample]
    sample_sums = abundance[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
    top_idx = sample_sums.nlargest(n_mag).index
    ab_sub = abundance.loc[top_idx, ["Genome"] + keep_samples].reset_index(drop=True)

    ko_synth = build_synthetic_ko_long(ab_sub)
    meta_sub = metadata[metadata["SampleID"].isin(keep_samples)].reset_index(drop=True)
    env_synth = build_synthetic_env(meta_sub)
    keep_genomes = set(ab_sub["Genome"])
    tax_sub = taxonomy[taxonomy["Genome"].isin(keep_genomes)].reset_index(drop=True)
    qual_sub = quality[quality["Name"].isin(keep_genomes)].reset_index(drop=True)

    print(f"  [cell] N_MAG={n_mag}, N_sample={n_sample}, "
          f"synth_KOs={len(ko_synth)}, env_cols=4")

    runs = []
    figures = {
        "cycle_diagram": lambda: cycle_diagram.analyze(
            ko_synth, taxonomy_df=tax_sub, abundance_df=ab_sub,
            env_df=env_synth, metadata_df=meta_sub),
        "mag_heatmap": lambda: mag_heatmap.analyze(
            ab_sub, taxonomy_df=tax_sub, metadata_df=meta_sub),
        "pathway": lambda: pathway.analyze(
            ko_synth, taxonomy_df=tax_sub, abundance_df=ab_sub),
    }
    for fname, fn in figures.items():
        times = []
        mems = []
        for _ in range(repeats):
            try:
                with MemSampler() as ms:
                    t0 = time.perf_counter()
                    fn()
                    t = time.perf_counter() - t0
                times.append(t)
                mems.append(ms.delta_mb)
            except Exception as e:
                print(f"    [fail] {fname}: {type(e).__name__}: {e}")
                break
        if times:
            r = {
                "dataset": f"liu_synth_dense_m{n_mag}s{n_sample}",
                "figure": fname,
                "n_mags": n_mag,
                "n_samples": n_sample,
                "n_groups": meta_sub["Group"].nunique(),
                "wall_s_median": float(np.median(times)),
                "wall_s_min": float(np.min(times)),
                "wall_s_max": float(np.max(times)),
                "peak_mem_mb_median": float(np.median(mems)),
                "peak_mem_mb_max": float(np.max(mems)),
                "n_success": len(times),
                "n_repeats": repeats,
                "status": "ok",
            }
            print(f"    {fname:<14} wall={r['wall_s_median']:7.2f}s "
                  f"peak_dRSS={r['peak_mem_mb_max']:6.1f}MB")
            runs.append(r)
    return runs


def main() -> int:
    print(f"[load] Liu 2023 base data...")
    abundance = pd.read_csv(LIU_BASE / "abundance.tsv", sep="\t")
    taxonomy = pd.read_csv(LIU_BASE / "mag_taxonomy_labels.tsv", sep="\t",
                            header=None, names=["Genome", "Taxonomy"])
    quality = pd.read_csv(LIU_BASE / "quality_report.tsv", sep="\t")
    metadata = pd.read_csv(LIU_BASE / "metadata.tsv", sep="\t")
    print(f"  {len(abundance)} MAGs x {abundance.shape[1]-1} samples available")

    cells = [(200, 30), (200, 87), (500, 30), (500, 60), (500, 87),
             (1000, 30), (1000, 60), (1000, 87)]
    print(f"[sweep] {len(cells)} synthetic-dense cells x 3 figures x 2 repeats")

    all_results = []
    for m, s in cells:
        all_results.extend(run_cell(m, s, abundance, taxonomy, quality, metadata))

    df = pd.DataFrame(all_results)
    out_path = OUT / "liu_synth_dense_runtime.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\n[done] -> {out_path}  ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
