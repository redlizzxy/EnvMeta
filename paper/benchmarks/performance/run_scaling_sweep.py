"""Run subsample scaling sweep on Liu 2023.

Tests cycle_diagram + mag_heatmap + pathway across a 3 MAG x 3 sample grid:
- N_MAG: 200 / 500 / 1000
- N_sample: 30 / 60 / 87
= 9 cells x 3 figures = 27 runs (+ 2 repeats each = ~80 total run-equivalents)

Output: paper/benchmarks/performance/results/liu_sweep_*.tsv files combined.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HARNESS = ROOT / "paper" / "benchmarks" / "performance" / "bench_harness.py"
PYTHON = sys.executable

# 9-cell sweep
N_MAG_GRID = [200, 500, 1000]
N_SAMPLE_GRID = [30, 60, 87]
FIGURES = "cycle_diagram,mag_heatmap,pathway"
REPEATS = 2  # keep tight; sweep is exploratory


def main() -> int:
    out_dir = ROOT / "paper" / "benchmarks" / "performance" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[sweep] {len(N_MAG_GRID)}x{len(N_SAMPLE_GRID)} = "
          f"{len(N_MAG_GRID)*len(N_SAMPLE_GRID)} cells x {FIGURES.count(',')+1} figures")
    cells = [(m, s) for m in N_MAG_GRID for s in N_SAMPLE_GRID]
    for i, (m, s) in enumerate(cells, 1):
        print(f"\n[cell {i}/{len(cells)}] N_MAG={m}, N_sample={s}")
        out = out_dir / f"liu_sweep_m{m}s{s}.tsv"
        cmd = [
            PYTHON, str(HARNESS),
            "--dataset", "liu",
            "--subsample", f"{m},{s}",
            "--figures", FIGURES,
            "--repeats", str(REPEATS),
            "--out", str(out),
        ]
        rc = subprocess.call(cmd)
        if rc != 0:
            print(f"  WARNING cell {i} returned rc={rc}, continuing")
    print("\n[sweep done]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
