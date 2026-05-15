"""
Threshold sensitivity analysis — auxiliary evidence that calibration outcomes
are robust to default-threshold choice (strong=0.75, suggestive=0.40).

Sweep strong_threshold ∈ {0.65, 0.70, 0.75, 0.80, 0.85} for each calibration
arm and stress run; record the label at each threshold. STRONG outcomes that
remain stable across this range are robust to the specific default choice.

Output:
    paper/benchmarks/external/threshold_sensitivity/
        threshold_sweep.tsv     (one row per dataset × threshold)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")

EXTERNAL = ROOT / "paper" / "benchmarks" / "external"
OUT_DIR = EXTERNAL / "threshold_sensitivity"


# All 4 calibration arms + 3 stress runs (v1)
RUNS = [
    {
        "key": "arm_a", "label": "Arm A (in-house) calibration",
        "yaml": str(ROOT / "paper" / "hypotheses" / "arsenic_steel_slag.yaml"),
        "data_loader": "sample_data",
        "compare_groups": ["CK", "A", "B"],
    },
    {
        "key": "arm_b_wei", "label": "Arm B (Wei 2024) calibration",
        "yaml": str(EXTERNAL / "wei_2024_paddy" / "wei2024_hypothesis.yaml"),
        "dataset": "wei_2024_paddy",
    },
    {
        "key": "arm_c1_liu", "label": "Arm C1 (Liu 2023) calibration",
        "yaml": str(EXTERNAL / "liu_2023_coldseep" / "liu2023_hypothesis.yaml"),
        "dataset": "liu_2023_coldseep",
    },
    {
        "key": "arm_c2a_grett", "label": "Arm C2-A (Grettenberger 2021) calibration",
        "yaml": str(EXTERNAL / "grettenberger_2021_amd_stream" / "grettenberger2021_hypothesis.yaml"),
        "dataset": "grettenberger_2021_amd_stream",
    },
    {
        "key": "arm_c2b_ayala", "label": "Arm C2-B (Ayala 2020) calibration",
        "yaml": str(EXTERNAL / "ayala_2020_pitlake" / "ayala2020_hypothesis.yaml"),
        "dataset": "ayala_2020_pitlake",
    },
    {
        "key": "stress_liu_v1", "label": "Liu 2023 stress (v1)",
        "yaml": str(EXTERNAL / "liu_2023_coldseep" / "liu2023_hypothesis_stress.yaml"),
        "dataset": "liu_2023_coldseep",
    },
    {
        "key": "stress_grett", "label": "Grettenberger 2021 stress",
        "yaml": str(EXTERNAL / "grettenberger_2021_amd_stream" / "grettenberger2021_hypothesis_stress.yaml"),
        "dataset": "grettenberger_2021_amd_stream",
    },
    {
        "key": "stress_ayala_v1", "label": "Ayala 2020 stress (v1)",
        "yaml": str(EXTERNAL / "ayala_2020_pitlake" / "ayala2020_hypothesis_stress.yaml"),
        "dataset": "ayala_2020_pitlake",
    },
]

THRESHOLDS = [0.65, 0.70, 0.75, 0.80, 0.85]


def load_inputs(data_dir: Path):
    metadata = pd.read_csv(data_dir / "metadata.tsv", sep="\t")
    env = pd.read_csv(data_dir / "env_factors.tsv", sep="\t")
    abund = pd.read_csv(data_dir / "abundance.tsv", sep="\t")
    ko_long = pd.read_csv(data_dir / "kegg_target_only.tsv", sep="\t")
    qual = pd.read_csv(data_dir / "quality_report.tsv", sep="\t")
    tax = pd.read_csv(
        data_dir / "mag_taxonomy_labels.tsv",
        sep="\t", header=None, names=["MAG", "classification"],
    )
    return metadata, env, abund, ko_long, qual, tax, None


def load_inputs_sample_data():
    dd = ROOT / "tests" / "sample_data"
    metadata = pd.read_csv(dd / "metadata.txt", sep="\t")
    env = pd.read_csv(dd / "env_factors.txt", sep="\t")
    abund = pd.read_csv(dd / "abundance.tsv", sep="\t")
    ko_long = pd.read_csv(dd / "kegg_target_only.tsv", sep="\t")
    qual = pd.read_csv(dd / "quality_report.tsv", sep="\t")
    tax = pd.read_csv(
        dd / "mag_taxonomy_labels.tsv",
        sep="\t", header=None, names=["MAG", "classification"],
    )
    keystone = pd.read_csv(dd / "keystone_species.txt", sep="\t")
    return metadata, env, abund, ko_long, qual, tax, keystone


def label_for(score: float, strong: float, suggestive: float) -> str:
    if score >= strong:
        return "strong"
    if score >= suggestive:
        return "suggestive"
    if score > 0:
        return "weak"
    return "insufficient"


def run_one(cfg: dict) -> list[dict]:
    """Score one YAML at multiple thresholds; return one row per threshold."""
    from envmeta.analysis import cycle_diagram, cycle_compare
    from envmeta.geocycle.hypothesis import load_hypothesis, score
    import yaml as yamllib

    yaml_path = Path(cfg["yaml"])
    if not yaml_path.exists():
        print(f"  [SKIP] missing yaml: {yaml_path}")
        return []

    if cfg.get("data_loader") == "sample_data":
        metadata, env, abund, ko_long, qual, tax, keystone = load_inputs_sample_data()
    else:
        data_dir = EXTERNAL / cfg["dataset"] / "input_data_local"
        if not data_dir.exists():
            print(f"  [SKIP] missing data: {data_dir}")
            return []
        metadata, env, abund, ko_long, qual, tax, keystone = load_inputs(data_dir)

    cycle_result = cycle_diagram.analyze(
        ko_annotation_df=ko_long, taxonomy_df=tax,
        abundance_df=abund, env_df=env, metadata_df=metadata,
        keystone_df=keystone,
    )
    compare_df = None
    if cfg.get("compare_groups"):
        compare_df = cycle_compare.compare_groups(
            ko_annotation_df=ko_long, taxonomy_df=tax,
            abundance_df=abund, env_df=env, metadata_df=metadata,
            keystone_df=keystone, groups=cfg["compare_groups"],
        )

    # Score once with default thresholds; we recompute label by hand at each
    # threshold from the same overall_score and required-veto state.
    with open(yaml_path, "r", encoding="utf-8") as f:
        raw = yamllib.safe_load(f)
    hyp = load_hypothesis(yaml_path)
    s = score(hyp, cycle_result.data, compare_df=compare_df,
              run_null=False, run_sensitivity=False)

    rows = []
    for thr in THRESHOLDS:
        sugg = thr - 0.35  # keep gap = 0.35 (matches default 0.75-0.40)
        if s.veto_reasons:
            label = "insufficient"
        else:
            label = label_for(s.overall_score, thr, sugg)
        rows.append({
            "run_key": cfg["key"],
            "run_label": cfg["label"],
            "overall_score": s.overall_score,
            "strong_threshold": thr,
            "suggestive_threshold": round(sugg, 2),
            "label_at_threshold": label,
            "veto_active": bool(s.veto_reasons),
            "n_satisfied": s.n_satisfied,
            "n_total": s.n_total,
            "n_skipped": s.n_skipped,
        })
    return rows


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    for cfg in RUNS:
        print(f"\n=== {cfg['label']} ===")
        rows = run_one(cfg)
        for r in rows:
            print(f"  strong={r['strong_threshold']:.2f} → label={r['label_at_threshold']}"
                  f" (overall={r['overall_score']:.3f}, veto={r['veto_active']})")
        all_rows.extend(rows)

    if not all_rows:
        print("No results.")
        return

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT_DIR / "threshold_sweep.tsv", sep="\t", index=False)

    print("\n" + "=" * 70 + "\n  THRESHOLD STABILITY MATRIX\n" + "=" * 70)
    pivot = df.pivot(index="run_label", columns="strong_threshold",
                     values="label_at_threshold")
    print(pivot.to_string())
    pivot.to_csv(OUT_DIR / "threshold_stability_matrix.tsv", sep="\t")

    print(f"\nResults written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
