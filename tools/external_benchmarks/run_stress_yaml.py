"""
通用 stress YAML runner — 跑任一 external dataset 的 stress hypothesis YAML。

用法:
    python tools/external_benchmarks/run_stress_yaml.py \\
        --dataset liu_2023_coldseep \\
        --yaml liu2023_hypothesis_stress.yaml

或同时跑多个:
    python tools/external_benchmarks/run_stress_yaml.py --all

输入: paper/benchmarks/external/{dataset}/input_data_local/
YAML: paper/benchmarks/external/{dataset}/{yaml}
输出: paper/benchmarks/external/{dataset}/envmeta_outputs/{yaml_stem}_score.md

Stress YAML 跑分逻辑：
  1. 读 input data (metadata/env/abund/ko_long/qual/tax)
  2. 跑 cycle_diagram.analyze 得 CycleData
  3. load_hypothesis(stress_yaml)
  4. score(hyp, cycle_data) — 不调阈值，不传 compare_df
  5. 输出 fig6_<yaml_stem>_score.md（与 calibration fig6_hypothesis_score.md 并存）

Pre-registration 纪律：
  - 不修改 stress YAML 解决意外结果
  - 不基于结果改 expected_label
  - 跑完后写 stress_test_results.md 做 diff
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import pandas as pd

EXTERNAL = ROOT / "paper" / "benchmarks" / "external"

# 默认 stress YAML 配置（--all 时遍历）
STRESS_TARGETS = [
    ("liu_2023_coldseep", "liu2023_hypothesis_stress.yaml"),
    ("grettenberger_2021_amd_stream", "grettenberger2021_hypothesis_stress.yaml"),
    # ayala 等 GhostKOALA 数据后再加
]


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
    return metadata, env, abund, ko_long, qual, tax


def run_cycle_and_score(dataset: str, yaml_name: str) -> dict:
    data_dir = EXTERNAL / dataset / "input_data_local"
    out_dir = EXTERNAL / dataset / "envmeta_outputs"
    yaml_path = EXTERNAL / dataset / yaml_name
    out_dir.mkdir(parents=True, exist_ok=True)

    if not yaml_path.exists():
        return {"dataset": dataset, "yaml": yaml_name, "error": f"YAML not found: {yaml_path}"}
    if not data_dir.exists():
        return {"dataset": dataset, "yaml": yaml_name, "error": f"Data dir missing: {data_dir}"}

    print(f"\n{'=' * 70}")
    print(f"  STRESS RUN -- {dataset} / {yaml_name}")
    print('=' * 70)

    print(f"Loading inputs from: {data_dir}")
    metadata, env, abund, ko_long, qual, tax = load_inputs(data_dir)
    print(f"  metadata: {len(metadata)}, env: {len(env)}, abund: {abund.shape}, "
          f"ko_long: {len(ko_long)}, qual: {len(qual)}, tax: {len(tax)}")

    # cycle diagram
    print("\n[1/2] Cycle diagram (recompute for stress YAML)")
    from envmeta.analysis import cycle_diagram
    t0 = time.time()
    cycle_result = cycle_diagram.analyze(
        ko_annotation_df=ko_long, taxonomy_df=tax,
        abundance_df=abund, env_df=env, metadata_df=metadata,
    )
    print(f"    runtime: {time.time() - t0:.2f} s")

    # hypothesis score (stress YAML)
    print(f"\n[2/2] Stress hypothesis scoring (YAML: {yaml_path.name})")
    from envmeta.geocycle.hypothesis import load_hypothesis, score
    t0 = time.time()
    hyp = load_hypothesis(yaml_path)
    s = score(hyp, cycle_result.data, run_null=True, null_n=999, run_sensitivity=True)
    dt = time.time() - t0

    yaml_stem = yaml_path.stem
    out_md = out_dir / f"fig6_{yaml_stem}_score.md"
    lines = [
        f"# {hyp.name} — Stress Hypothesis Score",
        "",
        f"- overall_score: **{s.overall_score:.3f}**",
        f"- label: **{s.label.upper()}**",
        f"- n_satisfied / n_total / n_skipped: {s.n_satisfied} / {s.n_total} / {s.n_skipped}",
        (f"- null_p: {s.null_p:.4f} (n={s.null_p_samples} permutations)"
         if s.null_p is not None else "- null_p: n/a"),
        (f"- weight_robust: {s.weight_robust} (OAT ±20%)"
         if s.weight_robust is not None else ""),
        (f"- veto_reasons: {s.veto_reasons}"
         if s.veto_reasons else "- veto_reasons: none"),
        "",
        "## Per-claim results",
        "",
        "| ID | Type | Status | Score | Weight |",
        "|---|---|---|---|---|",
    ]
    for cr in s.claim_results:
        lines.append(
            f"| {cr.claim_id} | {cr.claim_type} | {cr.status} | "
            f"{cr.score:.2f} | {cr.weight} |"
        )
    lines += ["", "## Evidence per claim", ""]
    for cr in s.claim_results:
        lines.append(f"### {cr.claim_id}")
        lines.append("")
        lines.append(f"- status={cr.status} score={cr.score:.2f} weight={cr.weight}")
        if cr.evidence:
            lines.append(f"- evidence: `{cr.evidence}`")
        if cr.explanation:
            lines.append(f"- evaluator note: {cr.explanation}")
        lines.append("")
    out_md.write_text("\n".join(lines), encoding="utf-8")

    print(f"  → {out_md.name}")
    print(f"    overall_score={s.overall_score:.3f} label={s.label}")
    print(f"    null_p={s.null_p} weight_robust={s.weight_robust}")
    print(f"    satisfied={s.n_satisfied}/{s.n_total} skipped={s.n_skipped}")
    print(f"    runtime: {dt:.2f} s")

    return {
        "dataset": dataset,
        "yaml": yaml_name,
        "overall_score": s.overall_score,
        "label": s.label,
        "n_satisfied": s.n_satisfied,
        "n_total": s.n_total,
        "n_skipped": s.n_skipped,
        "null_p": s.null_p,
        "weight_robust": s.weight_robust,
        "claim_results": [
            {"id": cr.claim_id, "status": cr.status, "score": cr.score,
             "weight": cr.weight, "type": cr.claim_type}
            for cr in s.claim_results
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", help="dataset folder name (under paper/benchmarks/external/)")
    p.add_argument("--yaml", help="YAML filename (relative to dataset folder)")
    p.add_argument("--all", action="store_true",
                   help="Run all configured stress targets")
    args = p.parse_args()

    if args.all:
        results = []
        for ds, yml in STRESS_TARGETS:
            results.append(run_cycle_and_score(ds, yml))
        print("\n" + "=" * 70)
        print("  SUMMARY")
        print("=" * 70)
        for r in results:
            if "error" in r:
                print(f"  [ERROR] {r['dataset']}: {r['error']}")
            else:
                print(f"  - {r['dataset']}: overall={r['overall_score']:.3f} "
                      f"label={r['label']} ({r['n_satisfied']}/{r['n_total']} satisfied)")
    elif args.dataset and args.yaml:
        run_cycle_and_score(args.dataset, args.yaml)
    else:
        p.error("Need either --all or both --dataset and --yaml")


if __name__ == "__main__":
    main()
