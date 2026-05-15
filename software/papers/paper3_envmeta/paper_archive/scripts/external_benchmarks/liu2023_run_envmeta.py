"""
跑 EnvMeta 5 张图 + 1 假说评分 on Liu et al. 2023 reshape 数据。

输入：paper/benchmarks/external/liu_2023_coldseep/input_data_local/
输出：paper/benchmarks/external/liu_2023_coldseep/envmeta_outputs/
YAML：paper/benchmarks/external/liu_2023_coldseep/liu2023_hypothesis.yaml
       (pre-registered before this run, commit 42168da)
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

DS = "liu_2023_coldseep"
DATA_DIR = ROOT / "paper" / "benchmarks" / "external" / DS / "input_data_local"
OUT_DIR = ROOT / "paper" / "benchmarks" / "external" / DS / "envmeta_outputs"
YAML_PATH = ROOT / "paper" / "benchmarks" / "external" / DS / "liu2023_hypothesis.yaml"


def load_inputs():
    metadata = pd.read_csv(DATA_DIR / "metadata.tsv", sep="\t")
    env = pd.read_csv(DATA_DIR / "env_factors.tsv", sep="\t")
    abund = pd.read_csv(DATA_DIR / "abundance.tsv", sep="\t")
    ko_long = pd.read_csv(DATA_DIR / "kegg_target_only.tsv", sep="\t")
    qual = pd.read_csv(DATA_DIR / "quality_report.tsv", sep="\t")
    tax_raw = pd.read_csv(
        DATA_DIR / "mag_taxonomy_labels.tsv",
        sep="\t", header=None, names=["MAG", "classification"],
    )
    return metadata, env, abund, ko_long, qual, tax_raw


def save_figure(fig, name: str):
    fig.savefig(OUT_DIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{name}.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"  → {name}.pdf")


def _build_ko_per_sample(abund: pd.DataFrame, ko_long: pd.DataFrame) -> pd.DataFrame:
    ko_clean = ko_long.copy()
    ko_clean["KEGG_ko"] = ko_clean["KEGG_ko"].str.replace("ko:", "", regex=False)
    ko_pivot = ko_clean.assign(present=1).pivot_table(
        index="MAG", columns="KEGG_ko", values="present",
        aggfunc="sum", fill_value=0,
    )
    abund_idx = abund.set_index("Genome")
    common_mags = abund_idx.index.intersection(ko_pivot.index)
    mat = abund_idx.loc[common_mags].T @ ko_pivot.loc[common_mags]
    out = mat.T.reset_index()
    out.columns.name = None
    out = out.rename(columns={out.columns[0]: "KEGG_ko"})
    return out


def run_cycle(metadata, env, abund, ko_long, tax):
    print("[1/5] Cycle diagram")
    from envmeta.analysis import cycle_diagram
    t0 = time.time()
    result = cycle_diagram.analyze(
        ko_annotation_df=ko_long,
        taxonomy_df=tax,
        abundance_df=abund,
        env_df=env,
        metadata_df=metadata,
    )
    dt = time.time() - t0
    save_figure(result.figure, "fig1_cycle_diagram")
    result.stats.to_csv(OUT_DIR / "fig1_cycle_diagram_stats.tsv",
                        sep="\t", index=False)
    print(f"    runtime: {dt:.2f} s, stats rows: {len(result.stats)}")
    return result


def run_mag_quality(qual, tax):
    print("[2/5] MAG quality")
    from envmeta.analysis import mag_quality
    t0 = time.time()
    try:
        r = mag_quality.analyze(quality_df=qual, taxonomy_df=tax)
        save_figure(r.figure, "fig2_mag_quality")
        print(f"    runtime: {time.time() - t0:.2f} s")
    except Exception as e:
        print(f"    [SKIP] {type(e).__name__}: {e}")


def run_gene_heatmap(metadata, abund, ko_long):
    print("[3/5] Gene heatmap")
    from envmeta.analysis import gene_heatmap
    t0 = time.time()
    try:
        ko_per_sample = _build_ko_per_sample(abund, ko_long)
        r = gene_heatmap.analyze(
            ko_abundance_df=ko_per_sample, metadata_df=metadata,
        )
        save_figure(r.figure, "fig3_gene_heatmap")
        print(f"    runtime: {time.time() - t0:.2f} s")
    except Exception as e:
        print(f"    [SKIP] {type(e).__name__}: {e}")


def run_log2fc(metadata, abund, ko_long):
    """Liu 是单组 'All'，log2fc 跳过（无 group A vs group B）"""
    print("[4/5] log2FC (skipped: Liu 是单组 pooled 分析，无 case/control)")


def run_mag_heatmap(metadata, abund, tax):
    print("[5/5] MAG abundance heatmap")
    from envmeta.analysis import mag_heatmap
    t0 = time.time()
    try:
        r = mag_heatmap.analyze(
            abundance_df=abund, metadata_df=metadata, taxonomy_df=tax,
        )
        save_figure(r.figure, "fig5_mag_heatmap")
        print(f"    runtime: {time.time() - t0:.2f} s")
    except Exception as e:
        print(f"    [SKIP] {type(e).__name__}: {e}")


def run_hypothesis_score(cycle_result):
    print(f"[6/6] Hypothesis scoring (YAML: {YAML_PATH.name})")
    from envmeta.geocycle.hypothesis import load_hypothesis, score
    t0 = time.time()
    hyp = load_hypothesis(YAML_PATH)
    s = score(hyp, cycle_result.data, run_null=True, null_n=999, run_sensitivity=True)
    dt = time.time() - t0

    out_md = OUT_DIR / "fig6_hypothesis_score.md"
    lines = [
        f"# {hyp.name} — Hypothesis Score",
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
    print(f"    overall_score={s.overall_score:.3f} label={s.label} "
          f"null_p={s.null_p} satisfied={s.n_satisfied}/{s.n_total}")
    print(f"    runtime: {dt:.2f} s")
    return s


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Loading inputs from: {DATA_DIR}")
    metadata, env, abund, ko_long, qual, tax = load_inputs()
    print(f"  metadata: {len(metadata)}, env: {len(env)}, abund: {abund.shape}, "
          f"ko_long: {len(ko_long)}, qual: {len(qual)}, tax: {len(tax)}")
    print(f"\nWriting outputs to: {OUT_DIR}\n")

    cycle_result = run_cycle(metadata, env, abund, ko_long, tax)
    run_mag_quality(qual, tax)
    run_gene_heatmap(metadata, abund, ko_long)
    run_log2fc(metadata, abund, ko_long)
    run_mag_heatmap(metadata, abund, tax)

    if YAML_PATH.exists():
        run_hypothesis_score(cycle_result)
    else:
        print(f"[6/6] YAML not found at {YAML_PATH}, skip scoring")

    print("\nDone.")


if __name__ == "__main__":
    main()
