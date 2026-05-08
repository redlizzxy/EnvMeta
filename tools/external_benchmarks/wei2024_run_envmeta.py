"""
跑 EnvMeta 5-6 张图 on Wei et al. 2024 reshape 数据。

输入：paper/benchmarks/external/wei_2024_paddy/input_data_local/（reshape 输出）
输出：paper/benchmarks/external/wei_2024_paddy/envmeta_outputs/

用法:
    python tools/external_benchmarks/wei2024_run_envmeta.py

调 EnvMeta analysis API（matplotlib non-interactive backend），无需启 streamlit。
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

DATA_DIR = ROOT / "paper" / "benchmarks" / "external" / "wei_2024_paddy" / "input_data_local"
OUT_DIR = ROOT / "paper" / "benchmarks" / "external" / "wei_2024_paddy" / "envmeta_outputs"


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
    out_pdf = OUT_DIR / f"{name}.pdf"
    out_png = OUT_DIR / f"{name}.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"  → {out_pdf.name}")


def run_cycle_diagram(metadata, env, abund, ko_long, tax):
    print("[1/5] Cycle diagram (4 元素 × 18 通路 自动推断)")
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
    print("[2/5] MAG quality scatter")
    from envmeta.analysis import mag_quality

    t0 = time.time()
    result = mag_quality.analyze(
        quality_df=qual,
        taxonomy_df=tax,
    )
    dt = time.time() - t0
    save_figure(result.figure, "fig2_mag_quality")
    print(f"    runtime: {dt:.2f} s")
    return result


def _build_ko_per_sample(abund: pd.DataFrame, ko_long: pd.DataFrame) -> pd.DataFrame:
    """KO long → KO × sample wide（去掉 'ko:' 前缀，列名 KEGG_ko）"""
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


def run_gene_heatmap(metadata, abund, ko_long):
    print("[3/5] Gene heatmap (精选 KO × group means)")
    from envmeta.analysis import gene_heatmap
    t0 = time.time()
    ko_per_sample = _build_ko_per_sample(abund, ko_long)
    try:
        result = gene_heatmap.analyze(
            ko_abundance_df=ko_per_sample,
            metadata_df=metadata,
        )
        dt = time.time() - t0
        save_figure(result.figure, "fig3_gene_heatmap")
        print(f"    runtime: {dt:.2f} s")
        return result
    except Exception as e:
        print(f"    [SKIP] {type(e).__name__}: {e}")
        return None


def run_log2fc(metadata, abund, ko_long):
    print("[4/5] log2FC (AsContam vs NoContam)")
    from envmeta.analysis import log2fc
    t0 = time.time()
    ko_per_sample = _build_ko_per_sample(abund, ko_long)
    try:
        result = log2fc.analyze(
            ko_abundance_df=ko_per_sample,
            metadata_df=metadata,
            params={"group_a": "AsContam", "group_b": "NoContam"},
        )
        dt = time.time() - t0
        save_figure(result.figure, "fig4_log2fc")
        print(f"    runtime: {dt:.2f} s")
        return result
    except Exception as e:
        print(f"    [SKIP] {type(e).__name__}: {e}")
        return None


def run_mag_heatmap(metadata, abund, tax):
    print("[5/5] MAG abundance heatmap")
    from envmeta.analysis import mag_heatmap
    t0 = time.time()
    try:
        result = mag_heatmap.analyze(
            abundance_df=abund,
            metadata_df=metadata,
            taxonomy_df=tax,
        )
        dt = time.time() - t0
        save_figure(result.figure, "fig5_mag_heatmap")
        print(f"    runtime: {dt:.2f} s")
        return result
    except Exception as e:
        print(f"    [SKIP] {type(e).__name__}: {e}")
        return None


def run_hypothesis_score(cycle_result, yaml_path: Path):
    """[6/6] YAML 假说评分（EnvMeta 独家卖点）— 复现 Wei 2024 As-N 耦合结论"""
    print(f"[6/6] Hypothesis scoring (YAML: {yaml_path.name})")
    from envmeta.geocycle.hypothesis import load_hypothesis, score
    t0 = time.time()
    hyp = load_hypothesis(yaml_path)
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
        (f"- weight_robust: {s.weight_robust} (OAT ±20% perturbation)"
         if s.weight_robust is not None else ""),
        (f"- veto_reasons: {s.veto_reasons}"
         if s.veto_reasons else "- veto_reasons: none (no required claim failed)"),
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

    cycle_result = run_cycle_diagram(metadata, env, abund, ko_long, tax)
    run_mag_quality(qual, tax)
    run_gene_heatmap(metadata, abund, ko_long)
    run_log2fc(metadata, abund, ko_long)
    run_mag_heatmap(metadata, abund, tax)

    yaml_path = ROOT / "paper" / "benchmarks" / "external" / "wei_2024_paddy" / "wei2024_hypothesis.yaml"
    if yaml_path.exists():
        run_hypothesis_score(cycle_result, yaml_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
