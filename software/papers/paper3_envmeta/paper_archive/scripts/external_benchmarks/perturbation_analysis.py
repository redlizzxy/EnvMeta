"""
Perturbation analysis — auxiliary evidence for Paper 3 calibration claims.

研究设计（mock review v0.9.2 Major #1 应对）：
    把 4 个 STRONG calibration YAML 的 pathway_active / keystone_in_pathway /
    env_correlation 三类 claim 的 `pathway` 参数随机替换为同一 KB element
    内的其他通路，重新评分。N=20 次置换 / dataset。

    若 STRONG outcome 依然机械产生 → KEGG-coverage 主导，calibration 失去意义。
    若大多数置换降级到 weak / suggestive → 原 STRONG 不是机械产物，对作者
    选择的 target pathway 敏感。

输出：
    paper/benchmarks/external/perturbation/
        ├── perturbation_results.tsv     # 每次 run 一行
        ├── perturbation_summary.tsv     # 每个 dataset 一行汇总
        └── perturbation_curve.{pdf,png,svg}  # original vs perturbed 分布
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from copy import deepcopy
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXTERNAL = ROOT / "paper" / "benchmarks" / "external"
OUT_DIR = EXTERNAL / "perturbation"


# 4 校准 dataset 配置
# `restrict_claim_types`：若设，仅对这些 claim type 做 perturbation；其他 claim
# 保持原 YAML 不动（Arm A 的 partial perturbation 用，遵循 mock review v0.9.3
# Major #1 reviewer 建议：仅扰动 3 个 pathway_active claim）。
# `data_loader`: optional callable returning (metadata, env, abund, ko_long, qual, tax)；
# 默认走 load_inputs（input_data_local/{file}.tsv 标准命名）。
# `compare_groups`: optional list — 若给出，对每次评分都计算 compare_df 用于 group_contrast claim。
DATASETS = [
    {
        "key": "arm_a_arsenic_steel_slag",
        "label": "Arm A (in-house)",
        "yaml": str(ROOT / "paper" / "hypotheses" / "arsenic_steel_slag.yaml"),
        "topic_element": "arsenic",
        "data_loader": "sample_data",   # 标记，runner 用专用 loader
        "compare_groups": ["CK", "A", "B"],
        "restrict_claim_types": ["pathway_active"],   # partial perturbation
    },
    {
        "key": "liu_2023_coldseep",
        "label": "Liu 2023",
        "yaml": "liu2023_hypothesis.yaml",
        "topic_element": "arsenic",
    },
    {
        "key": "grettenberger_2021_amd_stream",
        "label": "Grettenberger 2021",
        "yaml": "grettenberger2021_hypothesis.yaml",
        "topic_element": "sulfur",
    },
    {
        "key": "ayala_2020_pitlake",
        "label": "Ayala 2020",
        "yaml": "ayala2020_hypothesis.yaml",
        "topic_element": "sulfur",
    },
]

# Pathway 名 → element 的反向映射（基于 KB elements.json v2.0）
PATHWAY_TO_ELEMENT = {
    # arsenic
    "Arsenate reduction": "arsenic",
    "Arsenite oxidation": "arsenic",
    "Resp. arsenate red.": "arsenic",
    "As transport/detox": "arsenic",
    "As methylation": "arsenic",
    "As regulation": "arsenic",
    # nitrogen
    "Nitrate reduction": "nitrogen",
    "Nitrite reduction": "nitrogen",
    "NO reduction": "nitrogen",
    "N2O reduction": "nitrogen",
    "Ammonia oxidation": "nitrogen",
    "N fixation": "nitrogen",
    # sulfur
    "Assim. sulfate red.": "sulfur",
    "Dissim. sulfate red.": "sulfur",
    "Sulfide oxidation": "sulfur",
    "Thiosulfate metab.": "sulfur",
    # iron
    "Fe transport": "iron",
    "Fe uptake regulation": "iron",
}

ELEMENT_TO_PATHWAYS = {}
for pw, el in PATHWAY_TO_ELEMENT.items():
    ELEMENT_TO_PATHWAYS.setdefault(el, []).append(pw)


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


def load_inputs_sample_data():
    """tests/sample_data 加载器 — 文件名 .txt 而非 .tsv（Arm A 作者数据精简版）。

    返回 7 元组（多一个 keystone_df），与 load_inputs 区分。
    """
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


def perturb_yaml(
    orig_dict: dict, rng: random.Random, mode: str = "within_element",
    restrict_claim_types: list[str] | None = None,
) -> tuple[dict, list[dict]]:
    """随机替换所有含 `pathway` 参数的 claim 的 pathway 名。

    mode = "within_element"：从同 element 内随机选另一通路（保守测试，
        审稿人原始建议）。
    mode = "cross_element"：从**不同** element 中随机选通路（强阴性对照，
        预期数据不包含相应通路 → 大多 skipped → 退化到 0 分）。
    restrict_claim_types: 若给出，仅扰动这些 claim type；其他保持不变
        （Arm A partial perturbation 用，仅扰 3 个 pathway_active claim，
        保留 coupling_possible / env_correlation / group_contrast 不动）。

    Returns:
        (新 YAML dict, 替换记录 list[{claim_id, original, perturbed, ...}])
    """
    new = deepcopy(orig_dict)
    log: list[dict] = []
    for claim in new.get("claims", []):
        if restrict_claim_types is not None and claim.get("type") not in restrict_claim_types:
            continue
        params = claim.get("params") or {}
        pw = params.get("pathway")
        if not pw or pw not in PATHWAY_TO_ELEMENT:
            continue
        element = PATHWAY_TO_ELEMENT[pw]
        if mode == "within_element":
            candidates = [p for p in ELEMENT_TO_PATHWAYS[element] if p != pw]
        elif mode == "cross_element":
            candidates = [
                p for p, el in PATHWAY_TO_ELEMENT.items() if el != element
            ]
        else:
            raise ValueError(f"unknown perturbation mode: {mode}")
        if not candidates:
            continue
        new_pw = rng.choice(candidates)
        params["pathway"] = new_pw
        claim["params"] = params
        log.append({
            "claim_id": claim["id"],
            "claim_type": claim["type"],
            "element": element,
            "original_pathway": pw,
            "perturbed_pathway": new_pw,
            "mode": mode,
        })
    return new, log


def score_one(yaml_dict: dict, cycle_data, compare_df=None) -> dict:
    """评分一份 hypothesis dict，返回核心字段（不跑 null/sensitivity 节省时间）。"""
    from envmeta.geocycle.hypothesis import load_hypothesis, score
    hyp = load_hypothesis(yaml_dict)
    s = score(hyp, cycle_data, compare_df=compare_df, run_null=False, run_sensitivity=False)
    return {
        "overall_score": s.overall_score,
        "label": s.label,
        "n_satisfied": s.n_satisfied,
        "n_total": s.n_total,
        "n_skipped": s.n_skipped,
        "veto_reasons": ";".join(s.veto_reasons) if s.veto_reasons else "",
    }


def run_dataset(dataset_cfg: dict, n_perturb: int, seed_base: int = 0) -> list[dict]:
    """For one dataset: load data, run cycle once, score original + N perturbations
    (within-element + cross-element)."""
    key = dataset_cfg["key"]
    label = dataset_cfg["label"]
    raw_yaml = dataset_cfg["yaml"]
    yaml_path = Path(raw_yaml) if Path(raw_yaml).is_absolute() else EXTERNAL / key / raw_yaml
    restrict_claim_types = dataset_cfg.get("restrict_claim_types")
    compare_groups_arg = dataset_cfg.get("compare_groups")

    print(f"\n{'=' * 70}\n  PERTURBATION RUN -- {label}\n{'=' * 70}")
    if restrict_claim_types:
        print(f"  [partial] restricted to claim types: {restrict_claim_types}")

    if not yaml_path.exists():
        print(f"  [SKIP] missing yaml: {yaml_path}")
        return []

    keystone = None
    if dataset_cfg.get("data_loader") == "sample_data":
        metadata, env, abund, ko_long, qual, tax, keystone = load_inputs_sample_data()
        print(f"  [+keystone] {len(keystone)} entries loaded for keystone_in_pathway claim")
    else:
        data_dir = EXTERNAL / key / "input_data_local"
        if not data_dir.exists():
            print(f"  [SKIP] missing data dir: {data_dir}")
            return []
        metadata, env, abund, ko_long, qual, tax = load_inputs(data_dir)
    print(f"  inputs: {len(metadata)} samples / {len(qual)} MAGs / {abund.shape[0]} abund rows")

    print("  [1/N+1] Cycle diagram (1x cached)...")
    from envmeta.analysis import cycle_diagram
    t0 = time.time()
    cycle_result = cycle_diagram.analyze(
        ko_annotation_df=ko_long, taxonomy_df=tax,
        abundance_df=abund, env_df=env, metadata_df=metadata,
        keystone_df=keystone,
    )
    print(f"      cycle runtime: {time.time() - t0:.2f} s")

    compare_df = None
    if compare_groups_arg:
        print(f"  [+compare_groups] for group_contrast claims: {compare_groups_arg}")
        from envmeta.analysis import cycle_compare
        compare_df = cycle_compare.compare_groups(
            ko_annotation_df=ko_long, taxonomy_df=tax,
            abundance_df=abund, env_df=env, metadata_df=metadata,
            keystone_df=keystone,
            groups=compare_groups_arg,
        )

    with open(yaml_path, "r", encoding="utf-8") as f:
        orig_dict = yaml.safe_load(f)

    rows: list[dict] = []

    # Original
    print("  [original] scoring...")
    s_orig = score_one(orig_dict, cycle_result.data, compare_df=compare_df)
    rows.append({
        "dataset": label,
        "mode": "original",
        "run_idx": -1,
        "is_original": True,
        "seed": None,
        "n_perturbed_claims": 0,
        "perturbations_log": "",
        **s_orig,
    })
    print(f"      overall={s_orig['overall_score']:.3f} label={s_orig['label']}"
          f" satisfied={s_orig['n_satisfied']}/{s_orig['n_total']}")

    # Within-element perturbations
    for i in range(n_perturb):
        seed = seed_base + i
        rng = random.Random(seed)
        perturbed_dict, log = perturb_yaml(
            orig_dict, rng, mode="within_element",
            restrict_claim_types=restrict_claim_types,
        )
        s = score_one(perturbed_dict, cycle_result.data, compare_df=compare_df)
        rows.append({
            "dataset": label,
            "mode": "within_element",
            "run_idx": i,
            "is_original": False,
            "seed": seed,
            "n_perturbed_claims": len(log),
            "perturbations_log": json.dumps(log, ensure_ascii=False),
            **s,
        })
        print(f"      [within] {i:02d} seed={seed} overall={s['overall_score']:.3f}"
              f" label={s['label']} satisfied={s['n_satisfied']}/{s['n_total']}"
              f" (n_pert={len(log)})")

    # Cross-element perturbations (stronger negative control)
    for i in range(n_perturb):
        seed = seed_base + 1000 + i
        rng = random.Random(seed)
        perturbed_dict, log = perturb_yaml(
            orig_dict, rng, mode="cross_element",
            restrict_claim_types=restrict_claim_types,
        )
        s = score_one(perturbed_dict, cycle_result.data, compare_df=compare_df)
        rows.append({
            "dataset": label,
            "mode": "cross_element",
            "run_idx": i,
            "is_original": False,
            "seed": seed,
            "n_perturbed_claims": len(log),
            "perturbations_log": json.dumps(log, ensure_ascii=False),
            **s,
        })
        print(f"      [cross]  {i:02d} seed={seed} overall={s['overall_score']:.3f}"
              f" label={s['label']} satisfied={s['n_satisfied']}/{s['n_total']}"
              f" (n_pert={len(log)})")

    return rows


def make_summary(all_rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(all_rows)
    out = []
    for (ds, mode), group in df[df["mode"] != "original"].groupby(["dataset", "mode"]):
        orig = df[(df["dataset"] == ds) & (df["mode"] == "original")].iloc[0]
        perturbed = group
        out.append({
            "dataset": ds,
            "mode": mode,
            "original_score": orig["overall_score"],
            "original_label": orig["label"],
            "perturbed_n": len(perturbed),
            "perturbed_score_min": perturbed["overall_score"].min(),
            "perturbed_score_median": perturbed["overall_score"].median(),
            "perturbed_score_max": perturbed["overall_score"].max(),
            "perturbed_score_mean": perturbed["overall_score"].mean(),
            "perturbed_n_strong": int((perturbed["label"] == "strong").sum()),
            "perturbed_n_suggestive": int((perturbed["label"] == "suggestive").sum()),
            "perturbed_n_weak": int((perturbed["label"] == "weak").sum()),
            "perturbed_n_insufficient": int((perturbed["label"] == "insufficient").sum()),
            "fraction_label_strong": float((perturbed["label"] == "strong").mean()),
            "fraction_label_changed": float((perturbed["label"] != orig["label"]).mean()),
        })
    return pd.DataFrame(out).sort_values(["dataset", "mode"]).reset_index(drop=True)


def make_figure(all_rows: list[dict], out_path: Path):
    """Two-panel figure: within-element vs cross-element perturbation per dataset.

    Each dataset has a column with original (red star) + 20 within-element points
    (green-orange-red by label) + 20 cross-element points (offset right).
    """
    import numpy as np
    df = pd.DataFrame(all_rows)
    datasets = list(df[df["mode"] == "original"]["dataset"])
    n = len(datasets)

    fig, ax = plt.subplots(1, 1, figsize=(8.5, 4.5), dpi=300)

    label_colors = {
        "strong": "#2E8B57",
        "suggestive": "#F39C12",
        "weak": "#E67E22",
        "insufficient": "#C0392B",
    }

    # X-positions: each dataset gets two sub-columns (within / cross)
    bar_w = 0.32
    x_within = np.arange(n) - bar_w / 2
    x_cross = np.arange(n) + bar_w / 2

    for i, ds in enumerate(datasets):
        sub = df[df["dataset"] == ds]
        orig = sub[sub["mode"] == "original"].iloc[0]
        within = sub[sub["mode"] == "within_element"]
        cross = sub[sub["mode"] == "cross_element"]

        rng = np.random.default_rng(i)

        # within-element
        jit = rng.uniform(-0.10, 0.10, size=len(within))
        for j, (_, row) in enumerate(within.iterrows()):
            color = label_colors.get(row["label"], "#666")
            ax.scatter(x_within[i] + jit[j], row["overall_score"], s=22, alpha=0.7,
                       color=color, edgecolors="none", zorder=2)

        # cross-element
        jit = rng.uniform(-0.10, 0.10, size=len(cross))
        for j, (_, row) in enumerate(cross.iterrows()):
            color = label_colors.get(row["label"], "#666")
            ax.scatter(x_cross[i] + jit[j], row["overall_score"], s=22, alpha=0.7,
                       color=color, edgecolors="none", marker="s", zorder=2)

        # Median lines
        med_w = within["overall_score"].median()
        med_c = cross["overall_score"].median()
        ax.hlines(med_w, x_within[i] - 0.13, x_within[i] + 0.13,
                  colors="#222", linewidth=1.6, zorder=3, linestyles="-")
        ax.hlines(med_c, x_cross[i] - 0.13, x_cross[i] + 0.13,
                  colors="#222", linewidth=1.6, zorder=3, linestyles="-")

        # Original red star (between the two columns)
        ax.scatter(i, orig["overall_score"], s=200, marker="*",
                   color="#C0392B", edgecolors="black", linewidths=0.9, zorder=4,
                   label="Original (calibration)" if i == 0 else None)

    # Threshold lines
    ax.axhline(0.75, color="#2E8B57", linestyle=":", linewidth=0.8, alpha=0.6, zorder=1)
    ax.axhline(0.40, color="#F39C12", linestyle=":", linewidth=0.8, alpha=0.6, zorder=1)
    ax.text(n - 0.55, 0.76, "STRONG (0.75)", fontsize=8, color="#2E8B57", va="bottom")
    ax.text(n - 0.55, 0.41, "SUGGESTIVE (0.40)", fontsize=8, color="#F39C12", va="bottom")

    # Sub-column labels
    for i in range(n):
        ax.text(x_within[i], -0.05, "within", ha="center", fontsize=8,
                color="#444", transform=ax.get_xaxis_transform())
        ax.text(x_cross[i], -0.05, "cross", ha="center", fontsize=8,
                color="#444", transform=ax.get_xaxis_transform())

    ax.set_xticks(range(n))
    ax.set_xticklabels(datasets, fontsize=10)
    ax.tick_params(axis="x", which="major", pad=18)
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("overall score (after perturbation)", fontsize=11)
    ax.set_title("Pathway-target perturbation: within-element vs cross-element\n"
                 "(N=20 each / dataset; circle=within, square=cross, red⋆=original)",
                 fontsize=10, pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower left", fontsize=8, frameon=False)

    # Color legend (label colors)
    from matplotlib.lines import Line2D
    legend_lines = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=label_colors["strong"],
               markersize=7, label="strong"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=label_colors["suggestive"],
               markersize=7, label="suggestive"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=label_colors["insufficient"],
               markersize=7, label="insufficient (veto)"),
    ]
    ax.legend(handles=legend_lines, loc="upper left", fontsize=8, frameon=False,
              title="Perturbed label", title_fontsize=8)

    plt.tight_layout()
    for ext in ("pdf", "png", "svg"):
        fig.savefig(out_path.with_suffix(f".{ext}"))
    plt.close(fig)
    print(f"  figure: {out_path.with_suffix('.pdf').name} (+png/svg)")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=20, help="perturbations per dataset")
    p.add_argument("--seed", type=int, default=0, help="seed base")
    args = p.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []
    for cfg in DATASETS:
        rows = run_dataset(cfg, n_perturb=args.n, seed_base=args.seed)
        all_rows.extend(rows)

    if not all_rows:
        print("No results — abort.")
        return

    df_full = pd.DataFrame(all_rows)
    df_full.to_csv(OUT_DIR / "perturbation_results.tsv", sep="\t", index=False)

    df_summary = make_summary(all_rows)
    df_summary.to_csv(OUT_DIR / "perturbation_summary.tsv", sep="\t", index=False)

    print("\n" + "=" * 70 + "\n  SUMMARY\n" + "=" * 70)
    print(df_summary.to_string(index=False))

    make_figure(all_rows, OUT_DIR / "perturbation_curve")

    print(f"\nResults written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
