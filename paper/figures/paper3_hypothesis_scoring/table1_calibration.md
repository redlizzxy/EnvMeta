# Table 1. Four-Arm Calibration Results

> Paper 3 Results §X.1 — calibration of EnvMeta hypothesis scoring engine across
> four metagenomic datasets spanning a gradient of annotation breadth and study
> topics. All hypothesis YAMLs pre-registered (committed to git before EnvMeta
> runs); EnvMeta default thresholds (`min_completeness=30`, `strong=0.75`,
> `suggestive=0.40`) used unchanged.

| Arm | Dataset | Sample × MAG | Annotation | overall_score | Label | Claims (sat/total) | Topic |
|---|---|---|---|---|---|---|---|
| **A** | In-house steel-slag (CK/A/B) | 10 samples × 168 MAGs | KofamScan KEGG (full, 57 KOs / 4 elements) | **1.000** | **STRONG** | n/n | Arsenic remediation (high-As) |
| **B** | Wei 2024 *Microbiome* | 36 samples × 179 MAGs | Custom ROCker (14 genes only) | 0.63 | **INSUFFICIENT** | 3/5 (required veto) | Arsenic-contaminated paddy soils |
| **C1** | Liu 2023 *npj Biofilms Microbiomes* | 87 samples × 1084 MAGs | DRAM (KEGG-curated) | **1.000** | **STRONG** | 4/4 | Deep-sea cold seep (high-As) |
| **C2-A** | Grettenberger 2021 *Appl Environ Microbiol* | 1 site × 29 MAGs | METABOLIC step-level KEGG | **1.000** | **STRONG** | 4/4 | AMD stream — **non-arsenic** ⭐ |
| **C2-B** | Ayala 2020 *Microorganisms* | 1 sample × 13 MAGs | Pyrodigal + GhostKOALA (KEGG) | **1.000** | **STRONG** | 4/4 | Pit-lake deep layer — **non-arsenic** ⭐ |

⭐ = cross-topic dataset (no dominant arsenic in environment); see Table 2 for stress-test cross-topic discrimination.

## Notes

- **Arm B INSUFFICIENT**: required-veto activated because Wei's ROCker-only 14-gene set provides only 2/6 canonical KOs of EnvMeta's `Nitrate reduction` pathway (`napA + narG`), giving per-MAG completeness ≤ 33% and the As(III)↔NO₃⁻ coupling scored *partial* due to one terminus unsatisfied. null_p = 0.90 (n=999), weight_robust = True. This is calibration evidence (faithful annotation-coverage diagnostic), not engine malfunction.
- **Arm C2-B GhostKOALA pipeline**: 13 MAG genomes from BioProject PRJNA646106 → Pyrodigal ORF prediction (24,841 proteins) → GhostKOALA online annotation (45.2% hit rate, 11,234 proteins) → reshape to 6 EnvMeta inputs. End-to-end reproducible from raw genome FASTA.
- All four KEGG-curated calibration `STRONG` results use **claims targeting backbone biogeochemical pathways** near-universally expected in their respective environments (e.g. arsenate reduction in anoxic high-As; sulfide oxidation in AMD). This is *calibration*, not *discrimination* evidence.

## Pre-registration anchors (git timestamps)

| Arm | YAML pre-registered at | Run results at |
|---|---|---|
| A | (in-house dataset; see `paper/hypotheses/arsenic_steel_slag.yaml`) | — |
| B | `paper/benchmarks/external/wei_2024_paddy/wei2024_hypothesis.yaml` | — |
| C1 | commit `42168da` | commit `1e4571f` |
| C2-A | commit `44d7f5f` | commit `e8d5133` |
| C2-B | commit `76a4f77` | commit `60a5be4` |

## Source data

`paper/manuscript/scoring_validation_experiment_results.md` §1 + per-dataset
`fig6_*_score.md` files in each dataset's `envmeta_outputs/` folder.
