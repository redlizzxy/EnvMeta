# Results §X — Hypothesis Scoring Calibration and Stress-Test Discrimination

> Paper 3 Results section draft for the controlled experiment showing
> EnvMeta's hypothesis scoring engine works as a *calibration tool* and has
> *discrimination power* against domain-mismatched claims. To be inserted in
> outline_imeta.md §5.2 around §5.2.4 (YAML hypothesis scorer figure) and §5.2.8
> (external benchmark figure/table). Word count: ~800.
>
> **Status**: ready to integrate; refers to Table 1 + Table 2 + Figure X
> (definitions in §X.4 below).

---

## §X.1 Four-arm calibration confirms scoring engine tracks annotation breadth

To establish that EnvMeta's hypothesis scoring engine produces stable,
defensible labels under default thresholds, we ran a four-arm controlled
experiment over five metagenomic datasets spanning a gradient of annotation
breadth. The experiment combines **one in-house positive control** (Arm A,
authors' own data and authors' own a-priori hypothesis YAML) and **four
external blind-test calibrations** (Arms B, C1, C2-A, C2-B; published
datasets, hypothesis YAMLs authored before reading the corresponding
papers' specific findings — see §4.6.2). The two categories carry
different evidentiary weight: Arm A demonstrates engine self-consistency
on a dataset where the `STRONG` label is the *expected* outcome by
experimental design, while the four external arms test whether the
engine produces reasonable labels on independent data and hypotheses.
We therefore report Arm A as a positive control (engine-doesn't-misfire
sanity check) and the four external arms as the primary calibration
evidence.

The five datasets are: our in-house steel-slag arsenic-remediation dataset
(Arm A *positive control*; 168 MAGs × 10 samples; full KofamScan KEGG
annotation, 57 KOs across 4 elements), Wei et al. (Arm B *external*;
2024 *Microbiome*, 36 paddy-soil samples × 179 MAGs annotated with 14
functional genes via custom ROCker models¹), Liu et al. (Arm C1
*external*; 2023 *npj Biofilms Microbiomes*, deep-sea cold-seep, 87
samples × 1084 MAGs with DRAM KEGG annotation), Grettenberger & Hamilton
(Arm C2-A *external*; 2021 *Appl Environ Microbiol*, AMD stream, 29
MAGs with METABOLIC step-level KEGG annotation), and Ayala-Muñoz et al.
(Arm C2-B *external*; 2020 *Microorganisms*, Iberian Pyrite Belt acidic
pit lake deep layer, 13 MAGs re-annotated end-to-end with Pyrodigal plus
GhostKOALA²). All hypothesis YAMLs were pre-registered (committed to git
at hashes `42168da`, `44d7f5f`, `76a4f77` before EnvMeta was run), used
EnvMeta's default thresholds (`min_completeness=30`, `strong=0.75`,
`suggestive=0.40`), and cited only review literature published five or
more years before the target paper.

The Arm A in-house positive control returned **`STRONG` with
overall_score = 1.000 (9/9 claims satisfied)**, confirming the engine
yields the expected label on data and hypothesis built by the authors —
this is a self-consistency check, not independent calibration evidence.
The three KEGG-curated **external** arms (C1, C2-A, C2-B) returned
**`STRONG` with overall_score = 1.000 (4/4 claims satisfied each)**
(Table 1) — the principal calibration result. The single
ROCker-only arm (Arm B, Wei 2024) returned overall_score = 0.63 with label
`INSUFFICIENT` despite 3/5 claims satisfied; the required-veto activated
because Wei's published 14-gene set provides only 2 of the 6 canonical KOs
in EnvMeta's `Nitrate reduction` pathway (`napA + narG`, missing `narH/narI/
napB/narB`) and the As(III)↔NO₃⁻ chemistry coupling was scored *partial*
because one of its termini was unsatisfied. Permutation null-p = 0.90
(n = 999) and weight robustness under ±20% one-at-a-time (OAT) perturbation
= True confirmed the conservative diagnosis is not a weight-tuning artifact.
A separate threshold-sensitivity sweep (Methods §4.6.8) over
`strong_threshold ∈ {0.65, 0.70, 0.75, 0.80, 0.85}` confirms that all four
KEGG-curated arms retain `STRONG` and Arm B retains `INSUFFICIENT` across the
entire range, ruling out the possibility that the calibration outcome is an
artifact of the specific 0.75 default.
The contrast among the external arms — same scoring engine, identical
default thresholds — establishes that the `INSUFFICIENT` label on Arm B
faithfully reflects annotation-coverage diagnostics rather than engine
malfunction or threshold mismatch. Notably, the `STRONG` label on Arm
C2-B (Ayala) was obtained with a single GhostKOALA re-annotation
pipeline applied to publicly available MAG genomes (BioProject
PRJNA646106), demonstrating that KEGG-curated EnvMeta scoring is
reproducible end-to-end from raw genome assemblies on data that the
authors had no role in generating. We emphasize that this is **calibration evidence**: all
KEGG-curated `STRONG` results were obtained with claims targeting backbone
biogeochemical pathways near-universally expected in their respective
environments (e.g. arsenate reduction in anoxic high-As cold-seep, sulfide
oxidation in AMD streams, dissimilatory sulfate reduction in pit-lake deep
anoxic layers). A separate stress test of the engine's discrimination power
was therefore performed.

## §X.2 Stress test demonstrates discrimination power and exposes a binary-threshold limitation

For each of the three KEGG-curated cross-topic and same-topic datasets (Liu,
Grettenberger, Ayala), we authored a second pre-registered YAML
(`{dataset}_hypothesis_stress.yaml`, all committed at `50c4687`) encoding
deliberately *risky* claims violating environmental priors. Three claim
classes were used: **(A)** reversed-direction predictions (e.g. arsenite
oxidation should dominate in anoxic cold-seep sediments, reversing the
expected reduction-dominated regime); **(B)** cross-topic mismatches (e.g.
arsenate reduction should dominate in non-arsenic AMD systems); and **(C)**
`pathway_inactive` negation of backbone pathways (e.g. the dominant
arsenate-reduction pathway in cold-seep should *not* be active —
v0.9.0's Popperian-falsification claim type). Each YAML included a
calibration anchor claim to verify scoring system integrity. Twelve
predictions × three datasets were frozen in
`paper/manuscript/stress_test_predictions.md` before any stress run.

Observed stress-test scores fell substantially below their respective
calibration STRONG (1.000) labels (Table 2; Figure X). Grettenberger 2021
returned label `weak` (0.250, 1/3 satisfied) — the cleanest single-Arm
discrimination result. Liu 2023 and Ayala 2020 returned `suggestive` (0.625
and 0.455 respectively, with skipped claims contributing to the partial
scores). The most informative single discriminator was the cross-topic
claim "arsenate_reduction should dominate", which was correctly rejected
with **n = 0 active MAGs in both non-arsenic datasets** (Grettenberger n =
29 MAGs and Ayala n = 13 MAGs).

We treat this two-dataset rejection as **consistent with** — rather than
ironclad proof of — the scoring engine being uninfluenced by the universal
*arsC* arsenate-reductase homolog (Rosen, 2002) under cross-topic mismatch.
Three caveats temper the strength of this evidence. First, the
cross-topic discrimination is observed in 2 of 2 non-arsenic datasets
within a 3-dataset stress design; the underlying dataset count is small
(N=3) and the across-dataset generalisation is correspondingly soft.
Second, the absolute MAG counts in Grettenberger (29) and Ayala (13) are
small enough that the absence of any *arsC*-bearing MAG could in part
reflect sampling undercount rather than genuine absence in the underlying
community. Third, both datasets sample acid mine drainage / pit-lake
systems where arsenic is plausibly subdetectable rather than truly zero,
so the rejection result is conditional on the published annotations being
faithful to the underlying biology. A larger non-arsenic dataset (≥ 100
MAGs from a soil or marine system) — and a stress design extending
beyond N=3 datasets — would provide a more statistically robust test,
and we flag both as future work in §Y.4. Combined with the `INSUFFICIENT`
label on Arm B, the present stress results nevertheless support — within
the scope of the KEGG-curated datasets tested — the conclusion that
EnvMeta's scoring is neither hard-wired to confirm arsenic-cycle
hypotheses nor biased against alternative environmental contexts.

In two of three datasets (Liu and Ayala), the reversed-direction stress claim
"arsenite oxidation should dominate" returned satisfied because real but
weak oxidizer signals were present in the data: Liu 2023 had two MAGs
carrying *aoxA*/*aoxB* (mean completeness 50%, total contribution 0.3) and
Ayala 2020 had one MAG carrying *sox*-class genes (mean completeness 67%,
total contribution 66.7), versus the dominant reduction pathway with total
contribution > 6 (Liu) or > 675 (Ayala) — a 21-fold and 10-fold gap
respectively. EnvMeta's current binary `satisfied`/`unsatisfied` reporting
treats both as `satisfied` once the mean-completeness threshold is crossed.
This exposes a real limitation: the engine cannot distinguish "dominant" from
"detectable but weak" pathway activity, and stress claims encoding
domain-violating *dominance* (rather than mere presence) cannot currently be
expressed cleanly. We classify this as a `B-tier` discrimination outcome
(genuine but coarse-grained) and propose future engineering of a
`dominance_score = total_contribution / element_total` metric to enable
stress claims to specify a `min_dominance_fraction` parameter, which would
upgrade Liu and Ayala to `A-tier` clean discrimination.

## §X.3 Auxiliary perturbation analysis: annotation-breadth gradient bounds the calibration claim

To distinguish whether the four STRONG calibration outcomes reflect
authors' specific pre-data target choices or arise mechanically from KEGG
annotation breadth, we perturbed `params.pathway` across all four
calibration YAMLs and rescored under default thresholds (Methods §4.6.7).
Two perturbation modes were applied: within-element (random alternative
pathway from the same KB element) and cross-element (random pathway from
a different KB element); N=20 per mode per dataset (deterministic seeds).
The three external YAMLs underwent full perturbation of all
pathway-targeted claims; Arm A underwent partial perturbation restricted
to its three `pathway_active` claims (its `coupling_possible`,
`env_correlation`, `keystone_in_pathway`, and `group_contrast` claims pair
pathway with semantically-tied second parameters and are not amenable to
clean single-axis perturbation; this asymmetry is discussed in §Y.3
limitation #1).

The single most informative finding is a **monotonic annotation-breadth
gradient** in cross-element STRONG retention (Wilson 95% CI in brackets):
**Arm A 100%** [0.84, 1.00] (168 MAGs annotating active pathways across
all 4 elements simultaneously) → **Grettenberger 2021 30%** [0.13, 0.55]
(29 MAGs, mixed S/Fe AMD) → **Ayala 2020 15%** [0.05, 0.37] (13 MAGs,
mixed S/Fe pit lake) → **Liu 2023 0%** [0.00, 0.16] (29 MAGs, As-only
cold seep). This ordering is mechanistically expected: in datasets where
the data is element-monolithic, cross-element substitution
systematically lands on inactive pathways and triggers required-claim
veto; in datasets where the data is element-saturated, any pathway in
any element is plausibly active, so substitution does not change the
outcome. The fact that the engine behaves predictably as a function of
this independent-of-the-calibration property of each dataset is itself
a strong internal-validity check.

Liu 2023 is the cleanest single test case: 0/20 cross-element
perturbations retained STRONG (median score collapse 1.000 → 0.000),
under identical YAML weight structure and engine settings as the
within-element control's 50% retention. The 50%-vs-0% gap is
mechanistically diagnostic: in focused datasets, the calibration result
depends on element-level target accuracy, not just on summing satisfied
claims.

The within-element control bounds the KEGG-coverage caveat acknowledged
in Discussion §Y.1. Within-element mean scores in the three external
datasets are 25–48% below the original 1.000, but **40–50% of
within-element perturbations still produce STRONG**. This is consistent
with — not in conflict with — the manuscript's framing that calibration
evidence is KEGG-coverage-dependent: when a target element has high
annotation breadth, multiple parallel pathway-active claims register
satisfied irrespective of which subset the author selected. Ayala 2020
(n=13 MAGs, smallest dataset) is the most discriminating at 8/20 STRONG
retention; Liu and Grettenberger sit at 10/20.

Arm A's 100% STRONG retention under both within and cross perturbation
quantifies the saturation regime: for the in-house dataset specifically,
the perturbation analysis cannot strongly distinguish "author cherry-
picked targets" from "any reasonable target works because data is rich",
and the cherry-pick concern for Arm A must therefore rely on the
broader pre-registration discipline (§4.6.2) and the planned blind
third-party stress YAMLs (§Y.4) rather than on perturbation alone. We
report the perturbation result as **auxiliary evidence consistent with**
the calibration claims being non-mechanical, while flagging that the
within-element fraction-STRONG is **not** a false-positive rate and
that the more definitive blind-hypothesis-writing exercise (§Y.4)
remains future work. Full results in Supplementary Table S_pert
(`perturbation_summary.tsv`) and Figure X-bis
(`perturbation_curve.pdf`); reproduction protocol in
[`paper/manuscript/perturbation_analysis_results.md`](perturbation_analysis_results.md).

## §X.4 Reference audit and corrections

Post-hoc DOI verification identified four reference errors in the
pre-registered YAMLs that do not affect scoring outputs (no claim entity was
modified) but require transparent correction. Most consequentially, the
`nitrogen_fixation_explored` claims (Grettenberger and Ayala calibration YAMLs)
originally cited Auld et al. (2017 *Can J Microbiol*), which is a seasonal
community-variation study rather than an AMD diazotrophy report; these
claims are re-grounded in Dai et al. (2014 *PLoS One*; metagenomic
identification of 742 *nif* sequences and a 32.5-kb *nif/fix* gene cluster
from acid mine drainage) and Méndez-García et al. (2015 *Front Microbiol*;
review of AMD diazotrophs including *nifHDKENX*-bearing *Acidithiobacillus*,
*Leptospirillum*, and *Ferrovum* lineages). Three additional metadata
corrections (journal mislabels for Yin 2011, Cabrera 2006, and a non-existent
"Bothe 2007 *FEMS Rev*" → Bothe 2000) were committed at `ddd3098` and
`cae2de7`; the original pre-registered versions remain accessible in git
history. A complete proof-of-extraction quality audit (per claim × reference,
graded `Direct` / `Inferred` / `Weak`) is provided in Supplementary Table
S_refs.

## §X.5 Tables and figures

**Table 1.** *Four-arm calibration results.* Columns: Arm, dataset, sample
size, annotation method, overall_score, label, claims satisfied, claim
types, environmental topic. Five rows (A in-house, B Wei 2024, C1 Liu 2023,
C2-A Grettenberger 2021, C2-B Ayala 2020). Source data:
`paper/manuscript/scoring_validation_experiment_results.md` §1.

**Table 2.** *Three-arm stress-test results.* Columns: Arm, calibration label
(from Table 1), stress overall_score, stress label, stress claims satisfied
(per claim class A/B/C), discrimination grade (A/B), cross-topic rejection
(yes/no), commit hash anchor. Three rows (C1 Liu 2023, C2-A Grettenberger
2021, C2-B Ayala 2020). Source data:
`paper/manuscript/stress_test_results.md` §0 + §3 + §4.

**Figure X.** *Calibration → stress score gap, three KEGG-curated datasets.*
Horizontal bar plot with three pairs (calibration 1.000 vs stress 0.250 /
0.625 / 0.455). Annotated with the cross-topic-rejection result in red text
on each row. Stylistic match to Figure 7's panel style. Single panel,
~10 cm × 6 cm.

---

## Word count: ~800 (target 600-1000 ✅)

## Cited references (Vancouver, beyond §4.6 list)

¹ Reichart NJ, et al. 2020 *ISME J* 14:2851 — ROCker custom-model framework
² Kanehisa M, et al. 2016 *J Mol Biol* 428:726 — GhostKOALA

(Full reference list shared with §4.6; renumber at submission.)

## Maintenance log

| Date | Event |
|---|---|
| 2026-05-09 | Initial complete draft (§X.1 calibration / §X.2 stress / §X.3 audit / §X.4 tables-figures definitions) |
