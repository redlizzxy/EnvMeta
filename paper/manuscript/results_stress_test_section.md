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
experiment over four metagenomic datasets spanning a gradient of annotation
breadth: our in-house steel-slag arsenic-remediation dataset (Arm A; 168 MAGs ×
10 samples; full KofamScan KEGG annotation, 57 KOs across 4 elements), Wei
et al. (Arm B; 2024 *Microbiome*, 36 paddy-soil samples × 179 MAGs annotated
with 14 functional genes via custom ROCker models¹), Liu et al. (Arm C1; 2023
*npj Biofilms Microbiomes*, deep-sea cold-seep, 87 samples × 1084 MAGs with
DRAM KEGG annotation), Grettenberger & Hamilton (Arm C2-A; 2021 *Appl Environ
Microbiol*, AMD stream, 29 MAGs with METABOLIC step-level KEGG annotation),
and Ayala-Muñoz et al. (Arm C2-B; 2020 *Microorganisms*, Iberian Pyrite Belt
acidic pit lake deep layer, 13 MAGs re-annotated end-to-end with Pyrodigal
plus GhostKOALA²). All hypothesis YAMLs were pre-registered (committed to git
at hashes `42168da`, `44d7f5f`, `76a4f77` before EnvMeta was run), used
EnvMeta's default thresholds (`min_completeness=30`, `strong=0.75`,
`suggestive=0.40`), and cited only review literature published five or more
years before the target paper.

All four KEGG-curated arms returned **`STRONG` labels with overall_score = 1.000
and 4/4 claims satisfied** (Table 1; Arms A, C1, C2-A, C2-B). The single
ROCker-only arm (Arm B, Wei 2024) returned overall_score = 0.63 with label
`INSUFFICIENT` despite 3/5 claims satisfied; the required-veto activated
because Wei's published 14-gene set provides only 2 of the 6 canonical KOs
in EnvMeta's `Nitrate reduction` pathway (`napA + narG`, missing `narH/narI/
napB/narB`) and the As(III)↔NO₃⁻ chemistry coupling was scored *partial*
because one of its termini was unsatisfied. Permutation null-p = 0.90
(n = 999) and weight robustness under ±20% one-at-a-time (OAT) perturbation
= True confirmed the conservative diagnosis is not a weight-tuning artifact.
The contrast among the four arms — same scoring engine, identical default
thresholds — establishes that the `INSUFFICIENT` label on Arm B faithfully
reflects annotation-coverage diagnostics rather than engine malfunction or
threshold mismatch. Notably, the `STRONG` label on Arm C2-B (Ayala) was
obtained with a single GhostKOALA re-annotation pipeline applied to publicly
available MAG genomes (BioProject PRJNA646106), demonstrating that
KEGG-curated EnvMeta scoring is reproducible end-to-end from raw genome
assemblies. We emphasize that this is **calibration evidence**: all
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
scores). The single most informative discriminator was the cross-topic claim
"arsenate_reduction should dominate", which was correctly rejected with
**n = 0 active MAGs in both non-arsenic datasets** (Grettenberger and Ayala).
This double-rejection rules out the *a priori* concern that the universal
*arsC* arsenate-reductase homolog (Rosen, 2002) would inflate cross-topic
scores to satisfied: in environments without sustained arsenic pressure,
KEGG annotations of *arsC* and related operon members are absent or below
detection threshold, and EnvMeta's scoring engine reflects this faithfully.
Together with the `INSUFFICIENT` label on Arm B, the converging stress
results provide direct evidence that EnvMeta's hypothesis scoring is **domain-
neutral** — it is neither hard-wired to confirm arsenic-cycle hypotheses
(which would falsely "pass" the AMD and pit-lake stress tests) nor biased
against alternative environmental contexts.

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

## §X.3 Reference audit and corrections

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

## §X.4 Tables and figures

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
