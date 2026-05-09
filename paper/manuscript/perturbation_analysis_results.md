# Perturbation Analysis — Auxiliary Evidence for Calibration Specificity

> **Purpose**: Address Mock Review v0.9.2 Major Issue #1 — provide auxiliary
> evidence that the four STRONG calibration outcomes are not mechanical
> artifacts of KEGG annotation density / pre-data target choice. Run as the
> "1-day perturbation analysis" alternative to the deferred blind-hypothesis-
> writing exercise (Discussion §Y.4 future-work item).
>
> **Date**: 2026-05-09
> **Engine commit**: (current `master`)
> **Runner**: [`tools/external_benchmarks/perturbation_analysis.py`](../../tools/external_benchmarks/perturbation_analysis.py)
> **Outputs**: [`paper/benchmarks/external/perturbation/`](../benchmarks/external/perturbation/)

---

## 1. Design

### Hypothesis under test

If the four pre-registered calibration STRONG outcomes (作者 / Liu 2023 /
Grettenberger 2021 / Ayala 2020) arise *mechanically* from KEGG annotation
density rather than from authors' specific pre-data target choices, then
**replacing the target pathways with random alternatives should still yield
STRONG outcomes at high rates**. Conversely, if the original choices are
substantively informed, **target perturbation should degrade scores
substantially**.

### Two perturbation modes

For each of three external calibration YAMLs (Liu 2023, Grettenberger 2021,
Ayala 2020 — the author-data Arm A is excluded due to its joint
`coupling_possible` + `group_contrast` claims that are not pathway-targeted),
we perturb every claim whose `params.pathway` field is set:

| Mode | Replacement rule | Expected outcome under "calibration is real" |
|---|---|---|
| **within_element** | random different pathway from same KB element (e.g., `Arsenate reduction` → `As methylation`) | Moderate score degradation; score retention reflects KEGG-coverage breadth |
| **cross_element** | random pathway from **different** KB element (e.g., `Arsenate reduction` → `Sulfide oxidation`) | Strong score degradation; required-claim veto kicks in for absent pathways |

All other YAML fields (weight, required, min_completeness, env_factor for
`env_correlation`, expected_sign) are kept unchanged. **N=20** perturbations
per mode per dataset → 120 perturbed scorings + 3 originals = **123 runs**.
Random seeds are deterministic (within: 0–19; cross: 1000–1019) for
reproducibility.

### Engine settings

`run_null=False, run_sensitivity=False` (saves time; we only need
`overall_score` + `label`). Default thresholds (strong=0.75, suggestive=0.40)
and default `min_completeness=30` from each YAML.

---

## 2. Results

### Headline table

[`paper/benchmarks/external/perturbation/perturbation_summary.tsv`](../benchmarks/external/perturbation/perturbation_summary.tsv)

| Dataset | Mode | Median | Mean | Strong fraction | Label changed |
|---|---|---|---|---|---|
| Liu 2023 (As-focused, cold seep) | within_element | 0.770 | 0.787 | **10/20 (50%)** | 50% |
| Liu 2023 | **cross_element** | **0.000** | **0.000** | **0/20 (0%)** | **100%** |
| Grettenberger 2021 (AMD stream) | within_element | 0.772 | 0.753 | 10/20 (50%) | 50% |
| Grettenberger 2021 | cross_element | 0.455 | 0.471 | 6/20 (30%) | 70% |
| Ayala 2020 (AMD pit lake) | within_element | 0.500 | 0.521 | 8/20 (40%) | 60% |
| Ayala 2020 | cross_element | 0.500 | 0.464 | 3/20 (15%) | 85% |

(Original calibration scores: all three datasets gave 1.000 STRONG.)

### Three findings

1. **Cross-element perturbation is strongly discriminating.** When target
   pathways are replaced with pathways from a *different* KB element, label
   degradation is 70–100% across the three datasets. Liu 2023 is the most
   striking case: **0/20 cross-element perturbations retained STRONG**, with
   median score collapsing to 0.000. This is because Liu 2023 cold seep data
   is As-element-focused; replacing arsenic-pathway claims with N/S/Fe
   pathway claims triggers the required-claim veto (no active MAGs).

2. **Within-element perturbation degrades scores moderately but allows
   incidental STRONG.** Within-element mean scores are 25–48% below the
   original 1.000, but 40–50% of within-element perturbations *still*
   produce STRONG. This is consistent with the manuscript's
   "KEGG-coverage-dependent" framing (Discussion §Y.1): when a dataset has
   high KEGG annotation breadth across an entire element, multiple
   parallel pathway-active claims will register satisfied irrespective of
   which subset the author selected. **This bounds the calibration claim,
   it does not invalidate it**: a dataset's STRONG outcome reflects both
   author target choice (real signal) AND annotation breadth in the
   targeted element (KEGG-coverage effect).

3. **Smallest dataset is most discriminating.** Ayala 2020 (n=13 MAGs) has
   the lowest fraction of within-element STRONG (40%) and the lowest mean
   score (0.521). Liu and Grettenberger (n=29, n>1000) keep ~50% STRONG.
   This bounds within-element perturbation as a function of annotation
   breadth: small or sparsely-annotated datasets are more sensitive to
   target choice, large/densely-annotated datasets are less so.

### Why "0/20 STRONG for Liu cross-element" is the strongest single result

Liu 2023 is the cleanest test case because (a) the dataset is monolithically
As-focused (no functional N/S/Fe activity for cross-element replacements to
hit), and (b) the YAML has 2 of 4 claims marked `required=true`. Cross-element
replacement therefore systematically lands on inactive pathways, triggering
veto on the required claims, collapsing the label to insufficient. The
contrast between within-element 50% STRONG and cross-element 0% STRONG —
**under identical YAML weight structure and engine settings** — is
mechanistically diagnostic: the calibration result depends on element-level
target accuracy, not just on summing satisfied claims.

### Figure

`paper/benchmarks/external/perturbation/perturbation_curve.{pdf,png,svg}` —
3-column scatter (one per dataset) showing original (red star) + 20
within-element points (circles) + 20 cross-element points (squares),
colored by label, with within/cross median lines and STRONG/SUGGESTIVE
threshold reference lines.

---

## 3. Limitations

1. **Within-element saturation effect**: when a dataset has rich KEGG
   coverage across an entire element, within-element perturbation
   under-discriminates because most candidate pathways are also active.
   This is a *feature* (it correctly reflects KEGG-coverage dependence),
   not a bug. But it means "fraction of within-element perturbations
   retaining STRONG" should not be interpreted as a false-positive rate.

2. **Pre-registered claim count differs across datasets**. Liu has 4 claims,
   Grettenberger has 4 claims, Ayala has 4 claims — each with different
   `required`/non-required composition. We did not normalize for this; the
   raw fraction-STRONG numbers should be read as dataset-specific, not
   directly comparable.

3. **Author-data Arm A excluded**. The author calibration includes
   `coupling_possible` + `group_contrast` + `env_correlation` claims whose
   pathway-targeting is intertwined with species coupling and group-level
   comparison. A clean "perturb only the pathway field" analysis is not
   meaningful for those claims; including the author dataset would require
   designing a separate species/group perturbation, which is out of scope
   for this 1-day exercise.

4. **No coupling/keystone perturbation**. We perturb only the
   `params.pathway` field. The four calibration YAMLs do not heavily use
   `coupling_possible` or `keystone_in_pathway` claims; for any future
   coupling-rich YAML, a separate species perturbation analysis would be
   required.

5. **N=20 is small**. Confidence intervals on the "strong fraction" are
   broad (Wilson 95% for 10/20 ≈ [0.30, 0.70]). For publication-ready
   precision, N=200 or higher would tighten estimates. We treat the
   present result as auxiliary evidence rather than a primary statistical
   test.

---

## 4. Interpretation for the manuscript

The perturbation analysis provides **auxiliary evidence consistent with**
(but not ironclad proof of) the claim that the four pre-registered
calibration STRONG outcomes are not mechanical artifacts of KEGG annotation
density alone:

- **Cross-element perturbation result**: 70–100% label degradation (and 0%
  STRONG retention for Liu) shows that **element-level target accuracy
  matters**. A randomly chosen target from outside the data's biogeochemical
  domain does not reproduce the calibration.

- **Within-element perturbation result**: moderate score degradation
  (25–48% mean drop) with 40–50% incidental STRONG retention quantifies
  the "KEGG-coverage-dependent" caveat already acknowledged in Discussion
  §Y.1. We do not over-claim that within-element perturbation eliminates
  STRONG; we report it as a soft upper bound on annotation-density-driven
  inflation.

- **Limitation of this analysis**: the perturbation tests author-target-
  *choice* sensitivity, not author-*selection*-bias for the outcome. The
  blind-hypothesis-writing exercise (Discussion §Y.4) remains the gold
  standard for the latter and is preserved as future work. The
  perturbation analysis is offered as the auxiliary evidence that Mock
  Review v0.9.2 (Major #1) suggested, not as a substitute.

### Suggested location in manuscript

| Location | Treatment |
|---|---|
| **Methods §4.6.6** (NEW) | Brief subsection: "Auxiliary perturbation analysis" — protocol + reference to this doc |
| **Results §X.3** (NEW) | One paragraph: cross-element 0/20 STRONG for Liu as headline; within/cross summary table |
| **Discussion §Y.3** | Update Major #1 limitation: "perturbation analysis (Methods §4.6.6) provides auxiliary evidence; blind-writing remains future work for definitive bias assessment" |
| **Figure** | `perturbation_curve.pdf` as Figure X-bis or supplementary |
| **Table** | `perturbation_summary.tsv` as Table X-bis or supplementary table |

---

## 5. Reproduction

```powershell
conda activate envmeta
cd D:\workdata\envmeta
python tools/external_benchmarks/perturbation_analysis.py --n 20

# Outputs:
#   paper/benchmarks/external/perturbation/perturbation_results.tsv
#   paper/benchmarks/external/perturbation/perturbation_summary.tsv
#   paper/benchmarks/external/perturbation/perturbation_curve.{pdf,png,svg}
```

Runtime ≈ 30 s (3 cycle_diagram runs × ~1.2 s + 123 hypothesis scorings).
Deterministic — same `--seed` returns identical results.

---

## Maintenance log

| Date | Event |
|---|---|
| 2026-05-09 | Initial run; both modes; N=20. Cross-element Liu 0/20 STRONG = headline result. |
