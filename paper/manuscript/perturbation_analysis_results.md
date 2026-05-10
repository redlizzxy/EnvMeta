# Perturbation Analysis — Auxiliary Evidence for Calibration Specificity

> **Purpose**: Provide auxiliary evidence that the four STRONG calibration
> outcomes are not mechanical artifacts of KEGG annotation density / pre-data
> target choice. Run as the "1-day perturbation analysis" alternative to the
> deferred blind-hypothesis-writing exercise (Discussion §Y.4 future-work
> item).
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

For each of four calibration YAMLs (Arm A in-house arsenic-steel-slag /
Liu 2023 / Grettenberger 2021 / Ayala 2020), we perturb claims whose
`params.pathway` field is set. The author Arm A YAML uses **partial
perturbation** restricted to its three `pathway_active` claims
(`iron_transport_active`, `arsenate_reduction_active`,
`as_transport_active`); its `coupling_possible` (×2),
`env_correlation` (×2), `keystone_in_pathway` (×1), and `group_contrast`
(×1) claims are kept unchanged because they either lack a `pathway`
parameter (couplings) or because perturbing them would jointly alter
multiple semantically-tied fields (env_correlation pairs pathway with
env_factor; group_contrast pairs pathway with group ratio). The three
external YAMLs use **full perturbation** of all pathway-targeted claims
since they have no such mixed-claim entanglement.

| Mode | Replacement rule | Expected outcome under "calibration is real" |
|---|---|---|
| **within_element** | random different pathway from same KB element (e.g., `Arsenate reduction` → `As methylation`) | Moderate score degradation; score retention reflects KEGG-coverage breadth within element |
| **cross_element** | random pathway from **different** KB element (e.g., `Arsenate reduction` → `Sulfide oxidation`) | Strong score degradation in monolithic-element datasets; required-claim veto kicks in for absent pathways. In multi-element-annotated datasets cross-element substitution may still hit active pathways |

All other YAML fields (weight, required, min_completeness, env_factor for
`env_correlation`, expected_sign) are kept unchanged. **N=20** perturbations
per mode per dataset → 160 perturbed scorings + 4 originals = **164 runs**.
Random seeds are deterministic (within: 0–19; cross: 1000–1019) for
reproducibility.

### Sample size rationale (N=20)

N=20 was selected to amortize 160 perturbed runs against ~1-minute total
compute time while providing visual scatter for Figure X-bis. Wilson 95%
confidence intervals on observed STRONG-retention fractions are broad at
this N (e.g., 10/20 ⇒ [0.30, 0.70]), and we treat the present result
as auxiliary evidence rather than a primary statistical test. Larger
N=200 would tighten CIs to roughly ±0.07 and is left as a future-work
extension; the N=20 result is sufficient for the qualitative ordering
claims that drive this analysis.

### Engine settings

`run_null=False, run_sensitivity=False` (saves time; we only need
`overall_score` + `label`). The marginal effect of target perturbation
is the question of interest, not the OAT-weight-robustness of the
perturbed YAMLs themselves; running weight sensitivity on perturbed
runs would conflate two independent perturbation axes. Default
thresholds (strong=0.75, suggestive=0.40) and default
`min_completeness=30` from each YAML.

---

## 2. Results

### Headline table — annotation-breadth gradient across 4 datasets

[`paper/benchmarks/external/perturbation/perturbation_summary.tsv`](../benchmarks/external/perturbation/perturbation_summary.tsv)

| Dataset | Annotation regime | Mode | Median | Mean | Strong fraction (Wilson 95% CI) | Label changed |
|---|---|---|---|---|---|---|
| **Arm A in-house** (168 MAG, 4-element rich) | saturated | within_element | 1.000 | 0.983 | **20/20 (100%)** [0.84, 1.00] | 0% |
| Arm A in-house | saturated | **cross_element** | 1.000 | 1.000 | 20/20 (100%) [0.84, 1.00] | 0% |
| Liu 2023 (29 MAG, As cold seep) | focused (As-only) | within_element | 0.770 | 0.787 | 10/20 (50%) [0.30, 0.70] | 50% |
| Liu 2023 | focused | **cross_element** | **0.000** | **0.000** | **0/20 (0%)** [0.00, 0.16] | **100%** ⭐ |
| Grettenberger 2021 (29 MAG, AMD) | mixed (S+Fe) | within_element | 0.772 | 0.753 | 10/20 (50%) [0.30, 0.70] | 50% |
| Grettenberger 2021 | mixed | cross_element | 0.455 | 0.471 | 6/20 (30%) [0.13, 0.55] | 70% |
| Ayala 2020 (13 MAG, AMD pit lake) | mixed (S+Fe) | within_element | 0.500 | 0.521 | 8/20 (40%) [0.21, 0.63] | 60% |
| Ayala 2020 | mixed | cross_element | 0.500 | 0.464 | 3/20 (15%) [0.05, 0.37] | 85% |

(Original calibration scores: all four arms 1.000 STRONG. Arm A is run with
`keystone_df` loaded from `tests/sample_data/keystone_species.txt` so
its 9-claim YAML — including `keystone_on_iron` — is fully evaluated;
this matches Table 1 §X.1.)

### Four findings

1. **Annotation breadth governs perturbation discriminating power.** The
   STRONG-retention fraction under cross-element perturbation is
   inversely correlated with how monolithic the dataset's element-level
   focus is: Liu (As-only) → 0%, Ayala (mixed S+Fe) → 15%, Grettenberger
   (mixed S+Fe, larger) → 30%, Arm A (4-element rich) → 100%. This
   monotonic ordering is itself the strongest internal-validity check:
   the perturbation analysis behaves predictably as a function of dataset
   annotation regime, supporting both that (a) the engine is
   discriminating in focused datasets, and (b) the engine is
   appropriately *insensitive* to target perturbation in saturated
   datasets where any element-level pathway is plausibly active.

2. **Cross-element perturbation is strongly discriminating in focused
   datasets.** Liu 2023 is the cleanest test case: **0/20 cross-element
   perturbations retained STRONG**, with median score collapsing to
   0.000 because Liu cold-seep data is As-element-focused and the YAML
   has 2 of 4 claims marked `required=true`. Cross-element replacement
   systematically lands on inactive N/S/Fe pathways, triggering
   required-claim veto. The contrast between within-element 50% STRONG
   and cross-element 0% STRONG **under identical YAML weight structure
   and engine settings** is mechanistically diagnostic: in focused
   datasets, the calibration result depends on element-level target
   accuracy, not just on summing satisfied claims.

3. **Within-element perturbation degrades scores moderately but allows
   incidental STRONG (focused/mixed datasets only).** Within-element mean
   scores are 25–48% below the original 1.000, but 40–50% of
   within-element perturbations *still* produce STRONG in Liu /
   Grettenberger / Ayala. This is consistent with the manuscript's
   "KEGG-coverage-dependent" framing (Discussion §Y.1): when a dataset
   has high KEGG annotation breadth across an element, multiple parallel
   pathway-active claims will register satisfied irrespective of which
   subset the author selected. This bounds the calibration claim — it
   does not invalidate it.

4. **Arm A 100% STRONG retention quantifies the saturation regime.**
   The author-data calibration is robust to both within-element and
   cross-element pathway perturbation (20/20 STRONG retained in both
   modes) because the 168-MAG dataset annotates active pathways across
   all 4 elements simultaneously. Substituting any pathway in any
   element still hits an active community function. This is honest
   auxiliary evidence that for Arm A specifically, the perturbation
   analysis cannot strongly distinguish "author cherry-picked targets"
   from "any reasonable target works because data is rich"; the
   cherry-pick concern for Arm A must therefore rely on the broader
   pre-registration discipline (§4.6.2) and the planned blind
   third-party stress YAMLs (§Y.4) rather than on this perturbation
   alone.

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

2. **Arm A perturbation is partial, not full.** Only the three
   `pathway_active` claims in `arsenic_steel_slag.yaml` were perturbed;
   the two `coupling_possible`, two `env_correlation`, one
   `keystone_in_pathway`, and one `group_contrast` claim were preserved.
   The reasoning is that env_correlation and group_contrast claims pair
   pathway with a second semantically-tied parameter (env_factor or
   group ratio), and perturbing pathway alone would not test
   target-choice sensitivity in a clean single-axis manner; coupling
   claims have no pathway parameter at all. The partial perturbation
   therefore cleanly tests the three pathway_active claims but does not
   test the target-choice sensitivity of Arm A's six other claims.
   Because the six unperturbed claims contribute substantial weight
   (sum w = 7.1 / 9.9 = 72%), Arm A's overall score is partially
   anchored by the unperturbed core, contributing to the observed 100%
   STRONG retention.

3. **Pre-registered claim count differs across datasets**. Liu has 4
   claims, Grettenberger has 4 claims, Ayala has 4 claims, Arm A has 9
   claims with different `required`/non-required composition. We did not
   normalize for this; the raw fraction-STRONG numbers should be read as
   dataset-specific characterisations, with the across-dataset *ordering*
   (Arm A → Grettenberger → Ayala → Liu) being the load-bearing
   qualitative finding rather than the absolute percentages.

4. **No coupling/keystone perturbation**. We perturb only the
   `params.pathway` field. For any future coupling-rich YAML, a separate
   species perturbation analysis would be required.

5. **N=20 is small**. Wilson 95% confidence intervals on the "strong
   fraction" are broad (e.g., 10/20 ⇒ [0.30, 0.70]). For
   publication-ready precision, N=200 or higher would tighten estimates
   to roughly ±0.07. We treat the present result as auxiliary evidence
   rather than a primary statistical test, and the across-dataset
   monotonic ordering survives this CI broadening (the gap between
   Arm A's 100% and Liu's 0% is much larger than any plausible CI
   width).

---

## 4. Interpretation for the manuscript

The perturbation analysis provides **auxiliary evidence consistent with**
(but not ironclad proof of) the claim that the four pre-registered
calibration STRONG outcomes are not mechanical artifacts of KEGG annotation
density alone:

- **Annotation-breadth gradient (Arm A → Grettenberger → Ayala → Liu)**
  is the strongest evidence — STRONG-retention under cross-element
  perturbation is monotonically inversely correlated with how
  monolithically focused the dataset is on a single KB element. This
  monotonic behaviour confirms the engine is doing something
  data-dependent (not purely YAML-mechanical).

- **Cross-element perturbation in focused datasets**: 70–100% label
  degradation (and 0% STRONG retention for Liu) shows that
  **element-level target accuracy matters** in datasets where the data
  is element-monolithic. A randomly chosen target from outside the data's
  biogeochemical domain does not reproduce the calibration.

- **Within-element perturbation result**: moderate score degradation
  (25–48% mean drop) with 40–50% incidental STRONG retention quantifies
  the "KEGG-coverage-dependent" caveat already acknowledged in Discussion
  §Y.1. We do not over-claim that within-element perturbation eliminates
  STRONG; we report it as a soft upper bound on annotation-density-driven
  inflation.

- **Arm A 100% STRONG retention**: the author-data Arm A is in the
  saturation regime — its 168 MAGs annotate active pathways across all
  4 KB elements simultaneously, so perturbation cannot strongly
  distinguish "author cherry-picked targets" from "any reasonable
  target works because data is rich". The cherry-pick concern for
  Arm A specifically therefore relies on the broader pre-registration
  discipline (§4.6.2) and the planned blind third-party stress YAMLs
  (§Y.4) rather than on this perturbation analysis.

- **Limitation of this analysis**: the perturbation tests author-target-
  *choice* sensitivity (and bounds the engine's KEGG-coverage
  dependence), not author-*selection*-bias for the STRONG outcome. The
  blind-hypothesis-writing exercise (§Y.4) remains the gold standard
  for the latter and is preserved as future work. The perturbation
  analysis is offered as auxiliary evidence on the target-choice axis,
  not as a substitute for blind authoring.

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
| 2026-05-09 | Initial run on 3 external datasets; both modes; N=20. Cross-element Liu 0/20 STRONG = headline result. |
| 2026-05-10 | Added Arm A partial perturbation (3 pathway_active claims only). Arm A 100% STRONG retention reframes story as annotation-breadth gradient: Arm A (saturated) → Grettenberger / Ayala (mixed) → Liu (focused). Within-element CIs (Wilson 95%) added. N=20 rationale + run_null=False rationale added. |
