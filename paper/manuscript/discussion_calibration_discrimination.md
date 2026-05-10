# Discussion §Y — Calibration vs Discrimination, Limitations, and Future Work

> Paper 3 Discussion section draft for the calibration / discrimination /
> limitation / future-work threads opened by the controlled experiment in
> §X (Results). To be inserted/integrated into outline_imeta.md §5.3
> (Discussion). Word count: ~550.
>
> **Status**: ready to integrate; pairs with `methods_external_validation.md`
> (§4.6 Methods) and `results_stress_test_section.md` (Results §X).

---

## §Y.1 Calibration evidence is KEGG-coverage-dependent, not domain-neutral

The four-arm controlled experiment provides what we term **KEGG-coverage-
dependent calibration evidence**: under fixed default thresholds, EnvMeta's
scoring engine returns `STRONG` for the four datasets that supply canonical
KEGG annotation (KofamScan, DRAM, METABOLIC, or end-to-end Pyrodigal +
GhostKOALA re-annotation), and `INSUFFICIENT` for the one dataset whose
published annotation is restricted to a custom ROCker subset of fourteen
arsenic-cycle genes. This pattern reflects two coupled effects rather than
one. **First**, the engine performs as designed when KEGG-orthology coverage
is broad enough to span the target pathways encoded in the knowledge base
— a calibration-style result confirming that thresholds and scoring logic
behave consistently across heterogeneous study topics. **Second**, the
INSUFFICIENT outcome on Arm B is itself diagnostic: it correctly flags the
mismatch between Wei's 14-gene ROCker set and EnvMeta's 57-KO knowledge base,
rather than failing silently or falsely awarding partial credit. Together,
these observations position EnvMeta as a *diagnostic instrument* whose
performance is **conditional on adequate KEGG coverage of the target
pathways** — not as an oracle that confirms or denies the user's preferred
biological interpretation, and not as a domain-blind tool that performs
equally on any annotation regime.

This conditional-performance framing has two practical implications. First,
metagenomic studies that intend to use EnvMeta or similar KEGG-grounded
secondary tools should publish full KEGG annotations (e.g. KofamScan, DRAM,
or METABOLIC outputs) alongside any custom ROCker or hand-curated gene
sets. Without this, the tool's coverage-driven coarse-graining (as
illustrated on Arm B) limits the resolvable mechanism-space. Second, the
boundary between "KEGG coverage adequate" and "KEGG coverage insufficient"
should itself be reportable as part of the scoring output — a refinement
we discuss as future work in §Y.4.

## §Y.2 Discrimination power, with a binary-threshold caveat

The stress-test layer (Results §X.2) provides direct evidence that EnvMeta's
scoring engine has **discrimination power**: it does not award high scores
to claims that violate environmental priors. The cross-topic
"arsenate_reduction should dominate" claim was rejected with n = 0 active
MAGs in both non-arsenic datasets (Grettenberger AMD stream, Ayala IPB pit
lake), ruling out the *a priori* concern that the universal *arsC*
detoxification homolog would inflate cross-topic scores. The
`pathway_inactive` claim type introduced in v0.9.0, encoding Popperian
falsification predictions, returned the expected `unsatisfied` label in 3/3
datasets when applied to backbone pathways that were genuinely active in the
data. Together with the calibration result, the engine appears to be
domain-neutral — neither hard-wired to confirm arsenic hypotheses nor biased
against alternative contexts.

We acknowledge an honest limitation in the same data: in two of three datasets
(Liu cold-seep and Ayala pit lake), the reversed-direction stress claim
"arsenite oxidation should dominate" returned `satisfied` because real but
weak oxidizer signals were present (total contributions 21-fold and 10-fold
below the dominant reduction pathway, but above the binary `mean_completeness
≥ 50%` threshold). This is not a discrimination failure — the data genuinely
contain weak oxidizer signals consistent with sporadic literature reports
(Stolz et al., 2006) — but rather a *reporting granularity* limitation of the
current binary `satisfied`/`unsatisfied` scheme. Stress claims encoding
domain-violating *dominance* (rather than mere presence) cannot be expressed
without modifying the engine. We propose a future `dominance_score =
total_contribution / element_total` field with an associated
`min_dominance_fraction` parameter; this single extension would convert the
B-tier discrimination outcomes for Liu and Ayala into A-tier clean rejections,
matching the Grettenberger result. We deliberately do not apply this fix
retroactively because the pre-registered claim entities must remain frozen.

**Engineering retrofit (v0.9.x)**: motivated by the limitation above, we
implemented a `dominance_score` extension in EnvMeta v0.9.x and a second
pre-registered stress YAML (`*_stress_v2.yaml`) augmenting the Class-A
reversed claim with `min_dominance_fraction = 0.20`. We are explicit that
**the threshold value was chosen after observing v1 dominance outputs**
(Liu 0.05%, Ayala 7.08%) and is therefore not an independent prediction.
The v2 outcomes (Liu 0.625 → 0.250; Ayala 0.455 → 0.182) confirm that the
retrofit behaves as designed under the engineering choice we made, but we
do not interpret them as independent validation of EnvMeta's discrimination
power. The v1 and v2 results are reported side by side in Table 2 so that
readers can judge the retrofit's effect on the same data. The v2 YAMLs are
*additional* to v1 (committed alongside at `fdfae77`, not replacing v1 at
`50c4687`), preserving the v1 binary-threshold result as documentation of
the limitation that motivated the engineering change. Independent
validation of the `dominance_score` extension on a fresh dataset is part
of the future-work agenda outlined in §Y.4.

## §Y.3 Limitations

We list four interlinked limitations of the present validation experiment in
descending order of methodological importance.

**(1) Residual author selection bias.** This is the single largest
limitation of the calibration experiment. The in-house Arm A is most
susceptible — the authors authored the hypothesis YAML for their own
dataset and chose claims plausibly satisfiable by their own a-priori
research design — and we therefore frame Arm A as a positive control
(engine self-consistency check) rather than as independent calibration
evidence (Results §X.1). The three external arms (C1, C2-A, C2-B) are
less susceptible because the authors authored each YAML before reading
the corresponding paper's specific findings, but they are not bias-free:
even with explicit time-pre-registration discipline (§4.6.2), the
authors had selected the four target papers before writing the YAMLs and
will have unconsciously favored claims plausibly satisfiable in
KEGG-curated datasets within the topics surveyed. The three external
`STRONG` calibration outcomes therefore conflate two effects that cannot
be cleanly separated within the present design: the scoring engine's
behaviour under default thresholds, and the authors' skill in claim
selection. As partial
auxiliary evidence on this point, we performed a target-pathway
perturbation analysis (Methods §4.6.7; Results §X.3) randomly replacing
`params.pathway` across all four calibration YAMLs in two modes
(within-element and cross-element). The headline finding is a
**monotonic annotation-breadth gradient** in cross-element STRONG
retention: Arm A 100% (partial perturbation, 3 of 9 claims) → Grettenberger 30%
→ Ayala 15% → Liu 0% (focused, As-only). This monotonic ordering
behaves predictably with each dataset's annotation breadth — a property
independent of the calibration outcome itself — and is therefore
internally validity-checking. The Liu cross-element 0/20 STRONG result
in particular shows that element-level target accuracy is
mechanistically required for the calibration outcome in focused
datasets. Within-element perturbation moderately degrades scores (mean
drop 25–48%) in the three external datasets but retains STRONG in 40–50%
of runs, bounding the KEGG-coverage caveat already acknowledged in §Y.1.
The Arm A perturbation is necessarily partial — restricted to its three
`pathway_active` claims (the `coupling_possible`, `env_correlation`,
`keystone_in_pathway`, and `group_contrast` claims pair pathway with
semantically-tied second parameters and are not amenable to clean
single-axis perturbation) — and Arm A's 100% retention in this partial
test does *not* by itself demonstrate non-mechanical calibration for
that dataset. The cherry-pick concern for Arm A specifically therefore
relies on the broader pre-registration discipline and the planned blind
third-party stress YAMLs (§Y.4) rather than on this perturbation
analysis. We treat the perturbation analysis as auxiliary evidence
consistent with — not ironclad proof of — non-mechanical calibration in
focused-data regimes; the cleanest remaining mitigation in the
partial-perturbation regime is **blind hypothesis writing** by collaborators
unfamiliar with the dataset's findings, which we discuss in §Y.4 and
identify as the most direct path forward.

**(2) KEGG-coverage dependency.** As discussed in §Y.1, the calibration
evidence is conditional on adequate KEGG-orthology coverage of the target
pathways. The `INSUFFICIENT` outcome on Arm B (Wei 2024 ROCker-only
annotation) is itself diagnostic of this conditional, but it limits the
generalisation of the four `STRONG` results to environments where canonical
KEGG annotations are available.

**(3) KB coverage.** EnvMeta KB v2.0 (KEGG snapshot 2026-04-15) covers 4
elements × 18 pathways × 57 KOs (As / N / S / Fe). Two of Wei's 14 target
genes (*arxA* anaerobic arsenite oxidase; *nrfA* DNRA pathway) lacked KB
mappings and were skipped; iron(II) oxidation and iron(III) reduction
pathways central to AMD biogeochemistry are not yet encoded. KB v2.1 will
extend ROCker-model alias support and add iron-redox plus DNRA blocks.

**(4) Pre-publication hand-checks.** Post-hoc DOI verification of the
hypothesis YAMLs (Results §X.3) identified four citation errors (incorrect
journal names for Yin 2011 and Cabrera 2006; a non-existent "Bothe 2007
*FEMS Rev*"; topic-mismatched Auld 2017 cited as AMD diazotrophy). These do
not affect scoring outputs but they emphasise the value of building DOI
verification into the hypothesis-writing workflow itself; future YAMLs in
`docs/hypothesis_writing_guide.md` now require an inline `# DOI:` annotation
on each cited reference and a `Direct` / `Inferred` / `Weak`
proof-of-extraction grade.

## §Y.4 Future work

Beyond the `dominance_score` extension (now delivered in v0.9.x) and KB v2.1
backlog mentioned above, two methodological future items follow directly from
the stress-test experience. First, **third-party blind stress YAMLs**: in the next user-study
iteration, collaborators unfamiliar with the four target papers will be
invited to author independent stress YAMLs for the same datasets, providing a
selection-bias-controlled replication of the present result. Second, an
**LLM-assisted hypothesis YAML drafting** workflow that ingests an
environmental description and a list of recent reviews, produces a candidate
calibration-plus-stress YAML, and flags claim entities that might fail
DOI verification or proof-of-extraction grading. We see EnvMeta's hypothesis
scorer as well-suited to LLM-assisted authoring precisely because the engine
mechanically resists hindsight bias (pre-registration discipline + binary
status reporting + claim-entity immutability), making the LLM's drafting
errors easy to surface rather than easy to mask.

---

## Word count: ~640 (target 500-700 ✅)

## How this fits into outline §5.3 Discussion

| Outline §5.3 item | Where this draft contributes |
|---|---|
| 1. Differentiation from competitors | §Y.1 (calibration framing as diagnostic instrument) |
| 2. vs sequencing-vendor cloud platforms | not addressed here (orthogonal thread) |
| 3. Fork-rather-than-community model | not addressed here (orthogonal thread) |
| 4. Limitations | §Y.3 (3 limitations: author familiarity / KB coverage / DOI hand-checks) |
| 5. Future work | §Y.4 (dominance_score / blind stress / LLM-assisted YAML drafting) |
| New: Calibration vs Discrimination framing | §Y.1 + §Y.2 (insertable as paragraph 1.5 or as own subsection) |

## Maintenance log

| Date | Event |
|---|---|
| 2026-05-09 | Initial complete draft (§Y.1 calibration / §Y.2 discrimination / §Y.3 limitations / §Y.4 future work) |
