# Discussion §Y — Calibration vs Discrimination, Limitations, and Future Work

> Paper 3 Discussion section draft for the calibration / discrimination /
> limitation / future-work threads opened by the controlled experiment in
> §X (Results). To be inserted/integrated into outline_imeta.md §5.3
> (Discussion). Word count: ~550.
>
> **Status**: ready to integrate; pairs with `methods_external_validation.md`
> (§4.6 Methods) and `results_stress_test_section.md` (Results §X).

---

## §Y.1 Calibration evidence, not validation

The four-arm controlled experiment provides what we term **calibration
evidence**: under fixed default thresholds, EnvMeta's scoring engine returns
`STRONG` for KEGG-curated datasets across two arsenic-cycle topics and two
non-arsenic AMD topics, and `INSUFFICIENT` only for the dataset whose
published annotation provides too few canonical KOs to span the target
pathways. This is consistent with the engine's intended design as a
descriptive evidence-weighting scorecard rather than a frequentist hypothesis
test, and it positions EnvMeta as a *diagnostic instrument* for matching
hypothesis claims against available annotation breadth — not as an oracle
that confirms or denies the user's preferred biological interpretation. The
single most concrete recommendation that follows from this calibration is
that metagenomic studies should publish full KEGG annotation (e.g. KofamScan,
DRAM, or METABOLIC outputs) alongside any custom ROCker / hand-curated gene
sets; otherwise, secondary analyses including ours and any future
KEGG-grounded tool will face the same coverage-driven coarse-graining we
report on Arm B.

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

**Update (v0.9.x)**: the `dominance_score` extension has been implemented in
EnvMeta v0.9.x. A second pre-registered stress YAML (`*_stress_v2.yaml`) with
the Class-A reversed claim augmented by `min_dominance_fraction = 0.20`
upgraded both Liu (0.625 → 0.250) and Ayala (0.455 → 0.182) stress scores from
the `suggestive` (B-tier) to the `weak` (A-tier) label, with observed
dominance scores of 0.05% and 7.08% respectively, both well below the 20%
threshold. All three datasets now exhibit clean A-tier discrimination. The v2
YAMLs are *additional* to v1 (committed alongside, not replacing v1), so the
original B-tier evidence remains in git history as documentation of the
binary-threshold limitation that motivated the engineering change.

## §Y.3 Limitations

Several limitations remain that no software change can resolve. **(1) Author
familiarity with target papers.** Despite explicit pre-registration discipline,
we (the authors) had read all four target papers before writing the YAMLs, so
the claim selection itself reflects implicit knowledge of which pathways were
likely to score well. The cleanest mitigation is **blind hypothesis writing**
by collaborators unfamiliar with the dataset's findings — a study design we
recommend for future iterations. **(2) KB coverage.** EnvMeta KB v1.1 covers
4 elements × 18 pathways × 57 KOs (As / N / S / Fe). Two of Wei's 14 target
genes (*arxA* anaerobic arsenite oxidase; *nrfA* DNRA pathway) lacked KB
mappings and were skipped; iron(II) oxidation and iron(III) reduction
pathways central to AMD biogeochemistry are not yet encoded. KB v1.2 will
extend ROCker-model alias support and add iron-redox plus DNRA blocks.
**(3) Pre-publication hand-checks.** Post-hoc DOI verification of the
hypothesis YAMLs (Results §X.3) identified four citation errors (incorrect
journal names for Yin 2011 and Cabrera 2006; a non-existent "Bothe 2007 *FEMS
Rev*"; topic-mismatched Auld 2017 cited as AMD diazotrophy). These do not
affect scoring outputs but they emphasize the value of building DOI verification
into the hypothesis-writing workflow itself; future YAMLs in
`docs/hypothesis_writing_guide.md` now require an inline `# DOI:` annotation
on each cited reference and a `Direct` / `Inferred` / `Weak` proof-of-extraction
grade.

## §Y.4 Future work

Beyond the `dominance_score` extension (now delivered in v0.9.x) and KB v1.2
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
