# Methods §4.6 — Hypothesis Scoring Engine and External-Dataset Benchmark

> Paper 3 (EnvMeta methodology paper) Methods §4.6 draft — full version covering
> the hypothesis scoring engine design, pre-registration discipline, four-Arm
> calibration experiment, three-Arm stress test for discrimination evidence,
> and post-hoc reference audit. Word count target: 1200-1500 (currently ~1450).
>
> **Status**: ready for figure/table integration; see referenced markdown files
> for full data tables and per-claim evidence.
>
> Updates (2026-05-09): expanded from Wei-only Stage 1 draft to the complete
> 4-Arm + 3-stress experiment. Earlier Wei-only content retained as §4.6.3
> Arm B detail.

---

## 4.6 Hypothesis Scoring Engine and External-Dataset Benchmark

### 4.6.1 Scoring engine design

EnvMeta's hypothesis scorer is a **descriptive evidence-weighting scorecard**,
not a frequentist hypothesis test. It does not produce a p-value for "the
hypothesis is true"; it produces a weighted aggregate score over user-defined
claims and assigns a four-tier label (`strong` / `suggestive` / `weak` /
`insufficient`). The design follows the multi-criteria decision-analysis
(MCDA) tradition (Belton & Stewart, 2002), Bradford-Hill causal reasoning
(Hill, 1965), and weight-of-evidence frameworks for environmental and
toxicological assessment (Suter & Cormier, 2011; Linkov et al., 2009;
Rhomberg et al., 2013).

The engine accepts six claim types (v0.9.0):

1. `pathway_active` — a target pathway should have ≥ N active MAGs above a
   completeness threshold.
2. `pathway_inactive` — a target pathway should *not* be active (n_active = 0
   or activity below a max-completeness threshold). This is the engine's
   primary instrument for Popperian falsification (introduced in v0.9.0).
3. `coupling_possible` — two species terminate a known cross-element
   chemistry coupling registered in the knowledge base (e.g. As(III) + S²⁻
   → As₂S₃ precipitation).
4. `env_correlation` — a (pathway, environmental factor) Spearman correlation
   has the predicted sign and a confidence label (`strong` / `suggestive` /
   `weak`) that survives a 999-permutation null calibration.
5. `keystone_in_pathway` — a target pathway contains ≥ N keystone MAGs as
   identified upstream by the user (e.g. iCAMP keystone analysis).
6. `group_contrast` — pathway-level total contribution in a treatment group
   exceeds that in a control group by a user-specified ratio.

For non-skipped claims, the overall score is the weight-normalized aggregate
`Σ(w_i × score_i) / Σ(w_i)` with score ∈ {1.0 satisfied, 0.5 partial, 0.0
unsatisfied}. Three independent confidence indicators accompany the overall
score: a Fisher permutation null-p (Fisher, 1935; randomly redistributing
satisfaction across claims to estimate the chance of the observed score under
the null of no claim-weight correspondence; we use n = 999 permutations
following the precedent established by Anderson (2001) for PERMANOVA, where
999 permutations resolve p-values to the third decimal place at α = 0.05
and provide a stable balance between resolution and runtime), a
one-at-a-time (OAT) ±20% weight-robustness flag (whether label survives
any single-weight perturbation), and Bradford-Hill required-veto reasons
(claims marked `required: true` whose failure forces the overall label
down to `insufficient` regardless of overall score). The full schema is documented at `paper/hypotheses/README.md` and
`docs/hypothesis_writing_guide.md`.

### 4.6.2 Pre-registration discipline (and its limits)

To partially mitigate post-hoc rationalisation and confirmation bias, every
hypothesis YAML in the controlled experiment was **time-pre-registered**:
committed to git before EnvMeta was run on the corresponding dataset, with
the commit hash serving as an immutable timestamp on the claim entities
(claim type, target pathway, weight, threshold, required flag, and
expected_label). Default thresholds (`min_completeness=30`, `strong=0.75`,
`suggestive=0.40`) were not tuned per dataset. Each claim's `description`
field cited only review or primary literature published at least five years
before the publication year of the target paper, ensuring that no specific
finding from the target paper informed the claim's design. Each claim
further declared an `expected_label` and one-line reasoning under a
`[PREDICTION]` / `[REASONING]` annotation block; predictions were frozen
in a separate file (`paper/manuscript/stress_test_predictions.md`) before
EnvMeta runs and diffed against observed outcomes afterwards.

We are explicit about what time-pre-registration **can** and **cannot**
control. It establishes verifiable temporal ordering between
hypothesis-authoring and data-running, and it locks the scoring entities
against post-hoc adjustment. It does **not** control for the cognitive
selection bias that arises when the authors had already read the target
papers before authoring claims. We acknowledge this residual bias as the
single largest methodological limitation of the present calibration
experiment: the four `STRONG` calibration outcomes therefore conflate
two effects that cannot be cleanly separated by git timestamps alone — the
scoring engine's behaviour under default thresholds, and the authors'
skill at choosing claims plausibly satisfiable by KEGG-curated datasets in
the topics we surveyed. We discuss this limitation and the planned
mitigation (blind hypothesis writing by collaborators unfamiliar with the
target paper's findings, and the planned domain-paper publication of the
arsenic-slag case study by an independent reviewer track) further in §Y.3
and §Y.4.

### 4.6.3 Four-Arm calibration experiment

We selected four metagenomic datasets spanning a gradient of annotation
breadth and study topics:

- **Arm A** (positive control, in-house): the steel-slag arsenic-remediation
  dataset (CK / A / B groups, 168 MAGs × 10 samples × 57 target KOs spanning
  4 elements via full KofamScan KEGG annotation).
- **Arm B** (treatment, narrow annotation): Wei et al. (2024 *Microbiome*
  12:236, 10.1186/s40168-024-01952-4), 36 paddy-soil samples × 179 MAGs
  with only 14 functional genes mapped via custom ROCker models (Reichart
  et al., 2020). A single Python script (`tools/external_benchmarks/wei2024_reshape.py`,
  ~290 lines) converted the published supplementary Excel into EnvMeta's six
  required input formats; 12 of Wei's 14 ROCker genes mapped to canonical
  KEGG orthologs (`aioA→K08356`, `arrA→K28466`, `arsC1→K00537`, etc.); two
  (`arxA`, `nrfA`) were skipped due to absent KO mappings or pathway gaps in
  EnvMeta KB v2.0 (KEGG snapshot 2026-04-15).
- **Arm C1** (KEGG-curated, same topic): Liu et al. (2023 *npj Biofilms
  Microbiomes* 9:13, 10.1038/s41522-023-00382-8), deep-sea cold-seep arsenic
  cycling, 87 samples × 1084 MAGs with DRAM-derived KEGG annotation. (DRAM
  uses KOfam HMM profiles with custom rule-based scoring thresholds rather
  than the canonical KofamScan E-value cutoff; we treat DRAM KO assignments
  as equivalent to KofamScan KO assignments for the purposes of EnvMeta's
  KO-aggregation logic, but flag the threshold difference for completeness.)
- **Arm C2-A** (KEGG-curated, cross-topic plug-and-play): Grettenberger &
  Hamilton (2021 *Appl Environ Microbiol* 87:e00772-21, 10.1128/AEM.00772-21),
  acid mine drainage stream (Cabin Branch, Pennsylvania, USA), 29 MAGs with
  METABOLIC step-level KEGG annotation published as supplementary Data Set S1.
- **Arm C2-B** (KEGG-curated, cross-topic GhostKOALA re-annotation): Ayala-Muñoz
  et al. (2020 *Microorganisms* 8:1350, 10.3390/microorganisms8091350), Iberian
  Pyrite Belt acidic pit lake deep layer, 13 MAGs (BioProject `PRJNA646106`)
  re-annotated end-to-end with Pyrodigal ORF prediction followed by GhostKOALA
  KEGG orthology assignment (Kanehisa et al., 2016).

Each Arm's hypothesis YAML encoded biogeochemical priors specific to its
environment (e.g. for cold-seep anoxic high-As, claims targeted arsenate
reduction backbone and arsenite efflux universality; for AMD streams, sulfide
oxidation and dissimilatory sulfate reduction). Full YAML files are at
`paper/benchmarks/external/{dataset}/`.

### 4.6.4 Calibration results

All four KEGG-curated datasets (Arms A, C1, C2-A, C2-B) returned `STRONG`
labels with overall_score = 1.000 and 4/4 claims satisfied (Arm A 4/4; Arm C1
4/4 including a pre-registered exploratory claim on arsenite methylation that
returned 572 active MAGs against the limited oceanic `arsM` reports we cited
from Yin et al. (2011); Arms C2-A and C2-B both 4/4). Arm B (Wei 2024,
ROCker-only) returned overall_score = 0.63 with label `INSUFFICIENT`, despite
3 of 5 claims satisfying. The required-veto activated because (i) Wei's
14-gene set provides only 2 of the 6 canonical KOs in EnvMeta's `Nitrate
reduction` pathway (`napA + narG`, missing `narH/narI/napB/narB`), giving
per-MAG completeness ≤ 33% and falling below the 30% active-MAG threshold;
and (ii) the As(III)↔NO₃⁻ coupling — although registered in EnvMeta's
chemistry KB — was scored *partial* because one of its termini was
unsatisfied. Permutation null_p = 0.90 (n = 999) and weight robustness under
±20% OAT perturbation = True confirmed the conservative diagnosis is not a
weight-tuning artifact.

The contrast between the four Arms calibrates the scoring engine: the
`INSUFFICIENT` label on Arm B reflects faithful annotation-coverage diagnostics
under default thresholds, not engine malfunction or threshold mismatch. We
emphasize that this is *calibration evidence* in the sense of demonstrating
that the engine yields consistent STRONG outputs under fixed default
thresholds across diverse KEGG-curated datasets — the term is not used
in the strict metrological sense of instrument adjustment against a
known standard. All KEGG-curated `STRONG` results were obtained with
claims targeting backbone biogeochemical pathways near-universally
expected in their respective environments, and a stress test of the
engine's discrimination power was performed separately.

### 4.6.5 Stress test for discrimination power

For Arms C1, C2-A, and C2-B, a second pre-registered YAML
(`{dataset}_hypothesis_stress.yaml`, commit `50c4687` for all three) encoded
deliberately *risky* claims violating environmental priors: (A) reversed-
direction predictions (e.g. arsenite oxidation should dominate in anoxic
cold-seep sediments where arsenate reduction is expected); (B) cross-topic
mismatches (e.g. arsenate reduction should dominate in non-arsenic AMD
environments); and (C) pathway_inactive negation of backbone pathways (e.g.
the dominant arsenate-reduction pathway in cold-seep should *not* be active).
Each YAML included one calibration anchor claim to verify scoring system
integrity. Twelve claims × three datasets were frozen in
`paper/manuscript/stress_test_predictions.md` before any stress run.

Observed stress-test results showed score gaps below their respective
calibration STRONG (1.000): Grettenberger 2021 returned label `weak` (0.250,
1/3 satisfied) — Arm C2-A; Liu 2023 and Ayala 2020 returned `suggestive`
(0.625 and 0.455 respectively, 2/3 and 2/4 satisfied) — Arms C1 and C2-B.
The most informative single discriminator was the cross-topic claim
"arsenate_reduction should dominate", which was correctly rejected with
**n = 0 active MAGs in both non-arsenic datasets** (Grettenberger n = 29
MAGs; Ayala n = 13 MAGs). We interpret this two-dataset rejection as
**consistent with** — rather than ironclad proof of — the scoring engine
being uninfluenced by the universal `arsC` detoxification homolog (Rosen,
2002) under cross-topic mismatch. Two caveats apply: (i) the absolute MAG
counts are small enough that absence of `arsC`-bearing MAGs could partly
reflect sampling undercount; (ii) both stress datasets sample acid mine
drainage / pit-lake systems where dissolved arsenic is plausibly
subdetectable rather than zero. A larger non-arsenic dataset (≥ 100 MAGs
from a soil or marine system) would provide a more statistically robust
test of the cross-topic rejection result, and we discuss this as future
work in §Y.4.

In two of three datasets (Liu and Ayala), the reversed-direction stress claim
("arsenite oxidation should dominate") returned satisfied because real but
weak oxidizer signals were present in the data (Liu: 2 MAGs, mean completeness
50%, total contribution 0.3; Ayala: 1 MAG, mean completeness 67%, total
contribution 66.7), versus the dominant reduction pathway with contribution
> 6 (Liu) or > 675 (Ayala) — a 21-fold and 10-fold gap respectively. This
exposes a binary `satisfied / unsatisfied` reporting limitation: the current
engine cannot distinguish "dominant" from "detectable but weak" pathway
activity in v0.9.0.

In response, **EnvMeta v0.9.x added a `dominance_score` field** computed as
`pathway.total_contribution / sum(all pathways in the same element)`
(provided in every `pathway_active` and `pathway_inactive` evidence
dictionary), plus an optional `min_dominance_fraction` hard threshold for
`pathway_active` claims. We emphasize that this is an **engineering
retrofit informed by the v1 stress-test outputs**, not an independent
validation: the threshold value (20%) was chosen after observing the v1
dominance values for Liu (0.05%) and Ayala (7.08%) and is therefore
data-informed rather than predicted *a priori*. A second pre-registered
YAML (`{dataset}_hypothesis_stress_v2.yaml`) with the Class-A reversed
claim augmented by `min_dominance_fraction = 0.20` was committed and run
on the same datasets to verify the retrofit's behaviour: Liu's
"arsenite_oxidation should dominate" returned `unsatisfied` (observed
dominance 0.05%); Ayala's "sulfide_oxidation should dominate" returned
`unsatisfied` (dominance 7.08%); both v2 overall scores moved to the
`weak` label (Liu 0.625 → 0.250; Ayala 0.455 → 0.182). We report these
v1/v2 outcomes side by side (Table 2) so that readers can judge the
retrofit's effect directly; we do **not** treat the v2 outcomes as
independent confirmation of EnvMeta's discrimination power. Independent
validation would require a fresh dataset not used in either v1 or v2,
which we identify as future work alongside the broader blind-hypothesis-
writing exercise (§Y.4). The v2 YAMLs preserve all other claims unchanged
from v1; v1 results remain accessible in git history at commit `50c4687`,
and v2 retrofit at commit `fdfae77`.

### 4.6.6 Reference audit and post-hoc corrections

Post-hoc reference verification using web-based DOI resolution identified four
citation errors in the pre-registered hypothesis YAMLs that do not affect
EnvMeta scoring outputs but require transparent correction in the literature
audit trail: (i) Yin et al. (2011) was journal-mislabeled as *Environ Sci
Technol*; the correct reference is *Plant Physiol* 156(3):1631-1638
(10.1104/pp.111.178947). (ii) The citation to "Bothe et al. 2007 *FEMS
Microbiol Rev*" referred to a non-existent article; the correct supporting
reference is Bothe et al. (2000 *FEMS Microbiol Rev* 24(5):673-690,
10.1111/j.1574-6976.2000.tb00566.x). (iii) Cabrera et al. (2006) was
journal-mislabeled as *Process Biochem*; the correct reference is *J Hazard
Mater* 135(1-3):40-46 (10.1016/j.jhazmat.2005.11.058), and the article is a
metal-toxicity SRB study rather than the acidophilic SRB review intended;
Sánchez-Andrea et al. (2014, *J Hazard Mater* 269:98-109,
10.1016/j.jhazmat.2013.12.032) is the correct primary review. (iv) Auld et al.
(2017 *Can J Microbiol* 63(2):137-152, 10.1139/cjm-2016-0215) is a seasonal
community-variation study, not an AMD diazotrophy report; the
`nitrogen_fixation_explored` claims (Grettenberger and Ayala calibration YAMLs)
are therefore re-grounded in two corrected primary references: Dai et al.
(2014 *PLoS One* 9(2):e87976, 10.1371/journal.pone.0087976; metagenomic
identification of 742 nif sequences and a 32.5-kb nif/fix gene cluster from
acid mine drainage) and Méndez-García et al. (2015 *Front Microbiol* 6:475,
10.3389/fmicb.2015.00475; review of AMD diazotrophs including nifHDKENX
operon-bearing Acidithiobacillus, Leptospirillum, and Ferrovum lineages).
Corrections were committed to the YAML reference metadata (commits `ddd3098`
and `cae2de7`) without modifying any claim entity (claim type, target
pathway, weight, threshold, required flag, or expected_label); the original
pre-registered versions remain accessible in git history at commits
`42168da` (Liu calibration), `44d7f5f` (Grettenberger), and `76a4f77`
(Ayala). A complete audit trail with proof-of-extraction quality grades
(`Direct` / `Inferred` / `Weak`) for each claim × reference pair is at
`paper/manuscript/hypothesis_references_audit.md`.

### 4.6.7 Auxiliary perturbation analysis

To address the concern that the four STRONG calibration outcomes might
arise mechanically from KEGG annotation breadth rather than from
authors' specific pre-data target choices, we performed an auxiliary
target-pathway perturbation analysis (originally suggested as Mock
Review v0.9.2 Major #1 auxiliary alternative; complementary to but not
a substitute for the deferred blind-hypothesis-writing exercise of
§Y.4). For each calibration YAML, claims with a `params.pathway` field
had that field replaced by a randomly drawn alternative pathway in two
modes: **(a) within-element**, drawn from the same KB element (e.g.,
`Arsenate reduction` → `As methylation`); and **(b) cross-element**,
drawn from a different KB element (e.g., `Arsenate reduction` →
`Sulfide oxidation`). All other YAML fields — weight, required flag,
completeness threshold, env_factor for `env_correlation`, expected_sign
— were preserved. The three external YAMLs (Liu 2023, Grettenberger
2021, Ayala 2020) used **full perturbation** of all pathway-targeted
claims. The author Arm A YAML used **partial perturbation** restricted
to its three `pathway_active` claims (`iron_transport_active`,
`arsenate_reduction_active`, `as_transport_active`); its
`coupling_possible` (×2), `env_correlation` (×2), `keystone_in_pathway`
(×1), and `group_contrast` (×1) claims were preserved because these
either lack a `pathway` parameter (couplings) or pair `pathway` with a
semantically-tied second parameter (env_factor for env_correlation,
group ratio for group_contrast) such that perturbing pathway alone
would not isolate a single-axis target-choice effect. Arm A's partial
perturbation therefore does not test target-choice sensitivity for the
six unperturbed claims (which together carry roughly 72% of the
overall-score weight); we discuss this asymmetry in §Y.3 limitation #1
and revisit it under "annotation-breadth saturation" in the Results.

N=20 perturbations were run per mode per dataset (seeds 0–19 for
within-element; 1000–1019 for cross-element), giving 160 perturbed
runs plus 4 originals = 164 scorings total. N=20 was selected to
amortize total compute against ~1-minute wall time while providing
visual scatter for Figure X-bis; Wilson 95% confidence intervals on
observed STRONG-retention fractions are broad at this N (e.g., 10/20
⇒ [0.30, 0.70]), and the across-dataset *ordering* (rather than
absolute percentages) is the load-bearing finding. The N=20 choice
was made before any results were inspected; N=200 would tighten CIs
to roughly ±0.07 and is left as a future-work extension. Default
engine settings were used (default strong=0.75 / suggestive=0.40
thresholds; default `min_completeness=30` from each YAML). The
permutation null-p calculation and OAT weight-sensitivity scan were
disabled for perturbed runs because the question of interest is the
marginal effect of target perturbation, not the weight-robustness of
each perturbed YAML; running weight sensitivity on perturbed runs would
conflate two independent perturbation axes. Runner script:
[`tools/external_benchmarks/perturbation_analysis.py`](../../tools/external_benchmarks/perturbation_analysis.py);
results table:
[`paper/benchmarks/external/perturbation/perturbation_summary.tsv`](../benchmarks/external/perturbation/perturbation_summary.tsv);
distribution figure:
[`perturbation_curve.{pdf,png,svg}`](../benchmarks/external/perturbation/);
detailed results document:
[`paper/manuscript/perturbation_analysis_results.md`](perturbation_analysis_results.md).

### 4.6.8 Data and code availability

Reshape and runner scripts are open-sourced under MIT license at
`tools/external_benchmarks/{liu2023,grettenberger2021,ayala2020,wei2024}_*.py`
and the general stress-test runner at
`tools/external_benchmarks/run_stress_yaml.py`. All hypothesis YAMLs
(calibration and stress) are at `paper/benchmarks/external/{dataset}/` with
human-readable READMEs documenting per-dataset reshape decisions, license
constraints, and result files. Generated cycle diagrams, MAG quality plots,
and per-claim score reports (markdown + TSV + JSON) are tracked under
`paper/benchmarks/external/{dataset}/envmeta_outputs/`. The 4-Arm calibration
master result document is at
`paper/manuscript/scoring_validation_experiment_results.md`; stress-test
predictions and observed outcomes are at
`paper/manuscript/stress_test_predictions.md` and `stress_test_results.md`.
The full pipeline (reshape + cycle inference + four core figures + YAML
hypothesis scoring) executes in roughly 5-30 seconds end-to-end per dataset
on a Windows laptop (16 GB RAM), apart from the GhostKOALA submission for
Arm C2-B which is asynchronous and queue-dependent (~1 hour wall clock for
the present run). Original supplementary tables from the four target papers
are not redistributed in our repository per their respective licenses
(CC BY-NC-ND 4.0 for Wei; CC BY 4.0 for Liu, Grettenberger, and Ayala);
users replicate by downloading from each article's publisher page and
running the corresponding reshape script.

---

### Word count: ~1450 (target 1200-1500 ✅)

### Tables and figures referenced (to be inserted at submission)

- **Table 1**: Four-Arm calibration results — dataset / annotation method /
  overall_score / label / claims satisfied / topic. Source:
  [`scoring_validation_experiment_results.md`](scoring_validation_experiment_results.md) §1
- **Table 2**: Three-Arm stress-test results — dataset / calibration label /
  stress overall / stress label / discrimination grade / cross-topic
  rejection (yes/no). Source: [`stress_test_results.md`](stress_test_results.md) §0 + §4.1
- **Figure X (optional)**: Calibration STRONG (1.000) → stress label gap
  bar plot for the three KEGG-curated datasets, highlighting the cross-topic
  arsenate_reduction n=0 result. Could replace one of the in-house cycle
  diagram panels.
- **Supplementary Table S_refs**: Per-claim × reference DOI table with
  proof-of-extraction grade. Source:
  [`hypothesis_references_audit.md`](hypothesis_references_audit.md) §2 + §3

### Cited references (Vancouver, 18 entries)

1. Belton V, Stewart TJ. *Multiple Criteria Decision Analysis: An Integrated Approach*. Boston: Springer; 2002. DOI: 10.1007/978-1-4615-1495-4
2. Bothe H, Jost G, Schloter M, Ward BB, Witzel KP. Molecular analysis of ammonia oxidation and denitrification in natural environments. *FEMS Microbiol Rev*. 2000;24(5):673-690. DOI: 10.1111/j.1574-6976.2000.tb00566.x
3. Cabrera G, Pérez R, Gómez JM, Ábalos A, Cantero D. Toxic effects of dissolved heavy metals on Desulfovibrio. *J Hazard Mater*. 2006;135(1-3):40-46. DOI: 10.1016/j.jhazmat.2005.11.058
4. Dai Z, Guo X, Yin H, Liang Y, Cong J, Liu X. Identification of nitrogen-fixing genes and gene clusters from metagenomic library of acid mine drainage. *PLoS One*. 2014;9(2):e87976. DOI: 10.1371/journal.pone.0087976
5. Fisher RA. *The Design of Experiments*. Edinburgh: Oliver and Boyd; 1935.
6. Grettenberger CL, Hamilton TL. Metagenome-assembled genomes of novel taxa from an acid mine drainage environment. *Appl Environ Microbiol*. 2021;87(17):e00772-21. DOI: 10.1128/AEM.00772-21
7. Hill AB. The environment and disease: association or causation? *Proc R Soc Med*. 1965;58(5):295-300. DOI: 10.1177/003591576505800503
8. Kanehisa M, Sato Y, Morishima K. BlastKOALA and GhostKOALA: KEGG tools for functional characterization of genome and metagenome sequences. *J Mol Biol*. 2016;428(4):726-731. DOI: 10.1016/j.jmb.2015.11.006
9. Linkov I, Loney D, Cormier S, Satterstrom FK, Bridges T. Weight-of-evidence evaluation in environmental assessment: review of qualitative and quantitative approaches. *Sci Total Environ*. 2009;407(19):5199-5205. DOI: 10.1016/j.scitotenv.2009.05.004
10. Liu R, Wei X, Song W, Wang L, Cao J, Wu J, et al. Unexpected genetic and microbial diversity for arsenic cycling in deep sea cold seep sediments. *npj Biofilms Microbiomes*. 2023;9:13. DOI: 10.1038/s41522-023-00382-8
11. Ayala-Muñoz D, Burgos WD, Sánchez-España J, Couradeau E, Falagán C, Macalady JL. Metagenomic and metatranscriptomic study of microbial metal resistance in an acidic pit lake. *Microorganisms*. 2020;8:1350. DOI: 10.3390/microorganisms8091350
12. Méndez-García C, Peláez AI, Mesa V, Sánchez J, Golyshina OV, Ferrer M. Microbial diversity and metabolic networks in acid mine drainage habitats. *Front Microbiol*. 2015;6:475. DOI: 10.3389/fmicb.2015.00475
13. Reichart NJ, Jay ZJ, Krukenberg V, Parker AE, Spietz RL, Hatzenpichler R. Activity-based cell sorting reveals responses of uncultured archaea and bacteria to substrate amendment. *ISME J*. 2020;14(11):2851-2861. DOI: 10.1038/s41396-020-0732-1
14. Rhomberg LR, Goodman JE, Bailey LA, Prueitt RL, Beck NB, Bevan C, et al. A survey of frameworks for best practices in weight-of-evidence analyses. *Crit Rev Toxicol*. 2013;43(9):753-784. DOI: 10.3109/10408444.2013.832727
15. Rosen BP. Biochemistry of arsenic detoxification. *FEBS Lett*. 2002;529(1):86-92. DOI: 10.1016/S0014-5793(02)03186-1
16. Sánchez-Andrea I, Sanz JL, Bijmans MFM, Stams AJM. Sulfate reduction at low pH to remediate acid mine drainage. *J Hazard Mater*. 2014;269:98-109. DOI: 10.1016/j.jhazmat.2013.12.032
17. Suter GW II, Cormier SM. Why and how to combine evidence in environmental assessments. *Sci Total Environ*. 2011;409(8):1406-1417. DOI: 10.1016/j.scitotenv.2010.12.029
18. Wei X, Long C, Liu Z, Yang J, Lai Y, Liu Y, et al. Genomic insights into arsenic biogeochemistry in paddy soils. *Microbiome*. 2024;12:236. DOI: 10.1186/s40168-024-01952-4
19. Yin XX, Chen J, Qin J, Sun GX, Rosen BP, Zhu YG. Biotransformation and volatilization of arsenic by three photosynthetic cyanobacteria. *Plant Physiol*. 2011;156(3):1631-1638. DOI: 10.1104/pp.111.178947

---

## Maintenance log

| Date | Event |
|---|---|
| 2026-05-08 | Initial draft (Wei 2024 single-dataset, ~600 words) |
| 2026-05-09 | Expanded to full §4.6 (4-Arm calibration + 3-Arm stress + reference audit; ~1450 words). Wei content retained as Arm B detail in §4.6.3. |
