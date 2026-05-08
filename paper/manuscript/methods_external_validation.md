# Methods — External Dataset Validation

> 此段落为 Paper 3 (EnvMeta 方法学论文) Methods 章节"外部数据集验证"小节
> 的草稿。最终投稿前需根据 outline_imeta.md 第 4.6 节调整段落结构。

---

## 4.6 External-dataset benchmark

To assess whether EnvMeta generalizes beyond the in-house arsenic-slag /
steel-slag remediation dataset (CK / A / B groups, 168 MAGs × 10 samples),
we re-analyzed an independent published dataset focused on a related but
distinct system: arsenic-contaminated paddy soils across South China.

### 4.6.1 Dataset

We selected Wei et al. (2024, *Microbiome* 12:236, doi:10.1186/s40168-024-01952-4),
a 36-sample metagenomic study (12 sites × 3 replicates, pH 4.6-8.0) that
recovered 179 medium- and high-quality MAGs (completeness > 50%,
contamination < 10%) and reported coupling between As(III) oxidation and
denitrification genes. The study published its processed data as
`Additional file 2` (Excel workbook, CC BY-NC-ND 4.0), containing sample
metadata (Table S1), soil physicochemical parameters (Table S2), MAG
basic info with relative abundances (Table S7), and 14 target functional
genes mapped to MAGs by ROCker custom models (Table S8). Raw FASTQ data
in NCBI BioProject `PRJNA1068274` were not used; all EnvMeta inputs were
derived from the published supplementary tables.

### 4.6.2 Re-shaping pipeline

A single Python script (`tools/external_benchmarks/wei2024_reshape.py`,
~290 lines) converted the supplementary Excel into EnvMeta's six required
input formats:

- `metadata.tsv`: 36 samples × {SampleID, Group, Replicate}; Group derived
  from Wei's *Arsenic level* column (AsContam vs NoContam, 18 + 18)
- `env_factors.tsv`: 16 physicochemical parameters from Wei Table S2
  (pH, water content, eight metals/metalloids in mg kg⁻¹, five totals in g kg⁻¹)
- `mag_taxonomy_labels.tsv`: 179 MAGs × concatenated GTDB-Tk lineage strings
- `quality_report.tsv`: CheckM2-style schema populated from Wei Table S7
  completeness/contamination columns
- `abundance.tsv`: 179 MAGs × 36 samples; because Wei published only
  group-level mean abundances rather than per-sample mapping coverage, we
  replicated the AsContam / NoContam means across samples within group
  (a known dimensionality-degenerate substitute that disables PCoA, RDA,
  and per-sample ordination but permits group-aware figures)
- `kegg_target_only.tsv`: 540 long-format records mapping Wei's 14
  ROCker-detected genes to canonical KEGG orthologs (12 mapped: aioA→K08356,
  arrA→K28466, arsC1→K00537, arsC2→K03741, napA→K02567, narG→K00370,
  nirK→K00368, nirS→K15864, norB→K04561, nosZ→K00376, pmoA→K10944,
  pmoB→K10945; 2 unmapped: arxA — anaerobic arsenite oxidase ROCker
  custom model with no canonical KEGG KO; nrfA — DNRA pathway, not
  encoded in EnvMeta KB v1.1).

The full pipeline (re-shape + cycle inference + four core figures + YAML
hypothesis scoring) executed in **~24 s** end-to-end on a Windows laptop
(16 GB RAM), confirming the design target of 30-second cycle-figure
generation.

### 4.6.3 Reproduced figures

EnvMeta automatically rendered five publication-grade figures in addition
to the cycle diagram. Figure 3 (gene heatmap, 12 KO × group means with
z-score normalization) qualitatively reproduces Wei's Fig. S10–S11
gene-abundance comparisons. Figure 4 (log₂ fold-change AsContam vs
NoContam) parallels Wei's Fig. S7 / S14. The biogeochemical-cycle figure
(Figure 1) shows three As pathways and three N pathways active; sulfur
and iron quadrants remain empty because Wei's selected 14 genes do not
cover those elements (a data-availability artifact, not a tool failure).

### 4.6.4 Hypothesis scoring (key result)

A YAML hypothesis was authored encoding Wei's central claim: that
As(III)-oxidizing populations co-occur with denitrification machinery,
suggesting nitrate-mediated As(III) oxidation. Five claims captured the
component evidence (`pathway_active` for arsenite oxidation,
nitrate reduction, nitrite reduction, N₂O reduction; `coupling_possible`
for As(III) ↔ NO₃⁻), with required veto applied to the three Bradford
Hill "necessary conditions" (As(III) oxidation, nitrate reduction, and
the As(III)–NO₃⁻ coupling).

EnvMeta returned **overall_score = 0.63 with label INSUFFICIENT** despite
satisfying 3 of 5 claims. The required veto was triggered because (i)
Wei's published 14-gene set provides only two of the six canonical KOs
in EnvMeta's `Nitrate reduction` pathway (napA + narG, missing narH/narI/
napB/narB), giving per-MAG completeness ≤ 33 % and falling below the
30 % active-MAG threshold for several MAGs; and (ii) the As(III)↔NO₃⁻
coupling — although recognized in EnvMeta's chemistry KB as a redox
coupling — was scored *partial* because one of its termini (the nitrate-
reduction pathway) was unsatisfied. Permutation null_p = 0.90 (n = 999)
and weight robustness under ±20 % OAT perturbation = True confirmed that
the conservative diagnosis is not an artifact of weight tuning.

This outcome exemplifies EnvMeta's "non-endorsement" design principle:
the scorecard does not auto-confirm an investigator's preferred
interpretation. It surfaces a specific, actionable gap — *"the literature
hypothesis is plausible but the published gene set does not span enough
of the canonical pathway for EnvMeta to assign 'strong' support"* —
leaving causal interpretation with the researcher and pointing toward a
concrete remediation (extend the supplementary annotation to additional
narH / narI / napB / narB orthologs).

### 4.6.5 Limitations and KB-coverage feedback

Two of Wei's 14 target genes lack mappings in EnvMeta KB v1.1. **arxA**
is a ROCker custom model without a canonical KEGG KO and was skipped
(17 MAGs / 18 gene copies). **nrfA** belongs to the DNRA pathway, which
v1.1 does not encode (33 MAGs / 33 gene copies). These omissions are
recorded in `tools/external_benchmarks/wei2024_reshape.py:GENE_TO_KO`
with `None` values; KB v1.2 will incorporate a DNRA pathway block and
ROCker-model alias support to close these gaps.

The dimensionality-degenerate `abundance.tsv` (group-mean replication)
prevents PCoA / RDA / α-boxplot / LEfSe / co-occurrence-network analyses,
which is a property of Wei's published data layer rather than EnvMeta;
authors who can release per-sample MAG mapping coverage as standard
supplementary tables will enable downstream tools — including EnvMeta —
to provide the full analytic suite.

### 4.6.6 Data and code availability

Re-shape and runner scripts are open-sourced under MIT license at
`tools/external_benchmarks/wei2024_*.py`. The hypothesis YAML is at
`paper/benchmarks/external/wei_2024_paddy/wei2024_hypothesis.yaml`.
Generated figures (PNG + PDF, 600 dpi) are tracked in
`paper/benchmarks/external/wei_2024_paddy/envmeta_outputs/`.
Wei's original supplementary Excel is **not** redistributed in our
repository per its CC BY-NC-ND license; users replicate by downloading
from the *Microbiome* article page and running our reshape script.

### Word count: ~600 (target 300-500 — needs trimming for final version)
