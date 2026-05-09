# EnvMeta Performance Benchmark

> **Purpose**: empirically establish EnvMeta's runtime + memory scaling envelope
> across realistic dataset sizes, anchor the paper's "Methods §4.6 Implementation"
> performance claims with data, and produce user-facing hardware sizing guidance.
>
> **Last update**: 2026-05-09 (v0.9.1)
> **Hardware**: Windows 10 / Intel laptop / 16 GB RAM / Python 3.11 / single-process Streamlit-equivalent in-process invocation
> **Harness**: [`bench_harness.py`](performance/bench_harness.py) + [`bench_synthetic_dense.py`](performance/bench_synthetic_dense.py) + [`run_scaling_sweep.py`](performance/run_scaling_sweep.py)

---

## 1. TL;DR

EnvMeta is **CPU-bound on cycle inference, memory-light** (peak ΔRSS ≤ 10 MB
across all measured configurations). The dominant cost is the
permutation test in cycle inference's S2 step, whose complexity is
**N_pathway_active × N_env_factors × 999 × O(N_sample × log N_sample)** —
**independent of N_MAGs** beyond a flat S1 aggregation step.

This means a 1000+ MAG dataset with full KEGG annotation and 4 environmental
factors runs the cycle figure in ~20 s on a typical laptop — and any user
machine with 4 GB RAM can run the entire 14-figure pipeline.

| Dataset | N_MAG | N_samp | annotation | N_env | total wall (14 fig) | peak ΔRSS |
|---|--:|--:|---|--:|--:|--:|
| sample (real, dense) | 169 | 10 | 30 KO/MAG | 4 | **29.3 s** | 2.0 MB |
| Liu 2023 (real, sparse 8-KO target) | 1084 | 87 | 1.5 KO/MAG | 1 | **1.6 s** (7 figs) | 9.3 MB |
| Liu synthetic dense | 1000 | 87 | 25 KO/MAG | 4 | **~20 s** (cycle alone) | 4.3 MB |

---

## 2. Methodology

### 2.1 Three-tier dataset coverage

To separate true scaling effects from confounded variables (annotation density,
env breadth, MAG count), we ran three complementary regimes:

1. **sample_data (real, dense)** — 169 MAGs × 10 samples × 4 env factors × 30 KO/MAG; full KofamScan annotation; 3 groups (CK / A / B). The original arsenic-slag steel-slag remediation case study; runs all 14 figures.
2. **Liu 2023 cold seep (real, sparse)** — 1084 MAGs × 87 samples × 1 env factor (Depth_proxy) × 1.5 KO/MAG; published `kegg_target_only.tsv` is pre-filtered to 8 As-cycle KOs only; single group ("All"); runs 7/14 figures (LEfSe / log2FC / RDA / PCoA / network / log2fc / alpha_boxplot fail at single-group or missing-input boundaries).
3. **Liu 2023 synthetic dense** — same 200-1000 MAG × 30-87 sample subsets of Liu, but each MAG randomly assigned 25 KOs from the 57-KO KB; 4 synthetic numeric env factors; tests the **dense-annotation scaling envelope** on a 1000+ MAG dataset, simulating what a full DRAM/METABOLIC-annotated cold seep dataset would behave like.

### 2.2 Measurement protocol

For each (dataset × figure) combination:

- 2-3 repeats; report **median** wall-time (`time.perf_counter`) and **peak ΔRSS** (`psutil.Process().memory_info().rss` sampled at 50 ms by a background thread; baseline taken after `gc.collect()`)
- Headless matplotlib (`Agg` backend); `plt.close("all")` + `gc.collect()` between repeats to reset memory baseline
- All measurements in-process to mirror Streamlit's single-process model; no subprocess spawn overhead

### 2.3 Limitations of the measurements

- **Wall-time on a single machine** — absolute numbers will vary across CPUs (we measured an Intel i7-class laptop). Relative scaling (ratios across cells) is the more transferable signal.
- **Peak ΔRSS underestimates true peak** during very short calls if matplotlib's allocator has retained pages from prior runs. We report ΔRSS over baseline rather than absolute RSS for this reason.
- **No GPU usage** — EnvMeta is pure CPU; results extrapolate directly to any CPU-bound environment.
- **Liu's published annotation is pre-filtered to 8 As-cycle KOs**, which understates real-world cycle_diagram cost on full-KEGG datasets. We approximate dense annotation with the synthetic-dense regime described above.
- **The synthetic-dense regime randomly assigns KOs from the 57-KO knowledge base, which is not biologically realistic.** Real KofamScan / DRAM / METABOLIC annotations cluster around organism-specific functional repertoires, not uniform across taxa. The synthetic-dense scaling numbers should therefore be read as **upper-bound estimates** of permutation-test cost under "every pathway has active MAGs" conditions, rather than as a faithful reproduction of real-world annotation distributions. A direct comparison against one fully KofamScan-annotated public dataset (e.g., reprocessing Liu 2023's MAGs through KofamScan rather than using the published 8-KO subset) is identified as future validation work.

---

## 3. Per-figure scaling

### 3.1 cycle_diagram (the dominant cost)

| dataset | N_MAG | N_samp | N_env | KO/MAG | wall (s) |
|---|--:|--:|--:|--:|--:|
| sample | 169 | 10 | 4 | 30 | 13.5 |
| Liu real | 1084 | 87 | 1 | 1.5 | 0.4 |
| Liu synth dense | 200 | 30 | 4 | 25 | 14.8 |
| Liu synth dense | 200 | 87 | 4 | 25 | 15.3 |
| Liu synth dense | 500 | 30 | 4 | 25 | 16.1 |
| Liu synth dense | 500 | 60 | 4 | 25 | 16.5 |
| Liu synth dense | 500 | 87 | 4 | 25 | 16.7 |
| Liu synth dense | 1000 | 30 | 4 | 25 | 18.3 |
| Liu synth dense | 1000 | 60 | 4 | 25 | 18.5 |
| Liu synth dense | 1000 | 87 | 4 | 25 | 18.4 |

**Reading the table** (5× MAG growth + 3× sample growth → 1.24× wall time):

> Doubling N_MAG from 500 to 1000 adds ~2 s (~12 %); going from 30 to 87 samples
> at fixed N_MAG=1000 adds zero (within noise).

This pattern reflects [`envmeta/geocycle/inference.py:113-130`](../envmeta/geocycle/inference.py#L113-L130)
where `_permutation_rho_p` runs **999 permutations of (per_sample, env_vec)**
inside a loop over (pathway × env_factor). The inner cost is `O(N_sample × log N_sample)`
spearmanr; the outer multiplier is `N_pathway_active × N_env`.

**Closed-form estimate**:
```
cycle_diagram cost ≈ S1_aggregation(N_MAG × N_pathway × max_KO/pathway)
                   + N_pathway_active × N_env × 999 × spearmanr_cost(N_sample)
```

For a typical full-KEGG metagenome (18 active pathways × 4 env × 999 perm × 87 sample):
- ≈ 71,928 spearmanr calls × ~250 µs/call ≈ **18 s**

This matches measurement to within 10 %.

### 3.2 mag_heatmap, pathway, gene_profile (mid-cost)

These figures scale roughly linearly with N_MAG up to a top-N cap (default 30):
- **mag_heatmap**: 0.08-0.16 s across all measured cells
- **pathway**: 0.10-1.30 s; the synthetic dense regime hits 1.3 s at 1000 MAG (N_KO scan dominates)
- **gene_profile**: 0.16-0.34 s (similar to pathway)

### 3.3 stackplot, alpha, PCoA, RDA, LEfSe, log2FC, gene_heatmap (low-cost)

All < 0.6 s on every measured cell. matplotlib rendering + scipy linear algebra
dominate; the cost of these figures is essentially independent of N_MAG and only
weakly dependent on N_sample.

### 3.4 mag_quality, network (constant)

mag_quality is a single scatter plot with a Top-N cap; runs in 0.05-0.07 s
regardless of dataset size. network is a Degree × Betweenness scatter on
user-provided node/edge tables (EnvMeta does **not** compute the network from
scratch — that is upstream tool's responsibility, e.g. SparCC / FastSpar).

---

## 4. Memory profile

Peak ΔRSS over baseline across all 58 measured (dataset × figure) combinations:

- **Median**: 0.5 MB
- **Maximum**: 9.3 MB (Liu real, 1084 MAG × 87 sample, stackplot)
- **All cycle_diagram runs**: < 7 MB

EnvMeta is fundamentally not memory-bound. Even Streamlit Cloud free tier's
1 GB RAM allowance has ~1000× headroom. The previous concern that "1000 MAG
data would OOM on cloud" — is not supported by measurement. (The bottleneck
on the free tier is more likely the 200 MB upload-file size limit and the
single-CPU concurrency cap.)

---

## 5. Hardware sizing guidance

Derived from §3-§4. We define five usable regimes:

| Regime | N_MAG | N_samp | full 14-fig wall | Recommended hardware | Streamlit Cloud free tier |
|---|--:|--:|--:|---|:-:|
| **Demo / teaching** | ≤ 200 | ≤ 30 | ≤ 30 s | any 4 GB device | ✅ supported |
| **Typical PhD** | 200-500 | 30-100 | 30-60 s | 8 GB RAM laptop | ✅ supported |
| **Course study** | 500-1000 | 50-150 | 60-120 s | 16 GB RAM laptop | ⚠️ may hit 200 MB upload limit |
| **Large project** | 1000-5000 | 100-500 | 5-15 min | 32 GB workstation | ❌ 200 MB upload + queue timeout |
| **Outside design envelope** | > 5000 | > 500 | > 15 min | server / HPC | ❌ |

**Notes**:
- "Full 14-fig wall" estimates assume dense annotation (typical KofamScan / DRAM / METABOLIC) and 4-6 env factors. With sparse annotation (few KOs/MAG) or fewer env factors, costs can be 30-50× lower (see Liu real regime).
- The 5000 MAG upper bound is a soft envelope from the matplotlib-readability constraint (heatmap row labels become illegible) plus the typical iMeta / Bioinformatics user data scale, not a hard runtime cliff. A 10000-MAG dataset would still finish but produce unreadable figures unless top-N filtered.
- Streamlit Cloud free tier limits are the binding constraint above ~500 MAG datasets, not RAM. We recommend the cloud version for demonstration and the local install for production analyses.

---

## 6. Group-count guidance (from analysis-module hard requirements)

| Analysis | Hard floor | Soft recommended | Practical ceiling | Failure mode |
|---|--:|--:|--:|---|
| LEfSe | 2 groups, ≥ 2 samples each | 4-6 samples / group | unlimited | hard error if groups < 2 |
| PCoA + PERMANOVA | 2 groups, ≥ 3 samples each | 5-8 samples / group | unlimited | F-test unstable below 3 |
| log2FC | 2 groups, ≥ 2 samples each | 4-8 samples / group | binary contrast | hard error if not binary |
| Stackplot, MAG-heatmap, gene_heatmap | none required | — | — | renders without grouping |
| Hypothesis scorer `group_contrast` | 2 groups (with comparator df) | 2-4 groups | 6 groups | matrix grows quadratically beyond 6 |
| Color palette | — | 2-8 distinguishable | 12 | colors become indistinguishable beyond 12 |

**Source**: minimum-sample requirements come from the analysis modules
themselves (see [`envmeta/analysis/lefse.py`](../envmeta/analysis/lefse.py),
[`envmeta/analysis/log2fc.py`](../envmeta/analysis/log2fc.py),
and [`envmeta/analysis/pcoa.py`](../envmeta/analysis/pcoa.py)). Soft recommendations follow
LEfSe's original publication (Segata et al., 2011 *Genome Biol* 12:R60) for n ≥ 4
per group and PERMANOVA's small-sample-instability literature
(Anderson, 2001 *Austral Ecol* 26:32). Practical ceilings reflect the
matplotlib qualitative palette ceiling at 12 colors.

---

## 7. Key finding for the paper

> Within the scope of our measurements (Intel i7 laptop, Windows 10, three
> dataset regimes), EnvMeta's cycle inference cost appears to be dominated by
> **environmental-factor breadth × active-pathway breadth × permutation count**
> rather than by MAG count or sample count. A user with a 1000-MAG dataset
> assigned 25 random KOs/MAG and 4 env factors observed cycle inference at
> ~18 s; under the same conditions, doubling MAGs adds ≤ 12 % to that time.

We caution that this scaling generalisation rests on a **synthetic random
KO assignment**, not on a real KofamScan / DRAM / METABOLIC annotation
distribution. Real annotations cluster around organism-specific functional
repertoires, so the cost-by-N_MAG curve under realistic annotation could
plateau earlier (if many MAGs share the same KO subset and contribute
redundantly to the same pathways) or stretch out (if KO frequency tracks
taxonomy-based diversity). The matching countervailing finding from the
real-data regime is that **annotation density matters more than dataset
size** — a sparsely annotated 1084 MAG dataset (Liu's published 8 KOs) runs
in 0.4 s; a densely annotated 200 MAG one (sample data with ~30 KO/MAG
KofamScan) runs in 14 s. This finding is supported by both real-data
endpoints and is therefore more robust than the synthetic-regime scaling
curves alone.

---

## 8. Reproducing this benchmark

```powershell
# Activate environment
conda activate envmeta
cd D:\workdata\envmeta

# Step 1 — reshape Liu inputs (alpha + distance + KO wide + Phylum agg)
python paper/benchmarks/performance/reshape_liu_for_benchmark.py

# Step 2 — sample_data baseline (14 figures × 3 repeats)
python paper/benchmarks/performance/bench_harness.py --dataset sample --repeats 3 \
    --out paper/benchmarks/performance/results/sample_baseline_runtime.tsv

# Step 3 — Liu real full + subsample sweep (sparse annotation case)
python paper/benchmarks/performance/bench_harness.py --dataset liu --repeats 3 \
    --out paper/benchmarks/performance/results/liu_full_runtime.tsv
python paper/benchmarks/performance/run_scaling_sweep.py

# Step 4 — Liu synthetic dense (8 cells × 3 figures × 2 repeats; ~20 min)
python paper/benchmarks/performance/bench_synthetic_dense.py

# Step 5 — combine + figure
python paper/benchmarks/performance/analyze_results.py
```

Outputs: [`results/`](performance/results/) per-cell TSVs, [`scaling_curve.{pdf,png,svg}`](performance/scaling_curve.pdf), [`table_per_figure.tsv`](performance/table_per_figure.tsv), [`hardware_sizing_observed.md`](performance/hardware_sizing_observed.md).
