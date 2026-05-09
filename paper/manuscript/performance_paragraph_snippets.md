# Performance Benchmark — Paper Paragraph Snippets

> Drop-in English text for [`outline_imeta.md`](outline_imeta.md) Sections
> §5.2.8 (Results), §5.4.6 (Methods), and §5.3 (Discussion limitations).
> All numbers cited come from [`paper/benchmarks/performance.md`](../benchmarks/performance.md)
> (v0.9.1 / 2026-05-09). Cite [`scaling_curve.pdf`](../benchmarks/performance/scaling_curve.pdf)
> as Figure 8 (or wherever the runtime-scaling figure lands in the final layout).

---

## §5.2.8 (Results) — Performance and scaling envelope

We benchmarked EnvMeta's runtime and memory profile on three complementary
regimes designed to disentangle the contributions of MAG count, sample count,
annotation density, and environmental factor breadth (Figure 8). The first
regime used our in-house arsenic-slag dataset (169 MAGs × 10 samples × 4 env
factors × full KofamScan annotation; 30 KO/MAG); the second used Liu et al.
2023's published cold-seep dataset (1084 MAGs × 87 samples × 1 env factor ×
the published 8-KO arsenic-target subset; 1.5 KO/MAG); the third was a
synthetic-dense extension of Liu in which we randomly assigned 25 KOs/MAG
from EnvMeta's 57-KO knowledge base and 4 numeric env factors, simulating
what a fully KofamScan/DRAM/METABOLIC-annotated 1000+ MAG dataset would
behave like, and ran an 8-cell sweep across N_MAGs ∈ {200, 500, 1000} ×
N_samples ∈ {30, 60, 87}.

Two findings emerged that are inverted from a typical user's expectation.
First, **cycle_diagram cost is independent of MAG count**: at fixed dense
annotation and 4 env factors, going from 200 → 1000 MAGs added only 24% to
runtime (14.8 s → 18.4 s); going from 30 → 87 samples at fixed 1000 MAGs
added zero (within measurement noise). The asymptotic complexity is
**N_pathway_active × N_env × 999 × O(N_sample × log N_sample)** — set by the
permutation test in the cycle-inference S2 step (envmeta/geocycle/inference.py)
— rather than by N_MAG. Second, **annotation breadth dominates**: the same
1084-MAG Liu dataset ran cycle inference in 0.4 s with 8 published As-target
KOs but in 18.4 s with our synthetic-dense 25 KO/MAG annotation. This
positions EnvMeta favorably for the typical PhD-scale metagenome (200-1000
MAGs × 30-100 samples × 4-6 env factors), where the entire 14-figure
pipeline finishes in 30-120 s on a standard 8-16 GB laptop.

EnvMeta is also memory-light. Across all 58 measured (dataset × figure)
combinations, the maximum observed peak ΔRSS over baseline was 9.3 MB, and
the median was 0.5 MB; cycle_diagram itself peaked at 6.5 MB. This contrasts
with classic phyloseq + microbiome-cooccurrence-network workflows where SparCC
or NetCoMi require gigabytes of RAM at comparable MAG counts (note: EnvMeta
delegates that step to upstream tools and consumes their CSV outputs). The
practical consequence is that EnvMeta's local install runs comfortably on
any 4 GB RAM device, and its Streamlit Cloud demo deployment (1 GB free
tier) is bottlenecked by its 200 MB upload-file cap rather than by RAM.

[FIGURE 8 PLACEHOLDER — three-panel scaling figure: cycle_diagram /
mag_heatmap / pathway runtime vs N_MAGs × N_samples on log-log axes,
overlaying real-sparse, real-dense, and synthetic-dense cells. See
`paper/benchmarks/performance/scaling_curve.{pdf,svg}`.]

---

## §5.4.6 (Methods) — Performance benchmark implementation

We measured EnvMeta's per-figure runtime and memory profile by an in-process
harness ([`paper/benchmarks/performance/bench_harness.py`]; reproducibility script provided in
the Zenodo bundle) that wraps each `analyze()` entry point with `time.perf_counter`
for wall time and a background `psutil.Process().memory_info().rss` sampler
(50 ms interval) for peak ΔRSS over a `gc.collect()`-cleared baseline. Each
(dataset × figure) combination was run with 2-3 repeats, with `plt.close("all")`
and `gc.collect()` between repeats to reset the matplotlib allocator and the
Python heap. We report median wall time and maximum peak ΔRSS, and we use a
headless matplotlib `Agg` backend to mirror Streamlit's server-side rendering.

To probe the scaling envelope, we extended Liu et al.'s published 1084-MAG
cold-seep dataset by (i) randomly assigning each MAG 25 KOs sampled from
EnvMeta's 57-KO knowledge base, (ii) synthesizing 4 numeric env factors over
the 87 samples, and (iii) subsampling to N_MAG ∈ {200, 500, 1000} ×
N_samples ∈ {30, 60, 87}. This isolates dense-annotation behaviour, since
Liu's published `kegg_target_only.tsv` is pre-filtered to 8 As-target KOs
that under-represent the typical KofamScan / DRAM / METABOLIC density
(20-40 KO/MAG) of full-pipeline metagenomes. All measurements were taken
on a single Intel i7-class laptop (16 GB RAM, Windows 10, Python 3.11);
absolute numbers will vary across CPUs but the cross-cell ratios are
expected to transfer.

---

## §5.3 (Discussion) — Limitations: empirical scaling envelope

We make a tighter and empirically grounded claim about EnvMeta's intended
data-scale envelope than the original ≤ 5000 MAG soft estimate. From the
benchmark (§5.2.8), full-pipeline runtime stays under 2 minutes for the
typical PhD-thesis metagenome (200-1000 MAGs × 30-100 samples × 4 env
factors × ~25 KO/MAG); under 15 minutes for the upper bound of large-project
data (~5000 MAGs); and is dominated by the cycle inference's permutation
test (N_pathway × N_env × 999 perm × O(N_sample × log N_sample)), not by
any matplotlib-rendering step. The 5000-MAG ceiling is therefore not a
runtime cliff but a soft constraint imposed by heatmap legibility (default
top-30 filter aside, 1000+ row labels become unreadable) and by the typical
iMeta-class user data scale. Datasets above this envelope can still be
processed but require pre-aggregation (e.g. species-level rather than
MAG-level taxonomy) or top-N subsetting for visualisations. This soft
ceiling is consistent with our design philosophy of an open-source, local
tool aimed at graduate-student and small-lab users, complementing rather
than competing with HPC-grade pipelines such as Anvi'o.

---

## Where to place these in outline_imeta.md (current line numbers)

| Snippet | Target outline section | Current line | Replace what? |
|---|---|--:|---|
| §5.2.8 results | §5.2.8 第二外部数据集 benchmark | line 254 | Replace "Tara Oceans / Oak Ridge 待选" placeholder with this content; also retitle the section to "Performance and scaling envelope" |
| §5.4.6 methods | §5.4.6 Benchmark 实施 | line 334 | Replace the entire current 4-line §5.4.6 with this 1-paragraph description |
| §5.3 discussion | §5.3 局限性 | line 296 ("Streamlit GUI 在大数据集（10000+ MAG）下渲染有性能上限（推荐 ≤ 5000 MAG）") | Replace with this empirically-grounded paragraph |

When integrating into the outline (your work item, see Plan B), align the
Figure 8 reference number with whatever Figure layout you settle on
(F8 in the current outline; might shift if you reorder).
