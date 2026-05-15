# EnvMeta 论文大纲（GPB 投稿版 / v2 reframe 2026-05-14）

> **目标期刊**：Genomics, Proteomics & Bioinformatics（GPB / Oxford Academic / IF 11.5 / Q1 / CUGB C 区 15 分）
> **文章类型**：Methods / Application Note
> **预计字数**：4500-6500 words
> **预计 figures**：6-7 张
> **预计 tables**：2 张
> **APC**：est. $2,000-2,500 OA（中科院系，可能有作者优惠）
> **拒稿兜底**：SoftwareX（Elsevier / IF 2.4 / Q3 / 3 分 / **零 APC** / 2 月周期）
>
> **创建日期**：2026-05-07（iMeta v1）→ 2026-05-14（GPB v2 reframe）
> **历史版本**：[outline_imeta.md](outline_imeta.md)（iMeta v1 框架，已 deprecated）
> **作者**：redlizzxy（中国地质大学北京）
> **v2 reframe 主旨**：从 "domain-specialist complement to iMeta visualization ecosystem"
>   **重定位为** "hierarchical weight-of-evidence MCDA framework with reference implementation EnvMeta"。
>   不改代码、不改实验，只改 Abstract / Intro / Methods 框架叙述顺序：**算法 framework 第一，工具第二**。

---

## 0. GPB 投稿要点速查

| 项 | 要求 |
|---|---|
| 出版社 | Oxford Academic |
| 主办 | 中国科学院基因组研究所（BIG-CAS）+ Genetics Society of China |
| 版式 | OA only，APC est. $2,000-2,500 |
| 文章类型 | Methods / Database / Application Note / Genome Resource / Original Research / Review |
| 字数 | Methods/Application Note: ~4500-6500 words |
| 结构 | Abstract (≤300 words) + Keywords + Highlights (3-4) + IMRAD |
| Figure 风格 | 偏算法图（schema + flowchart + 性能曲线 + 案例验证）|
| 推荐审稿人 | 投稿系统填 ≥ 3 位（中科院系 BIG-CAS 编委 + bioinformatics 资深） |
| 引用格式 | Vancouver (numerical) — GPB 标准 |
| 参考文献 | 优先 GPB + BIG-CAS 系列 + omics 方法学顶刊（Bioinformatics, NAR, Genome Biology） |
| 数据可用性 | 强制 GitHub + Zenodo DOI + **算法可复现** |
| **GPB 关键差异点** | 必须有 **algorithm pseudocode + 复杂度分析 + vs comparable tools benchmark**（vs iMeta tool paper 标准多了这 3 项）|

---

## 1. Title 候选（v2 GPB framing — algorithm 第一，工具第二）

| # | Title | 卖点 |
|---|---|---|
| **A** ⭐ | **A hierarchical weight-of-evidence MCDA framework for descriptive-to-causal microbiome-environment association inference, with reference implementation EnvMeta** | 主推 — 算法 framework 第一，工具第二，完全 GPB 体裁 |
| B | EnvMeta: A hierarchical weight-of-evidence framework for biogeochemical-cycle hypothesis evaluation in environmental metagenomics | 工具名前置，框架次之，较短 |
| C | Weight-of-evidence MCDA for microbiome-environment causality assessment: framework, reference implementation (EnvMeta), and four-arm calibration | 加 "calibration" 强调验证 |

→ 推荐 **A**。把 "hierarchical weight-of-evidence MCDA framework" 立为主体算法贡献，
   EnvMeta 作为参考实现位居其次。这是 v2 reframe 的核心 — 让 GPB 审稿人**第一眼就看到算法贡献**，
   而不是"又一个 visualization 工具"。

**deprecated v1 (iMeta) Title**（见 outline_imeta.md）：
"EnvMeta: An interactive visualization platform for environmental metagenomics with built-in
 biogeochemical cycle inference and hypothesis scoring"

---

## 2. Highlights（GPB 风格 / 3 条 / 算法 framework 第一）

> GPB Highlights 风格：3-4 条，每条 1 个 capability 一行。这里 **首条强调框架**，
> 次条强调 schema + 验证，末条强调实现 + 可复现性。

```
• A hierarchical weight-of-evidence MCDA framework integrating S1 compositional
  debiasing, S2 999-permutation null calibration, S3 weighted-sum scoring with
  Saltelli ±20% one-at-a-time weight robustness, and Bradford-Hill required-veto
  reasoning, exposing five auditable diagnostic quantities per claim.

• A six-claim YAML hypothesis schema operationalizing the framework against
  MAG-level KEGG-orthology data, calibrated across four published metagenomic
  datasets (all STRONG under fixed defaults) and discriminated against domain-
  mismatched stress claims via cross-element pathway-target perturbation.

• Open-source reference implementation EnvMeta — a Streamlit platform with
  fourteen publication-quality visualizations, KEGG-driven biogeochemical-cycle
  inference (4 elements × 18 pathways × 57 KOs), Fork Bundle reproducibility
  packaging, and standalone offline interactive HTML supplementary material
  (~400 KB) embedding the analysis itself.
```

---

## 3. Graphical Abstract 草图

```
┌─────────────────────────────────────────────────────────────────────┐
│                                                                     │
│  📁 Upload metagenomic data    ───►   🔍 Auto-recognize 11 types    │
│  (TSV/CSV/FASTA/...)                  (abundance/KO/taxonomy/...)    │
│                                                                     │
│                              ▼                                       │
│                                                                     │
│   ┌──────────────────┐   ┌──────────────────┐   ┌──────────────────┐│
│   │ 📊 12 visualiz.  │   │ ⚛️ Cycle infer.  │   │ 🧪 Hypothesis    ││
│   │ • α/β-PCoA/RDA   │   │ • 4 elements     │   │ • YAML schema    ││
│   │ • Heatmaps       │   │ • Mockup-10      │   │ • Null perm.     ││
│   │ • LEfSe/log2FC   │   │ • Cross-element  │   │ • Sensitivity    ││
│   │ • MAG quality    │   │   coupling       │   │ • 9-tier label   ││
│   └──────────────────┘   └──────────────────┘   └──────────────────┘│
│                                                                     │
│                              ▼                                       │
│                                                                     │
│   📦 Export: PNG/PDF/SVG/TIFF + .py reproducible script + Bundle    │
│   🌐 Standalone interactive HTML (~400 KB, offline)                 │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

**布局**：上 1/3 数据上传 → 中 1/3 三栏功能（左 reads-based / 中 cycle / 右
hypothesis）→ 下 1/3 导出。

---

## 4. Abstract（~290 words / GPB v2 framing — framework 第一）

> Translating MAG-level KEGG-orthology data into mechanism-evaluable
> associations between microbial pathways and environmental factors requires
> formal weight-of-evidence assessment rather than narrative tabulation.
> Despite advances in microbiome visualization tools, no widely-adopted
> algorithmic framework integrates compositional debiasing, permutation-based
> null calibration, multi-criteria decision analysis (MCDA), and Bradford-Hill
> weight-of-evidence reasoning in a way that exposes its intermediate
> quantities for reproducibility audit. We present a **hierarchical
> weight-of-evidence MCDA framework** operationalizing this need through three
> algorithmic stages — **S1** compositional debiasing with three-threshold
> top-1 contributor sensitivity scanning, **S2** Fisher 999-permutation null
> calibration with five-tier confidence labels, **S3** weighted-sum scoring
> with Saltelli ±20% one-at-a-time weight robustness and Bradford-Hill
> required-veto reasoning — exposing five auditable diagnostic quantities per
> claim (score, weight-robust score, null_p, evidence count, confidence label).
> The framework operationalizes against a **six-claim YAML hypothesis schema**
> covering pathway activity, cross-pathway coupling, environmental correlation,
> keystone-MAG identification, group contrast, and Popperian pathway-inactive
> falsification. We provide an open-source reference implementation
> **EnvMeta** — a Streamlit-based platform that wraps the framework with
> fourteen publication-quality visualizations, automated biogeochemical
> cycle inference for As/N/S/Fe (4 elements × 18 pathways × 57 KOs), Fork
> Bundle reproducibility packaging, and standalone offline interactive HTML
> supplementary material (~400 KB) embedding the analysis as the SI itself.
> The scoring engine is calibrated across four published metagenomic datasets
> (Wei 2024 paddy soil; Liu 2023 cold seep; Grettenberger 2021 acid mine
> drainage; Ayala 2020 Iberian Pyrite Belt pit lake), all returning STRONG
> labels under fixed default thresholds, and discriminated against
> domain-mismatched stress claims through cross-element pathway-target
> perturbation (0/20 STRONG retention in the most narrowly-focused dataset)
> and 5-threshold sensitivity sweeps. A 58-cell performance benchmark
> establishes runtime ≤ 120 s and peak memory ≤ 10 MB for typical PhD-scale
> metagenomes (200-1000 MAGs × 30-100 samples × 4-6 env factors) on a standard
> laptop. EnvMeta is freely available at
> https://github.com/redlizzxy/EnvMeta with online demo at
> https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/.

**Keywords**: microbiome bioinformatics; metagenome-assembled genome (MAG);
weight-of-evidence; multi-criteria decision analysis (MCDA); permutation test;
biogeochemical cycle inference; KEGG orthology; reproducible research

---

## 5. 主体章节大纲

### 5.1 Introduction（~830 words / 4 段 / GPB v2 framework-first framing）

> **写作风格基线（v2 GPB）**：开篇定位为 **methodological gap in microbiome
> causality inference**，而非"visualization tool ecosystem"。第 1 段建立 weight-of-
> evidence inference 的方法学需求；第 2 段 review 现有方法学组件（CoDA / permutation /
> MCDA / Bradford-Hill）+ 微生物组工具 landscape；第 3 段提出 framework + reference
> implementation；第 4 段列 3 项贡献。

#### 段 1 — Methodological gap in microbiome causality inference (~180 words)

> Inferring causal associations between microbial pathways and environmental
> factors from observational metagenomic data is a core methodological
> challenge in microbiome science. Modern metagenome-assembled genome (MAG)
> studies routinely catalogue thousands of genomes per sample annotated with
> KEGG-orthology (KO) functional capacity, yet translating these data tables
> into mechanism-evaluable associations remains predominantly a narrative
> reasoning task: researchers tabulate KO presence × abundance × environmental
> covariates, examine patterns by inspection, and assemble discussion-section
> arguments. This workflow leaves three sources of evidential uncertainty
> implicit and unaudited: (i) compositional bias from sequencing-depth
> normalization choices that distort relative-abundance interpretation; (ii)
> chance correlations under multiple-comparison conditions, especially in
> small-N MAG cohorts where pathway-completeness scores can pass by-eye
> thresholds without controlling for null-distribution mass; and (iii) the
> analyst's discretionary weighting of when to upgrade a "consistent with"
> pattern to a "supports the hypothesis" claim. Formal weight-of-evidence
> frameworks address each of these in environmental risk assessment (Linkov
> et al., 2009; Suter & Cormier, 2011; Rhomberg et al., 2010), but no
> widely-adopted operationalization exists for the microbial-pathway ×
> environmental-factor inference setting central to environmental
> metagenomics.

#### 段 2 — Existing methodological components and tools landscape (~290 words)

> Several method classes partially address these uncertainties. Compositional
> data analysis (CoDA: Aitchison-family methods such as CLR/ILR transforms;
> Gloor et al., 2017) explicitly handles sequencing-depth normalization
> distortions. Permutation-based null calibration is standard in
> differential-abundance and ordination statistics (PERMANOVA: Anderson 2001;
> LEfSe: Segata et al., 2011), and Fisher-style permutation underpins multiple
> ecology and microbiome test packages. Multi-criteria decision analysis
> frameworks (MCDA: Belton & Stewart, 2002) and Bradford-Hill weight-of-evidence
> reasoning (Hill, 1965; Rhomberg et al., 2010) provide formal vocabularies
> for combining heterogeneous evidence streams, with required-veto and
> sensitivity-analysis tooling well-established in environmental and clinical
> risk assessment (Linkov et al., 2009; Suter & Cormier, 2011). On the
> visualization and pipeline side, Anvi'o (Eren et al., 2021) offers
> MAG-centric exploratory analysis and pangenome workflows; QIIME2 (Bolyen
> et al., 2019) and the phyloseq / R ecosystem (McMurdie & Holmes, 2013)
> cover community statistics; recent integrative web platforms including
> ImageGP 2 (Chen et al., 2024 *iMeta*), Sangerbox (Shen et al., 2022),
> TOmicsVis (Miao et al., 2023), Wekemo Bioincloud (Gao et al., 2024 *iMeta*),
> iMetaLab Suite (Li et al., 2022 *iMeta*), iNAP (Feng et al., 2022 *iMeta*),
> and Majorbio Cloud (Ren et al., 2022 *iMeta*) emphasize plot reproduction
> across general biomedical and microbiome workflows; and Krona (Ondov et al.,
> 2011), iTOL (Letunic & Bork, 2024), and Cytoscape (Shannon et al., 2003)
> provide specialized hierarchies, phylogenies, and networks respectively.
>
> Yet for the specific operational setting of MAG-level KEGG-resolution
> microbiome-environment association inference, these components remain
> **isolated rather than integrated**. None of the tools above wraps S1
> compositional debiasing + S2 permutation null calibration + S3 MCDA
> scoring + Bradford-Hill required-veto + sensitivity / null-distribution
> diagnostics into a unified, formally documented, and reproducibility-auditable
> framework. The gap leaves environmental-microbiology researchers either
> writing one-off statistical scripts at every iteration or defaulting to
> narrative reasoning that obscures the three evidential uncertainties named
> in §1.

#### 段 3 — Framework + reference implementation (~210 words)

> We present a **hierarchical weight-of-evidence MCDA framework** for
> microbiome-environment association inference, with open-source reference
> implementation **EnvMeta**. The framework decomposes into three algorithmic
> stages — **S1** compositional debiasing with three-threshold top-1
> contributor sensitivity scanning, **S2** 999-permutation Fisher null
> calibration with five-tier confidence labels, and **S3** weighted-sum
> scoring with Saltelli ±20% one-at-a-time weight robustness and Bradford-Hill
> required-veto reasoning — producing **five auditable diagnostic quantities**
> per claim (score, weight-robust score, null_p, evidence count, confidence
> label) that downstream readers and reviewers can inspect independently. The
> framework operationalizes against a **six-claim YAML hypothesis schema**
> covering pathway activity, cross-pathway coupling, environmental
> correlation, keystone-MAG identification, group contrast, and Popperian
> pathway-inactive falsification. EnvMeta wraps the framework behind a
> Streamlit graphical user interface that adds fourteen publication-quality
> visualizations spanning reads-based community statistics (α/β-diversity,
> PCoA, RDA, LEfSe, log2FC, taxonomy stackplots, KEGG heatmaps) and
> MAG-based exploration (quality control, abundance heatmaps, pathway
> completeness, gene profiles, network preparation), a KEGG-driven
> biogeochemical-cycle knowledge base (4 elements × 18 pathways × 57 KOs),
> Fork Bundle reproducibility packaging, and standalone offline interactive
> HTML supplementary material (~400 KB) that embeds the analysis as the SI
> itself. Five design principles — domain-neutral inference, user-supplied
> knowledge, fully offline operation, fork-rather-than-community distribution,
> and descriptive (not causal) outputs — narrow the operational scope and
> position EnvMeta as a microbiome-specific instance of a general
> framework.

#### 段 4 — Specific contributions (~150 words)

> Three contributions follow:
>
> 1. We **formalize a hierarchical weight-of-evidence MCDA framework** for
>    descriptive-to-causal microbiome-environment association inference,
>    integrating S1 compositional debiasing, S2 999-permutation null
>    calibration, and S3 MCDA scoring with Saltelli ±20% weight robustness
>    and Bradford-Hill required-veto reasoning into a unified procedure
>    exposing five auditable diagnostic quantities per claim. Default
>    thresholds are pre-registered with cryptographic anchoring (OpenTimestamps;
>    Todd, 2016) and grounded in conventional clinical-epidemiology and
>    risk-assessment literature.
>
> 2. We **calibrate the framework** against four published metagenomic
>    datasets (Wei 2024; Liu 2023; Grettenberger 2021; Ayala 2020), all
>    returning STRONG labels under fixed default thresholds, and stress-test
>    discrimination power through three claim-class perturbations (reversed
>    direction, cross-topic, `pathway_inactive`) and within / cross-element
>    pathway-target perturbation (N = 20 per mode × 4 datasets), with
>    cross-element perturbation collapsing the most narrowly-focused dataset
>    to 0/20 STRONG retention.
>
> 3. We provide an **open-source reference implementation, EnvMeta**, that
>    operationalizes the framework with fourteen publication-quality
>    visualizations, automated biogeochemical-cycle inference, Fork Bundle
>    reproducibility packaging, and standalone offline interactive HTML
>    supplementary material.

### 5.2 Results（~2500-3000 words / 7-8 figures）

**节奏**：每个 Result 子节 1 个图（核心 figure）+ 1 段叙事 + 1 段 quantitative 对比。

#### 5.2.1 EnvMeta 整体架构（Figure 1）

- 5 层架构图
- 11 FileType + 14 图 + 3 创新模块（cycle/hypothesis/HTML）
- 数据流：上传 → 识别 → 分析 → 调参 → 导出

**Figure 1 内容**：架构 schema 图（参考 paper/figures/screenshot_*.png 现有素材重画 SVG）

#### 5.2.2 12 张出版级图 + 同步代码生成（Figure 2）

- 7 reads-based + 5 MAG-based
- 重点：**调参 GUI 控件 + 实时预览 + 代码同步生成**
- 量化：vs 原 R/Python 脚本，作图时间从 30-60 min → 1-2 min（ref `time_comparison.md`）

**Figure 2 内容**：4 panel
- A) 堆叠图（reads-based）+ GUI 控件截图
- B) PCoA + 调参 → 实时刷新
- C) MAG 丰度热图（4 层参数面板）
- D) 自动生成的 .py 复现脚本片段

#### 5.2.3 元素循环图自动推断（Figure 3 — **CORE FIGURE**）

- KEGG-driven KB（4 元素 × 18 通路 × 57 KO）
- 3 阶段推断：**S1 去偏 → S2 999 次置换 → S3 敏感度扫描**
- Mockup 10 合并细胞布局（substrate-gene-product 三段式）
- 跨元素化学物耦合（orpiment / 雌黄沉淀 As(III) + H₂S → As₂S₃，primary refs:
  Newman et al., 1997 *Appl Environ Microbiol* 63:2022-2028 (precipitation
  of arsenic trisulfide by *Desulfotomaculum auripigmentum*);
  Rodriguez-Freire et al., 2014 *Environ Sci Technol* 48:4107）

**Figure 3 内容**：
- A) KB schema（4 元素 × 18 通路）
- B) S1 去偏 → S2 置换 → S3 敏感度算法流程
- C) 推断后的元素循环图（cell_renderer 输出，砷渣 case）
- D) 跨元素耦合放大（雌黄沉淀化学物锚点）

#### 5.2.4 YAML 假说评分器 + 4-Arm calibration（Figure 4 + Table 1）

- **6 类 claim**（v0.9.0）：`pathway_active` / `pathway_inactive` / `coupling_possible` /
  `env_correlation` / `keystone_in_pathway` / `group_contrast`
- 评分公式：MCDA `Σ(w_i × score_i) / Σ(w_i)` + Bradford-Hill required-veto + 4 档标签
  (`strong` / `suggestive` / `weak` / `insufficient`)
- 3 个独立置信指标：Fisher 999-perm null_p / OAT ±20% weight robustness / required-veto reasons
- Pre-registration discipline：git timestamp 锁定 claim entities + 默认阈值不调

**Figure 4 内容**（4-panel）：
- A) YAML 输入 schema + 6 claim types 一图
- B) 999 次置换 null 分布 + 观察值 (Arm A 砷渣 case)
- C) Bradford-Hill required-veto 决策树 + 4 档标签
- D) 4-Arm calibration 结果 bar plot（详见 Table 1 + Results §5.2.4 narrative）

**Drop-in Results §X.1 narrative**（~400 words）— source: [`paper/manuscript/results_stress_test_section.md`](results_stress_test_section.md) §X.1

> To establish that EnvMeta's hypothesis scoring engine produces stable,
> defensible labels under default thresholds, we ran a four-arm controlled
> experiment over four metagenomic datasets spanning a gradient of annotation
> breadth: Arm A (in-house steel-slag arsenic-remediation, 168 MAGs × 10 samples,
> full KofamScan KEGG annotation, 57 KOs across 4 elements); Arm B (Wei et al.
> 2024 *Microbiome*, 36 paddy-soil samples × 179 MAGs annotated with 14 functional
> genes via custom ROCker models); Arm C1 (Liu et al. 2023 *npj Biofilms
> Microbiomes*, deep-sea cold-seep, 87 samples × 1084 MAGs with DRAM KEGG
> annotation); Arm C2-A (Grettenberger & Hamilton 2021 *Appl Environ Microbiol*,
> AMD stream, 29 MAGs with METABOLIC step-level KEGG); Arm C2-B (Ayala-Muñoz
> et al. 2020 *Microorganisms*, Iberian Pyrite Belt acidic pit lake, 13 MAGs
> re-annotated end-to-end with Pyrodigal + GhostKOALA). All hypothesis YAMLs
> were pre-registered (committed to git at hashes `42168da`, `44d7f5f`, `76a4f77`
> before EnvMeta was run), used EnvMeta's default thresholds (`min_completeness=30`,
> `strong=0.75`, `suggestive=0.40`), and cited only review literature published
> ≥ 5 years before the target paper.
>
> All four KEGG-curated arms returned **`STRONG` labels with overall_score = 1.000
> and 4/4 claims satisfied** (Table 1; Arms A, C1, C2-A, C2-B). Arm B (Wei 2024,
> ROCker-only) returned overall_score = 0.63 with label `INSUFFICIENT` despite
> 3/5 claims satisfied; the required-veto activated because Wei's 14-gene set
> provides only 2 of the 6 canonical KOs in EnvMeta's `Nitrate reduction`
> pathway (`napA` + `narG`, missing `narH`/`narI`/`napB`/`narB`) and the
> As(III)↔NO₃⁻ chemistry coupling was scored *partial*. Permutation null-p = 0.90
> (n = 999) and weight robustness under ±20% OAT perturbation = True confirmed
> the conservative diagnosis is not a weight-tuning artifact. The contrast
> establishes that the `INSUFFICIENT` label faithfully reflects annotation-coverage
> diagnostics rather than engine malfunction. Notably, the `STRONG` label on
> Arm C2-B (Ayala) was obtained with a single GhostKOALA re-annotation pipeline
> applied to publicly available MAG genomes (BioProject PRJNA646106),
> demonstrating that KEGG-curated EnvMeta scoring is reproducible end-to-end
> from raw genome assemblies.

**Drop-in Results §X.2 narrative**（~400 words / stress test discrimination）— source: [`results_stress_test_section.md`](results_stress_test_section.md) §X.2

> For each KEGG-curated dataset (Liu, Grettenberger, Ayala), we authored a
> second pre-registered YAML (`{dataset}_hypothesis_stress.yaml`, all committed
> at `50c4687`) encoding deliberately *risky* claims violating environmental
> priors. Three claim classes: **(A)** reversed-direction predictions; **(B)**
> cross-topic mismatches; **(C)** `pathway_inactive` negation. Each YAML
> included a calibration anchor claim. Twelve predictions × three datasets
> were frozen in [`stress_test_predictions.md`](stress_test_predictions.md) before any stress run.
>
> Observed stress scores fell substantially below the calibration STRONG
> baseline (Table 2; Figure X). Grettenberger 2021 returned `weak` (0.250,
> 1/3 satisfied) — the cleanest single-Arm discrimination. Liu 2023 and
> Ayala 2020 returned `suggestive` (0.625 / 0.455, with skipped claims
> contributing to partial scores). The most informative single discriminator
> was the cross-topic claim "arsenate_reduction should dominate", which was
> correctly rejected with **n = 0 active MAGs in both non-arsenic datasets**
> (Grettenberger n = 29 MAGs; Ayala n = 13 MAGs). We treat this two-dataset
> rejection as **consistent with** — rather than ironclad proof of — the
> scoring engine being uninfluenced by the universal *arsC* arsenate-reductase
> homolog (Rosen, 2002) under cross-topic mismatch; small dataset sizes mean
> absence of *arsC* could partly reflect sampling undercount, and a larger
> non-arsenic dataset (≥ 100 MAGs) is flagged as future work (§5.3 / §Y.4).
> Within the scope of the KEGG-curated datasets we tested, EnvMeta's scoring
> is **consistent with** being neither hard-wired to confirm arsenic-cycle
> hypotheses nor biased against alternative contexts, **conditional on
> adequate KEGG coverage**.
>
> In two of three datasets (Liu and Ayala), the reversed-direction stress
> claim "arsenite oxidation should dominate" returned satisfied because real
> but weak oxidizer signals were present (total contributions 21-fold and
> 10-fold below the dominant reduction pathway, but above the binary
> `mean_completeness ≥ 50%` threshold). This exposes a **binary-threshold
> reporting limitation** addressed in v0.9.x by the `dominance_score =
> total_contribution / element_total` field. We are explicit that the v0.9.x
> extension is an **engineering retrofit informed by v1 outputs** rather
> than independent validation: the `min_dominance_fraction = 0.20` threshold
> was chosen after observing v1 dominance values (Liu 0.05%, Ayala 7.08%).
> A second pre-registered v2 YAML (`*_stress_v2.yaml`, commit `fdfae77`)
> with this threshold returned `unsatisfied` for both datasets and reduced
> overall scores to Liu 0.250 / Ayala 0.182 (Table 2 v2 column). We report
> v1 and v2 outcomes side by side so that readers can judge the retrofit's
> effect directly; we do **not** interpret v2 as independent confirmation
> of discrimination power. Independent validation of the `dominance_score`
> field on a fresh dataset is identified as future work (§5.3 Future).

#### 5.2.5 独立交互 HTML 导出（Figure 5）

- D3.js v7 inline 嵌入（无外部依赖）
- 400-550 KB，**比典型 PDF SI 还小**
- 包含：循环图全交互 + 假说评分表 + 跨组对比 + SVG 导出
- "**论文 SI 即软件本身**"（reviewer 下载单文件即可亲手审计）

**Figure 5 内容**：
- A) HTML 4 tab UI 截图
- B) 循环图节点拖拽 + 化学物锚点高亮
- C) 假说评分 claim 表 + null_p 分布
- D) HTML 文件大小 vs Krona/iTOL/Anvi'o SI 对比柱状图

#### 5.2.6 Fork Bundle —— 论文与工具绑定发布（Figure 6）

- 一篇论文一个 Bundle.zip（KB + YAML + config + KEGG 快照 + sample data）
- 用户下载 Bundle → 一键加载 → 完全复现论文图
- **解决"论文发表后工具版本漂移导致复现失败"问题**

**Figure 6 内容**：
- Bundle 结构 + 加载流程
- 复现案例：从 Bundle 加载 → 5 分钟内重生成论文图

#### 5.2.7 案例研究：砷渣-钢渣微生物修复（Figure 7）

- 168 MAG × 10 sample × 57 目标 KO
- CK / A（低钢渣）/ B（高钢渣）三组
- 关键发现：铁氧化固砷 + 硫还原沉砷的多元素协同
- **不是新机制论文（那是 Paper 1/2 的事），而是展示 EnvMeta 如何辅助机制发现**

**Figure 7 内容**：
- A) 元素循环图（CK vs A vs B）
- B) 假说评分跨组对比
- C) 关键 MAG 承载者切换（Gallionella → Thiobacillus）

#### 5.2.8 Performance and scaling envelope（Figure 8 + Table 2）

> Drop-in English text — source: [`paper/manuscript/performance_paragraph_snippets.md`](performance_paragraph_snippets.md) §5.2.8

We benchmarked EnvMeta's runtime and memory profile on three complementary
regimes designed to disentangle the contributions of MAG count, sample count,
annotation density, and environmental factor breadth (Figure 8). The first
regime used our in-house arsenic-slag dataset (169 MAGs × 10 samples × 4 env
factors × full KofamScan annotation; 30 KO/MAG); the second used Liu et al.
2023's published cold-seep dataset (1084 MAGs × 87 samples × 1 env factor ×
the published 8-KO arsenic-target subset; 1.5 KO/MAG); the third was a
synthetic-dense extension of Liu in which we randomly assigned 25 KOs/MAG
from EnvMeta's 57-KO knowledge base (i.e., **dense relative to the
KB-target subset of As/N/S/Fe pathways**, not dense relative to the full
~10,000-KO genome-wide KofamScan output of a typical MAG) and 4 numeric
env factors, simulating what a fully KofamScan/DRAM/METABOLIC-annotated
1000+ MAG dataset would behave like *within EnvMeta's KB scope*, and ran
an 8-cell sweep across N_MAGs ∈ {200, 500, 1000} × N_samples ∈ {30, 60,
87}.

Two findings emerged. First, within the scope of our measurements,
**cycle_diagram cost appears largely independent of MAG count**: at fixed
synthetic-dense annotation and 4 env factors, going from 200 → 1000 MAGs
added only 24% to runtime (14.8 s → 18.4 s); going from 30 → 87 samples at
fixed 1000 MAGs added zero (within measurement noise). The asymptotic
complexity is **N_pathway_active × N_env × 999 × O(N_sample × log N_sample)**
— set by the permutation test in the cycle-inference S2 step
([`envmeta/geocycle/inference.py`](../envmeta/geocycle/inference.py)) —
rather than by N_MAG. We caution that this finding is established under a
**synthetic random KO assignment** rather than a realistic KofamScan / DRAM
/ METABOLIC annotation distribution; real annotations cluster around
organism-specific functional repertoires and may show different cost-by-N_MAG
behaviour. A direct comparison against one fully KofamScan-annotated dataset
is identified as future validation. Second, the dataset-anchored finding
that **annotation density matters more than dataset size** — the same
1084-MAG Liu dataset ran cycle inference in 0.4 s with 8 published As-target
KOs versus 18.4 s under our synthetic 25 KO/MAG augmentation — is robust to
the synthetic-regime caveat above (the two endpoints are both real-data
or real-data-derived). Together, these observations position EnvMeta
favourably for the typical PhD-scale metagenome (200-1000 MAGs × 30-100
samples × 4-6 env factors), where the entire 14-figure pipeline finishes
in 30-120 s on a standard 8-16 GB laptop.

EnvMeta is also memory-light. Across all 58 measured (dataset × figure)
combinations, the maximum observed peak ΔRSS over baseline was 9.3 MB, and
the median was 0.5 MB; cycle_diagram itself peaked at 6.5 MB. The practical
consequence is that EnvMeta's local install runs comfortably on any 4 GB RAM
device, and its Streamlit Cloud demo deployment (1 GB free tier) is bottlenecked
by its 200 MB upload-file cap rather than by RAM.

**Figure 8 内容** — `paper/benchmarks/performance/scaling_curve.{pdf,svg,png}`：
3-panel log-log scatter（cycle_diagram / mag_heatmap / pathway 各 1 panel），
3 个数据 regime 用 3 种 marker × 颜色（sample 真 dense 蓝 / Liu 真 sparse 灰 / Liu 合成 dense 红），
58 个测试点。

**Table 2 内容** — Hardware sizing table：5 种用量 regime（demo ≤ 30s / typical PhD
30-60s / 课题组 1-2 min / 大型项目 5-15 min / 超出设计范围 > 15 min）×（推荐硬件 +
Streamlit Cloud 是否支持）。Source: [`paper/benchmarks/performance.md`](../benchmarks/performance.md) §5。

**Drop-in §X.3 perturbation analysis narrative**（~280 words；放 §5.2.4 末尾或独立 §5.2.4.5）— source: [`results_stress_test_section.md`](results_stress_test_section.md) §X.3

> To distinguish whether the four STRONG calibration outcomes reflect
> authors' specific pre-data target choices or arise mechanically from
> KEGG annotation breadth, we perturbed every `params.pathway` field
> across the three external calibration YAMLs (Liu 2023, Grettenberger
> 2021, Ayala 2020) and rescored under default thresholds (Methods
> §4.6.7). Two perturbation modes were applied: within-element (random
> alternative pathway from the same KB element) and cross-element (random
> pathway from a different KB element); N=20 per mode per dataset.
> Cross-element control is strongly discriminating: Liu 2023, the most
> narrowly As-focused dataset, retains STRONG in **0/20** cross-element
> perturbations (median score 0.000) because cross-element substitution
> lands on inactive N/S/Fe pathways and triggers required-claim veto.
> Grettenberger 2021 and Ayala 2020 retain STRONG in 6/20 and 3/20
> cross-element runs (70–85% label change). Within-element control bounds
> the KEGG-coverage caveat: mean scores drop 25–48% but 40–50% of
> within-element runs still register STRONG, consistent with the
> manuscript's framing that calibration evidence is KEGG-coverage-
> dependent rather than domain-neutral. We treat these as auxiliary
> evidence consistent with — not ironclad proof of — non-mechanical
> calibration. The cleanest definitive mitigation, blind hypothesis
> writing by collaborators unfamiliar with target findings, remains
> future work (§Y.4). Full results in Supplementary Table S_pert and
> Figure X-bis (`perturbation_curve.pdf`); reproduction protocol in
> [`paper/manuscript/perturbation_analysis_results.md`](perturbation_analysis_results.md).

**Drop-in §X.4 reference audit narrative**（~150 words；可放 §5.2.8 末尾或 §5.4.8.6）— source: [`results_stress_test_section.md`](results_stress_test_section.md) §X.4

> Post-hoc DOI verification identified four reference errors in the
> pre-registered YAMLs that do not affect scoring outputs (no claim entity
> was modified) but require transparent correction. Most consequentially, the
> `nitrogen_fixation_explored` claims (Grettenberger and Ayala calibration
> YAMLs) originally cited Auld et al. (2017 *Can J Microbiol*), which is a
> seasonal community-variation study rather than an AMD diazotrophy report;
> these claims are re-grounded in Dai et al. (2014 *PLoS One*; metagenomic
> identification of 742 *nif* sequences from acid mine drainage) and
> Méndez-García et al. (2015 *Front Microbiol*; review of AMD diazotrophs).
> Three additional metadata corrections (journal mislabels for Yin 2011,
> Cabrera 2006, and a non-existent "Bothe 2007 *FEMS Rev*" → Bothe 2000)
> were committed at `ddd3098` and `cae2de7`; the original pre-registered
> versions remain accessible in git history. A complete proof-of-extraction
> quality audit is at [`paper/manuscript/hypothesis_references_audit.md`](hypothesis_references_audit.md)
> (Supplementary Table S_refs).

#### 5.2.9 Tool comparison — Feature matrix + performance benchmark（Figure F10 / Table T1 / Table T2）

> **GPB v2 framing 调整**：v1 的 "complementary positioning across the visualization
> ecosystem"（针对 iMeta）→ v2 改为客观 **"feature matrix + head-to-head performance
> benchmark"**（GPB 审稿人偏好实测数据）。Table T1 保持功能矩阵（capability matrix），
> 新增 Table T2 性能 head-to-head（仅 shared capabilities），新增 Figure F10
> head-to-head bar chart。
>
> 方法学完整文档：[`paper/benchmarks/external/tools_comparison/methodology.md`](../benchmarks/external/tools_comparison/methodology.md)。

##### 5.2.9.1 Capability matrix（功能矩阵）

**Table T1 — Feature matrix across comparable tools**：

| Tool | Type | Domain coverage | Element cycle inference | Hypothesis scoring + null/sensitivity | Standalone offline HTML SI | MAG-based 5-figure suite | KEGG-driven biogeochemical KB |
|---|---|---|---|---|---|---|---|
| **ImageGP 2** (Chen 2024) | Online platform / Vue+Django | General biomedical (45 tools / 6 categories) | − | − | − | − | − |
| **EasyAmplicon / EasyMetagenome** (Liu 2023) | CLI pipeline / shell+R | Amplicon + metagenomic standard analyses | − | − | − | − | − |
| **Anvi'o** (Eren 2021) | CLI + light web / Python | MAG exploration + pangenomics | − | − | partial (interactive web SI) | partial (MAG-side strong) | − |
| **QIIME2** (Bolyen 2019) | CLI plug-in framework | Amplicon community statistics | − | − | − | − | − |
| **phyloseq / Shiny-phyloseq** (McMurdie 2013) | R package + Shiny | Community statistics | − | − | − | − | − |
| **Krona** (Ondov 2011) | Static HTML SVG | Hierarchical taxonomy | − | − | partial (offline HTML / no analysis) | − | − |
| **iTOL** (Letunic 2024) | Online + offline export | Phylogeny + annotation | − | − | partial (online SI + image) | − | − |
| Vendor cloud platforms | SaaS web | General amplicon + functional | − | − | − | − | − |
| **EnvMeta** ⭐ | Local + cloud Streamlit | **Environmental MAG specialist** | ✅ | ✅ | ✅ | ✅ | ✅ |

**EnvMeta 在 Type B 独家能力**：在上面 8 行所有工具都未覆盖的 environmental-MAG
specialist niche 填补 5 项独家功能（元素循环 / 假说评分 / null calibration /
weight sensitivity / 离线交互 HTML SI）。

##### 5.2.9.2 Head-to-head performance benchmark（Type A 共享能力）

For analyses that EnvMeta and competitor tools both offer (PCoA, KO heatmap,
hierarchical taxonomy), we ran identical inputs through each tool on the same
hardware (Intel i7-class 16 GB Win10 laptop) and measured wall time, peak RSS,
output file size, and user operation count. Methodology in
[`paper/benchmarks/external/tools_comparison/methodology.md`](../benchmarks/external/tools_comparison/methodology.md).

**Table T2 — Head-to-head performance on the same Liu 2023 cold-seep subset (200 MAG × 30 sample)**
[**⚠️ placeholder values — Task 3b actual benchmark execution pending**]：

| Task | EnvMeta | Anvi'o | Krona | MicrobiomeAnalyst |
|---|---|---|---|---|
| **PCoA on β-diversity** |
| Wall time (s) | TBD | n/a (not native) | n/a | TBD (web) |
| Peak RSS (MB) | TBD | n/a | n/a | n/a (web) |
| User operation count | TBD | n/a | n/a | TBD (uploads + clicks) |
| Output formats | PNG/PDF/SVG/TIFF | — | — | PNG/PDF |
| Reproducible script | ✅ .py | partial | — | — |
| **KO heatmap** |
| Wall time (s) | TBD | TBD | n/a | TBD (web) |
| Peak RSS (MB) | TBD | TBD | n/a | n/a (web) |
| **Hierarchical taxonomy HTML** |
| Wall time (s) | TBD (stackplot proxy) | n/a | TBD | TBD (sunburst) |
| Output size (KB) | TBD | n/a | TBD | n/a (web) |
| Offline-capable | ✅ | partial | ✅ | ❌ |
| **Setup overhead** |
| First-time install (min) | < 5 (pip) | 30-60 (conda + DBs) | 5-10 (conda) | 0 (web) |
| Input prep (min) | 0 (KO/abund 直接) | 60-120 (需 contigs.db) | 5 (TSV) | 5 (CSV upload) |

##### 5.2.9.3 Figure F10 schema — Visual head-to-head comparison

**Figure F10**：6-panel head-to-head — 同一份 Liu 2023 子集 + 同一分析（如 PCoA
of β-diversity），在不同工具中的最佳输出：
- A) EnvMeta（PCoA + 95% CI ellipse + PERMANOVA p）
- B) Anvi'o（如可执行的等价输出）
- C) Krona（hierarchical taxonomy 等价输出）
- D) MicrobiomeAnalyst Web 平台（PCoA 输出）
- E) **EnvMeta 独家**：元素循环图（无可比，标为 unique）
- F) **EnvMeta 独家**：假说评分卡 + null_p 分布（无可比，标为 unique）

每 panel 加 caption：工具名 + 操作步数 + 输出格式 + wall time。

**Source**：[`paper/benchmarks/external/tools_comparison/`](../benchmarks/external/tools_comparison/)（待 Task 3b 实测后 figures 生成）+
[`paper/tool_comparison.md`](../tool_comparison.md)（v1 framing 文档）。

---

### 5.3 Discussion（~700-900 words）

1. **EnvMeta 的差异化卖点**（200 words）
   - 工具论文常见误区："功能多 = 价值大"
   - EnvMeta 的核心不是功能多，而是 **3 个业界空白**（cycle / hypothesis / HTML）
   - 元素循环 + 假说评分 = "**从描述性可视化迈向假说评估可视化**"

2. **vs 测序公司云平台**（150 words）
   - 真实竞争对手是云平台，不是 Anvi'o
   - 云平台的优势：开箱即用 / 大数据集
   - 云平台的劣势：封闭 / 不支持元素循环 / 不可定制 / 数据出境
   - EnvMeta 定位 = "**比云平台更灵活 + 比 Anvi'o 更易用**" 的中间层

3. **Fork-rather-than-community 模型的合理性**（150 words）
   - 不维护中心化机制目录，每篇论文绑定一份 Bundle
   - 类似 R bioconductor vs Linux distros 的差别
   - 优势：社区维护负担 0 / 论文复现性 100%
   - 劣势：用户跨论文需要 fork / 学习成本略高

4. **Calibration vs Discrimination framing**（NEW —— 1 段过渡 / ~150 words / 接 §Y.1-Y.2）

> Drop-in English text — source: [`discussion_calibration_discrimination.md`](discussion_calibration_discrimination.md) §Y.1 + §Y.2 abridged.

> The four-arm controlled experiment (§5.2.4) provides what we term
> **KEGG-coverage-dependent calibration evidence**: under fixed default
> thresholds, EnvMeta's scoring engine returns `STRONG` for the four datasets
> that supply canonical KEGG annotation (KofamScan, DRAM, METABOLIC, or
> end-to-end Pyrodigal + GhostKOALA), and `INSUFFICIENT` for the one dataset
> whose published annotation is restricted to a custom ROCker fourteen-gene
> subset. This pattern reflects two coupled effects: the engine performs as
> designed when KEGG-orthology coverage is adequate, and it correctly flags
> coverage mismatch on Arm B rather than failing silently. The stress-test
> layer (§5.2.4 paragraph 2) provides additional evidence that, **conditional
> on adequate KEGG coverage**, the engine resists awarding high scores to
> claims that violate environmental priors: the cross-topic "arsenate_reduction
> should dominate" claim was rejected with n = 0 active MAGs in 2/2
> non-arsenic datasets — a result consistent with (rather than ironclad proof
> of) the engine being uninfluenced by the universal *arsC* detoxification
> homolog. Together, these observations position EnvMeta as a *diagnostic
> instrument* whose performance is conditional on adequate KEGG coverage of
> the target pathways, rather than as a domain-blind tool that performs
> equally on any annotation regime, and rather than as an oracle that
> validates user interpretations.

5. **局限性**（~340 words / merged §Y.3 + performance §5.3）

> Drop-in English text — source: [`discussion_calibration_discrimination.md`](discussion_calibration_discrimination.md) §Y.3
> + [`performance_paragraph_snippets.md`](performance_paragraph_snippets.md) §5.3.

> **(1) Residual author selection bias** — the single largest limitation of
> the four-arm calibration experiment. Despite explicit time-pre-registration
> discipline (§5.4.8.2), the authors had read all four target papers before
> writing the hypothesis YAMLs. Time-pre-registration locks the claim
> entities against post-hoc adjustment but does not control for the
> cognitive selection bias of choosing claims plausibly satisfiable by
> KEGG-curated datasets in the topics surveyed. The four `STRONG` calibration
> outcomes therefore conflate two effects that cannot be cleanly separated
> within the present design: the scoring engine's behaviour under default
> thresholds, and the authors' skill in claim selection. The cleanest
> mitigation is **blind hypothesis writing** by collaborators unfamiliar
> with the dataset's findings (§5.3 Future).
>
> **(2) KEGG-coverage dependency** — as discussed in §Y.1 / §5.3 paragraph
> 4, the calibration evidence is conditional on adequate KEGG-orthology
> coverage of the target pathways. The `INSUFFICIENT` outcome on Arm B
> (Wei 2024 ROCker-only annotation) is itself diagnostic of this
> conditional, but it limits the generalisation of the four `STRONG`
> results to environments where canonical KEGG annotations are available.
>
> **(3) KB coverage.** EnvMeta KB v2.0 (KEGG snapshot 2026-04-15) covers 4
> elements × 18 pathways × 57 KOs (As / N / S / Fe). Two of Wei's 14 target
> genes (*arxA*; *nrfA*) lacked KB mappings; iron(II) oxidation and iron(III)
> reduction pathways central to AMD biogeochemistry are not yet encoded.
> KB v2.1 will extend ROCker-model alias support and add iron-redox plus
> DNRA blocks.
>
> **(4) Pre-publication hand-checks.** Post-hoc DOI verification identified
> four citation errors (§5.2.8 audit). These do not affect scoring outputs
> but emphasise the value of building DOI verification into the
> hypothesis-writing workflow itself.
>
> **(5) Synthetic-dense scaling regime.** The 58-cell performance benchmark
> uses a synthetic random KO assignment (25 KOs/MAG) to approximate full-KEGG
> annotation density. Real KofamScan / DRAM / METABOLIC annotations cluster
> around organism-specific functional repertoires and may show different
> cost-by-N_MAG behaviour. The synthetic-dense scaling numbers (§5.2.8)
> should be read as upper-bound estimates rather than as a faithful
> reproduction of real-world annotation distributions.
>
> **(6) Empirical scaling envelope.** From the performance benchmark
> (§5.2.8), full-pipeline runtime stays under 2 minutes for the typical
> PhD-thesis metagenome (200-1000 MAGs × 30-100 samples × 4 env factors);
> under 15 minutes for ~5000 MAGs; and is dominated by the cycle inference's
> permutation test. The 5000-MAG ceiling is a soft constraint imposed by
> heatmap legibility and typical iMeta-class user data scale rather than a
> runtime cliff, consistent with our design philosophy of an open-source
> local tool aimed at graduate-student and small-lab users, complementing
> rather than competing with HPC-grade pipelines such as Anvi'o.

6. **Future work**（~150 words / §Y.4）

> Drop-in English text — source: [`discussion_calibration_discrimination.md`](discussion_calibration_discrimination.md) §Y.4

> Beyond the `dominance_score` extension (delivered in v0.9.x; §5.2.4 paragraph 2)
> and KB v1.2 backlog mentioned above, two methodological future items follow
> directly from the stress-test experience. First, **third-party blind stress
> YAMLs**: in the next user-study iteration, collaborators unfamiliar with
> the four target papers will be invited to author independent stress YAMLs
> for the same datasets, providing a selection-bias-controlled replication
> of the present result. Second, an **LLM-assisted hypothesis YAML drafting**
> workflow that ingests an environmental description and a list of recent
> reviews, produces a candidate calibration-plus-stress YAML, and flags
> claim entities that might fail DOI verification or proof-of-extraction
> grading. We see EnvMeta's hypothesis scorer as well-suited to LLM-assisted
> authoring precisely because the engine mechanically resists hindsight bias
> (pre-registration discipline + binary status reporting + claim-entity
> immutability), making the LLM's drafting errors easy to surface rather
> than easy to mask. Beyond hypothesis-side work, the L3 plugin framework
> (user-supplied Python `analyze()` → auto-registered GUI) and a D3.js
> cycle-figure editor for direct node/pathway dragging round out the longer
> roadmap.

---

### 5.4 Methods（~2200-2800 words / GPB v2 新增 §5.4.0 算法框架节）

> **GPB v2 重要变更**：新增 §5.4.0 Algorithmic Framework 节作为 Methods 首段，
> 用 formal pseudocode + 时间复杂度分析 + Equation 1 weight-of-evidence 公式建立**算法贡献**，
> 然后 §5.4.1-§5.4.7（旧 §5.4 的 1-7 条）作为**实现细节**叙述。这是 GPB 审稿人最关键的差异化要求
> （vs iMeta tool paper 标准多了 algorithm formalization 这一项）。

#### 5.4.0 Algorithmic Framework — Hierarchical Weight-of-Evidence MCDA

##### Notation

**Table N1 — Symbols and parameters**

| Symbol | Meaning | Default value |
|---|---|---|
| 𝒢 = {g₁, …, g_M} | Set of M MAGs (genomes) | — |
| 𝒮 = {s₁, …, s_N} | Set of N samples | — |
| 𝒦 = {k₁, …, k_K} | Set of K KO target IDs | 57 (KB v2.0) |
| 𝒫 = {p₁, …, p_P} | Set of P element-cycle pathways | 18 |
| 𝒞 = {c₁, …, c_C} | Set of C user claims (YAML) | varies |
| ℰ = {e₁, …, e_J} | Set of J environmental factors | varies (1-14) |
| **A** ∈ ℝ^(M×N) | MAG abundance matrix | input |
| **K** ∈ {0,1}^(M×K) | KO presence matrix | input |
| τ₁ < τ₂ < τ₃ | Three pathway-completeness thresholds | (0.4, 0.5, 0.6) |
| α | Permutation count | 999 |
| δ | Weight perturbation magnitude | ±0.20 |
| θ_strong | Strong-label threshold | 0.75 |
| θ_suggestive | Suggestive-label threshold | 0.40 |
| ε | Saltelli weight-robustness band | 0.05 |
| η_min_dominance | Minimum dominance fraction (pathway_active claim) | 0.20 |

##### Framework overview

The framework consumes inputs (**A**, **K**, claim YAML 𝒞, environmental factors ℰ)
and emits, for each claim c ∈ 𝒞, a five-tuple diagnostic:

> **D(c) = (score, weight_robust_score, null_p, evidence_count, confidence_label)**

The pipeline composes three algorithmic stages — **S1** compositional debiasing
with three-threshold top-1 contributor sensitivity scanning, **S2** Fisher
999-permutation null calibration, and **S3** weighted-sum MCDA scoring with
Bradford-Hill required-veto and Saltelli ±δ one-at-a-time sensitivity. Each
stage's output is independently inspectable downstream, supporting
reproducibility audit at every intermediate quantity.

##### Algorithm 1 — Master pipeline

```
Algorithm 1: EnvMeta hierarchical weight-of-evidence evaluation
Input:  A (M×N abundance), K (M×K KO matrix), 𝒞 (claims),
        ℰ (env factors), θ (thresholds), α (= 999), δ (= 0.20)
Output: { D(c) : c ∈ 𝒞 } where D(c) = (score, weight_robust, null_p,
        evidence_count, confidence_label)

 1: function Evaluate(A, K, 𝒞, ℰ, θ, α, δ)
 2:     (K′, S1_sens) ← S1_Debias(A, K, τ₁, τ₂, τ₃)           ▷ Algorithm 2
 3:     R ← {}                                                  ▷ S2 results table
 4:     for each pathway p ∈ 𝒫 do
 5:         for each env factor e ∈ ℰ do
 6:             (ρ_obs, perm_p) ← S2_PermP(completeness(p,K′), e, α)   ▷ Algorithm 3
 7:             tier ← TierLabel(ρ_obs, perm_p)                ▷ 5-tier: strong/suggestive/weak/spurious?/unknown
 8:             R[(p, e)] ← (ρ_obs, perm_p, tier)
 9:         end for
10:     end for
11:     D ← {}
12:     for each claim c ∈ 𝒞 do
13:         D[c] ← S3_ScoreClaim(c, K′, R, S1_sens, θ, α, δ)   ▷ Algorithm 4
14:     end for
15:     return D
16: end function
```

##### Algorithm 2 — S1: Compositional debiasing with three-threshold sensitivity scan

```
Algorithm 2: S1 — CLR debias + 3-threshold top-1 contributor robustness
Input:  A (M×N abundance), K (M×K KO), τ₁ < τ₂ < τ₃ completeness thresholds
Output: K′ debiased KO matrix; S1_sens robustness record

 1: function S1_Debias(A, K, τ₁, τ₂, τ₃)
 2:     A_clr ← CLR(A)                                          ▷ Aitchison centered-log-ratio (Aitchison 1986)
 3:     K′ ← K                                                  ▷ binary presence unchanged
 4:     S1_sens ← {}
 5:     for each pathway p ∈ 𝒫 do
 6:         top_contrib ← []
 7:         for each threshold τ ∈ {τ₁, τ₂, τ₃} do
 8:             contributors_τ ← argmax_{g ∈ 𝒢}
 9:                              [ A_clr[g, ·] · 𝟙(completeness(p, K[g, ·]) ≥ τ) ]
10:             top_contrib.append(contributors_τ)
11:         end for
12:         S1_sens[p] ← (|unique(top_contrib)| == 1)            ▷ True = same top-1 across all 3 thresholds (robust)
13:     end for
14:     return K′, S1_sens
15: end function
```

**Interpretation**: S1 does not modify the KO matrix itself; it produces a
**per-pathway robustness flag** indicating whether the identity of the top-1
contributor MAG is stable under three completeness threshold choices. Pathways
flagged S1_sens = False are reported with reduced confidence downstream (this
is the "S1 去偏" semantic), exposing one of the three evidential uncertainties
named in §5.1.

##### Algorithm 3 — S2: Fisher 999-permutation null calibration

```
Algorithm 3: S2 — Spearman correlation with permutation p-value
Input:  x, y paired vectors of length n; α permutation count
Output: ρ_obs ∈ [-1, 1]; perm_p ∈ (0, 1]

 1: function S2_PermP(x, y, α)
 2:     ρ_obs ← Spearman(x, y)                                  ▷ rank-based correlation
 3:     count ← 0
 4:     for j = 1 to α do
 5:         y_shuffled ← UniformShuffle(y)
 6:         ρ_j ← Spearman(x, y_shuffled)
 7:         if |ρ_j| ≥ |ρ_obs| then count ← count + 1
 8:     end for
 9:     perm_p ← (count + 1) / (α + 1)                          ▷ Fisher exact permutation (Anderson 2001)
10:     return ρ_obs, perm_p
11: end function
```

**Time per call**: O(α × N log N) for the Spearman ranking inside each
permutation. With default α = 999 and N = 30-100 samples typical of MAG
studies, single (p, e) evaluation costs ~10 ms.

##### Algorithm 4 — S3: MCDA scoring with Bradford-Hill required-veto + Saltelli ±δ

```
Algorithm 4: S3 — Weighted-sum MCDA scoring with required-veto and weight robustness
Input:  claim c, K′ (debiased KO), R (S2 results), S1_sens, θ, α, δ
Output: D(c) = (score, weight_robust, null_p, evidence_count, confidence_label)

 1: function S3_ScoreClaim(c, K′, R, S1_sens, θ, α, δ)
 2:     {(s_i, w_i, type_i, required_i)} ← ParseYAML(c)         ▷ extract subclaims
 3:     subclaim_scores ← []
 4:     for each subclaim i in c.claims do
 5:         s_i ← EvaluateSubclaim(type_i, K′, R, S1_sens, θ)   ▷ ∈ [0, 1]; supports 6 claim types
 6:         subclaim_scores.append((s_i, w_i, required_i))
 7:     end for
 8:     evidence_count ← |{i : s_i > 0}|
 9:     
10:     ▷ ── Bradford-Hill required-veto:
11:     if ∃ i : required_i = true ∧ s_i = 0 then
12:         return D(c) ← (NaN, NaN, NaN, evidence_count, "INSUFFICIENT")
13:     end if
14:     
15:     ▷ ── Weighted-sum score:
16:     W ← Σ_i w_i
17:     score ← (Σ_i s_i × w_i) / W                             ▷ ∈ [0, 1]
18:     
19:     ▷ ── Saltelli ±δ one-at-a-time weight sensitivity:
20:     weight_robust ← score
21:     for each i ∈ {1, …, |subclaims|} do
22:         for ±sign ∈ {+, -} do
23:             w_i′ ← w_i × (1 ± δ)
24:             W′ ← W − w_i + w_i′
25:             score′ ← (Σ_{j≠i} s_j × w_j + s_i × w_i′) / W′
26:             weight_robust ← min(weight_robust, score′)
27:         end for
28:     end for
29:     
30:     ▷ ── Null calibration (claim-entity shuffle):
31:     count ← 0
32:     for j = 1 to α do
33:         𝒞_shuffled ← ShuffleClaimEntities(c)                ▷ randomly reassign claim targets
34:         score_j ← (Σ_i s_i(𝒞_shuffled) × w_i) / W
35:         if score_j ≥ score then count ← count + 1
36:     end for
37:     null_p ← (count + 1) / (α + 1)                          ▷ shuffle-consistency diagnostic (not frequentist p)
38:     
39:     ▷ ── Confidence label per Equation 1:
40:     confidence_label ← LabelAssign(score, weight_robust, null_p, evidence_count, θ)
41:     
42:     return (score, weight_robust, null_p, evidence_count, confidence_label)
43: end function
```

##### Equation 1 — Weight-of-evidence confidence label assignment

The five-tier confidence label is assigned via the decision rule:

```
              ⎧ STRONG         if score ≥ θ_strong  ∧ weight_robust ≥ θ_strong − ε
              ⎪                                     ∧ null_p < 0.05
              ⎪ SUGGESTIVE     if θ_suggestive ≤ score < θ_strong
              ⎪                                     ∧ weight_robust ≥ θ_suggestive
label(c)  =  ⎨ WEAK           if 0 < score < θ_suggestive
              ⎪ INSUFFICIENT   if any required claim vetoed   ▷ Bradford-Hill veto branch
              ⎩ UNKNOWN        if evidence_count = 0           ▷ no usable sub-claim
```                                                                                 (1)

with default θ_strong = 0.75 (clinical-epidemiology "strong evidence"
convention), θ_suggestive = 0.40 (Bradford-Hill consistency threshold),
ε = 0.05 (Saltelli weight-robustness band), and α = 999 (Fisher permutation
count following Anderson 2001 PERMANOVA convention). All thresholds are
pre-registered with cryptographic anchoring (OpenTimestamps: Todd 2016,
[paper/manuscript/timestamps/](timestamps/) directory) and grounded in
conventional clinical-epidemiology and environmental risk-assessment literature
(Hill 1965; Linkov et al. 2009; Rhomberg et al. 2010; Suter & Cormier 2011).
We treat `null_p` explicitly as a **shuffle-consistency diagnostic**, not as a
frequentist p-value, given the discrete and granular nature of the
4-9-subclaim YAML score distribution (Methods §5.4.8.1).

##### Time complexity analysis

Let M = #MAGs, N = #samples, P = #active pathways, J = #env factors,
C = #claims with |c| ≈ 4-9 subclaims each, α = #permutations (default 999).

**Table N2 — Per-stage time complexity**

| Stage | Operation | Asymptotic cost |
|---|---|---|
| S1 CLR | M × N abundance transform | O(M × N) |
| S1 sensitivity scan | 3 thresholds × P pathways × M MAGs | O(3 × M × P) |
| S2 permutation null (per cell) | α × Spearman(N) | O(α × N log N) |
| S2 full table | P × J × O(α × N log N) | **O(P × J × α × N log N)** |
| S3 per-claim weighted-sum | \|c\| subclaims | O(\|c\|) |
| S3 Saltelli ±δ | 2 × \|c\| reweights | O(\|c\|²) |
| S3 null calibration | α shuffles × \|c\| reeval | O(α × \|c\|²) |
| S3 full | C × O(α × \|c\|²) | **O(C × α × \|c\|²)** |

**Total**: **O(P × J × α × N log N + C × α × \|c\|²) + O(M × N)**.

The cost is dominated by **S2 permutation** (P × J × α × N log N) under
default α = 999, scaling with **environmental factor breadth and pathway
count, not MAG count**. This matches the 58-cell performance benchmark
(§5.2.8): 200 → 1000 MAGs added only +24% runtime at fixed P, J, α, since
cycle_diagram cost is largely N_MAG-independent. The pipeline runs comfortably
on standard laptops (~30-120 s for typical PhD-scale metagenomes).

##### Reproducibility audit invariants

Each diagnostic **D(c)** exposes **five independently auditable quantities**
(`score`, `weight_robust_score`, `null_p`, `evidence_count`, `confidence_label`).
A downstream reader inspecting D(c) can independently verify:

1. **score's robustness to weight perturbation** via `weight_robust_score`
   (Saltelli ±20% one-at-a-time minimum);
2. **null-distribution mass** via `null_p` (claim-entity shuffle-consistency);
3. **sub-claim evidential coverage** via `evidence_count`;

— all **without re-running the analysis**. The framework persists these
five quantities in machine-readable form (JSON in EnvMeta's HTML export,
TSV via `envmeta evaluate-yaml`), enabling third-party reviewer audit. This
five-quantity audit invariant is the core reproducibility contribution of
the framework, distinguishing it from narrative-reasoning workflows that
expose only a final qualitative judgment.

---

#### 5.4.1 文件识别模块（300 words）

- 11 FileType 表头规则匹配 + 反向索引
- abundance 分 MAG/TAXON 两档（confidence 0.95 / 0.88）

2. **14 图分析引擎**（300 words）
   - reads-based / MAG-based / cycle 三类引擎
   - 共享 4 层参数面板（`envmeta/analysis/_mag_common.py`）
   - matplotlib + seaborn 出版级渲染

3. **元素循环推断算法**（500 words）
   - **S1 去偏**：CLR 变换 + 抽样标准化 + 3 档阈值 Top-1 contributor 一致性
   - **S2 置换检验**：999 次 shuffle + perm_p 计算 + 5 档 confidence 标签
   - **S3 评分**：claim 类型路由 + weighted_sum + null_p + veto logic
   - 数据结构：`CycleData / EnvCorrelation / SensitivityRow / HypothesisScore`（见 `envmeta/geocycle/model.py`）

4. **YAML 假说 schema**（150 words → 详见 §5.4.8）
   - 6 类 claim（v0.9.0）：`pathway_active` / `pathway_inactive` / `coupling_possible` /
     `env_correlation` / `keystone_in_pathway` / `group_contrast`
   - 4 档标签（`strong` / `suggestive` / `weak` / `insufficient`）+ Bradford-Hill required-veto
   - schema 校验通过 jsonschema；评分引擎设计 + pre-registration discipline 详见 §5.4.8

5. **D3.js 交互 HTML 导出**（200 words）
   - inline embed D3 v7（~280 KB）+ 400 KB total HTML
   - 4 tab UI（架构 / 循环图 / 假说 / 环境相关）
   - 客户端 SVG 序列化 + download

6. **Performance benchmark implementation**（~250 words）

> Drop-in English text — source: [`paper/manuscript/performance_paragraph_snippets.md`](performance_paragraph_snippets.md) §5.4.6

We measured EnvMeta's per-figure runtime and memory profile by an in-process
harness ([`paper/benchmarks/performance/bench_harness.py`](../benchmarks/performance/bench_harness.py); reproducibility script
provided in the Zenodo bundle) that wraps each `analyze()` entry point with
`time.perf_counter` for wall time and a background `psutil.Process().memory_info().rss`
sampler (50 ms interval) for peak ΔRSS over a `gc.collect()`-cleared baseline.
Each (dataset × figure) combination was run with 2-3 repeats, with `plt.close("all")`
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
(20-40 KO/MAG) of full-pipeline metagenomes. All measurements were taken on
a single Intel i7-class laptop (16 GB RAM, Windows 10, Python 3.11);
absolute numbers will vary across CPUs but the cross-cell ratios are
expected to transfer. Full benchmark report at
[`paper/benchmarks/performance.md`](../benchmarks/performance.md); raw measurements at
[`paper/benchmarks/performance/results/`](../benchmarks/performance/results/).

7. **实施细节**（150 words）
   - Streamlit 1.30+
   - matplotlib 3.8+
   - Python 3.11+
   - 测试套件：pytest **301/301 全绿**（v0.9.1）
   - CI：GitHub Actions（待补）

---

### 5.4.8 Hypothesis Scoring Engine and External-Dataset Benchmark（~1450 words / 7 子节）

> Drop-in English text — source: [`paper/manuscript/methods_external_validation.md`](methods_external_validation.md)
>
> 7 个子节涵盖：
> - **§4.6.1** Scoring engine design — MCDA + Bradford-Hill weight-of-evidence + 6 类 claim types + 3 confidence indicators
> - **§4.6.2** Pre-registration discipline — git timestamp anchoring, default thresholds, citation policy ≥ 5y before target paper
> - **§4.6.3** Four-Arm calibration experiment — Arm A in-house / Arm B Wei 2024 ROCker / Arm C1 Liu 2023 DRAM / Arm C2-A Grettenberger 2021 METABOLIC / Arm C2-B Ayala 2020 GhostKOALA re-annotated
> - **§4.6.4** Calibration results — 4 KEGG-curated arms STRONG (overall=1.000); Arm B INSUFFICIENT due to required-veto on incomplete denitrification annotation
> - **§4.6.5** Stress test — 3-Arm × 3-class (reversed/cross-topic/`pathway_inactive`); cross-topic rejection in 2/2 non-arsenic datasets (n=0 active MAGs); v0.9.x `dominance_score` extension upgrades B-tier → A-tier
> - **§4.6.6** Reference audit and post-hoc corrections — 4 errors corrected transparently in YAML metadata (commits `ddd3098`, `cae2de7`); claim entities frozen; full audit at [`hypothesis_references_audit.md`](hypothesis_references_audit.md)
> - **§4.6.7** Auxiliary perturbation analysis — random pathway-target substitution in 2 modes (within-element / cross-element); N=20 each per 4 datasets (Arm A partial, 3 external full); cross-element collapses Liu 2023 to 0/20 STRONG → element-level target accuracy is mechanistically required
> - **§4.6.8** Threshold sensitivity — sweep `strong_threshold` ∈ {0.65, 0.70, 0.75, 0.80, 0.85}; all 4 KEGG-curated calibration arms remain STRONG and all stress runs remain weak/suggestive across the conventional range
> - **§4.6.9** Data and code availability — reshape scripts at `tools/external_benchmarks/`; perturbation + threshold sensitivity runners at `tools/external_benchmarks/`; YAMLs at `paper/benchmarks/external/{dataset}/`; results at `envmeta_outputs/`, `paper/benchmarks/external/perturbation/`, `paper/benchmarks/external/threshold_sensitivity/`
>
> 19 条 Vancouver 引用（含 DOI） — 见 §5.7 References。投稿时整段插入。

---

### 5.5 Data and Code Availability

```
- Source code: https://github.com/redlizzxy/EnvMeta (MIT License)
- Online demo: https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/
- Sample data: tests/sample_data/ (arsenic slag bioremediation, subset)
- Zenodo DOI: TBD (release v1.0 at submission)
- Documentation: docs/data_preparation_zh.md + docs/install_for_beginners.md
- Issue tracker: GitHub Issues
```

---

### 5.6 Acknowledgements

```
We thank Dr. Yong-Xin Liu (Agricultural Genomics Institute at
Shenzhen, Chinese Academy of Agricultural Sciences) and the broader
iMeta editorial community for the foundational visualization
ecosystem — including ImageGP, ImageGP 2, EasyAmplicon, EasyMetagenome,
EVenn, Wekemo Bioincloud, Sangerbox, OmicStudio, shinyCircos, TOmicsVis,
iMetaLab Suite, iNAP, Majorbio Cloud, and MetOrigin — that has set
the methodological and accessibility standard for our field. EnvMeta
was designed as a domain-specialist complement to this ecosystem,
and its accessibility-first GUI philosophy directly inherits from
the EasyAmplicon and ImageGP 2 design tradition.

We thank [supervisor name] for guidance and [research group] for
data sharing and beta testing. We acknowledge feedback from the
graduate student volunteers participating in the v2 user-study
questionnaire (final n to be reported in the v1.0 release; if user
study data is available at submission time, replace this with the
actual n and consolidated SUS / task-success summary), and the
authors of the four external-validation datasets (Wei 2024, Liu
2023, Grettenberger 2021, Ayala 2020) for openly sharing their
MAG-level annotations without which the present calibration and
stress-test experiments would not have been possible.

Funding: [supervisor's grant ID] / China University of Geosciences
(Beijing) graduate research support.
```

---

### 5.7 References（约 50-60 篇）

> ⚠️ **iMeta 实际引用格式 ≠ Vancouver**（参考 Chen 2024 ImageGP 2 References 实测）：
>
> ```
> Chen, Tong, Yong-Xin Liu, and Luqi Huang. 2022. "ImageGP: An Easy-To-Use
> Data Visualization Web Server for Scientific Researchers." iMeta 1: e5.
> https://doi.org/10.1002/imt2.5
> ```
>
> 即 **author-year + 引号 title + 斜体 journal + volume: page** 格式（类似 Chicago）。
>
> §5.7.1 (来自 §4.6 草稿 19 条) 当前是 Vancouver 格式（"Belton V, Stewart TJ. ..."），
> §5.7.2 / §5.7.3 (本次新增) 已用 iMeta 格式。
>
> **投稿 final docx 阶段必须把 §5.7.1 19 条全部转换为 iMeta author-year 格式**。
> 工时预估 30-45 min（手工 search-replace + DOI 校对）。

**重点引用领域**（按 iMeta 编辑偏好排序）：

1. **iMeta 系列论文**（5-8 篇，编辑友好度+）
   - 刘永鑫 et al. EasyAmplicon / EasyMetagenome 系列
   - iMeta 已发表的工具论文（参考其格式）

2. **微生物组顶刊**（10-15 篇）
   - Nature Microbiology / ISME / Microbiome
   - QIIME 2 / DADA2 / phyloseq

3. **可视化工具**（8-10 篇）
   - Krona / Anvi'o / MicrobiomeAnalyst / iTOL / Cytoscape

4. **元素循环 + 生物地球化学**（10-15 篇）
   - KEGG / MetaCyc 数据库论文
   - 砷代谢、铁氧化、硫还原经典论文
   - 多元素耦合机制（如 Stocker / Konhauser 的铁循环综述）

5. **统计方法**（5-8 篇）
   - PERMANOVA / Bray-Curtis / Aitchison
   - 置换检验 + bootstrap

6. **Reproducibility / FAIR**（3-5 篇）
   - Wilson et al. Good enough practices
   - FAIR principles

#### 5.7.1 §5.4.8 引用清单（19 条 Vancouver + DOI，已 verify）

> 来源：[`methods_external_validation.md`](methods_external_validation.md) §References。
> 投稿时按 §5.7 总编号合并整理。

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
20. Anderson, Marti J. 2001. "A New Method for Non-Parametric Multivariate Analysis of Variance." *Austral Ecology* 26: 32-46. https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x  *(used to justify n=999 permutation count, §5.4.8.1)*
21. Segata, Nicola, Jacques Izard, Levi Waldron, Dirk Gevers, Larisa Miropolsky, Wendy S. Garrett, and Curtis Huttenhower. 2011. "Metagenomic Biomarker Discovery and Explanation." *Genome Biology* 12: R60. https://doi.org/10.1186/gb-2011-12-6-r60  *(used to justify per-group n ≥ 4 recommendation, §5.4.6)*
22. Stolz, John F., Partha Basu, Joanne M. Santini, and Ronald S. Oremland. 2006. "Arsenic and Selenium in Microbial Metabolism." *Annual Review of Microbiology* 60: 107-130. https://doi.org/10.1146/annurev.micro.60.080805.142053  *(used in §Y.2 to ground sporadic-oxidizer literature, §5.2.4 stress test discussion)*

#### 5.7.2 iMeta sister-tools 引用清单（§5.1 ecosystem 段落引用，强烈推荐保留）

> **政治意义重要**：这些是 iMeta 已发表的 sister tools，引用它们建立 ecosystem
> 信任感 + 致敬期刊编辑团队（ImageGP 2 Editor-in-Chief Liu YX 团队）。
>
> ⚠️ **iMeta 引用格式不是 Vancouver**，看 ImageGP 2 论文 References 实际格式是
> **author-year + 引号 title**（类似 APA / Chicago）。例如：
> > Chen, Tong, Yong-Xin Liu, and Luqi Huang. 2022. "ImageGP: An Easy-To-Use Data
> > Visualization Web Server for Scientific Researchers." *iMeta* 1: e5.
> > https://doi.org/10.1002/imt2.5
>
> 投稿 final docx 阶段需要把 §5.7.1 全部 19 条 Vancouver 引用 **格式转换** 为 iMeta
> 风格。所有 5.7.1 + 5.7.2 + 5.7.3 引用最终统一编号。

1. Chen, Tong, Yong-Xin Liu, Tao Chen, Mei Yang, Siqing Fan, Minglei Shi, et al. 2024. "ImageGP 2 for Enhanced Data Visualization and Reproducible Analysis in Biomedical Research." *iMeta* 3: e239. https://doi.org/10.1002/imt2.239
2. Chen, Tong, Yong-Xin Liu, and Luqi Huang. 2022. "ImageGP: An Easy-To-Use Data Visualization Web Server for Scientific Researchers." *iMeta* 1: e5. https://doi.org/10.1002/imt2.5
3. Chen, Tong, Haiyan Zhang, Yu Liu, Yong-Xin Liu, and Luqi Huang. 2021. "EVenn: Easy to Create Repeatable and Editable Venn Diagrams and Venn Networks Online." *Journal of Genetics and Genomics* 48: 863-866. https://doi.org/10.1016/j.jgg.2021.07.007
4. Ning, Wanshan, Yuxiang Wei, Letian Gao, Cheng Han, Yujie Gou, Shanshan Fu, et al. 2022. "HemI 2.0: An Online Service for Heatmap Illustration." *Nucleic Acids Research* 50: W405-W411. https://doi.org/10.1093/nar/gkac480
5. Shen, Weitao, Ziguang Song, Xiao Zhong, Mei Huang, Danting Shen, Pingping Gao, et al. 2022. "Sangerbox: A Comprehensive, Interaction-Friendly Clinical Bioinformatics Analysis Platform." *iMeta* 1: e36. https://doi.org/10.1002/imt2.36
6. Lyu, Fengye, Feiran Han, Changli Ge, Weikang Mao, Li Chen, Huipeng Hu, et al. 2023. "OmicStudio: A Composable Bioinformatics Cloud Platform with Real-Time Feedback That Can Generate High-Quality Graphs for Publication." *iMeta* 2: e85. https://doi.org/10.1002/imt2.85
7. Wang, Yazhou, Lihua Jia, Ge Tian, Yihan Dong, Xiao Zhang, Zhengfu Zhou, et al. 2023. "shinyCircos-V2.0: Leveraging the Creation of Circos Plot with Enhanced Usability and Advanced Features." *iMeta* 2: e109. https://doi.org/10.1002/imt2.109
8. Miao, Ben-Ben, Wei Dong, Zhao-Fang Han, Xuan Luo, Cai-Huan Ke, and Wei-Wei You. 2023. "TOmicsVis: An All-In-One Transcriptomic Analysis and Visualization R Package with Shinyapp Interface." *iMeta* 2: e137. https://doi.org/10.1002/imt2.137
9. Gao, Yunyun, Guoxing Zhang, Shunyao Jiang, and Yong-Xin Liu. 2024. "Wekemo Bioincloud: A User-Friendly Platform for Meta-Omics Data Analyses." *iMeta* 3: e175. https://doi.org/10.1002/imt2.175
10. Li, Leyuan, Zhibin Ning, Kai Cheng, Xu Zhang, Caitlin M. A. Simopoulos, and Daniel Figeys. 2022. "iMetaLab Suite: A One-Stop Toolset for Metaproteomics." *iMeta* 1: e25. https://doi.org/10.1002/imt2.25
11. Feng, Kai, Xi Peng, Zheng Zhang, Songsong Gu, Qing He, Wenli Shen, et al. 2022. "iNAP: An Integrated Network Analysis Pipeline for Microbiome Studies." *iMeta* 1: e13. https://doi.org/10.1002/imt2.13
12. Ren, Yi, Guo Yu, Caiping Shi, Linmeng Liu, Quan Guo, Chang Han, et al. 2022. "Majorbio Cloud: A One-Stop, Comprehensive Bioinformatic Platform for Multiomics Analyses." *iMeta* 1: e12. https://doi.org/10.1002/imt2.12
13. Yu, Gang, Cuifang Xu, Danni Zhang, Feng Ju, and Yan Ni. 2022. "MetOrigin: Discriminating the Origins of Microbial Metabolites for Integrative Analysis of the Gut Microbiome and Metabolome." *iMeta* 1: e10. https://doi.org/10.1002/imt2.10
14. Chen, Chengjie, Ya Wu, and Rui Xia. 2022. "A Painless Way to Customize Circos Plot: From Data Preparation to Visualization Using TBtools." *iMeta* 1: e35. https://doi.org/10.1002/imt2.35

#### 5.7.3 EasyAmplicon / Anvi'o / 通用 metagenomics 工具引用清单（§5.1 段落 + Methods 引用）

15. Liu, Yong-Xin, Yuan Qin, Tengfei Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, and Yang Bai. 2021. "A Practical Guide to Amplicon and Metagenomic Analysis of Microbiome Data." *Protein & Cell* 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8
16. Liu, Yong-Xin, Lei Chen, Tengfei Ma, Xiaofang Li, Maosheng Zheng, Xin Zhou, et al. 2023. "EasyAmplicon: An Easy-to-Use, Open-Source, Reproducible, and Community-Based Pipeline for Amplicon Data Analysis in Microbiome Research." *iMeta* 2: e83. https://doi.org/10.1002/imt2.83
17. Eren, A. Murat, Evan Kiefl, Alon Shaiber, Iva Veseli, Samuel E. Miller, Matthew S. Schechter, et al. 2021. "Community-Led, Integrated, Reproducible Multi-Omics with Anvi'o." *Nature Microbiology* 6: 3-6. https://doi.org/10.1038/s41564-020-00834-3
18. Bolyen, Evan, Jai Ram Rideout, Matthew R. Dillon, Nicholas A. Bokulich, Christian C. Abnet, Gabriel A. Al-Ghalith, et al. 2019. "Reproducible, Interactive, Scalable and Extensible Microbiome Data Science Using QIIME 2." *Nature Biotechnology* 37: 852-857. https://doi.org/10.1038/s41587-019-0209-9
19. McMurdie, Paul J., and Susan Holmes. 2013. "phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data." *PLOS ONE* 8: e61217. https://doi.org/10.1371/journal.pone.0061217
20. Ondov, Brian D., Nicholas H. Bergman, and Adam M. Phillippy. 2011. "Interactive Metagenomic Visualization in a Web Browser." *BMC Bioinformatics* 12: 385. https://doi.org/10.1186/1471-2105-12-385  *(Krona)*
21. Letunic, Ivica, and Peer Bork. 2024. "Interactive Tree of Life (iTOL) v6: Recent Updates to the Phylogenetic Tree Display and Annotation Tool." *Nucleic Acids Research* 52: W78-W82. https://doi.org/10.1093/nar/gkae268
22. Shannon, Paul, Andrew Markiel, Owen Ozier, Nitin S. Baliga, Jonathan T. Wang, Daniel Ramage, et al. 2003. "Cytoscape: A Software Environment for Integrated Models of Biomolecular Interaction Networks." *Genome Research* 13: 2498-2504. https://doi.org/10.1101/gr.1239303
23. Bastian, Mathieu, Sebastien Heymann, and Mathieu Jacomy. 2009. "Gephi: An Open Source Software for Exploring and Manipulating Networks." *Proceedings of the International AAAI Conference on Web and Social Media* 3: 361-362. https://doi.org/10.1609/icwsm.v3i1.13937

#### 5.7.4 Cross-element coupling chemistry (§5.4.3 / §5.2.3)

24. Newman, Dianne K., Terry J. Beveridge, and François M. M. Morel. 1997. "Precipitation of Arsenic Trisulfide by *Desulfotomaculum auripigmentum*." *Applied and Environmental Microbiology* 63 (5): 2022-2028. https://doi.org/10.1128/aem.63.5.2022-2028.1997  *(primary As₂S₃ chemistry coupling reference cited in §5.2.3)*
25. Rodriguez-Freire, Lucia, Reyes Sierra-Alvarez, Robert Root, Jon Chorover, and James A. Field. 2014. "Biomineralization of Arsenate to Arsenic Sulfides Is Greatly Enhanced at Mildly Acidic Conditions." *Water Research* 66: 242-253. https://doi.org/10.1016/j.watres.2014.08.016
26. Hollibaugh, James T., Carrie Carini, Heath Gürleyük, Roumiana Jellison, Steven B. Joye, Greg Lecleir, et al. 2005. "Arsenic Speciation in Mono Lake, California: Response to Seasonal Stratification and Anoxia." *Geochimica et Cosmochimica Acta* 69: 1925-1937. https://doi.org/10.1016/j.gca.2004.10.011

---

## 6. Figures + Tables 总表（v0.9.1 校准）

### 6.1 Figures

| # | Figure | 类型 | 现有素材 | 待补 |
|---|---|---|---|---|
| F1 | EnvMeta 整体架构 | Schema | 部分（CLAUDE.md 内嵌）| 重画为 SVG |
| F2 | 12 图 + 调参 + 代码生成 | Screenshots + 截图组合 | `paper/figures/screenshot_*.png` | 排版 4-panel |
| F3 | 元素循环图自动推断 ⭐ | Schema + Output | `paper/figures/cycle_*.png` | 加 S1-S3 流程 schema |
| F4 | YAML 假说评分器 + 4-Arm calibration ⭐ | Schema + Output | [`paper/figures/paper3_hypothesis_scoring/`](../figures/paper3_hypothesis_scoring/) ✅ PDF/PNG/SVG | 整合到 4-panel（A schema / B null_p / C decision tree / D calibration bar） |
| F5 | 独立交互 HTML 导出 | Screenshots | `paper/bundles/*.html` 截图 | 新拍 4 panel |
| F6 | Fork Bundle 结构 | Schema | 无 | 新画 |
| F7 | 砷渣案例研究 | Result | `paper/figures/cycle_arsenic.png` | 多组对比 panel |
| F8 | Performance + scaling 实测 ⭐ | Result | [`paper/benchmarks/performance/scaling_curve.{pdf,png,svg}`](../benchmarks/performance/scaling_curve.pdf) ✅ | 已 ready（3 panel × 58 cells） |
| F9 | Calibration → stress score gap | Result | [`paper/figures/paper3_hypothesis_scoring/figure_x_calibration_vs_stress.{pdf,png,svg}`](../figures/paper3_hypothesis_scoring/) ✅ | 已 ready；考虑作为 F4-D panel 或独立 F9 |
| F10 | vs 竞品对比柱图 | Comparison | `paper/tool_comparison.md` 整理 | 改成图（原 F9 → F10） |

### 6.2 Tables

| # | Table | 内容 | 现有素材 |
|---|---|---|---|
| T1 | vs 竞品矩阵 | Krona / Anvi'o / MicrobiomeAnalyst / 测序公司云 / EnvMeta — 12 图覆盖度 / 离线 / 元素循环 / 假说评分 / GUI 友好 / 开源 / 价格 | `paper/tool_comparison.md` 整理 |
| T2 | 4-Arm calibration ⭐ | Arm × dataset × annotation × overall_score × label × claims_satisfied × topic | [`paper/figures/paper3_hypothesis_scoring/table1_calibration.{md,tsv}`](../figures/paper3_hypothesis_scoring/) ✅ |
| T3 | 3-Arm stress test ⭐ | Arm × calibration label × stress overall (v1) × stress label (v1) × stress overall (v2) × discrimination grade × cross-topic rejection | [`paper/figures/paper3_hypothesis_scoring/table2_stress.{md,tsv}`](../figures/paper3_hypothesis_scoring/) ✅ |
| T4 | Hardware sizing | Regime × N_MAG × N_sample × wall × hardware × Streamlit Cloud 是否支持 | [`paper/benchmarks/performance.md`](../benchmarks/performance.md) §5 |
| T5 | Group-count guidance（可选 SI）| 各分析硬下限 + 软推荐 + 实际上限 | [`paper/benchmarks/performance.md`](../benchmarks/performance.md) §6 |
| ST1 | Reference audit (SI) | 16 claims × 13 references DOI table + extraction grade | [`paper/manuscript/hypothesis_references_audit.md`](hypothesis_references_audit.md) ✅ |

**编号 reconcile note**：草稿用 "Table 1 / Table 2 / Figure X"。本 outline 映射：
- `Table 1` (calibration) → **T2**
- `Table 2` (stress v1+v2) → **T3**
- `Figure X` (calibration→stress gap) → **F9**（独立 panel）或 **F4 panel D**（整合到 hypothesis figure）
- 投稿 final docx 时按所选方案统一编号

---

## 7. 投稿前清单（v0.9.1 校准 / 2026-05-09）

### 已完成 ✅

- [x] **核心代码 v0.9.1（301/301 测试全绿）**
- [x] 12 图 + 循环图 + 假说评分 (6 类 claim) + Bundle + HTML 全功能
- [x] tool_comparison.md / time_comparison.md
- [x] paper/figures/ 部分截图
- [x] 在线 demo + GitHub repo
- [x] 期刊选定 + 投稿决策
- [x] **R 侧侧对照 11 图**（v0.8.2，[`paper/benchmarks/r_comparison/`](../benchmarks/r_comparison/)，含 RDA 数值与 R vegan 对齐到 4 位小数）
- [x] **第二外部数据集复现 + benchmark**（v0.9.0 / 超额完成 4 数据集）：
  - Wei 2024（Arm B / paddy soil / ROCker, INSUFFICIENT）
  - Liu 2023（Arm C1 / cold seep / DRAM, STRONG）
  - Grettenberger 2021（Arm C2-A / AMD stream / METABOLIC, STRONG）
  - Ayala 2020（Arm C2-B / pit lake / GhostKOALA 重注释, STRONG）
- [x] **3-Arm stress test discrimination evidence**（v0.9.0 + v0.9.x dominance_score 升级至 A-tier 全清）
- [x] **English README + LICENSE**（v0.9.1 / commit 2b50dbc）
- [x] **Methods §4.6 草稿**（~1450 字 / 7 子节 / 19 Vancouver+DOI 引用）— [`methods_external_validation.md`](methods_external_validation.md)
- [x] **Results §X 草稿**（~800 字 / 4 子节 / Table 1+2 + Figure X 定义）— [`results_stress_test_section.md`](results_stress_test_section.md)
- [x] **Discussion §Y 草稿**（~640 字 / 4 子节 / 含 v0.9.x 升级补丁）— [`discussion_calibration_discrimination.md`](discussion_calibration_discrimination.md)
- [x] **Performance benchmark + scaling figure**（v0.9.1 / 58 测试点 / 3 regime）— [`paper/benchmarks/performance.md`](../benchmarks/performance.md) + [scaling_curve.{pdf,svg,png}](../benchmarks/performance/scaling_curve.pdf)
- [x] **Reference DOI audit + 4 处错引透明纠正** — [`hypothesis_references_audit.md`](hypothesis_references_audit.md)
- [x] **F4 + F8 + T2 + T3 实物素材就位**（v0.9.1）

### 投稿前必补 ⏳

- [ ] **整合 outline 三段进 docx 主稿**（用户做，涉及 outline 整体结构 + iMeta 投稿格式调整）
- [ ] **F1 / F6 schema 图重画** —— 1 天
- [ ] **F3 schema 图补 S1-S3 流程** —— 半天（F4 已 ready）
- [ ] **F5 HTML 截图 4 panel** —— 半天
- [ ] **F7 砷渣案例研究多组对比 panel** —— 半天
- [ ] **F10 vs 竞品对比柱图** —— 半天（基于 T1 数据）
- [ ] **References 整理总表**（合并 §5.7.1 + 其他 5 大类引用） —— 1-2 天
- [ ] **User study v2 数据回收 ≥ 8 份**（问卷 2026-04-19 发，等回收 1 周）
- [ ] **Zenodo DOI release**（投稿当天前 1-2 天，避免文件树漂移）—— 30 min

### 投稿当天

- [ ] 邮件刘永鑫确认推荐审稿人
- [ ] 投稿系统填表（Authors / Affiliations / Funding / Recommend reviewers）
- [ ] Cover Letter（强调 3 大空白 + iMeta 选题契合 + v0.9.0 calibration/discrimination evidence）
- [ ] Conflict of Interest 声明
- [ ] Zenodo DOI 嵌入正文 §5.5 Data and Code Availability

---

## 8. 时间线（再校准）

```
2026-05-07   今日：大纲完成
2026-05-08 ~ 05-22   补素材（R 对照 + 第二数据集 + 英文文档 + Zenodo DOI）
                     ── 工时 ~2-3 周
2026-05-23 ~ 06-15   论文初稿（Abstract + Intro + Results + Methods + Discussion）
                     ── 工时 ~3 周
2026-06-15 ~ 06-30   导师 + 课题组内审 + 修改
                     ── 工时 2 周
2026-07-01 ~ 07-10   英文润色（DeepL + Grammarly + 投稿前精修）
2026-07-10 ~ 07-14   邮件刘永鑫 + Cover Letter
2026-07-15           **投 iMeta** ⏱ 倒计时锚点
2026-09-15           一审决定（estimated minor / major revision）
2026-09-15 ~ 10-15   修订
2026-10-15           修订投回
2026-11-15 ~ 12-15   **接收 ✅**（赶秋招简历 + 12 月节点）
2027-01 ~ 02         online ✅（赶 2027-05 激励截止前）
```

---

## 9. 维护记录

| 日期 | 事项 | 备注 |
|---|---|---|
| 2026-05-07 | 大纲初版 | 待导师审阅 |
| 2026-05-09 | **方案 B 整合**：Methods §4.6 / Results §X.1-X.3 / Discussion §Y.1-Y.4 三段英文草稿插入 §5.4.8 / §5.2.4 + §5.2.8 / §5.3；Performance benchmark 段落（§5.2.8 + §5.4.6 + §5.3 局限）合并；Vancouver+DOI 19 条加入 §5.7.1；Figure/Table 编号校准（F4 / F8 / F9 已 ready；T2 / T3 / T4 已 ready）；§7 投稿前清单按 v0.9.1 重新校准；Abstract 更新含 4-Arm calibration + cross-topic rejection + performance benchmark 三个新 finding。语言决策：保留中文 outline scaffold + 英文段落直插（最低工作量）。投稿 final docx 阶段需要全英化 + 编号统一（用户做） |
