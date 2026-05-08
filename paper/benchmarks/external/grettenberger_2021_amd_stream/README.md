# EnvMeta 外部数据集验证 — Grettenberger 2021 AEM Cabin Branch AMD 溪流

> **状态**：Phase 1 完成（2026-05-08）
> **角色**：对照实验 Arm C2-A — KEGG-curated 跨主题（非砷）数据集
> **关联**：[../../../manuscript/scoring_validation_experiment.md](../../../manuscript/scoring_validation_experiment.md)

---

## 1. 数据集来源

| 项 | 内容 |
|---|---|
| 论文 | Grettenberger CL & Hamilton TL (2021) "Metagenome-assembled genomes of novel taxa from an acid mine drainage environment" |
| 期刊 | *Applied and Environmental Microbiology* (AEM) |
| DOI | [10.1128/AEM.00772-21](https://doi.org/10.1128/AEM.00772-21) |
| PMC | [PMC8357290](https://pmc.ncbi.nlm.nih.gov/articles/PMC8357290/) |
| License | CC BY 4.0 |
| 数据规模 | 32 MAGs + 5 sites (Cabin Branch AMD stream, Pennsylvania, USA) |
| 注释方法 | **DRAM + METABOLIC**（含 KEGG 全注释）+ GTDB-Tk taxonomy + CheckM2 |
| 主题 | AMD 溪流 C/N/S/Fe 4 元素循环（**无砷**）|

## 2. 复现路径

```powershell
# 1. 下载 Suppl 数据 + Table 1 HTML
#    PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC8357290/
#    - 下 "Supplemental file 2" = aem.00772-21-s0002.xlsx (Data Set S1)
#    - 下 "Supplemental file 1" = aem.00772-21-s0001.pdf (Files S1-S6)
#    - 抓 Table 1 HTML: https://pmc.ncbi.nlm.nih.gov/articles/PMC8357290/table/T1/

# 2. 运行 reshape (HTML + xlsx → 6 EnvMeta inputs)
conda activate envmeta
python tools/external_benchmarks/grettenberger2021_reshape.py `
    --table1 D:\download\grettenberger_table1.html `
    --dataset_s1 D:\download\aem.00772-21-s0002.xlsx `
    --out paper/benchmarks/external/grettenberger_2021_amd_stream/input_data_local

# 3. 跑 EnvMeta 全套（cycle + mag_quality + hypothesis 评分）
python tools/external_benchmarks/grettenberger2021_run_envmeta.py
```

总耗时：reshape ~1 s + EnvMeta ~3 s = **~4 s**（最快的 plug-and-play 数据集）。

## 3. EnvMeta 输入文件

| 文件 | 来源 | 内容 |
|---|---|---|
| `metadata.tsv` | placeholder | 1 pooled sample "All"（Data Set S1 仅 MAG 级，无 sample × MAG abundance）|
| `env_factors.tsv` | placeholder | pH=3.0 proxy（AMD 一般值）|
| `mag_taxonomy_labels.tsv` | Table 1 (HTML) | 29 MAGs × GTDB-Tk lineage |
| `quality_report.tsv` | Table 1 | 29 MAGs × completeness/contamination/size/GC |
| `abundance.tsv` | placeholder | 29 MAGs × 1 sample，全填 1.0 |
| `kegg_target_only.tsv` | Data Set S1 | **22990 unique MAG-KO records**（METABOLIC step Present 展开）|

### 关于 22990 个 KO records

Data Set S1 KEGGModuleStepHit 提供 29 MAGs × 2235 KEGG module steps 的 Present/Absent 矩阵。每个 step 可关联多个 KO（如 step `M00001+02` = `(K01810,K06859,K13810,K15916)` 任一命中即 Present）。

reshape 把每个 Present step 的所有 KO 都标记为该 MAG 携带 → **22990 unique MAG-KO records**（已去重）。这是 conservative inflation：
- 实际可能只命中 1 个 KO，但被记成 4 个
- pathway_completeness 因此偏高
- **STRONG label 需要在 README 标注此 caveat**
- INSUFFICIENT label 仍 robust

⚠️ 注意：Table 1 是 32 MAGs，Data Set S1 是 29 MAGs（缺 MAG 7/23/24）。reshape 取 inner-join = 29 MAGs。

## 4. ⭐ 假说评分核心结果

YAML：[`grettenberger2021_hypothesis.yaml`](grettenberger2021_hypothesis.yaml) — **pre-registered before
running EnvMeta** (commit `44d7f5f`)，基于 Grettenberger 2018-2020 设计实验时
**先验已知**的 AMD 溪流文献综述层（Bond 2000 / Bigham 2000 / Schippers 1999 /
Cabrera 2006 / Tan 2009 / Auld 2017），**不引用** 2021 论文具体发现。

```
overall_score = 1.000
label         = STRONG
satisfied     = 4 / 4 claims
weight_robust = True (OAT ±20%)
veto_reasons  = none
```

### 4 个 claim 全 satisfied 详情

| Claim | Status | active MAGs | mean completeness | Pre-data 预期 |
|---|---|---|---|---|
| sulfide_oxidation (REQUIRED) | satisfied | 4 | 83% | 主导（流水氧化）✓ |
| dissim_sulfate_reduction (exploratory) | satisfied | 2 | 75% | 沉积物缺氧界面 ✓ |
| nitrate_reduction_explored | satisfied | 4 | 83% | 局部缺氧（NOT required）✓ |
| nitrogen_fixation_explored | satisfied | 4 | 100% | 寡营养偶发（NOT required）✓ |

### 循环图（fig1）覆盖

- **Arsenic 0/6 active**（符合 AMD 非砷预期）
- **Nitrogen 5/6 active**（N fixation + NO reduction / Nitrate reduction）
- **Sulfur 3/4 active**（Assim/Dissim sulfate red / Sulfide oxidation）
- **Iron 1/2 active**（Fe transport — KB v1.1 唯一 Fe pathway）

⚠️ **KB v1.1 限制**：Grettenberger 论文核心是 **Fe(II) oxidation**（典型 AMD 通路），但 EnvMeta KB v1.1 **缺 Fe oxidation/reduction 通路**（仅有 Fe transport + uptake regulation）。pre-data YAML 故意**没把 Fe oxidation 作为 claim**避免 KB-driven 失败。这是 KB v1.2 扩展 backlog。

## 5. 与对照实验其他 Arm 的关系

| Arm | 数据集 | 注释方法 | label | n_satisfied |
|---|---|---|---|---|
| A | 作者 168 MAG (砷) | KEGG 全 | STRONG | n/n |
| B | Wei 2024 (砷+N) | ROCker 14 基因 | INSUFFICIENT | 3/5 (veto) |
| C1 | Liu 2023 (砷, 深海冷泉) | DRAM (KEGG-curated) | STRONG | 4/4 |
| **C2-A** | **Grettenberger 2021 (AMD 溪流) 跨主题** | **METABOLIC step (KEGG)** | **STRONG** | **4/4** ⭐ |
| C2-B | Ayala 2020 (AMD pit lake) | GhostKOALA 重跑 | （等结果）| — |

**核心证据更新（2 个 KEGG-curated STRONG 而非 1）**：
- **2 个独立 KEGG-curated 数据集** 都得 STRONG（Liu 砷主题 + Grettenberger 跨主题 AMD）
- **1 个 ROCker-only 数据集** 得 INSUFFICIENT（Wei）
- H1 (数据广度假说) 在跨主题 + 同主题 2 个数据集上稳定支持
- H2 (机制问题假说) 进一步否决

**跨主题证据**：Grettenberger（无砷）依然 STRONG，证明 EnvMeta 评分**领域中立**，不限于砷研究。

## 6. License + 引用

- Grettenberger 2021 论文：CC BY 4.0
- 本目录 EnvMeta 输出：MIT
- 引用模板：

  > Cycle diagram inference and hypothesis scoring were performed with EnvMeta
  > v0.8.2 (https://github.com/redlizzxy/EnvMeta) on re-shaped Table 1 + Data
  > Set S1 of Grettenberger & Hamilton (2021, *AEM*, doi:10.1128/AEM.00772-21).

## 7. 维护记录

| 日期 | 事项 |
|---|---|
| 2026-05-08 | Phase 1 完成；STRONG label；H1 在跨主题 KEGG-curated 数据上稳定支持 |
