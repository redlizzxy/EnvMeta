# EnvMeta 外部数据集验证 — Ayala-Muñoz 2020 Microorganisms IPB Pit Lake 深层

> **状态**：Phase 1 完成（2026-05-09 凌晨）
> **角色**：对照实验 Arm C2-B — KEGG-curated 跨主题（IPB 金属耐性 pit lake 深层）数据集
> **关联**：[../../../manuscript/scoring_validation_experiment.md](../../../manuscript/scoring_validation_experiment.md)

---

## 1. 数据集来源

| 项 | 内容 |
|---|---|
| 论文 | Ayala-Muñoz D, Burgos WD, Sánchez-España J, Couradeau E, Falagán C, Macalady JL (2020) "Metagenomic and Metatranscriptomic Study of Microbial Metal Resistance in an Acidic Pit Lake" |
| 期刊 | *Microorganisms* 8(9):1350 |
| DOI | [10.3390/microorganisms8091350](https://doi.org/10.3390/microorganisms8091350) |
| License | CC BY 4.0 |
| 数据规模 | 13 MAGs × 1 pooled sample (Cueva de la Mora pit lake, IPB Spain, 35m deep) |
| 注释方法 | **GhostKOALA 重跑**（论文未公开 KEGG 注释，重做 ORF + GhostKOALA）|
| 主题 | AMD pit lake 缺氧深层金属耐性（**无砷主导**，主背景 Fe/Cu/Zn/SO4）|

## 2. 复现路径

```powershell
# 1. 下载 13 MAG genomes from BioProject PRJNA646106 (NCBI FTP)
python tools/external_benchmarks/ayala2020_download_and_predict.py
# → 输出: input_data_local/all_mags_proteins.faa (24841 蛋白序列)

# 2. 上传 all_mags_proteins.faa 到 GhostKOALA (KEGG)
#    选项: Eukaryotes + Prokaryotes + Viruses
#    使用机构邮箱 (cugb.edu.cn 等)
#    异步等待 4-12h（实际 ~1h，11234/24841 = 45.2% hit）
#    下载 user_ko.txt → 保存到 input_data_local/

# 3. 运行 reshape (user_ko.txt + MAG list → 6 EnvMeta inputs)
conda activate envmeta
python tools/external_benchmarks/ayala2020_reshape.py `
    --user_ko paper/benchmarks/external/ayala_2020_pitlake/input_data_local/user_ko.txt `
    --out paper/benchmarks/external/ayala_2020_pitlake/input_data_local

# 4. 跑 EnvMeta 假说评分（calibration + stress 双 YAML）
python tools/external_benchmarks/run_stress_yaml.py `
    --dataset ayala_2020_pitlake --yaml ayala2020_hypothesis.yaml
python tools/external_benchmarks/run_stress_yaml.py `
    --dataset ayala_2020_pitlake --yaml ayala2020_hypothesis_stress.yaml
```

总耗时：MAG 下载 ~5 min + ORF prediction ~2 min + GhostKOALA 异步 ~1h + reshape ~1s + EnvMeta ~3s = **~70 min** wall clock。

## 3. EnvMeta 输入文件

| 文件 | 来源 | 内容 |
|---|---|---|
| `metadata.tsv` | placeholder | 1 pooled sample "All"（Ayala deep layer 35m，单 sample）|
| `env_factors.tsv` | placeholder | pH=2.5 proxy（IPB pit lake 一般值）|
| `mag_taxonomy_labels.tsv` | 文件名推断 | 13 MAGs × phylum-level taxonomy（A_CRE/A_EUR/B_ACI/...）|
| `quality_report.tsv` | placeholder | 13 MAGs × 70/5 placeholder（待 Suppl Table 替换）|
| `abundance.tsv` | placeholder | 13 MAGs × 1 sample，全填 1.0 |
| `kegg_target_only.tsv` | GhostKOALA user_ko.txt | **11243 unique MAG-KO records** |

⚠️ **MAG 列表**（5 Archaea + 8 Bacteria）:
- Archaea: A_CRE_07, A_EUR_01, A_EUR_06, A_MIC_10, A_NAN_12
- Bacteria: B_ACI_09, B_ACT_02, B_ACT_11, B_CHL_03, B_DOR_08, B_NIT_04, B_PAT_13, B_PRO_05

⚠️ **Taxonomy + Quality 是 placeholder**（用文件名推 phylum；quality 70/5 默认）—
完整 GTDB-Tk + CheckM2 结果待 Ayala 2020 Suppl Table S1 提供后替换；hypothesis
scoring 结果**不依赖**这两个字段（cycle_diagram 和 score 主要用 ko_long）。

## 4. ⭐ 假说评分核心结果

### 4.1 Calibration YAML — STRONG (4/4)

YAML：[`ayala2020_hypothesis.yaml`](ayala2020_hypothesis.yaml)（**pre-registered** commit `76a4f77`）

```
overall_score = 1.000
label         = STRONG
satisfied     = 4 / 4 claims
weight_robust = True (OAT ±20%)
veto_reasons  = none
```

| Claim | Status | active MAGs | mean completeness | Pre-data 预期 |
|---|---|---|---|---|
| dissim_sulfate_reduction (REQUIRED) | satisfied | 8 | 84% | 缺氧高 SO4 backbone ✓ |
| sulfide_oxidation_microaerophilic | satisfied | 1 | 67% | 微氧界面 ✓ |
| nitrate_reduction_explored | satisfied | 3 | 50% | 寡营养偶发 ✓ |
| nitrogen_fixation_explored | satisfied | 2 | 100% | AMD diazotrophy ✓ |

### 4.2 Stress YAML — SUGGESTIVE (2/4 satisfied)

YAML：[`ayala2020_hypothesis_stress.yaml`](ayala2020_hypothesis_stress.yaml)（**pre-registered** commit `50c4687`）

```
overall_score = 0.455
label         = SUGGESTIVE
satisfied     = 2 / 4 claims
weight_robust = True
```

| Claim | Type | Expected | Observed | 评估 |
|---|---|---|---|---|
| A. sulfide_oxidation_should_dominate (反向) | pathway_active | unsatisfied | **satisfied** ⚠️ | 数据真有 1 MAG 带 sox/sqr (contrib=66.7)，二元阈值放大成 satisfied |
| B. arsenate_reduction_should_dominate (cross-topic) | pathway_active | unsatisfied | **unsatisfied** ✅ | n=0 active MAGs；与 Grettenberger 一致推翻 universal arsC 担忧 |
| C. dissim_sulfate_reduction_should_NOT_dominate (negative) | pathway_inactive | unsatisfied | **unsatisfied** ✅ | n=8, mean_comp=84%, 明显违反 negative |
| D. nitrate_reduction_calibration_anchor | pathway_active | satisfied | **satisfied** ✅ | 跑分系统正常 |

### 4.3 Discrimination 等级

**B 级**（与 Liu 类似）：
- ✅ Cross-topic arsenate_reduction unsatisfied —— 关键正面证据
- ⚠️ Reversed sulfide oxidation 意外 satisfied —— 暴露二元阈值 limit（同 Liu）
- ✅ Negative claim 工作正常

**核心发现**：
1. **领域中立性证据加固**：第 2 个无砷 dataset (Ayala) cross-topic arsenate_reduction 实测 n=0，与 Grettenberger 一致 → **n=2 KEGG-curated cross-topic stress test 双双 reject 砷代谢主导**，universal arsC 担忧被双重反驳
2. **二元阈值 limit 在 Ayala 重现**：Sulfide oxidation 1 MAG contrib=66.7 vs Sulfate reduction 8 MAG contrib=675（10× 弱），但都 satisfied → **加强了 dominance_score 改进的必要性**（v0.9 / Paper 4 future work）

## 5. 与对照实验其他 Arm 的关系

| Arm | 数据集 | 注释方法 | Calibration | Stress | Discrimination |
|---|---|---|---|---|---|
| A | 作者 168 MAG (砷) | KEGG 全 | STRONG | — | (positive control) |
| B | Wei 2024 (砷+N) | ROCker 14 基因 | INSUFFICIENT | — | (无法 stress) |
| C1 | Liu 2023 (砷, 深海冷泉) | DRAM | STRONG (1.000) | suggestive (0.625) | B 级 |
| C2-A | Grettenberger 2021 (AMD 溪流) 跨主题 | METABOLIC | STRONG (1.000) | **weak (0.250)** | **A 级** ⭐ |
| **C2-B** | **Ayala 2020 (pit lake 深层) 跨主题** | **GhostKOALA** | **STRONG (1.000)** | **suggestive (0.455)** | **B 级** |

**对照实验 final**：
- **n=4 KEGG-curated 数据集 calibration 全 STRONG**
- **n=3 stress test 全 score 显著低于 calibration**（0.455-0.625 vs 1.000）
- **2/2 无砷数据集 cross-topic arsenate_reduction 都 unsatisfied** → 领域中立性铁证
- **2/3 stress test 暴露二元阈值 limit**（Liu A + Ayala A） → dominance_score future work

## 6. License + 引用

- Ayala 2020 论文：CC BY 4.0
- 13 MAG 序列：BioProject PRJNA646106 公开
- 本目录 EnvMeta 输出：MIT
- 引用模板：

  > Cycle diagram inference and hypothesis scoring were performed with EnvMeta
  > v0.8.2 (https://github.com/redlizzxy/EnvMeta) on 13 MAGs from
  > BioProject PRJNA646106 (Ayala-Muñoz et al. 2020, *Microorganisms*,
  > doi:10.3390/microorganisms8091350) re-annotated with GhostKOALA (Kanehisa
  > Labs, KEGG).

## 7. 维护记录

| 日期 | 事项 |
|---|---|
| 2026-05-08（早晚）| 数据准备：13 MAG 下载 + ORF 预测，提交 GhostKOALA 异步 |
| 2026-05-09（凌晨）| GhostKOALA 完成 (1h) → reshape + 跑 EnvMeta calibration + stress |
| 2026-05-09 | Calibration STRONG (1.000) + Stress SUGGESTIVE (0.455)；B 级 discrimination；cross-topic arsenate_reduction reject 加固领域中立性 |
