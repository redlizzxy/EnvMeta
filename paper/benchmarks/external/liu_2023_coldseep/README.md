# EnvMeta 外部数据集验证 — Liu et al. 2023 npj Biofilms 深海冷泉

> **状态**：Phase 1 完成（2026-05-08）
> **用途**：对照实验 Arm C1 — KEGG-curated 砷数据集，验证 Wei INSUFFICIENT
> 是否由数据广度而非工具机制导致
> **关联**：[../../../manuscript/scoring_validation_experiment.md](../../../manuscript/scoring_validation_experiment.md)

---

## 1. 数据集来源

| 项 | 内容 |
|---|---|
| 论文 | Liu R et al. (2023) "Unexpected genetic and microbial diversity for arsenic cycling in deep sea cold seep sediments" |
| 期刊 | *npj Biofilms and Microbiomes* 9:13 |
| DOI | [10.1038/s41522-023-00382-8](https://doi.org/10.1038/s41522-023-00382-8) |
| License | CC BY 4.0 |
| 数据规模 | 87 metagenomes（13 全球冷泉点）+ 1741 MAGs (MOESM4) + 1084 As-cycling MAGs (MOESM5) + 2140 As gene hits (MOESM6) |
| 注释方法 | DRAM v1.3.5（含 KEGG 全注释）+ METABOLIC pipeline + GTDB-Tk r207 |
| 主题 | 深海冷泉砷循环（缺氧 + H₂S 富集 + 砷酸盐还原主导）|
| 数据入口 | NCBI BioProject `PRJNA831433` (raw + 22 arsenotrophic MAGs); 论文 Suppl `MOESM2/4/5/6_ESM.xlsx` |

## 2. 复现路径（用户操作）

```powershell
# 1. 下载 4 个 Suppl xlsx（MOESM2/4/5/6）从论文页：
#    https://www.nature.com/articles/s41522-023-00382-8#Sec19
#    保存到 D:\download\

# 2. 跑 reshape (xlsx → 6 EnvMeta inputs)
conda activate envmeta
python tools/external_benchmarks/liu2023_reshape.py `
    --moesm2 "D:\download\41522_2023_382_MOESM2_ESM.xlsx" `
    --moesm4 "D:\download\41522_2023_382_MOESM4_ESM.xlsx" `
    --moesm5 "D:\download\41522_2023_382_MOESM5_ESM.xlsx" `
    --moesm6 "D:\download\41522_2023_382_MOESM6_ESM.xlsx" `
    --out paper/benchmarks/external/liu_2023_coldseep/input_data_local

# 3. 跑 EnvMeta 全套
python tools/external_benchmarks/liu2023_run_envmeta.py
```

总耗时：reshape ~5 s + EnvMeta 全套 ~3 s = **~8 s**（远快于 Wei 24 s，因 Liu 数据已整理好）。

## 3. EnvMeta 输入文件

| EnvMeta 输入 | Liu 来源 | 内容 |
|---|---|---|
| `metadata.tsv` | MOESM5 sample 列名 | 87 samples × {SampleID, Group="All", Replicate} |
| `env_factors.tsv` | placeholder | 深海冷泉 MOESM2 没发布 pH/化学元素，env 留 placeholder |
| `mag_taxonomy_labels.tsv` | MOESM5 Taxonomy 列 | 1084 As-cycling MAGs × GTDB-Tk lineage |
| `quality_report.tsv` | MOESM4 | 1083 MAGs × completeness/contamination |
| `abundance.tsv` | MOESM5 | **1084 MAGs × 87 samples（per-sample，无维度退化）** ⭐ 比 Wei 强一档 |
| `kegg_target_only.tsv` | MOESM6 | 1678 MAG-KO records (8 mapped genes; arxA/arsP/arsI 跳过, KB v1.1 unmapped) |

### KO 映射（Liu 12 As 基因 → EnvMeta KB v1.1）

| Liu | KO | KB pathway |
|---|---|---|
| aioA | K08356 | Arsenite oxidation |
| arrA | K28466 | Resp. arsenate red. |
| arsC1 | K00537 | Arsenate reduction |
| arsC2 | K03741 | Arsenate reduction |
| arsM | K07755 | As methylation |
| acr3 | K03325 | As transport/detox |
| arsB / arsB_935 | K03893 | As transport/detox |
| arsH | K11811 | As transport/detox |
| **arxA** | None | **KB v1.1 未收录**（与 Wei 一致） |
| **arsP** | None | **KB v1.1 未收录** |
| **arsI** | None | **KB v1.1 未收录** |

跳过 462 records（arxA 3 / arsP 449 / arsI 10）。**这与 Wei reshape 跳过 arxA + nrfA 的策略一致**。

## 4. EnvMeta 输出

| 文件 | 状态 |
|---|---|
| `fig1_cycle_diagram.{pdf,png}` | ✅ As 5/6 pathway active；N/S/Fe 全空（数据未发布） |
| `fig1_cycle_diagram_stats.tsv` | 36 rows（仅 As 通路 + sensitivity）|
| `fig2_mag_quality.{pdf,png}` | ✅ 1084 MAGs scatter |
| `fig3_gene_heatmap.{pdf,png}` | ✅ 8 KO × group means |
| `fig4_log2fc` | ⏸ 跳过（Liu 单组 pooled，无 case/control）|
| `fig5_mag_heatmap.{pdf,png}` | ✅ top MAGs × 87 samples |
| `fig6_hypothesis_score.md` | ✅ ⭐ **STRONG** label |

## 5. ⭐ 假说评分核心结果

YAML：[`liu2023_hypothesis.yaml`](liu2023_hypothesis.yaml) — **pre-registered before
running EnvMeta** (commit `42168da`)，基于 Liu 2022 设计实验时**先验已知**的
生物地球化学原理（不引用 Liu 论文具体发现，避免循环论证）。

```
overall_score = 1.000
label         = STRONG
satisfied     = 4 / 4 claims
weight_robust = True (OAT ±20%)
veto_reasons  = none
```

### 4 个 claim 全 satisfied 详情

| Claim | Status | active MAGs | Pre-data 预期 | 是否符合 |
|---|---|---|---|---|
| arsenate_reduction (required) | satisfied | 18 | 应主导（缺氧+高As）| ✅ |
| as_transport/detox (required) | satisfied | 1（边界）| 应活跃 | ✅ |
| as_methylation (exploratory) | satisfied | **572** | 海洋 arsM 报道有限 → 不确定 | ✅ **远超预期** |
| respiratory_as_reduction (exploratory) | satisfied | 17 | 某些海洋沉积物报告过 | ✅ |

⭐ **意外发现**：arsM (砷甲基化) 在 572 个 MAG 中 active —— 远超 pre-data 预期
(基于 Yin 2011 ES&T "海洋 arsM 报道有限")。这恰好对应 Liu 论文里关于
"unexpected diversity" 的核心叙述（Asgardarchaeota + 多新门的 arsM 携带）。

### 与 Wei 的对照

| 数据集 | 注释方法 | label | Veto? | 结论 |
|---|---|---|---|---|
| 作者 168 MAG | KEGG 全 | STRONG | none | Positive control ✓ |
| **Wei 2024** | ROCker-only 14 基因 | **INSUFFICIENT** | nitrate_reduction + as_n_coupling | 数据广度不足 |
| **Liu 2023** | DRAM (KEGG) | **STRONG** | none | KEGG-curated 数据正确评分 |

**核心证据**：在 KEGG-curated 数据上（Liu），EnvMeta 评分 = STRONG。
在 ROCker-only 数据上（Wei），EnvMeta 评分 = INSUFFICIENT。
**机制本身工作正常 — Wei INSUFFICIENT 是数据广度限制，不是工具机制问题。**

## 6. 时间记录

| 步骤 | 耗时 |
|---|---|
| 下载 4 个 Suppl xlsx | < 1 min |
| reshape | ~5 s |
| Cycle diagram | 0.41 s |
| 其他 4 图 | < 3 s |
| Hypothesis scoring | 0.02 s |
| **总（不含人工）** | **~9 s** |

## 7. License + 引用

- Liu 2023 论文：CC BY 4.0（比 Wei CC BY-NC-ND 4.0 更宽松，理论上可 redistribute
  derivative，但本仓库**仍按统一对照实验原则**不入 input_data_local/）
- 本目录 EnvMeta 输出：MIT
- 引用模板：

  > Cycle diagram inference and hypothesis scoring were performed with EnvMeta
  > v0.8.2 (https://github.com/redlizzxy/EnvMeta) on re-shaped supplementary
  > data (MOESM2, MOESM4, MOESM5, MOESM6) of Liu et al. (2023, *npj Biofilms
  > Microbiomes* 9:13, doi:10.1038/s41522-023-00382-8).

## 8. 维护记录

| 日期 | 事项 |
|---|---|
| 2026-05-08 | Phase 1 完成；STRONG label；对照实验 H1 假说强支持 |
