# EnvMeta 外部数据集验证 — Wei et al. 2024 Microbiome 砷稻田

> **状态**：Phase 1 完成（2026-05-08）
> **目标**：在独立第 2 砷场地数据集上跑通 EnvMeta，证明工具不只能跑作者自己的
> 砷渣-钢渣数据。
> **关联**：[../../../manuscript/external_validation_plan.md](../../../manuscript/external_validation_plan.md)

---

## 1. 数据集来源

| 项 | 内容 |
|---|---|
| 论文 | Wei et al. (2024) "Various microbial taxa couple arsenic transformation to nitrogen and carbon cycling in paddy soils" |
| 期刊 | *Microbiome* 12:236 |
| DOI | [10.1186/s40168-024-01952-4](https://doi.org/10.1186/s40168-024-01952-4) |
| License | CC BY-NC-ND 4.0 |
| 数据范围 | 36 paddy soil metagenomes（华南 12 个采样点 × 3 重复，pH 4.6-8.0）|
| 分组 | AsContam (n=18) vs NoContam (n=18) |
| 数据入口 | NCBI BioProject `PRJNA1068274` (raw FASTQ); 论文 Suppl `MOESM2_ESM.xlsx`（处理后表）|

## 2. 复现路径（用户操作）

> ⚠️ **License 限制**：Wei 论文 Suppl 是 CC BY-NC-ND 4.0，禁止 redistribute
> derivative works，因此 reshape 后的 `input_data_local/` 不入 git 仓库。
> 复现步骤：

```powershell
# 1. 下载 Wei 论文 Suppl（公开免费）
#    https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-024-01952-4
#    → "Supplementary Information" 区 → 下 MOESM2_ESM.xlsx
#    保存到 D:\download\40168_2024_1952_MOESM2_ESM (1).xlsx

# 2. 跑 reshape 脚本（xlsx → 6 个 EnvMeta input file）
conda activate envmeta
python tools/external_benchmarks/wei2024_reshape.py `
    --xlsx "D:\download\40168_2024_1952_MOESM2_ESM (1).xlsx" `
    --out paper/benchmarks/external/wei_2024_paddy/input_data_local

# 3. 跑 EnvMeta 全套（cycle + 4 plots + hypothesis）
python tools/external_benchmarks/wei2024_run_envmeta.py
```

总耗时：reshape ~3 s + EnvMeta 全套 ~24 s = **~30 s**（Win 16 GB 笔记本）。

## 3. EnvMeta 输入文件（reshape 输出）

reshape 脚本从 Wei `MOESM2_ESM.xlsx` 提取并转换：

| EnvMeta 输入 | Wei 来源 | 内容 |
|---|---|---|
| `metadata.tsv` | Table S1 | 36 samples × {SampleID, Group=AsContam/NoContam, Replicate} |
| `env_factors.tsv` | Table S2 | 36 samples × {pH, WH2O%, Cr, Sb, As, Cd, Pb, Cu, NH₄⁺, NO₃⁻, TotalFe, OM, TC, TN, TP, TS} |
| `mag_taxonomy_labels.tsv` | Table S7 | 179 MAGs × GTDB-Tk taxonomy |
| `quality_report.tsv` | Table S7 | 179 MAGs × CheckM-style schema (completeness/contamination) |
| `abundance.tsv` | Table S7 | 179 MAGs × 36 samples（**group-mean 复制**，见局限性） |
| `kegg_target_only.tsv` | Table S8 | 540 MAG-KO records (12 mapped genes × MAG presence; arxA + nrfA 跳过) |

### 局限性透明记录

1. **abundance.tsv 维度退化**：Wei Suppl 仅提供 group-level 平均丰度（AsContam
   / NoContam 两列），不提供 per-sample MAG mapping coverage。reshape 把
   group mean 复制到该组所有样本，所以同组 36 列内的值全部相同。后果：
   - **PCoA / RDA / α-boxplot / LEfSe / 共现网络**：跑不通或维度退化（同组样
     本完全重合到 1 点），**未尝试**
   - **环境相关性 panel**：因 abundance 全 group-level，所有通路 vs 所有 env
     因子产生**伪相关 ρ=0.82**（reflect 不了真实样本间变异）；图里这部分仅
     供布局展示，**数值不应解读**
   - **可跑且数值合理**：循环图 / MAG 质量 / gene heatmap / log2FC / MAG
     abundance heatmap（这些图按 group means 计算，维度限制不影响）

2. **KO 注释精选不全**：Wei 用 ROCker 自定义模型，仅注释 14 个目标基因
   （aioA/arxA/arrA/arsC1/arsC2/napA/narG/nirK/nirS/norB/nosZ/nrfA/pmoA/pmoB），
   非完整 KEGG KO 集。这影响：
   - **arxA**：ROCker 自定义模型，KEGG 无对应 KO，EnvMeta KB v1.1 未收录 → 跳过
     （17 MAG / 18 gene copies skipped）
   - **nrfA**：DNRA pathway，EnvMeta KB v1.1 未收录 → 跳过
     （33 MAG / 33 gene copies skipped）
   - **Nitrate reduction**：KB 含 6 KO（narG/H/I + napA/B + narB），Wei 仅提供 2
     个（napA + narG），单 MAG 完整度 ≤ 33%，部分 MAG 只带 1 个（17%）低于
     30% 阈值，导致 pathway-level 显示"无活跃 MAG"
   - **S / Fe 元素**：Wei 未注释，循环图显示 sulfur (0/4 active) / iron (0/2
     active)，是数据限制，非工具缺陷

## 4. EnvMeta 输出（5 图 + 1 假说评分）

`envmeta_outputs/` 含：

| 文件 | 说明 | 状态 |
|---|---|---|
| `fig1_cycle_diagram.{pdf,png}` | 4 元素 × 18 通路循环图（179 MAG, 6 通路 active）| ✅ EnvMeta 独家 |
| `fig1_cycle_diagram_stats.tsv` | 通路活性 + env coupling 长表（180 rows） | |
| `fig2_mag_quality.{pdf,png}` | MAG completeness × contamination 散点（179 MAGs）| ✅ |
| `fig3_gene_heatmap.{pdf,png}` | 12 mapped KO × {AsContam, NoContam} z-score heatmap | ✅ |
| `fig4_log2fc.{pdf,png}` | 12 mapped KO × log₂(AsContam / NoContam) | ✅ |
| `fig5_mag_heatmap.{pdf,png}` | top MAGs × 36 samples abundance heatmap | ✅ |
| `fig6_hypothesis_score.md` | YAML 假说评分（5 claims）| ✅ EnvMeta 独家 |

## 5. 假说评分关键结果（fig6）

YAML：[`wei2024_hypothesis.yaml`](wei2024_hypothesis.yaml) — 复现作者核心论断
"As(III) oxidation 与 denitrification 通过 NO₃⁻ 作为电子受体共发生"。

```
overall_score = 0.627
label         = INSUFFICIENT  ← 被 required veto 否决
null_p        = 0.9000 (n=999 permutations)
weight_robust = True (OAT ±20%)
satisfied     = 3 / 5 claims
```

**veto reasons**:
1. `nitrate_reduction_active: unsatisfied (required=true)` — Wei 仅提供 2/6 KO
2. `as_n_coupling_arsenite_nitrate: partial (required=true)` — KB 有 As(III)↔NO3-
   redox coupling 记录，但因 Nitrate reduction 通路 unsatisfied，only one end observed

**这是 EnvMeta 设计原则的诚实体现**：

> "描述而非断言。因果判断是用户职责；工具避免附和假说"
> ([../../../../CLAUDE.md](../../../../CLAUDE.md) §产品定位与核心设计决策)

EnvMeta 没有盲目说"作者假说成立 ✓"，而是输出了一个**结构化诊断**：

- ✅ As(III) 氧化通路 active（9 active MAG, 50% completeness）
- ❌ Nitrate reduction 通路 KO 覆盖不足（KB 6-KO 中仅 2 个被 Wei 注释）
- ✅ Nitrite reduction 通路 active（10 active MAG）
- ✅ N₂O reduction 通路 active（52 active MAG, 100% completeness — 仅需 nosZ）
- ⚠️ As-N coupling KB-recognized 但通路 unsatisfied → partial

**审稿人价值**：工具的"不附和假说"姿态（不强行打 strong 标签）证明 EnvMeta
是**研究助手**而非"假说背书机"。Wei 论文的 As-N 耦合结论没问题，但**复现性
受限于作者发布的 KO 注释广度**——这是科学社区可改进的方向（推 KEGG 全注释 +
扩 EnvMeta KB v1.2 加 DNRA pathway）。

## 6. EnvMeta 在此数据集上的"额外价值"

Wei 原文已展示的图（Fig 1-7 + S1-S18）覆盖了：物种系统发育树 / 地理分布 /
基因丰度 box plot / NMDS / Pearson 相关。EnvMeta **新增**：

| 新增能力 | 论文里没有 | EnvMeta 自动给出 |
|---|---|---|
| 4 元素循环图自动推断 | ❌ | ✅ fig1_cycle_diagram |
| YAML 假说评分 + null_p + weight sensitivity | ❌ | ✅ fig6_hypothesis_score |
| Bradford Hill required veto 机制 | ❌ | ✅ "INSUFFICIENT" 透明诊断 |
| 30 秒一键出 5 张 publication-grade PDF | 论文用 R/Python 自定义脚本 | ✅ fig1-5 |

## 7. License + 引用

- Wei et al. 2024 论文：CC BY-NC-ND 4.0 — 必须 cite，禁止 redistribute Suppl 原文件
- 本目录 EnvMeta 输出图：新生成产物，按 EnvMeta 仓库 [LICENSE](../../../../LICENSE) (MIT) 分发
- 引用模板：

  > Cycle diagram inference and hypothesis scoring on this paddy-soil dataset
  > were performed with EnvMeta v0.8.2 (https://github.com/redlizzxy/EnvMeta).
  > Input data were re-shaped from Supplementary Tables S1, S2, S7, and S8 of
  > Wei et al. (2024, *Microbiome* 12:236, doi:10.1186/s40168-024-01952-4).

## 8. 时间记录

| 步骤 | 耗时 | 备注 |
|---|---|---|
| 下载 Suppl xlsx | < 1 min | 浏览器一键下 |
| reshape 脚本运行 | ~3 s | xlsx → 6 input file |
| Cycle diagram | 22.79 s | 含 999-permutation null |
| MAG quality | 0.04 s | |
| Gene heatmap | 0.05 s | |
| log2FC | 0.06 s | |
| MAG heatmap | 0.09 s | |
| Hypothesis score | 0.02 s | 5 claims + null_p + sensitivity |
| **总（不含人工）** | **~24 s** | Win 16 GB 笔记本 |

## 9. 维护记录

| 日期 | 事项 |
|---|---|
| 2026-05-08 | Phase 1 完成（reshape + 5 图 + 1 假说） |
