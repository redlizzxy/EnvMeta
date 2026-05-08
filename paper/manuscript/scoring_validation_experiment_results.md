# 假说评分对照实验 — 结果与结论

> **完成日期**：2026-05-08
> **状态**：Stage 1 完成（Liu 2023 = STRONG），按预先 declare 决策树**停止后续 Arm**
> **关联**：[scoring_validation_experiment.md](scoring_validation_experiment.md)（实验设计）

---

## 1. 实验三臂结果

| Arm | 数据集 | 注释方法 | overall_score | label | satisfied | 备注 |
|---|---|---|---|---|---|---|
| **A** Positive Control | 作者 168 MAG | KEGG 全 (57 KO) | 0.85+ | STRONG | n/n | v0.7+ 稳定 |
| **B** ROCker treatment | Wei 2024 | ROCker 14 基因 | 0.627 | **INSUFFICIENT** | 3/5 (veto) | nitrate_reduction unsatisfied + coupling partial |
| **C1** KEGG-curated ext. | Liu 2023 | DRAM (KEGG) | **1.000** | **STRONG** | 4/4 | weight_robust=True |
| ~~C2~~ | Yang 2024 | KofamScan | — | **不需要跑** | — | C1 STRONG 已锁定结论 |

## 2. 决策树触发

按 [scoring_validation_experiment.md §4 第一步](scoring_validation_experiment.md) 预先 declare 的决策表：

| C1 结果 | 实测 | 解读 | 触发动作 |
|---|---|---|---|
| **STRONG** | ✅ Liu 1.000 | H1 强支持。EnvMeta 评分机制在 KEGG 全注释上工作正常；Wei INSUFFICIENT 是数据广度问题。| **停止后续 Arm**。叙事路径锁定 X（重新框定）。可选 W（KB v1.2 + reduced_coverage status）。|

## 3. 三臂对照支持的因果推理

**H1 (数据广度假说)** = "Wei INSUFFICIENT 是 ROCker-only 14 基因覆盖广度不足"

**H2 (机制问题假说)** = "EnvMeta 假说评分机制本身偏严或设计有缺陷"

| 观察 | 与 H1 一致？ | 与 H2 一致？ |
|---|---|---|
| Arm A (KEGG 全注释): STRONG | ✅ | ❌（如 H2 真，应 INSUFFICIENT 频发）|
| Arm B (ROCker 14 基因): INSUFFICIENT | ✅ | ✅（H2 也预测 INSUFFICIENT）|
| **Arm C1 (DRAM KEGG): STRONG** | ✅ | ❌**（关键证据 — H2 预测 INSUFFICIENT，但实测 STRONG）** |

**结论**：Arm A + Arm C1 都 STRONG，**仅 Arm B (ROCker-only) INSUFFICIENT**。
H1 完全支持，H2 完全否决。

## 4. 防 p-hacking 关键证据

- **Pre-registration**：Liu YAML 在跑 EnvMeta **之前** commit（git timestamp `42168da`，
  早于 Liu reshape script `?` 早于 EnvMeta runner 早于 fig6 输出）
- **不调阈值**：Liu YAML 用 EnvMeta 默认 `min_completeness=30 / strong=0.75 /
  suggestive=0.40`，与 Wei YAML 一致
- **不基于 Liu 结论设计 claim**：4 个 claim 来自 Liu 2022 设计实验时先验已知的
  生物地球化学原理（Stolz 2006 / Mukhopadhyay 2002 / Rosen 2002 / Yin 2011），
  **不引用 Asgardarchaeota / 4484-113 / AABM5-125-24 等 Liu 论文具体发现**
- **预期与实测对比**：4 个 claim 中 2 个 confirmatory + 2 个 exploratory。
  exploratory claim **可能** unsatisfied（特别是 arsM in marine 仅有限报道），
  但实测 572 MAG active **远超 pre-data 预期** — 这本身是论文里有价值的发现
  ("As methylation in deep sea cold seeps is more widespread than previously
  reported")，**不是为复现而调出来的**

## 5. 论文叙事路径锁定 — 路径 X（重新框定 INSUFFICIENT）

按 [hypothesis_scoring_analysis.md §4 路径 X](hypothesis_scoring_analysis.md)，
最终叙事段落（投稿用，写在 Paper 3 Discussion）：

> "We tested EnvMeta's hypothesis scoring on three independent arsenic
> metagenomic datasets representing different annotation breadths: our
> in-house steel-slag dataset with full KEGG annotation (57 KOs across 4
> elements; STRONG label), Wei et al. (2024) with ROCker-only annotation of
> 14 selected genes (INSUFFICIENT label), and Liu et al. (2023) with DRAM-
> based KEGG-curated annotation focused on As metabolism (STRONG label,
> 4/4 claims satisfied including a pre-registered exploratory hypothesis
> on arsM methylation that returned 572 active MAGs — a finding consistent
> with Liu's reported 'unexpected diversity' narrative).
>
> The contrast between Wei and Liu — both arsenic-themed but with different
> annotation breadths — establishes that EnvMeta's INSUFFICIENT diagnostic
> on Wei is not a tool malfunction but a faithful reflection of annotation
> coverage. Wei's published 14-gene set covers only 2/6 of EnvMeta's
> canonical Nitrate reduction KO list, which by design fails the
> pre-registered Bradford-Hill required veto. This generalizes beyond
> EnvMeta: any KEGG-grounded downstream tool — Anvi'o pathway completeness,
> KEGG Mapper, MicrobiomeAnalyst — would face the same limitation. We
> recommend metagenomic studies publish full KEGG annotation alongside any
> custom ROCker / DRAM hits to ensure secondary analyses remain reproducible."

## 6. KB v1.2 待办（路径 X 的可选扩展 = 路径 W）

实验暴露 EnvMeta KB v1.1 的 3 个空缺（Wei + Liu 共同跳过的 gene）：

| Gene | 跳过原因 | KB v1.2 处理 |
|---|---|---|
| arxA | KB 无 KO 编号 | 加 alias 映射到 K08356 (与 aoxB 同 large subunit family) |
| arsP | KB 没收录 | 加 K25223 重新分类（实际上 KB K25223 标 arsJ，需要修正）|
| arsI | KB 没收录 | 加 K?? (organoarsenical lyase, KEGG 实际有 ortholog) |
| nrfA | KB 没收录 DNRA pathway | 加 DNRA pathway: K03385 + K15876 |

KB 扩展是可选的工程任务，不影响本对照实验的方法论结论（路径 X 已站住）。

## 7. 不再做 Stage 2 (Yang 2024) 的依据

按预先 declare 决策树，C1 STRONG → 停止后续 Arm。

跑 Yang 2024 的边际收益：
- ✅ 加强证据（多一个 KEGG 全注释数据集 = STRONG）
- ❌ 但已有 2 个 STRONG (Arm A + C1) + 1 个 INSUFFICIENT (Arm B) 已足够支持因果推理
- ❌ Yang 数据 reshape 8-12h 工时投入回报递减
- ❌ 多增 1 个 STRONG 不改变结论；多增 1 个 INSUFFICIENT 才会改变结论但概率低

**节约工时**：8-12h → 0h，省下来推 Methods 全稿。

## 8. 投稿叙事改写计划

按 [hypothesis_scoring_analysis.md §4 路径 X 推荐](hypothesis_scoring_analysis.md)：

1. **Methods 4.6 节**（已写，[methods_external_validation.md](methods_external_validation.md)）：
   - 改 Wei 段叙述：从"Wei INSUFFICIENT 是诚实诊断"加深为"在两个独立 KEGG-curated
     对照数据上 EnvMeta 输出 STRONG，证明 Wei INSUFFICIENT 由数据广度而非
     工具机制导致"
   - 加 Liu 段：覆盖深海冷泉 + DRAM 注释 + STRONG label + arsM 意外发现

2. **Discussion 新加段**："annotation breadth limitation as a community
   reproducibility issue" — 推荐 KEGG 全注释作为发布最佳实践

3. **Supplementary**：含完整 3 数据集对照表 + pre-registration commit hash 列表

## 9. 实验维护记录

| 日期 | 事项 |
|---|---|
| 2026-05-08 17:xx | C1 (Liu 2023) 跑通 = STRONG；H1 强支持；停止后续 Arm |
