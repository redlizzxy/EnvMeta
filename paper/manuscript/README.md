# Paper 3 论文稿件 + 方法学文档

> EnvMeta 方法学论文（目标期刊：iMeta / Bioinformatics / Front Microbiol）的草稿
> 与素材文档。**本目录是论文工作区**，不是稿件本身（最终 docx 在论文撰写阶段产出）。

---

## 文件导览

### 论文骨架

| 文件 | 角色 |
|---|---|
| [outline_imeta.md](outline_imeta.md) | iMeta 投稿大纲 |
| [paper3_pre_submission_checklist.md](paper3_pre_submission_checklist.md) | 投稿前 4 大任务进度（R 对照 / 第二数据集 / 英文 README / Methods 起草）|
| [zenodo_doi_steps.md](zenodo_doi_steps.md) | Zenodo DOI 5 步操作流程（投稿时执行）|

### Methods 章节素材

| 文件 | 角色 |
|---|---|
| [methods_external_validation.md](methods_external_validation.md) | Methods 章节 — 外部数据集验证设计 |
| [external_validation_plan.md](external_validation_plan.md) | 外部数据集 Phase 1.0 计划（Wei/Liu 选型）|

### 假说评分对照实验（核心卖点）

文档分层结构（推荐按顺序读）：

| 文件 | 角色 |
|---|---|
| [scoring_validation_experiment.md](scoring_validation_experiment.md) | **实验设计** — 4 Arm 顺序自适应设计 + 决策树 |
| [hypothesis_scoring_analysis.md](hypothesis_scoring_analysis.md) | **4 路径 narrative 分析** + 路径 X 推荐 |
| [scoring_validation_experiment_results.md](scoring_validation_experiment_results.md) | ⭐ **主结果文档** — 4 Arm + Stress test 完整结果 + 维护记录 |
| [scoring_validation_self_critique.md](scoring_validation_self_critique.md) | 严肃自检 — calibration vs discrimination + 局限承认 |
| [stress_test_predictions.md](stress_test_predictions.md) | **冻结的盲预测**（pre-registration anchor `50c4687`）|
| [stress_test_results.md](stress_test_results.md) | ⭐ **Stress test 实测对照** + 论文叙事段落（直接拷贝）|
| [hypothesis_references_audit.md](hypothesis_references_audit.md) | **引用审计** — DOI 表 + 提炼度评估 + 4 处错引纠正 |

### 配套外部资源（不在本目录但相关）

- 6 个 hypothesis YAML：[`paper/benchmarks/external/*/`](../benchmarks/external/) （Liu / Grettenberger / Ayala 各 2 个 calibration + stress）
- 假说写作指南：[`docs/hypothesis_writing_guide.md`](../../docs/hypothesis_writing_guide.md)
- 假说设计原则：[`paper/hypotheses/HYPOTHESIS_DESIGN_PRINCIPLES.md`](../hypotheses/HYPOTHESIS_DESIGN_PRINCIPLES.md)
- Stress runner：[`tools/external_benchmarks/run_stress_yaml.py`](../../tools/external_benchmarks/run_stress_yaml.py)

---

## 当前状态（2026-05-09）

**对照实验完成度**：✅ **4 Arm 全部完成 + Stress test 3/3 完成**

| Arm | 数据集 | Calibration | Stress | Discrimination |
|---|---|---|---|---|
| A | 作者 168 MAG (砷渣) | STRONG | — | (positive control) |
| B | Wei 2024 (砷+N) | INSUFFICIENT | — | — |
| C1 | Liu 2023 (冷泉砷) | STRONG (1.000) | suggestive (0.625) | B 级 |
| C2-A | Grettenberger 2021 (AMD 跨主题) | STRONG (1.000) | weak (0.250) | **A 级** ⭐ |
| C2-B | Ayala 2020 (pit lake 跨主题) | STRONG (1.000) | suggestive (0.455) | B 级 |

**最强单条证据**：cross-topic arsenate_reduction 在 2/2 无砷数据集 (Grettenberger + Ayala) 都 n=0 active MAGs → EnvMeta 领域中立性铁证。

**待完成**（按优先级）：
1. 论文 Methods 4.6 假说评分章节 + Results stress test 章节起草
2. ~~Verify Korehi 2014 / Mendez-Garcia 2015~~ ✅ 完成（2026-05-09）：Korehi 2014 主题错（同 Auld 错引模式），改引 **Dai 2014 PLoS One** (10.1371/journal.pone.0087976) + **Méndez-García 2015 Front Microbiol** (10.3389/fmicb.2015.00475)
3. 实现 dominance_score 字段解决二元阈值 limit（v0.9 future work）
4. 英文 README + LICENSE + Zenodo DOI（投稿硬指标）
