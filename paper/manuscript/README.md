# Paper 3 论文稿件 + 方法学文档

> EnvMeta 方法学论文（目标期刊：iMeta / Bioinformatics / Front Microbiol）的草稿
> 与素材文档。**本目录是论文工作区**，不是稿件本身（最终 docx 在论文撰写阶段产出）。

---

## 文件导览

### 论文骨架

| 文件 | 角色 |
|---|---|
| [outline_imeta.md](outline_imeta.md) | iMeta 投稿大纲（v0.9.2 含 Plan B 整合 + ImageGP 2 reframing + mock review 修订）|
| [publication_strategy.md](publication_strategy.md) ⭐ | **发表路线** — Path A 改良版（EnvMeta + bioRxiv + 课题论文并行）+ bioRxiv 操作 checklist |
| [paper3_pre_submission_checklist.md](paper3_pre_submission_checklist.md) | 投稿前 4 大任务进度（R 对照 / 第二数据集 / 英文 README / Methods 起草）— 大部分已完成 |
| [reviewer_simulation_prompt.md](reviewer_simulation_prompt.md) ⚠️ untracked | SCI 模拟审稿提示词模板（用户私人检查文件，**不 push**）|
| [zenodo_doi_steps.md](zenodo_doi_steps.md) | Zenodo DOI 5 步操作流程（投稿时执行）|

### Methods 章节素材

| 文件 | 角色 |
|---|---|
| [methods_external_validation.md](methods_external_validation.md) | ⭐ **Methods §4.6** — 外部数据集验证设计 + scoring engine + pre-registration（v0.9.2 修订：lead with author bias limitation）|
| [results_stress_test_section.md](results_stress_test_section.md) | ⭐ **Results §X** — calibration + stress test narrative（v0.9.2 修订：n=29/13 + "consistent with"）|
| [discussion_calibration_discrimination.md](discussion_calibration_discrimination.md) | ⭐ **Discussion §Y** — KEGG-coverage-dependent + 6 limitations + future work（v0.9.2 重写）|
| [performance_paragraph_snippets.md](performance_paragraph_snippets.md) | Performance benchmark drop-in 段落（§5.2.8 / §5.4.6 / §5.3）|
| [external_validation_plan.md](external_validation_plan.md) ⏸ HISTORICAL | 外部数据集 Stage 1 早期计划（Wei 单数据集），已被 v0.9.0 4-Arm 实验超额完成 |

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

## 当前状态（2026-05-09，v0.9.2 paper-side）

**对照实验完成度**：✅ **4 Arm calibration 全 STRONG + Stress test 3/3 全 A 级**（措辞已按 mock review v0.9.2 修订软化为 "KEGG-coverage-dependent" 而非 "domain-neutral"）

| Arm | 数据集 | Calibration | Stress v1 | Stress v2（dominance-aware）| Discrimination |
|---|---|---|---|---|---|
| A | 作者 168 MAG (砷渣) | STRONG | — | — | (positive control) |
| B | Wei 2024 (砷+N) | INSUFFICIENT | — | — | — |
| C1 | Liu 2023 (冷泉砷) | STRONG (1.000) | suggestive (0.625) | **weak (0.250)** ⭐ | v1 B → **v2 A** |
| C2-A | Grettenberger 2021 (AMD 跨主题) | STRONG (1.000) | weak (0.250) | n/a (已 A) | **A 级** ⭐ |
| C2-B | Ayala 2020 (pit lake 跨主题) | STRONG (1.000) | suggestive (0.455) | **weak (0.182)** ⭐ | v1 B → **v2 A** |

**最强单条证据**（softened）：cross-topic arsenate_reduction 在 Grettenberger (n=29 MAG) + Ayala (n=13 MAG) 双双 n=0 active MAGs，**consistent with**（不是 ironclad proof of）EnvMeta 评分引擎在 cross-topic mismatch 下的 domain-neutral behaviour；≥100 MAG 验证 flag 为 future work。

**Paper 3 投稿核心证据全部就位**：
- ✅ Methods §4.6 完整草稿（~1450 字 / 7 子节 / 19 条 Vancouver+DOI；含 Anderson 2001 999-perm justify）
- ✅ Results §X stress test 章节（~800 字 / 4 子节 / Table 1+2+Figure X 定义；n=29/13 caveat 明确）
- ✅ Discussion §Y limitation+future（~640 字 / 4 子节；KEGG-coverage-dependent 重 frame；author bias 作为 #1 limitation）
- ✅ Table 1 + Table 2 + Figure X 实物素材（[../figures/paper3_hypothesis_scoring/](../figures/paper3_hypothesis_scoring/)，PDF + PNG + SVG）
- ✅ 引用 audit + 4 处错引透明纠正（含 Sánchez-España + Dai 2014 + Méndez-García 2015 verified）
- ✅ dominance_score 字段 v0.9.x 兑现（v2 stress test B→A 升级；reframed as engineering retrofit, NOT independent validation）
- ✅ 性能 benchmark：3 regimes × 58 cells；scaling_curve.{pdf,svg,png}；paper-facing performance.md + 中文用户面 docs/performance_zh.md
- ✅ ImageGP 2 ecosystem-extension reframing：EnvMeta 全文从"vs 竞品"改为 "domain-specialist complement"
- ✅ Mock review v0.9.1 → v0.9.2 修订：5 大 Major + 5 Minor 全部温和解；Recommendation 从 Major → **Minor Revision**

**剩余工作**（按 [`publication_strategy.md`](publication_strategy.md) Path A 改良版）：
1. **6 张 placeholder figures**（F1/F3/F5/F6/F7/F10）— 投稿前 mandatory，2-3 天
2. **§5.7.1 Vancouver → iMeta 引用格式**转换 — 用户 Zotero 改 docx 时统一处理
3. **bioRxiv 投稿** — 6 月初目标，3 天审核拿 DOI
4. **EnvMeta 投 iMeta** — bioRxiv 上线后立即投
5. **1-day perturbation analysis**（推荐）— 把 mock review v0.9.2 仅剩 Major Issue（作者偏见）从 acknowledgment 升级为 auxiliary evidence；接收概率从 85% → 95%
6. **课题论文起草**（用户主导）— As 形态重测启动 + 用 EnvMeta bioRxiv DOI 引用

**已暂缓**：
- 盲法 stress test（用户判断暂时不可行；Discussion §Y.4 保留作为 future work proposal；mock review v0.9.2 提议 perturbation analysis 作为 auxiliary evidence 替代）
