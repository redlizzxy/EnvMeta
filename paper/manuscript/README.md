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
| `mock_review_v0.9.{1,2,3}_*.md` ⚠️ untracked | 历次模拟审稿报告（v0.9.1 Major Rev → v0.9.2 Minor Rev → v0.9.3 Minor Rev acceptance-track；本地检查文件，**不 push**）|
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
| [perturbation_analysis_results.md](perturbation_analysis_results.md) | ⭐ **Perturbation analysis 结果** — within + cross 双模式；Liu cross 0/20 STRONG headline；mock review v0.9.2 Major #1 auxiliary evidence |

### 配套外部资源（不在本目录但相关）

- 6 个 hypothesis YAML：[`paper/benchmarks/external/*/`](../benchmarks/external/) （Liu / Grettenberger / Ayala 各 2 个 calibration + stress）
- 假说写作指南：[`docs/hypothesis_writing_guide.md`](../../docs/hypothesis_writing_guide.md)
- 假说设计原则：[`paper/hypotheses/HYPOTHESIS_DESIGN_PRINCIPLES.md`](../hypotheses/HYPOTHESIS_DESIGN_PRINCIPLES.md)
- Stress runner：[`tools/external_benchmarks/run_stress_yaml.py`](../../tools/external_benchmarks/run_stress_yaml.py)

---

## 当前状态（2026-05-10，v0.9.6 paper-side — handling-editor 4 cheap Major 修完，stop mock review，待 figures + bioRxiv）

**对照实验 + auxiliary evidence 完成度**：✅ **4 Arm calibration 全 STRONG + Stress test 3/3 全 A 级 + Perturbation analysis 双模式 N=20/each (4 datasets，Arm A partial)**

| Arm | 数据集 | Calibration | Stress v2 | Discrimination | Perturbation (within / cross STRONG retention) |
|---|---|---|---|---|---|
| A | 作者 168 MAG (砷渣) | STRONG (1.000, 9/9) | — | (positive control) | 100% / **100%** ⭐ partial-perturbation robustness (3 of 9 claims) |
| B | Wei 2024 (砷+N) | INSUFFICIENT | — | — | — |
| C1 | Liu 2023 (冷泉砷) | STRONG (1.000) | weak (0.250) ⭐ | A 级 | 50% / **0%** ⭐ headline |
| C2-A | Grettenberger 2021 (AMD) | STRONG (1.000) | A 级 ⭐ | A 级 | 50% / 30% |
| C2-B | Ayala 2020 (pit lake) | STRONG (1.000) | weak (0.182) ⭐ | A 级 | 40% / 15% |

**Annotation-breadth gradient**（mock review v0.9.4 头号证据）：
Arm A (partial perturbation, 3/9 claims) 100% → Grettenberger (mixed) 30% → Ayala (mixed, smaller) 15% → Liu (focused) 0% — 3 个 external 数据集的 ordering 是 load-bearing pattern（Arm A 100% 反映 partial perturbation 设计 + 6 unperturbed 锚点，不是数据 saturation）。

**两条 headline 证据**：
1. cross-topic arsenate_reduction 在 Grettenberger (n=29) + Ayala (n=13) 双双 n=0 active MAGs（**consistent with** cross-topic mismatch 下 domain-neutral，softened 表述）
2. **Liu 2023 cross-element perturbation 0/20 STRONG**（median 0.000）→ 元素级 target accuracy 是机制刚需；within-element 40-50% retention 与 §Y.1 KEGG-coverage-dependent 一致

**Paper 3 投稿核心证据**：
- ✅ Methods §4.6 完整草稿（8 子节 / ~1700 字含 §4.6.7 perturbation；19 条 Vancouver+DOI；含 Anderson 2001 999-perm justify）
- ✅ Results §X 章节（5 子节 / ~1080 字含 §X.3 perturbation；Table 1+2+Figure X+X-bis 定义；n=29/13 caveat 明确）
- ✅ Discussion §Y limitation+future（~700 字 / 4 子节；KEGG-coverage-dependent 重 frame；§Y.3 limitation #1 整合 perturbation auxiliary evidence）
- ✅ Table 1 + Table 2 + Figure X 实物素材（[../figures/paper3_hypothesis_scoring/](../figures/paper3_hypothesis_scoring/)）+ Figure X-bis perturbation_curve（[../benchmarks/external/perturbation/](../benchmarks/external/perturbation/)）
- ✅ 引用 audit + 4 处错引透明纠正
- ✅ dominance_score 字段 v0.9.x 兑现（reframed as engineering retrofit）
- ✅ 性能 benchmark：3 regimes × 58 cells；scaling_curve；paper/benchmarks/performance.md + docs/performance_zh.md
- ✅ ImageGP 2 ecosystem-extension reframing
- ✅ Mock review v0.9.1 → v0.9.2 → v0.9.3：Major Rev → Minor Rev → **Minor Rev (acceptance-track)**

**剩余工作**（已 stop mock review；v0.9.5 4 cheap Major 修完）：
1. **6 张 placeholder figures**（F1/F3/F5/F6/F7/F10）— 投稿前 mandatory，2-3 天
2. **bioRxiv 投稿** — 6 月初目标，3 天审核拿 DOI
3. **EnvMeta 投 iMeta** — bioRxiv 上线后立即投
4. **课题论文起草**（用户主导）— As 形态重测启动 + 用 EnvMeta bioRxiv DOI 引用

**carry 到 revision 阶段处理（不影响投稿）**：
- 2-D band gap threshold sweep（v0.9.5 Major #2）— 1-2h
- 引用列表扩到 ≥40（v0.9.5 Major #5）— 1-2h
- 9 个 v0.9.5 Minor + 10 个 v0.9.4 Minor — 投稿后 reviewer 反馈再迭代

**已暂缓**：
- 盲法 stress test（Discussion §Y.4 仍作为 future work；perturbation analysis 已作为 auxiliary evidence 替代 acknowledgment-only state）
