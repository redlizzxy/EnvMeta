# EnvMeta 假说设计原则（论文方法学）

> **本文档**：用于 Paper 3 / Paper 4 Methods 章节引用 EnvMeta hypothesis scoring 的设计哲学
> 与方法学边界。**用户教程 + schema 速查**见 [docs/hypothesis_writing_guide.md](../../docs/hypothesis_writing_guide.md)
> 和 [paper/hypotheses/README.md](README.md)。

---

## 1. 设计哲学

EnvMeta hypothesis scorer 的核心设计哲学源自三条原则：

### 1.1 描述性而非因果性

scorer 输出 `overall_score` 和 4 档 label（strong / suggestive / weak / insufficient）
反映**数据对预设 claim 的支持度**。它不下"因果证实"结论。因果判断是研究者的科学
判断，不是工具的职责（cf. Bradford-Hill 1965; Hill 1965 *Proc R Soc Med*）。

### 1.2 Pre-registration 强制

YAML 必须在跑 EnvMeta 之前 commit 到 git；`description` 字段引用早于目标论文 5
年以上的综述/原创论文；不调阈值。这与 medical / psych 领域的 OSF pre-registration
精神一致（cf. Nosek et al. 2018 *PNAS*）。

### 1.3 Calibration vs Discrimination 的分层

单数据集 scoring 实验存在两类目标：
- **Calibration**：评分引擎是否在 backbone 通路 claim 上工作正常
- **Discrimination**：评分引擎是否能 reject 反向预测 / 领域错配 claim

EnvMeta v0.9 起强制要求双层假说（calibration + stress claims）来同时报告两类
能力。仅写 calibration claims 的旧 YAML 仅产出 calibration 证据，不构成 discrimination
power 证明。

---

## 2. Bradford-Hill 准则的 partial 对应

EnvMeta scorer 的设计松散对应 Bradford-Hill (1965) 因果推理 9 准则的工程实现：

| Bradford-Hill 准则 | EnvMeta 实现 | 局限 |
|---|---|---|
| **Strength** | `overall_score`（0-1 加权和）| 单数据集；不做 effect size 标定 |
| **Consistency** | 跨数据集 calibration（多 dataset 都 STRONG）| 需要用户主动多数据集复现 |
| **Specificity** | `pathway_inactive` claim 检测领域错配 | v0.9 加入；旧 YAML 无 |
| **Temporality** | YAML pre-registration（git timestamp）| 仅时间锁，不证因果 |
| **Biological gradient** | `min_completeness` / `max_completeness` 阈值 | 二元 satisfied/partial/unsatisfied，不报连续 dose-response |
| **Plausibility** | `description` 字段引用早期文献 | 引用质量靠用户自查 |
| **Coherence** | 多 claim_type 互证（pathway+coupling+env_correlation）| 不做联合 effect estimation |
| **Experiment** | OAT weight sensitivity (`weight_robust`) | 仅扰动权重，不扰动数据本身 |
| **Analogy** | 跨数据集复现（4-arm 对照实验）| 需要研究者主动设计对照 |

我们**不声称** EnvMeta 实现了完整 Bradford-Hill 框架。它实现了其中**可工程化、
可自动评分**的子集。剩余准则（如 dose-response 量化、effect size confidence interval）
需要研究者另行做。

---

## 3. 6 类 claim 的方法学定位

| Claim type | 方法学定位 | 适用层 |
|---|---|---|
| `pathway_active` | 期望某通路活跃 — 默认 calibration claim | calibration |
| `pathway_inactive` | 期望某通路 NOT 活跃 — Popperian falsifiability 主力 | **stress** |
| `coupling_possible` | 跨元素化学物耦合 — 高阶 mechanism | calibration |
| `env_correlation` | 通路 × 环境因子相关性符号 | calibration / stress（用 expected_sign）|
| `keystone_in_pathway` | 通路含 ≥N keystone — 物种层面证据 | calibration（物种身份风险）|
| `group_contrast` | 跨组比值 — 处理效应 | calibration（实验设计依赖）|

---

## 4. Stress test 三种失败模式

| 类型 | 含义 | 代表例 | 期望 |
|---|---|---|---|
| **A. 反向预测** | 把先验颠倒 | 冷泉缺氧 → "as_oxidation 应主导" | unsatisfied |
| **B. 领域错配** | 把 X 主题机制套到 Y 主题环境 | AMD 无砷 → "arsenate_reduction 应主导" | unsatisfied 或 inactive satisfied |
| **C. Negative claim** | 用 pathway_inactive 否定本应主导的通路 | 砷渣 → "arsenate_reduction 应 inactive" | unsatisfied |

推荐覆盖 A + B 两类。C 类容易"过严"（任何环境都有微弱背景信号，inactive 难
satisfied），用作**辅助而非主要 stress claim**。

---

## 5. 方法学局限

EnvMeta hypothesis scoring 实验有 4 类**无法用工具消除**的局限：

### 5.1 作者读过目标论文（Selection Bias）

即使 explanation 不引用论文 specific finding，作者**潜意识**仍可能挑"作者会发现的
通路"做 calibration claim。**唯一根治**是 **盲法 pre-registration**：请未读过目标
论文的同事独立写 YAML。EnvMeta 工具无法强制这一点；研究者应在 Methods 中
report 是否做了盲法。

### 5.2 KB 覆盖率不全（KB Coverage Limit）

KB v1.1 含 4 元素 × 18 通路 × 57 KO（envmeta/geocycle/knowledge_base/elements.json）。
KB 缺失的通路（如 Fe(II)/Fe(III) redox、arsP/arsI、nrfA DNRA）会导致相关 claim
被 skipped 或评分偏低。这不是 scorer 缺陷，是 KB 工程问题。研究者应在 Methods
明确 KB version 和已知缺失。

### 5.3 评分阈值的人为标定（Threshold Anchoring）

`strong=0.75 / suggestive=0.40` 是经验值，没有 frequentist null 校准。EnvMeta 的
`null_p`（permutation test）部分缓解但不消除该问题。Reviewer 可能问"为什么 0.75
是 strong"，回答只能是"经验值，与社区 SUS / NPS 分数文化类似"。

### 5.4 单数据集 vs Meta-analysis

每次 score() 是单数据集证据。跨数据集 calibration 是 narrative 层面的论证，工具
不做正式 meta-analysis（如 random-effect model）。研究者要做正式跨数据集统计，
应导出 per-dataset overall_score / null_p 进 meta-analysis 软件（R `meta` 包等）。

---

## 6. 论文表述建议

报告 EnvMeta hypothesis scoring 实验时，优先使用以下措辞：

✅ **推荐**：
- "We **calibrated** EnvMeta's scoring on dataset X"
- "stress test claims targeting reversed/cross-topic predictions"
- "pre-registered hypothesis YAML committed to git before EnvMeta was run"
- "calibration evidence (not discrimination evidence)"
- "EnvMeta returned label X under default thresholds"

❌ **避免**：
- "EnvMeta correctly identifies ..."（暗示 ground truth）
- "validated / proved / demonstrated the scoring engine"（overclaim）
- "EnvMeta's accuracy is ..."（accuracy 暗示 binary classification）

---

## 7. 引用建议

如果 Methods 引用 EnvMeta hypothesis scoring，建议引用形式：

> Hypothesis scoring was performed with EnvMeta v0.X (https://github.com/redlizzxy/EnvMeta).
> The scoring engine implements weighted evidence aggregation across user-defined
> claims (pathway_active / pathway_inactive / coupling_possible / env_correlation /
> keystone_in_pathway / group_contrast types) with permutation-based null
> calibration and one-at-a-time weight sensitivity scanning. Hypothesis YAMLs
> were pre-registered (committed to git) before running the scoring engine.
> Default thresholds (strong=0.75, suggestive=0.40, min_completeness=30) were
> not tuned per dataset.

---

## 8. 维护

| 日期 | 事项 |
|---|---|
| 2026-05-08 | 初版 — 双层假说设计哲学 + Bradford-Hill 对应 + 6 类 claim + 4 类局限 |
