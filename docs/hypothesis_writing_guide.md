# EnvMeta 假说 YAML 写作指南

> **文档角色**：EnvMeta 假说评分器（hypothesis scorer）的**方法论 + 模板**指南。
> Schema 速查表请见 [paper/hypotheses/README.md](../paper/hypotheses/README.md)。
> **本指南给"如何写出方法论稳健的假说 YAML"提供框架**——尤其防止 reviewer
> 诘问"倒果为因 / hindsight bias / 假说选择性偏差"。

---

## 1. 核心原则：双层假说设计

写假说 YAML 不是**只**写"我希望被证实的东西"。一份**方法学稳健**的假说要包含
两层 claim：

| 层级 | claim 子集 | 目的 | 期望结果 |
|---|---|---|---|
| **Calibration** | 基于该环境的 backbone biogeochemistry 先验，几乎必然存在的通路 | 校准评分引擎是否在该数据上工作正常 | satisfied |
| **Stress** | 故意违背先验 / 领域错配 / 用 pathway_inactive 否定主导 | 测试评分引擎的 discrimination power | unsatisfied |

**只写 Calibration 层**：评分会偏向 STRONG，有"挑容易题"嫌疑。
**只写 Stress 层**：评分会偏向 INSUFFICIENT，但不证明引擎工作正常。
**两层都写**：能区分"通路不存在 vs 评分机制偏严"，是真正的 discrimination 证据。

> 详见 [scoring_validation_self_critique.md](../paper/manuscript/scoring_validation_self_critique.md)
> 中"Claim 设计代表性"批评 —— 只写 backbone 通路 claim 是 calibration 不是 discrimination。

---

## 2. Pre-registration：防 p-hacking 的工程纪律

**核心承诺**：YAML 必须在跑 EnvMeta **之前** commit 到 git。git timestamp 是
不可逆证据。

### 2.1 必须做的（hard requirement）

1. ✅ **YAML commit 在跑结果 commit 之前**——reviewer 可独立核验 git log
2. ✅ **不调阈值**——用 EnvMeta 默认 `min_completeness=30 / strong=0.75 / suggestive=0.40`
3. ✅ **每条 claim 的 `explanation` 字段**只引用**早于目标论文发表年份 5 年以上**的综述/原创论文
4. ✅ **commit message 标记 pre-registration**：建议 `feat(benchmarks): <dataset> pre-data 假设 YAML — pre-registration commit`

### 2.2 不能做的（hard violation）

| 违规 | 后果 | 例子 |
|---|---|---|
| explanation 引用目标论文 specific finding | 倒果为因 | "Asgardarchaeota 是 keystone（Liu 2023）" → 不能写 |
| 看了 EnvMeta 结果再调 weight 或 threshold | p-hacking | 跑出 0.74 → 把 strong_threshold 调成 0.70 |
| 看了结果再删/加 claim | post-hoc cherry-picking | 某 claim unsatisfied → 删掉 |
| 用论文已知物种身份做 species_a/species_b | 循环论证 | 把 Liu 论文的 keystone 物种写进 coupling claim |

### 2.3 容易忽视的 leakage 风险

⚠️ **作者读过目标论文这事本身**就是 selection bias 来源——即使 explanation 不引用，潜意识仍可能挑"作者会发现的通路"。**唯一有效的 100% 防护**是 **盲法**：
请未读过目标论文的师弟/师妹独立写 YAML。

如果做不到盲法，**至少**包含 stress claim（见 §3.3）作为 partial mitigation。

---

## 3. Pre-prediction：防 hindsight bias 的工程纪律

**新规则（v0.9 起）**：每条 claim 必须在 `description` 字段里**预先**声明
**expected_label** 和 **reasoning**。这样跑完结果 diff 预测时，hindsight bias
被结构性削弱。

### 3.1 模板

```yaml
- id: arsenate_reduction_should_dominate
  type: pathway_active
  weight: 1.5
  required: true
  description: |
    [PREDICTION] expected_label=satisfied
    [REASONING] 冷泉缺氧+高 As(V) → arsC1/arsC2 主导是 backbone 预测
    （Stolz 2006）。如未 satisfied，说明数据 KEGG 注释缺失或样本量不足。
  params:
    pathway: "Arsenate reduction"
    min_completeness: 30
```

跑完后在 `*_predictions.md` 文件里 diff：

| claim_id | expected | observed | match |
|---|---|---|---|
| arsenate_reduction_should_dominate | satisfied | satisfied | ✅ |
| nitrification_should_NOT_dominate (stress) | satisfied | unsatisfied ⚠️ | ❌ — 意外结果 |

❌ 行就是科学发现的来源——意外结果触发数据复查或先验修正，**不要**通过改 YAML 把它消除。

### 3.2 expected_label 的合法值

- `satisfied` — 期望该 claim 在数据中被证实
- `unsatisfied` — 期望该 claim 在数据中被否定（**stress claim 必用**）
- `partial` — 期望弱信号（边界情况）
- `skipped` — 期望数据/KB 缺失（不能反驳/支持）
- `unknown` — 完全探索性，无明确预期（应使用稀少，每假说不超 20%）

---

## 4. 6 类 claim 选择指南

> 完整 schema 见 [paper/hypotheses/README.md](../paper/hypotheses/README.md)。

| Claim type | 用途 | calibration 适用 | stress 适用 |
|---|---|---|---|
| `pathway_active` | 期望通路活跃 | ✅ 默认 | ❌ |
| **`pathway_inactive`** ⭐ | **期望通路 NOT 活跃**（v0.9 加） | ❌ | ✅ **stress 主力** |
| `coupling_possible` | 跨元素耦合存在 | ✅ | ❌ |
| `env_correlation` | 通路 × 环境因子相关性符号 | ✅ | ✅（用 expected_sign 反向）|
| `keystone_in_pathway` | 通路含 ≥N keystone | ✅ | ⚠️ 受限（KB 依赖）|
| `group_contrast` | 跨组 ratio ≥ min | ✅ | ⚠️ 受限（要 compare_df）|

### 4.1 pathway_inactive — Popperian falsifiability 主力

**新加的 negative claim type**（v0.9 / Paper 4 stress test 引入）：

```yaml
# 例 1: AMD 溪流 stress claim — 砷代谢机制不应主导（领域错配测试）
- id: arsenate_reduction_should_NOT_dominate_in_amd
  type: pathway_inactive
  weight: 1.0
  required: false
  description: |
    [PREDICTION] expected_label=satisfied (即 pathway 应 inactive)
    [REASONING] AMD 溪流是无砷环境（Cabin Branch 系统主背景金属是 Fe/Cu/Zn），
    arsenate reduction 应不存在或仅微弱信号。如 satisfied (= inactive)，
    支持 EnvMeta 领域中立性；如 unsatisfied (= 通路意外活跃)，说明
    数据里有通用 arsC 解毒型基因（不是 arsenic-specific），需重新解读。
  params:
    pathway: "Arsenate reduction"
    max_completeness: 50      # 平均完整度 < 50% 仍算 partial 违反
```

**评估状态**：
- `satisfied` (1.0)：n_active_mags == 0（完全不活跃，符合 negative 预期）
- `partial` (0.5)：有活跃 MAG 但 mean_completeness < max_completeness（弱违反）
- `unsatisfied` (0.0)：有活跃 MAG 且 mean_completeness ≥ max_completeness（明显违反）
- `skipped`：数据里找不到通路（KB 缺失，不能反驳/支持）

### 4.2 stress claim 的三种失败模式

写 stress claim 时建议覆盖**三类失败模式**之中的至少 1-2 种：

| 类型 | 含义 | 例子 |
|---|---|---|
| **A. 反向预测** | 把先验颠倒 | 冷泉（缺氧）→ "as_oxidation 应主导"（实际 reduction 主导）|
| **B. 领域错配** | 把 X 主题机制套到 Y 主题环境 | AMD（无砷）→ "arsenate_reduction 应主导" |
| **C. 程度反转** | 通用 backbone 改成 negative 形式 | 用 pathway_inactive 否定本应主导的通路 |

**推荐 A + B**。C 虽然简单但容易"过严"（背景信号无法剔除）。

---

## 5. 完整模板：calibration + stress 双层假说

```yaml
# Pre-registered: 2026-XX-YY commit XXXXX
name: "<dataset>_dual_layer_hypothesis"
description: |
  <dataset> pre-data 假说，含 calibration (backbone 通路) + stress (反向/错配)
  双层 claims，用于评估 EnvMeta scorer 的 calibration accuracy + discrimination
  power。所有 claim 在跑 EnvMeta 之前 commit；不调阈值；不引用目标论文 specific finding。

version: "1.0"
strong_threshold: 0.75
suggestive_threshold: 0.40

claims:
  # ═════════════════════════════════════════════════════════
  #  CALIBRATION LAYER — backbone 通路，期望 satisfied
  # ═════════════════════════════════════════════════════════
  - id: backbone_pathway_should_dominate
    type: pathway_active
    weight: 1.5
    required: true
    description: |
      [PREDICTION] expected_label=satisfied
      [REASONING] <环境类型> 的 backbone 通路。先验依据 <早于目标论文 5 年的综述>。
    params:
      pathway: "<pathway_name>"
      min_completeness: 30

  - id: secondary_pathway_explored
    type: pathway_active
    weight: 0.5
    required: false
    description: |
      [PREDICTION] expected_label=partial
      [REASONING] 探索性 claim，<comment>。
    params:
      pathway: "<pathway_name>"
      min_completeness: 30

  # ═════════════════════════════════════════════════════════
  #  STRESS LAYER — risky claims，期望 unsatisfied
  # ═════════════════════════════════════════════════════════
  - id: reversed_prediction_should_fail
    type: pathway_active
    weight: 1.0
    required: false
    description: |
      [PREDICTION] expected_label=unsatisfied (stress / 反向预测)
      [REASONING] 把先验颠倒。<环境>是 X，故意预测 Y 应主导。期望 EnvMeta reject。
    params:
      pathway: "<oppositely_predicted_pathway>"
      min_completeness: 30

  - id: cross_topic_should_NOT_dominate
    type: pathway_inactive
    weight: 1.0
    required: false
    description: |
      [PREDICTION] expected_label=satisfied (即 inactive，stress / 领域错配)
      [REASONING] 把 <X 主题机制> 套到 <Y 主题环境>。期望 inactive。
    params:
      pathway: "<cross_topic_pathway>"
      max_completeness: 50
```

---

## 6. Bradford-Hill 对应（论文 Methods 可引用）

EnvMeta scorer 的设计松散对应 Bradford-Hill (1965) 因果推理 9 准则：

| Bradford-Hill 准则 | EnvMeta 实现 |
|---|---|
| Strength | `overall_score` |
| Consistency | 跨数据集 calibration（多个 dataset 都 STRONG）|
| Specificity | `pathway_inactive` claim 检测"领域错配" |
| Temporality | YAML pre-registration（git timestamp）|
| Biological gradient | `min_completeness` 阈值；`max_completeness` 反向 |
| Plausibility | `description` 引用早期文献 |
| Coherence | 多 claim_type 互证（pathway + coupling + env_correlation）|
| Experiment | OAT weight sensitivity (`weight_robust`) |
| Analogy | 跨数据集复现（4-arm 对照实验）|

EnvMeta **不下因果结论**——这是用户职责。scorer 的输出是"数据支持度"，不是
"因果证实"。详见 [hypothesis_scoring_analysis.md](../paper/manuscript/hypothesis_scoring_analysis.md)。

---

## 7. 常见错误自查清单

写完 YAML，在 commit 之前对照过一遍：

- [ ] explanation/description 字段是否引用了**早于目标论文 5 年以上**的综述？
- [ ] 是否完全没有提论文 specific finding（具体物种身份、site-specific 差异、新发现 phylum）？
- [ ] 阈值是否保留 EnvMeta 默认值 (30/0.75/0.40)？
- [ ] 是否有至少 1 个 stress claim？（pathway_inactive 或 type A 反向预测）
- [ ] 每条 claim 的 description 是否含 [PREDICTION] expected_label= 和 [REASONING] ?
- [ ] required 字段：veto 应只用于该环境的 backbone 通路（如冷泉的 arsenate reduction）
- [ ] 是否有未读过目标论文的同事独立审过一遍？（盲法 best practice）

---

## 8. 跑分后的工作流

1. **不要修改 YAML** —— pre-registration 一旦 commit 即冻结
2. **写 `*_predictions_vs_observed.md`** —— 表格 diff [PREDICTION] expected_label vs 实际 label
3. **意外结果如实报告** —— ❌ 行是科学发现的来源；分析为什么
4. **Calibration / Discrimination 拆分报告**：
   - Calibration claims 全 satisfied → 评分引擎工作正常
   - Stress claims 多 unsatisfied → discrimination power 强
   - Stress claims 也 satisfied → 引擎可能偏宽松，或 stress 设计不够 risky（**不要**重写 stress claim 解决，要在 Discussion 里分析）

---

## 9. 模板与示例

| 文件 | 角色 |
|---|---|
| [paper/hypotheses/arsenic_steel_slag.yaml](../paper/hypotheses/arsenic_steel_slag.yaml) | 单层 calibration 假说示例（v1，仅 calibration）|
| [paper/benchmarks/external/liu_2023_coldseep/liu2023_hypothesis.yaml](../paper/benchmarks/external/liu_2023_coldseep/liu2023_hypothesis.yaml) | 单层 calibration（pre-registered）|
| `paper/benchmarks/external/*/<*>_hypothesis_stress.yaml` | 双层 calibration+stress（v0.9 引入）|

---

## 10. 论文表述建议

报告 hypothesis scoring 实验时（论文 Methods / Discussion），优先使用以下措辞：

- ✅ "We **calibrated** EnvMeta's scoring on dataset X"（不是 "demonstrated" 或 "validated"）
- ✅ "stress test claims targeting reversed/cross-topic predictions"
- ✅ "pre-registered hypothesis YAML committed to git before EnvMeta was run"
- ❌ 不要说 "EnvMeta correctly identifies ..."（这暗示 ground truth，是 overclaim）
- ❌ 不要说 "validated the scoring engine"（应说 "characterized the scoring engine's behavior"）

---

## 11. 维护

| 日期 | 事项 |
|---|---|
| 2026-05-08 | 初版 — 含双层 calibration+stress 设计、pathway_inactive、pre-prediction 模板 |
