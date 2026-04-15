# 机制假说 YAML 评分器（S3 / L2 层）

## 这是什么

EnvMeta 的循环图推断是**描述性、无假说输入**的（给什么数据就推出什么通路 /
相关，不对"对不对、有没有 drive"下结论）。但实际写论文的流程是：

> "我们**假说**钢渣促进铁氧化，Fe(III) 表面吸附砷酸盐。数据能支持吗？"

本模块用一份 **YAML 假说清单** + 现有 CycleData 输出 → 证据链评分。

- ✅ 解耦：假说评估与推断不互相污染
- ✅ 描述性：评分反映"数据支持度"，不等于"因果证实"
- ✅ 用户可控：阈值、权重、claim 组合都在 YAML 里，一看就知道评分怎么算来的

## 快速开始

1. 复制 [arsenic_steel_slag.yaml](arsenic_steel_slag.yaml) 到你自己的项目
2. 按你的假说改 claim（见下方 schema）
3. 在 EnvMeta 循环图页底部展开「🧪 假说评分 (可选)」
4. 上传 YAML → 点击「评分」→ 看结果表 + 标签 + 下载 TSV

## YAML Schema

```yaml
name: "假说名（显示用，必填）"
description: "自由文本（可选）"
author: "..."
version: "1.0"

strong_threshold: 0.75         # 默认 0.75，≥ 此值 → strong 标签
suggestive_threshold: 0.40     # 默认 0.40，≥ 此值 → suggestive

claims:                        # 至少 1 条
  - id: 唯一标识
    type: claim 类型（4 选 1）
    weight: 1.0                # 权重（任意正数）
    description: "可选说明"
    params:
      # 每种 type 的 params 不同，见下表
      ...
```

## 四类 Claim

### 1. `pathway_active` — 某通路活跃

**params**：
| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `pathway` | str | — | 通路名（匹配 display_name 或 pathway_id，模糊匹配） |
| `min_completeness` | float | 50 | 平均完整度 % 阈值 |
| `min_contribution` | float | 0 | total_contribution 阈值 |

**评分**：
- 活跃 + 两阈值达标 → `satisfied` (1.0)
- 活跃但阈值差一档 → `partial` (0.5)
- 无活跃 MAG → `unsatisfied` (0)
- 通路名匹配不到 → `skipped`（从分母里除掉，不扣分）

### 2. `coupling_possible` — 跨元素耦合

**params**：
| 字段 | 类型 | 说明 |
|---|---|---|
| `species_a` | str | KB 里登记的物种名（如 `As(V)`、`Fe(III)`、`S-2`、`NO3-`） |
| `species_b` | str | 另一物种 |

**评分**：
- KB 里有 `species_a ↔ species_b` 耦合 AND 两端物种都被数据观测到 → `satisfied`
- KB 有但只一端观测到 → `partial`
- KB 有但两端都未观测 → `unsatisfied`
- KB 里根本没这对 → `skipped`（用户写错物种名也不扣分，而是提示 skip）

**目前 KB 支持的耦合**（详见 `envmeta/geocycle/knowledge_base/elements.json` →
`couplings`）：
- `As(III) + S-2 → As2S3`（precipitation）
- `As(V) + Fe(III) → Fe-As_surface`（adsorption）
- `Fe(II) + S-2 → FeS`（precipitation）
- `NO3- + As(III) → As(V)`（redox）

### 3. `env_correlation` — 通路 × 环境因子相关

**params**：
| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `pathway` | str | — | 通路名 |
| `env_factor` | str | — | 环境因子列名（如 `Eh`、`Total_As`、`pH`） |
| `expected_sign` | `positive`/`negative`/`any` | `any` | 期望符号 |
| `min_confidence` | `strong`/`suggestive`/`weak` | `suggestive` | 最低置信度要求 |

**评分**：
- 符号对 + confidence 达标 → `satisfied`
- 符号对但 confidence 低一档 → `partial`
- 符号不对 → `unsatisfied`
- 找不到 (pathway, env_factor) 记录 → `skipped`

**置信度阶梯**（来自 S2 置换检验）：
`strong > suggestive > weak > spurious? > none`
`spurious?` 永远不算满足（即使符号对），因为置换检验已判其"看起来相关但不可信"。

### 4. `keystone_in_pathway` — 通路含 keystone MAG

**params**：
| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `pathway` | str | — | 通路名 |
| `min_keystones` | int | 1 | 最少 keystone 数 |

**评分**：
- ≥ min_keystones → `satisfied`
- 1 ~ min_keystones-1 个 → `partial`
- 0 个 → `unsatisfied`
- 通路匹配不到 → `skipped`

## 评分聚合

```
overall_score = Σ(weight_i × score_i) / Σ(weight_i)        # 仅对非 skipped 的 claim
```

**标签映射**：
| overall | 标签 | 含义 |
|---|---|---|
| ≥ strong_threshold | `strong` | 主要 claim 都被数据支持 |
| ≥ suggestive_threshold | `suggestive` | 部分支持，但有明显缺口 |
| > 0 | `weak` | 零星支持 |
| = 0 或 total_weight = 0 | `insufficient` | 数据不足以评估 |

## 撰写 checklist

写假说 YAML 时建议：

- [ ] 每条 claim 都写 `description`，让读者（审稿人）看得懂为什么你列这条
- [ ] 权重分配合理：**核心机制 2.0**，辅助证据 1.0，可选 0.5
- [ ] 把"预期正 / 负相关"明确写在 `expected_sign`，而不是 `any`（更诚实）
- [ ] `min_confidence` 先用 `suggestive`；若期刊严格可紧到 `strong`
- [ ] 如果你的假说涉及 group 对比（A vs B 增量），当前 v1 还不支持
  `group_contrast`，workaround：分别对 group=A / group=B 跑评分，
  比较两次的 overall_score

## 与传统统计检验的关系

假说评分**不是**单一 p-value；它是**证据加权**。好处：

1. 一条 claim 被打脸不会让整个假说塌掉（还有其他 claim 支撑）
2. 反过来：多条 weak 证据可叠加出 strong 标签，比"找 1 条显著相关"更稳健
3. 可审计：每条 claim 的 evidence 都在 TSV 里，审稿人能逐条核对

## 局限（v1）

- ❌ 不支持跨 group 对比 claim（`group_contrast`）—— v2 加
- ❌ 不支持分层假说（sub-mechanisms）—— 扁平 claim list
- ❌ 不会帮用户"生成 YAML"—— 这个是用户的科学判断，不是工具的职责
- ❌ 不做 multiple testing correction —— 本身是证据权衡不是 hypothesis testing

## 示例输出

跑 `arsenic_steel_slag.yaml` 对 B 组数据的预期（仅示意）：

```
overall_score: 0.80  [strong]

claim_id                        status      score  explanation
iron_transport_active           satisfied   1.00   Fe transport: 6 MAG 活跃...
keystone_on_iron                satisfied   1.00   Fe transport: 含 1 个 keystone (Sulfuricaulis sp. Mx_141)
arsenate_reduction_active       satisfied   1.00   Arsenate reduction: 活跃...
fe_as_adsorption_coupling       satisfied   1.00   As(V) + Fe(III) → Fe-As_surface: 两端物种都观测到
s_as_precipitation_coupling     satisfied   1.00   As(III) + S-2 → As2S3: ...
ammonia_ox_eh_link              satisfied   1.00   Ammonia ox ↔ Eh: ρ=+0.84 [strong]
arsenate_total_as_link          satisfied   1.00   Arsenate red ↔ Total_As: ρ=+0.76 [suggestive]
as_transport_active             partial     0.50   As transport 活跃但完整度低...
```
