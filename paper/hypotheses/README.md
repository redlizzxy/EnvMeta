# 机制假说 YAML 评分器（S3 / L2 层）

EnvMeta 的循环图推断是**描述性、无假说输入**的。本模块让用户**自带假说 YAML**，
对现有 CycleData 输出打"证据加权分"，返回 strong / suggestive / weak / insufficient
标签。与推断架构解耦。

> 📖 **理论依据 + 完整注释示例 + 参考文献（带 DOI）**：
> 请打开 [`arsenic_steel_slag.yaml`](arsenic_steel_slag.yaml) —— 顶部 80 行
> 是**教程 + 理论框架**，底部是 **6 篇标准 Vancouver 格式文献**（自查 DOI）。
>
> 本 README 仅做 **schema 速查表**。

---

## 快速开始

1. 复制 [arsenic_steel_slag.yaml](arsenic_steel_slag.yaml) 到你的项目，
   **先读顶部注释段**（不要跳过）
2. 按你的假说改 claims 列表
3. 在 EnvMeta 循环图页底部展开「🧪 假说评分 (可选)」
4. 上传 YAML → 点击「评分」→ 看结果表 + 彩色标签 + 下载 TSV/JSON

---

## Schema 速查

```yaml
name: "假说名（必填）"
description: "自由文本（可选）"
author: "..."
version: "1.0"

strong_threshold: 0.75         # overall ≥ 此值 → strong
suggestive_threshold: 0.40     # overall ≥ 此值 → suggestive

claims:                        # 至少 1 条
  - id: 唯一标识
    type: claim 类型（4 选 1）
    weight: 1.0                # 推荐比例 核心=2.0 / 主辅=1.0 / 可选=0.5
    description: "..."
    params:
      # 每种 type 的 params 不同，见下方四类 claim 表
```

---

## 四类 Claim 参数

### 1. `pathway_active` — 某通路活跃

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `pathway` | str | — | 通路名（匹配 display_name 或 pathway_id，模糊匹配）|
| `min_completeness` | float | 50 | 平均完整度 % 阈值 |
| `min_contribution` | float | 0 | total_contribution 阈值 |

**评分**：活跃+两阈值达标 → `satisfied`(1.0)；活跃但阈值差 → `partial`(0.5)；
无活跃 MAG → `unsatisfied`(0)；通路名匹配不到 → `skipped`（不扣分）

### 2. `coupling_possible` — 跨元素化学耦合

| 字段 | 类型 | 说明 |
|---|---|---|
| `species_a` | str | KB 里登记的物种名（如 `As(V)`、`Fe(III)`、`S-2`、`NO3-`）|
| `species_b` | str | 另一物种 |

**评分**：KB 有 + 两端观测到 → `satisfied`；仅一端 → `partial`；
两端都没有 → `unsatisfied`；KB 里没这对 → `skipped`

**目前 KB 支持**（`envmeta/geocycle/knowledge_base/elements.json` → `couplings`）：
- `As(III) + S-2 → As2S3` (precipitation)
- `As(V) + Fe(III) → Fe-As_surface` (adsorption)
- `Fe(II) + S-2 → FeS` (precipitation)
- `NO3- + As(III) → As(V)` (redox)

### 3. `env_correlation` — 通路 × 环境因子相关性

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `pathway` | str | — | 通路名 |
| `env_factor` | str | — | 环境因子列名（`Eh`/`Total_As`/`pH` 等）|
| `expected_sign` | `positive`/`negative`/`any` | `any` | 期望符号 |
| `min_confidence` | `strong`/`suggestive`/`weak` | `suggestive` | 最低置信度要求 |

**评分**：符号对+confidence 达标 → `satisfied`；符号对但 confidence 低一档 → `partial`；
符号反 → `unsatisfied`；找不到记录 → `skipped`

**置信度阶梯**（来自 S2 置换检验）：
`strong > suggestive > weak > spurious? > none`
`spurious?` 永不满足（即使符号对，置换检验判其不可信）。

### 4. `keystone_in_pathway` — 通路含 keystone MAG

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `pathway` | str | — | 通路名 |
| `min_keystones` | int | 1 | 最少 keystone 数 |

**评分**：≥ min → `satisfied`；1~min-1 → `partial`；0 → `unsatisfied`；
通路匹配不到 → `skipped`

---

## 评分聚合

```
overall_score = Σ(weight_i × score_i) / Σ(weight_i)   # 仅对非 skipped
```

| overall | 标签 |
|---|---|
| ≥ strong_threshold | `strong` |
| ≥ suggestive_threshold | `suggestive` |
| > 0 | `weak` |
| = 0 或 all skipped | `insufficient` |

---

## 局限（v1）

- ❌ 不支持 `group_contrast` claim（跨组比较）—— v2 加
- ❌ 不支持分层假说（sub-mechanisms）
- ❌ 不做 multiple testing correction（本身是证据权衡，不是 hypothesis testing；
  理论依据见 YAML 顶部 § 1）
- ❌ 不会帮用户"生成 YAML"—— 假说是用户的科学判断，不是工具的职责

---

## 示例输出

跑 `arsenic_steel_slag.yaml` 对 B 组数据（实测）：

```
overall_score: 1.000  [strong]

claim_id                        status      score  explanation
iron_transport_active           satisfied   1.00   Fe transport: 活跃...
keystone_on_iron                satisfied   1.00   含 1 keystone (Sulfuricaulis sp. Mx_141)
arsenate_reduction_active       satisfied   1.00   Arsenate reduction: 活跃...
fe_as_adsorption_coupling       satisfied   1.00   As(V) + Fe(III) → Fe-As_surface
s_as_precipitation_coupling     satisfied   1.00   As(III) + S-2 → As2S3
ammonia_ox_eh_link              satisfied   1.00   ρ=+0.84 [strong]
arsenate_total_as_link          satisfied   1.00   ρ=+0.76 [suggestive]
```

CK/A/B 三组实跑结果见 [`paper/benchmarks/validation/hypothesis/`](../benchmarks/validation/hypothesis/)。
