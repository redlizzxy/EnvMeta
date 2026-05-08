# 假说评分机制 — Wei 2024 复盘 + 论文叙事决策

> **创建日期**：2026-05-08
> **目的**：在 Wei 2024 复现得到 INSUFFICIENT 后，系统梳理评分机制运行原理、
> 不一致原因、对论文的指导价值、4 条叙事路径，供决策。
> **关联**：
> - 复盘对象：[../benchmarks/external/wei_2024_paddy/envmeta_outputs/fig6_hypothesis_score.md](../benchmarks/external/wei_2024_paddy/envmeta_outputs/fig6_hypothesis_score.md)
> - 评分器实现：[../../envmeta/geocycle/hypothesis.py](../../envmeta/geocycle/hypothesis.py)
> - YAML 模板：[../../paper/hypotheses/arsenic_steel_slag.yaml](../../paper/hypotheses/arsenic_steel_slag.yaml)

---

## TL;DR（30 秒读完）

EnvMeta 在 Wei 数据集上输出 INSUFFICIENT，根源是**三层叠加**：(1) EnvMeta KB v1.1 缺 DNRA + arxA；(2) Wei 仅发布 14 个 ROCker 精选基因，不是 KEGG 全注释；(3) EnvMeta YAML 默认阈值偏严 + 3 个 required veto。

**机制不是"失败"** —— 它正确捕捉到了"作者发布的注释广度限制了下游可复现性"这一更深层的社区数据实践问题。

**叙事风险**：单纯调阈值复现 STRONG 标签 = p-hacking 嫌疑。**最稳的路径是把 INSUFFICIENT 升级为对社区数据实践的方法学反思**，写在 Discussion 而非 Methods。

---

## 1. 运行原理（结构化总结）

EnvMeta 假说评分是一个 **MCDA evidence-weighting scorecard**（多准则决策分析的证据加权评分卡），不是 statistical hypothesis testing。

### 1.1 输入

用户提供一份 YAML，包含若干 claim（论断），每个 claim 含：

```yaml
- id: arsenite_oxidation_active     # claim 名
  type: pathway_active               # 5 类之一
  weight: 1.5                        # MCDA 权重
  required: true                     # Bradford Hill 必要条件
  params:                            # 类型特异参数
    pathway: "Arsenite oxidation"
    min_completeness: 30
```

### 1.2 5 类 claim 评估器

| Claim type | 评估什么 | 输出 score |
|---|---|---|
| `pathway_active` | 通路里有多少 MAG 完整度 ≥ 阈值 | active MAG 数 / 期望阈值 → [0, 1] |
| `coupling_possible` | KB 是否记录该化学物对的 coupling + 数据是否两端都观测到 | 0 / 0.5 (partial) / 1 |
| `env_correlation` | 通路 vs env 因子的 ρ 是否符合 expected_sign + min_confidence | 0 / 1 |
| `keystone_in_pathway` | 通路是否含 ≥ N 个 keystone MAG | 0 / 1 |
| `group_contrast` | 通路在指定组间是否有差异 | 0 / 1 |

每个 claim 输出 `(status, score, weight, evidence)`，其中 status ∈ {satisfied, partial, unsatisfied, skipped}。

### 1.3 三层聚合

```
overall_score = Σ(weight × score) / Σ(weight)         # MCDA 加权平均
                                                       # 仅 status≠skipped 的 claim 参与
label_base = strong   if overall ≥ strong_threshold
             suggestive if overall ≥ suggestive_threshold
             weak       if overall > 0
             insufficient if overall = 0

# Bradford Hill required veto（硬否决）
if any required claim has status ≠ satisfied:
    label = INSUFFICIENT  ← 不管 overall 多高
    veto_reasons = [那些失败的 required claim id]
```

### 1.4 三个 robustness 指标（S3.5）

- **null_p**（permutation test）：把 score 在 claim 间**随机重洗 999 次**，看
  random_overall ≥ observed_overall 的频率。null_p < 0.05 = 不像随机巧合。
- **weight_robust**（OAT sensitivity）：每个权重 ±20% 扰动，看 label 是否变化。
  weight_robust=True = label 在权重扰动下稳定。
- **9 档解读标签**：strong / suggestive / weak / insufficient × {robust, sensitive, opposite_under_null} 等组合。

### 1.5 设计原则（CLAUDE.md 明文写）

> **"描述而非断言。因果判断是用户职责；工具避免附和假说。"**

具体表现：
- 工具不输出 p-value，不做多重检验校正
- 工具不说"假说成立 / 不成立"，只输出 4 档 label
- required veto 机制直接对应 Bradford Hill "必要条件"语义

---

## 2. 为什么 Wei 数据集会出现 INSUFFICIENT —— 三层叠加

### Wei 实测结果

```
overall_score = 0.627  →  按阈值应为 suggestive (0.40 ≤ x < 0.75)
label         = INSUFFICIENT  ←  被 required veto 否决
veto_reasons  = [nitrate_reduction_active: unsatisfied (required=true),
                 as_n_coupling_arsenite_nitrate: partial (required=true)]
satisfied     = 3 / 5 claims
null_p        = 0.9000 (n=999)         ← overall 容易被随机重洗达到
weight_robust = True (OAT ±20%)         ← label 稳定
```

### 不一致的三层成因

| 层 | 谁的问题 | 具体表现 | 是否可补救 |
|---|---|---|---|
| **L1** EnvMeta KB v1.1 不全 | EnvMeta 缺陷 | DNRA pathway (nrfA = K03385) + arxA 自定义 ROCker model 未收录 → 33 + 17 MAG 50 records 直接跳过 | ✅ KB v1.2 加 DNRA pathway，~30 min |
| **L2** Wei 注释广度有限 | Wei 数据局限 | Wei 仅发布 ROCker 14 个目标基因，不是 KEGG 全注释。Nitrate reduction 通路 KB 含 6 KO（narG/H/I + napA/B + narB），Wei 只覆盖 napA + narG = 2/6 = 33% | ❌ 改不了，是别人发表数据 |
| **L3** 默认 YAML 阈值偏严 | 用户设计选择 | min_completeness = 30%（单 MAG 含 1/6 = 17% 不够，含 2/6 = 33% 边缘）+ 3 个 required veto（任一失败硬否决整体） | ✅ 阈值可调（但有 p-hacking 风险，见 §4） |

### 关键观察

**INSUFFICIENT 不等于"Wei 结论错"。** 它说的是"在当前 KB + 注释广度 + 阈值组合下，EnvMeta 不能给 Wei 假说 ≥ suggestive 标签"。Wei 的 As-N 耦合结论本身基于他们自己 ROCker hits + 系统发育分析 + 文献综合，**EnvMeta 没否定它**，只是在 KB-anchored 框架下凑不够证据强度。

### 反过来想：作者自己 168 MAG 数据为什么没出这问题？

作者自己数据用 **KEGG 全注释**（57 KO 覆盖 4 元素 18 通路），所以 Nitrate reduction 6 KO 全有，pathway completeness 容易达 50%+，假说 YAML 默认阈值跑出 strong/suggestive 是合理的。**Wei 的 INSUFFICIENT 是数据广度差异的产物，不是工具失效**。

---

## 3. 机制还有没有对论文的指导意义？—— 有，且更高一档

INSUFFICIENT 表面看是"工具失败"，但站到方法学论文的高度看，**揭示了 4 个有价值的发现**：

### 价值 1：揭示了一个被业界忽视的可复现性问题

任何 KEGG-grounded 下游工具（Anvi'o pathway completeness / KEGG Mapper /
MicrobiomeAnalyst / EnvMeta）都依赖完整 KO 注释。但学术界很多论文（包括 Wei 在内）只发布 author-selected gene set，导致：

- 二次分析（meta-analysis / cross-study comparison）难做
- 任何"自动化 pathway analysis"工具都会在这类数据上输出 reduced confidence

EnvMeta 的 INSUFFICIENT 标签**正确捕捉**了这个问题，并把它显示给读者，**这是工具的功能不是 bug**。

### 价值 2：证明 EnvMeta 的"非附和"设计原则

如果 EnvMeta 在 Wei 数据上盲目输出"假说成立"，说明工具有"过度附和"倾向 — 这是 MCDA 类工具的常见问题（Saltelli & Tarantola 2002 警告过）。EnvMeta 没出现这个问题，**正好是设计原则的实战验证**。

### 价值 3：给 EnvMeta KB v1.2 提供了具体扩展方向

Wei 数据暴露了 EnvMeta KB v1.1 的两个空缺：
- DNRA pathway（nrfA + nrfH）需补
- ROCker custom model alias 需要支持机制

这是一份**自带 backlog** 的 case study。

### 价值 4：可以引出"社区数据发布最佳实践"建议

论文 Discussion 可以写一段：

> "We recommend metagenomic studies publish full KEGG annotation alongside
> any author-selected ROCker / DRAM hits, to enable reproducible secondary
> analysis by KEGG-grounded tools."

iMeta 这种方法论期刊喜欢这种"工具论文 + 社区实践批评"的双层贡献。

---

## 4. 如何调整叙述 —— 4 条路径对比

### 路径 Z（高风险）：单纯调阈值让结果变好 ❌

**做法**：看到 INSUFFICIENT → 把 `min_completeness` 从 30 降到 15，去掉 nitrate
required veto → 跑出 SUGGESTIVE/STRONG。

**问题**：调整理由依赖"我已知 Wei 结论是对的"。审稿人质疑：

> "你怎么知道 15% 阈值是合适的？是因为试出来 15 能 reproduce Wei，还是因为 15 在科学上有道理？"

如果答不出后半个问题 → 循环论证 → 拒稿。

**这条不能走。**

### 路径 Y（中度可接受）：发明 preregistered scaling rule，独立于 Wei

**做法**：在 Methods 里**先 declare 一条规则**（不只针对 Wei），再应用：

> "**Adjustment rule**: 当某数据集发布的 gene set 仅覆盖 EnvMeta KB 通路 KO
> 列表的子集，pathway-completeness 阈值按覆盖比例缩放：
> `effective_threshold = default × (covered_KOs / canonical_KOs)`
> 这条规则保证不同注释广度的研究在 EnvMeta 评分时**可比较**。"

应用到 Wei：Nitrate reduction 2/6 → effective threshold = 30 × 0.33 = 10% →
单 MAG 含 1 KO (17%) > 10% → satisfied。

**可信度风险**：
1. 仍可能被质疑"是不是为 Wei 量身定做"
2. 必须在**至少 2 个独立数据集**上应用同一条 rule（不能只 Wei）
3. 必须在论文 Methods 里**充分说明 rule 的科学动机**（不依赖 Wei 结论）

**降低风险的关键**：
- 同一条 rule 也应用回作者 168 MAG（应该影响小，因为是 KEGG 全注释）
- 同一条 rule 也应用到 Tara Oceans（如果做异方向数据集）
- 看 rule 是否在 3 个数据集上**差异化输出**（不是无脑全 satisfied）

### 路径 X（最诚实）：保留 INSUFFICIENT，重新框定叙事 ⭐ **推荐**

**做法**：不调任何阈值。把 INSUFFICIENT 写在 Discussion，叙事升级：

| 原叙事（弱） | 新叙事（强）|
|---|---|
| "EnvMeta 没复现 Wei 结论" | "EnvMeta 揭示了 Wei 等论文发布的 ROCker-only annotation 限制了下游 KEGG-grounded 工具的可复现性 —— 这是社区数据实践问题" |

具体段落（Discussion 用）：

> "Our re-analysis of Wei et al. (2024) returned an INSUFFICIENT label not
> because the original biological conclusion is unsupported, but because the
> published 14-gene ROCker hits cover only a subset of canonical KEGG
> pathways (e.g., 2 of 6 KOs in Nitrate reduction). This finding generalizes
> beyond EnvMeta: any KEGG-grounded downstream tool — Anvi'o pathway
> completeness, KEGG Mapper, MicrobiomeAnalyst — would face the same
> limitation when confronted with author-selected gene sets. We therefore
> recommend that future metagenomic studies publish full KEGG annotation
> alongside any custom ROCker / DRAM hits, to ensure secondary analyses
> remain reproducible. EnvMeta's design — surfacing this coverage gap as a
> structured diagnostic rather than silently averaging over missing KOs —
> embodies the 'transparent insufficiency' principle adopted from
> epidemiological evidence frameworks (Bradford Hill 1965; Suter & Cormier
> 2011)."

**优点**：
- ✅ **零 p-hacking 风险**
- ✅ 方法学贡献更高一档（不只是工具，还是对社区实践的批评）
- ✅ iMeta 偏向方法论批评，比"工具复现"更喜欢
- ✅ 省 2 小时
- ✅ 不引入额外 YAML 维护负担

**潜在风险**：
- 审稿人可能仍会问"为什么不调一下阈值看看？" → 回答："调阈值会引入 p-hacking 风险，本研究坚持默认 strict mode 以保证 INSUFFICIENT 诊断的方法学纯净"

### 路径 W（最强）：X + 工具升级（KB v1.2 加 reduced_coverage status）

**做法**：X 的所有内容 + 在 EnvMeta 评分器里加一个新 status `reduced_coverage`：

> "**reduced_coverage** (区别于 unsatisfied): 当某 pathway 在数据集中的 KO 覆
> 盖率 < 50% 时，pathway_active claim 不输出 satisfied 也不输出 unsatisfied，
> 而输出 reduced_coverage(0.5)，明确告诉用户'数据广度不足，不能判定'。"

这是把 X 的"数据广度问题"概念**真的实现成工具的常驻能力**，而不只在论文里口头提一句。

**实施成本**：1-2 小时（hypothesis.py 加新 status + 测试）

**额外回报**：
- KB v1.2 + 新 status 可以单独写入 Methods 4.4 节（"YAML schema enhancements"）
- EnvMeta v0.9.0 release 可以 highlight 这个新功能
- 论文里可以说"as illustrated by the Wei et al. case study, EnvMeta v0.9 added a `reduced_coverage` status to make annotation-broadth limitations explicit"

---

## 5. 我的推荐 + 决策表

| 路径 | 方法论合法性 | 工时 | 论文叙事强度 | 风险 |
|---|---|---|---|---|
| **X** 重新框定 | ⭐⭐⭐ 完全干净 | 0.5 h（写 Discussion 段）| ⭐⭐⭐⭐ 高 | 几乎为零 |
| **W** X + 工具升级 | ⭐⭐⭐ 完全干净 | 2-3 h | ⭐⭐⭐⭐⭐ 最高 | 几乎为零 |
| **Y** preregistered rule | ⭐⭐ 可接受但要严防 | 4-6 h（含跨 dataset 验证）| ⭐⭐⭐ 中 | 中（要加做 Tara 才能站住）|
| **Z** 单纯调阈值 | ❌ 不可接受 | 1 h | ⭐⭐ 低（被识破就崩）| 高（拒稿）|

**推荐顺序**：W > X > Y > Z（绝不做）

---

## 6. 等待用户决策

3 个具体决策点：

1. **路径选 X / W / Y 哪个？**（我推荐 W；最稳保守是 X）
2. **是否需要做 Tara Oceans 异方向数据集？**（如果选 Y 必须做；选 X/W 可选）
3. **commit 当前 Wei Phase 1 的方式**（已 commit `7eb29d0`；后续路径决定后再补 Discussion 段 / KB v1.2 / 跨 dataset 验证）

请总览后告诉我。
