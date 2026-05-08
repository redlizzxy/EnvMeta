# 假说评分对照实验 — 严肃自检（"倒果为因"风险评估）

> **创建日期**：2026-05-08
> **触发**：用户要求自检 3 份 pre-data YAML 是否存在倒果为因痕迹、是否符合常规科研假说设置、与论文结论一致性如何
> **关联**：[scoring_validation_experiment_results.md](scoring_validation_experiment_results.md)
> **目的**：投稿前先自找 reviewer 会问的硬问题，不要给自己镀金

---

## 1. 三个维度的诚实评估

### 维度 A — 形式上的倒果为因风险（git / 阈值 / 文献）

| 防护措施 | Liu (C1) | Grettenberger (C2-A) | Ayala (C2-B) | 评估 |
|---|---|---|---|---|
| YAML commit 早于跑结果 commit | ✅ `42168da` < `1e4571f` | ✅ `44d7f5f` < `e8d5133` | ✅ `76a4f77` < (待跑) | 通过 |
| 不调阈值（用 EnvMeta default 30/0.75/0.40）| ✅ | ✅ | ✅ | 通过 |
| explanation 仅引用早于目标论文的综述 | ✅ Stolz 2006 / Rosen 2002 / Yin 2011 | ✅ Bond 2000 / Schippers 1999 / Cabrera 2006 | ✅ Johnson & Hallberg 2005 / Falagán 2014 | 通过 |
| 不出现目标论文 specific finding 名词 | ✅ 未提 Asgardarchaeota / 4484-113 | ✅ 未提具体 phylum / site 差异 | ✅ 未提 Coccomyxa / Ca. Acidulodesulfobacterium | 通过 |

**结论**：形式风险**低**。git timestamp 是不可逆证据，reviewer 可独立核验。

### 维度 B — 人类知识泄漏（无法被 git 防护）

⚠️ **这是最难洗的风险点**，我必须诚实承认：

- 我（作者）在写 YAML 之前**已经读过**这 3 篇论文，否则不会选它们做 dataset。
- 即便 explanation 只引用早期综述，潜意识里可能已经倾向选"作者会发现的通路"作为 claim。
- 这种 bias 无法用"我承诺没引用"消除——它是 **dataset selection bias** 的延伸。

**reviewer 会怎么问**：
> "您怎么证明您选这些 claim 不是因为预先知道 EnvMeta 会满足它们？"

**最强的回答**：**Wei 2024 反例**。
- Wei 2024 数据集**作者也读过**，但仍 INSUFFICIENT，因为 claim 写完后 ROCker-only 数据本身覆盖广度不够。
- 这证明"作者读过 → STRONG"不是机械成立的——还要数据本身能 carry 这些 claim。
- 但这只是**部分缓解**，不是完全反驳。最干净的做法是 **由不知情第三方写 YAML**，下次复盘可以请师弟独立写一份对照。

### 维度 C — Claim 设计的代表性（最严肃的批评）

**这一条我必须严厉批评自己**：

| 数据集 | REQUIRED claim | 通路名 | 在该环境**先验上**几乎必然存在吗？ |
|---|---|---|---|
| Liu | arsenate_reduction | arsC1/arsC2 cytosolic As 还原 | ✅ 砷渣/冷泉/任何高 As 几乎必有 |
| Liu | as_transport_detox | arsB/Acr3 efflux | ✅ universal As 抗性，几乎必有 |
| Grettenberger | sulfide_oxidation | sqr/sox AB | ✅ AMD 溪流 backbone，几乎必有 |
| Ayala | dissim_sulfate_reduction | dsrAB/aprAB | ✅ 缺氧高 SO₄² backbone，几乎必有 |

**问题**：所有 REQUIRED claim 都选了**通用 backbone 通路** —— 即在该环境类型里几乎必然检测到的通路。
这意味着只要 KEGG 全注释+任何相关环境，几乎一定 STRONG。

**这降低了"3 STRONG"的诊断意义**。

→ 真正能给评分机制 stress test 的应该是 **risky claim**，比如：
- "深海冷泉**砷甲基化通路应主导**"（实际 As 还原才主导）
- "AMD 溪流**N fixation 应在大多数 MAG 活跃**"（实际偶发）
- "pit lake 深层**砷代谢应活跃**"（实际无砷）

→ 这些都会被 EnvMeta INSUFFICIENT，但我**故意没写**这种 claim，因为我已经知道结果。
→ 如果只看 STRONG/INSUFFICIENT 比例，确实有"挑容易题"嫌疑。

---

## 2. 是否符合一般科研过程的假说设置？

### 合规面 ✅

- 1 REQUIRED + 2-3 exploratory 结构：符合 confirmatory + exploratory 常规科研框架
- 每 claim 引用早期综述：符合"假说基于先验知识"
- weight 分级（1.5 / 1.0 / 0.5 / 0.3）：体现假说强弱区分

### 不符合面 ⚠️

- **真实科研假说通常具有"区分力"**——能区分两个互斥假说（例如"砷固定主要靠 sulfidogenesis 还是 ferrolysis"）。
- 我的 claim 大多是"backbone 通路应活跃"，缺乏区分两个假说的张力。
- 这接近**生物地球化学常识 checklist**，不是"如果 H 真则 Y 应被观察到"的 falsifiable 假说。

→ 这不是 EnvMeta 评分机制的缺陷，而是我**写假说时偷懒**——挑容易满足的来证明评分能给 STRONG。

---

## 3. 测试结果与论文结论一致性

### Liu 2023 npj Biofilms

**论文核心 claim**：1741 As-cycling MAGs；arsenate reduction 主导；arsM widespread (572 active 是论文 highlight)；Asgardarchaeota 为新发现。

| 我的 claim | EnvMeta 结果 | 论文实际结论 | 一致性 |
|---|---|---|---|
| arsenate_reduction REQUIRED | satisfied | arsC1/arsC2 主导 | ✅ 高度一致 |
| as_transport_detox REQUIRED | satisfied | arsB/Acr3 几乎全 MAG | ✅ 高度一致 |
| as_methylation_explored | satisfied (572 active) | 论文 highlight | ✅ 一致；**但**预期是"有限"，结果远超预期 |
| respiratory_as_reduction_explored | satisfied (3 active) | 论文未重点提 | ⚠️ EnvMeta 检测到弱信号，论文未强调 |

**亮点**：as_methylation 我**预期低**（基于 Yin 2011），实测**高**（572 active）→ **数据驱动的真实发现**，不是 fit。

### Grettenberger 2021 AEM

**论文核心 claim**：32 novel MAGs 多样性描述；包含 acidophilic chemolithotrophs；新颖 phyla candidate。

| 我的 claim | EnvMeta 结果 | 论文实际结论 | 一致性 |
|---|---|---|---|
| sulfide_oxidation REQUIRED | satisfied (4 active 83%) | 包含 sulfide oxidizers | ✅ 一致 |
| dissim_sulfate_reduction | satisfied (2 active 75%) | 论文未重点 | ⚠️ EnvMeta 检测到，论文未强调 |
| nitrate_reduction_explored | satisfied (4 active 83%) | 论文未重点 | ⚠️ 同上 |
| nitrogen_fixation_explored | satisfied (4 active 100%) | 论文未重点 | ⚠️ EnvMeta 给出强信号，论文未涉及 |

**关键观察**：Grettenberger 论文核心是**物种多样性描述**，不是**机制验证**。所以我的 claim（机制层面）与论文 claim（物种层面）**不在同一层**。EnvMeta STRONG 不等于"复现论文结论"，更准确地说是"在该数据集上检测到通用先验预期的通路活跃信号"。

### Ayala 2020 Microorganisms（待 GhostKOALA）

**论文核心 claim**（已知）：13 MAGs deep layer；含 dsrAB / aprAB carriers；"硫化物沉淀 not limited by genetic potential"。

| 我的 claim | 预期结果 |
|---|---|
| dissim_sulfate_reduction REQUIRED | 大概率 satisfied（论文核心机制）|
| sulfide_oxidation | 可能 satisfied（IPB 微氧界面）|
| nitrate_reduction_explored | 不确定 |
| nitrogen_fixation_explored | 偶发，不确定 |

---

## 4. 我对论文叙事的修正建议

### 当前 narrative 的脆弱处

> "We tested EnvMeta's hypothesis scoring on **four** independent metagenomic
> datasets ... The contrast among the four ... establishes that EnvMeta's
> INSUFFICIENT label on Wei reflects faithful annotation-coverage diagnostics,
> not tool malfunction."

**reviewer 可能的诘问**：
1. "您所有 STRONG dataset 都用 backbone 通路 claim，这是 confirmation 容易场景，无法证伪"
2. "您 selection 这些 dataset 时已读过论文，怎么排除 selection bias？"
3. "您 STRONG 的诊断意义是什么？是 EnvMeta 工作正常，还是 KEGG 全注释 + 通用 claim 必然 STRONG？"

### 建议的修正叙事（更诚实更稳健）

> We **calibrated** EnvMeta's hypothesis scoring on four metagenomic datasets
> spanning a gradient of annotation breadth. With KEGG-grounded backbone
> pathway claims (deliberately chosen to reflect biogeochemical priors well
> established before the publication years of these papers), EnvMeta returned
> STRONG on all three KEGG-curated datasets and INSUFFICIENT on the
> ROCker-only dataset. This is **calibration evidence** — the scoring engine
> tracks annotation coverage as designed — rather than discrimination evidence.
> A stress test using risky/falsifiable claims (e.g. predicting nitrification
> dominance in AMD) would be a stronger demonstration; we leave such **blind**
> hypothesis-writing experiments to future user studies, where third parties
> unaware of dataset findings would author claims independently.

**核心改动**：
1. "demonstration" → **"calibration"**
2. 承认 STRONG 不等于 "discrimination power"
3. 主动暴露未做的 stress test，将其作为 future work
4. 把 Wei INSUFFICIENT 重新定位为 **calibration anchor**（而非 "证明 EnvMeta 工作正常"）

---

## 5. 改进路线（如果时间允许）

### 立刻可做（~1 h）

更新 [scoring_validation_experiment_results.md](scoring_validation_experiment_results.md) 第 5 节叙事段落，
把 "establishes ... not tool malfunction" 改为 "calibrates ... future stress
tests with blind/risky claims would strengthen this evidence"。

### 中期（建议）

请实验室内**未读过这 3 篇论文的师弟/师妹**独立为 1-2 个数据集写 pre-data
YAML，跑 EnvMeta，看是否仍是 STRONG。这才是真正的盲法 pre-registration。

### 长期（v0.9 / Paper 4 候选）

构造**故意 risky claim 测试集**：在每个数据集上**额外**写 1 条 claim 故意
违背先验（如 "深海冷泉 N fixation 应主导"），看 EnvMeta 是否给 INSUFFICIENT。
这是评分机制 **discrimination power** 的硬证据。

---

## 6. 给自己的诚实总结

| 维度 | 我的 grade | 理由 |
|---|---|---|
| 形式合规（git/阈值/引用）| **A** | 不可逆证据齐全 |
| 人类知识泄漏防护 | **C+** | 作者读过论文，唯一缓解是 Wei 反例和未来盲法 |
| Claim 设计代表性 | **C** | 都挑了 backbone 通路，避开 risky claim |
| 与论文结论一致性 | **B** | Liu 高度一致；Grettenberger 不在同一层；Ayala 待验 |
| 整体方法学贡献 | **B-** | 是 calibration 证据，不是 discrimination 证据，不应过度声称 |

**核心反思**：我把对照实验做成了"pre-registration + KEGG-curated 都 STRONG → 评分机制工作正常"的 confirmation。这在形式上无懈可击，但**实质上是 calibration 不是 stress test**。reviewer 大概率会指出这点；与其被动应答，不如主动在 Discussion 中坦诚承认局限并把 risky-claim 实验列为 future work。

**这种诚实自评本身**就是论文方法学严谨性的证据，比硬撑"4-arm 全胜利" narrative 更稳。
