# 生物地球化学循环图 v2 — 合并细胞级联 + 描述性推断

## v2 升级（2026-04-15, Mockup 10 落地）

原 v1 `envmeta_cycle_v1.pdf` 为 2×2 元素抽象柱图；v2 用"合并细胞 + 基因级联"
替代：每元素象限挑选 Top-K (通路, MAG) 对，每个细胞显示一个 MAG 在该通路上
实际持有的 KO 链（底物 → 酶 → 中间产物 → … → 最终产物），读取
KB v1.1 新增的 substrate/product 字段驱动。

- `envmeta_cycle_v2.pdf` — 新版级联渲染（默认输出）
- `envmeta_cycle_v2.png` — 新版位图预览
- `envmeta_cycle_v2_stats.tsv` — stats 长表
- `envmeta_cycle_v1.pdf` — 保留 v1 供对比，亦可用 ``params={"cell_mode": "bars"}`` 回退

**S2.5-3 完成**（2026-04-15）：化学物耦合线已接入。本数据实际画出 3 条：
- As(III) + S-2 → As₂S₃（紫，precipitation）
- As(V) + Fe(III) → Fe-As_surface（棕橙，adsorption）
- NO3- + As(III) → As(V)（绿，redox）

Fe(II) + S-2 → FeS（黑蓝）未出现，因当前数据 Fe(II) 不作为任一细胞产物；
若用户加载含 dissim Fe 还原（Fe(III)→Fe(II)）的 KO 数据应会显示。

**Phase 3 S2.5 后续**：S2.5-4 组选择下拉（CK/A/B 单组模式）；S2.5-5 SVG 导出。

---

## MAG selection criteria（S2.5-9，`contributor_ranking`）

"对循环最关键"的 MAG 没有唯一定义。本工具提供 4 档可切换判据：

| ranking | 回答的问题 | 排序 key | 过滤 |
|---|---|---|---|
| `abundance`（默认） | 谁在占主导（丰度视角） | `completeness × log1p(ab)` | — |
| `completeness` | 谁对这条通路覆盖最完整 | `completeness` | — |
| `keystone_only` | 网络枢纽视角（忽略高丰度非 keystone） | 同 abundance | `is_keystone==True` |
| `keystone_priority` | 优先看 keystone，但不排除高丰度非 keystone | `base × (10 if keystone else 1)` | — |

Streamlit 侧栏 radio 切换。非默认判据时 subtitle 会显示 `ranking=X` 标签。

### 论文写作建议
主图用 `abundance` 展示"谁占主导"；SI 附 `keystone_only` 版本做对比，叙述：
> "Under the abundance-weighted view (main Fig X), Arsenate reduction was dominated
> by *Gallionella*. However, when restricted to keystone species (Fig S1), the same
> pathway is shown to be carried by *SPCO01* (a co-occurrence hub), suggesting a
> functional dependency that quantity-based analysis would miss."

---

## MAG 显示标签格式（S2.5-13）

为了让用户一眼看出不同 MAG 是"同属不同菌株"而非"同种未合并"，标签遵循：

| 情形 | 格式 | 例 |
|---|---|---|
| GTDB `g__Genus;s__species` | `Genus species` | `Sulfuricaulis marinus` |
| GTDB `g__Genus;s__` 空（多数 MAG） | `Genus sp. Mx_XX` | `Sulfuricaulis sp. Mx_141` |
| 无分类 | `MAG_id` | `Mx_All_141` |

这样 `Sulfuricaulis sp. Mx_141` 和 `Sulfuricaulis sp. Mx_110` 清楚呈现为"同属 Sulfuricaulis 的两个不同菌株"，而不是"两个同种未合并"。

---

## 为什么 Fe uptake regulation 没有底物 / 产物（S2.5-13）

`fur` (K03711) 和 `tonB` (K03832) 是**转录调控因子**（transcriptional regulators），不催化"底物→产物"的化学反应。它们通过结合 DNA 和调节其他基因的表达来控制铁摄取——**没有化学催化活性**。因此 KEGG 官方数据里 substrate / product 都是空，我们的 KB 也相应标记为 null。在循环图上表现为"细胞膜内只有酶椭圆、没有穿膜底物产物箭头"，这是**正确的生物学表达**。

如果这种调控型 cell 干扰阅读，可在 Streamlit 勾选"🧪 隐藏纯调控型 cell"，工具会跳过所有基因都是调控型的 cell（如单独的 fur+tonB 组合），使画面聚焦于催化型通路。保留"同 MAG 承载 Fe uptake regulation + Fe transport"这类混合 cell（因为 Fe transport 提供了 Fe(III)→Fe_internal 的催化内容）。

---

## 如何叙述 A/CK/B 三组循环图的差异（S2.5-7d）

用户常见误解：看三张单组图看起来 "通路潜力一样"，觉得图没用。
**这是 KB 覆盖率偏高导致的视觉同构**，真实差别藏在三个维度：

### 维度 1：承载丰度（primary）
同一通路在不同组的 `total_contribution` 反映**谁真的占据这条通路的生态位**。
例如 `Arsenate reduction` 在 A 组 total_contribution=1200，CK 组=400 →
加钢渣组砷还原能力是对照组的 3 倍，尽管两组都"有"这条通路。

在 Streamlit 循环图页底部"跨组对比表"一键生成，TSV 可直接进论文 Table S1。

### 维度 2：承载者身份（secondary）
对比 `top_mag_genus` 跨组变化：
- 同一通路在 CK 组由 A 物种承载，在 A 组换成 B 物种 → 说明处理改变了微生物承载结构
- `top_mag_is_keystone` 标记进一步强化：keystone 级的承载者说服力更强

### 维度 3：环境耦合（组内不独立）
单组循环图的 env 面板多为空，因为组内 n≤4 相关性极不稳定。
**用全组（group=All）那张图读 env 耦合**，配合 S2 置换零假设的可信度标签
（strong / suggestive / spurious?）叙述哪些相关经得起统计考验。

### 论文写作建议
> "While the qualitative pathway potentials were largely conserved across treatments
> (Fig. X), abundance-weighted contribution analysis revealed quantitative shifts:
> Arsenate reduction activity in high-slag treatment (A) was 3× that of control (CK)
> (Table S1). Moreover, the top-contributing MAG for As reduction shifted from
> *Genus_A* (CK) to *Gallionella* (A), a keystone species identified by
> co-occurrence network analysis (Section Z)."

---

# 生物地球化学循环图 v1 — 自动描述性推断

## 重要声明（Interpretation Guidance）

**EnvMeta 输出是描述性的，不是因果性的。**

本模块报告：
- 哪些通路的完整度 ≥ 阈值（事实）
- 每条通路的最高完整度/最高丰度 MAG 是谁（事实）
- 哪些环境因子与通路活性在统计上有相关（事实）

本模块**不报告**：
- 谁"驱动"了什么（因果）
- 哪个机制成立（假说评估）
- 该研究应该得出什么结论（科学判断）

**因果解读需要领域专家基于这些描述性证据自行判断。** 用户和审稿人**务必**将本图
视为"假说生成器"而非"假说验证器"。若需要对特定机制做评分，请使用 Session 3 的
"机制 YAML 评分器"（用户自带假说 YAML + 显式证据链）。

---

## 输入
- `tests/sample_data/kegg_target_only.tsv`（MAG × KO 长表）
- `tests/sample_data/mag_taxonomy_labels.tsv`
- `tests/sample_data/keystone_species.txt`
- `tests/sample_data/abundance.tsv`（MAG × 样本丰度，用于通路贡献加权）
- `tests/sample_data/env_factors.txt`（pH / Eh / TOC / Total_As）
- `tests/sample_data/metadata.txt`

## EnvMeta 输出
- `envmeta_cycle_v1.pdf` — 2×2 元素象限（As/N/S/Fe）+ 环境耦合面板
- `envmeta_cycle_v1_stats.tsv` — 活跃通路长表 + env-pathway 相关性（过滤后）
- `envmeta_cycle_v1_full_corr.tsv` — **完整**环境相关矩阵（不过滤，供用户自主判读）
- `envmeta_cycle_v1_sensitivity.tsv` — 完整度阈值扫描（30/50/70%）的 Top-1 稳定性

## 推断结果（基于论文真实数据）— 描述性报告

### 元素循环活跃通路（completeness ≥ 50%）
| 元素 | 活跃/总 | 最活跃通路（按贡献）| Top-completeness contributor |
|---|---|---|---|
| As | 5/6 | As regulation (133 MAG) | Fen-1038 |
| N | 6/6 | NO reduction (64 MAG) | Gallionella |
| S | 4/4 | Thiosulfate metab. (113 MAG) | DASXPG01；Sulfuricaulis 是 Sulfide oxidation 的 top contributor |
| Fe | 2/2 | Fe uptake regulation (153 MAG) | DASXPG01 |

> **Top-completeness contributor** = 在该通路里 `completeness × log1p(abundance)` 最高的 MAG，
> **不等于**"驱动该通路的物种"。完整度高只说明该 MAG 拥有该通路的大部分 KO；
> 丰度高说明它在样本中多。两者都不构成**生态驱动性**的因果证据。

### 可信度分层（S2 置换零假设验证）

对每条环境-通路相关跑 999 次置换零假设检验（打乱 env 重算 ρ），打上可信度标签：

| 标签 | 判据 | 数量 |
|---|---|---|
| **strong** | \|ρ\|>0.7 + perm_p<0.01 + 敏感度稳健 | **1** |
| **suggestive** | 0.5<\|ρ\|≤0.7 + perm_p<0.05 | 10 |
| **weak** | 0.3<\|ρ\|≤0.5 + perm_p<0.1 | — |
| **spurious?** | \|ρ\|≥0.5 但 perm_p>0.05（超阈值但置换不支持）| **5** |

**唯一的 strong 发现**：Ammonia oxidation ↔ Eh（ρ=0.835, perm_p<0.01, robust sensitivity）

**5 条 spurious?**：S1 讨论里怀疑的"ρ=0.764 共变"假设**被置换检验证实** — 5 条通路
共享同一 ρ=0.764 vs Total_As 的相关不经得起置换检验，提示它们**不是各自独立响应
As 压力**，而是追踪某个共同样本级信号。

### 环境-通路 Spearman 相关（已过滤：|ρ| ≥ 0.5，p < 0.05）
| 通路 | 环境因子 | Spearman ρ | p |
|---|---|---|---|
| Arsenate reduction | Total_As | 0.854 | 0.002 ** |
| Ammonia oxidation | Eh | 0.835 | 0.003 ** |
| As transport/detox | Total_As | 0.809 | 0.005 ** |
| As methylation | Total_As | 0.764 | 0.010 * |
| N fixation | Eh | 0.762 | 0.010 * |

> **相关 ≠ 因果**。ρ=0.85 可以来自：(a) 通路确实响应 As 压力；
> (b) 两者共同受其他因子驱动；(c) 小样本（n=10）偶然相关。
> **Session 2 将加入置换零假设检验 + 可信度标签**，把 (c) 排除。

### 完整相关矩阵（供用户自主判读，不过滤）
见 `envmeta_cycle_v1_full_corr.tsv`。**注意观察同一通路与多个环境因子同时相关的情况**，
这常提示潜在共变（confounding）。

### 阈值敏感度
见 `envmeta_cycle_v1_sensitivity.tsv`。每条通路的 Top-1 contributor 在 30/50/70%
三档完整度阈值下是否稳定：
- `robust`：三档 Top-1 一致 → 该发现对阈值不敏感
- `threshold-sensitive`：Top-1 随阈值变化 → 该发现不稳健，慎重引用

## 该数据集的观察（用户可参考但请自行判断）

描述层面，本样本的证据与"铁氧化固砷、硫/氮循环调控 Eh"类机制**相容**，原因：
- Gallionella（Fe 氧化代表属）同时在 NO reduction 和 Arsenate reduction 中
  完整度较高 + 丰度较高
- Sulfuricaulis（S 氧化代表属）是 Sulfide oxidation 的 top contributor
- Ammonia oxidation 与 Eh 正相关（ρ=0.84, p=0.003）
- 5/5 砷代谢通路与 Total_As 正相关

**但这些只是相容证据，不是机制证明**。举例，Arsenate reduction ↔ Total_As 的
高 ρ 也可能是因为"高 As 环境选择了有 arsC 的 MAG 定殖"，而不是"arsC 表达改变了 As 形态"。
机制判断需要：
- 转录组 / 蛋白组数据（而非仅 KO 基因存在）
- 时序数据（而非截面）
- 富集实验证伪（如敲除 Gallionella 观察 As 变化）

**Session 3 的机制 YAML 评分器**将让用户把具体假说写成结构化证据规则，EnvMeta 对
每条规则独立评分，输出`supported / partial / refuted`；这比当前的描述性相关更可靠。

## 算法

1. 每个 MAG 按 KB 18 通路计算完整度（复用 `pathway.analyze`）
2. 通路"活跃"判定：完整度 ≥ 阈值（默认 50%，可调）
3. 活跃 MAG 排序：`completeness × log1p(abundance_mean)` 降序
4. 通路贡献：`Σ (completeness × abundance_mean)`
5. env-pathway 相关性：
   - 样本级通路活性 = `Σ_MAG (abundance_MAG × |MAG 的 KO ∩ 通路 KO|)`
   - Spearman 相关于每个环境因子
   - 过滤 `|ρ| ≥ 0.5` 且 `p ≤ 0.05`
6. **完整相关矩阵**（不过滤）同时输出，供用户判读 confounding
7. **敏感度扫描**：默认 thresholds=[30, 50, 70]%，记录每通路 Top-1 是否稳定
8. **置换零假设检验**（S2 新增）：每条相关跑 999 次置换，empirical p 替代 scipy p
9. **可信度标签**（S2 新增）：综合 |ρ| + perm_p + 敏感度打
   `strong / suggestive / weak / spurious?` 标签，论文可按标签选择引用强度

## 文件结构
- `envmeta/geocycle/model.py`：dataclass
- `envmeta/geocycle/inference.py`：推断引擎
- `envmeta/geocycle/renderer.py`：matplotlib 静态渲染
- `envmeta/analysis/cycle_diagram.py`：对齐 analysis 接口的薄壳

Phase 4 将把 `CycleData` 导出为 JSON → D3.js 交互编辑。
