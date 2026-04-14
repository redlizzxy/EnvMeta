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

## 文件结构
- `envmeta/geocycle/model.py`：dataclass
- `envmeta/geocycle/inference.py`：推断引擎
- `envmeta/geocycle/renderer.py`：matplotlib 静态渲染
- `envmeta/analysis/cycle_diagram.py`：对齐 analysis 接口的薄壳

Phase 4 将把 `CycleData` 导出为 JSON → D3.js 交互编辑。
