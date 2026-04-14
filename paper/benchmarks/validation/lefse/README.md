# LEfSe 差异特征验证

## 输入
- `tests/sample_data/Species.txt`（113 种 × 10 样本）
- `tests/sample_data/Genus.txt`（84 属 × 10 样本）
- `tests/sample_data/metadata.txt`

## EnvMeta 输出
- `envmeta_lefse_species.pdf` + `envmeta_lefse_species_stats.tsv` — 种水平
- `envmeta_lefse_genus.pdf` + `envmeta_lefse_genus_stats.tsv` — 属水平

## 统计方法（对标 Segata 2011）
原 R 脚本 `scripts/R/04_LEfSe.R` 只绘制 Galaxy LEfSe **预跑**的 `input.res` 结果；
本模块直接从丰度表做完整 LEfSe 流程（无需 Galaxy 依赖）：

1. 逐 taxon 做 **Kruskal-Wallis** 秩和检验（所有组），保留 p < `alpha_kw`（默认 0.05）
2. 对每个显著 taxon：
   - **Enriched 组** = 组均值最大的组
   - **LDA 效应量** = `log10(1 + lda_scale × |max_mean − min_mean|)`
     （`lda_scale` 默认 1e6，把百分比尺度的丰度差异映射到 ~2–6 的典型 LDA 输出范围）
3. 按 `lda_threshold`（默认 2.0）过滤
4. 可选按 `tax_levels` 限制层级（Genus / Species / …）

与 Galaxy LEfSe 的差异：
- **无 subclass 分层**（子组 Wilcoxon 步骤暂未实现，样本量小时影响有限）
- LDA 效应量用组均值差 log10 近似，而非完整 LDA 判别分析的特征值分解
（简化但保留"显著差异 + 最大丰度组"的核心语义；完全对齐 Galaxy 留作 Phase 2 精细化）

## 本数据集结果（样本 n=10）
- **种水平**：KW α=0.1 + LDA≥2.0 → 39 个显著种
- **属水平**：KW α=0.1 + LDA≥2.0 → 32 个显著属

> 样本量小（CK=3, A=4, B=3），KW 显著性门槛放宽到 0.1；真实研究（n≥10/组）
> 应回到 α=0.05 + LDA≥3.0 的论文标准。

## 与原 R 脚本对比
| 项 | R 04_LEfSe.R | EnvMeta lefse.py |
|---|---|---|
| LDA 计算 | 依赖 Galaxy LEfSe（外部） | 内置（KW + 均值差 log10） |
| 输入 | `input.res` 预处理结果 | 原始丰度表 + metadata |
| 分类层级识别 | `str_count(Feature, ".")` 计点 | `k__/p__/.../s__` 前缀正则 |
| 出图 | ggplot2 水平柱 | matplotlib 水平柱 |
| 多版本输出 | main / full / genus_species 三版 | 单图 + `tax_levels` / `max_features` 参数控制 |
| 一键复现 | 需 Galaxy + R + 多包 | `pip install envmeta` 点 3 下 |

装 R + Galaxy LEfSe 后可做完整算法对齐验证。
