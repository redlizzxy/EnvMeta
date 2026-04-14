# RDA 排序图验证

## 输入
- `tests/sample_data/Species.txt`（113 种 × 10 样本）
- `tests/sample_data/env_factors.txt`（4 个环境因子：pH / Eh / TOC / Total_As）
- `tests/sample_data/metadata.txt`

## EnvMeta 输出
- `envmeta_rda.pdf` — RDA 排序图（样本散点 + 环境因子箭头 + F 检验显著性）
- `envmeta_rda_stats.tsv` — RDA1/RDA2 解释度（约束方差比例）+ 逐因子 ANOVA-terms F/p + Mantel r/p + 箭头坐标

## 统计方法（对标 R vegan::rda）
- Hellinger 变换丰度矩阵
- `skbio.stats.ordination.rda()` 做约束排序
- **解释度**：默认 `explained_ref="constrained"`，RDA1/RDA2 占约束方差的比例（与 R `summary(rda)$cont$importance[2,]` 一致）
- **逐因子显著性**（主指标）：每个环境因子做**边际置换 F 检验**（999 次置换，随机打乱该因子列后重拟合，对标 R `anova.cca(by="terms")` 思路，数值上等效于 marginal 版；边际 vs 序贯的 p 值对小样本可能不同）
- **Mantel**（辅助指标）：Bray-Curtis 距离对每个环境因子做 999 次置换的距离相关检验

## 本数据集结果
- RDA1 解释 36.93%，RDA2 解释 32.12%（约束方差；R 原脚本在完整数据集上分别是 32.96% / 24.94%）
- 约束方差 / 总方差 = 1.875 / 3.330 = 56.3%
- ANOVA-terms（边际，999 置换）：
  - Total_As F=1.18, p=0.064
  - TOC F=1.25, p=0.057
  - Eh F=0.98, p=0.487
  - pH F=1.12, p=0.216
- Mantel（辅助）：
  - Total_As r=0.63, p=0.003（**）— 最强距离关联
  - Eh r=0.46, p=0.012（*）
  - TOC r=0.29, p=0.033（*）
  - pH r=0.05, p=0.814（ns）

> n=10 / 4 个因子 置换下 F 检验的统计功效有限；原 R 脚本使用完整论文数据（更多物种 + 相同 10 样本）得出 pH** 等较强信号。

## 样本标签策略
- `use_alias_labels=True`（默认）：从 metadata 的 `Group+Replicate` 合成 `CK_1 / A_1 / B_1` 风格标签（对标 R 脚本 `id_map` 的 `2_1 → CK_1`）。
- env_factors 用别名 SampleID（`CK_1`），abundance 用原始 ID（`2_1`），模块内部按
  (Group, 出现顺序) 自动对齐。

## 与原 R 脚本对比
| 项 | R vegan | EnvMeta（skbio） |
|---|---|---|
| RDA 实现 | `vegan::rda(sp_hell ~ ., env_scaled)` | `skbio.stats.ordination.rda(Y, X)` |
| 解释度分母 | 约束方差总和（`cont$importance`） | 同（`explained_ref="constrained"`） |
| 箭头方向 | bp 得分 | biplot_scores（符号约定可能镜像） |
| 因子显著性 | `anova.cca(by="terms")` 序贯 | 边际置换 F 检验（small n 下数值可能不同） |
| Mantel | `vegan::mantel()` | `skbio.stats.distance.mantel()` |
| 标签 | `id_map` 手工映射 | `Group+Replicate` 自动合成 |

装 R 后做侧侧 PDF 对比，补 F 统计量和 p 值对照。
