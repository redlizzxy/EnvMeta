# RDA 排序图验证

## 输入
- `tests/sample_data/Species.txt`（113 种 × 10 样本）
- `tests/sample_data/env_factors.txt`（4 个环境因子：pH / Eh / TOC / Total_As）
- `tests/sample_data/metadata.txt`

## EnvMeta 输出
- `envmeta_rda.pdf` — RDA 排序图（样本散点 + 环境因子箭头）
- `envmeta_rda_stats.tsv` — RDA1/RDA2 解释度 + 每个环境因子的 Mantel r/p + 箭头坐标

## 统计
- Hellinger 变换丰度矩阵
- `skbio.stats.ordination.rda()` 做约束排序
- Bray-Curtis 距离 + `skbio.stats.distance.mantel()` 逐因子做 999 次置换检验
- 本数据集结果：
  - RDA1 解释 20.8%，RDA2 解释 18.1%
  - Total_As Mantel r=0.63（p=0.005，**）— 最强关联
  - Eh Mantel r=0.46（p=0.010，*）
  - TOC Mantel r=0.29（p=0.040，*）
  - pH Mantel r=0.05（p=0.845，ns）

## 样本 ID 对齐策略
env_factors 用别名 SampleID（`CK_1`），abundance 用原始 ID（`2_1`），模块内部按
(Group, 出现顺序) 自动对齐，兼容 metadata → env 的二级命名。

## 与原 R 脚本对比（待补）
原 `scripts/R/03_RDA.R` 使用 `vegan::rda()` + `vegan::mantel()`。
skbio 和 vegan 的 RDA 结果轴方向可能镜像，但解释度和 Mantel 结果数值应一致。
下次装 R 做侧侧对比。
