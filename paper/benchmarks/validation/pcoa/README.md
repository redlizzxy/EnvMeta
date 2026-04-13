# PCoA + PERMANOVA 验证

## 输入
- `tests/sample_data/beta_bray.txt`（10×10 Bray-Curtis 距离矩阵）
- `tests/sample_data/metadata.txt`（CK/A/B 三组）

## EnvMeta 输出
- `envmeta_pcoa.pdf` — PCoA 散点图（PC1 vs PC2，含 PERMANOVA 注释）
- `envmeta_pcoa_stats.tsv` — 解释方差 + PERMANOVA 全局结果
- `envmeta_pcoa_pairwise.tsv` — CK vs A / CK vs B / A vs B 两两对比

## 验证要点
- 前两轴解释方差之和 > 30%（测试 `test_pcoa_variance_explained_reasonable` 已断言）
- 三组整体 PERMANOVA p < 0.1（`test_pcoa_permanova_significant`）
- 近似 R² 由 pseudo-F、组数、样本数推算
