# 堆叠图验证

## 输入
- `tests/sample_data/Phylum.txt` (10 门 × 10 样本，CK/A/B 三组)
- `tests/sample_data/metadata.txt`

## EnvMeta 输出
- `envmeta_phylum_sample.pdf` — 横轴=样本（10 柱）
- `envmeta_phylum_group.pdf` — 横轴=分组（3 柱）
- `envmeta_phylum_{style}_pct.tsv` — 百分比矩阵（每列和为 100）

## 与原 R 脚本对比（待补）
原脚本 `scripts/R/01_tax_stackplot.R` 需要 R 环境（vegan/cowplot）。
本机暂未装 R，迭代 2 装上 R 后再补 `R_original_phylum.pdf` 做侧侧对比。

## 结果一致性检查
两种样式下，所有列百分比之和精确等于 100（test_stackplot.py 已验证）。
Top-10 + Others（如有）的处理符合 R 脚本的 `tax_stackplot(sorted="abundance")` 行为。
