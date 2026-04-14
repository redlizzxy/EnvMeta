# α 多样性箱线图验证

## 输入
- `tests/sample_data/alpha.txt`（5 个指数 × 10 样本：observed_species, shannon, simpson, invsimpson, Pielou_evenness）
- `tests/sample_data/metadata.txt`

## EnvMeta 输出
- `envmeta_alpha.pdf` — 3 列 × 2 行布局的 5 个箱线图
- `envmeta_alpha_stats.tsv` — Kruskal-Wallis + 两两 Mann-Whitney U + BH 校正

## 统计结果说明
- 每个指数：15 行总（5 指数 × 3 组对 CK/A, CK/B, A/B）
- 小样本（n=3-4）导致本数据集上多数 padj > 0.05
- `padj` 列是 BH 校正后的 p，`significance` 列是标签（`***/**/*/空`）

## 与原 R 脚本对比（待补）
`scripts/R/02_alpha_diversity.R` 用 Wilcoxon 且不做全局 BH 校正，做到多指数组合后 BH 校正是对齐论文写作习惯的做法。
装上 R 后跑 PDF 做侧侧对比（迭代 4 阻塞项）。
