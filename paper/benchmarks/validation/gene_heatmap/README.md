# 元素循环基因热图验证

## 输入
- `tests/sample_data/ko_tpm.spf`（51 个知识库 KO × 10 样本 TPM）
- `tests/sample_data/metadata.txt`

## EnvMeta 输出
- `envmeta_heatmap.pdf` — 4 元素合一热图（KO × Group，左侧色块标注元素+通路）
- `envmeta_heatmap_zscore.tsv` — 行 Z-score 归一后的值
- `envmeta_heatmap_rawmeans.tsv` — 每组 TPM 均值原始表

## 数据覆盖
- 知识库共 57 KO（As 17 / N 17 / S 15 / Fe 8）
- eggnog 表里仅 51 个 KO 有数据，其余 6 个缺失（在热图中自动跳过）
- 3 组 CK / A / B 的组均值

## 与原脚本对比
原 `scripts/python/05_gene_heatmap_log2fc.py` 采用 2×2 分元素子图布局（Fig2-8），
EnvMeta 迭代 2 做的是合并单图版本。log2FC（Fig2-9）将在迭代 3 实现。
