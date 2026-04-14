# log2FC 差异柱图验证

## 输入
- `tests/sample_data/ko_tpm.spf`（51 个知识库 KO × 10 样本，TPM）
- `tests/sample_data/metadata.txt`（CK / A / B 三组）

## EnvMeta 输出
- `envmeta_A_vs_CK.pdf`（以及 B_vs_CK、A_vs_B）
- 对应的 `_stats.tsv`（含 ko, gene, pathway, element, mean_a, mean_b, log2fc, t, p, padj, significant, significance）

## 统计设计
- Welch's t-test（`scipy.stats.ttest_ind(equal_var=False)`）
- log2FC = log2((mean_a + 0.5) / (mean_b + 0.5))（加 pseudocount 避免除零）
- BH 校正：所有 51 个 KO 统一 BH FDR 校正
- 显著性：padj < 0.05 且 |log2FC| ≥ 1

## 当前结果
小样本（3-4 per group）使所有 padj > 0.05，没有达到双重显著阈值（padj + log2FC）。
在真实 n=10-20 的研究数据上应能观察到若干显著差异基因。

## 与原 Python 脚本对比
原脚本 `scripts/python/05_gene_heatmap_log2fc.py` 后半段同样使用 Welch's + BH，所以 padj 和 log2FC 数值本应一致。下次可跑原脚本做侧侧对比。
