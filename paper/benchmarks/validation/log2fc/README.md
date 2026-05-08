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

---

## 2026-05-08 Python 脚本对照执行结果 ✅

原脚本同时生成 Fig2-8（热图）和 Fig2-9（log2FC），输出在
`paper/benchmarks/validation/gene_heatmap/` 下：
- `Fig2-9_log2fc_combined.pdf` — 原脚本 log2FC 柱图（3 组对比合一）
- `gene_log2fc_all_comparisons.txt` — 全 51 KO × 3 比较的 log2FC + p + padj

### 数值对照（抽样：K11811 arsH）

| 比较 | 原脚本 | EnvMeta |
|---|---|---|
| `CK_vs_A` | log2FC=-0.6448, p=0.1357, padj=1.0 | `A_vs_CK`: +0.6448 |
| Welch's t-test | ✅ | ✅ |
| BH FDR 校正 | ✅ | ✅ |
| pseudocount +0.5 | ✅ | ✅ |

📌 **方向由用户决定**（不是 bug）：
- 用户传入 `group_a, group_b` → EnvMeta 算 `log2(mean_a / mean_b)`
- xlabel 自动写 `log2({group_a}/{group_b})` 明确方向
- 原脚本固定 CK_vs_A 即 `log2(CK/A)`；EnvMeta 用 group_a=A, group_b=CK 时即 `log2(A/CK)`
- log2FC 绝对值一致，符号取决于参数顺序

**约定建议**：在论文中保持 baseline (e.g. CK) 做 `group_b`（分母），treatment 做 `group_a`（分子），符号语义"treatment vs baseline"。

### 复现命令

见 `gene_heatmap/README.md`（Fig2-8 + Fig2-9 同一脚本）。

### 维护记录

| 日期 | 事项 |
|---|---|
| 2026-04-14 | 初版（仅 EnvMeta 3 组对比 PDF + stats）|
| 2026-05-08 | 原脚本 log2fc 输出 + 数值对照（小样本 padj > 0.05 双方一致；命名方向反转已记录）|
