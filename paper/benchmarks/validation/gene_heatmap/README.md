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

---

## 2026-05-08 Python 脚本对照执行结果 ✅

原 Python 脚本输出已生成：
- `Fig2-8_heatmap_combined.pdf` / `.png` — 2×2 元素子图布局
- `Fig2-9_log2fc_combined.pdf` / `.png` — 3 组对比 log2FC（CK_vs_A / A_vs_B / CK_vs_B）
- `FigS_arsenic_overview.pdf` / `.png` — 砷代谢补充图
- `FigS_iron_overview.pdf` / `.png` — 铁代谢补充图
- `FigS_nitrogen_overview.pdf` / `.png` — 氮代谢补充图
- `FigS_sulfur_overview.pdf` / `.png` — 硫代谢补充图
- `gene_log2fc_all_comparisons.txt` — 完整 log2FC + p + padj 表

### 数据覆盖一致性

- ✅ 51 KO（57 - 6 全 0）双方一致
- ✅ 3 组 CK / A / B 双方一致
- ✅ 元素分类（As/N/S/Fe）双方一致

### log2FC stats 对照（抽样）

| 指标 | 原脚本 | EnvMeta | 一致性 |
|---|---|---|---|
| K11811 arsH (CK→A) log2FC | -0.6448 | (envmeta A_vs_CK) +0.6448 ⚠️方向 | 数值一致，符号反转（CK_vs_A vs A_vs_CK 命名差异）|
| Welch's t-test | ✅ | ✅ | 双方使用同一统计 |
| BH FDR 校正 | ✅ all 51 KO | ✅ all 51 KO | 一致 |
| pseudocount +0.5 | ✅ | ✅ | 一致 |

⚠️ **方向命名约定差异**：
- 原脚本：列 `Comparison=CK_vs_A` 表示 `log2(mean_CK / mean_A)`
- EnvMeta：文件名 `A_vs_CK` 表示 `log2(mean_A / mean_CK)`

数值绝对值应完全一致，符号反转。论文里需要在 figure caption 明确指明
比较方向，避免审稿人误读。

### Bonus：4 个 FigS 元素补充图

EnvMeta 当前未提供按元素拆分的补充图（FigS_arsenic / iron / nitrogen / sulfur
overview）。这是原 Python 脚本的额外产物，可作为：
- ✅ 论文 SI 直接引用（已生成）
- ⚠️ 或加入 EnvMeta 作为可选导出（"按元素分组导出 4 张图"功能）

### 复现命令

```powershell
$py = "C:\Users\REDLIZZ\.conda\envs\envmeta\python.exe"
& $py "d:\workdata\envmeta_thesis\scripts\python\05_gene_heatmap_log2fc.py" `
  -i "data\raw\eggnog.KEGG_ko.TPM.spf" `
  --metadata "data\raw\metadata.txt" `
  -o "paper\benchmarks\validation\gene_heatmap" `
  --stat_dir "paper\benchmarks\validation\gene_heatmap" --style paper
```

### 维护记录

| 日期 | 事项 |
|---|---|
| 2026-04-14 | 初版（仅 EnvMeta 输出，缺 log2FC）|
| 2026-05-08 | 原脚本输出 + 4 FigS + log2fc 完整对照；方向命名差异已识别 |
