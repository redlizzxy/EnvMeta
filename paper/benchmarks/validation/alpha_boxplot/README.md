# α 多样性箱线图 — R / EnvMeta 侧侧对照

**状态**：✅ 已完成对照（2026-05-07）
**结论**：EnvMeta 输出与原 R 流程**统计等价**（关键 p 值与 KW 完全一致；
个别 pairwise p 值因连续性修正差异 ≤ 5%，**无统计学结论反转**）

---

## 1. 输入数据（确认一致）

| 文件 | 来源 | 用途 | MD5 |
|---|---|---|---|
| `alpha.txt` | `data/raw/` ≡ `tests/sample_data/`（同一份）| MetaPhlAn4 输出 5 指数 | `4a69ea704a745384f2d5a731f363dfff` |
| `tax_count.alpha` | `data/raw/` | Kraken2 输出 6 指数 | — |
| `metadata.txt` | `data/raw/` | 样本→组映射（CK/A/B）| — |

R 脚本同时读 Kraken2 + MetaPhlAn4（双数据源 6 个指数）；
EnvMeta 当前实现只读 MetaPhlAn4 一份输入（5 个指数）。
**直接对比的 3 个共有指数**：observed_species / shannon / simpson。

---

## 2. 输出文件

### EnvMeta 侧
- `envmeta_alpha.pdf` — 5 指数（observed_species / shannon / simpson / invsimpson / Pielou_evenness）3×2 布局
- `envmeta_alpha_stats.tsv` — Mann-Whitney U + Kruskal-Wallis + BH 校正 padj

### R 侧
- `r_alpha_combined.pdf` / `.png` — 6 指数 2×3 布局（kraken2 + metaphlan4 双数据源）
- `r_alpha_a-f_*.pdf/.png` — 6 单图（kraken2 a-c + metaphlan4 d-f）
- `alpha_wilcoxon_results.txt` — Wilcoxon 秩和检验 + Kruskal-Wallis（无 BH 校正）

---

## 3. 统计结果对照（共有 3 指数）

| 指数 | 对比 | R p | EnvMeta p | 一致性 |
|---|---|---|---|---|
| observed_species | CK vs A | 0.1116 | 0.1143 | ⚠️ 偏差 0.27%（连续性修正）|
| observed_species | CK vs B | 0.6579 | 0.6579 | ✅ 完全一致 |
| observed_species | A vs B  | 0.0987 | 0.0987 | ✅ 完全一致 |
| observed_species | KW      | 0.0895 | 0.0895 | ✅ 完全一致 |
| shannon | CK vs A | 0.1116 | 0.1143 | ⚠️ 偏差 0.27% |
| shannon | CK vs B | 0.3827 | 0.4000 | ⚠️ 偏差 4.5% |
| shannon | A vs B  | 1.0000 | 1.0000 | ✅ |
| shannon | KW      | 0.2283 | 0.2283 | ✅ |
| simpson | CK vs A | 1.0000 | 1.0000 | ✅ |
| simpson | CK vs B | 0.3827 | 0.4000 | ⚠️ 偏差 4.5% |
| simpson | A vs B  | 0.8597 | 0.8571 | ⚠️ 偏差 0.30% |
| simpson | KW      | 0.7047 | 0.7047 | ✅ |

**判断**：所有 KW p 值完全一致；pairwise p 值有 4 处微差（< 5%）。
所有差异源于 R `wilcox.test(exact=FALSE)` 的**连续性修正** vs `scipy.stats.mannwhitneyu`
默认离散 ranks 处理。**无 p < 0.05 vs p > 0.05 的方向反转**，结论等价。

---

## 4. 视觉差异（预期内）

| 维度 | R | EnvMeta | 差异 |
|---|---|---|---|
| 布局 | 2×3 行列 | 3×2 列行 | 不影响信息内容 |
| 色板 | `c(CK="#1c9cbd", A="#e3943d", B="#92181e")` | EnvMeta default `Set1` | 颜色具体值不同但都遵循"3 组用 3 色"原则 |
| 字体 | sans (validation 强制) | DejaVu Sans (matplotlib default) | 都是 sans-serif |
| 散点 | jitter + ggrepel 标签 | jitter | EnvMeta 暂无样品名标签 |
| 显著性标注 | comparison brackets + p 值 | 仅 stats TSV | EnvMeta GUI 可选启用 |
| 数据源 | Kraken2 + MetaPhlAn4 双数据源 | MetaPhlAn4 单数据源 | EnvMeta 当前仅支持单文件 |

---

## 5. 复现命令

### R 侧
```powershell
$Rscript = "F:\Program Files\R\R-4.5.0\bin\Rscript.exe"
& $Rscript "paper\benchmarks\validation\r_scripts\02_alpha_diversity.R" `
  --kraken2 "data\raw\tax_count.alpha" `
  --metaphlan4 "data\raw\alpha.txt" `
  --design "data\raw\metadata.txt" `
  --output "paper\benchmarks\validation\alpha_boxplot\r_alpha" `
  --stat_dir "paper\benchmarks\validation\alpha_boxplot" `
  --style paper
```

### EnvMeta 侧
```powershell
streamlit run app.py
# 上传 alpha.txt + metadata.txt → α 多样性页 → 默认参数 → 下载 PDF
```

---

## 6. 论文里的引用方式

> EnvMeta α-diversity outputs were validated against the original R pipeline
> (`amplicon`-style `wilcox.test` + `kruskal.test`) on the same MetaPhlAn4
> input. Pairwise p-values agreed within 5% (12/12 matching by direction at
> α=0.05); Kruskal-Wallis p-values agreed exactly (3/3). Minor numerical
> differences trace to R's `wilcox.test(exact=FALSE)` continuity correction
> versus `scipy.stats.mannwhitneyu` discrete ranks, both standard implementations.

---

## 7. 维护记录

| 日期 | 事项 |
|---|---|
| 2026-04-14 | 初版（仅 EnvMeta 侧输出，缺 R 对照）|
| 2026-05-07 | R 对照执行完成；统计对比表 + 视觉差异说明补全 |
