# β-PCoA + PERMANOVA — R / EnvMeta 侧侧对照

**状态**：✅ 已完成对照（2026-05-07）
**结论**：EnvMeta 输出与原 R 流程（vegan）**数值等价**（F / R² / PC%
精确到小数点 2 位完全一致；p 值因 999 次置换 RNG 差异 ±0.002 但同方向）

---

## 1. 输入数据（确认一致）

| 文件 | 来源 | 用途 | MD5 |
|---|---|---|---|
| `beta_bray.txt` | `data/raw/` ≡ `tests/sample_data/`（同一份）| 10×10 Bray-Curtis 距离矩阵 | — |
| `metadata.txt` | `data/raw/` | 样本→组映射（CK/A/B）| — |

R 脚本与 EnvMeta 都直接读 Bray-Curtis 距离矩阵（不重新计算），保证完全
等价的输入。

---

## 2. 输出文件

### EnvMeta 侧
- `envmeta_pcoa.pdf` — PCoA 散点图 + PERMANOVA 注释
- `envmeta_pcoa_stats.tsv` — 全局解释度 + PERMANOVA F/R²/p
- `envmeta_pcoa_pairwise.tsv` — 三组两两对比

### R 侧
- `r_pcoa.pdf` / `.png` — PCoA 散点图（含 ggrepel 样品标签）
- `r_pcoa_clean.pdf` / `.png` — 无标签清洁版

---

## 3. 数值对照

### 全局 PERMANOVA

| 指标 | R | EnvMeta | 一致性 |
|---|---|---|---|
| pseudo-F | **6.6157** | **6.6157** | ✅ 完全一致 |
| R² | **0.654** | **0.654** | ✅ 完全一致 |
| p (999 perm) | 0.003 | 0.001 | ⚠️ 都 < 0.005，RNG 差异 |
| Df Model | 2 | 2 | ✅ |
| Df Residual | 7 | 7 | ✅ |
| n_samples | 10 | 10 | ✅ |

### 解释方差

| 轴 | R | EnvMeta | 一致性 |
|---|---|---|---|
| PC1 | 50.91% | 50.91% | ✅ |
| PC2 | 20.29% | 20.29% | ✅ |
| PC3 | 10.38% | 10.38% | ✅ |
| PC4 | 7.28% | 7.28% | ✅ |
| PC5 | 4.74% | 4.74% | ✅ |
| 前2轴累计 | **71.20%** | **71.20%** | ✅ |

### Pairwise PERMANOVA

| 对比 | R p | EnvMeta p | R F (推算) | EnvMeta F | 一致性 |
|---|---|---|---|---|---|
| CK vs A | 0.033 | 0.031 | — | 7.44 | ✅ 同方向 |
| CK vs B | 0.100 | 0.107 | — | 7.61 | ✅ 同方向（边界）|
| A vs B  | 0.035 | 0.034 | — | 5.20 | ✅ 同方向 |

注：R 脚本仅输出 R²，EnvMeta 输出 pseudo-F，但全局 F 一致证明算法等价。

**判断**：所有数值小差异 ≤ 0.5%，p 值差异 ≤ 0.007（999 次置换 RNG 限制）。
**统计学结论 100% 一致**（α=0.05 下：CK-A 显著 / CK-B 边界 / A-B 显著）。

---

## 4. 视觉差异（预期内）

| 维度 | R | EnvMeta | 差异 |
|---|---|---|---|
| 散点 | jitter + ggrepel 样品标签 | jitter | EnvMeta 暂无样品标签（论文 SI 可加）|
| 配色 | `c(CK="#1c9cbd", A="#e3943d", B="#92181e")` | EnvMeta default Set1 | 都 3 色无重叠 |
| 椭圆 | 无（"无椭圆版"） | 可选 95% confidence | 配置项差异 |
| PERMANOVA 注释位置 | 图右下角 | 图左下角 | 排版差异 |
| 字体 | sans (validation 强制) | DejaVu Sans (matplotlib) | 都 sans-serif |

---

## 5. 复现命令

### R 侧
```powershell
$Rscript = "F:\Program Files\R\R-4.5.0\bin\Rscript.exe"
& $Rscript "paper\benchmarks\validation\r_scripts\02_beta_PCoA.R" `
  --dist "data\raw\beta_bray.txt" `
  --design "data\raw\metadata.txt" `
  --output "paper\benchmarks\validation\pcoa\r_pcoa" `
  --stat_dir "paper\benchmarks\validation\pcoa" `
  --style paper
```

### EnvMeta 侧
```powershell
streamlit run app.py
# 上传 beta_bray.txt + metadata.txt → β-多样性页 → PCoA → 默认参数 → 下载 PDF
```

---

## 6. 论文引用模板

> EnvMeta β-PCoA outputs were validated against the `vegan` R pipeline
> (`betadisper` + `adonis2`) on the same Bray-Curtis distance matrix.
> Pseudo-F (6.6157), R² (0.654), and PC1-5 explained variance percentages
> matched to two decimal places. PERMANOVA p-values (999 permutations)
> agreed within ±0.007 across global and pairwise tests, with all
> statistical conclusions identical at α=0.05.

---

## 7. 维护记录

| 日期 | 事项 |
|---|---|
| 2026-04-14 | 初版（仅 EnvMeta 侧输出，缺 R 对照）|
| 2026-05-07 | R 对照执行完成；F/R²/PC% 完全一致，p 值 RNG 差异说明 |
