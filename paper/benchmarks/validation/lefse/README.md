# LEfSe 差异特征验证

## 输入
- `tests/sample_data/Species.txt`（113 种 × 10 样本）
- `tests/sample_data/Genus.txt`（84 属 × 10 样本）
- `tests/sample_data/metadata.txt`

## EnvMeta 输出
- `envmeta_lefse_species.pdf` + `envmeta_lefse_species_stats.tsv` — 种水平
- `envmeta_lefse_genus.pdf` + `envmeta_lefse_genus_stats.tsv` — 属水平

## 统计方法（对标 Segata 2011）
原 R 脚本 `scripts/R/04_LEfSe.R` 只绘制 Galaxy LEfSe **预跑**的 `input.res` 结果；
本模块直接从丰度表做完整 LEfSe 流程（无需 Galaxy 依赖）：

1. 逐 taxon 做 **Kruskal-Wallis** 秩和检验（所有组），保留 p < `alpha_kw`（默认 0.05）
2. 对每个显著 taxon：
   - **Enriched 组** = 组均值最大的组
   - **LDA 效应量** = `log10(1 + lda_scale × |max_mean − min_mean|)`
     （`lda_scale` 默认 1e6，把百分比尺度的丰度差异映射到 ~2–6 的典型 LDA 输出范围）
3. 按 `lda_threshold`（默认 2.0）过滤
4. 可选按 `tax_levels` 限制层级（Genus / Species / …）

与 Galaxy LEfSe 的差异：
- **无 subclass 分层**（子组 Wilcoxon 步骤暂未实现，样本量小时影响有限）
- LDA 效应量用组均值差 log10 近似，而非完整 LDA 判别分析的特征值分解
（简化但保留"显著差异 + 最大丰度组"的核心语义；完全对齐 Galaxy 留作 Phase 2 精细化）

## 本数据集结果（样本 n=10）
- **种水平**：KW α=0.1 + LDA≥2.0 → 39 个显著种
- **属水平**：KW α=0.1 + LDA≥2.0 → 32 个显著属

> 样本量小（CK=3, A=4, B=3），KW 显著性门槛放宽到 0.1；真实研究（n≥10/组）
> 应回到 α=0.05 + LDA≥3.0 的论文标准。

## 与原 R 脚本对比
| 项 | R 04_LEfSe.R | EnvMeta lefse.py |
|---|---|---|
| LDA 计算 | 依赖 Galaxy LEfSe（外部） | 内置（KW + 均值差 log10） |
| 输入 | `input.res` 预处理结果 | 原始丰度表 + metadata |
| 分类层级识别 | `str_count(Feature, ".")` 计点 | `k__/p__/.../s__` 前缀正则 |
| 出图 | ggplot2 水平柱 | matplotlib 水平柱 |
| 多版本输出 | main / full / genus_species 三版 | 单图 + `tax_levels` / `max_features` 参数控制 |
| 一键复现 | 需 Galaxy + R + 多包 | `pip install envmeta` 点 3 下 |

---

## 2026-05-07 R 重绘对照执行结果

R 侧 PDF 输出已生成（基于 `data/raw/input.res` 中 Galaxy LEfSe 的 91 个显著特征）：
- `r_lefse_LDA_main.pdf` — 正文版（属+种，LDA > 4，14 个特征）
- `r_lefse_LDA_genus_species.pdf` — 推荐版（属+种，LDA > 3，37 个特征）
- `r_lefse_LDA_full.pdf` — 完整版（全层级，LDA > 2，91 个特征）
- `LEfSe_significant_features.txt` — 91 显著特征清单（CK 富集 54 / A 富集 22 / B 富集 15）

⚠️ **算法层数值对照（vs Galaxy input.res 的 91 个特征）尚未做** —— 留作下次 session：
- EnvMeta 跑同一个 Species.txt 输出 ~39 个显著种（α=0.1, LDA≥2）
- Galaxy 输出 91 个特征但跨多个 tax level
- 需筛 Galaxy 的 Genus/Species 部分（约 42 个）后逐 feature 比对：Group 标注一致性 + LDA 数值偏差

### 复现命令（R 侧）
```powershell
$Rscript = "F:\Program Files\R\R-4.5.0\bin\Rscript.exe"
& $Rscript "paper\benchmarks\validation\r_scripts\04_LEfSe.R" `
  --input "data\raw\input.res" `
  --output "paper\benchmarks\validation\lefse\r_lefse" `
  --stat_dir "paper\benchmarks\validation\lefse" --style paper
```
