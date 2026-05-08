# RDA — R / EnvMeta 侧侧对照

**状态**：✅ 已完成对照（2026-05-08，bug 修复后）
**结论**：EnvMeta 输出与 R vegan **数值精确一致**。

| 指标 | R vegan | EnvMeta | 一致性 |
|---|---|---|---|
| Total inertia | 0.1617 | 0.16168 | ✅ 0% |
| Constrained inertia | 0.1125 | 0.11248 | ✅ 0% |
| Constrained / Total | 69.6% | 69.5% | ✅ |
| RDA1 explained | 32.96% | 32.96% | ✅ 0% |
| RDA2 explained | 24.94% | 24.94% | ✅ 0% |
| pH ANOVA F | 3.905 | 3.905 | ✅ 0% |
| Eh ANOVA F | 3.877 | 3.877 | ✅ 0% |
| TOC ANOVA F | 1.945 | 1.945 | ✅ 0% |
| Total_As ANOVA F | 1.708 | 1.708 | ✅ 0% |
| pH Mantel r | 0.347 | 0.347 | ✅ 0% |
| Eh Mantel r | 0.288 | 0.288 | ✅ 0% |
| TOC Mantel r | 0.620 | 0.620 | ✅ 0% |
| Total_As Mantel r | 0.705 | 0.705 | ✅ 0% |

p 值因 999 次置换 RNG 随机性有微差（< 0.015），但**所有因子 α=0.05 显著性方向一致**：
pH/Eh 显著，TOC 边际，Total_As 显著。

---

## 1. 输入数据

| 文件 | 来源 | 用途 |
|---|---|---|
| `Species.txt` | `data/raw/`（113 种 × 10 样本）| 物种丰度矩阵 |
| `env_factors.txt` | `data/raw/`（4 因子：pH/Eh/TOC/Total_As）| 环境因子表 |
| `metadata.txt` | `data/raw/` | SampleID→组（CK/A/B）|

样本命名差异（已自动处理）：
- `Species.txt` 列名 `2_1, 2_5, 8_1...`（原始测序 ID）
- `env_factors.txt` SampleID `CK_1, CK_5, A_1...`（重命名后）

R 用 hard-coded id_map；EnvMeta `_align_env_to_abundance()` 按 (Group, 出现顺序) 匹配 — 双方数值一致证明对齐正确。

---

## 2. 修复历史

### Bug 发现（2026-05-07）

EnvMeta 用 `skbio.stats.ordination.rda()` 直接拿 eigvals，与 R `vegan::rda()`
归一化方式不同，导致 inertia 数值差 16-20×、ANOVA F 差 3.5×、p 值方向反转。

### 修复方案 B（2026-05-08）

保留 `skbio.rda()` 输出 ordination axes（site / biplot scores 视觉 OK），
但**自己按 vegan 公式重算 inertia + ANOVA**：

```python
# Vegan-equivalent inertia
def _ss(arr): return np.sum(arr ** 2) / (n_samples - 1)
Y_centered = Y - Y.mean(axis=0)
total_inertia = _ss(Y_centered)

# Constrained via lstsq
beta, _, _, _ = np.linalg.lstsq(X, Y_centered, rcond=None)
Y_fitted = X @ beta
constrained_inertia = _ss(Y_fitted)

# RDA axes eigvals via SVD on Y_fitted
_, s_vals, _ = np.linalg.svd(Y_fitted, full_matrices=False)
constrained_eigvals = (s_vals ** 2) / (n_samples - 1)

# Sequential ANOVA (Type I, vegan by="terms")
def _seq_contrib(X_in, j):
    before = _ss(X_in[:, :j] @ lstsq_fit) if j > 0 else 0
    after = _ss(X_in[:, :j+1] @ lstsq_fit)
    return after - before
```

详见 `envmeta/analysis/rda.py` 第 160-243 行。

**关键 default 调整**：`explained_ref` 默认从 `"constrained"` 改为 `"total"`，
对标 R vegan `summary()$cont$importance[2,]` 默认行为。

---

## 3. 输出文件

### EnvMeta 侧
- `envmeta_rda.pdf`
- `envmeta_rda_stats.tsv`（修复后重新生成）

### R 侧
- `r_rda.pdf` / `.png`（带样品标签）
- `r_rda_detail.pdf` / `.png`（详细版）
- `RDA_results.txt`（完整 vegan 输出 + Mantel）

---

## 4. 复现命令

### R 侧
```powershell
$Rscript = "F:\Program Files\R\R-4.5.0\bin\Rscript.exe"
& $Rscript "paper\benchmarks\validation\r_scripts\03_RDA.R" `
  --species "data\raw\Species.txt" `
  --env "data\raw\env_factors.txt" `
  --output "paper\benchmarks\validation\rda\r_rda" `
  --stat_dir "paper\benchmarks\validation\rda" --style paper
```

### EnvMeta 侧
```powershell
streamlit run app.py
# 上传 Species.txt + env_factors.txt + metadata.txt → RDA 页 → 默认参数 → 下载 PDF
```

---

## 5. 论文引用模板

> EnvMeta RDA outputs were validated against vegan's `rda()` with Hellinger
> transformation and Type I sequential ANOVA on identical input matrices.
> Constrained and total inertia, per-axis explained variance, and per-factor
> ANOVA F statistics matched vegan to four decimal places (< 0.01% deviation).
> Permutation p-values agreed within ±0.015 (999 permutations, with all four
> factors maintaining the same α=0.05 significance direction). Mantel
> correlations between Hellinger-Bray-Curtis species distance and per-factor
> Euclidean distance matched vegan exactly.

---

## 6. 维护记录

| 日期 | 事项 |
|---|---|
| 2026-04-14 | 初版（仅 EnvMeta 侧输出，缺 R 对照）|
| 2026-05-07 | R 对照完成，发现 inertia 16-20× 差 + ANOVA p 反转 |
| 2026-05-08 | 修复方案 B：SS-based 重算 inertia + ANOVA + explained_ref 默认改 "total"；F/r/解释度 4 位精度对齐 R vegan；293/293 测试全绿 |
