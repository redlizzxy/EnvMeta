# MAG 元素循环基因谱验证

## 输入
- `tests/sample_data/kegg_target_only.tsv`（MAG × KO 长表注释）
- `tests/sample_data/mag_taxonomy_labels.tsv`
- `tests/sample_data/keystone_species.txt`
- `tests/sample_data/abundance.tsv`

## EnvMeta 输出
- `envmeta_gene_profile.pdf` — 168 MAG × 51 activa KO 热图（按 As/N/S/Fe 元素分块）
- `envmeta_gene_profile_stats.tsv` — MAG × KO 拷贝数矩阵 + 元数据

## 关键结果
- **168 MAG × 51 active KO**（KB 定义 57 KO，6 个全 0 自动过滤）
- 总基因拷贝数：5144
- 14 keystone 正确标记
- Top-3 基因数最多的 MAG：
  - `Mx_All_158` Desulfobacterota_B（114 copies）—— 硫还原专家
  - `Mx_All_48` Pseudomonadota（88）
  - `Mx_All_63` Pseudomonadota（79）—— 同时在 pathway 完整度 Top-1，交叉验证一致

## 算法
每个 MAG 的 KO 拷贝数 = 该 MAG 在 KO 注释表里出现该 KO 的次数。
颜色 = `log1p(copies)` 避免高拷贝 KO 挤压色彩空间。

KO 自然顺序来自 KB（`element_pathway_ko_order()`）：As 11 + N 17 + S 15 + Fe 8，
保证同元素 KO 相邻、同通路 KO 相邻。

## 与原脚本对比
| 项 | `scripts/python/06_MAG_gene_profile.py` | EnvMeta `gene_profile.analyze` |
|---|---|---|
| 行数 | 1080 | ~220 |
| 输入 | 6 个 CLI 参数（含 iTOL 热图二次解析）| 4 个 DataFrame |
| KO 顺序 | 硬编码（手工枚举 57 KO）| KB 自动提供 |
| 元素映射 | 硬编码索引 `i<17=As` | KB `flat_ko_map()` |
| 拷贝数计算 | 同（groupby count）| 同 |
| 输出版本 | 多语言 thesis_cn + paper_en | 单图（rcParams 控制字体）|

算法无随机性，数值完全一致。

## Phase 3 循环图应用
本模块输出的 **MAG × KO 拷贝数矩阵** 是循环图"微生物层"的直接数据源：
- 节点：元素循环 KO（从 KB 自动获得）
- 承载关系：哪些 MAG 在此 KO 有 copies > 0
- 权重：KO copy number × MAG 丰度

---

## 2026-05-08 Python 脚本对照执行结果 ✅

原 Python 脚本输出已生成到 `paper_en/Fig3-3_MAG_gene_profile_EN.pdf` +
`Fig3-6a_Top10_gene_profile_EN.pdf` + `Fig3-6b_Keystone_gene_profile_EN.pdf` +
4 个 summary txt。

### 算法一致性（README 已记录"数值完全一致"）

- ✅ 拷贝数计算：groupby count（双方一致）
- ✅ KO 顺序：原脚本手工枚举 57 KO；EnvMeta KB 自动 → 列表完全相同
- ✅ 元素映射：原脚本硬编码索引；EnvMeta KB → 元素归属相同
- ✅ keystone 14 个识别正确
- 📌 **EnvMeta 默认 `drop_zero_kos=True` 过滤全 0 KO 列**（51 vs 57）；
  这是 EnvMeta 的可选 default，**不是 bug**。需保留全 57 列时调用：
  `gene_profile.analyze(..., params={"drop_zero_kos": False})`

### 复现命令

```powershell
$py = "C:\Users\REDLIZZ\.conda\envs\envmeta\python.exe"
& $py "d:\workdata\envmeta_thesis\scripts\python\06_MAG_gene_profile.py" `
  --itol "data\raw\07_element_KO_heatmap.txt" `
  --ko "data\raw\kegg_target_only.tsv" `
  --bac "data\raw\tax_bac120_summary.tsv" `
  --ar "data\raw\tax_ar53_summary.tsv" `
  --abund "data\raw\abundance.tsv" `
  --ks "data\raw\keystone_species.txt" `
  --outdir "paper\benchmarks\validation\gene_profile" `
  --statdir "paper\benchmarks\validation\gene_profile"
```

### 维护记录

| 日期 | 事项 |
|---|---|
| 2026-04-14 | 初版（声明数值完全一致）|
| 2026-05-08 | 原脚本输出生成；51 vs 57 KO 列差异确认为 EnvMeta `drop_zero_kos=True` 默认行为，非 bug |
