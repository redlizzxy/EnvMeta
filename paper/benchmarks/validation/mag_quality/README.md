# MAG 质量散点图验证

## 输入
- `tests/sample_data/quality_report.tsv`（CheckM2 输出，168 个 MAG）
- `tests/sample_data/mag_taxonomy_labels.tsv`（MAG + GTDB classification）
- `tests/sample_data/keystone_species.txt`（14 个 keystone MAG）

## EnvMeta 输出
- `envmeta_mag_quality.pdf`
- `envmeta_mag_quality_stats.tsv`（summary 3 行 + 168 个 MAG 明细）

## 质量分布
| 质量 | 数量 | Completeness 均值 | Contamination 均值 | Genome Size 均值 (Mb) |
|---|---|---|---|---|
| High（≥90% / <5%）| 35 | 95.4% | 2.33% | 3.67 |
| Medium（≥50% / <10%）| 111 | 76.2% | 4.60% | 3.33 |
| Low | 22 | 62.0% | 8.62% | 3.66 |

168 个 MAG 中 14 个被标为 keystone 物种（菱形黑边高亮 + 标签）。

## 与原脚本对比
| 项 | `scripts/python/06_MAG_quality.py` | EnvMeta `mag_quality.analyze` |
|---|---|---|
| 输入 | 4 个 CLI 参数（checkm2 + bac120 + ar53 + keystone）| 3 个 DataFrame（quality + taxonomy + keystone，后两者可选）|
| GTDB 处理 | 内部合并 bac120 + ar53 | 交给用户（外部可合并为一个 taxonomy 表）|
| 门归类 | 少于 3 个的归 Other + 不在预定义颜色中归 Other | 同，参数 `min_phylum_count` 可调 |
| 质量阈值 | 写死 90/5, 50/10 | 参数 `high_completeness` / `high_contamination` / `med_completeness` / `med_contamination` 可调 |
| 输出 stats | 两个 txt（summary + detail）| 单个长表 `stats_df`（type=summary / detail）|

算法无随机性，数值应与原脚本完全一致。

---

## 2026-05-08 Python 脚本对照执行结果 ✅

原 Python 脚本输出已生成在本目录（`Fig3-1_MAG_quality.pdf` + `MAG_quality_summary.txt`）。

### 数值对照（完全一致）

| Quality | Count (orig) | Count (EnvMeta) | Completeness mean (orig) | Completeness mean (EnvMeta) | Contamination mean (orig) | Contamination mean (EnvMeta) | Genome Size Mb (orig) | Genome Size Mb (EnvMeta) |
|---|---|---|---|---|---|---|---|---|
| High | 35 | 35.0 | 95.4 | 95.38 | 2.33 | 2.33 | 3.7 | 3.67 |
| Medium | 111 | 111.0 | 76.2 | 76.24 | 4.60 | 4.60 | 3.3 | 3.33 |
| Low | 22 | 22.0 | 62.0 | 62.04 | 8.62 | 8.62 | 3.7 | 3.66 |

⭐ **三档分类计数完全一致；指标均值差异 ≤ 0.04%（仅小数四舍五入差异）**。

### 复现命令

```powershell
$py = "C:\Users\REDLIZZ\.conda\envs\envmeta\python.exe"
& $py "d:\workdata\envmeta_thesis\scripts\python\06_MAG_quality.py" `
  --checkm2 "data\raw\quality_report.tsv" `
  --bac120 "data\raw\tax_bac120_summary.tsv" `
  --ar53 "data\raw\tax_ar53_summary.tsv" `
  --keystone "data\raw\keystone_species.txt" `
  --output "paper\benchmarks\validation\mag_quality" `
  --stat_dir "paper\benchmarks\validation\mag_quality" --style paper
```

### 论文引用模板

> EnvMeta MAG quality assessment outputs were validated against the
> original CheckM2-based Python pipeline. Quality tier counts (35 High /
> 111 Medium / 22 Low at thresholds 90/5, 50/10) and per-tier mean
> completeness/contamination/genome-size matched exactly.

### 维护记录

| 日期 | 事项 |
|---|---|
| 2026-04-14 | 初版（仅 EnvMeta 输出）|
| 2026-05-08 | 原 Python 脚本输出 + 数值对照完成 — 完全一致 |
