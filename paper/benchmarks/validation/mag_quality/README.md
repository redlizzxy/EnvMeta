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
