# mag_heatmap validation (S6)

Real sample data: 168 MAGs × 10 samples (3 CK / 4 A / 3 B), GTDB classification
+ 14 keystone species list.

## Inputs

- `tests/sample_data/abundance.tsv` — MAG × sample 相对丰度 (%)
- `tests/sample_data/mag_taxonomy_labels.tsv` — MAG + GTDB classification
- `tests/sample_data/keystone_species.txt` — 14 keystone MAGs
- `tests/sample_data/metadata.txt` — SampleID + Group

## Outputs

| 文件 | 说明 |
|---|---|
| `top30_EN.pdf` | Top-30 MAG (mean 选择) 热图；门内聚类 + 三段配色 + keystone ★ |
| `top30_stats.tsv` | Top-30 MAG × sample 百分比 + Phylum + is_keystone + row_order |
| `top30_variance_EN.pdf` | Top-30 MAG (variance 选择) — 突出组间差异大的 MAG |
| `top30_variance_stats.tsv` | 同上，variance 分数 |
| `_run.py` | 复现脚本 |

## 结果摘要（mean 模式）

- **Top-30 中 keystone 5/14**（35.7%）—— 高丰度 ∩ 关键物种的交集
- **8 个门** 在 Top-30 中出现；Pseudomonadota (11) / Chloroflexota (8) /
  Acidobacteriota (6) 占主导
- **Top-1 MAG `Mx_All_102`**（Chloroflexota）mean=0.698%，variance 模式
  同样 top-1（分数 0.854）—— 说明它是高丰度 **且** 组间分化大的 MAG，
  值得作为后续机制验证候选

## variance 模式补充

同数据用 `selection_by="variance"` → Top-3: `Mx_All_102 / Mx_All_6 / Mx_All_98`。
variance 前 30 的 keystone 重叠情况与 mean 模式不同，互补使用可识别：
- **mean 独有**：广义"核心菌群"MAG
- **variance 独有**：对处理响应敏感的 MAG（A/B vs CK 差异信号）

## 与原脚本对照

原 `scripts/python/07_MAG_abundance_heatmap.py`：
- 544 行；硬编码 `SAMPLE_ORDER` / `SAMPLE_GROUP`；专为论文 3-组 10-样本设计
- 图例 / 配色 / 分组分隔线、门色带、keystone 菱形 + rank 数字全部内联

本模块 `envmeta/analysis/mag_heatmap.py`：
- ~220 行；无硬编码样本/组；metadata 驱动组彩条顺序
- 保留三段非线性配色核心算法（Blues → YlGn → YlOrRd）
- 保留门内聚类、keystone 前缀 ★、组彩条
- 简化图例（交给 figure colorbar）；后续可加专用图例区
