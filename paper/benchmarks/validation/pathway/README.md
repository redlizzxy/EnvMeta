# MAG 通路完整度验证

## 输入
- `tests/sample_data/kegg_target_only.tsv`（MAG × KO 长表注释）
- `tests/sample_data/mag_taxonomy_labels.tsv`（MAG + GTDB classification）
- `tests/sample_data/keystone_species.txt`（14 keystone）
- `tests/sample_data/abundance.tsv`（MAG × sample 丰度）

## EnvMeta 输出
- `envmeta_pathway_heatmap.pdf` — 168 MAG × 18 通路完整度热图（元素彩条 + 门彩条 + keystone ★）
- `envmeta_pathway_bubble_top30.pdf` — Top30 丰度 MAG 气泡图（颜色=完整度，大小=丰度）
- `envmeta_pathway_stats.tsv` — 完整度矩阵 + 元数据

## 算法
每个 MAG 的通路完整度 = `|MAG 的 KO ∩ 通路的 KO| / |通路的 KO| × 100%`

通路定义源：`envmeta/geocycle/knowledge_base/elements.json`（4 元素 × 18 通路 × 57 KO）

与原 R 脚本 `scripts/python/08_pathway_completeness.py` 定义完全一致。

## 关键结果
- 168 MAG 全部可分析
- 14 keystone 识别正确（从 `keystone_species.txt` 输入自动标记）
- Top-3 总完整度：
  - `Mx_All_63` Pseudomonadota（1142%，即 18 条通路平均每条 ~63%）
  - `Mx_All_48` Pseudomonadota（1042%）
  - `Mx_All_76` Pseudomonadota（1042%）

**Phase 3 循环图推断引擎**可直接消费这份完整度矩阵：通路完整度 > 50%
的 MAG 视为"该通路活跃承载者"。

## 与原脚本对比
| 项 | `08_pathway_completeness.py` | EnvMeta `pathway.analyze` |
|---|---|---|
| 通路定义 | 硬编码 17 条通路 | 从 KB（`elements.json`）动态读，18 条 |
| 元素映射 | 硬编码 `PATHWAY_ELEMENT` | KB 自动提供 `pathway_element_map()` |
| 输入 | 6 个 CLI 参数 | 4 个 DataFrame（taxonomy/keystone/abundance 可选）|
| 输出 | 3 张独立图 + 3 个 txt | 1 张图（style 参数切换 heatmap/bubble）+ 单个长表 |
| element 过滤 | 无 | `element_filter=['arsenic',...]` |
| MAG 过滤 | 硬编码 Top10/Top30 | `max_mags` + `sort_by` 参数 |

算法无随机性，和原脚本完整度数值应完全一致。
