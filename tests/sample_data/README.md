# 测试样例数据

来自砷渣-钢渣修复宏基因组研究（CK/A/B 三组），已精简以保持仓库轻量。

## 文件索引

| 文件 | 用途 | 适用模块 |
|------|------|---------|
| `metadata.txt` | 样本分组（CK/A/B） | 所有模块 |
| `Phylum.txt` | 门水平丰度表 | 堆叠图 |
| `Genus.txt` | 属水平丰度表 | 堆叠图、LEfSe |
| `Species.txt` | 种水平丰度表 | 堆叠图、RDA |
| `alpha.txt` | α多样性指数（Shannon/Simpson/Chao1） | α多样性 |
| `beta_bray.txt` | Bray-Curtis 距离矩阵 | PCoA |
| `env_factors.txt` | 环境因子（理化指标） | RDA、理化图 |
| `07_element_KO_heatmap.txt` | 元素循环 KO 丰度表 | 基因热图、log2FC |
| `kegg_target_only.tsv` | 全量 KEGG KO 注释 | 通路完整度、循环图 |
| `quality_report.tsv` | CheckM MAG 质量报告 | MAG 质量图 |
| `abundance.tsv` | MAG 丰度矩阵 | MAG 丰度热图 |
| `mag_taxonomy_labels.tsv` | MAG 分类注释 | MAG 丰度热图、基因谱 |
| `gephi_nodes.csv` | 共现网络节点 | 网络图 |
| `gephi_edges.csv` | 共现网络边 | 网络图 |
| `keystone_species.txt` | Keystone 物种列表 | 网络图、循环图 |

## 实验设计

| 组 | 样本 | 钢渣配比 | 配色 |
|----|------|---------|------|
| CK | 2_1, 2_5, 2_7 | 0 | `#4DAF4A` 绿 |
| A  | 8_1, 8_2, 8_3, 8_4 | 低 | `#377EB8` 蓝 |
| B  | 9_1, 9_2, 9_3 | 高 | `#E41A1C` 红 |
