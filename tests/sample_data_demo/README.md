# Sample data — demo subset (30 MAGs)

> **本目录是在线 demo 用的轻量化测试样本**，由 `tests/sample_data/` 完整 168 MAG
> 数据集抽取 30 MAG 而成。
>
> 抽样规则：(1) 保留全部 keystone MAGs；(2) 保证 4 元素（As/N/S/Fe）每元素 ≥ 5
> 个活跃 MAG；(3) 按总丰度补足到 30 个。

## 用途

- **在线 demo**（Streamlit Cloud）：app.py 默认从这里加载示例数据，减少
  并发崩溃风险（30 MAG 内存占用 ≈ 168 MAG 的 1/5-1/6）。
- **本目录不用于**：单元测试（pytest 跑 `tests/sample_data/`）/ 论文 Arm A 校准
  实验 / 性能 benchmark — 这些场景都用完整 168 MAG 数据集。

## 仅作功能演示，不是科学论断的数据

本目录数据经过 MAG 子采样，**不反映原研究的科学结论**。仅供用户在线点几下
"加载示例数据 → 生成图" 看看 EnvMeta 各功能 UI 是否正常。要复现砷渣修复研究
的结果，请用完整 168 MAG 数据集（`tests/sample_data/`）或下载论文 Fork Bundle。

## 生成方式

运行：

```
python tests/sample_data_demo/_build_demo_subset.py
```

会从 `tests/sample_data/` 重新抽取并覆盖本目录所有文件。

## 文件清单（与 tests/sample_data/ 对应）

- MAG-keyed（已抽样）：abundance.tsv / quality_report.tsv /
  mag_taxonomy_labels.tsv / kegg_target_only.tsv / keystone_species.txt /
  gephi_nodes.csv / gephi_edges.csv
- Sample-keyed（原样拷贝）：metadata.txt / env_factors.txt / alpha.txt /
  beta_bray.txt / ko_tpm.spf / 07_element_KO_heatmap.txt
- Taxon-aggregated（原样拷贝）：Phylum.txt / Genus.txt / Species.txt
- 假说 YAML 模板（原样拷贝）：sample_hypothesis.yaml

## 被抽样的 30 个 MAG

- Mx_All_106
- Mx_All_110
- Mx_All_141
- Mx_All_150
- Mx_All_151
- Mx_All_152
- Mx_All_153
- Mx_All_168
- Mx_All_173
- Mx_All_32
- Mx_All_35
- Mx_All_53
- Mx_All_90
- Mx_All_98
- unmapped
- Mx_All_102
- Mx_All_6
- Mx_All_94
- Mx_All_66
- Mx_All_49
- Mx_All_128
- Mx_All_82
- Mx_All_88
- Mx_All_72
- Mx_All_36
- Mx_All_112
- Mx_All_116
- Mx_All_127
- Mx_All_45
- Mx_All_34
