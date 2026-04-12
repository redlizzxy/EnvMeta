# 钢渣-脱硫石膏修复砷选冶渣的微生物群落宏基因组研究

## Metagenomic Study of Microbial Communities in Arsenic Smelting Slag Remediated with Steel Slag–Desulfurization Gypsum

---

## 项目结构

```
thesis_arsenic_metagenomics/
│
├── config/                          # 统一配置
│   ├── plot_config.py               #   Python绘图配置（配色/字体/路径/尺寸）
│   ├── theme_thesis.R               #   R绘图主题配置
│   └── __init__.py
│
├── scripts/                         # 分析脚本
│   ├── python/                      #   Python脚本（按图表编号排序）
│   │   ├── 01_physicochemical.py    #     Fig1-2 理化指标柱状图
│   │   ├── 02_alpha_diversity.py    #     Fig2-4 α多样性箱线图
│   │   ├── 03_gene_heatmap.py       #     Fig2-8 元素循环基因热图
│   │   ├── 04_gene_log2fc.py        #     Fig2-9 基因log2FC对比图
│   │   ├── 05_MAG_quality.py        #     Fig3-1 MAG质量散点图
│   │   ├── 06_MAG_heatmap.py        #     Fig3-3 MAG丰度热图
│   │   ├── 07_pathway_completeness.py #   Fig3-4 通路完整度
│   │   └── 08_species_contribution.py #   Fig3-5 贡献物种气泡图
│   │
│   ├── R/                           #   R脚本
│   │   ├── 01_tax_stackplot.R       #     Fig2-1~2-3 物种组成堆叠图
│   │   ├── 02_beta_PCoA.R          #     Fig2-5 PCoA + PERMANOVA
│   │   ├── 03_RDA_Mantel.R         #     Fig2-6 RDA排序图
│   │   ├── 04_LEfSe_analysis.R     #     Fig2-7 LEfSe差异分析
│   │   ├── 05_phylogenetic_tree.R  #     Fig3-2 MAG系统发育树
│   │   └── 06_network_prep.R       #     Fig3-7 网络分析预处理
│   │
│   └── shell/                       #   自动化脚本
│       ├── run_all.sh               #     一键运行全部分析 + 打包
│       ├── sync_to_win.sh           #     快速同步到Windows D盘
│       └── pack_copyright.sh        #     软件著作权代码打包
│
├── data/                            # 数据文件
│   ├── raw/                         #   原始数据
│   │   ├── metadata.txt             #     样本分组信息
│   │   ├── Phylum_fixed.txt         #     门水平丰度表
│   │   ├── Genus_fixed.txt          #     属水平丰度表
│   │   ├── Species_fixed.txt        #     种水平丰度表
│   │   ├── alpha.txt                #     α多样性指数
│   │   ├── beta_bray.txt            #     β多样性距离矩阵
│   │   └── ko_list.txt              #     KEGG KO注释列表
│   │
│   └── processed/                   #   处理后数据（脚本生成）
│       ├── species_percentage.txt   #     物种百分比统计
│       ├── log2fc_results.csv       #     基因差异分析结果
│       └── ...
│
├── figures/                         # 图表输出
│   ├── thesis_cn/                   #   毕业论文版（中文标注，300dpi PDF）
│   ├── paper_en/                    #   小论文版（英文标注，600dpi TIFF/PDF）
│   └── supplementary/               #   补充图
│
├── tables/                          # 表格（CSV/TSV）
│
├── docs/                            # 文档
│   ├── copyright/                   #   软件著作权材料
│   │   ├── 前30页源代码.txt
│   │   ├── 后30页源代码.txt
│   │   └── 代码统计.txt
│   └── 图表对照清单.md             #   图号 ↔ 脚本 ↔ 数据 对照表
│
├── logs/                            # 运行日志
├── requirements.txt                 # Python依赖
├── README.md                        # 本文件
└── LICENSE                          # 许可证
```

## 快速开始

### 1. 环境准备

```bash
# Python 依赖
pip install -r requirements.txt --break-system-packages

# R 依赖
Rscript scripts/R/install_packages.R
```

### 2. 放入数据

将原始数据文件放入 `data/raw/` 目录。

### 3. 运行分析

```bash
# 一键运行全部分析 + 自动打包到 D盘
bash scripts/shell/run_all.sh

# 只跑第二章
bash scripts/shell/run_all.sh --chapter 2

# 只打包（已有结果时用）
bash scripts/shell/run_all.sh --skip-analysis

# 日常快速同步到D盘
bash scripts/shell/sync_to_win.sh
```

### 4. 软著打包

```bash
bash scripts/shell/pack_copyright.sh
```

## 实验设计

| 组别 | 样本 | 钢渣配比 | 说明 |
|------|------|---------|------|
| CK   | 2_1, 2_5, 2_7 | 0 | 对照组 |
| A    | 8_1, 8_2, 8_3, 8_4 | 低 | 低钢渣比 |
| B    | 9_1, 9_2, 9_3 | 高 | 高钢渣比 |

## 图表清单

| 图号 | 描述 | 脚本 | 数据 |
|------|------|------|------|
| Fig1-2 | 理化指标柱状图 | `python/01_physicochemical.py` | metadata.txt |
| Fig2-1~3 | 物种组成堆叠图 | `R/01_tax_stackplot.R` | Phylum/Genus/Species_fixed.txt |
| Fig2-4 | α多样性箱线图 | `python/02_alpha_diversity.py` | alpha.txt |
| Fig2-5 | PCoA分析图 | `R/02_beta_PCoA.R` | beta_bray.txt |
| Fig2-6 | RDA排序图 | `R/03_RDA_Mantel.R` | metadata.txt + Species |
| Fig2-7 | LEfSe差异分析 | `R/04_LEfSe_analysis.R` | Genus_fixed.txt |
| Fig2-8 | 元素循环基因热图 | `python/03_gene_heatmap.py` | ko_list.txt |
| Fig2-9 | 基因log2FC对比 | `python/04_gene_log2fc.py` | ko_list.txt |
| Fig3-1 | MAG质量散点图 | `python/05_MAG_quality.py` | checkm_results.txt |
| Fig3-2 | MAG系统发育树 | `R/05_phylogenetic_tree.R` | tree_metadata.txt |
| Fig3-3 | MAG丰度热图 | `python/06_MAG_heatmap.py` | group_mean_abundance.txt |
| Fig3-7 | 共现网络图 | Gephi + `R/06_network_prep.R` | network_*.txt |
| Fig3-9 | 生物地球化学模型 | AI/PPT手绘 | — |

## 配色方案

- CK组: `#4DAF4A` (绿) — A组: `#377EB8` (蓝) — B组: `#E41A1C` (红)

## 依赖

- Python ≥ 3.8: matplotlib, numpy, pandas, scipy, seaborn
- R ≥ 4.0: ggplot2, vegan, amplicon, ggtree, pheatmap, RColorBrewer
- 其他: Gephi (网络图), Adobe Illustrator/PPT (机制图)

## 作者

[你的姓名] — [学校/导师] — 2026
