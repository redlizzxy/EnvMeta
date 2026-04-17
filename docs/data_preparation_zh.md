# EnvMeta 数据准备指南（中文）

> 目标读者：**拿到测序数据但不知道哪个文件能塞进 EnvMeta**的环境微生物研究生。
> 本文回答「我手上的文件 → EnvMeta 哪个输入口」这条路径。

## 0. 3 分钟快速导航

| 你的情景 | 直接看哪节 |
|---|---|
| 测序公司给了压缩包，不知道哪个能用 | §1 测序公司包解压指南 |
| 自己跑了 CoverM / HUMAnN / eggNOG-mapper 等 | §2 上游工具 → EnvMeta 映射表 |
| 想知道 EnvMeta 每种输入具体长什么样 | §3 EnvMeta 标准格式样板 |
| 文件就是读不进去 / 格式报错 | §4 常见问题 FAQ |

---

## 1. 测序公司包解压指南

**国内常见公司**（2026 年行情）：MajorBio（美吉）/ 诺禾致源 / BGI / 基迪奥 /
Wekemo / Biomarker（百迈客）。宏基因组基础套餐结果通常是一个 `xxx_Analysis.zip`
（2-10 GB），里面按流程分目录。

### 1.1 典型目录结构（以 MajorBio 为例）

```
MetagenomeAnalysis/
├── 01.QC/                    # 质控原始 fastq（EnvMeta 用不上）
├── 02.Assembly/              # megahit 组装结果（.fa）
├── 03.GenePrediction/        # prodigal 基因预测
├── 04.Annotation/
│   ├── KEGG/                 # ⭐ KO 注释长表 → KO_ANNOTATION_LONG
│   ├── COG/
│   └── CAZy/
├── 05.Binning/
│   ├── MAGs/                 # MAG fasta
│   ├── CheckM2/              # ⭐ quality_report.tsv → CHECKM_QUALITY
│   └── GTDB-Tk/              # ⭐ summary.tsv → MAG_TAXONOMY
├── 06.Abundance/
│   ├── CoverM/               # ⭐ abundance.tsv → ABUNDANCE_WIDE
│   └── MetaPhlAn4/           # 备选物种丰度
├── 07.Diversity/             # ⭐ alpha / beta → ALPHA_DIVERSITY / DISTANCE_MATRIX
└── 08.Statistics/
    └── LEfSe/                # 备选（EnvMeta 内置，不需要）
```

**⭐ = EnvMeta 直接用的文件**。其他文件虽然测序公司给了，EnvMeta 不依赖。

### 1.2 不同公司的命名差异

| 标准分析 | MajorBio 叫法 | 诺禾叫法 | BGI 叫法 |
|---|---|---|---|
| KO 注释长表 | `all.gene.kegg.tsv` | `ko.annotation.xls` | `KEGG_annotation.txt` |
| MAG 质量表 | `quality_report.tsv` | `bin_stats.tsv` | `mag.quality.csv` |
| MAG 分类 | `gtdbtk.bac120.summary.tsv` | `gtdb.classification.tsv` | `mag_taxonomy.tsv` |
| 丰度矩阵 | `coverm_relative_abundance.tsv` | `taxa.abundance.xls` | `abundance.xls` |

命名虽然不同，**内容字段高度一致**：拿到后直接上传 EnvMeta「文件管理」页，
置信度 > 0.7 的就是识别对了。

### 1.3 公司常「漏给」而 EnvMeta 需要的文件

- **keystone species 列表**：共现网络分析的关键输入。一般需要用 iCAMP /
  SpiecEasi / NetCoMi 从样本丰度矩阵**自己算**，或用 Cytoscape / Gephi 的"网络
  中心性"功能手动筛（Degree ≥ 10 且 Betweenness ≥ 200 是论文常用阈值）。
- **环境因子表**：需要自己整理（采样时记录的 pH / Eh / Total_C / Total_N 等实验室
  测定值）。每列一个因子，每行一个样本。

---

## 2. 上游工具 → EnvMeta 映射表

### 2.1 完整映射矩阵（按工具字母序）

| 上游工具 | 典型输出文件 | EnvMeta FileType | 是否需预处理 |
|---|---|---|---|
| **Bowtie2 / BWA** + samtools | sorted BAM → `coverage.tsv` | `ABUNDANCE_WIDE` | 提取覆盖度矩阵 |
| **CheckM2** | `quality_report.tsv` | `CHECKM_QUALITY` | ❌ 无需 |
| **CoverM** | `abundance.tsv` | `ABUNDANCE_WIDE` | ❌ 无需 |
| **DRAM** | `annotations.tsv` 或 `distill/metabolism.tsv` | `KO_ANNOTATION_LONG` | ❌ 无需（DRAM 已含 KO 列） |
| **eggNOG-mapper** | `*.emapper.annotations` | `KO_ANNOTATION_LONG` | ⚠️ 取 `KEGG_ko` 列，去掉 "ko:" 前缀 |
| **GTDB-Tk** | `gtdbtk.bac120.summary.tsv` | `MAG_TAXONOMY` | ❌ 无需 |
| **HUMAnN3** | `genefamilies.tsv` | `KO_ABUNDANCE_WIDE` | ⚠️ `humann_regroup_table --groups uniref90_ko` |
| **KofamScan** | `*.mapper.tsv` | `KO_ANNOTATION_LONG` | ❌ 无需 |
| **Kraken2 + Bracken** | `*_bracken.tsv` | `ABUNDANCE_WIDE` | ⚠️ 合并多样本为宽表 |
| **MetaPhlAn4** | `*_bugs_list.tsv` | `ABUNDANCE_WIDE` | ❌ 无需 |
| **QIIME2** (diversity alpha) | `shannon_entropy.tsv` | `ALPHA_DIVERSITY` | ⚠️ 多指标合并一张表 |
| **QIIME2** (diversity beta) | `*-distance-matrix.tsv` | `DISTANCE_MATRIX` | ❌ 无需 |
| **iCAMP / SpiecEasi** | `network.edges.tsv` | `GEPHI_EDGES` | ⚠️ 先在 Gephi/R 算 Degree+Betweenness → `GEPHI_NODES` |

### 2.2 工具详解：最常见 5 个

#### CoverM（MAG/物种丰度计算）

```bash
# 标准命令
coverm genome -d mags/ -1 R1.fq.gz -2 R2.fq.gz \
    --methods relative_abundance --min-read-percent-identity 0.95 \
    -o abundance.tsv
```

输出**直接可用**：首列 MAG ID（或 Taxonomy），后续列每列一个样本。

#### HUMAnN3（KO 通路丰度）

```bash
# Step 1: 运行 HUMAnN3
humann --input input.fq.gz --output output/

# Step 2: 把 UniRef90 聚合成 KO（EnvMeta 要 KO）
humann_regroup_table --input output/input_genefamilies.tsv \
    --groups uniref90_ko --output input_ko_abundance.tsv

# Step 3: 多样本合并成宽表
humann_join_tables --input output/ --output ko_abundance_all.tsv \
    --file_name ko_abundance
```

最终 `ko_abundance_all.tsv` 首列 KEGG_ko（K00001 格式），其他列样本 → 对应
EnvMeta `KO_ABUNDANCE_WIDE`。

#### eggNOG-mapper（MAG 基因 KO 注释）

```bash
emapper.py -i predicted_genes.faa --output mag_annotation \
    --cpu 16 --dbmem
```

输出 `mag_annotation.emapper.annotations` 有 20+ 列。**只保留两列**：

```bash
# 提取 MAG + KEGG_ko（去掉 "ko:" 前缀）
awk -F'\t' 'NR>4 && $12!="-" {
    split($12, a, ","); for(i in a) {
        gsub("ko:", "", a[i]); print $1 "\t" a[i]
    }
}' mag_annotation.emapper.annotations > ko_annotation_long.tsv
```

这个 `ko_annotation_long.tsv` 是 EnvMeta `KO_ANNOTATION_LONG`（长表格式）。

> **注意**：EnvMeta 要的是 **MAG × KO 长表**（每行 1 MAG + 1 KO），不是宽表。
> 如果 eggNOG 输出是基因级别，需要先按 MAG 聚合（如果 gene_id 包含 MAG 前缀，
> 提取前缀）。

#### GTDB-Tk（MAG 分类）

```bash
gtdbtk classify_wf --cpus 32 --extension fa \
    --genome_dir mags/ --out_dir gtdb_out/
```

输出 `gtdb_out/gtdbtk.bac120.summary.tsv` 含 `user_genome` + `classification`
（`d__Bacteria;p__Pseudomonadota;c__...;o__...;f__...;g__...;s__...`）→ 直接对应
`MAG_TAXONOMY`。

#### CheckM2（MAG 质量）

```bash
checkm2 predict --threads 16 --input mags/ --output-directory checkm2_out/
```

输出 `checkm2_out/quality_report.tsv` 含 `Name / Completeness / Contamination /
Genome_Size` → 直接对应 `CHECKM_QUALITY`。

### 2.3 共现网络两文件（特殊说明）

EnvMeta 的共现网络分析是 **Gephi 辅助工具**，需要两个 Gephi 格式 CSV：

**gephi_nodes.csv**：
```
Id,Label,Degree,Betweenness,Module,Phylum
Mx_All_63,,18,342,1,Pseudomonadota
Mx_All_102,,24,456,2,Chloroflexota
...
```

**gephi_edges.csv**：
```
Source,Target,Weight,Type
Mx_All_63,Mx_All_102,0.82,Undirected
...
```

**如何生成**：在 R 里用 SpiecEasi 或 iCAMP 做相关性网络（|Spearman r| > 0.9
& p < 0.05），把节点和边的表格导出；Degree / Betweenness 用 Gephi 的"统计"
面板计算（或 igraph `degree()` / `betweenness()`）。

---

## 3. EnvMeta 标准格式样板

### 3.1 `METADATA` 样本分组表

```
SampleID    Group    Replicate
CK_1    CK    1
CK_2    CK    2
CK_3    CK    3
A_1    A    1
A_2    A    2
A_3    A    3
B_1    B    1
B_2    B    2
B_3    B    3
```

**必需列**：`SampleID` + `Group`。其他列（Replicate / 时间 / pH / 等）都是
bonus，不影响识别。

### 3.2 `ABUNDANCE_WIDE` 物种丰度（宽表）

```
Taxonomy    CK_1    CK_2    CK_3    A_1    A_2    A_3    B_1    B_2    B_3
k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;...    0.123    0.145    0.098    ...
k__Bacteria;p__Chloroflexota;c__...    0.056    0.078    0.062    ...
```

**规则**：
- 首列是分类信息（`k__Domain;p__Phylum;...;s__Species`）**或** MAG ID
- 其他列每列一个样本（**列名要和 METADATA 的 SampleID 完全匹配**）
- 数值可以是**相对丰度 %** 或 **TPM** 或 **读长计数**（EnvMeta 会自动归一化）

### 3.3 `DISTANCE_MATRIX` 距离矩阵

```
    CK_1    CK_2    CK_3    A_1    ...
CK_1    0.000    0.234    0.256    ...
CK_2    0.234    0.000    0.198    ...
CK_3    0.256    0.198    0.000    ...
A_1    ...    ...    ...    0.000
```

对称方阵，对角为 0。来源可以是 Bray-Curtis / Unifrac / Jaccard 等。

### 3.4 `ALPHA_DIVERSITY` α 多样性指数表

```
SampleID    Shannon    Simpson    Chao1    Observed_species
CK_1    5.21    0.974    156    128
CK_2    5.18    0.972    151    124
...
```

**识别关键列名**：`shannon` / `simpson` / `chao1` / `observed_species` / `invsimpson`
/ `ace` / `richness` / `evenness`（大小写不敏感）。至少需要一个。

### 3.5 `CHECKM_QUALITY` CheckM 质量表

```
Name    Completeness    Contamination    Genome_Size
Mx_All_1    95.4    1.2    4521000
Mx_All_2    87.3    3.1    3876000
...
```

**识别关键列**：`Completeness` + `Contamination`（CheckM / CheckM2 都用这两列）。

### 3.6 `ENV_FACTORS` 环境因子表

```
SampleID    Group    pH    Eh    Total_C    Total_N    Total_As
CK_1    CK    6.8    124    12.3    1.45    2.3
...
```

**识别关键**：`SampleID` + `Group` + ≥ 2 个数值环境因子列。

### 3.7 `KO_ABUNDANCE_WIDE` KO 丰度宽表

```
KEGG_ko    CK_1    CK_2    CK_3    A_1    ...
K00001    0.00123    0.00134    ...
K00002    0.00045    ...
...
```

**识别关键**：首列符合 `K##### ` 五位数字模式 + 多个数值样本列。

### 3.8 `KO_ANNOTATION_LONG` KO 注释长表

```
MAG    KEGG_ko
Mx_All_1    K00001
Mx_All_1    K00002
Mx_All_1    K00023
Mx_All_2    K00001
...
```

**识别关键**：必须有 MAG 相关列名（`MAG` / `Genome` / `Name` / `Bin` 等）+
KEGG_ko 列。

### 3.9 `KEYSTONE_SPECIES` keystone 物种列表

```
MAG    Degree    Betweenness    Module
Mx_All_63    18    342    1
Mx_All_102    24    456    2
...
```

**识别关键**：MAG 列 + `Degree` 或 `Betweenness` 中的至少一列。

### 3.10 `MAG_TAXONOMY` MAG 分类表

**格式 1**（GTDB-Tk 原生）：
```
user_genome    classification
Mx_All_1    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__...;f__...;g__...;s__...
Mx_All_2    d__Bacteria;p__Chloroflexota;c__...
```

**格式 2**（简化两列无表头）：
```
Mx_All_1    d__Bacteria;p__Pseudomonadota;...
Mx_All_2    d__Bacteria;p__Chloroflexota;...
```

EnvMeta 两种都支持（2 列无表头会自动补 `MAG / Taxonomy` 列名）。

### 3.11 `GEPHI_NODES` / `GEPHI_EDGES` Gephi 网络

见 §2.3。

---

## 4. 常见问题 FAQ

### Q1. 编码乱码？GBK 还是 UTF-8？

EnvMeta 自动嗅探编码（chardet），`.txt` / `.tsv` 优先 UTF-8，GBK 也能读。
如果出现乱码：

```bash
# Windows 检查
file -i xxx.tsv
# 强制转 UTF-8
iconv -f GBK -t UTF-8 xxx.tsv > xxx_utf8.tsv
```

### Q2. 分隔符是 Tab 还是逗号？

EnvMeta 自动检测（比较首行 `\t` 和 `,` 的出现数）。两种都支持，但 **不要混用**
（一个文件内 Tab 和逗号都出现会报错）。推荐 Tab（`.tsv`）因为物种名里常含逗号。

### Q3. 样本 ID 跨文件必须完全一致？

**必须**。如果 metadata 里写 `CK_1`，丰度表列名就必须是 `CK_1`，不能是 `CK1`
或 `CK-1`。常见坑：

- Excel 把 `CK_1` 读成 `CK 1`（下划线变空格）
- 某些测序公司给的文件用原始样本号（如 `2_1`）而 metadata 用别名（`CK_1`）
  → 需统一（推荐在 metadata 里加一列 `OriginalID`，用 `SampleID` 列对应别名）

### Q4. 缺失值怎么表示？

EnvMeta 把 `""` / `NA` / `null` / `-` / `nan` / `NaN` 都视为缺失。数值列的缺失
自动按 0 处理（物种丰度）或排除该行（环境因子 PCA）。

### Q5. 上传后文件类型识别错了怎么办？

文件管理页每个文件卡片中间有**「手动修正类型」**下拉框，选对后 EnvMeta 立即
重新分类。常见误判：

- `alpha.txt` 列名只有 `Shannon` 一项 → 有时会被判成 `ABUNDANCE_WIDE`。
  手动改成 `alpha_diversity`
- `quality_report.tsv` 如果改过列名（改成中文）→ 可能识别不到 `Completeness`
  列。建议保留英文原名

### Q6. 文件太大（> 100 MB）跑不动？

EnvMeta 大多数分析在 < 50 MB 文件下 < 10 秒跑完。如果遇到：

- 丰度表超过 500 物种 / 超过 50 样本 → 建议先在本地用 `head -1000` 精简
- KO 注释长表超过 100 万行 → 检查是否没按 MAG 聚合（eggNOG 原始是基因级）
- 文件内含大量 0 值行 → 先过滤掉（`awk '{if(NR>1){s=0;for(i=2;i<=NF;i++)s+=$i; if(s>0) print; else next}else print}'`）

### Q7. Windows 路径里中文 / 空格导致读不出来？

EnvMeta 通过 Streamlit `file_uploader` 上传，不走原始路径，所以路径里有中文
/ 空格**不影响**。但如果你在命令行用 `scripts/` 里的代码直接读文件，用
`Path(r"D:\my path\xxx.tsv")` 或带引号的字符串。

### Q8. 我要分析的文件不在支持列表里？

EnvMeta 支持「**手动修正类型**」+ 下游分析自动适配。只要你的数据结构大致符合
§3 里某个 FileType 的特征，把它手动标成那个类型，就能跑对应分析。

如果你的文件结构完全是新的（如转录组 TPM 矩阵），请在 GitHub 上提 Issue，
我们考虑加规则：https://github.com/redlizzxy/EnvMeta/issues

---

## 5. 附录：scripts/to_gephi.R 共现网络导出示例

```r
library(igraph)
library(Hmisc)

# 输入：样本 × 物种丰度矩阵（行是样本）
abund <- read.table("abundance.tsv", sep="\t", header=TRUE, row.names=1)
abund_t <- t(abund)  # 物种 × 样本

# Spearman 相关
cor_res <- rcorr(t(abund_t), type="spearman")
r_mat <- cor_res$r
p_mat <- cor_res$P

# 筛阈值：|r| > 0.9 且 p < 0.05
r_mat[abs(r_mat) < 0.9 | p_mat > 0.05 | is.na(r_mat)] <- 0
diag(r_mat) <- 0

# 构建 igraph
g <- graph_from_adjacency_matrix(r_mat, mode="undirected", weighted=TRUE)
g <- delete.vertices(g, degree(g) == 0)

# 计算 Degree / Betweenness
V(g)$Degree <- degree(g)
V(g)$Betweenness <- betweenness(g, normalized=FALSE)

# 导出 Gephi 节点 CSV
nodes_df <- data.frame(
    Id = V(g)$name,
    Label = "",  # EnvMeta prepare_gephi_csv 会补
    Degree = V(g)$Degree,
    Betweenness = round(V(g)$Betweenness, 2)
)
write.csv(nodes_df, "gephi_nodes.csv", row.names=FALSE, quote=FALSE)

# 导出 Gephi 边 CSV
edges_df <- as.data.frame(get.edgelist(g))
colnames(edges_df) <- c("Source", "Target")
edges_df$Weight <- E(g)$weight
edges_df$Type <- "Undirected"
write.csv(edges_df, "gephi_edges.csv", row.names=FALSE, quote=FALSE)
```

导出的两个 CSV 直接上传 EnvMeta「共现网络图」页即可。

---

## 结语

如果你从测序公司收到数据包、或自己跑完上游分析后看完这份指南仍不知道哪个文件
能用，欢迎：

- 在 EnvMeta「文件管理」页直接拖上去试 —— 识别不到的文件会标为 `unknown` +
  给出识别失败原因
- 开 Issue：https://github.com/redlizzxy/EnvMeta/issues

**EnvMeta 的设计哲学**：上游工具千千万，下游分析只要格式合法就能跑。不要被
上游差异吓退。
