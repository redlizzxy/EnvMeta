# 第二外部数据集 — Phase 0 数据可用性检查日志

> **状态**：进行中（2026-05-08 二次调研后切换异方向数据集）

## ⚠️ 2026-05-08 重大决定：切换异方向数据集

### 原计划（已废弃）

原推荐：MUCC v2 Wilson 2024 永冻土，下载 Zenodo 14532347。

### 实测发现（已确认）

Zenodo 14532347 网页 description **宣称**含 4745 MAGs / KO 注释 / 丰度矩阵，但实际 API 返回的 `files[]` **只有 5 个 Methanoregula 单属子集文件，共 62 MB**。description 里宣传的全量数据**根本没上传到这条 Zenodo 记录**，且这条记录与 Wilson 2024 *Nat Microbiol* 不是同一篇论文。

Wilson 2024 论文的真实数据散落在：
- Zenodo 10426238 (MAG fasta + metadata, 12.4 GB)
- Zenodo 7587534 (DRAM 注释 644 MB，**含 KO**)
- Nature Source Data ZIP（**MAG×sample 丰度矩阵可能在这里**，但不一定能挖到）

四件套（MAG / KO / 丰度 / metadata）需要从 4 处拼装，且丰度矩阵**没有**直接发布的处理后版本。

### 新选择：Tara Oceans MAGs (Meren Lab curated) ⭐

| 项 | 内容 |
|---|---|
| 论文 | Delmont et al. 2018 *Nat Microbiol* + 后续多篇 |
| 入口 | <https://merenlab.org/data/tara-oceans-mags/> |
| 数据 | 1888 MAGs + Anvi'o profile DB（含 mapping coverage 矩阵）+ KEGG/COG 注释 + 样本 metadata |
| 大小 | ~15 GB（曾被国内多组镜像过）|
| 主题 | 海洋 surface/DCM/mesopelagic + 64 站点的 N/S 循环 |
| 优势 | **完整四件套 + 文档化最完整 + 中国可稳下 + 海洋 N/S 命中 EnvMeta 4-元素 KB（N+S）+ 论文审稿人最熟** |

为什么改选这个：
- ✅ 唯一**已发布完整 KO 矩阵 + MAG 丰度矩阵**的非砷数据集
- ✅ Anvi'o profile DB 含 mapping coverage（直接用于 MAG 丰度热图）
- ✅ 海洋 N/S 循环 = EnvMeta 4-元素 KB 中 2 个元素的命中
- ✅ Meren Lab 文档化是行业标杆，复现叙事强

### 备选（若 Tara 失败）

**SMAG (Soil MAGs)** — Ma et al. 2023 *Nat Commun* "A genomic catalogue of soil microbiomes"：40,039 MAGs + KEGG + GTDB-Tk + abundance table。中科院南土所发布，国内访问稳定。需切片到生境子集（如 grassland/forest），压到 ~20 GB。


> **目的**：在投入大量解码和处理时间前，先确认两个候选数据集的实际可下载性 +
> 文件清单 + EnvMeta 输入兼容性，以便决定走 Phase 1 / Phase 2 哪条路径。
> **关联**：[../../../paper/manuscript/external_validation_plan.md](../../../paper/manuscript/external_validation_plan.md)

---

## 1. ~~MUCC v2 永冻土~~ → 改用 Tara Oceans MAGs (Meren Lab curated)

### Phase 0 实测：MUCC 已废弃（见上方决定记录）

### 新数据源（Tara Oceans）

- Meren Lab 主页：<https://merenlab.org/data/tara-oceans-mags/>
- 含 MAG fasta、Anvi'o profile DB、KEGG/COG 注释、样本 metadata 全套

### 下载状态
- [ ] Meren Lab 页面已打开浏览
- [ ] 1888 MAG profile DB 下载（约 5 GB）
- [ ] KEGG/COG 注释 tsv 下载
- [ ] 样本 metadata tsv 下载

### EnvMeta 兼容性核对（Tara Oceans）

需确认 5 项：

| 必需输入 | 是否存在 | 文件名 | 备注 |
|---|---|---|---|
| MAG 丰度表（MAG × sample） | ⬜ | | Anvi'o profile DB 含 mapping coverage |
| KO 注释（MAG × KO 或 sample × KO） | ⬜ | | KEGG modules 注释由 anvi-run-kegg-kofams 输出 |
| 样本 metadata（surface/DCM + 站点 + 海洋区） | ⬜ | | |
| GTDB-Tk taxonomy / Meren 自带 taxonomy | ⬜ | | Meren Lab 用 anvi'o 而非 GTDB-Tk，需转换 |
| 环境因子（温度 / 盐度 / nitrate 等，可选） | ⬜ | | Tara Oceans 公开多种 env factor |
| MAG 质量（completeness/contamination） | ⬜ | | Anvi'o stats 输出 |

### 决策

- ✅ **5 项齐全** → Phase 1 直接 reshape + 跑图
- ⚠️ **Anvi'o profile DB 需要装 anvi'o 才能解** → 评估是否值得装（4-6 h 装）
- ❌ **缺 KO 矩阵** → 退备选 SMAG（中科院南土所）

实际决策：⬜（待填）

---

## 2. Wei 2024 砷稻田（*Microbiome*）

### 数据源
- 论文 PDF: <https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-024-01952-4>
- BioProject: <https://www.ebi.ac.uk/ena/browser/view/PRJNA1068274>
- 论文 Suppl Tables: 通过论文页面 "Supplementary Information" 下载

### 下载状态
- [ ] 论文 PDF 已下载
- [ ] Supplementary Tables 已下载（特别留意 Table S2 / S3 / S4 等可能含 KO 矩阵的表）
- [ ] BioProject 样本列表已浏览（确认 36 samples 是否齐全）

### Suppl Tables 清单（2026-05-08 实测）

**MOESM1_ESM.docx**（仅图）：Fig. S1-S18 — 系统发育/地理位置/物化参数对比/亚社区分类等

**MOESM2_ESM.xlsx**（数据表，关键）：

| Table | 内容 | EnvMeta 用途 |
|---|---|---|
| **S1** | Sample metadata（n=36） | metadata.tsv（样本分组：As-contaminated vs noncontaminated）|
| **S2** | Physicochemical parameters of paddy soil | env.tsv（pH / EC / As / NO₃⁻ / NH₄⁺ 等环境因子，喂 RDA + 假说评分）|
| S3 | Metagenomic dataset statistics | （读图统计，跑图非必需）|
| S4 | Metatranscriptomic dataset statistics | （RNA 层，可选）|
| S5 | UniProt IDs for ROCker references | （注释方法学）|
| S6 | As genes-carrying reads in transcriptomes | （RNA 层）|
| **S7** | **179 MAGs with completeness > 50% & contamination < 10% — Taxa + relative abundance** | **mag_abundance.tsv + mag_taxonomy.tsv + checkm.tsv 三合一** ⭐ |
| **S8** | **Functional genes (As ox/red + denitrification + DNRA + methane oxidation) in MAGs** | **ko_long.tsv（精选 KO，5 大功能基因×179 MAG）** ⭐ |
| **S9** | **20 MAGs carrying As(III) ox + denitrification/DNRA gene** | **作者核心假说证据 → 喂 EnvMeta YAML 评分器** ⭐ |
| S10 | Transcriptional activity of functional genes | （可选 RNA 层对照）|

### EnvMeta 兼容性核对（实测）

| 必需输入 | Suppl 中是否提供 | 来源 | 备注 |
|---|---|---|---|
| 样本 metadata + pH 等 | ✅ | S1 + S2 | 36 paddy samples, As-contaminated vs noncontaminated |
| MAG 丰度（MAG × sample） | ✅ | S7 | 179 高质量 MAG (completeness>50%) |
| KO 注释（MAG × KO） | ⚠️ **精选版** | S8 | **不是 KEGG 全 KO**，是论文核心 5 类基因（aioA/arxA/arrA/arsC1/arsC2 + denit + DNRA + mcrA）|
| GTDB-Tk taxonomy | ✅ | S7 | 同行 |
| 环境因子 (env.tsv) | ✅ | S2 | pH / EC / As / NO₃⁻ / NH₄⁺ 等齐全 |
| MAG 质量 | ✅ | S7 | completeness + contamination 列 |

### 决策

✅ **路径 A'（修正版）确立**：

- Wei Suppl 提供"精选 KO"（不是 KEGG 全集），但**正中 EnvMeta 砷+N 循环图卖点**
- 1-2 天可跑通：直接 xlsx → reshape → EnvMeta
- 跑得通的图（≥ 8 张）：
  1. 物种堆叠图（179 MAG taxonomy）
  2. α 多样性（As-cont vs noncont 两组）
  3. β-PCoA + PERMANOVA（基于 MAG 丰度）
  4. RDA（基于 S2 物化参数 vs MAG 丰度，pH 主轴）
  5. LEfSe（As-cont vs noncont 差异 MAG）
  6. 基因热图（精选 KO × MAG）
  7. MAG 丰度热图
  8. **🔑 循环图（As + N 元素，从 S8 精选 KO 推断）— EnvMeta 独家**
  9. **🔑 YAML 假说评分**：直接写"As(III) oxidation 与 denitrification 共发生"假说（S9 是金标准） → 喂评分器 → 复现作者结论
- 跑不通/简化的图（≤ 6 张）：MAG 质量散点（S7 给了 comp/cont 但缺 N50 等）/ 通路完整度（精选 KO 不全 KEGG module）/ MAG 基因谱（精选）/ 共现网络（缺 sample × MAG 共现）

**License 注意**：Wei 论文为 CC BY-NC-ND 4.0。Suppl 原始 xlsx **不入 EnvMeta git repo**，仅本地保留 + reshape 脚本入仓库；EnvMeta 跑出的图作为新生成产物可入仓库 + 论文使用，并 cite Wei 2024。

实际决策：✅ **路径 A'（2026-05-08 确立）**

---

## 3. Phase 0 决策（2026-05-08 收尾）

### 最终路径

| 数据集 | 决策 | 理由 |
|---|---|---|
| **Wei 2024 砷稻田** | ✅ **首选 — 路径 A' 确立** | Suppl Table S1/S2/S7/S8/S9 全提供，1-2 天跑通；正中 EnvMeta 假说评分卖点 |
| **Tara Oceans (Meren)** | ⏸ **次选，待 Wei 完成后视情况启动** | ~15 GB + 装 anvi'o 麻烦；Wei 完成后若叙事已强可省略 |
| ~~MUCC v2~~ | ❌ **废弃** | 数据散落 4 处 + 丰度矩阵无独立发布 |

### 节奏

**Phase 1 — Wei 复现**（2026-05-08 至 2026-05-10，2-3 天）：
- Day 1：reshape Wei MOESM2.xlsx → EnvMeta 6 个 input file
- Day 2：跑 9 张图 + 写 YAML 假说 + 截图归档
- Day 3：写 README + compare_to_original.md + Methods 段落初稿

**Phase 1.5 决策点**（2026-05-10）：
- 如 Wei 复现叙事已强 → 跳 Tara Oceans，进 Methods 全稿
- 如需"跨主题"加强 → 启动 Tara Oceans Phase 1（再加 3-5 天）

---

## 维护记录

| 日期 | 事项 | 操作人 |
|---|---|---|
| 2026-05-08 | 模板创建，Phase 0 启动 | Claude |
