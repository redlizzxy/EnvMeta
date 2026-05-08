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

### Suppl Tables 清单（实测）

> 用户翻完 Suppl 后填入：

```
[此处列每张 Suppl 表的标题 + 内容简述]
例：
Table S1 — Sample metadata (n=36, pH/EC/As 等)
Table S2 — MAG list (高质量 20 个，含 GTDB-Tk taxonomy)
Table S3 — ROCker hits (As-related KO/gene)
Table S4 — KO matrix (MAG × KO)？？？  ← 关键，需确认
```

### EnvMeta 兼容性核对

| 必需输入 | Suppl 中是否提供 | 备注 |
|---|---|---|
| 样本 metadata + pH | ⬜ | |
| MAG 丰度（MAG × sample） | ⬜ | |
| KO 注释（MAG × KO） | ⬜ | **决定 A/B/C 路径的关键** |
| GTDB-Tk taxonomy | ⬜ | |

### 决策（A/B/C 路径）

- **路径 A**：Suppl 含 KO 矩阵 + MAG 丰度 → 直接 reshape，1-2 天
- **路径 B**：Suppl 缺 KO 矩阵但有 MAG fasta → ENA 下 fasta + 跑 KofamScan，3-5 天
- **路径 C**：Suppl 仅有 ROCker hits（不含 full KO）→ 子采样 8 样本重跑全流程，1-2 周
- **路径 D（搁置）**：Suppl 完全不够 → 暂搁置 Wei，只做 MUCC

实际决策：⬜（待填）

---

## 3. Phase 0 完成后的下一步

填完上述两份清单后，根据决策结果选节奏：

| MUCC | Wei | 节奏 | 预计完成 |
|---|---|---|---|
| ✅ 路径 1 | A | **节奏 2 激进**：MUCC + Wei 双管齐下 | 5-7 天 |
| ✅ 路径 1 | B | 节奏 1.5：MUCC 优先 + Wei 中等耗时 | 7-10 天 |
| ✅ 路径 1 | C/D | **节奏 1 保守**：仅 MUCC | 3-5 天 |
| ❌ MUCC 退备选 | — | 重新评估，可能换 Tara Oceans | TBD |

---

## 维护记录

| 日期 | 事项 | 操作人 |
|---|---|---|
| 2026-05-08 | 模板创建，Phase 0 启动 | Claude |
