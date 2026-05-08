# 第二外部数据集 — Phase 0 数据可用性检查日志

> **状态**：进行中
> **目的**：在投入大量解码和处理时间前，先确认两个候选数据集的实际可下载性 +
> 文件清单 + EnvMeta 输入兼容性，以便决定走 Phase 1 / Phase 2 哪条路径。
> **关联**：[../../../paper/manuscript/external_validation_plan.md](../../../paper/manuscript/external_validation_plan.md)

---

## 1. MUCC v2 永冻土（Wilson et al. 2024 *Nat Microbiol*）

### 数据源
- Zenodo: <https://zenodo.org/records/14532347>
- DOI: 10.5281/zenodo.14532347（待确认）

### 下载状态
- [ ] zip 已下载到本地（路径：）
- [ ] 已解压
- [ ] 文件清单已抄录（见下）

### 文件清单（解压后实测）

> 用户解压后填入：

```
[此处粘贴 ls -la 输出，列出所有文件名 + 大小]
```

### EnvMeta 兼容性核对（Phase 0 关键点）

需确认 5 项（任一缺失则路径需调整）：

| 必需输入 | 是否存在 | 文件名 | 备注 |
|---|---|---|---|
| MAG 丰度表（MAG × sample） | ⬜ | | |
| KO 注释（MAG × KO 或 sample × KO） | ⬜ | | |
| 样本 metadata（生境 + 深度） | ⬜ | | |
| GTDB-Tk taxonomy | ⬜ | | |
| 环境因子（pH / 温度 / 深度，可选） | ⬜ | | |
| MAG 质量（CheckM completeness/contamination） | ⬜ | | |

### 决策

- ✅ **5 项齐全** → Phase 1 直接 reshape + 跑图（最优路径）
- ⚠️ **缺 KO 注释但有 MAG fasta** → 需要本地跑 KofamScan，加 1-2 天
- ❌ **缺 MAG 丰度** → 需要从 raw reads 重新 mapping，加 3-5 天
- ❌❌ **缺多项** → 退备选数据集（Tara Oceans）

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
