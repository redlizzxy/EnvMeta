# Performance & Feature Benchmark: EnvMeta vs Comparable Tools

> **创建日期**：2026-05-14
> **状态**：Task 3a 方法学设计完成 / Task 3b 实测填表 pending
> **关联**：[outline_gpb.md §5.2.9 / §5.2.10](../../../manuscript/outline_gpb.md)
> **目标**：为 GPB 投稿提供与同类工具的客观性能 + 功能对比（GPB Reviewers 必查项）

---

## 1. Benchmark scope

EnvMeta 提供 14 项分析 + 元素循环图 + 假说评分 + 独立交互 HTML。同类工具
在不同维度有各自强项 / 短板。本 benchmark 分两类对比：

### Type A — 共享能力（可直接对比）

EnvMeta 与某工具均提供同一类分析时，跑同一份数据，比 wall time / RSS peak /
output 维度。

| 分析任务 | EnvMeta | Anvi'o | Krona | MicrobiomeAnalyst | QIIME2 (ref) |
|---|---|---|---|---|---|
| 物种组成堆叠图 | ✅ stackplot | partial (anvi-display-pan) | ❌ | ✅ Taxa visualization | ✅ taxa-bar-plot |
| α-diversity | ✅ alpha_boxplot | partial | ❌ | ✅ Alpha diversity | ✅ alpha-rarefaction |
| β-diversity / PCoA | ✅ pcoa | ❌ | ❌ | ✅ Beta diversity | ✅ beta-diversity |
| KO 热图 | ✅ gene_heatmap | partial | ❌ | partial | ❌ |
| MAG 丰度热图 | ✅ mag_heatmap | ✅ interactive-binning | ❌ | ❌ | ❌ |
| 通路完整度 | ✅ pathway | partial (KEGG metabolism) | ❌ | partial (KEGG enrichment) | ❌ |
| 层级分类饼/扇 | partial (stackplot) | ❌ | ✅ HTML | ✅ Sunburst | ❌ |
| 网络（Gephi prep）| ✅ network | ❌ | ❌ | ✅ Correlation network | ❌ |

### Type B — 独家能力（仅功能矩阵）

EnvMeta 独有的 5 项功能，列功能矩阵显示同类工具均**不提供**。

| 功能 | EnvMeta | Anvi'o | Krona | MicrobiomeAnalyst | iTOL | iMeta sister tools (ImageGP 2 等) |
|---|---|---|---|---|---|---|
| 自动元素循环推断（4 元素 × 18 通路 × 57 KO）| ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| YAML 假说评分（6 claim types + null calibration）| ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| 独立离线交互 HTML SI（含分析本体）| ✅ ~400 KB | partial (在线 web) | partial (HTML 但无分析) | ❌ web-only | ❌ web-only | ❌ |
| 跨元素化学物耦合（orpiment 等）| ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| 假说权重 ±20% 敏感度 + 置换 null calibration | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |

---

## 2. 工具选取与理由

### 2.1 Anvi'o (v8+)

- **理由**：MAG-centric 工具的事实标准；与 EnvMeta 在 MAG 探索 / pangenome
  分析有重叠
- **限制**：主要面向交互式 MAG 探索，不是图表批量生成；与 EnvMeta 的"出版级
  图 + 调参 GUI"workflow 哲学不同
- **可比分析**：MAG quality stats / abundance heatmap / pathway completeness
- **install**：`conda create -n anvio-8 -c conda-forge -c bioconda python=3.10 anvio==8`
  （~5-10 GB 环境）
- **input prep**：需要 contigs.db + profile.db；从 KO/abundance 表起步不直接

### 2.2 Krona (v2.8+)

- **理由**：HTML 输出独立离线，与 EnvMeta 的 HTML 导出有比较点
- **限制**：仅做层级分类饼，不支持 PCoA / 热图 / 网络
- **可比分析**：物种层级 HTML（vs EnvMeta stackplot HTML）
- **install**：`conda install -c bioconda krona`
- **input prep**：tab-delimited TaxID + abundance

### 2.3 MicrobiomeAnalyst (web, v2.0)

- **理由**：用户友好 GUI 同侪；微生物组下游分析 web 平台代表
- **限制**：纯 web，手动 upload + 点击；无 CLI/API；session 不持久
- **可比分析**：Taxa visualization / Alpha-Beta diversity / Heatmap / KEGG enrichment
- **install**：无（web 访问 https://www.microbiomeanalyst.ca/）
- **测量限制**：无法精确测 wall time / RSS（在浏览器中），只能报"上传 → 出图"用户
  操作时长 + 输出维度

### 2.4 QIIME2 (v2024+, 参考对照)

- **理由**：amplicon 标准管线；不是 EnvMeta 直接竞品，作 reference baseline
- **不深度对比**，仅在 Table 提到能力覆盖
- **install**：`conda install -c qiime2 qiime2`

### 2.5 排除的工具

- **iMeta sister tools (ImageGP 2 等)**：web-only / 通用 biomedical / 与 EnvMeta
  scope 重叠少，不在 Type A 直接对比；在 §5.2.9 ecosystem matrix 作 reference
- **phyloseq / R**：CLI / R 包，与 GUI 工具不同 paradigm
- **iTOL**：phylogeny 专精，不覆盖 EnvMeta 主要分析

---

## 3. 测试数据集

统一使用 **Liu 2023 cold seep dataset 子集**作为 benchmark 数据（合理 N_MAG / N_sample 规模）：

- N_MAG = 200（从 1084 抽样，保持元素覆盖）
- N_sample = 30
- 已注释 KO：8 个（published） + 25 synthetic-dense (扩展)
- env factors = 4 numeric

每个工具用**同样的数据子集**（按工具能接受的格式预转换：abundance.tsv + ko_long.tsv +
metadata.tsv），跑相同的 3 个核心分析（PCoA + heatmap + 通路完整度）。

测试硬件：**单台 Intel i7-class 16 GB Win10 笔记本 / Python 3.11**（与 §5.2.8 一致）。

---

## 4. 测量指标

### 4.1 性能指标

| 指标 | 测量方式 |
|---|---|
| Wall time | `time.perf_counter()` 包裹工具入口；3 次重复取 median |
| Peak RSS | `psutil.Process().memory_info().rss` 50 ms 间隔后台采样 |
| 用户操作步数 | "从打开工具到出图"的鼠标点击 / 命令数（手数） |
| 输出文件大小 | bytes（PNG / SVG / PDF / HTML） |

### 4.2 体验指标（半量化）

| 指标 | 描述 |
|---|---|
| Setup overhead | 首次跑通的安装 + 配置时间（分钟）|
| Input format prep | 工具能否直接吃用户拿到的标准格式（KO 表 / 丰度表）|
| 调参成本 | 切换主题色 / 重画一次出图需多少点击 / 代码修改 |
| 复现脚本 | 工具是否生成可复现 .py / .R / .sh / .yaml |
| 离线可用性 | 是否完全离线（无网络）/ 部分需联网 / web-only |

### 4.3 输出维度

| 指标 | 描述 |
|---|---|
| 视觉变体 | 同一分析能输出多少种图（如 PCoA 可否 PC1 vs PC2 / vs PC3 切换） |
| 出版级格式 | PNG / PDF / SVG / TIFF / EPS 支持多少种 |
| 交互性 | 静态 / hover / 点击 / 拖拽 |
| SI 嵌入度 | 输出是否能直接作 SI 文件（vs 需要后期 PowerPoint 拼接）|

---

## 5. 测量协议（per 工具）

```
Step 1: 准备输入数据（标准化到工具期望的格式）
Step 2: 启动工具（CLI / GUI / web）
Step 3: 跑每个共享分析 3 次（重复测）
Step 4: 测 wall time（time.perf_counter）+ RSS peak（psutil）
Step 5: 截图最终输出 + 记录文件大小 + 记录用户操作步数
Step 6: 填入对应行的 Table T2 / T3
```

测量脚本 stub：[`bench_vs_tools.py`](bench_vs_tools.py)（仅 EnvMeta 侧已就绪 + 占位）。

---

## 6. 预期结果（基于已有 §5.2.8 EnvMeta benchmark + 工具公开 benchmark / 用户经验）

### Table T2 — Type A 共享能力性能对比（预期填表，待 Task 3b 实测）

| 任务 | EnvMeta | Anvi'o | Krona | MicrobiomeAnalyst |
|---|---|---|---|---|
| PCoA | TBD | (不适用) | (不适用) | TBD |
| KO heatmap | TBD | TBD | (不适用) | partial |
| 层级分类 HTML | TBD | (不适用) | TBD | TBD |
| Setup 时间 (min) | <5 (pip install) | 30-60 | 5-10 | 0 (web) |
| 跑通一次完整 pipeline | TBD | TBD | TBD | TBD |
| 输出格式数 | 4（PNG/PDF/SVG/TIFF）| 1-2 | 1 (HTML) | 2 (PNG/PDF) |
| 离线可用 | ✅ | partial | ✅ | ❌ |
| 复现脚本 | ✅ .py | partial | ❌ | ❌ |

### Table T3 — Type B 独家能力矩阵（基本已确定）

见 §1 Type B 表。

---

## 7. Task 3b 执行状态

### 7.1 Windows env 实测结果（2026-05-14）

| 工具 | 状态 | 实际原因 |
|---|---|---|
| **Anvi'o** | 🔴 **不可行** | 需要 raw contigs FASTA + read mapping BAMs；Liu/Grett/Ayala 公开数据集只有 abundance + KO 表，无 contigs；如果重新组装 ~5-10 GB / 数据集 × 4 数据集 + 装 Anvi'o 5-10 GB env = 2-3 天专门工时，不在本周可行范围内 |
| **Krona** | 🔴 **Windows install 失败** | `conda install -c bioconda krona` 在 Windows 失败（`[Errno 22] Invalid argument: 'ktClassifyBLAST'` — Krona Perl bin scripts 在 Windows NTFS 不兼容）；需 WSL2 或 Linux/Mac 重跑（~30 min WSL2 setup + 1-2h 实测）|
| **MicrobiomeAnalyst** | 🟡 **Web 待手动** | 无 CLI/API；需用户手动 Web 上传 + 截图测量；与 GPB 投稿 timeline 异步可行 |

### 7.2 Pivoted 决策

**bioRxiv 优先级 > Task 3b 全完成**。bioRxiv 接受 caveated TBD with methodology
documentation。Table T2 / §5.2.9.2 内显式标注"待 WSL2 / Web 后续实测"，文献中常见做法。

**WSL2 + Linux Krona 回归后做的事**（acceptable post-bioRxiv timeline）：

| 步骤 | 工时 |
|---|---|
| Setup WSL2 Ubuntu + miniconda Linux 端 | 30 min |
| conda install -c bioconda krona（Linux WSL2 中应成功）| 5 min |
| 准备 Liu 2023 mag_taxonomy_labels.tsv → Krona TaxID + abundance input | 30 min |
| 跑 ktImportText 3 次重复 + 测 wall time + 输出 HTML 大小 | 30 min |
| MicrobiomeAnalyst Web 上传 Liu 子集 + 跑 3 个共享分析（PCoA / heatmap / KEGG）| 1-2 h |
| 填 Table T2 / 生成 Figure F10 / 写 Results §2.9.2 narrative | 2-3 h |

**Anvi'o 仍跳过**（无 raw FASTA 实质阻塞，与 timeline 无关）。

### 7.3 GPB 审稿应对

GPB Reviewers 可能问"Why no Anvi'o benchmark?" — 回应：(a) Anvi'o 与 EnvMeta
**workflow paradigm 不同**（探索式 MAG 浏览 vs 出版级图批量生成），analogous
analyses 不直接可比；(b) 公开数据集格式（abundance + KO 表）与 Anvi'o 期望
（contigs.db + profile.db）不兼容，**重新组装 4 个数据集**会使 Anvi'o benchmark
失去 head-to-head 意义；(c) Feature matrix Table T1 已记录 Anvi'o 优势区
（MAG exploration）vs EnvMeta 独家区（cycle inference + hypothesis scoring）。

---

## 8. 参考文献占位

- Anvi'o: Eren et al., 2021 *Nature Microbiol* 6:3-6
- Krona: Ondov et al., 2011 *BMC Bioinformatics* 12:385
- MicrobiomeAnalyst: Chong et al., 2020 *Nucleic Acids Res* 48:W374
- QIIME2: Bolyen et al., 2019 *Nat Biotechnol* 37:852
- ImageGP 2: Chen et al., 2024 *iMeta* 3:e239
