# Paper 3 任务 2 — 第二外部数据集 Benchmark 方案 [HISTORICAL ARCHIVE]

> ⚠️ **本文件已被 v0.9.0/0.9.1 超额完成**（2026-05-09）。Paper 3 实际跑了 **4
> 个外部数据集**（Wei 2024 + Liu 2023 + Grettenberger 2021 + Ayala 2020）+ **3 个
> stress test discrimination test**（含 dominance_score v2 B→A 升级）。当前进度
> 见 [scoring_validation_experiment_results.md](scoring_validation_experiment_results.md)
> 主结果文档 + [stress_test_results.md](stress_test_results.md)。
>
> 本文件保留作为 Stage 1 (Wei 单数据集) 早期计划的 **historical archive**，记录
> 决策过程；**不要据此规划新工作**，请优先看上述主结果文档。
>
> ---
>
> **原创建日期**：2026-05-08
> **原状态**：Phase 0 进行中（已被超越）
> **关联**：[paper3_pre_submission_checklist.md](paper3_pre_submission_checklist.md) 任务 2
> **目标**：提供"EnvMeta 不只能跑作者自己数据"的泛化性证据，回答 iMeta 审稿
> 人 90% 概率会问的"only works on your own data?"

---

## 数据集选定

### 🔵 砷方向（同方向）：Wei et al. 2024 *Microbiome*

| 项 | 内容 |
|---|---|
| 论文 | "Various microbial taxa couple arsenic transformation to nitrogen and carbon cycling in paddy soils" *Microbiome* 12:236 (2024) |
| DOI | [10.1186/s40168-024-01952-4](https://doi.org/10.1186/s40168-024-01952-4) |
| 数据 | NCBI BioProject `PRJNA1068274`（36 paddy samples × 20 高质量 MAGs） |
| 下载入口 | EBI ENA：<https://www.ebi.ac.uk/ena/browser/view/PRJNA1068274> |
| 主题 | 中国广东稻田，pH 4.6-8.0 梯度，As-N-C 三元素耦合 |

**选这个的理由**：
- As + N 是 EnvMeta 4-元素 KB 的命中区
- 36 样本 + pH 梯度完美喂 RDA + 假说评分
- 论文有"As(III) 氧化 + 反硝化共发生"的明确假说，可写为 YAML 喂 S3.5 评分器
- 2024 *Microbiome* 高引（IF 12+）

### 🟢 异方向（跨主题）：Tara Oceans MAGs (Meren Lab curated) ⭐

> **2026-05-08 调研后修正**：原推荐 MUCC v2 (Zenodo 14532347) 实测仅含 Methanoregula 单属子集；Wilson 2024 真实数据散落 4 处且缺丰度矩阵。改选 Tara Oceans，原因：唯一已发布完整四件套的非砷大数据集 + Meren Lab 文档化最完整 + 中国可稳定下载 + 海洋 N/S 命中 EnvMeta 4-元素 KB（N+S）。

| 项 | 内容 |
|---|---|
| 论文 | Delmont et al. 2018 *Nat Microbiol* "Nitrogen-fixing populations of Planctomycetes and Proteobacteria are abundant in surface ocean metagenomes" + 后续多篇 |
| 入口 | <https://merenlab.org/data/tara-oceans-mags/> |
| 规模 | 1888 MAGs × 64 站点（surface / DCM / mesopelagic 三深度）|
| 数据 | MAG fasta + **Anvi'o profile DB（含 mapping coverage 即丰度）** + KEGG/COG 注释 + 样本 metadata 全套 |
| 大小 | ~15 GB（曾被国内多组镜像）|
| 主题 | 海洋 N 固定 + S 循环（与陆地砷修复完全异主题）|

**选这个的理由**：
- ✅ **唯一已发布完整四件套**（MAG / KO / 丰度 / metadata）的非砷数据集
- ✅ Anvi'o profile DB 直接给 MAG mapping coverage（不必从 raw reads 重 mapping）
- ✅ 海洋 N/S 循环 = EnvMeta 4-元素 KB 中 2 个元素的命中
- ✅ Meren Lab 文档化是行业标杆，复现叙事强
- ✅ 64 站点 + 3 深度 = 多组对照天然存在

### 备选（若 Tara 失败）：SMAG 土壤 MAGs（中科院南土所）

Ma et al. 2023 *Nat Commun* "A genomic catalogue of soil microbiomes"：
- DOI 10.1038/s41467-023-43000-z
- 40,039 MAGs + KEGG + GTDB-Tk + abundance（已发表处理后表）
- 中国课题组发布，国内访问稳定
- 切片版（grassland 子集）压到 ~20 GB

---

## 验证方案 5 阶段（5-10 天）

### Phase 0 — 数据下载与可用性验证（1 天）

**Day 1 上午**：
- 下载 MUCC Zenodo zip（~10-20 GB）→ 解压 → 查清单：必须确认含
  - MAG 丰度表（MAG × sample）
  - KO 注释（KO × MAG 或 KO × sample）
  - 样本 metadata（生境 + 深度）
  - GTDB-Tk taxonomy
- 下载 Wei 2024 Supplementary Tables（先尝试免重跑流程的捷径）

**Day 1 下午**：
- 评估 Wei Suppl 是否含 KO 矩阵：
  - **如有** → 进 Phase 2，省 1-2 周
  - **如无** → 决定是否子采样 8 样本跑 KofamScan（3-5 天），或暂搁置

**checkpoint**：写 `paper/benchmarks/external/dataset_inspection_log.md` 记录两份数据的实际可用程度。

### Phase 1 — MUCC 永冻土全流程（2-3 天，优先做）

**Day 2 — reshape**：
将 MUCC 数据 reshape 成 EnvMeta 11 种 FileType 格式：
- `metadata.tsv`（生境 + 深度）
- `mag_abundance.tsv`（MAG × sample）
- `ko_long.tsv` 或 `ko_wide.tsv`（KO 矩阵）
- `mag_taxonomy.tsv`（GTDB-Tk）
- `checkm.tsv`（completeness/contamination）
- 可选 `env.tsv`（pH / depth / 温度）

上传 EnvMeta → 文件管理器自动识别（**这本身就是验证点 1**）

**Day 3 — 跑核心 8 张图**（不必全 14）：
1. 物种堆叠图 — 复现 palsa→bog→fen 物种演替
2. PCoA + PERMANOVA — 复现生境差异
3. MAG 丰度热图 — methanogen 丰度梯度
4. 通路完整度 — N + CH₄ 循环通路三生境差异
5. 基因谱 — 关键 KO（pmoA/mcrA/nifH）的 MAG 分布
6. 🔑 循环图（N 元素）— 自动推 N 循环耦合（EnvMeta 独家）
7. YAML 假说评分 — 写"palsa→bog 过渡含 methanogen 增多"假说，喂评分器
8. 共现网络 — keystone MAG 识别

**Day 4 — 输出归档**：
```
paper/benchmarks/external/wilson_2024_permafrost/
├── README.md              # 数据集 + 预处理 + EnvMeta 参数 + 8 图截图缩略
├── input_data/            # reshape 后的 EnvMeta 输入（≤ 50 MB）
├── envmeta_outputs/       # 8 图 PDF + .py 复现脚本 + stats TSV
├── compare_to_original.md # EnvMeta 出图 vs 原 paper 图 1-1 对照
└── benchmark_table.tsv    # runtime / memory
```

### Phase 2 — Wei 砷场地（路径依 Phase 0 决定）

**路径 A（捷径，1-2 天）— Suppl 含 KO 矩阵**：
- 同 Phase 1 的 reshape + 跑图流程
- 输出 `paper/benchmarks/external/wei_2024_paddy/`

**路径 B（中等，3-5 天）— 子采样 8 样本重跑 KofamScan**：
- 选 4 个 acid + 4 个 neutral pH 样本（覆盖梯度）
- ENA aspera 下载 ~25 GB raw FASTQ
- MEGAHIT + Prokka + KofamScan
- 之后同路径 A

**路径 C（搁置）— 只做 MUCC**：
- 异方向证据已强（4745 MAG / *Nat Microbiol*）
- Methods 章节叙事改为"用 MUCC 永冻土证明跨主题泛化"
- 砷方向泛化性留给"今后工作"或预印本补充

### Phase 3 — 与原论文图一一对照（1 天）

每数据集都答这 5 个验证问题：

| 验证点 | MUCC | Wei |
|---|---|---|
| 1. EnvMeta 文件管理器能识别新数据吗？ | 11 FileType / N 命中 | 同 |
| 2. 能复现原作者展示的 1-2 个核心结果？ | palsa→bog→fen 演替 + N 循环 | As-N 耦合 + pH 梯度 |
| 3. EnvMeta 输出比原论文多了什么？ | 循环图 + 假说评分 | 同 |
| 4. 跑全套耗时？ | 目标 < 30 min | < 30 min |
| 5. 出图 publication-grade？ | 是/否 + 截图 | 同 |

### Phase 4 — Methods 段落初稿（半天）

写 200-300 字 `paper/manuscript/methods_external_validation.md` 给 Methods 第二数据集小节用。

---

## 输出标准（按 paper3 checklist 任务 2 规定）

```
paper/benchmarks/external/<dataset_name>/
├── README.md                  # 数据集来源 + 预处理 + EnvMeta 跑全套 + 结果总结
├── input_data/                # 子采样后输入数据（≤ 100 MB）
├── envmeta_outputs/           # 14 图（或 8 核心图）PDF + .py 复现脚本
├── benchmark_table.tsv        # runtime / memory / 出图质量对比
└── compare_to_original.md     # 对比原 publication 图（如有）
```

---

## 节奏决定点

**节奏 1（保守，5-7 天）**：先 Phase 0 + Phase 1（MUCC），再视 Wei Suppl 决定。
保证任务 2 在 1 周内有完整产出，砷方向留给 review 补。

**节奏 2（激进，10-15 天）**：Phase 0 + Phase 1 + Phase 2 全做。
砷方向 + 跨主题双管齐下，最差 2 周。

---

## 备选数据集（若 1A 或 2A 出问题）

### 砷方向备选 — Yang et al. 2024 *ISME J* 中国 7 矿区尾矿
- DOI 10.1093/ismejo/wrae110；BioProject PRJNA989741；17 metagenomes × 7 矿区
- 全部需重跑流程

### 异方向备选 — Tara Oceans MAG 子集
- Meren Lab 已 curated：<https://merenlab.org/data/tara-oceans-mags/>
- 957 MAGs，~50 GB 子集
- 海洋 N/S 循环（与陆地循环 KB 略有不同）

---

## 维护记录

| 日期 | 事项 |
|---|---|
| 2026-05-08 | 方案存档创建，Phase 0 启动 |
