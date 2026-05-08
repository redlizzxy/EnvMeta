# EnvMeta 论文大纲（iMeta 投稿版）

> **目标期刊**：iMeta（Wiley，IF 24，2022 创刊）
> **文章类型**：Tools and Resources / Methods（方法学论文）
> **预计字数**：5000-7000 words（不含 Methods 详细附录）
> **预计 figures**：8-9 张（iMeta 偏好图丰富）
> **预计 tables**：2 张
>
> **创建日期**：2026-05-07
> **作者**：redlizzxy（中国地质大学北京）
> **状态**：大纲（待导师审阅）

---

## 0. iMeta 投稿要点速查

| 项 | 要求 |
|---|---|
| 出版社 | Wiley |
| 版式 | OA only，APC ~¥10000-15000 |
| 文章类型 | Article / Methods / Tools and Resources / Reviews |
| 字数 | Article: ~5000-7000；Methods 章节可放附录 |
| 结构 | Abstract + Highlights + Graphical Abstract + IMRAD |
| Figure 风格 | 彩色多图 + 大量 schema 图 + 案例展示 |
| 推荐审稿人 | 投稿系统可填 ≥ 3 位（**这里挂刘永鑫**）|
| 引用格式 | Wiley iMeta 模板（参考 https://onlinelibrary.wiley.com/journal/2770596x） |
| 参考文献 | 优先 iMeta + 中科院 BIG 系列 + 微生物组顶刊 |
| 数据可用性 | 强制要求 GitHub + Zenodo DOI + 在线 demo |

---

## 1. Title 候选（3 个）

| # | Title | 卖点 |
|---|---|---|
| **A** | **EnvMeta: An interactive visualization platform for environmental metagenomics with built-in biogeochemical cycle inference and hypothesis scoring** | 完整 + 工具名 + 元素循环 + 假说评分 |
| B | EnvMeta: From metagenomic data to publication-ready figures with biogeochemical cycle automation | 强调"出版级"（吸引研究生用户）|
| C | Bridging environmental microbiology and visualization: A GUI-driven platform with element cycle inference (EnvMeta) | 强调"桥接"（吸引非生信用户）|

→ 推荐 **A**，最完整且关键词命中。副标题可考虑 "for non-bioinformaticians"。

---

## 2. Highlights（iMeta 强制 4 条）

```
- A GUI-driven Streamlit platform integrating 12 publication-ready
  visualizations for environmental metagenomics, with synchronized
  Python/R script generation for full reproducibility.

- Automated biogeochemical cycle inference for As/N/S/Fe at
  KEGG-driven pathway resolution (4 elements × 18 pathways × 57 KOs)
  with permutation-based confidence labeling.

- A YAML-based hypothesis scoring engine that converts user-supplied
  mechanistic claims into evidence-weighted scores with null
  distribution testing.

- Standalone interactive HTML export (~400 KB, offline-ready) embeds
  the entire analysis as the supplementary information itself,
  enabling reviewer-level reproducibility.
```

---

## 3. Graphical Abstract 草图

```
┌─────────────────────────────────────────────────────────────────────┐
│                                                                     │
│  📁 Upload metagenomic data    ───►   🔍 Auto-recognize 11 types    │
│  (TSV/CSV/FASTA/...)                  (abundance/KO/taxonomy/...)    │
│                                                                     │
│                              ▼                                       │
│                                                                     │
│   ┌──────────────────┐   ┌──────────────────┐   ┌──────────────────┐│
│   │ 📊 12 visualiz.  │   │ ⚛️ Cycle infer.  │   │ 🧪 Hypothesis    ││
│   │ • α/β-PCoA/RDA   │   │ • 4 elements     │   │ • YAML schema    ││
│   │ • Heatmaps       │   │ • Mockup-10      │   │ • Null perm.     ││
│   │ • LEfSe/log2FC   │   │ • Cross-element  │   │ • Sensitivity    ││
│   │ • MAG quality    │   │   coupling       │   │ • 9-tier label   ││
│   └──────────────────┘   └──────────────────┘   └──────────────────┘│
│                                                                     │
│                              ▼                                       │
│                                                                     │
│   📦 Export: PNG/PDF/SVG/TIFF + .py reproducible script + Bundle    │
│   🌐 Standalone interactive HTML (~400 KB, offline)                 │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

**布局**：上 1/3 数据上传 → 中 1/3 三栏功能（左 reads-based / 中 cycle / 右
hypothesis）→ 下 1/3 导出。

---

## 4. Abstract（200-250 words 草稿）

> Environmental metagenomics has matured technically but the downstream
> analysis-to-publication pipeline remains a bottleneck for graduate
> students and non-bioinformatics researchers. Existing tools either
> require command-line proficiency (Anvi'o, MicrobiomeAnalyst), focus
> on a single visualization type (Krona), or are constrained by
> commercial cloud platforms with limited customization. We present
> **EnvMeta**, an open-source Streamlit-based platform that unifies
> 12 publication-ready visualizations with three differentiating
> capabilities: (1) automated biogeochemical cycle inference for
> As/N/S/Fe at KEGG-driven pathway resolution (4 elements × 18
> pathways × 57 KOs), with permutation-based confidence labels and
> sensitivity scanning; (2) a YAML-based hypothesis scoring engine
> that evaluates user-supplied mechanistic claims against the data
> with null distribution testing and weight robustness checks;
> (3) standalone interactive HTML export (~400 KB, offline-ready) that
> embeds the analysis itself as supplementary information, enabling
> reviewer-level reproducibility. EnvMeta is designed under a
> "domain-neutral, user-supplied knowledge, fully offline, fork-rather-
> than-community" philosophy. We demonstrate EnvMeta on an arsenic
> slag-steel slag bioremediation case study (168 MAGs × 10 samples)
> and validate its scalability on a second public dataset
> (Tara Oceans subset / Oak Ridge uranium / TBD). EnvMeta is freely
> available at https://github.com/redlizzxy/EnvMeta with an online
> demo at https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/.

---

## 5. 主体章节大纲

### 5.1 Introduction（~700-900 words）

**段落结构**：

1. **背景与痛点**（200 words）
   - 环境宏基因组学技术成熟（Illumina + Nanopore + MAG 流程标准化）
   - 但下游 "数据 → 出版图" 这条最后一公里成为研究生瓶颈
   - 痛点 1：需要写大量定制 R/Python 脚本（matplotlib/ggplot2）
   - 痛点 2：测序公司云平台预设固定，无法自定义
   - 痛点 3：元素循环图必须手画 PPT（Adobe Illustrator）
   - 痛点 4：假说论证靠 reviewer 主观，缺统计标尺
   - 痛点 5：复现性差，论文 SI 是静态 PDF

2. **现有工具的局限**（250 words）
   - Krona（饼图分类，不能定量元素循环）
   - Anvi'o（强大但学习曲线陡，CLI 优先）
   - MicrobiomeAnalyst（Web 平台，定制有限，需联网）
   - phyloseq / Shiny-phyloseq（R 生态，对零基础用户门槛高）
   - 测序公司云平台（封闭，不支持元素循环图 / 假说评分）
   - **空白**：GUI + 元素循环 + 假说评分 + 离线 + 开源 → 业界无任何工具同时具备

3. **EnvMeta 的设计哲学**（200 words）
   - 5 大原则（领域中立 / 用户自带知识 / 完全离线 / Fork 而非社区 / 描述而非断言）
   - 5 层架构（L1 通用循环图 → L2 假说评分 → L3 插件框架 → L4 Bundle → L5 KB 工具）
   - 目标用户群（环境微生物方向的硕博研究生 + 非生信背景研究人员）

4. **本文贡献**（150 words）
   ```
   1. We present EnvMeta, the first GUI platform integrating 12
      publication-ready visualizations + automated biogeochemical
      cycle inference + YAML hypothesis scoring + standalone
      interactive HTML.
   2. We propose a permutation-based, sensitivity-aware framework
      for inferring element cycle activity from metagenomic data,
      moving beyond descriptive Krona-style charts toward
      hypothesis-evaluable visualizations.
   3. We demonstrate EnvMeta on As/Fe bioremediation case study
      and a second public dataset, with full reproducibility via
      Fork Bundles and DOI-tagged Zenodo releases.
   ```

### 5.2 Results（~2500-3000 words / 7-8 figures）

**节奏**：每个 Result 子节 1 个图（核心 figure）+ 1 段叙事 + 1 段 quantitative 对比。

#### 5.2.1 EnvMeta 整体架构（Figure 1）

- 5 层架构图
- 11 FileType + 14 图 + 3 创新模块（cycle/hypothesis/HTML）
- 数据流：上传 → 识别 → 分析 → 调参 → 导出

**Figure 1 内容**：架构 schema 图（参考 paper/figures/screenshot_*.png 现有素材重画 SVG）

#### 5.2.2 12 张出版级图 + 同步代码生成（Figure 2）

- 7 reads-based + 5 MAG-based
- 重点：**调参 GUI 控件 + 实时预览 + 代码同步生成**
- 量化：vs 原 R/Python 脚本，作图时间从 30-60 min → 1-2 min（ref `time_comparison.md`）

**Figure 2 内容**：4 panel
- A) 堆叠图（reads-based）+ GUI 控件截图
- B) PCoA + 调参 → 实时刷新
- C) MAG 丰度热图（4 层参数面板）
- D) 自动生成的 .py 复现脚本片段

#### 5.2.3 元素循环图自动推断（Figure 3 — **CORE FIGURE**）

- KEGG-driven KB（4 元素 × 18 通路 × 57 KO）
- 3 阶段推断：**S1 去偏 → S2 999 次置换 → S3 敏感度扫描**
- Mockup 10 合并细胞布局（substrate-gene-product 三段式）
- 跨元素化学物耦合（雌黄沉淀 As(III)↔H₂S → As₂S₃）

**Figure 3 内容**：
- A) KB schema（4 元素 × 18 通路）
- B) S1 去偏 → S2 置换 → S3 敏感度算法流程
- C) 推断后的元素循环图（cell_renderer 输出，砷渣 case）
- D) 跨元素耦合放大（雌黄沉淀化学物锚点）

#### 5.2.4 YAML 假说评分器（Figure 4）

- 5 类 claim（taxon_anchored / pathway_pair / env_correlation / cross_element_coupling / population_structure）
- 评分公式：weighted_sum × null_p × veto_logic
- 9 档解读标签（strong / suggestive / spurious? / unknown / ...）
- 权重敏感度（用户改 claim 权重，看 overall 是否稳健）

**Figure 4 内容**：
- A) YAML 输入 → 评分输出 schema
- B) 999 次置换 null 分布 + 观察值
- C) 9 档解读标签 + 决策树
- D) 权重敏感度热图（claim weight × overall score）

#### 5.2.5 独立交互 HTML 导出（Figure 5）

- D3.js v7 inline 嵌入（无外部依赖）
- 400-550 KB，**比典型 PDF SI 还小**
- 包含：循环图全交互 + 假说评分表 + 跨组对比 + SVG 导出
- "**论文 SI 即软件本身**"（reviewer 下载单文件即可亲手审计）

**Figure 5 内容**：
- A) HTML 4 tab UI 截图
- B) 循环图节点拖拽 + 化学物锚点高亮
- C) 假说评分 claim 表 + null_p 分布
- D) HTML 文件大小 vs Krona/iTOL/Anvi'o SI 对比柱状图

#### 5.2.6 Fork Bundle —— 论文与工具绑定发布（Figure 6）

- 一篇论文一个 Bundle.zip（KB + YAML + config + KEGG 快照 + sample data）
- 用户下载 Bundle → 一键加载 → 完全复现论文图
- **解决"论文发表后工具版本漂移导致复现失败"问题**

**Figure 6 内容**：
- Bundle 结构 + 加载流程
- 复现案例：从 Bundle 加载 → 5 分钟内重生成论文图

#### 5.2.7 案例研究：砷渣-钢渣微生物修复（Figure 7）

- 168 MAG × 10 sample × 57 目标 KO
- CK / A（低钢渣）/ B（高钢渣）三组
- 关键发现：铁氧化固砷 + 硫还原沉砷的多元素协同
- **不是新机制论文（那是 Paper 1/2 的事），而是展示 EnvMeta 如何辅助机制发现**

**Figure 7 内容**：
- A) 元素循环图（CK vs A vs B）
- B) 假说评分跨组对比
- C) 关键 MAG 承载者切换（Gallionella → Thiobacillus）

#### 5.2.8 第二外部数据集 benchmark（Figure 8 / Table 2）

- 数据集候选：Tara Oceans 子集（500+ MAG × 100+ sample，海洋）/ Oak Ridge 铀污染（陆地）
- **跑 14 图全套 + 元素循环 + 假说评分**
- 量化：runtime / memory / 出图质量（vs 原始 publication）
- 证明：EnvMeta 不是只能跑作者自家砷渣数据，泛化能力可演示

**Figure 8 / Table 2 内容**：
- Figure 8：第二数据集元素循环图（演示 N/S 循环为主，区别于砷渣 case）
- Table 2：runtime / memory / scaling 对比表（小 / 中 / 大三档数据集）

#### 5.2.9 vs 竞品对比（Figure 9 / Table 1）

- vs Krona / Anvi'o / MicrobiomeAnalyst / 测序公司云平台
- 对比维度：12 图覆盖度 / 离线可用 / 元素循环 / 假说评分 / GUI 友好度 / 开源 / 价格

**Table 1 内容**：参考现有 `paper/tool_comparison.md` 整理为正式表格。

---

### 5.3 Discussion（~700-900 words）

1. **EnvMeta 的差异化卖点**（200 words）
   - 工具论文常见误区："功能多 = 价值大"
   - EnvMeta 的核心不是功能多，而是 **3 个业界空白**（cycle / hypothesis / HTML）
   - 元素循环 + 假说评分 = "**从描述性可视化迈向假说评估可视化**"

2. **vs 测序公司云平台**（150 words）
   - 真实竞争对手是云平台，不是 Anvi'o
   - 云平台的优势：开箱即用 / 大数据集
   - 云平台的劣势：封闭 / 不支持元素循环 / 不可定制 / 数据出境
   - EnvMeta 定位 = "**比云平台更灵活 + 比 Anvi'o 更易用**" 的中间层

3. **Fork-rather-than-community 模型的合理性**（150 words）
   - 不维护中心化机制目录，每篇论文绑定一份 Bundle
   - 类似 R bioconductor vs Linux distros 的差别
   - 优势：社区维护负担 0 / 论文复现性 100%
   - 劣势：用户跨论文需要 fork / 学习成本略高

4. **局限性**（200 words，**诚实**）
   - User study 样本有限（n=2 内测 + v2 在跑）；以 case study 风格诚实报告
   - KB 仅限 As/N/S/Fe（不含 C/P/H）—— 投稿后可由用户社区扩展
   - Streamlit GUI 在大数据集（10000+ MAG）下渲染有性能上限（推荐 ≤ 5000 MAG）
   - 假说评分依赖用户提供机制 YAML，错误 YAML 可能产生误导分数（schema 校验是软约束）

5. **Future work**（150 words）
   - L3 插件框架（用户上传 Python `analyze()` → 自动注册 GUI）
   - KB 扩展到 C/P/H 循环
   - D3.js 编辑器（用户直接拖拽编辑节点 / 通路）
   - LLM-assisted hypothesis YAML 起草

---

### 5.4 Methods（~1500-2000 words / 可放附录）

1. **文件识别模块**（300 words）
   - 11 FileType 表头规则匹配 + 反向索引
   - abundance 分 MAG/TAXON 两档（confidence 0.95 / 0.88）

2. **14 图分析引擎**（300 words）
   - reads-based / MAG-based / cycle 三类引擎
   - 共享 4 层参数面板（`envmeta/analysis/_mag_common.py`）
   - matplotlib + seaborn 出版级渲染

3. **元素循环推断算法**（500 words）
   - **S1 去偏**：CLR 变换 + 抽样标准化 + 3 档阈值 Top-1 contributor 一致性
   - **S2 置换检验**：999 次 shuffle + perm_p 计算 + 5 档 confidence 标签
   - **S3 评分**：claim 类型路由 + weighted_sum + null_p + veto logic
   - 数据结构：`CycleData / EnvCorrelation / SensitivityRow / HypothesisScore`（见 `envmeta/geocycle/model.py`）

4. **YAML 假说 schema**（200 words）
   - 5 类 claim：taxon_anchored / pathway_pair / env_correlation / cross_element_coupling / population_structure
   - veto_reasons 9 类
   - schema 校验通过 jsonschema

5. **D3.js 交互 HTML 导出**（200 words）
   - inline embed D3 v7（~280 KB）+ 400 KB total HTML
   - 4 tab UI（架构 / 循环图 / 假说 / 环境相关）
   - 客户端 SVG 序列化 + download

6. **Benchmark 实施**（300 words）
   - 数据集：砷渣 168 MAG × 10 sample（自有）+ 第二数据集 ≥ 500 MAG（待选）
   - 测试硬件：4 核 CPU / 16 GB RAM / Python 3.11
   - 对照工具：Krona 2.x / Anvi'o 8.x / MicrobiomeAnalyst Web 版本快照
   - 量化指标：runtime / memory / 用户操作步数

7. **实施细节**（150 words）
   - Streamlit 1.30+
   - matplotlib 3.8+
   - Python 3.11+
   - 测试套件：pytest 293/293 全绿
   - CI：GitHub Actions（待补）

---

### 5.5 Data and Code Availability

```
- Source code: https://github.com/redlizzxy/EnvMeta (MIT License)
- Online demo: https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/
- Sample data: tests/sample_data/ (arsenic slag bioremediation, subset)
- Zenodo DOI: TBD (release v1.0 at submission)
- Documentation: docs/data_preparation_zh.md + docs/install_for_beginners.md
- Issue tracker: GitHub Issues
```

---

### 5.6 Acknowledgements

```
We thank Dr. Yong-Xin Liu (CAS Microbiology Institute, founder of
iMeta) for foundational training in metagenomic upstream/downstream
analyses, whose published pipelines and instructional materials
inspired several of EnvMeta's design choices for accessibility.

We thank [supervisor name] for guidance and [research group] for
data sharing and beta testing. We acknowledge feedback from
n=8+ graduate student users in the v2 user study.

Funding: [supervisor's grant ID] / China Geological University
(Beijing) graduate research support.
```

---

### 5.7 References（约 50-60 篇）

**重点引用领域**（按 iMeta 编辑偏好排序）：

1. **iMeta 系列论文**（5-8 篇，编辑友好度+）
   - 刘永鑫 et al. EasyAmplicon / EasyMetagenome 系列
   - iMeta 已发表的工具论文（参考其格式）

2. **微生物组顶刊**（10-15 篇）
   - Nature Microbiology / ISME / Microbiome
   - QIIME 2 / DADA2 / phyloseq

3. **可视化工具**（8-10 篇）
   - Krona / Anvi'o / MicrobiomeAnalyst / iTOL / Cytoscape

4. **元素循环 + 生物地球化学**（10-15 篇）
   - KEGG / MetaCyc 数据库论文
   - 砷代谢、铁氧化、硫还原经典论文
   - 多元素耦合机制（如 Stocker / Konhauser 的铁循环综述）

5. **统计方法**（5-8 篇）
   - PERMANOVA / Bray-Curtis / Aitchison
   - 置换检验 + bootstrap

6. **Reproducibility / FAIR**（3-5 篇）
   - Wilson et al. Good enough practices
   - FAIR principles

---

## 6. Figures 总表

| # | Figure | 类型 | 现有素材 | 待补 |
|---|---|---|---|---|
| F1 | EnvMeta 整体架构 | Schema | 部分（CLAUDE.md 内嵌）| 重画为 SVG |
| F2 | 12 图 + 调参 + 代码生成 | Screenshots + 截图组合 | `paper/figures/screenshot_*.png` | 排版 4-panel |
| F3 | 元素循环图自动推断 ⭐ | Schema + Output | `paper/figures/cycle_*.png` | 加 S1-S3 流程 schema |
| F4 | YAML 假说评分器 | Schema + Output | `paper/figures/hypothesis_*.png` | 加 null_p 分布图 |
| F5 | 独立交互 HTML 导出 | Screenshots | `paper/bundles/*.html` 截图 | 新拍 4 panel |
| F6 | Fork Bundle 结构 | Schema | 无 | 新画 |
| F7 | 砷渣案例研究 | Result | `paper/figures/cycle_arsenic.png` | 多组对比 panel |
| F8 | 第二数据集 benchmark | Result | 无 ❌ | **必补**（待选数据集后跑） |
| F9 | vs 竞品对比柱图 | Comparison | `paper/tool_comparison.md` 整理 | 改成图 |

**T1 = vs 竞品矩阵表**；**T2 = 性能 benchmark 表**。

---

## 7. 投稿前清单（按完成度）

### 已完成 ✅

- [x] 核心代码 v0.8.1（293/293 测试全绿）
- [x] 12 图 + 循环图 + 假说评分 + Bundle + HTML 全功能
- [x] tool_comparison.md / time_comparison.md
- [x] paper/figures/ 部分截图
- [x] 在线 demo + GitHub repo
- [x] 期刊选定 + 投稿决策

### 投稿前必补 ⏳

- [ ] **R 侧侧对照 11 图**（仅 log2fc 1 张）—— 1-2 周
- [ ] **第二外部数据集**复现 + benchmark —— 2-3 天
  - 推荐数据集：先用 Tara Oceans 子集（海洋 N/S 循环演示）+ Oak Ridge 铀污染（陆地金属循环）
  - 用户已确认服务器有剩 → 公开数据下载 + EnvMeta 跑全套
- [ ] **English README + LICENSE + Zenodo DOI** —— 4-6h
- [ ] **F1 / F6 schema 图重画** —— 1 天
- [ ] **F3 / F4 schema 图补 S1-S3 流程 + null_p 分布** —— 半天
- [ ] **F5 HTML 截图 4 panel** —— 半天
- [ ] **F8 第二数据集结果图** —— 跟数据集复现一起做
- [ ] **F9 / T1 / T2 整理为图表** —— 1 天
- [ ] **References 整理** —— 1-2 天
- [ ] **User study v2 数据回收 ≥ 8 份**（目前已收 ~2-3 份，等回收）

### 投稿当天

- [ ] 邮件刘永鑫确认推荐审稿人
- [ ] 投稿系统填表（Authors / Affiliations / Funding / Recommend reviewers）
- [ ] Cover Letter（强调 3 大空白 + iMeta 选题契合）
- [ ] Conflict of Interest 声明

---

## 8. 时间线（再校准）

```
2026-05-07   今日：大纲完成
2026-05-08 ~ 05-22   补素材（R 对照 + 第二数据集 + 英文文档 + Zenodo DOI）
                     ── 工时 ~2-3 周
2026-05-23 ~ 06-15   论文初稿（Abstract + Intro + Results + Methods + Discussion）
                     ── 工时 ~3 周
2026-06-15 ~ 06-30   导师 + 课题组内审 + 修改
                     ── 工时 2 周
2026-07-01 ~ 07-10   英文润色（DeepL + Grammarly + 投稿前精修）
2026-07-10 ~ 07-14   邮件刘永鑫 + Cover Letter
2026-07-15           **投 iMeta** ⏱ 倒计时锚点
2026-09-15           一审决定（estimated minor / major revision）
2026-09-15 ~ 10-15   修订
2026-10-15           修订投回
2026-11-15 ~ 12-15   **接收 ✅**（赶秋招简历 + 12 月节点）
2027-01 ~ 02         online ✅（赶 2027-05 激励截止前）
```

---

## 9. 维护记录

| 日期 | 事项 | 备注 |
|---|---|---|
| 2026-05-07 | 大纲初版 | 待导师审阅 |
