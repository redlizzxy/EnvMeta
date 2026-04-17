# EnvMeta — 环境微生物宏基因组可视化分析平台

## 项目概述

EnvMeta 是一个面向环境微生物研究的宏基因组下游可视化分析桌面工具。核心解决"从数据文件到发表级图形"的最后一公里问题，特色功能是自动生成微生物-环境互作的生物地球化学循环图。

**目标用户**：环境微生物方向的硕博研究生和科研人员。

**核心差异化**：现有工具（Shiny-phyloseq、plotmicrobiome、Anvi'o 等）无一同时满足：GUI 交互调参 + 代码同步生成 + 生物地球化学循环图自动生成。EnvMeta 填补这个空白。

## 技术栈

- **前端 UI**：Streamlit
- **交互图表**：Plotly（预览用）
- **高质量输出**：matplotlib（最终渲染出版级 PDF/TIFF）
- **循环图编辑器**：D3.js（通过 st.components 嵌入）
- **后端逻辑**：Python（pandas / scipy / networkx / scikit-bio）
- **打包分发**：pip install + streamlit run（暂不做 exe）

## 项目结构（目标）

```
envmeta/
├── app.py                    # Streamlit 入口
├── CLAUDE.md                 # 本文件（Claude Code 自动读取）
├── envmeta/                  # 核心包
│   ├── __init__.py
│   ├── file_manager/         # 模块A：文件识别
│   │   ├── detector.py       #   文件类型检测（表头规则匹配）
│   │   ├── validator.py      #   schema 校验
│   │   └── schemas/          #   各文件类型的 JSON schema
│   ├── analysis/             # 模块B：分析引擎
│   │   ├── base.py           #   统一接口 analyze(files, params) → figure + stats
│   │   ├── stackplot.py      #   物种组成堆叠图
│   │   ├── pcoa.py           #   PCoA + PERMANOVA
│   │   ├── alpha.py          #   α多样性箱线图
│   │   ├── rda.py            #   RDA 排序图
│   │   ├── lefse.py          #   LEfSe 差异分析
│   │   ├── heatmap.py        #   元素循环基因热图
│   │   ├── log2fc.py         #   基因 log2FC 对比
│   │   ├── mag_quality.py    #   MAG 质量散点图
│   │   ├── mag_heatmap.py    #   MAG 丰度热图
│   │   ├── pathway.py        #   通路完整度气泡图
│   │   ├── gene_profile.py   #   MAG 元素循环基因谱
│   │   └── network.py        #   共现网络图
│   ├── params/               # 模块C：参数管理
│   │   ├── param_registry.py #   参数注册与管理
│   │   └── presets.py        #   预设方案（毕业论文版/SCI版）
│   ├── export/               # 模块D：导出
│   │   ├── figure_export.py  #   图形导出（PDF/SVG/TIFF/PNG）
│   │   └── code_generator.py #   调参后生成可运行的 .py 脚本
│   └── geocycle/             # 模块E：循环图引擎
│       ├── knowledge_base/   #   元素循环知识库（JSON）
│       │   └── elements.json #   As/S/Fe/N/C/P 的 KO→通路→反应映射
│       ├── inference.py      #   通路推断（扫描 KO 表→活跃通路）
│       ├── model.py          #   循环图数据模型（节点/边/层）
│       └── renderer.py       #   可视化渲染
├── tests/                    # 测试
│   ├── test_detector.py
│   └── sample_data/          # 测试用样例文件（从论文真实数据精简）
├── docs/                     # 文档
├── requirements.txt
└── README.md
```

## 现有代码基础

本项目基于一个已完成的硕士论文可视化代码库。现有 `scripts/` 目录下有：

### Python 脚本（scripts/python/）
- `01_physicochemical.py` — 理化指标柱状图
- `05_gene_heatmap_log2fc.py` — 元素循环基因热图 + log2FC（包含完整的 As/S/Fe/N 四大元素 KO→基因→通路映射字典）
- `06_MAG_quality.py` — MAG 质量散点图
- `06_MAG_gene_profile.py` — MAG 元素循环基因谱热图
- `07_MAG_abundance_heatmap.py` — Top30 MAG 丰度热图（544行，含聚类、三段非线性配色）
- `08_pathway_completeness.py` — 通路完整度热图 + 气泡图
- `09_cooccurrence_network.py` — 共现网络图（力导向布局，Keystone 标注）

### R 脚本（scripts/R/）
- `01_tax_stackplot.R` — 门/属/种物种组成堆叠图
- `02_alpha_diversity.R` — α多样性箱线图
- `02_beta_PCoA.R` — PCoA + PERMANOVA
- `03_RDA.R` — RDA 排序图
- `04_LEfSe.R` — LEfSe LDA 分析

### 关键数据结构
所有脚本使用 argparse 命令行接口。需要改造为函数接口 `generate_figure(params_dict) → Figure`。

## 元素循环知识库

已有的 KO 字典（在 `05_gene_heatmap_log2fc.py` 中定义），覆盖 51 个 KO：

- **砷代谢（11个KO）**：arsC-grx/trx/myc、aoxA/B、arsB、ACR3、arsA、arsH、AS3MT、arsR
- **氮代谢（17个KO）**：narG/H/I、napA/B、narB、nirK/S/A、CYP55、norB、nosZ、amoA/B/C、nifD/K
- **硫代谢（15个KO）**：cysJ/I/H、sir、aprA/B、dsrA/B、sqr、soxA/X/B/Y/Z、TST
- **铁代谢（8个KO）**：ABC.FEV.A/P/S、TC.FEV.OM、afuA、FTR、fur、tonB

## 研究背景（用于循环图推断）

本工具源自一个砷渣微生物修复的宏基因组研究：
- 三组实验：CK（对照）、A（低钢渣配比）、B（高钢渣配比）
- 配色：CK=#4DAF4A（绿）、A=#377EB8（蓝）、B=#E41A1C（红）
- 核心机制："铁氧化固砷 → 硫/氮循环调控Eh → 局部硫化沉砷"
- 关键物种：Sulfuricaulis（硫氧化）、Gallionella（铁氧化）、Thiobacillus（硫氧化）等

## 五大功能模块

### 模块A：智能文件管理器
上传文件 → 表头规则匹配 → 自动识别类型（metadata/丰度表/KO注释/CheckM/距离矩阵/α指数/环境因子）→ 推荐可用分析

### 模块B：可视化分析引擎
内置 12 种图表，按 Reads-based（堆叠图/α/β/RDA/LEfSe/热图/log2FC）和 MAG-based（质量图/丰度热图/通路图/基因谱/网络图）分层

### 模块C：交互式图形调参
Streamlit 控件实时调参（颜色/字体/大小/图例等）→ 参数字典 → 刷新预览。先实现"调完点刷新"，后续优化为实时预览。

### 模块D：高质量导出
图形（PDF/SVG/TIFF/PNG）+ 可运行代码（.py）+ 统计结果。预设"毕业论文版"和"SCI投稿版"。

### 模块E：生物地球化学循环图生成器 ⭐
核心创新。双层结构（环境层+微生物层），从 KO 注释自动推断活跃通路，用户可交互编辑（添加/删除节点和通路），导出矢量图。

## 开发规范

- **语言**：Python 3.11+，代码注释用中文
- **Git commit 前缀**：feat（新功能）/ fix（修bug）/ refactor（重构）/ docs（文档）/ test（测试）
- **函数**：写 docstring，参数类型标注
- **分支**：main（稳定版）→ dev（日常开发）→ feature/xxx（功能分支）
- **测试数据**：放在 tests/sample_data/，用论文真实数据的精简版

## 日常开发循环

每次开始开发时：

```powershell
conda activate envmeta
cd D:\workdata\envmeta
streamlit run app.py       # 浏览器 http://localhost:8501
```

提交代码：

```powershell
git add <改动的文件>       # 避免 git add -A（会误入 .env 等敏感文件）
git commit -m "feat: xxx"  # 前缀：feat/fix/refactor/docs/test
git push                    # 推送到 origin/master
```

**远程仓库**：`https://github.com/redlizzxy/EnvMeta`（Private）

## 开发路线图

- **Phase 0（已完成）**：项目骨架 + git init + conda 环境 + 知识库 v1
- **Phase 1（当前）**：文件识别 + 堆叠图/PCoA/热图 + 基础调参 + 导出
- **Phase 2**：全部 12 种图表 + 代码生成器
- **Phase 3**：循环图 v1（推断引擎 + 静态渲染）
- **Phase 4**：循环图交互编辑 + 产品打磨 → v1.0 发布

## 产品定位与核心设计决策（2026-04-17 确认）

### 定位

EnvMeta **不是**预装所有环境机制的权威工具，而是一个：

> **可定制的环境生信分析框架 + 通用循环图推断引擎 + 论文-EnvMeta 绑定发布协议**

### 核心设计原则

| 原则 | 含义 | 理由 |
|---|---|---|
| **领域中立** | 核心代码不内嵌任何具体机制/研究主题目录 | 科学问题空间无穷，不可枚举；避免确认偏差 |
| **用户自带知识** | 机制 YAML / 知识库扩展 / 分析插件都由用户提供 | 维护人不承担社区知识更新负担 |
| **完全离线** | 所有功能本地可用，无网络依赖 | 符合学术界数据保密 + 降低维护复杂度 |
| **Fork 而非社区** | 每篇论文绑定一份定制 EnvMeta 发布 | 分布式发布，无中心化维护瓶颈 |
| **描述而非断言** | 循环图输出"谁承载什么 + 什么与什么相关"，不下因果结论 | 因果判断是用户职责；工具避免附和假说 |

### 拒绝的设计

- ❌ 自动检测研究主题（env 列名 + KO 频度**推不出**"砷修复"等研究意图）
- ❌ 下拉选研究领域（无法穷举科学问题）
- ❌ 预装竞争机制目录（违反"不维护社区"原则）
- ❌ 嵌入 KEGG 快照（过期 + 违反离线原则）
- ❌ 完整 GUI KB 编辑器（5% 用户会用，投入回报低）

### 采用的架构（5 层）

| 层 | 内容 | 状态 |
|---|---|---|
| L1 通用循环图内核 | 描述性推断 + 去偏 + 敏感度 | Phase 3 v1 已有，S1 去偏待做 |
| L2 机制 YAML 评分器 | 用户上传假说 YAML → 证据评分 | S3 计划 |
| L3 插件框架 | 用户上传 Python `analyze()` → 自动注册 GUI | S8 计划（可选）|
| L4 Fork Bundle | 打包配置 + 插件 + YAML + KB + 样例数据 | S4 计划 |
| L5 KB 工具 | schema 校验 + template + diff/merge CLI | S5 计划 |

### 市场现实（诚实校准）

- 预期峰值：**50-200 活跃用户 / 1-2 年**（不是上千）
- 真实采用：自己课题组（确定）+ 同行复现（每篇 7-35 人）+ 小众外部（5-50 实验室）
- 论文目标：**iMeta / Bioinformatics / Front Microbiol 方法学论文**；博士论文章节完美匹配

### KB 制作难度分层

| 场景 | 工作量 | 难度 |
|---|---|---|
| 给现有元素加通路 | 30 min | ⭐ |
| 加新通路 + KO | 1-2 h | ⭐⭐ |
| 加新元素循环（如 Hg / C）| 4-8 h | ⭐⭐⭐⭐ |
| 换研究领域（全碳循环 100+ KO）| 1-2 天 | ⭐⭐⭐⭐⭐ |

**关键理解**：KB 成本 = 一次性 + 与论文研究副产物重叠（占论文周期 <1%）。EnvMeta
降低**机械**成本 ~50%，但**智力**成本（哪些 KO 属于哪通路）无法降低。

### 当前执行路线：B（27h / 8 session）

S1 v1 去偏 → S2 置换零假设 + 可信度标签 → S3 机制 YAML 评分器 → S4 Fork Bundle →
S5 KB 工具 → S6 `mag_heatmap` → S7 `network`(Gephi-prep) → S9 论文起草

S8 插件框架推迟到论文接收后再做。完整计划见 `C:\Users\REDLIZZ\.claude\plans\logical-drifting-metcalfe.md`。

## Backlog（累积 TODO，每次迭代前筛选）

**Phase 1 迭代 2 候选**：
- PCoA + PERMANOVA 分析器（移植 `scripts/R/02_beta_PCoA.R`）
- 元素循环基因热图（移植 `scripts/python/05_gene_heatmap_log2fc.py`）
- 模块 C 参数面板结构化（ParamRegistry + 预设方案）
- **堆叠图排序选项**（按均值 / 按最大值 / 按中位数 / 反转顺序，用户 2026-04-13 提出）
- 模块 A 扩展：识别 alpha_diversity、distance_matrix、checkm_quality 等
- combined 样式（sample + group 并排）
- 装 R 后做 EnvMeta vs R 的侧侧 PDF 对比

**UX 待优化**：
- "文件管理"首页同名覆盖策略（当前只处理新 name）
- Streamlit 里的中文字体支持（matplotlib 渲染图里中文字显示方块时）

## 论文发表数据积累

本工具计划作为方法学论文发表（目标期刊：iMeta / Bioinformatics / Frontiers in Microbiology）。
开发过程中需要持续积累以下四类数据：

### 积累规则

1. **功能验证**：每完成一个分析模块，用 `tests/sample_data/` 中的论文真实数据跑一遍，将输出保存到 `paper/benchmarks/validation/`，与原始 R/Python 脚本的输出对比，确认结果一致。
2. **效率对比**：每完成一个模块，记录"传统方式（写代码）需要多少行代码/多少步骤/多少时间" vs "EnvMeta 需要多少次点击/多少时间"，写入 `paper/benchmarks/time_comparison.md`。
3. **截图留档**：每个模块完成后截一张 EnvMeta 的界面截图，保存到 `paper/figures/screenshot_模块名.png`。
4. **开发日志量化**：每次开发日志除了记录做了什么，加一行量化数据（代码行数、验证结果、耗时对比）。

## 下次 session 计划（2026-04-17 末次更新）

**当前进度**（2026-04-17 收工，S6+S7 系列 10 commit）：
- Phase 1 Reads-based **7/7 完成**
- Phase 2 MAG-based **5/5 完成**（mag_quality / pathway / gene_profile /
  mag_heatmap / **network**）—— 全部按**统一 4 层参数**重构
- Phase 3 循环图 **S1→S4 全部完成**
- 测试 **245/245 全绿**
- 分析图表累计 **12 种**（7 Reads + 4 MAG + 1 Cycle）+ 网络 Gephi 辅助
- **5 张 MAG 图完全视觉统一**：Genus 标签规则（含 Family fallback）/ 门彩条
  位置 / 门图例 / 共享参数（filter_mode / row_order / top_n_by / max_mags）
- **v0.5 内部测试版就绪**

---

### 备选路径（按优先级）

1. **🟥 S8-ux 新手落地包（~10h）⭐ 最高优先** — 数据准备指南（覆盖 MajorBio/
   诺禾/BGI + CoverM/HUMAnN/eggNOG/GTDB-Tk/CheckM2/QIIME2 等上游工具） +
   "我该用哪张图" 研究问题向导 + 每图"如何解读"expander + 文件→分析反向
   索引（文件管理页一键提示可跑什么分析）。**详见 2026-04-17 阶段性自检 v2**
2. **S4.5 HTML 交互导出（~10-15h, 论文 SI 杀手锏）** — D3.js 独立 HTML 嵌入 cycle_data + hypothesis_score JSON，审稿人可交互
3. **S8-ui 导出中心统一重构（~3-4h）** — 把所有导出（bundle / 图表 / .py 脚本）汇集到"导出中心"页
4. **装 R 做 EnvMeta vs 原脚本侧侧对比** — 论文 Methods 关键证据，独立时间做（需装 R 环境）
5. **English README** — iMeta 要求，~2h
6. **S9 论文 Methods 起草** — 素材已齐（S1 去偏 + S2 置换 + S3+S3.5 评分器 + S4 Bundle + S6/S7 全 12 图）

### 已推迟（明确）

- **KB version 字段升 v2.0**：`elements.json` 里 `"version": "1.1"` 与实际含 v2.0 字段不匹配；独立小 session 修
- **S3.5-doc DOI 验证**：WebFetch 当时被 firewall 挡，写进 YAML 的 6 个 DOI 待用户手动过一遍
- **S8 插件框架**：论文接收后再做（CLAUDE.md 明确承诺）

### 里程碑

- ✅ Phase 1 全部完成（7/7 Reads-based）
- ✅ Phase 2 全部完成（5/5 MAG-based）
- ✅ Phase 3 核心功能全部完成（循环图 + 假说评分 + Bundle）
- ✅ **v0.5 内部测试版就绪**
- 📄 论文 Methods 可起草

## 开发日志

> 每次 session 结束前更新此区块。新对话开始时 Claude Code 自动读取，了解当前进度。

### 2026-04-17（**阶段性自检 v2** — 新手用户视角 + 开源免费定位 + 上游工具映射）

S7 后第二次全局自检。上次（2026-04-19 第一版自检）聚焦功能/竞争/审稿人
三角，本次新增 **新手用户视角**（CLAUDE.md 原承诺但从未具体审视）+ 用户
主动提的**开源免费定位** + **上游工具输出映射清单**，作为投稿前的补漏。

#### 1. 功能 + 竞争 + 审稿人（与 v1 比较）

- 功能：Phase 2 5/5 全完 + 12 图 + 循环图 + 假说评分 + Bundle；🟥 投稿
  阻塞项未变（R 对照 / English README / Methods / Zenodo / 二数据集）
- 竞争：三件套矩阵**升级为四件套**（新加"开源免费"列）—— 见下文
- 审稿人：接收率估计未变（现在 30-40% → 补 🟥 60-75% → 加 🟨 75-85%）

#### 2. 🆕 开源免费 vs 测序公司（用户新提）

环境微生物课题组常见痛点：**测序公司小众分析要加钱 / 有些根本不给做**。
价格现实（2026 年行情）：

| 分析 | 诺禾 / 美吉 / BGI / 基迪奥 | EnvMeta |
|---|---|---|
| 宏基因组测序（基础）| ~4500 元/样本 | — |
| α/β 多样性 / 堆叠图 / LEfSe 基础 | 含在基础套餐 | 免费 |
| MAG 组装 + 注释 | 加钱 5k-30k 元/项目 | 免费（接入上游输出）|
| 共现网络 + keystone 识别 | 个性化加价 or 不做 | 免费（Gephi 辅助）|
| 通路完整度 / 元素循环推断 | **几乎无公司提供** | **独有** |
| 假说评分 + Bundle 复现 | **无** | **独有** |
| 任意参数调整 | ❌ 参数锁死 | ✅ 全可调 |

**EnvMeta 四件套矩阵**：

| 工具 | Cycle | 假说评分 | Bundle | **开源免费** |
|---|---|---|---|---|
| **EnvMeta** | ✅ | ✅ | ✅ | ✅✅ |
| 测序公司云平台 | ❌ | ❌ | ❌ | ❌（付费）|
| Shiny-phyloseq / plotmicrobiome | ❌ | ❌ | ❌ | ✅ |
| Anvi'o | 部分 | ❌ | ❌ | ✅ |

> **关键洞察**：真实对手不是 Shiny-phyloseq，是**测序公司云平台**（MajorBio /
> 诺禾 / BGI / Wekemo）。云平台新手最好用但**结果质量不可控 + 小众分析要钱**。
> EnvMeta 定位 = "**比云平台更灵活 + 比 Anvi'o 更易用**" 的中间层。

#### 3. 🆕 新手用户视角（iMeta 论文核心叙事）

**真实场景**：环境微生物博士生，送样到测序公司，收到一堆表格，打开 EnvMeta，
然后……

**✅ 当前给新手的帮助**（已落地）：
- 自动文件识别（11 FileType + Gephi 刚加）
- 默认参数即可出图（5 张 MAG 图共享）
- 一键复现 `.py` 脚本 / Gephi 使用指南 / Help 悬停文字 / Genus 标签

**❌ 新手还缺的**（按优先级）：

- 🔴 **Tier 1 — 阻塞新手上手**：
  1. **上游工具输出 → EnvMeta 格式映射指南**（覆盖测序公司包 + 自跑上游）
  2. **"我应该用哪张图？" 决策向导**（研究问题 → 推荐分析）
  3. **每图"如何解读"expander**（参考假说评分器的 9 档解读模式）
- 🟡 **Tier 2**：样例数据一键加载 / 参数推荐对话（小样本时放宽阈值）/
  中文术语词汇表 / 常见错误 FAQ
- 🟢 **Tier 3**：云平台对比说明 / 预设方案（毕业论文版 / SCI 版）/ 视频教程

#### 4. 🆕 上游工具输出映射清单（扩展 Tier 1-1）

原来只想覆盖测序公司输出，实际要覆盖**两类用户**：
- 类型 A（非生信）：测序公司压缩包 → 不知道哪个文件对应 EnvMeta 输入
- 类型 B（半生信）：自己跑完上游分析 → 一堆输出不知道哪个接哪个

→ `docs/data_preparation_zh.md` 应含：

| 上游工具 | 典型输出 | EnvMeta 输入 |
|---|---|---|
| CoverM / Bowtie2 | abundance.tsv | `ABUNDANCE_WIDE` |
| MetaPhlAn4 | `*_bugs_list.tsv` | `ABUNDANCE_WIDE` |
| HUMAnN3 | `genefamilies.tsv` / `pathabundance.tsv` | `KO_ABUNDANCE_WIDE` |
| eggNOG-mapper | `*.emapper.annotations` | `KO_ANNOTATION_LONG` |
| DRAM | `annotations.tsv` / `distill/metabolism.tsv` | `KO_ANNOTATION_LONG` |
| GTDB-Tk | `gtdbtk.bac120.summary.tsv` | `MAG_TAXONOMY` |
| CheckM2 | `quality_report.tsv` | `CHECKM_QUALITY` |
| KofamScan | `*.mapper.tsv` | `KO_ANNOTATION_LONG` |
| QIIME2 | `alpha.tsv` / `beta-distance.tsv` | 多种 |
| Kraken2 + Bracken | `*_bracken.tsv` | `ABUNDANCE_WIDE` |
| iCAMP / SpiecEasi | 共现矩阵 | `GEPHI_NODES/EDGES` |

#### 5. 🆕 产品核心叙事（用户一句话点题）

> **"文件识别功能 = 解决 '文件到手不知道能干什么' 的痛点"**

这句话是 EnvMeta **跨越生信门槛的核心卖点**。建议写进：
- README 首段
- 论文 Abstract / Introduction
- 官网 hero slogan（以后有的话）

对应 UX 增强：文件管理页标完文件类型后，**反向提示可跑哪些分析**（~1h 工作量）：
```
📄 abundance.tsv — 📊 abundance (wide)
   → 可用于：堆叠图 / PCoA / RDA / MAG 丰度热图 / 共现网络
```

#### 6. 🆕 推荐路线（S8-ux 作为下一个 session 优先）

| 优先级 | 任务 | 工时 | 受众 |
|---|---|---|---|
| 🟥 **S8-ux 新手落地包** | 数据准备指南 + 问题向导 + 解读 expander | ~10h | 课题组非生信 |
| 🟥 文件→分析反向索引 | 文件管理页一键提示"可跑什么分析" | 1h | 所有用户 |
| 🟡 S4.5 HTML 交互导出 | SI 杀手锏 | 10-15h | 审稿人 |
| 🟡 R 侧侧对比 + English README | 投稿 🟥 清单 | ~1.5 天 | 审稿人 |
| 🟢 S9 论文 Methods 起草 | 素材已齐 | 1-2 周 | 审稿人 |

**核心判断**：**S8-ux 应插队到 S4.5 之前**。10h 就能把 EnvMeta 从
"开源工具"升级到"**非生信课题组真正能用的开源工具**"。这个定位差距比任何
新功能都值钱 —— 是 iMeta 论文 Results 里 user study 的先决条件。

---

### 2026-04-17（S7 — 共现网络 Gephi-prep 辅助 ⭐ Phase 2 5/5 全完）

- **目标**：Phase 2 最后一张"共现网络图"。不在 matplotlib 里画网络（Gephi
  更专业），EnvMeta = **Gephi 辅助工具**
- **三项交付**：
  1. **Degree vs Betweenness 散点图**（`envmeta/analysis/network.py` ~180 行）：
     keystone 深蓝 #1B3A5C + Genus 标注 + 可调阈值线
  2. **Gephi CSV 预处理**（`envmeta/tools/gephi_prep.py` ~120 行）：
     - `label_mode="keystone_only"` → 非 keystone Label 自动清空
       （省去用户在 Gephi 里一个个删标签）
     - `validate_gephi_format()` 校验必需列 / 引用完整性 / 重复边
  3. **Gephi 推荐参数指南**（app.py expander）：
     - 从用户 `.gephi` 项目提取：Fruchterman Reingold 区=10000/重力=10/速度=1
     - 节点大小 Degree 1-4 / 颜色 is_keystone 分区 / 标签 Arial Italic 32
- **S7-fix 4 轮**：
  - TYPE_BADGES 缺 GEPHI_NODES/EDGES → KeyError
  - detector.py 加 Gephi 文件识别规则（排在 KEYSTONE 之后避免误匹配）
  - gephi_prep 输出列名 `_x/_y` 冲突 → 直接用 nodes 自带 Genus 列
  - `annotate_taxonomy` 根治：merge 前 drop 已有列（一处改 5 模块受益）
  - 散点图标签截断 18→30 字符
- **测试**：235 → **245 全绿**（+8 network + 2 code_generator auto-parametrize）
- **里程碑**：**Phase 2 MAG-based 5/5 全完 → v0.5 内部测试版就绪**

---

### 2026-04-17（S6 系列收口 — S6-fix + S6-fix2 + S6-fix3 三轮用户反馈打磨）

S6 mag_heatmap 落地后用户浏览器实测 → 提出 13 点细节问题，一天内
3 轮迭代全部收口。**测试 215 → 235 全绿**。

#### 三轮 commit 时间线

| # | commit | 触发反馈 | 交付 |
|---|---|---|---|
| 1 | `fda28d5` **S6-fix** | alpha.txt 被误选当 MAG 丰度 / 上传 tax 仍显 Mx ID / 门彩条重叠无图例 / bubble 无 colorbar / gene_profile 缺"只看 keystone"过滤 / phylum_then_count 含义晦涩 | 文件自识别 (`_first_of`) / Genus 标签 / gridspec 布局 / Phylum 图例 / bubble colorbar / gene_profile filter_mode |
| 2 | `48e7f97` **S6-fix2** | 剩余 MAG 图（pathway / gene_profile / mag_quality）**未自动识别文件** + 通路完整度**无 Genus** + 基因谱**单色配色不明显** + phylum_then_XX 意义弱 + 想要**所有参数统一** | 抽出 `_mag_common.py` 共享模块（~260 行）；4 张 MAG 图对齐到**统一 4 层参数**（Layer 1-3 完全相同 + Layer 4 特有）；gene_profile 默认 cmap 改 **viridis** |
| 3 | `6ddcde7` **S6-fix3** | Mx_All_153 未显 Genus 原因？/ max_mags 语义？/ 基因谱 0 值占大部分视觉 | Family fallback（`JACPRB01 sp. Mx_153`）/ max_mags help 文案改"最终展示上限"/ gene_profile 加 **blank_zeros**（0 留白）+ **sort_ko_by_coverage** 选项 |

#### 架构成果：`envmeta/analysis/_mag_common.py` 新增（S6-fix2）

把 4 张 MAG 图的共性抽出来：

- `PHYLUM_COLORS` / `GROUP_COLORS` 共享常量
- `annotate_taxonomy()`: GTDB classification →
  Phylum / Family / Genus / Species / label 五列（label 用
  4 档 fallback：`Genus species` → `Genus sp. Mx_XX` → `Family sp. Mx_XX` → `Mx_XX`）
- `apply_filter_mode()`: **all / top_n / keystone_only / top_plus_keystone**
  统一子集过滤；score 来源可以是矩阵（mag_heatmap 丰度）或列
  （pathway abundance_mean / gene_profile gene_count）
- `order_rows()`: **phylum_cluster / metric_desc / abundance** 统一行排序；
  替代各图自己的 `phylum_then_total` / `phylum_then_count` 等晦涩名字
- `draw_phylum_bar()` / `draw_phylum_legend()` / `draw_keystone_note()`:
  统一左侧彩条 + 右侧图例区

#### 统一参数面板（4 张 MAG 图完全一致）

| 层 | 参数 | 说明 |
|---|---|---|
| L1 子集过滤 | filter_mode / top_n_by / top_n_count / max_mags | 所有 4 图同一套 widget |
| L2 视觉 | highlight_keystones / show_phylum_bar / show_phylum_legend + 画布 | 同上 |
| L3 行排序 | row_order / linkage_method | mag_heatmap / pathway / gene_profile 共享；mag_quality 不用（散点图） |
| L4 特有 | 各图独有（质量阈值 / 三段配色 / bubble 样式 / blank_zeros 等）| 按图变化 |

#### 向后兼容（所有老 bundle / 老调用不破）

hypothesis.py / bundle / cycle / 历史脚本写的老参数名（`sort_by=phylum_then_total`、
`top_n=30`、`selection_by=mean`、`cluster_rows=True`、`top_abundance_n`）
全部通过 `_normalize_deprecated_params()` 自动迁移到新名字，
旧测试与 paper/bundles 里的老 zip 仍能加载。

#### 真实数据验证（`paper/benchmarks/validation/mag_heatmap/`）

同一组数据（168 MAG × 10 sample + 14 keystone）：
- Top-30 mean 模式 Top-1 **Mx_All_102 (Chloroflexota)** mean=0.698%
- Top-30 variance 模式 Top-1 同样 **Mx_All_102** → 高丰度+组间分化大
- 5/14 keystone 进 Top-30 = 35.7% 重叠
- Mx_All_153（只到 f__JACPRB01）现在标签 **`JACPRB01 sp. Mx_153`**，
  保留科级信息而非纯 MAG id

#### 关键学习

1. **UI 大重构必须靠"先答疑 + 出参数表 + AskUser 定关键决策"三段式**：
   S6-fix2 用户明确要求"先为我设计每个功能侧边的可调参数"，不立刻写
   代码；AskUserQuestion 问 4 个关键决策后用户全部选推荐方案 → 一口气
   完成 4 图统一而不用反复改
2. **0 值的视觉处理比 cmap 选型更重要**：用户第一次反馈 cmap 单色
   不明显，但换成 viridis 后用户又发现"深紫满屏"——根本问题是
   viridis 的 0 值 ≠ 白色。`cmap.set_bad('white') + masked array`
   让"缺失 vs 低表达"一眼分
3. **向后兼容是架构自由度的代价**：为了保留 hypothesis.py 导入
   `pathway.PHYLUM_COLORS` + paper/bundles 里的老 YAML 参数名，
   每个模块都加了 `_normalize_deprecated_params()`；小成本换大收益
   （所有老 bundle 可加载 + 所有老测试不破）

---

### 2026-04-20（S6 — MAG 丰度热图 ⭐ Phase 2 达 4/5 + 12 图完整）

- **目标**：把 `scripts/python/07_MAG_abundance_heatmap.py`（544 行硬编码
  样本/组/配色）收敛为领域中立的 `envmeta/analysis/mag_heatmap.py`，
  让 EnvMeta 达到「12 图完整」投稿就绪里程碑
- **核心交付**：
  - `envmeta/analysis/mag_heatmap.py`（~220 行）：
    - `analyze(abundance_df, taxonomy_df, keystone_df, metadata_df, params)` → AnalysisResult
    - 13 个 DEFAULTS：top_n / selection_by(mean|sum|variance) /
      log_transform / cluster_rows / cluster_cols / cluster_within_phylum /
      linkage_method / color_breakpoints / show_phylum_bar / show_group_bar /
      highlight_keystones / width_mm / height_mm
    - `_build_tricolor_cmap()` — 保留原脚本的 Blues→YlGn→YlOrRd 三段非线性
      配色（适配少数极高 + 多数低值的长尾分布）
    - metadata 驱动样本排序（按 CK/A/B 组分组再组内保序）—— 替代原脚本
      的 `SAMPLE_ORDER` 硬编码
    - 门内聚类（`cluster_within_phylum=True`）保持门分块视觉一致性
  - `tests/test_mag_heatmap.py` **+12 case**（超额 6 case 目标）：
    smoke / top_n / selection_by 降序 / cluster 确定性 / no-tax / no-md /
    phylum+group bar 关闭 / color_breakpoints 自定义 / 极端 breakpoints /
    empty / missing sample_cols / PNG+PDF 导出
  - `app.py` 新页面「MAG 丰度热图」：4 文件上传（ab 必填 + tax/ks/md 可选）
    + 13 个参数滑块 / 下拉 / 勾选框 + 4 键下载（PNG/PDF/TSV + SVG/TIFF
    矢量）+ `_reproduce_button` 一键生成 .py 复现脚本
  - `envmeta/export/code_generator.py` SUPPORTED **10→11**，加
    `_tpl_mag_heatmap`（自动获得 +2 parametrize case = 37 code_generator 测试）
  - `paper/benchmarks/validation/mag_heatmap/`：`_run.py` + `top30_EN.pdf` +
    `top30_stats.tsv` + variance 模式副图 + README
  - `paper/benchmarks/time_comparison.md` 补行（544 行 / 90 min vs 3 点击 / 5 s）
- **真实数据亮点**（168 MAG × 10 样本 + 14 keystone）：
  - Top-30（mean 模式）含 **5/14 keystone**（35.7%）= 高丰度 ∩ 关键物种交集
  - **8 门** 分布；Pseudomonadota (11) / Chloroflexota (8) / Acidobacteriota (6)
  - **Top-1 `Mx_All_102`**（Chloroflexota）在 mean + variance 两模式均第一
    → 高丰度**且**组间分化大，机制验证首选候选
  - variance 模式与 mean 模式的 Top-30 不完全重叠 → 互补：mean=核心菌群，
    variance=处理响应敏感菌群
- **架构要点**：
  - 复用 `pathway.PHYLUM_COLORS` / `_extract_phylum` / `_mag_col`（避免重复）
  - 组色带颜色继承 `cycle_diagram` 的 CK/A/B 标准色（#4DAF4A/#377EB8/#E41A1C）
    —— 跨模块视觉一致性
  - 默认 `cluster_within_phylum=True`：先按门分组，门内聚类 → 行视觉上
    仍按门分块，但门内按相似性排序（原脚本同策略）
  - `color_breakpoints=(0.2, 0.5)` 元组化（原脚本硬编码），UI 可滑块调
- **测试**：215 → **227 全绿**（+12 mag_heatmap + code_generator 自动 +2 parametrize）
- **里程碑**：
  - Phase 2 MAG-based **3/5 → 4/5**
  - 分析图表 **11 → 12 种**（7 Reads + 4 MAG + 1 Cycle）——「12 图完整」✅
  - 投稿就绪清单 🟥 减一条（S6/S7 二合一只剩 S7）
- **下一步**：**S7 network Gephi-prep（~4h）** 收尾 Phase 2 → **v0.5 内部测试版**

---

### 2026-04-19（S4 — Fork Bundle 导出 ⭐ L4 层落地）

- **目标**：把 CLAUDE.md 承诺的"论文-EnvMeta 绑定发布协议"做成实际工具。
  读者/审稿人 fork 一个 zip 即可完整复现一篇论文的 EnvMeta 分析配置
- **Bundle 结构**：manifest + kb (elements.json + 可选 kegg_snapshot) +
  hypotheses (可多个 YAML) + config + 可选 README
- **不打包原始数据**：manifest 里留 paper_doi 指针
- **交付**：
  - `envmeta/tools/bundle.py` (~370 行): `create_bundle / load_bundle /
    inspect_bundle` + `BundleContents` dataclass + manifest schema 校验
  - 2 个新 CLI: `envmeta bundle-create` + `envmeta bundle-inspect`
  - app.py 循环图页 "📦 Fork Bundle" expander：左列加载 + 右列导出当前状态
  - `paper/bundles/arsenic_steel_slag_bundle.zip` (14.6 KB) —— 真实示例
  - `paper/bundles/README.md` —— "如何为你的论文打 bundle" 手册
  - +10 tests（round-trip / 多 hypothesis / 注释保真 / manifest 校验 /
    inspect / version mismatch / 错误路径）
- **版本兼容策略**：记录 envmeta_version；加载时 mismatch 只 warning 不 fail
  （向后容忍让老 bundle 可用）
- **注释保真**：hypothesis YAML 用 `zip.write()` 保留原文（含 S3.5-doc 写入的
  理论说明 + DOI 文献），不是 parse 后重 dump
- **与 S3/S3.5 衔接**：bundle 里的 YAML 可以直接拿到假说评分器里跑；
  "paper-bundled EnvMeta" 叙事完整闭环（L1-L2-L4）
- **测试**：204 → **214 全绿**（+10 bundle）
- **下一步**：S6 mag_heatmap（~3h 补 Phase 2）或 S4.5 HTML 交互导出（~10-15h
  SI 杀手锏）

---

### 2026-04-19（S3.5 — 假说评分器 v2：null 检验 + weight 敏感度 + required + validate CLI + group_contrast）

- **目标**：回应 reviewer 两大核心质疑（加权和是 ad-hoc / 阈值武断），把 S3 从
  "scorecard v1" 升到 "scorecard v2 with robustness indicators"；顺便补 group_contrast
- **5 项产出**（一次 session 全部封死）：
  1. **Null comparison**（排列检验 Fisher 1935 风格）：`HypothesisScore.null_p` +
     `null_p_samples`；n_non_skipped<3 或 weight/score 同质时退化为 None
  2. **Weight sensitivity** OAT ±20%：`weight_robust: bool | None` +
     `weight_sensitivity_rows: list[dict]`
  3. **required claim 硬否决**：YAML `required: true` 字段；任一必要条件失败 →
     label 强制 `insufficient` + `veto_reasons` 列原因。overall_score 仍透明显示
     + params 保留 `base_label_before_veto`
  4. **`envmeta hypothesis-validate` CLI**：检查 YAML 对 KB 的引用是否存在
     （pathway 名 / coupling species），给出相近候选建议
  5. **第 5 类 claim `group_contrast`**：`high_group/low_group` total_contribution
     比值 ≥ min_ratio → satisfied。`score(..., compare_df=)` 可选注入；
     app.py 检测到 group_contrast claim 自动调 compare_groups
- **架构**：hypothesis.py 不 import cycle_compare（避免循环依赖），compare_df
  作为参数注入；inference.py 完全不动
- **理论基础文献写入 YAML 顶部**（Fisher 1935 新加一条）：`required` 对应
  Bradford Hill 必要条件；null_p 对应 Fisher 排列检验；权重比例对应 MCDA
- **真实数据结果**（arsenic_steel_slag.yaml 9 claims vs 三组）：
  - **CK**: 0.901 strong | null_p=0.296（随机样）| robust ✅ | veto 0
  - **A**:  0.890 strong | null_p=0.764（最随机）| robust ✅ | veto 0
  - **B**:  **1.000 strong** | null_p=N/A（退化）| robust ✅ | veto 0
  - null_p 让 "B 组独特性" 看得见：CK/A null_p 接近随机，说明 YAML 里的
    假说对它们不特异；B 组全通过（退化为 None）+ label strong + robust = 真正支持
- **测试 189 → 202（+13）**：
  - null 3 (basic / too few / disable)
  - sensitivity 2 (present / stable)
  - required 3 (satisfied no-veto / unsatisfied veto / schema reject non-bool)
  - group_contrast 3 (satisfied / skipped no-compare_df / schema accept)
  - validate CLI 2 (accept good / warn unknown pathway)
- **commit**（待）：feat(hypothesis) S3.5 + docs(yaml) required/group_contrast
- **论文 Methods 可引用**：
  > "We report three complementary robustness indicators alongside the headline
  > score: (1) a permutation null p-value [Fisher 1935] quantifying the
  > probability of reaching this score by chance; (2) one-at-a-time weight
  > sensitivity analysis [Saltelli 2008] reporting label stability under ±20%
  > perturbation; (3) a Bradford-Hill-inspired `required` veto mechanism [Hill 1965]
  > for user-designated necessary conditions."
- **下一步**：S4 Fork Bundle (~4h) 或 S6 mag_heatmap (~3h, Phase 2 补齐)

---

### 2026-04-19（**阶段性自检** — 产品定位 / 竞争分析 / 投稿就绪度 三角审视）

在 S5/S6 启动前做一次全局自检，避免闭门造车。**不动代码**，结论写入此处供后续
session 随时引用。

#### 1. 功能 vs 初衷对齐

**✅ 已精进（超出预期）**：
- L1 循环图：从"2×2 抽象柱图"升级到 Mockup 10 发表级级联 + KEGG-driven KB + 跨元素耦合
- L2 假说评分器：原承诺"用户自带 YAML → strong/none"→ 现在 5 类 claim + 3 可信度指标
  (null_p / weight_robust / required-veto) + 9 档 UI 解读 + 跨组表 + CLI 校验 + 7 篇理论文献
- L4 Fork Bundle：论文-EnvMeta 绑定发布协议如期落地
- 架构解耦：hypothesis.py 不依赖 inference.py / bundle.py 不依赖 cycle_compare —— 单向依赖图
- 描述性不断言：S1 去偏反思彻底，评分器明确声明 "NOT hypothesis testing"

**⚠️ 未实现 / 落后**：
| 项 | 阻塞等级 | 工时 |
|---|---|---|
| S6 mag_heatmap + S7 network（12/12 图）| 🟥 阻塞投稿 | 7h |
| R 侧侧对比验证 | 🟥 阻塞投稿 | 1 天 |
| English README | 🟥 iMeta 要求 | 2h |
| tool_comparison.md 填表 | 🟥 审稿必挑 | 4h |
| 示例数据 Zenodo DOI | 🟥 投稿必备 | 2h |
| LICENSE 文件 | 🟥 合规 | 10 min |
| Methods + Results + Discussion 起草 | 🟥 论文本体 | 1-2 周 |
| User study（5-10 人 SUS）| 🟨 推荐 | 1 周 |
| S4.5 HTML 交互导出（SI）| 🟨 推荐 | 10-15h |
| 第二个数据集复现 | 🟨 推荐 | 2 天 |
| L3 插件框架 / L5 KB 工具 / Docker | 🟩 推迟 | — |

**🟡 可能冗余（已审视，保留为主）**：
- CLAUDE.md 日志超长：对 AI 上下文恢复必要，对人类阅读不友好 → 可抽"精简版发版日志"（不紧急）
- `scripts/` 下原论文 R/Python：保留作 validation reference，可整理入 `scripts/legacy/`

**结论：无实质偏离初衷**。L1/L2/L4 三层超预期，缺口全在"投稿前行政事务"。

---

#### 2. 竞争分析 — 三件套矩阵

| 工具 | Cycle inference | 假说评分 | Bundle / 论文绑定 |
|---|---|---|---|
| **EnvMeta** | ✅ **独** | ✅ **独** | ✅ **独** |
| Shiny-phyloseq | ❌ | ❌ | ❌ |
| plotmicrobiome | ❌ | ❌ | ❌ |
| Anvi'o | 部分 | ❌ | ❌ |
| iPath3 | 静态 | ❌ | ❌ |
| BioCyc（商用）| ✅ | ❌ | ❌ |

三件套同时具备 = **业界空白**。三个卖点相互强化（bundle 装的就是 KB+YAML，闭环）。

**相对劣势**：无 R 集成 / 用户基数小 / KB 当前聚焦 As/N/S/Fe / 依赖 GTDB+KEGG 上游 /
无在线 demo / 中文优先（国际化滞后）/ 单数据集验证。

---

#### 3. 审稿人视角 + 投稿就绪度

**会表扬的 8 点**：差异化独份；理论扎实；架构透明；诚实校准（不假装 hypothesis
testing）；null_p 真有区分力；测试 215；Fork Bundle 方法学创新；S1 去偏反思彻底。

**会质疑的 10 点**（按杀伤力排序）：
1. 只在单数据集验证 → 补公开数据集复现
2. 无 R 对照 → 装 R 做侧侧 PDF 对比
3. 无 user study 证据 → 5-10 人 SUS
4. 12 图不全 → S6 + S7
5. KEGG license 声明 → 加 NOTICE
6. 权重主观性 → 已引 MCDA 回应
7. 大数据集性能 benchmark 缺 → 跑 500+ MAG
8. YAML 写法门槛 → 已有 validate CLI
9. 其他元素循环通用性 → 补 C 或 Hg KB 示范
10. LICENSE 未定 → MIT / Apache-2.0

**接收率估计**：
- 现在投稿：**30-40%**（差异化强但证据单薄）
- 补完 🟥 清单后：**60-75%**（三大硬伤解决）
- 加 🟨 清单：**75-85%**（方法学论文稳妥）

---

#### 4. 推荐冲刺路线（4-6 周到投稿）

```
Week 1:  S6 mag_heatmap + S7 network + LICENSE + tool_comparison.md
Week 2:  R 对照（6 图）+ 第二数据集复现
Week 3:  English README + Zenodo DOI + 示例数据包 + Figures 整理
Week 4-5: Methods + Results + Discussion 起草
Week 6:  User study + S4.5 HTML SI + Cover Letter
```

**核心判断**：工具本身已达投稿级别，剩余都是"科研行政事务"而非技术创新。

**即刻下一步**：**先 S6**（~3h）收紧 11/12 → 12/12，再开其他。

---

### 2026-04-19（**Today's Session 总结** — 全日 8 commit / L2 假说评分器 v2 + L4 Fork Bundle 落地）

今日一日跑完 **8 个 commit / 测试 176 → 215（+39）**，L1-L2-L4 三层架构闭环，
循环图 + 假说评分器 + 绑定发布协议全部可发表级。

#### 时间线（commit 顺序）

| # | commit | 主题 | 测试变化 |
|---|---|---|---|
| 1 | `27c32d6` | fix(cycle) S2.5-14 multi-chain 耦合锚点回归（As(V)↔Fe(III) 棕线消失）| +1 |
| 2 | `3701f28` | docs: S2.5-14 session 日志 | — |
| 3 | `ff4b1b5` | feat(hypothesis) **S3** — 机制假说 YAML 评分器（L2 层首次落地）| +12 |
| 4 | `d35b7d2` | docs: 2026-04-19 阶段性总结（v1）| — |
| 5 | `13a1681` | docs(hypothesis) **S3.5-doc** — YAML 自带理论依据 + 6 篇 DOI 文献 | — |
| 6 | `8039f13` | feat(hypothesis) **S3.5** — null 检验 + weight 敏感度 + required + validate CLI + group_contrast | +13 |
| 7 | `2b0241c` | feat(app) **S3.5-ui** — 评分解读一句话 + 三指标参考卡 | — |
| 8 | `e075c6a` | feat(hypothesis) 跨组评分对比表（论文 Results 主轴一键呈现）| +2 |
| 9 | `1f85db5` | fix(app) 假说评分 UI 3 点打磨（not-allowed 遮挡 / 普适表述 / 减冗余）| — |
| 10 | `d495c3b` | feat(bundle) **S4** — Fork Bundle 导出（L4 层论文绑定发布协议）| +10 |
| 11 | `447546b` | fix(bundle) app UI 导出补齐 KEGG 快照 + README + 清空字段 | +1 |

#### 四大里程碑

**S3 / S3.5（L2 层 — 机制假说 YAML 评分器 v2）**：
- 5 类 claim：pathway_active / coupling_possible / env_correlation /
  keystone_in_pathway / group_contrast
- 3 个可信度指标：null_p（Fisher 1935 排列）/ weight_robust（OAT ±20%）/
  required-veto（Bradford Hill）
- UI 综合解读（9 档判断）+ 跨组对比表（CK/A/B 一屏呈现）
- `envmeta hypothesis-validate` CLI 校验 YAML 对 KB 的引用
- **hypothesis.py 不 import inference.py**，架构性解耦 —— 论文可
  诚实声明 "hypothesis-agnostic inference + separate evaluator"

**S3.5-doc（YAML 理论武装）**：
- `arsenic_steel_slag.yaml` 顶部 80 行教程 + 7 篇带 DOI 文献（MCDA
  Keeney&Raiffa 1993 / Belton&Stewart 2002 / Bradford Hill 1965 /
  WoE Suter&Cormier 2011 / Linkov 2009 / Rhomberg 2013 / Fisher 1935）
- 每条 claim 带 Bradford Hill 维度标注（strength / specificity /
  biological plausibility / coherence / gradient）
- README 精简，YAML 成自包含教程

**S4（L4 层 — Fork Bundle 绑定发布协议）**：
- `envmeta/tools/bundle.py` ~370 行：create/load/inspect
- `envmeta bundle-create` + `envmeta bundle-inspect` 双 CLI
- app.py 循环图页双列 UI（加载 / 导出）
- `paper/bundles/arsenic_steel_slag_bundle.zip` 14.6 KB 官方示例
- 读者/审稿人 fork zip → 一键复现论文分析

**S2.5-14（收尾）**：
- multi-chain cell 的 substrate_pos_map / product_pos_map 修正
- 4/4 KB 耦合线（含 As(V)↔Fe(III) 棕色）CK/A/B 全部回归

#### 真实数据亮点（论文 Results 一张表）

`arsenic_steel_slag.yaml` 9 claims vs 三组实测（S3.5 版）：

| 组 | overall | label | null_p | weight_robust | 判读 |
|---|---|---|---|---|---|
| CK | 0.901 | strong | 0.297 | ✅ | strong 但不特异（通过率运气主导）|
| A  | 0.890 | strong | 0.770 | ✅ | strong 但不特异（最随机）|
| **B** | **1.000** | **strong** | **N/A** | ✅ | **最强支持（退化式 N/A + robust = 真 strong）**|

null_p 让 "B 组独特性"看得见：CK/A 的 strong label 其实是"通过率幸运"，
只有 B 组是权重设计与数据贯穿的真 strong。

#### 测试累计

**176 → 215（+39）**：
- +1 multi-chain 耦合回归（S2.5-14）
- +12 hypothesis S3 单元（load/4 类 claim/overall 聚合）
- +13 S3.5 可信度指标（null×3 / sensitivity×2 / required×3 /
  group_contrast×3 / validate×2）
- +2 跨组评分 `score_by_groups`
- +10 bundle round-trip / manifest / inspect / 错误路径
- +1 manifest 空字段省略

#### 交付文件（今日新增全清单）

| 类别 | 文件 |
|---|---|
| 源码 | `envmeta/geocycle/hypothesis.py`（+430 行）/ `envmeta/tools/bundle.py`（~370 行）/ `envmeta/tools/hypothesis_validator.py`（~170 行）/ `envmeta/analysis/hypothesis_compare.py`（~150 行）|
| CLI | `envmeta bundle-create` / `envmeta bundle-inspect` / `envmeta hypothesis-validate` |
| UI | app.py 「📦 Fork Bundle」expander + 「🧪 假说评分」9 档解读 + 跨组对比表 + 三指标速查 |
| 依赖 | `requirements.txt` +pyyaml |
| 测试 | `tests/test_hypothesis.py`（+27 case）/ `tests/test_bundle.py`（+11 case）|
| 论文示例 | `paper/hypotheses/arsenic_steel_slag.yaml`（9 claims + 7 文献）+ `paper/bundles/arsenic_steel_slag_bundle.zip` |
| 论文文档 | `paper/hypotheses/README.md` + `paper/bundles/README.md` |
| 归档 | `paper/benchmarks/validation/hypothesis/` (CK/A/B 三组 TSV+JSON+README) |

#### 论文 Methods 素材

> "EnvMeta couples hypothesis-agnostic cycle inference with an
> evidence-weighting scorecard [MCDA: Keeney & Raiffa 1993; Bradford Hill
> 1965]. The scorer reports three complementary robustness indicators:
> (i) a permutation null p-value [Fisher 1935] quantifying the probability
> of reaching the observed overall score by chance; (ii) one-at-a-time
> weight sensitivity analysis reporting label stability under ±20%
> weight perturbation; (iii) a Bradford-Hill-inspired `required` veto
> mechanism for user-designated necessary conditions. Analysis
> configurations (knowledge base + hypothesis YAML + parameters) can be
> packaged into a Fork Bundle (zip) for reviewer-reproducible delivery."

#### 今日关键学习

1. **UI 解读 ≫ 原始指标**：null_p=0.30 单看没感觉，配上"通过率运气主导"
   的一句话解读，审稿人立刻懂。S3.5-ui 是本日投入 / 产出比最高的改动
2. **架构解耦的最佳时机是"不写"**：hypothesis.py 不依赖 inference.py、
   bundle.py 不依赖 cycle_compare.py，单向依赖维持系统可演化
3. **Bundle 的价值在"不打包什么"**：不打包原始数据（体积/敏感），留 DOI
   指针；不加密 / 不签名（违反分布式 fork 原则）
4. **"required claim 硬否决" 的语义独立性**：与软惩罚（加权重）的关键
   区别是它有独立的语义位面，避免用户用权重调参蒙混过关
5. **退化式 N/A vs 无信号 N/A**：null_p=None 在 overall=1.0 时是"好"，
   在 claim<3 时是"没信号"—— 必须在 UI 层区分，不能只展示原始值

---

### 2026-04-19（S3 — 机制假说 YAML 评分器 ⭐ L2 层落地）

- **目标**：EnvMeta 推断是描述性无假说的，用户想"评估我的假说数据支不支持"需要
  独立入口。S3 填 L2 层（"用户自带假说 YAML → 证据评分"），与 inference 架构解耦
- **4 类 claim**（v1 范围）：
  - `pathway_active`：通路活跃 + completeness/contribution 阈值
  - `coupling_possible`：KB 有 species_a↔species_b 配对 AND 两端物种都被观测
  - `env_correlation`：(pathway, env) 符号方向 + confidence 达 min 阈
  - `keystone_in_pathway`：通路含 ≥N keystone contributor
- **聚合**：`overall = Σ(w·score)/Σ(w)`（仅非 skipped），4 档标签
  strong(≥0.75) / suggestive(≥0.40) / weak(>0) / insufficient
- **skipped 不扣分**：KB 没有的耦合 / 数据里没跑到的通路 → 从分母剔除，避免
  "用户写得不完整 → 分数被压低"的误解
- **交付**（commit 待）：
  - `envmeta/geocycle/hypothesis.py`（~430 行）：dataclass + load_hypothesis +
    score + 4 evaluator + 辅助（pathway 模糊匹配、_collect_observed_species）
  - `paper/hypotheses/arsenic_steel_slag.yaml`：8 claim 示例（对应论文研究假说）
  - `paper/hypotheses/README.md`：schema 说明 + 4 类 claim 语义 + 撰写 checklist
  - `tests/sample_data/sample_hypothesis.yaml` + `tests/test_hypothesis.py` (+12 case)
  - `app.py` 循环图页底部 expander「🧪 假说评分 (可选)」 + 示例下载 +
    上传 + 评分按钮 + TSV/JSON 下载 + 彩色标签徽章
  - `paper/benchmarks/validation/hypothesis/`：CK/A/B 三组实跑结果
- **真实数据结果**（arsenic_steel_slag.yaml vs 样本数据）：
  - **CK**: 0.879 / strong (5/6 sat, 2 skipped) — 背景铁砷活性已存在
  - **A**:  0.868 / strong (6/7 sat, 1 skipped) — 中等处理
  - **B**:  **1.000 / strong** (6/6 sat, 2 skipped) — 高钢渣组满分支持
  - 三组 skipped 主要是 CK/A 没同时跑出 Fe(III) + S-2 耦合锚点（生物学合理），
    skipped 不扣分正确体现
- **论文可写**："hypothesis-agnostic cycle inference + hypothesis-testing YAML
  evaluator separates evidence from interpretation"
- **测试**：176 → **189 全绿**（+13：12 hypothesis 新 case + test_overall_label
  edge + existing regression 不受影响）
- **下一步**：S4 Fork Bundle 或 S4.5 HTML 交互导出（见 CLAUDE.md 路线）

---

### 2026-04-19（S2.5-14 — 回归修复：multi-chain 细胞丢失化学物耦合锚点）

- **症状**：用户发现 B 组图中 `As(V)↔Fe(III)` 的棕色耦合线（v29 之前有）在
  S2.5-11 合并同-MAG 多通路后消失。另外 `Fe(II)↔S-2`（FeS 前驱）也未画出
- **根因**：MAG-merge（commit 4cf2c89）后 Fe 细胞合并 `fur/tonB`（substrate/product
  都是 None）+ `Fe transport`（Fe(III)→Fe_internal）。单值 `substrate_species`
  锚点取 `merged_genes[0].substrate` → None（fur 排在最前）。KB 的 `As(V)↔Fe(III)`
  耦合永远找不到 Fe 侧锚点 → 静默 drop
- **修复**（commit `27c32d6`，3 文件 / +103 / −17 行）：
  - `substrate_species` / `product_species` 改列表（收集所有 merged_genes 的物种）
  - 多链 cell 额外导出 `substrate_pos_map` / `product_pos_map`，按化学物名精确定位
  - 单链外部 substrate/product 跳过 None（调控基因前置场景）
  - `_anchor_fig_point` 加 `species` 参数；`_find_best_pair` 支持列表成员匹配
  - +1 回归测试 `test_multichain_cell_exposes_fe_iii_for_coupling`
- **验证**：CK/A/B × 2 ranking 6 组合，**4/4 KB 耦合线全部回归**；
  `hide_regulator_only_cells` 仍正常过滤；**176/176 全绿**
- **自检结论**：S2.5-1 ~ S2.5-13 全部 13 项功能在线；唯一漏网的就是这条
  耦合回归，本次修复锁死。可以进入 S3

---

### 2026-04-18（**Session 总结 — S2.5 系列全部收口（13 个子 session）⭐ 循环图可发表级**）

本阶段是 Phase 3 循环图从 v1 抽象柱图 → **Mockup 10 发表级合并细胞级联图**
的完整落地，并把 KB 从手写升级为 **KEGG-driven**。累计 **18 个 commit**，
全测试 **101 → 176 全绿（+75 新 case）**。

#### 交付（按 commit 顺序）

| # | commit | 内容 | +测试 |
|---|---|---|---|
| S2.5-1 | 5bb0438 | KB schema 扩展（substrate/product/schematic/couplings） | +3 |
| 2a | 716df23 | MAGContribution.genes KO 级 cascade 字段 | +3 |
| 2b | e8f70d9 | cell_renderer primitive（合并细胞绘制） | +7 |
| 2c | 2afa9df | 主 renderer v2 级联改造 | +3 |
| S2.5-3 | 151d104 | 化学物耦合连线（跨元素虚线 + 产物节点） | +4 |
| S2.5-4 | 40e2771 | 组选择下拉（CK/A/B 单组模式） | +4 |
| S2.5-5 | bd9011e | SVG + TIFF 导出（SCI / 毕业论文双预设） | +9 |
| S2.5-6 | 24d3cbd | 文件识别扩展（ko_long/keystone/mag_taxonomy）+ 自动匹配 + 字符修复 v1 | +5 |
| S2.5-7 | 59617ef | 字符排版 v2 + 多亚基并列 + 耦合线去叠 + 跨组对比 cycle_compare | +10 |
| S2.5-8 | 5c5a4b2 | 跨组最活 ★ + keystone ✦ 标注 | +3 |
| S2.5-9 | 8521f1f | MAG 选择判据 4 档可切换（abundance/completeness/keystone_only/keystone_priority） | +4 |
| S2.5-10a | d05d0e6 | KEGG 快照构建器 + 种子 KO 表 + 57-KO snapshot JSON |  |
| S2.5-10b/c | 58560e0 | `envmeta kb-build` CLI + KB v2.0 迁移（29 KO 获 complex 字段） | +7 |
| S2.5-10d | 64ef2fe | renderer 段内复合体分段 + 缺失通路显示 | +4 |
| post-10d | 895c598 | parallel_complex 标注 complex id + 断链 label 垂直堆叠 |  |
| bundle | 0523fb4 | 同-complex 基因合并成单大椭圆（统一布局） | +2 |
| chain | 00bab82 | 断链改多行垂直布局 + bundle label 不截断 | +2 |
| ext-h | 4153a0a | 多链 cell_h 自动扩大 + 外置 substrate/product |  |
| mag-merge | 4cf2c89 | 同 MAG 多通路合并为单 cell（多行渲染） |  |
| S2.5-13 | ee5f307 | Genus+species 标签 + 可选隐藏调控型 cell（fur/tonB） | +2 |

#### 架构里程碑

**数据层（KEGG-driven KB, S2.5-10）**：
- `scripts/build_kegg_snapshot.py` 联网抓 KEGG REST → `envmeta/geocycle/kegg_snapshot.json`
- `envmeta/tools/kb_builder.py` + `envmeta kb-build` CLI 一键生成 elements.json
- 手写 KB ⭐⭐⭐ → KEGG-driven ⭐⭐⭐⭐⭐ 权威性升级
- `complex` 字段（KEGG MODULE id）解决 sqr+Sox 混淆、narG/H/I/napA/B 显式标注

**视觉层（Mockup 10 落地）**：
- `cell_renderer.py` 三分枝：parallel_complex / all_same_intermediate / segmented
- 同-complex 多基因 → 单大 bundle 椭圆（含多行 wrap 不截断）
- 断链场景 → 多行垂直布局，每行独立外置 substrate/product
- 化学式全部 mathtext 粗体（SO₄²⁻、NO₃⁻ 正确下上标）
- 跨元素化学耦合：紫/棕/绿虚线 + 产物节点 + 去叠
- redox 类型产物节点：`As(III)→As(V) (ox)` 明示转化方向

**分析层（多视角 / 跨组）**：
- 4 档 MAG 排序（S2.5-9）：abundance / completeness / keystone_only / keystone_priority
- 跨组对比表（S2.5-7d）`compare_groups()`：回答"A/CK/B 差别在哪"
- 跨组最活 ★ + keystone ✦ 双标注（S2.5-8）
- 标签 Genus + species 消除同属不同种视觉歧义（S2.5-13）
- 可选 `hide_regulator_only_cells` 过滤 fur/tonB 纯调控（S2.5-13）

**产出层**：
- PNG / PDF / SVG（SCI 投稿） / TIFF 600dpi（毕业论文）4 格式
- 每分析页下载 `.py` 复现脚本（既有）
- 11 张 Mockup（01-10）+ 大量 validation PDF 归档到 paper/

#### 真实数据科学发现（可直接写论文）

1. **加钢渣（B 组）显著抬升砷代谢活性**：
   - Arsenate reduction total_contribution: CK=189 / A=237 / B=361（1.9×）
   - Ammonia oxidation：CK=0 / B=24.3（**22×**）
2. **Keystone species 在处理组成为功能主力**：
   - A 组 Nitrate reduction top MAG = Sulfuricaulis（keystone）
   - B 组 Arsenate reduction top MAG = SPCO01（keystone）
   - **CK 组 keystone 不承担任何 top 通路** —— 说明钢渣处理让关键物种从网络边缘走向功能中心
3. **KEGG 打包 vs 生物学真相被正确呈现**：
   - Sulfide oxidation 细胞现在清楚分成两行：sqr 路径 + Sox 复合体路径
   - Nitrate reduction 6 基因显示为 `narG/narH/narI/napA/napB/narB (M00529 complex, 6 subunits)` 单椭圆

#### 关键用户交互修正（反映在架构中）

按用户反馈反复迭代的关键设计点：
- **领域中立**：拒绝"自动检测研究主题"假命题（S1）
- **用户自带 KB**：KEGG-driven 减少 50%+ 扩展成本（S2.5-10）
- **MAG ID 而非 genus 分组**：同属不同菌株是**生物学上不同** MAG，显示加 `sp. Mx_XX` 消除视觉误解（S2.5-13）
- **描述性不断言**：cycle diagram 不用"drives"因果词；置换零假设区分 strong/suggestive/spurious?（S2）
- **发表级排版**：化学式粗体 + 下标上标；基因椭圆按名长自适应；长 bundle 双行不截断

#### 测试进展

- **101 → 176**（+75 case），全绿
- 新增覆盖：KB loader、cell_renderer 分支、chain 拆分、bundle label、
  KEGG snapshot、kb_builder、cycle_compare、figure_export 多格式、
  cross-group 标注、4 档 ranking、label Species 格式、regulator 过滤

#### CLI 新增

- `python -m envmeta kb-build --elements ... --preserve-from ... --output ...`
- `python scripts/build_kegg_snapshot.py --seed ... --output ...`

#### 归档产物

- `paper/benchmarks/validation/cycle_diagram/`：
  v1.pdf (原 S1 柱图) / v2.pdf (Mockup 10 级联) / v2_kegg.pdf / v2_bundle_{A,B,CK}.pdf /
  v2_chain_{A,B,CK}.pdf / v2_magmerge_{A,B,CK}.pdf / v2_extheight_{A,B,CK}.pdf /
  v2_labeled_B.pdf + v2_labeled_B_noreg.pdf
- `envmeta/geocycle/kegg_snapshot.json`（24 KB，57 KO + 12 module）
- `envmeta/geocycle/knowledge_base/elements.json` v2.0（手工 + KEGG 合并）
- `paper/benchmarks/validation/cycle_diagram/README.md` 扩展段：
  MAG selection criteria / 跨组差异叙述 / MAG 标签格式 / Fe uptake 无 I/O 解释

#### 下次 Session 起点

**S2.5 循环图全部收口**。剩余 Phase 3 增强可选路径：
1. **S3 机制 YAML 评分器**（5h，路线 B+ 核心）：用户自带假说 YAML
   → CycleData 证据评分 → strong/suggestive/none 标签
2. **S4 Fork Bundle 导出**（4h）：打包 elements.json + 用户自定义 +
   paper-specific 配置为 zip（"论文-EnvMeta 绑定发布"协议）
3. **S4.5 HTML 交互导出**（10-15h，路线 B+ 大卖点）：D3.js 独立 HTML
   嵌入 cycle_data JSON + 拖拽/悬停/缩放/SVG 导出
4. **S6 mag_heatmap**（3h）：补 Phase 2 最后 2 张图之一
5. **S7 network Gephi-prep**（4h）：共现网络 Gephi 兼容导出
6. **S9 论文 Methods 起草**（3h）：iMeta / Bioinformatics 方法学

建议：**先 S3**（机制 YAML），它能回答用户"如何评估我的砷修复假说"这类
具体问题，与已有 CycleData 推断直接对接；再 S4.5（HTML 导出）作为论文
SI 杀手锏。S6 / S7 可并入 Phase 2 的尾声。

---

### 2026-04-17（**Today's Session 总结 — 战略重构 + S1 + S2 完成**）

今日是 **Phase 3 v1 完成后的战略反思 + 深度重构日**。重点不在代码量，在**方向校准**：

#### 战略层面的 5 轮深度讨论

1. **确认偏差担忧**：前一版 README 写的"Gallionella drives Arsenate reduction"
   是我在解读数据时附和了用户的原论文假说。用户敏锐识别出这个风险，引发讨论：
   - EnvMeta 核心算法**实际上是无偏的**（无假说输入通道）
   - 但**叙述层面**（README、label）容易被上下文带偏
   - **结论**：工具要刻意"描述性而非断言" → S1 去偏

2. **竞争机制生成器（Option C）讨论** → **否决**：
   - 初案：预装机制目录，用户研究问题自动匹配
   - 用户质疑："env 列名 + KO 频度推不出砷修复"——**完全正确**
   - 讨论用户真正愿景：EnvMeta 应是**可定制框架**，不是预装机制库
   - **产品定位彻底重写** → 领域中立 + 用户自带 KB + Fork Bundle 发布

3. **市场/可行性诚实评估**：
   - 预期用户峰值：50-200 / 1-2 年（**不是上千**）
   - 主要受众：自己课题组 + 同行复现（每篇 7-35 人/次）+ 小众外部（5-50）
   - 博士论文章节 / 方法学 SCI 论文（iMeta / Bioinformatics）— **完美匹配**
   - 成为领域标准工具 — **不现实，不追求**

4. **KB 制作难度分层评估**：
   - 加通路：30 min ⭐
   - 加新元素循环：4-8h ⭐⭐⭐⭐
   - 换研究领域：1-2 天 ⭐⭐⭐⭐⭐
   - 机械成本可降 ~50%，**智力成本无法降低**（哪些 KO 属哪通路是专业判断）
   - EnvMeta 提供：schema 校验 + template + diff/merge CLI

5. **循环图可视化大改造**（Mockup 01 → 10 连续 5 轮迭代）：
   - M01-03：早期 2×2 柱图 → 用户："不够直观"
   - M04：FeOOH 中心嵌入细胞 → 用户："不能预设中心结构"
   - M06-07：数据驱动 2×2 网格 → 用户："化学流程没经过细胞"
   - M08-09：化学物穿越细胞 → 用户："多基因同一 MAG 要合并"
   - **M10 定型**：合并细胞 + 细胞内级联（SO₄²⁻→cysJ→SO₃²⁻→cysI→S-cys）+
     元素耦合连具体化学物（As(III)─H₂S→As₂S₃）+ 数据驱动布局
   - 用户评价："方向是对的，够指导手绘图"

6. **交互方案 A/B/C 评估**：
   - A：Plotly 替换 matplotlib（6-8h，in-app 悬停）
   - **B：独立 HTML 导出（10-15h，论文 SI 杀手锏）⭐ 采纳**
   - C：D3.js 真编辑器（25-35h，推迟到 Phase 4）

#### 代码层面交付

| S# | 内容 | commit |
|---|---|---|
| S1 | 去偏（完整矩阵 + 敏感度 + 描述语言）| a3cefa5 |
| S2 | 置换零假设 + 可信度标签 | 49062e8 |

#### 真实数据的硬证据（S2 置换检验结果）

对你论文数据的 16 条超阈值 env-pathway 相关：
- 🟢 **1 strong**：Ammonia oxidation ↔ Eh（ρ=0.835, perm_p<0.01, robust）
- 🟡 **10 suggestive**：Arsenate reduction / As transport 等
- 🔴 **5 spurious?**：S1 讨论时**预测**的"5 条共享 ρ=0.764 是共变信号"，
  **被置换检验证实**。论文 Methods 可直接引用这个去伪存真过程

#### Mockup 归档

`paper/figures/mockups/01-10*.png` 10 张 + 对应 3 个 gen_*.py 生成器脚本

#### CLAUDE.md 结构性更新

在"开发路线图"和"Backlog"之间新增 **"产品定位与核心设计决策"**章节：
- 定位：可定制环境生信分析框架 + 通用循环图推断引擎 + 论文-EnvMeta 绑定发布协议
- 5 项核心设计原则
- 5 种拒绝的设计（列明理由）
- 5 层架构（L1-L5）
- 市场校准 + KB 难度分层
- 当前执行路线：B+（27h 原版 → 62-72h 升级版）

#### 截图归档

- screenshot_10_cycle_diagram.jpeg（S1 版）
- screenshot_11_cycle_s2_confidence.jpeg（S2 版，带可信度徽章）

#### 测试进度

Session 开始：98/98 → S1 结束：101/101 → S2 结束：**105/105 全绿**

#### 下次 Session 起点

**路线 B+ S2.5**（循环图 v2 可视化大升级）：
- KB 扩展（substrate/product/couplings/schematic 字段，~6h）
- 合并细胞 renderer + 级联布局（~8h）
- 化学物耦合连线（~4h）
- 组选择下拉（~2h）
- SVG 导出（~3h）
- **预计 2-3 个 session 完成**（每次 3-5h）

完整计划：`C:\Users\REDLIZZ\.claude\plans\logical-drifting-metcalfe.md`

#### 今日关键学习

1. **工具的算法无偏 ≠ 叙述无偏**：再好的算法，如果 README 用"drives / 主导"这类
   因果词，就等于把上下文偏见嵌进了工具本身。S1 的修正是小改动大影响
2. **用户永远比工具更懂自己的研究**：当我提议"自动检测研究主题"时，
   用户一句"推不出砷修复"就戳破了伪命题。尊重这个认知层级
3. **产品定位比功能列表重要**：从"预装权威目录"到"可定制框架"的定位
   转向，看起来是削减功能，实际上是把工具的天花板打开了
4. **Mockup 驱动设计**：5 轮 mockup 迭代比 50 句自然语言描述更高效
5. **置换检验的价值**：不是"多一个严谨性卖点"，是**真正捕捉到共变伪信号**

---

### 2026-04-17（S2 — 循环图 v1.2 置换零假设 + 可信度标签）
- **阶段**：路线 B+ Session 2 完成。为每条 env-pathway 相关打可信度等级
- **核心改造**：
  - `inference.py` 新增 `_permutation_rho_p()`：999 次置换 + empirical p 替代 scipy Spearman p
  - 新增 `_confidence_label()`：综合 |ρ| + perm_p + sensitivity robust 打
    `strong / suggestive / weak / spurious? / unknown` 标签
  - `_env_correlations()` 接受 perm_n/perm_seed 参数
  - `infer()` 拿到 sensitivity 后**后处理**每条相关打 confidence 标签
  - `EnvCorrelation` dataclass 加 `perm_p` 和 `confidence` 字段
  - `renderer._draw_env_panel()` 显示彩色可信度徽章（绿=strong/橙=suggestive/红=spurious）
  - `CycleData.meta` 加 n_confidence_strong/suggestive/spurious 计数
  - `to_flat_stats()` 把 confidence 放在 top_mag 列，perm_p 放在 n_active_mags 列
- **真实数据结果**（高价值验证）：
  - **1 strong**：Ammonia oxidation ↔ Eh（ρ=0.835, perm_p<0.01, robust）
  - **10 suggestive**：Arsenate red / As transport 等（相关性可信但非最强）
  - **5 spurious?**：S1 讨论里怀疑的"ρ=0.764 共变"**被置换检验证实**
    —— 5 条通路共享同一 ρ=0.764 vs Total_As，不经 999 次置换考验
  - 这正是 **S1 去偏 + S2 置换**的价值：告诉用户哪些是真信号、哪些是共变伪信号
- **测试**：`tests/test_cycle_diagram.py` +4 case（perm 单元、confidence 逻辑边界、
  字段存在、分布），**105/105 全绿 41.7 s**（加速：fast mode perm_n=99）
- **论文可写**："confidence labels differentiate robust findings (n=1 strong) from
  threshold-dependent or confounding-inflated correlations (n=5 spurious) that
  would not survive permutation-based null hypothesis testing."
- **下一步**：S2.5 循环图升级（合并细胞 + 级联 + 化学物耦合，~20-25h）
- **量化**：
  - 修改：model +3 字段、inference +75 行（perm + confidence）、renderer +25 行（徽章）
  - 测试：101 → 105 case（+4）
  - 合计：~120 行净新增

### 2026-04-17（S1 — 循环图 v1 去偏：描述性中立语言 + 完整矩阵 + 敏感度扫描）
- **阶段**：路线 B Session 1 完成。响应 2026-04-17 讨论的"避免确认偏差"要求
- **核心改造**：
  - `paper/benchmarks/validation/cycle_diagram/README.md` 重写：
    - "drives / 主导" → "top-completeness contributor"
    - 新增"Interpretation Guidance"警告段：描述性 ≠ 因果性
    - 指出 5/5 As↔Total_As 高相关的替代解释（"高 As 选择有 arsC 的 MAG"）
    - 承诺 S2 置换检验 + S3 机制 YAML 评分器
  - `envmeta/geocycle/model.py`：`CycleData` 新增 `full_corr_matrix` + `sensitivity` 字段；
    新 dataclass `SensitivityRow`
  - `envmeta/geocycle/inference.py`：
    - `_env_correlations()` 返回 (filtered, full) 双列表
    - 新增 `_sensitivity_scan()`：扫描 `sensitivity_thresholds=[30, 50, 70]` 三档，
      每通路 Top-1 contributor 一致性 → `robust` / `threshold-sensitive` 标签
    - `DEFAULTS` 加 `sensitivity_thresholds`
  - `app.py` 循环图页：
    - 顶部增加"⚠️ 输出是描述性的，不是因果性的"警示条
    - 3 个展开区：主结果 / 完整相关矩阵 / 阈值敏感度
  - 验证产物拆分为 3 份 TSV（主 / full_corr / sensitivity）
- **实测数据**：
  - 18 通路里 **12 robust / 6 threshold-sensitive**（告诉用户哪些结论稳健）
  - **完整相关矩阵 68 条 vs 过滤后 16 条** —— 用户可看到"Arsenate reduction 还和 Eh 相关 0.42"
    这类潜在 confounding
- **测试**：`tests/test_cycle_diagram.py` +3 case（full matrix、sensitivity、新 stats types），
  **101/101 全绿 26.7 s**
- **下一步**：S2 置换零假设检验 + 可信度标签（strong/suggestive/weak/spurious?）
- **量化**：
  - 修改：model +22 行、inference +65 行、app.py +20 行、README 重写 65→95 行
  - 测试：98 → 101 case（+3）
  - 合计：~150 行净新增

### 2026-04-17（**Phase 3 v1 — 生物地球化学循环图自动推断 ⭐ 核心卖点**）
- **阶段**：Phase 3 v1 完成 — 论文 hero figure 自动生成
- **架构**：
  - `envmeta/geocycle/model.py`（~90 行）：dataclass（CycleData / ElementCycle /
    PathwayActivity / EnvCorrelation / MAGContribution），Phase 4 D3 JSON 导出
    的直接来源
  - `envmeta/geocycle/inference.py`（~220 行）：3 步推断引擎
    1. MAG 基础表（Phylum/Genus/is_keystone/abundance_mean）
    2. 逐通路活跃度（按 KB 18 通路，completeness ≥ 50% → 活跃 MAG；
       `completeness × log1p(abundance)` 排序贡献；`Σ completeness × abundance`
       = 通路总贡献）
    3. env-pathway Spearman 相关（样本级通路活性 = `Σ_MAG abundance × |KO∩pw|`）
  - `envmeta/geocycle/renderer.py`（~160 行）：matplotlib 静态渲染，2×2 元素象限
    + 底部 env 耦合面板
  - `envmeta/analysis/cycle_diagram.py`（~50 行）：对齐 analyze() 签名的薄壳
  - app.py "生物地球化学循环图" 页面（6 文件上传，5 个可选 + 4 参数 + 3 键下载）
- **已完成**：
  - 测试 `tests/test_cycle_diagram.py` 6 case，**98/98 全绿 24 s**
  - 论文积累 `paper/benchmarks/validation/cycle_diagram/` 含 PDF + stats + README
  - `time_comparison.md` 补循环图行（**无现成脚本**，首次实现）
- **推断结果复现用户研究假设**（重大里程碑）：
  - **Gallionella**（Fe 氧化专家）同时主导 NO reduction 和 Arsenate reduction
    → 铁氧化-氮-砷耦合 ✓
  - **Sulfuricaulis**（S 氧化专家）主导 Sulfide oxidation → 硫调控 Eh ✓
  - Ammonia oxidation ↔ Eh 强正相关（ρ=0.84, p=0.003）→ N 循环影响 Eh ✓
  - 5/5 砷代谢通路 ↔ Total_As 显著正相关 → As 压力-应答响应 ✓
  - **论文 Methods 可直接引用"算法自动推断支持原机制假设"**
- **下一步**：
  - 路线 A：回补 Phase 2 `mag_heatmap` + `network` (Gephi-prep)
  - 路线 B：进 Phase 4 — D3.js 交互编辑 + JSON 导出 + 环状布局
- **量化**：
  - 新增代码：model 90 + inference 220 + renderer 160 + cycle_diagram 50 +
    test 80 + app.py +100 = ~700 行
  - 测试：92 → 98 case（+6）
  - **累计分析图表：7 Reads + 3 MAG + 1 Cycle = 11 种**
  - **差异化卖点落地**：唯一支持从原始数据自动推断生物地球化学循环图的可视化工具

### 2026-04-17（Phase 2 迭代 3 — MAG 元素循环基因谱 / **Phase 3 前置 2/2 完成**）
- **阶段**：**Phase 3 循环图推断引擎的全部硬前置到齐**
- **已完成**：
  - 模块 B：`envmeta/analysis/gene_profile.py`（~220 行，vs 原脚本 1080 行）
    - MAG × KO 拷贝数矩阵（按元素 KB 顺序排列：As 11 / N 17 / S 15 / Fe 8 = 57 KO，6 个
      全 0 自动过滤 → 51 active）
    - 颜色 `log1p(copies)`，顶部元素色带 + 左侧门色带 + keystone ★
    - `max_mags` / `sort_by` / `element_filter` / `show_gene_names` 可调
    - 复用 pathway 模块的 `_parse_ko_annotation` / `PHYLUM_COLORS` / `_mag_col`
  - 模块 D：code_generator SUPPORTED 9 → 10
  - app.py：MAG 元素循环基因谱 页面（4 文件上传 + 参数 + 4 键下载）
  - 测试：`tests/test_gene_profile.py` 6 个 case，**92/92 全绿 20.3 s**
  - 论文积累：
    - `paper/benchmarks/validation/gene_profile/` 含 PDF + stats TSV + README
    - 168 MAG × 51 active KO，总拷贝 5144，Top-1 `Mx_All_158`（硫还原专家，Desulfobacterota_B）
    - 与 pathway 模块 Top-1 `Mx_All_63` 交叉验证一致
    - `time_comparison.md` 补基因谱行（1080 行 / 150 min vs 3 点击 / 3 s）
- **Phase 3 前置状态**：
  - ✅ pathway — 通路完整度矩阵
  - ✅ gene_profile — MAG × KO 拷贝数矩阵
  - ✅ gene_heatmap — 样本级 KO 丰度
  - ✅ KB — 18 通路 × 57 KO 定义
  - ✅ RDA + log2FC — 环境-功能关联
  - **全部到齐，可启动 Phase 3 geocycle/inference.py**
- **下一步**：**Phase 3 v1 — 循环图推断 + 静态渲染**（核心卖点）
  或按路线：回补 `mag_heatmap`（3h）+ `network`(Gephi-prep, 3h），再进 Phase 3
- **量化**：
  - 新增代码：gene_profile 220 + test 65 + app.py +85 + code_gen +12 = ~380 行
  - 测试：84 → 92 case（+8：6 新 + 2 参数化）
  - 累计分析图表：7 Reads + 3 MAG = **10 种**

### 2026-04-17（Phase 2 迭代 2 — MAG 通路完整度）
- **阶段**：Phase 2 第 2 图 — Phase 3 循环图硬前置 1/2 完成
- **已完成**：
  - 模块 B：`envmeta/analysis/pathway.py`（~320 行）
    - 对标 `scripts/python/08_pathway_completeness.py`（970 行），收敛到 ~320 行
    - 通路定义从**知识库** `elements.json` 动态读取（4 元素 × 18 通路 × 57 KO），
      不再硬编码；自动继承 KB 更新
    - 支持 heatmap / bubble 双样式；`max_mags` / `sort_by` / `element_filter` 参数
    - 元素彩条 + 门彩条 + keystone ★ 标记
  - 模块 E（知识库）：`envmeta/geocycle/knowledge_base/__init__.py` 新增
    `pathway_ko_sets()` / `pathway_element_map()` 两个索引工具
  - 模块 D：`code_generator` SUPPORTED 8 → 9
  - app.py：MAG-based 下 "代谢通路完整度" 页面（4 文件上传，3 个可选 + heatmap/bubble 切换）
  - 测试：`tests/test_pathway.py` 6 个 case，**84/84 全绿 16.7 s**
  - 论文积累：
    - `paper/benchmarks/validation/pathway/` 含 heatmap + bubble 双 PDF + stats TSV
    - 168 MAG × 18 通路完整度计算；14 keystone 正确标注
    - `time_comparison.md` 补通路完整度行（py 970 行 / 120 min vs EnvMeta 3 点击 / 3 s）
- **Phase 3 前置进度**：✅ pathway 完成，⏳ gene_profile 下一步
- **下一步**：`gene_profile`（Phase 3 前置 2/2）→ 然后 **Phase 3 v1 循环图**
- **量化**：
  - 新增代码：pathway 320 + kb helpers +18 + test 70 + app.py +85 + code_gen +12 = ~505 行
  - 测试：76 → 84 case（+8：6 新 + 2 参数化）
  - 累计分析图表：7 Reads + 2 MAG = 9 种

### 2026-04-17（Phase 2 迭代 1 — MAG 质量散点图）
- **阶段**：进入 Phase 2（MAG-based）。首个 MAG 图完成
- **已完成**：
  - 模块 B：`envmeta/analysis/mag_quality.py`（~260 行）
    - 对标 `scripts/python/06_MAG_quality.py`，简化输入：CheckM2 质量 + GTDB 分类
      + Keystone 列表（后两者可选）
    - Phylum 颜色 12 门 + Other；少于 `min_phylum_count` 合并为 Other
    - 质量 3 级（High/Medium/Low）阈值可调（`high_completeness/high_contamination/…`）
    - Keystone 菱形 + 标签高亮
    - stats 长表：summary（3 行）+ detail（每 MAG 质量分类）
  - 模块 D：`envmeta/export/code_generator.py` 新增 `_tpl_mag_quality`（SUPPORTED 7 → 8）
  - app.py：MAG-based 页面下 "MAG 质量评估"（3 文件上传，后两可选）+ 4 键下载
  - 测试：`tests/test_mag_quality.py` 5 个 case，**76/76 全绿 13.6 s**
  - 论文积累：
    - `paper/benchmarks/validation/mag_quality/` 含 PDF + stats TSV + README
    - 168 MAG：High 35 / Medium 111 / Low 22；14 个 keystone 标注
    - `time_comparison.md` 补 MAG 质量行（py 287 行 / 45 min vs EnvMeta 3 点击 / 2 s）
- **下一步**：Phase 2 剩余 4 图
  - `mag_heatmap` ← `07_MAG_abundance_heatmap.py`（Top30 MAG + 聚类 + 三段配色 ~3h）
  - `pathway` ← `08_pathway_completeness.py`（通路完整度气泡图 ~2h）
  - `gene_profile` ← `06_MAG_gene_profile.py`（MAG 元素循环基因谱 ~2h）
  - `network` ← `09_cooccurrence_network.py`（共现网络 ~3h）
- **量化**：
  - 新增代码：mag_quality 260 + test 55 + app.py +75 + code_gen +12 = ~402 行
  - 测试：69 → 76 case（+7：5 新 + 2 参数化）
  - 累计分析图表：7 Reads-based + 1 MAG-based = 8 种

### 2026-04-16（Phase 1 迭代 5 — LEfSe 收口 Phase 1 Reads-based）
- **阶段**：迭代 5 完成 — Phase 1 Reads-based 全部 7 张图闭环
- **前置修复**（基于 2026-04-15 R 对照）：`fix(rda)` commit e7a1b96 — 约束方差归一 +
  逐因子边际置换 F 检验 + `use_alias_labels` CK_1 风格标签。讨论 A/B/C/D 方案后
  选 A 保持现状（算法层对齐，工程/出版美化留 Phase 2）
- **已完成**：
  - 模块 B：`envmeta/analysis/lefse.py`（~180 行）
    - 内置 LEfSe 流程（不依赖 Galaxy 外部工具）：KW 筛选 + 组均值最大组 = enriched
      + `log10(1 + 1e6 × |max_mean - min_mean|)` 近似 LDA 效应量
    - 自动识别分类层级（`k__/p__/.../s__` 前缀正则）
    - 参数：`alpha_kw` / `lda_threshold` / `tax_levels` / `max_features` / `group_order` / `palette`
    - 水平条形图按组上色，组内按 LDA 降序
  - 模块 D：`envmeta/export/code_generator.py` 新增 `_tpl_lefse`（SUPPORTED 6 → 7 种）
  - app.py：新 LEfSe 子页面（丰度表 + metadata 两输入 + 4 键下载 + 参数侧边栏）
  - 测试：`tests/test_lefse.py` 6 个 case，**69/69 全绿 11.7 s**（含自动获得的 2 个 code_generator 参数化 case）
  - 论文积累：
    - `paper/benchmarks/validation/lefse/` 含 species + genus 两版 PDF + stats TSV + README
    - `time_comparison.md` 补 LEfSe 行（R 206 行 + Galaxy 外部 / 75 min vs EnvMeta 3 点击 / 3 s）
- **遇到的问题**：
  - 原 R 脚本 `04_LEfSe.R` 只绘制 Galaxy LEfSe 预跑结果 → EnvMeta 须内置算法，避免用户装 Galaxy
  - 样本量小（n=3-4/组），KW 显著门槛须放宽到 α=0.1 才有结果；文档注明
- **下一步**：**Phase 1 Reads-based 全部闭环 → 可进入 Phase 2（MAG-based 5 图）**
  - 优先级：`mag_quality` → `mag_heatmap` → `pathway` → `gene_profile` → `network`
  - 阻塞项：装 R + Galaxy LEfSe 做侧侧对比（论文 Methods 证据）
- **量化**：
  - 新增代码：lefse 180 + test_lefse 55 + app.py +70 + code_generator +12 = ~317 行
  - 测试：61 → 69 case（+8：6 新 + 2 参数化）
  - Reads-based 分析：6 → 7 种
  - **Phase 1 基本闭环**：7/7 Reads-based 图 + 代码生成器 + combined 样式 + 文件识别 7 类

### 2026-04-15（Phase 1 迭代 4 — RDA + 代码生成器 + combined 堆叠图）
- **阶段**：迭代 4 完成（代码生成器 + combined 堆叠图 + RDA 排序图）
- **已完成**：
  - 模块 D 深化：`envmeta/export/code_generator.py`（~175 行）
    - 支持 6 种分析模板（stackplot / pcoa / gene_heatmap / alpha_boxplot / log2fc / rda）
    - `generate(analysis_id, file_paths, params, output_base) -> str` 返回独立可运行的 Python 源码
    - 自动剔除 `_` 开头的私有键（中间 DataFrame）
  - 模块 B 扩展：
    - `envmeta/analysis/stackplot.py` 新 style="combined"（sample/group 并排 58%:42%），stats 用 MultiIndex 列
    - `envmeta/analysis/rda.py`（~220 行）：Hellinger + skbio RDA + Mantel 逐因子检验 + 样本 ID 自动对齐（两级命名兜底）
  - app.py：
    - 5 个分析页下方加「⬇️ 复现脚本（.py）」按钮（`_reproduce_button` 辅助函数）
    - 堆叠图 radio 加 combined 选项
    - 新 RDA 页面（3 文件上传 + 参数侧边栏 + 4 键下载）
  - 测试：`test_code_generator.py`（15）+ `test_stackplot.py`（+2）+ `test_rda.py`（4），**61/61 全绿 8.35 s**
  - 论文积累：`paper/benchmarks/validation/rda/` 含 PDF + stats + README；time_comparison.md 补 RDA 行（264 行 R / 60 min vs EnvMeta 3 点击 / 2 s）
  - Git：3 个 commit（code_generator / combined / RDA）
- **遇到的问题**：
  - env_factors.txt 的 SampleID 用别名（`CK_1`）而 Species.txt 列用原始 ID（`2_1`）→ RDA 模块内部用 (Group, 出现顺序) 兜底对齐
  - skbio 和 R vegan 的 RDA 轴方向可能镜像（数值绝对值一致）→ 不影响统计意义，留论文写作时统一
- **下一步**：Phase 1 迭代 5 = LEfSe（收口 Phase 1），之后 Phase 2 做 MAG-based 5 图
- **量化**：
  - 新增代码：code_generator 175 + rda 220 + stackplot 重构 +70 + app.py +140 + 测试 +80 = ~685 行
  - 测试：40 → 61 case（+21）
  - Reads-based 分析图表：5 → 6 种 + combined 样式
  - 代码生成器受益面：6 种分析自动支持下载 .py 复现

### 2026-04-14（Phase 1 迭代 3）
- **阶段**：迭代 3 完成（α 箱线图 + log2FC 差异柱图）
- **已完成**：
  - 模块 B：
    - `envmeta/analysis/alpha_boxplot.py`（~145 行）：多指数箱线图 + Kruskal-Wallis + 两两 Mann-Whitney U + BH 校正，子图网格 + 散点叠加 + 显著性横线
    - `envmeta/analysis/log2fc.py`（~175 行）：Welch's t-test + BH 校正 + log2FC（pseudocount 避零），4 元素 2×2 水平柱图 + 星号标注
  - app.py：
    - α 多样性页面（+75 行）：指数多选 / 子图列数 / p 阈值 / 画布尺寸 + 三键下载
    - log2FC 页面（+80 行）：组 A/B 双下拉 / 元素过滤 / padj 阈值 / |log2FC| 阈值 + 三键下载
  - 测试：`tests/test_alpha_boxplot.py`（4 case）+ `tests/test_log2fc.py`（5 case），**40/40 全绿 9.07 s**
  - 论文积累：
    - `paper/benchmarks/validation/alpha_boxplot/` 含 PDF + stats TSV + README
    - `paper/benchmarks/validation/log2fc/` 含 3 组对 PDF + stats TSV + README
    - `time_comparison.md` 补 α（R 260 行 / 35 min vs EnvMeta 3 点击 / 5 s）和 log2FC（py 230 行 / 50 min vs 4 点击 / 3 s）两行
  - Git：2 个 commit（α + log2FC）
- **遇到的问题**：
  - 首版 log2fc `element_filter` 只过滤绘图不过滤 stats → 测试失败 → 改为同时过滤（TSV 导出也一致）
  - 小样本（n=3-4/组）导致所有 padj > 0.05 —— 真实研究数据（n≥10）应能观察到显著差异
- **下一步**：Phase 1 迭代 4 — 代码生成器雏形 + RDA + combined 样式 + LEfSe（见顶部计划）
- **量化**：
  - 新增代码：alpha_boxplot 145 + log2fc 175 + 测试 120 + app.py +155 = ~595 行
  - 测试：31 → 40 case（+9）
  - 累计分析图表：3 → 5 种（堆叠图/PCoA/元素循环热图/α 多样性/log2FC）

### 2026-04-13（晚 — Phase 1 迭代 2）
- **阶段**：迭代 2 完成（PCoA + 元素循环热图 + 文件识别扩展 + 模块 C 雏形）
- **已完成**：
  - 模块 A：detector 扩展 5 类文件 — DISTANCE_MATRIX、ALPHA_DIVERSITY、CHECKM_QUALITY、ENV_FACTORS、KO_ABUNDANCE_WIDE（规则按优先级排序 + 首个匹配即采纳）
  - 模块 B：
    - `envmeta/analysis/pcoa.py`（skbio PCoA + PERMANOVA + adjustText 标签防重叠，含组间两两对比）
    - `envmeta/analysis/gene_heatmap.py`（4 元素合并单图版，KO × Group，左侧通路/元素色块，含元素过滤）
    - `envmeta/analysis/stackplot.py` 追加 `sort_by` / `reverse_stack` 参数（用户 2026-04-13 反馈）
  - 模块 C 轻量：`envmeta/params/common.py`（render_figure_size / render_font_controls / render_dpi_selector）
  - 知识库加载器：`envmeta/geocycle/knowledge_base/__init__.py`（lru_cache + flat_ko_map / element_colors / pathway_display）
  - 测试数据：`scripts/python/gen_ko_sample_data.py` 过滤 eggnog 24592 行 → 51 个知识库 KO，生成 `tests/sample_data/ko_tpm.spf` (5.8 KB)
  - app.py：3 个分析页面（堆叠图/PCoA/热图）端到端可用，参数面板复用 common.py 组件
  - 测试：新增 test_pcoa (4)、test_gene_heatmap (4)、test_kb_loader (4)、扩展 test_detector (+4)、test_stackplot (+2)，**累计 31/31 全绿 5.5 s**
  - 论文积累：
    - `paper/benchmarks/validation/{pcoa,gene_heatmap}/` 各含 PDF + stats TSV + README
    - `paper/benchmarks/time_comparison.md` 补 PCoA 和基因热图两行
- **遇到的问题**：
  - 规则冲突：env_factors（SampleID+Group+多数值列）和 metadata（SampleID+Group）同时匹配，metadata conf 更高反而吃掉 env_factors → 改为按优先级顺序首个匹配即返回
  - `skbio.DistanceMatrix` 要求 C-contiguous 数组，`np.ix_` 切片结果不连续 → `np.ascontiguousarray`
- **遗留**：
  - log2FC 差异柱图（Fig2-9）、combined 样式堆叠图、α 多样性箱线图、RDA — 留到迭代 3
  - 装 R 后做 EnvMeta vs R 侧侧 PDF 对比
- **下一步**：Phase 1 迭代 3（剩余分析图 + 代码生成器）或进入 Phase 2
- **量化**：
  - 新增代码：detector +98 / stackplot +30 / pcoa 207 / gene_heatmap 260 / common 46 / kb_loader 61 / app.py +140 / 测试 +240 = ~1080 行
  - 依赖 +1：adjustText
  - 测试：12 → 31 个 case（+19）
  - 文件识别：2 → 7 种类型（覆盖 tests/sample_data/ 里 10/12 种常见文件）

### 2026-04-13（下午 — Phase 1 迭代 1）
- **阶段**：Phase 1 迭代 1 完成（端到端打通堆叠图）
- **已完成**：
  - 模块 A 雏形：`envmeta/file_manager/detector.py`（识别 metadata + abundance_wide，含编码/分隔符嗅探，规则可扩展）
  - 模块 B 接口：`envmeta/analysis/base.py`（AnalysisResult dataclass）
  - 第一图：`envmeta/analysis/stackplot.py`（sample / group 双样式，Top-N + Others，matplotlib 渲染）
  - 模块 D 雏形：`envmeta/export/figure_export.py`（PNG / PDF，含 export_to_bytes 供 Streamlit 用）
  - app.py 改造：文件管理页面（多文件上传 + 自动识别 + 预览 + 手动修正）+ 堆叠图页面（参数侧边栏 + 实时生成 + 三键下载 PNG/PDF/TSV）
  - 测试：`tests/test_detector.py` (8) + `tests/test_stackplot.py` (4)，**12/12 全绿**
  - 论文积累：`paper/benchmarks/validation/stackplot/` 含 sample/group 两版 PDF + 百分比 TSV + README，`time_comparison.md` 填入堆叠图行
- **遗留 TODO**：
  - VS Code Python 解释器需手动选 envmeta 环境（Ctrl+Shift+P → Python: Select Interpreter）
  - 待装 R 后做 EnvMeta vs R 的侧侧 PDF 对比
  - 截图 `paper/figures/screenshot_stackplot.png` 待手动截
- **下一步**：Phase 1 迭代 2 — 加 PCoA + 热图 + 模块 C（参数面板结构化）
- **量化**：
  - 新增代码：detector 149 + base 14 + stackplot 175 + figure_export 44 + 测试 80 + app.py +180 ≈ 640 行
  - 测试覆盖：12 个 case，1.32 秒跑完
  - 效率对比（堆叠图）：写 R 脚本 177 行 / ~30 min，vs EnvMeta 3 次点击 / ~10 s

### 2026-04-13
- **阶段**：Phase 0 全部完成 → 准备进入 Phase 1
- **已完成**：
  - 项目骨架：envmeta 包结构（5 子模块）+ tests/ + paper/
  - 开发环境：Miniconda + conda envmeta（Python 3.11）+ 13 个依赖全部验证通过
  - app.py Streamlit 框架（侧边栏导航 + 5 大模块占位页面，浏览器实测通过）
  - 测试数据：15 个论文真实数据精简文件放入 `tests/sample_data/`，带 README 索引
  - README.md 改写为 EnvMeta 工具说明（替换旧论文项目内容）
  - INSTALL.md 环境安装指南
  - 元素循环知识库 v1：`envmeta/geocycle/knowledge_base/elements.json`
    4 元素 × 18 通路 × 57 KO（As 17 / N 17 / S 15 / Fe 8），含元素间耦合关系
  - 论文积累框架 paper/：benchmarks/validation、time_comparison（含"一键复现"列）、
    performance、figures、user_study（SUS 量表 + task_protocol）、manuscript、tool_comparison
  - requirements.txt 干净环境复现验证通过（临时 envmeta_verify 环境）
  - GitHub 私有仓库 `redlizzxy/EnvMeta` 创建并完成首次推送
- **Git 历史**（3 个 commit）：
  - `6cab58f` feat: Phase 0 项目初始化
  - `d1372a4` docs: 初始化论文数据积累框架
  - `7778cc4` feat: Phase 0 收尾（测试数据 + README + 知识库）
- **下一步**：Phase 1 — 模块 A（文件识别引擎）+ 模块 B（堆叠图/PCoA/热图）+ 模块 C（基础调参）+ 模块 D（基础导出）
- **遇到的问题**：
  - Windows Store Python 占位符（exit code 49）→ 安装 Miniconda 解决
  - conda 未加入 PATH → `conda init powershell` 解决
  - conda 首次使用需接受 Terms of Service → `conda tos accept` 解决
  - PowerShell 转义规则与 bash 不同 → 验证命令在 bash 跑、不在 PowerShell 里跑
  - gh CLI 在 bash session 读不到 keyring → 推送改在 PowerShell 里执行
  - `gh repo create ... --push` 时 git 未配置 credential helper → `gh auth setup-git` 后再 push
- **量化**：
  - 代码/文档：60 + 18 = 78 个文件入库，累计 ~15,000 行（含原论文脚本参考）
  - 知识库：51 → 57 KO（源码实际值，超出 CLAUDE.md 原声明 6 个）
  - 环境复现：从零装 conda → app 跑通全程 ~15 分钟
  - 测试数据：586 KB，覆盖 12 种图表全部输入类型
