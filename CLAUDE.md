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

## 论文发表数据积累

本工具计划作为方法学论文发表（目标期刊：iMeta / Bioinformatics / Frontiers in Microbiology）。
开发过程中需要持续积累以下四类数据：

### 积累规则

1. **功能验证**：每完成一个分析模块，用 `tests/sample_data/` 中的论文真实数据跑一遍，将输出保存到 `paper/benchmarks/validation/`，与原始 R/Python 脚本的输出对比，确认结果一致。
2. **效率对比**：每完成一个模块，记录"传统方式（写代码）需要多少行代码/多少步骤/多少时间" vs "EnvMeta 需要多少次点击/多少时间"，写入 `paper/benchmarks/time_comparison.md`。
3. **截图留档**：每个模块完成后截一张 EnvMeta 的界面截图，保存到 `paper/figures/screenshot_模块名.png`。
4. **开发日志量化**：每次开发日志除了记录做了什么，加一行量化数据（代码行数、验证结果、耗时对比）。

## 开发日志

> 每次 session 结束前更新此区块。新对话开始时 Claude Code 自动读取，了解当前进度。

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
