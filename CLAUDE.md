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

## 下次 session 计划（2026-04-16 起）

**当前进度**：Phase 0 + Phase 1 迭代 1/2/3/4 完成。Phase 1 Reads-based 全部 6 种分析（堆叠图 / PCoA / 基因热图 / α 多样性 / log2FC / RDA）+ 代码生成器 + combined 样式端到端可用。文件识别覆盖 7 类，测试 61/61 全绿。

**Phase 1 基本收口** — 剩余只有 LEfSe（迭代 5 或并入 Phase 2）。

**建议下次做**（按优先级）：

1. **LEfSe 差异分析**（移植 `scripts/R/04_LEfSe.R`，~2 小时）
   - LDA + KW 分析，输出柱图（LDA score 由大到小）
   - Python 端用 `scipy.stats.kruskal` 做总检验 + `sklearn.discriminant_analysis.LinearDiscriminantAnalysis` 算 LDA score
   - 新建 `envmeta/analysis/lefse.py` + app.py 页面 + 测试
   - 完成后 Phase 1 Reads-based 全部 7 图闭环

2. **进入 Phase 2：MAG-based 图表**（按优先级排序）
   - `envmeta/analysis/mag_quality.py` ← `scripts/python/06_MAG_quality.py`（MAG 质量散点图，最简单 ~1.5h）
   - `envmeta/analysis/mag_heatmap.py` ← `scripts/python/07_MAG_abundance_heatmap.py`（Top30 MAG 丰度热图，含聚类 ~3h）
   - `envmeta/analysis/pathway.py` ← `scripts/python/08_pathway_completeness.py`（通路完整度气泡图 ~2h）
   - `envmeta/analysis/gene_profile.py` ← `scripts/python/06_MAG_gene_profile.py`（MAG 元素循环基因谱 ~2h）
   - `envmeta/analysis/network.py` ← `scripts/python/09_cooccurrence_network.py`（共现网络 ~3h）

3. **SVG / TIFF 导出**（~0.5 小时）
   - 修改 `envmeta/export/figure_export.py` 增加 SVG + TIFF 格式支持
   - 每页导出区追加 2 个下载按钮（对应"SCI 投稿"和"毕业论文"预设）

**建议顺序**：先 1（小但收尾 Phase 1），再按 MAG 图表的依赖与复杂度递增开 Phase 2。

**阻塞项**：装 R 做 EnvMeta vs 原脚本的 PDF 侧侧对比（论文关键验证数据，目前所有 6 张 Reads-based 都等这一步）。

**里程碑**：
- Phase 1 基本完成 → 开始准备论文 Methods 草稿（代码生成器 + 效率对比表 + validation PDF 都已齐全）
- Phase 2：补齐 5 张 MAG 图 → v0.5 内部测试版
- Phase 3：循环图生成器（核心创新 → v1.0 发布）

## 开发日志

> 每次 session 结束前更新此区块。新对话开始时 Claude Code 自动读取，了解当前进度。

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
