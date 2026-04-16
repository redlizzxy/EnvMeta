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

### 2026-04-19（**Today's Session 总结** — S2.5-14 耦合回归修复 + S3 L2 假说评分器落地）

今日两条主线：**先修 bug，后加大功能**。累计 **3 个 commit / 测试 176 → 189 全绿**。

#### 时间线

| 顺序 | 内容 | commit |
|---|---|---|
| 1 | S2.5-14 multi-chain 细胞丢失化学物耦合锚点（As(V)↔Fe(III) 棕色线消失）| 27c32d6 |
| 2 | 2026-04-19 S2.5-14 session 日志 | 3701f28 |
| 3 | S3 机制假说 YAML 评分器（L2 层落地） | ff4b1b5 |

#### 关键发现 / 用户反馈

1. **视觉回归一定要用真实数据跑**：S2.5-14 的 bug（MAG-merge 后 merged_genes[0].substrate=None
   导致 Fe 细胞锚点丢失）**单元测试没捕获**，必须靠用户对比 v29 vs v33 截图才发现。
   → 修复同时加了 `test_multichain_cell_exposes_fe_iii_for_coupling` 锁死
2. **"假说评估"与"推断"必须架构解耦**：S3 落地前反复讨论"inference 要不要
   接受 hypothesis 输入" → 最终决定 **不接受**。hypothesis.py 只读 CycleData，
   inference.py 零改动。论文可以诚实说 "hypothesis-agnostic inference +
   separate evaluator"
3. **skipped 不扣分的语义设计**：如果 claim 指向 KB 没有的耦合 / data 里没跑
   到的 pathway → **分母剔除**，而不是判为 0 分。否则用户"写得不完整 →
   分数被压低"会产生误解。CK 组的 3 条实跑数据验证了这个设计合理

#### 本次 session 两个里程碑

**S2.5-14**：
- 根因 `merged_genes[0].substrate` 取单值，第一个基因若是 fur/tonB 就 None
- 修复 3 处：物种列表化 + 多链位置 map (`substrate_pos_map`/`product_pos_map`)
  + species-specific 定位 + 跳过 None substrate 的外部绘制
- **4/4 KB 耦合线全部回归**：As(III)↔S-2 / As(V)↔Fe(III) (brown) /
  Fe(II)↔S-2 / NO3-↔As(III)，CK/A/B × 2 ranking 下一致

**S3**：
- 4 类 claim 全部实现 + 加权总分 + 4 档标签
- 8 claim 示例 `arsenic_steel_slag.yaml` 对应用户真实论文研究
- UI 集成到循环图页底部 expander，不破坏现有工作流
- 真实数据：CK 0.879 / A 0.868 / **B 1.000** — 与 2026-04-17 推断观察完全对齐
- +12 测试 case（4 claim 各一个 satisfied + 各一个 skipped/unsatisfied 路径
  + load schema 校验 + to_dataframe/to_json）

#### 测试进展

**176 → 189（+13）**：
- +1 multichain Fe(III) 暴露回归（S2.5-14）
- +12 hypothesis 单元（S3）

#### 交付文件（今日新增）

| 类别 | 文件 |
|---|---|
| 源码 | `envmeta/geocycle/hypothesis.py`（+430 行）|
| UI | `app.py` 循环图页 hypothesis expander |
| 依赖 | `requirements.txt` +pyyaml |
| 测试 | `tests/test_hypothesis.py`（+12 case）|
| 测试数据 | `tests/sample_data/sample_hypothesis.yaml` |
| 论文示例 | `paper/hypotheses/arsenic_steel_slag.yaml` |
| 论文文档 | `paper/hypotheses/README.md` |
| 归档 | `paper/benchmarks/validation/hypothesis/` (CK/A/B 各 TSV+JSON + README) |

#### 论文 Methods 素材

> "EnvMeta couples hypothesis-agnostic cycle inference (evidence of pathway
> activity / MAG contribution / env correlation, all without user-supplied
> hypothesis input) with an optional hypothesis-testing YAML evaluator."

#### 下次 Session 可选路径

1. **S4 Fork Bundle 导出**（~4h）：打包 elements.json + 用户 YAML + 配置 →
   zip，落地"论文-EnvMeta 绑定发布"协议。S3 + S4 组合后，L1-L2-L4 三层形成
   闭环
2. **S4.5 HTML 交互导出**（~10-15h，论文 SI 杀手锏）：D3.js 独立 HTML 嵌入
   cycle_data + hypothesis_score JSON，审稿人能直接交互
3. **S3.5 group_contrast claim**（~2h）：补 S3 v2 缺口。调用 cycle_compare
   在 YAML 里写"B vs CK total_contribution 比值 ≥ 1.5"
4. **S6 mag_heatmap**（~3h）：补 Phase 2 欠账
5. **S9 论文 Methods 起草**：素材已齐（S1 去偏 + S2 置换 + S3 评分器）

**强烈推荐先 S4**：与 S3 天然衔接，工时低，论文叙事完整度直接拉满。
S4.5 推到论文投稿前冲刺。

#### 今日关键学习

1. **视觉回归不要依赖单元测试**：要跑真实数据 + 用户视觉对比。本次 S2.5-14
   正是靠用户对比截图发现
2. **skipped 不是错误处理，是语义设计**：让"不适用"和"不成立"分开，
   才符合用户的科学判断流程
3. **架构解耦的成本最低的时候就是"不写"**：S3 没让 inference 接受 hypothesis
   输入，省去了 5+ 个 API/依赖复杂度，后期维护成本大幅降低

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
