# EnvMeta — 环境微生物宏基因组可视化分析平台

> **本文件的角色**：项目核心规范、架构决策、未实现功能路线。新 session 开始时
> Claude Code 自动读取。**不记历史日志** —— 那些在 [DEBUG_NOTES.md](DEBUG_NOTES.md)。

## 项目概述

EnvMeta 是面向环境微生物研究的宏基因组下游可视化分析工具。核心解决"从数据文件
到发表级图形"的最后一公里问题，独家能力是自动生成微生物-环境互作的生物地球化学
循环图。

**目标用户**：环境微生物方向的硕博研究生和科研人员（非生信背景也能用）。

**核心差异化**（vs Shiny-phyloseq / plotmicrobiome / Anvi'o / 测序公司云平台）：
GUI 交互调参 + 代码同步生成 + 元素循环图自动推断 + 假说 YAML 评分器 +
独立交互 HTML 导出 —— 五件套同时具备是业界空白。

**在线版**：<https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/>
**仓库**：<https://github.com/redlizzxy/EnvMeta>

## 技术栈

- **前端 UI**：Streamlit
- **出版级渲染**：matplotlib + seaborn
- **统计**：pandas / scipy / scikit-bio / statsmodels / networkx
- **交互 HTML**：D3.js v7（inline 嵌入，无外部依赖）
- **循环图推断**：自研（`envmeta/geocycle/`）— KEGG-driven KB + 置换检验 + 敏感度扫描
- **Python**：3.11+（所有开发按 3.11）

## 项目结构

```
envmeta/
├── app.py                         # Streamlit 入口
├── CLAUDE.md / DEBUG_NOTES.md     # 项目规范 + 开发日志（拆分）
├── envmeta/                       # 核心包
│   ├── file_manager/              # 模块 A — 文件识别（11 FileType）
│   ├── analysis/                  # 模块 B — 14 图分析引擎（7 Reads + 5 MAG + 1 Cycle + 1 hypothesis）
│   ├── params/                    # 模块 C — 参数管理 + 共享 widget
│   ├── export/                    # 模块 D — 导出（PNG/PDF/SVG/TIFF + code_generator）
│   ├── geocycle/                  # 模块 E — 循环图
│   │   ├── knowledge_base/        #   KB JSON（4 元素 × 18 通路 × 57 KO）
│   │   ├── inference.py           #   推断引擎（去偏 + 置换 + 敏感度）
│   │   ├── model.py               #   dataclass（CycleData / HypothesisScore / ...）
│   │   ├── hypothesis.py          #   YAML 假说评分器（5 类 claim）
│   │   ├── cell_renderer.py       #   matplotlib 细胞渲染（Mockup 10 布局）
│   │   ├── renderer.py            #   matplotlib 主渲染
│   │   ├── html_exporter.py       #   独立交互 HTML 导出
│   │   └── templates/             #   cycle_interactive.html + d3.v7.min.js
│   ├── help/                      # 新手落地包（向导 / 解读 / 反向索引）
│   └── tools/                     # Fork Bundle / KB builder / Gephi prep / hypothesis validator
├── docs/
│   ├── data_preparation_zh.md     # 上游工具 → EnvMeta 输入映射
│   └── install_for_beginners.md   # 小白安装指南
├── paper/
│   ├── bundles/                   # 论文 Fork Bundle 示例
│   ├── benchmarks/                # 验证数据 + 效率对比
│   ├── hypotheses/                # 假说 YAML 示例
│   ├── figures/                   # 截图 + mockup
│   └── user_study/                # 评测问卷设计 + 海报 + 部署指南
├── tests/
│   ├── sample_data/               # 论文精简数据（首页一键加载）
│   └── test_*.py                  # 291 case 全绿
├── scripts/                       # 原论文脚本（参考实现，不改）
└── requirements.txt
```

## 元素循环知识库

`envmeta/geocycle/knowledge_base/elements.json`（v2.0，KEGG-driven）：
**4 元素 × 18 通路 × 57 KO**

- **砷代谢（11 KO）**：arsC-grx/trx/myc、aoxA/B、arsB、ACR3、arsA、arsH、AS3MT、arsR
- **氮代谢（17 KO）**：narG/H/I、napA/B、narB、nirK/S/A、CYP55、norB、nosZ、amoA/B/C、nifD/K
- **硫代谢（15 KO）**：cysJ/I/H、sir、aprA/B、dsrA/B、sqr、soxA/X/B/Y/Z、TST
- **铁代谢（8 KO）**：ABC.FEV.A/P/S、TC.FEV.OM、afuA、FTR、fur、tonB

KB 可通过 `envmeta kb-build` CLI 从 KEGG snapshot 重建，或手工编辑 JSON。

## 研究背景（用于循环图推断的参考数据集）

本工具源自砷渣-钢渣微生物修复宏基因组研究：
- 三组实验：CK（对照）、A（低钢渣配比）、B（高钢渣配比）
- 标准配色：CK=#4DAF4A（绿）、A=#377EB8（蓝）、B=#E41A1C（红）
- 核心机制（作者假说）："铁氧化固砷 → 硫/氮循环调控 Eh → 局部硫化沉砷"
- 168 MAG × 10 样本 × 57 目标 KO
- 关键物种：Sulfuricaulis（硫氧化）、Gallionella（铁氧化）、Thiobacillus（硫氧化）

## 五大功能模块

### 模块 A：智能文件管理器
上传文件 → 表头规则匹配 → 自动识别 11 种类型 → 反向索引提示可跑分析。
**abundance 分 MAG 级 (conf=0.95) / TAXON 级 (conf=0.88) 两档**。

### 模块 B：可视化分析引擎
14 图：
- **Reads-based**（7）：堆叠图 / α / β-PCoA / RDA / LEfSe / 基因热图 / log2FC
- **MAG-based**（5）：MAG 质量 / 丰度热图 / 通路完整度 / 基因谱 / 共现网络 Gephi-prep
- **循环图**（1）：4 元素循环 + 跨元素耦合
- **假说评分**（1）：YAML evaluator

### 模块 C：交互式图形调参
Streamlit 控件实时调参（颜色 / 字体 / 大小 / 图例）→ 参数字典 → 刷新预览。
MAG-based 4 图共享统一 4 层参数面板（`envmeta/analysis/_mag_common.py`）。

### 模块 D：高质量导出
图形 PNG / PDF / SVG / TIFF 600dpi + 可运行 `.py` 复现脚本 + 统计 TSV。
**导出中心**统一入口（图表 / Bundle / 脚本 / 文档 4 tab + 批量 ZIP）。

### 模块 E：生物地球化学循环图生成器 ⭐
核心创新。Mockup 10 合并细胞 + 级联布局 + 跨元素化学物耦合 + 假说评分 +
独立交互 HTML 导出（400-550 KB 离线可用）。**EnvMeta 差异化独家能力。**

## 产品定位与核心设计决策

### 定位

EnvMeta **不是**预装所有环境机制的权威工具，而是：

> **可定制的环境生信分析框架 + 通用循环图推断引擎 + 论文-EnvMeta 绑定发布协议**

### 核心设计原则

| 原则 | 含义 | 理由 |
|---|---|---|
| **领域中立** | 核心代码不内嵌具体机制/研究主题 | 科学问题空间无穷；避免确认偏差 |
| **用户自带知识** | 机制 YAML / KB 扩展 / 插件都由用户提供 | 不承担社区知识更新负担 |
| **完全离线** | 所有功能本地可用，无网络依赖 | 数据保密 + 降低维护复杂度 |
| **Fork 而非社区** | 每篇论文绑定一份定制 EnvMeta 发布 | 分布式发布，无中心化维护瓶颈 |
| **描述而非断言** | 输出"谁承载什么 + 什么与什么相关"，不下因果结论 | 因果判断是用户职责；工具避免附和假说 |

### 拒绝的设计

- ❌ 自动检测研究主题（env 列名 + KO 频度**推不出**"砷修复"等研究意图）
- ❌ 下拉选研究领域（无法穷举科学问题）
- ❌ 预装竞争机制目录（违反"不维护社区"原则）
- ❌ 嵌入 KEGG 快照（过期 + 违反离线原则）
- ❌ 完整 GUI KB 编辑器（5% 用户会用，投入回报低）

### 5 层架构

| 层 | 内容 | 状态 |
|---|---|---|
| L1 通用循环图内核 | 描述性推断 + 去偏（S1）+ 置换（S2）+ 敏感度 | ✅ 完成 |
| L2 机制 YAML 评分器 | 用户上传假说 YAML → 证据评分 + null_p | ✅ 完成（S3/S3.5）|
| L3 插件框架 | 用户上传 Python `analyze()` → 自动注册 GUI | ⏸ 论文接收后 |
| L4 Fork Bundle | 打包 KB + YAML + config + KEGG 快照为 zip | ✅ 完成（S4）|
| L5 KB 工具 | schema 校验 + template + diff/merge CLI | ✅ 完成（kb-build CLI）|

### 市场现实（诚实校准）

- 预期峰值：**50-200 活跃用户 / 1-2 年**（不追求上千）
- 真实采用：自己课题组（确定）+ 同行复现（每篇 7-35 人）+ 小众外部（5-50 实验室）
- 论文目标：**iMeta / Bioinformatics / Front Microbiol 方法学论文**；博士论文章节完美匹配

### 真实竞争对手：测序公司云平台（非 Anvi'o）

环境微生物课题组真实痛点不是"找不到工具"，而是"**测序公司给了数据不知道怎么
做进阶分析**"。EnvMeta 定位 = "**比云平台更灵活 + 比 Anvi'o 更易用**" 的中间层。

## 开发规范

- **语言**：Python 3.11+，代码注释用中文
- **Git commit 前缀**：`feat`（新功能）/ `fix`（修 bug）/ `refactor` / `docs` / `test` / `chore`
- **函数**：写 docstring，参数类型标注
- **分支**：`master`（稳定版）→ `feature/xxx`（功能分支，推送前合并）
- **测试数据**：`tests/sample_data/`，用论文真实数据的精简版
- **.gitignore**：`data/` 本地数据不入库；`.claude/` 本地 Claude 状态不入库

### 日常开发循环

```powershell
conda activate envmeta
cd D:\workdata\envmeta
streamlit run app.py            # http://localhost:8501
pytest tests/ -q                # 跑测试
```

### 提交代码

```powershell
git add <改动的文件>            # 避免 git add -A
git commit -m "feat: xxx"
git push origin master          # Streamlit Cloud 自动 2-3 min 内重新部署
```

## 论文发表数据积累规则

本工具计划作为方法学论文发表（目标：iMeta / Bioinformatics / Frontiers in
Microbiology）。开发时持续积累：

1. **功能验证**：每完成一个分析模块，用 `tests/sample_data/` 跑一遍 → 输出存
   `paper/benchmarks/validation/`，与原 R/Python 脚本对比
2. **效率对比**：记录"传统代码行数/时间" vs "EnvMeta 点击/时间" →
   `paper/benchmarks/time_comparison.md`
3. **截图留档**：每模块截 `paper/figures/screenshot_模块名.png`
4. **开发日志量化**：每次日志加一行量化数据（代码行数、验证结果、耗时对比），
   写入 [DEBUG_NOTES.md](DEBUG_NOTES.md)

## 当前进度（2026-04-21，v0.8.1）

**v0.8.1 内测版（在线可用）**：
- Phase 0-3 全部完成
- 测试 **293/293 全绿**
- **12 种分析图表 + 独立交互 HTML**
- 更新日志见 [CHANGELOG.md](CHANGELOG.md)

### 已完成（按架构层分组）

- ✅ **Phase 0** 项目骨架 + 环境 + 知识库
- ✅ **Phase 1** 7 张 Reads-based 图（堆叠图/α/PCoA/RDA/LEfSe/基因热图/log2FC）+ 代码生成器
- ✅ **Phase 2** 5 张 MAG-based 图（MAG 质量/丰度热图/通路完整度/基因谱/网络 Gephi-prep）
- ✅ **Phase 3 核心** 循环图 v2（Mockup 10 合并细胞）+ 去偏/置换/敏感度 + YAML 评分器 + Fork Bundle + KEGG-driven KB
- ✅ **S8-ux 新手落地包** 数据准备指南 + 图表选择向导 + 14 图解读 expander + 样例数据一键加载
- ✅ **T1 导出中心** 4-tab 统一入口 + 批量 ZIP
- ✅ **T2 HTML 交互导出** 400-550 KB 独立 HTML（D3 嵌入）+ per-group 切换 + 化学物精确锚点 + 上下标
- ✅ **在线部署** Streamlit Cloud + 本地 Windows/Linux 行为一致性修复
- ✅ **内测素材** 腾讯问卷 + 海报 + 部署指南 + 小白安装
- ✅ **v0.8.1 Mac 端修复** conda ToS / pip protobuf resolver / macOS YAML 评分 ENAMETOOLONG / HTML 切组后 chem-link 失效

## Backlog（投稿前 + Phase 4）

### 🟥 iMeta/Bioinformatics 投稿硬指标

| 任务 | 工时 | 理由 |
|---|---|---|
| **English README + LICENSE + Zenodo DOI** | ~4-6h | iMeta 投稿硬性要求；LICENSE 先选 MIT |
| **R 侧侧对比验证** | 1-2 天 | 论文 Methods 关键证据，需装 R 环境跑 6 图 |
| **第二数据集复现** | 2-3 天 | 审稿人必问；候选：Oak Ridge 铀污染 / 秘鲁砷湖 / 阿拉斯加冻土 |
| **论文 Methods + Results + Discussion 起草** | 1-2 周 | 素材已齐（S1 去偏 + S2 置换 + S3/S3.5 评分器 + S4 Bundle + 全 12 图 + S8-ux + T1 + T2 HTML） |
| **User study 数据回收分析** | 1 周 | 问卷已发（2026-04-19），1 周内回收，汇总 `paper/user_study/results.md` |

### 🟡 加分项（推荐做）

| 任务 | 工时 | 理由 |
|---|---|---|
| Q7 新手 UX 精调 | ~1h | HTML 控件分组 / SVG-JSON 下拉 / 3 步 onboarding / tab 重排 / 参数 debug mode |
| 大数据集性能 benchmark | 1 天 | 500+ MAG × 50+ sample，证明可扩展性 |
| tool_comparison.md 细化 | 4h | 表格对比 Krona / Anvi'o / MicrobiomeAnalyst / 测序公司云平台 |

### 🟢 论文接收后再做（Phase 4）

- **L3 插件框架**：用户上传 Python `analyze()` → 自动注册 GUI（CLAUDE.md 明确承诺延后）
- **KB version 字段升 v2.0**：`elements.json` 里 `"version": "1.1"` 与实际含 v2.0 字段不匹配
- **S3.5-doc DOI 验证**：YAML 里 6 个理论文献 DOI 待用户手动过一遍
- **循环图 D3 编辑器**：用户直接拖拽编辑节点/通路（Phase 4 主线）

### 已知小 bug（待收口）

- "文件管理"首页同名覆盖策略（当前只处理新 name）
- matplotlib 中文字体在部分 Windows 环境显示方块（需要打包字体）

## 下次 session 建议起点

按优先级：
1. **English README + LICENSE** — 最小工时收益最大（4-6h → 扫清投稿硬指标 2/4）
2. **R 对照** — 要装 R 环境，时间成本高但论文必需
3. **论文 Methods 起草** — 素材已齐，可先写 Methods 3000 字 demo
4. **User study 回收** — 等 1 周问卷，并行做前三项

---

## 历史日志 / 自检报告

详见 [DEBUG_NOTES.md](DEBUG_NOTES.md) — 2026-04-13 至 2026-04-18 全部 session 日志、
阶段性自检（v1/v2）、产品定位讨论、技术决策记录、用户反馈处理。
