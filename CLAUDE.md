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
│   ├── install_for_beginners.md   # 小白安装指南
│   └── legacy/                    # 早期 EnvMeta 文档（INSTALL / 项目细化方案 v1）
├── paper/                         # EnvMeta 方法学论文素材（开源）
│   ├── manuscript/                # 论文稿件
│   ├── benchmarks/                # 验证数据 + 效率对比
│   ├── bundles/                   # 论文 Fork Bundle 示例
│   ├── hypotheses/                # 假说 YAML 示例
│   ├── figures/                   # EnvMeta app 截图 + mockup
│   ├── tool_comparison.md
│   └── user_study/                # 评测问卷设计 + 海报 + 部署指南
├── tests/
│   ├── sample_data/               # 论文精简数据（首页一键加载）
│   └── test_*.py                  # 293 case 全绿
├── scripts/                       # EnvMeta KB CLI 工具
│   ├── build_kegg_snapshot.py     # KB 重建 CLI（envmeta kb-build）
│   └── seed_ko_list.json          # KB 种子 KO 列表
└── requirements.txt
```

**学位论文（课题）资产**位于本地 `d:\workdata\envmeta_thesis\`（非 git repo），
与本仓库严格分离：

```
d:\workdata\envmeta_thesis\        # 本地工作区，不入 GitHub
├── figures/                       # 论文图集（paper_en + thesis_cn + supplementary）
├── tables/                        # 论文表格（实验设计 / 测序 QC / 拼接统计）
├── config/                        # 学位论文绘图配置（matplotlib + ggplot2）
└── scripts/                       # 原研究分析脚本（R + python + shell + legacy）
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

## 当前进度（2026-05-15，v0.10 GPB-targeted draft + bioRxiv ready）

**v0.10 — Paper 3 GPB-targeted bioRxiv 初稿就绪**（2026-05-14 ~ 2026-05-15 session）：
- **4 篇论文 v2 投稿计划锁定**（详见 [`software/master/4_papers_master_plan_2026-05.md`](software/master/4_papers_master_plan_2026-05.md)）：
  P1 BITEB→**Waste Management** (Q2/5) / P2 ES&T→**ISME J** (Q1/15) /
  P3 iMeta→**GPB + SoftwareX 兜底** (Q1/15→Q3/3) / P4 RESB→**Microbiological Research** (Q2/5)
- **Paper 3 v0.10 GPB-targeted draft 完整**（commit `35730e8`）：
  - GPB reframing：Title / Highlights / Abstract / Intro v2 重写为 "hierarchical
    weight-of-evidence MCDA framework with reference implementation EnvMeta"（[`outline_gpb.md`](paper/manuscript/outline_gpb.md) 1375 行）
  - Methods §4.0 Algorithmic Framework：4 算法伪代码 + Notation N1 + Equation 1
    weight-of-evidence label assignment + Time complexity Table N2
  - bioRxiv-ready manuscript：[`manuscript_v0.10_gpb_draft.md`](paper/manuscript/manuscript_v0.10_gpb_draft.md)
    （1389 行 / 10370 字） + [`manuscript_v0.10_gpb_draft.docx`](paper/manuscript/manuscript_v0.10_gpb_draft.docx)（1.94 MB / TOC + 图嵌入）
  - 9 figures：4 placeholders 已生成（F1/F5/F6/F10 in [`paper/figures/bioRxiv_v0.10/`](paper/figures/bioRxiv_v0.10/)）+
    4 现有 figures 复用 proxy（F2 mockup_03 / F3 mockup_10 / F4 cal_vs_stress / F7 liu cycle）+ F8 scaling_curve + S_pert
  - References 48 Vancouver 条目（4 calibration DOI verified：Wei 2024 / Liu 2023 /
    Grettenberger 2021 / Ayala 2020）
  - Tools comparison methodology：[`paper/benchmarks/external/tools_comparison/methodology.md`](paper/benchmarks/external/tools_comparison/methodology.md)
    Type A 共享 + Type B 独家矩阵；Task 3b Windows blocker（Krona NTFS / Anvi'o
    无 FASTA）caveat 已 honest 文档化，WSL2 后续补
- **Paper 1 v5.7.2 outline 完整** + audit/ledger 双文件 + F 核实清单（并行 session
  完成；详见 `software/papers/paper1_reads_arsenic/paper1_progress.md`）
- 测试 **301/301 全绿**
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
- ✅ **v0.8.2 R 侧侧对照** 11 图全部完成 + RDA 数值与 R vegan 对齐
- ✅ **v0.9.0 假说评分对照实验** Paper 3 核心证据：
  - 4 Arm KEGG-curated calibration 全 STRONG（作者 / Liu / Grettenberger / Ayala）
  - n=3 stress test discrimination evidence（Grettenberger A 级 + Liu/Ayala B 级）
  - **领域中立性铁证**：cross-topic arsenate_reduction 在 2/2 无砷数据集双双 n=0 active MAGs ⭐
  - 第 6 类 claim `pathway_inactive`（Popperian falsifiability 主力）
  - 双层假说写作教程（[`docs/hypothesis_writing_guide.md`](docs/hypothesis_writing_guide.md) + [`paper/hypotheses/HYPOTHESIS_DESIGN_PRINCIPLES.md`](paper/hypotheses/HYPOTHESIS_DESIGN_PRINCIPLES.md)）
  - 引用审计 + 6 YAML 引用 metadata 修订（4 处错引透明纠正 + Sánchez-España 2008 verified + Auld 2017 → Dai 2014 + Méndez-García 2015 verified）
- ✅ **v0.9.1 dominance_score 字段** + Paper 3 写作素材：
  - `min_dominance_fraction` 硬阈值参数（pathway_active claim）+ evidence 字段（两 evaluator 信息透明）
  - **Stress test v2 B → A 级升级**：Liu 0.625 → 0.250；Ayala 0.455 → 0.182；3/3 现全 A 级 clean discrimination ⭐
  - Paper 3 三段完整草稿：Methods §4.6（~1450 字）+ Results §X（~800 字）+ Discussion §Y（~640 字）
  - Table 1 + Table 2 + Figure X 实物素材（PDF/PNG/SVG）](paper/figures/paper3_hypothesis_scoring/)
- ✅ **v0.9.1 outline 全套整合 + ImageGP 2 reframing + 性能 benchmark**（2026-05-09 paper-side session）：
  - Plan B 整合：3 段 Paper 3 草稿 + Performance 段落 → outline_imeta.md §5.4.8 / §5.2.4 / §5.2.8 / §5.3 / §5.7.1（commit `e1da1d9`）
  - 3-tier 性能 benchmark：sample real-dense + Liu real-sparse + Liu synthetic-dense（58 cells × 14 figures）；scaling_curve.{pdf,svg,png} + paper/benchmarks/performance.md + docs/performance_zh.md（commit `e1c69a6`）
  - **关键 finding**：cycle_diagram cost 由 N_pathway × N_env × 999 perm × spearmanr(N_sample) 主导，**与 N_MAG 几乎无关**；annotation breadth 才是真正主导
  - **ImageGP 2 ecosystem-extension reframing**（commit `3c9662b`）：发现 Chen 2024 *iMeta* 3:e239（Liu YX 是 EIC，COI 已声明）→ EnvMeta 全文从"vs 竞品"改为 "domain-specialist complement to iMeta visualization ecosystem"；致敬 14 个 sister tools；Abstract / Highlights / §5.1 Intro / §5.2.9 / §5.6 Acknowledgments / §5.7 References 全部重写
- ✅ **mock review v0.9.1 → v0.9.2 修订**（commit `b9c5699`）：
  - 5 大 Major Issues 全部温和解修订（4 Resolved + 1 honestly acknowledged）
  - 5 个 Minor Issues 修复：As₂S₃ 化学引用 / KB v1.1→v2.0 / 999-perm Anderson 2001 justify / Stolz 2006 编号 / user-study 措辞 tempering
  - Recommendation: Major Revision → **Minor Revision**（投稿接收概率 60% → 85%）
- ✅ **publication strategy archived**：[`paper/manuscript/publication_strategy.md`](paper/manuscript/publication_strategy.md)
  - 选定 Path A 改良版：EnvMeta 先发（bioRxiv preprint + iMeta 投稿）+ 课题论文并行起草
  - 时间表满足 2026 内 ≥ 1 篇接收 + 2027-05 前所有文章见刊的硬约束
  - bioRxiv 投稿 checklist + 风险评估完整
- ✅ **1-day perturbation analysis**（mock review v0.9.2 Major #1 auxiliary evidence；commit `5a21e5f`）：
  - Runner: [`tools/external_benchmarks/perturbation_analysis.py`](tools/external_benchmarks/perturbation_analysis.py)（within-element + cross-element 双模式，N=20 / mode / dataset）
  - **Cross-element headline finding**：Liu 2023 collapse 到 0/20 STRONG（median 0.000）→ 元素级 target 准确性是机制刚需
  - Within-element bound：mean 跌 25-48%，但 40-50% 仍 STRONG（与 §Y.1 KEGG-coverage-dependent 一致）
  - 集成 Methods §4.6.7 / Results §X.3 / Discussion §Y.3 limitation #1
  - 写作素材：[`paper/manuscript/perturbation_analysis_results.md`](paper/manuscript/perturbation_analysis_results.md)；figure: `paper/benchmarks/external/perturbation/perturbation_curve.{pdf,png,svg}`
- ✅ **Mock review v0.9.3 (post-perturbation)** 跑完（untracked，本地检查文件）：
  - Recommendation: **Minor Revision (acceptance-track)**（v0.9.2 Minor → v0.9.3 仍 Minor，但加了 acceptance-track 标记）
  - **v0.9.2 Major #1 → Resolved with auxiliary evidence**（perturbation 兑换成功）
  - **v0.9.3 新 Major**：Arm A perturbation asymmetry — 作者数据排除于 perturbation 控制，可补 partial Arm A perturbation（≈30 min compute）封口
  - 余 1 Major (3-dataset stress caveat 仍未修) + 8 Minor（5 carried + 3 new）= 2-4 天到 acceptance-ready
  - 详见 [`paper/manuscript/mock_review_v0.9.3_post_perturbation.md`](paper/manuscript/mock_review_v0.9.3_post_perturbation.md)（不 push）
- ✅ **Deployment hardening + 在线 demo 轻量化 + FAQ 文档**（2026-05-11 session）：
  - **Streamlit Cloud 崩溃修复**：用户 push v0.9.5 后 Cloud 触发 redeploy，'Oh no. Error running app'；根因 Cloud 默认 Python 升到 3.13 导致 numpy 2.x / pandas 3.x / scikit-bio 0.7 / scipy 1.17 wheel 不可用 → 新增 `runtime.txt` pin Python 3.11；验证 `/healthz` 返回 `{"status":"ok"}` 恢复 (commit `56d296a`)
  - **30 MAG 轻量化测试样例**：新建 [`tests/sample_data_demo/`](tests/sample_data_demo/)（含 [_build_demo_subset.py](tests/sample_data_demo/_build_demo_subset.py) 重建脚本）；抽样规则 = 14 keystone 全留 + 4 元素每个 ≥5 MAG + 按总丰度补足；abundance 168→30 行 / ko_long 5008→937 行 / 内存约 1/5，Streamlit Cloud 并发承载 4-5 倍提升
  - **app.py 切换 + UI 提示**：[app.py:317](app.py#L317) `_load_sample_data()` 优先 demo / 回退 full；[app.py:654](app.py#L654) expander 显示 "30 MAG 子集 测试样例" caption（本地装零影响）
  - **新增 [docs/FAQ.md](docs/FAQ.md)**：12 条常见问题 — Q1 等候时间 / Q2 'Oh no' 怎么办 / Q3 本地装链 / **Q4 并发限制（1-3 OK / 8+ OOM）** / Q5 demo 样例 vs full 区别 / Q6 session 刷新数据丢 / Q7-9 本地装坑 / Q10-12 数据识别 / 循环图慢 / INSUFFICIENT 解读
  - **README + README_CN**：顶部 nav 加 FAQ 链接；'Online demo' 段加 30 MAG 说明 + 并发限制 + Oh-no 处理引导；'Sample data' 章节展开为 full vs demo 两层介绍
  - **导师版大纲 docx**：[software/planning/Paper3_Outline_for_Supervisor.docx](software/planning/Paper3_Outline_for_Supervisor.docx)（45.8 KB / 中文宋体 + 英文 TNR 混合字体 / 13 章 / 含 5 个待决问题）+ 构建脚本 [build_paper3_outline_docx.py](software/planning/build_paper3_outline_docx.py)；software/ 在 .gitignore，本地存档
  - 测试 **301/301 全绿**；2 commit push: `56d296a` (runtime.txt) + `8851910` (demo subset + FAQ)
- ✅ **Mock review v0.9.5 (handling editor) 4 cheap Major 修订**（2026-05-10 third session）：
  - **Major #3 OpenTimestamps 区块链锚点已落地**：4 anchor commits (`42168da`, `44d7f5f`, `76a4f77`, `50c4687`) × 3 公共 calendar (alice / bob / finney) = 12 anchor 响应 archived 在 [`paper/manuscript/timestamps/`](paper/manuscript/timestamps/)；§4.6.2 现引用 OpenTimestamps anchoring 而非"future commit"——pre-registration 担忧从 institutional-trust-based 升级为 cryptographically witnessed
  - **Major #1 null_p 重 frame**：§4.6.1 明确把 null_p 标为 "shuffle-consistency diagnostic"，**explicitly NOT a frequentist p-value**；解释 4-9 claim YAML + 离散 satisfaction 让 null distribution 必然粗糙；加 weight-of-evidence 软件 relationship 段（MCDA / Bradford-Hill / Suter & Cormier / Linkov / Rhomberg refs）
  - **Major #4 Arm A reframe**："saturation regime" → "**partial-perturbation robustness check on the three pathway_active claims**"；明确说明 6 unperturbed claims 占总权重 ≈72%，arithmetic 锚定 score 在 0.75 以上；不再 over-claim 是 dataset saturation
  - **Minor #4 distance-to-boundary**：Results §X.2 加 stress score 到最近 label boundary 距离的明确数字（Grett −0.150；Liu −0.125 / +0.225；Ayala +0.055 near-boundary 与 §4.6.8 阈值 sensitivity transition 一致）；Table 2 列描述更新
  - 测试 **301/301 全绿**；预期 mock review v0.9.6 全部 6 Major 中 4 个 Closed (剩 #2 2-D band gap sweep 和 #5 引用扩 + #6 figures)
- ✅ **Mock review v0.9.4 (independent first-time reviewer) 5 Major 修订**（2026-05-10 second session）：
  - **Major #1 mock-review metadata 泄漏**：删除 manuscript 4 个文件全部 "Mock Review v0.9.x" 标签（§4.6.7 + perturbation_analysis_results §1/§4 + outline §5.7.1 reference annotations + outline §5.3 author notes）
  - **Major #2 Arm A baseline 数值不一致**（Table 1 = 1.000 vs perturbation = 0.919）：扩展 perturbation runner 支持 keystone_df，sample_data + keystone 跑出 1.000 STRONG 9/9 satisfied；perturbation 文件 + README 同步；Results §X.1 "4/4" 改为 "Arm A 9/9, C1/C2-A/C2-B 4/4"
  - **Major #3 pre-registration 第三方 witness 缺**：§4.6.2 加一段 "institutional-trust-based" 局限性诚实承认 + 3 条 mitigation（GitHub 公开 commit log / anchor commit margin / perturbation 独立性）+ 未来 OpenTimestamps 承诺
  - **Major #4 默认阈值 0.75/0.40 缺 justification**：新增 [`tools/external_benchmarks/threshold_sensitivity.py`](tools/external_benchmarks/threshold_sensitivity.py) 跑 5 阈值 × 8 dataset；新增 Methods §4.6.8 Threshold sensitivity 段；4 KEGG-curated arm 全 STRONG 跨 0.65-0.85 全稳；Arm B 全 INSUFFICIENT；stress 跨范围全 weak/suggestive；Results §X.1 加引用句
  - **Major #5 Arm A 循环论证 reframe**：Results §X.1 完全重写——明确把 Arm A 标为 "in-house positive control (engine self-consistency check)"，3 个 external arm 标为 "primary calibration evidence"；Discussion §Y.3 limitation #1 同步 reframe（Arm A 最 susceptible，3 external 较少 susceptible 但仍非 bias-free）
  - 测试 **301/301 全绿**；预期 mock review v0.9.5 全部 Major Resolved
- ✅ **Mock review v0.9.3 全部 Major + Minor 修订**（2026-05-10 session）：
  - **v0.9.3 Major #1 (Arm A asymmetry) → Closed**：补 Arm A partial perturbation（仅 3 个 pathway_active claim）
    - **Arm A 100% STRONG retention**（within + cross 都 20/20）→ 揭示 saturation regime；168 MAG 4-element 丰标注让任何通路替换都仍命中 active；honest framing 为该数据集 cherry-pick 防御依然依赖 §4.6.2 pre-registration + §Y.4 blind future work
    - 跨数据集形成 monotonic gradient：**Arm A 100% (saturated) → Grettenberger 30% → Ayala 15% → Liu 0% (focused)**，Wilson 95% CI 标注；这是 internal validity check
  - **v0.9.2 Major #2 (3-dataset caveat) → Closed**：§X.2 加 "2-of-2 within 3-dataset stress design" 一句话
  - **v0.9.3 新 Minor 3 个全 Closed**：§4.6.7 加 N=20 rationale + Wilson CI + run_null=False 辩护
  - **v0.9.2 Minor 5 个全 Closed**：calibration evidence disambig (§4.6.4) / consistent with 重复软化 / Newman 1997→1998 reconcile (§5.7.1 #24 改 Newman 1997 AEM) / DRAM-KEGG 区别 (§4.6.3) / dense annotation (§5.2.8) clarify "dense relative to KB target pathways"
  - **§5.7.1 numbering 平铺 22-entry**（v0.9.2 Minor #3 carried）
  - 测试 **301/301 全绿**；预期 mock review v0.9.4 全部 Major Resolved → bioRxiv ready

## Backlog（投稿前 + Phase 4）

### 🟥 投稿前必做（按 [`publication_strategy.md`](paper/manuscript/publication_strategy.md) Path A 改良版）

| 任务 | 工时 | 备注 |
|---|---|---|
| **6 张 placeholder figures**（F1 架构 / F3 S1-S3 流程 / F5 HTML 4-panel / F6 Bundle / F7 砷渣案例 / F10 vs Tools 截图比对）| 2-3 天 | 投稿前 mandatory；用户配截图 + 我可辅助排版 |
| **§5.7.1 Vancouver → iMeta 引用格式转换** | 30-45 min | Zotero 改 docx 时统一处理；不急 |
| **bioRxiv 投稿 dry-run + 实际投稿** | 30 min + 3 天审核 | 见 publication_strategy.md §"bioRxiv 投稿操作 checklist" |
| **EnvMeta 投 iMeta** | 2-3 h 投稿系统填表 | bioRxiv DOI 上线后即可 |
| **Verify Korehi 2014 / Mendez-Garcia 2015 真正 AMD diazotrophy 引用** | 10 min | 已 Verified Dai 2014 + Méndez-García 2015 → §5.7.1 #4 + #12（见 [`hypothesis_references_audit.md`](paper/manuscript/hypothesis_references_audit.md)）|
| **User study 数据回收分析**（条件性）| 1 周 | 问卷已发 2026-04-19；如投稿前 n ≥ 8 回收成功则加进 §5.6 Acknowledgments；否则 Acknowledgments 已 tempering 处理 |

### 🟧 推荐做（stop mock review → 6 figures → bioRxiv）

| 任务 | 工时 | 理由 |
|---|---|---|
| ~~v0.9.3 全部 Major + Minor~~ | ✅ 完成 | Arm A partial / 3-dataset caveat / 8 Minor 全清 |
| ~~v0.9.4 独立审稿人 5 Major~~ | ✅ 完成 | mock metadata 泄漏 / Arm A 0.919↔1.000 / pre-reg / threshold sensitivity / Arm A reframe |
| ~~v0.9.5 handling editor 4 cheap Major~~ | ✅ 完成 | OpenTimestamps / null_p reframe / Arm A partial 重新框架 / distance-to-boundary |
| ~~Streamlit Cloud 崩溃修复~~ | ✅ 完成 (2026-05-11) | runtime.txt pin Python 3.11；30 MAG demo subset + FAQ.md 双管齐下防并发崩 |
| **导师审阅大纲 docx** | ⏳ 用户操作 | [Paper3_Outline_for_Supervisor.docx](software/planning/Paper3_Outline_for_Supervisor.docx) 已生成；用户发给导师；等反馈 |
| **2-D band gap threshold sweep**（v0.9.5 Major #2，carry 到 revision）| 1-2h | 在 revision 阶段处理；不影响投稿 |
| **引用列表扩到 ≥40**（v0.9.5 Major #5，carry 到 revision）| 1-2h | iMeta 方法学论文典型密度；revision 时补 |
| **19 个 Minor**（v0.9.4 / v0.9.5）| ~2h | 投稿后 reviewer 反馈一并回应 |

### 🟡 加分项（投稿后 / 课题论文阶段做）

| 任务 | 工时 | 理由 |
|---|---|---|
| ~~`dominance_score` 字段~~ | ✅ 完成 v0.9.1 | engineering retrofit, 不是 independent validation（mock review v0.9.2 已重 frame）|
| ~~性能 benchmark~~ | ✅ 完成 v0.9.1 | 58 cells × 3 regimes，scaling figure ready |
| ~~盲法 stress test~~ | ⏸ 暂缓（不可行）| Discussion §Y.4 仍作为 future work |
| Q7 新手 UX 精调 | ~1h | HTML 控件分组 / SVG-JSON 下拉 / 3 步 onboarding |
| tool_comparison.md 细化 | 4h | 表格对比 ImageGP 2 / Anvi'o / Krona / 测序公司云（已在 §5.2.9 outline，需扩成正式表）|

### 🟢 论文接收后再做（Phase 4）

- **L3 插件框架**：用户上传 Python `analyze()` → 自动注册 GUI
- ~~**KB version 字段升 v2.0**~~：✅ 完成 mock review v0.9.2（`elements.json` v1.1 → v2.0）
- **S3.5-doc DOI 验证**：YAML 里 6 个理论文献 DOI 待用户手动过一遍
- **循环图 D3 编辑器**：用户直接拖拽编辑节点/通路

### 已知小 bug（待收口）

- "文件管理"首页同名覆盖策略（当前只处理新 name）
- matplotlib 中文字体在部分 Windows 环境显示方块（需要打包字体）

## 下次 session 建议起点

按优先级（5 轮 mock review iteration stop；Streamlit Cloud 修复 + demo 轻量化已落地）：

1. **等导师反馈** —— 用户已生成 [Paper3_Outline_for_Supervisor.docx](software/planning/Paper3_Outline_for_Supervisor.docx)，含 5 个待决问题
2. **6 张 placeholder figures**（F1/F3/F5/F6/F7/F10）—— 投稿 mandatory；用户配截图 + 我辅助排版；2-3 天
3. **bioRxiv 投稿** —— 6 月初前完成；3 天审核拿 DOI
4. **EnvMeta 投 iMeta** —— bioRxiv DOI 上线后立即投
5. **课题论文起草并行** —— As 形态重测启动（用户主导）+ EnvMeta bioRxiv DOI 引用
6. （可选 carry 到 revision 阶段）2-D band gap sweep / 引用扩到 ≥40 / 19 个 Minor

---

## 历史日志 / 自检报告

详见 [DEBUG_NOTES.md](DEBUG_NOTES.md) — 2026-04-13 至 2026-04-18 全部 session 日志、
阶段性自检（v1/v2）、产品定位讨论、技术决策记录、用户反馈处理。
