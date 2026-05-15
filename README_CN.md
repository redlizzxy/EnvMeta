# EnvMeta

**环境微生物宏基因组可视化分析平台**

[English README](README.md) · [在线体验](https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/) · [常见问题 FAQ](docs/FAQ.md) · [License](LICENSE)

EnvMeta 解决环境微生物博士生的核心痛点：**测序公司给了一堆表格，不知道哪个文件能做什么分析**。文件一键识别 + 14 种发表级图表 + 元素循环图自动推断 + 假说评分器 + 独立交互 HTML 导出。**全部开源免费、离线可用。**

> 🌐 **在线体验（零安装）**：<https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/>
>
> 首页点「📦 加载砷渣修复示例数据」→ 3 秒跑通全部 14 个分析。

---

## 📣 内测体验官招募中（有奖）

EnvMeta 正在收集第二轮内测反馈用于方法学论文（目标：iMeta / Bioinformatics）。**双轨设计，控制在 30 分钟内**：

| 你的背景 | 时长 | 内容 |
|---|---|---|
| **完全无生物背景**（其他专业本科 / 研究生） | **12-15 min** | 4 个简短任务（点击 / 拖拽 / 导出），不要求理解结果 |
| **有生物或生信基础** | **22-25 min** | 多 1 页深度任务（PCoA 解读 + 假说 YAML 评分） |

**🎁 你将获得**：
- 🧧 完成填写可领红包（详情扫码后看）
- 🏷 名字写入论文致谢（可匿名）
- 📊 一周内收到完整反馈汇总

**👇 直接参与**：

🔗 **问卷链接**：<https://wj.qq.com/s2/26579245/b28b/>

（点击上方链接直接打开问卷；二维码海报随论文 supplement 私下分发）

> 之前填过 v1 长版问卷的同学，**这次仍欢迎再填一次** v2 —— 问卷结构变了，新数据无法和旧版直接拼接。v2 真的快很多 🙇

---

## 为什么做这个

| 真实痛点 | 传统方案 | EnvMeta |
|---|---|---|
| 测序公司给了一堆 .tsv，不知道能做什么 | 自己写代码识别 | **拖入自动识别** 11 种文件类型 |
| 想画堆叠图 / PCoA / LEfSe，但不会写 R | 付费升级高级分析 / 啃代码 | **3 次点击** 出发表级 PDF |
| 想画元素循环图论证铁-砷-硫耦合 | 业界无公开工具，只能 PPT 手画 | **从 KO 注释自动推断** 4 元素 18 通路 |
| 评估"我的机制假说数据支不支持" | 主观陈述 | **YAML 评分器** + 置换 null_p + 权重敏感度 |
| 审稿人要复现 SI | PDF 附件 + 数据表 + 邮件来回 | **400 KB 独立 HTML**，审稿人浏览器交互操作 |

## 核心卖点（竞品对比）

| 能力 | Krona | Anvi'o | MicrobiomeAnalyst | 测序公司云平台 | **EnvMeta** |
|---|:-:|:-:|:-:|:-:|:-:|
| 元素循环图自动推断 | ❌ | ❌ | ❌ | ❌ | **✅ 独有** |
| 假说评分 + null_p + 权重敏感度 | ❌ | ❌ | ❌ | ❌ | **✅ 独有** |
| 独立离线交互 HTML（SI 杀手锏） | ❌ | ⚠️ static | ❌ web-only | ❌ | **✅ 独有** |
| Fork Bundle 论文-工具绑定 | ❌ | ❌ | ❌ | ❌ | **✅ 独有** |
| 跨元素耦合（As↔H₂S→As₂S₃） | ❌ | ❌ | ❌ | ❌ | **✅ 独有** |
| 任意参数可调 | ⚠️ | ✅ | ⚠️ | ❌ 参数锁死 | ✅ |
| 完全开源免费 | ✅ | ✅ | ✅ | ❌ 加钱 | ✅ |

## 功能矩阵（v0.9.0 / 假说 stress test 完成）

| 模块 | 支持内容 |
|------|---------|
| 📁 文件识别 | metadata / abundance（MAG/Taxon 分层） / distance / alpha / CheckM / env / KO 宽/长表 / Keystone / MAG taxonomy / Gephi nodes+edges — **11 种** |
| 📊 Reads-based（7 图） | 物种组成堆叠图 / α 多样性 / β 多样性 PCoA / RDA 排序 / LEfSe 差异分析 / 元素循环基因热图 / 基因 log2FC |
| 🧬 MAG-based（5 图） | MAG 质量散点 / MAG 丰度热图 / 通路完整度 / 元素循环基因谱 / 共现网络 Gephi-prep |
| 🔄 生物地球化学循环图 ⭐ | 4 元素（As/N/S/Fe）× 18 通路自动推断 + 跨元素化学物耦合 + keystone ★ 标注 |
| 🧪 机制假说 YAML 评分器 | **6 类 claim**（pathway_active / **`pathway_inactive`** [v0.9 ⭐ Popperian 否证] / coupling_possible / env_correlation / keystone_in_pathway / group_contrast）+ 3 可信度指标（Fisher 置换 p / Saltelli 权重敏感度 / Bradford Hill required veto）+ 9 档解读 |
| 📐 假说写作教程（v0.9 新增） | 双层模板（calibration + stress claims）+ pre-registration 纪律 + pre-prediction 模板 + 6 类 claim 选择指南 + Bradford-Hill 对应。见 [`docs/hypothesis_writing_guide.md`](docs/hypothesis_writing_guide.md) |
| 📦 Fork Bundle | 打包 KB + YAML + config + KEGG 快照 → zip，审稿人一键复现 |
| 🌐 独立交互 HTML | 400 KB 单文件 D3.js 嵌入，4 象限力导向 + 点击穿透 + SVG/JSON 导出，离线可用 |
| 💾 导出中心 | PNG / PDF / SVG / TIFF 600dpi / TSV / 可运行 `.py` 复现脚本，批量 ZIP |
| 🧭 新手落地包 | 数据准备指南（11 上游工具 → EnvMeta 映射）+ 图表选择向导（问题 → 推荐分析）+ 14 图「如何解读」expander + 样例数据一键加载 |

## 目标用户

- 🎓 环境微生物 / 土壤修复 / 生态方向的**硕士、博士研究生**
- 📚 送样到测序公司的课题组，想自己跑下游分析
- 📝 发论文需要自动推断元素循环图作为 mechanism figure
- 🔬 不要求编程经验；有 R/Python 基础可获得更深的定制能力

## 三种使用方式

### 1. 在线零安装（最简单）

<https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/>

- 首次访问若显示「App is sleeping」，点「Wake up」等 **30-90 秒**。
- 在线 demo 默认载入 **30 MAG 测试样例**（`tests/sample_data_demo/`），
  非论文完整 168 MAG 数据集——够你点点看体验功能，但**不用于复现论文**。
  严肃使用请本地装。
- Streamlit Cloud 免费版只有 1 GB 内存：1-3 人并发 OK，8 人以上容易 OOM。
  详见 [FAQ Q4](docs/FAQ.md#q4几个人能同时用在线版会被挤崩吗)。
- 若网页显示「😟 Oh no. Error running app」→ 见 [FAQ Q2](docs/FAQ.md#q2网页显示-oh-no-error-running-app-怎么办)。

### 2. 本地安装（新手版）

见 [docs/install_for_beginners.md](docs/install_for_beginners.md)，零基础 20 分钟装完（装 Miniconda + 克隆 repo + 一键装依赖）。

### 3. 本地安装（生信基础）

```bash
git clone https://github.com/redlizzxy/EnvMeta.git
cd EnvMeta
conda create -n envmeta python=3.11 -y
conda activate envmeta
pip install -r requirements.txt
streamlit run app.py
```

浏览器自动打开 `http://localhost:8501`。

## 📜 更新日志（近期）

内测期频繁修 bug / 加功能，完整列表见 **[CHANGELOG.md](CHANGELOG.md)**。

### v0.9.0 — 2026-05-09（假说评分对照实验完成 + Stress test discrimination 证据）⭐

**Paper 3 投稿核心证据全部就位**。在 4 个 KEGG-curated 宏基因组数据集（作者数据 + Liu 2023 冷泉 + Grettenberger 2021 AMD 溪流 + Ayala 2020 pit lake）上的对照实验 —— 全部 hypothesis YAML 在跑 EnvMeta 之前 commit（git timestamp 锚定）。

- ✨ **新增 `pathway_inactive` claim type** — 第 6 类 claim，Popperian falsifiability 主力。`n_active_mags == 0 → satisfied`（符合"不应活跃"预期）；不破坏现有 YAML（向后兼容）。
- ✨ **双层假说写作教程** — [`docs/hypothesis_writing_guide.md`](docs/hypothesis_writing_guide.md) 用户教程（calibration + stress 双层模板、pre-registration 纪律、pre-prediction 防 hindsight bias、6 类 claim 选择指南、Bradford-Hill 对应）。配套的论文方法学版（`HYPOTHESIS_DESIGN_PRINCIPLES.md`）与跨数据集 stress runner 随论文 supplement 私下归档，不进公开仓库。
- 📊 **4 个 KEGG-curated 数据集 calibration 全 STRONG** + **3 个 stress test 分数显著低于 calibration**（Grettenberger weak 0.250 / Liu suggestive 0.625 / Ayala suggestive 0.455）。Cross-topic `arsenate_reduction_should_dominate` 在 **2/2 无砷数据集**（Grettenberger + Ayala）双双 reject（n=0 active MAGs）—— EnvMeta 评分领域中立性铁证。
- 📚 **引用 DOI 审计** — verify 16 claim × 13 文献 DOI；透明纠正 4 处错引（Yin 2011 期刊错；Bothe "2007 FEMS Rev" 不存在 → 应为 Bothe 2000；Cabrera 2006 期刊+主题错；Auld 2017 主题错 → 替换为 Dai 2014 PLoS One [primary AMD nifHDK metagenomic 实证] + Méndez-García 2015 Front Microbiol review）。
- 🐛 6 个 hypothesis YAML 引用 metadata 修订（不动 claim 实体；pre-registration audit trail 在 git history 完整保留）。
- 🧪 pytest **297/297 全绿**（+4 个 `pathway_inactive` 测试）。

### v0.8.2 — 2026-05-08（RDA 数值对齐 R vegan + 11 图侧侧对照完成）

- 🐛 修 RDA 数值与 R vegan 不一致（skbio 归一化差异致 inertia 16-20× 偏差、ANOVA F/p 反转）
  - 改用 SS-based 公式（vegan-equivalent），修复后 F / r / 解释度 4 位精度对齐 R
- 📚 R/Python 11 图侧侧对照工作完成（验证产物随论文 supplement 私下归档）
  - 5 图数值精确一致 + 6 图算法等价 + 11 个 README + 论文引用模板
- 📚 Paper 3（EnvMeta 方法学论文）投稿前 4 大任务清单（随论文私下归档）
- 🧪 pytest 293/293 全绿（无回归）

### v0.8.1 — 2026-04-21（Mac 端首批内测反馈修复）

- 🐛 修 HTML 交互导出切组后化学物-通路连线失效（拖拽不跟随、hover 不高亮）
- 🐛 修 macOS 上传 YAML 假说评分报 `[Errno 63] File name too long`
- 📚 小白安装指南补 Mac 3 大坑 FAQ（conda ToS / protobuf resolver / Streamlit 欢迎邮箱）
- 📚 README 新增"更新到新版本"章节

### v0.8.0 — 2026-04-19（v0.8 内测版 "Sunday Sprint"）

- ✨ 新手落地包（数据准备指南 / 图表向导 / 14 图解读 / 一键加载样例）
- ✨ 导出中心 4-tab 统一入口 + 批量 ZIP
- ✨ HTML 交互导出 v1.3（D3 inline 嵌入 + 独立离线）
- ✨ Streamlit Cloud 在线部署

---

## 🔄 更新到新版本（内测阶段）

内测阶段 EnvMeta 会**频繁修 bug / 加功能**，建议每周 pull 一次拿最新代码。

**在线版用户**无需任何操作，云端会自动部署。

**本地安装用户**，在 Terminal（Mac）或 Anaconda Prompt（Windows）里执行：

```bash
# 1) 如果 streamlit 正在跑，先 Ctrl+C 停掉
# 2) 进入 EnvMeta 目录
cd ~/Desktop/EnvMeta              # Mac
# cd %USERPROFILE%\Desktop\EnvMeta  # Windows

# 3) 激活环境
conda activate envmeta

# 4) 拉最新代码（这就是核心更新命令）
git pull origin master

# 5) 仅在终端提示 requirements 变化 / 启动报 ModuleNotFoundError 时才跑
pip install -r requirements.txt

# 6) 重新启动
streamlit run app.py
```

### 常见场景

**场景 A：`git pull` 提示 `Already up to date.`**
已经是最新，直接启动即可。

**场景 B：`git pull` 报 `Your local changes ... would be overwritten`**
说明你在本地改过文件。先暂存再拉：
```bash
git stash              # 暂存本地改动
git pull origin master
git stash pop          # 如果想恢复本地改动
```

**场景 C：更新后启动报错 `No module named 'xxx'`**
跑一次 `pip install -r requirements.txt` 补齐新依赖。

### 怎么知道有没有更新

- 看 [GitHub 仓库首页](https://github.com/redlizzxy/EnvMeta)的 commits 列表，有新 commit 就是有更新
- 或 `git fetch && git log HEAD..origin/master --oneline` 查看本地落后的 commit 数
- 遇到 bug 前建议先 pull 一次，可能早已修好

## 数据格式

读 [docs/data_preparation_zh.md](docs/data_preparation_zh.md)（也可在 app 内「数据准备指南」页面浏览）：覆盖 CoverM / HUMAnN3 / eggNOG / DRAM / GTDB-Tk / CheckM2 / KofamScan / QIIME2 / Kraken2+Bracken / MetaPhlAn4 / iCAMP 等 11 种上游工具的输出 → EnvMeta 输入格式映射。

## 样例数据

两层样例：

- **`tests/sample_data/`** — 砷渣-钢渣微生物修复研究精简数据，**168 MAG ×
  10 样本 × 57 目标 KO**，CK/A/B 三组。用于 pytest 单元测试、论文 Arm A
  positive-control 校准、扰动分析 runner。
- **`tests/sample_data_demo/`** — 上述的 **30 MAG 子集**，保留全部 14 个
  keystone MAG + 每元素 ≥5 个活跃 MAG。在线版「📦 加载砷渣修复测试样例」按钮
  加载这个子集，降低 Streamlit Cloud 并发崩溃风险。**仅用于功能演示，不反映
  原研究的科学结论**。

本地装默认优先用 demo 子集；若 demo 目录不存在则回退到完整版。
重建 demo 子集：`python tests/sample_data_demo/_build_demo_subset.py`。

## 技术栈

- **前端 UI**：Streamlit
- **出版级渲染**：matplotlib + seaborn
- **统计**：pandas / scipy / scikit-bio / statsmodels
- **交互 HTML**：D3.js v7（inline 嵌入，无外部依赖）
- **循环图推断**：自研（`envmeta/geocycle/`）— KEGG-driven KB + 置换检验 + 敏感度扫描

## 项目结构

```
envmeta/
├── app.py                         # Streamlit 入口
├── envmeta/                       # 核心包
│   ├── file_manager/              # 文件识别（11 种类型）
│   ├── analysis/                  # 14 图分析引擎
│   ├── geocycle/                  # 循环图推断 + 假说评分 + HTML 导出
│   │   ├── knowledge_base/        # 4 元素 × 18 通路 × 57 KO（KEGG-driven）
│   │   ├── hypothesis.py          # YAML 评分器
│   │   └── html_exporter.py       # 独立交互 HTML
│   ├── help/                      # 新手落地包（向导 / 解读 / 反向索引）
│   ├── tools/                     # Fork Bundle / KB builder / Gephi prep
│   └── export/                    # PNG/PDF/SVG/TIFF + .py 复现脚本
├── docs/
│   ├── data_preparation_zh.md     # 上游工具 → EnvMeta 映射
│   └── install_for_beginners.md   # 小白安装指南
├── tests/
│   ├── sample_data/               # 论文精简数据（一键加载）
│   └── test_*.py                  # 291 case 全绿
└── requirements.txt
```

## 开发进度

- [x] **Phase 0** — 项目骨架 + 环境 + 知识库
- [x] **Phase 1** — 7 张 Reads-based 图表 + 基础调参 + 导出
- [x] **Phase 2** — 5 张 MAG-based 图表 + 代码生成器
- [x] **Phase 3** — 循环图推断 + 假说评分 + Fork Bundle + 独立交互 HTML
- [x] **v0.8 Sunday Sprint** — 新手落地包 + 导出中心统一 + HTML v1.3
- [x] **English README + LICENSE**（本次发布）
- [ ] Zenodo DOI（v1.0 release 时挂载）
- [ ] 论文 Methods + Results + Discussion 起草
- [ ] Phase 4 — 插件框架（论文接收后）

## 论文引用

EnvMeta 计划作为方法学论文发表。

发表前引用本仓库 URL 即可；发表后会在此处追加 DOI。

## 致谢

EnvMeta 源自一项砷渣-钢渣微生物修复宏基因组研究。感谢内测阶段的同学贡献反馈（名单在论文发表后见致谢页）。

## License

[MIT License](LICENSE) — Copyright (c) 2026 redlizzxy and EnvMeta contributors.

## 联系（建议优先邮件联系）

- GitHub Issues：<https://github.com/redlizzxy/EnvMeta/issues>
- Bug 反馈 / 功能请求：18872605913@163.com

