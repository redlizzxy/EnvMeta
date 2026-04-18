# EnvMeta

**环境微生物宏基因组可视化分析平台**

EnvMeta 解决环境微生物博士生的核心痛点：**测序公司给了一堆表格，不知道哪个文件能做什么分析**。文件一键识别 + 14 种发表级图表 + 元素循环图自动推断 + 假说评分器 + 独立交互 HTML 导出。**全部开源免费、离线可用。**

> 🌐 **在线体验（零安装）**：<https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/>
>
> 首页点「📦 加载砷渣修复示例数据」→ 3 秒跑通全部 14 个分析。

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

## 功能矩阵（v0.8 / Phase 1+2+3 全通）

| 模块 | 支持内容 |
|------|---------|
| 📁 文件识别 | metadata / abundance（MAG/Taxon 分层） / distance / alpha / CheckM / env / KO 宽/长表 / Keystone / MAG taxonomy / Gephi nodes+edges — **11 种** |
| 📊 Reads-based（7 图） | 物种组成堆叠图 / α 多样性 / β 多样性 PCoA / RDA 排序 / LEfSe 差异分析 / 元素循环基因热图 / 基因 log2FC |
| 🧬 MAG-based（5 图） | MAG 质量散点 / MAG 丰度热图 / 通路完整度 / 元素循环基因谱 / 共现网络 Gephi-prep |
| 🔄 生物地球化学循环图 ⭐ | 4 元素（As/N/S/Fe）× 18 通路自动推断 + 跨元素化学物耦合 + keystone ★ 标注 |
| 🧪 机制假说 YAML 评分器 | 5 类 claim（pathway_active / coupling_possible / env_correlation / keystone_in_pathway / group_contrast）+ 3 可信度指标（Fisher 置换 p / Saltelli 权重敏感度 / Bradford Hill required veto）+ 9 档解读 |
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

首次访问若显示「App is sleeping」，点「Wake up」等 30 秒。

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

## 数据格式

读 [docs/data_preparation_zh.md](docs/data_preparation_zh.md)（也可在 app 内「数据准备指南」页面浏览）：覆盖 CoverM / HUMAnN3 / eggNOG / DRAM / GTDB-Tk / CheckM2 / KofamScan / QIIME2 / Kraken2+Bracken / MetaPhlAn4 / iCAMP 等 11 种上游工具的输出 → EnvMeta 输入格式映射。

## 样例数据

`tests/sample_data/` 提供砷渣-钢渣微生物修复研究（CK/A/B 三组，168 MAG × 10 样本 × 57 目标 KO）的精简数据。首页按钮一键加载即可跑通全部 14 分析。

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
├── paper/
│   ├── bundles/                   # 论文 Fork Bundle 示例
│   ├── benchmarks/                # 验证数据 + 效率对比
│   └── user_study/                # 评测问卷设计
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
- [ ] 论文 Methods + Results + Discussion 起草
- [ ] English README + LICENSE + Zenodo DOI（iMeta 投稿硬指标）
- [ ] Phase 4 — 插件框架（论文接收后）

## 论文引用

EnvMeta 计划作为方法学论文发表。目标期刊：**iMeta / Bioinformatics / Frontiers in Microbiology**。

发表前引用本仓库 URL 即可；发表后会在此处追加 DOI。

## 致谢

EnvMeta 源自一项砷渣-钢渣微生物修复宏基因组研究。感谢内测阶段的师弟师妹贡献反馈（名单见致谢页）。

## License

TBD（投稿前确定；预计 MIT 或 Apache-2.0）。

## 联系

- GitHub Issues：<https://github.com/redlizzxy/EnvMeta/issues>
- Bug 反馈 / 功能请求：欢迎 PR
