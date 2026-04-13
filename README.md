# EnvMeta

**环境微生物宏基因组可视化分析平台**

EnvMeta 是一个面向环境微生物研究（污染修复、元素循环）的宏基因组下游可视化分析工具，解决"从数据文件到发表级图形"的最后一公里问题。核心创新是**一键生成生物地球化学循环图**。

---

## 特色

| 维度 | 现有工具 | EnvMeta |
|------|---------|---------|
| 使用形式 | R Shiny / 命令行 | 本地 Web 应用（streamlit run） |
| 分析范围 | 通用微生物组 | **环境微生物专用**（污染修复、元素循环） |
| 可视化定制 | 代码级调参或有限 GUI | **GUI 实时调参 + 自动生成可复现代码** |
| 机制图 | 无 | **自动生成生物地球化学循环图** |

## 功能模块

- **A · 文件管理器** — 拖入文件自动识别类型（丰度表 / KO / CheckM / 距离矩阵 / metadata 等）
- **B · 可视化引擎** — 12 种图表：堆叠图 / α·β多样性 / PCoA / RDA / LEfSe / 基因热图 / log2FC / MAG质量 / 丰度热图 / 通路完整度 / 基因谱 / 共现网络
- **C · 交互调参** — 颜色、字体、图例、布局实时调整
- **D · 高质量导出** — PDF / SVG / TIFF / PNG + 可运行 Python 脚本
- **E · 循环图生成器** ⭐ — 从 KO 注释推断活跃通路，生成可编辑的元素循环机制图

## 安装

需要 Python 3.11+。推荐使用 conda。

```bash
conda create -n envmeta python=3.11 -y
conda activate envmeta
pip install -r requirements.txt
```

详细安装步骤见 [INSTALL.md](INSTALL.md)。

## 启动

```bash
conda activate envmeta
streamlit run app.py
```

浏览器打开 http://localhost:8501 即可使用。

## 目录结构

```
envmeta/
├── app.py                   # Streamlit 入口
├── envmeta/                 # 核心包
│   ├── file_manager/        # 模块 A：文件识别
│   ├── analysis/            # 模块 B：分析引擎
│   ├── params/              # 模块 C：参数管理
│   ├── export/              # 模块 D：导出
│   └── geocycle/            # 模块 E：循环图引擎
├── tests/
│   └── sample_data/         # 论文精简样例数据
├── paper/                   # 论文数据积累（验证/效率对比/用户反馈）
├── scripts/                 # 原始论文分析脚本（参考实现）
├── data/                    # 本地数据（不入库）
├── docs/
├── requirements.txt
├── CLAUDE.md                # 项目说明（Claude Code 读取）
└── INSTALL.md
```

## 技术栈

- 前端 UI：Streamlit
- 交互预览：Plotly
- 出版级渲染：matplotlib
- 后端：pandas / scipy / networkx / scikit-bio
- 循环图编辑器：D3.js（嵌入 st.components）

## 开发路线图

- [x] Phase 0 — 项目骨架 + 环境 + 知识库框架
- [ ] Phase 1 — 文件识别 + 堆叠图/PCoA/热图 + 基础调参 + 导出
- [ ] Phase 2 — 全部 12 种图表 + 代码生成器
- [ ] Phase 3 — 循环图 v1（推断引擎 + 静态渲染）
- [ ] Phase 4 — 循环图交互编辑 + v1.0 发布

## 样例数据

`tests/sample_data/` 提供来自砷渣-钢渣微生物修复研究（CK/A/B 三组）的精简数据，用于功能验证和演示。详见 [tests/sample_data/README.md](tests/sample_data/README.md)。

## 引用

EnvMeta 计划作为方法学论文发表，目标期刊：iMeta / Bioinformatics / Frontiers in Microbiology。

## License

TBD
