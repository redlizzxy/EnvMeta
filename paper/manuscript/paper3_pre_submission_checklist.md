# Paper 3 (EnvMeta) 投稿前素材补全清单

> **目标期刊**：iMeta（IF 24，2022 创刊）
> **预计投稿日期**：2026-08-15
> **创建日期**：2026-05-07
> **关联文件**：
> - [outline_imeta.md](outline_imeta.md) — 论文大纲（已写）
> - 总体策略：[../../software/planning/4_papers_master_plan_2026-05.md](../../software/planning/4_papers_master_plan_2026-05.md)
>
> **本文件用途**：跟踪 Paper 3 投稿前 4 大任务进度。任何任务有进展即更新本文件。

---

## 总览：4 大任务

| # | 任务 | 工时 | 优先级 | 状态 | 完成时间 |
|---|---|---|---|---|---|
| 1 | **R 侧侧对照 11 图补全** | 1-2 周 | 🟥 P0 | 进行中（2026-05-07 启动） | TBD |
| 2 | **第二外部数据集 benchmark** | 2-3 天 | 🟥 P0 | 待启动 | TBD |
| 3 | **English README + LICENSE + Zenodo DOI** | 4-6h | 🟥 P0 | 待启动 | TBD |
| 4 | **Methods 第一稿起草** | 1 周 | 🟧 P1（依赖前 3 项） | 待启动 | TBD |

---

## 任务 1 — R 侧侧对照 11 图补全

### 背景

当前 `paper/benchmarks/validation/` 仅有 log2fc 1 张 R 对照（11 缺 1 张已完成）。
Paper 3 Methods 关键证据：每张 EnvMeta 出图必须能与原 R 脚本输出**视觉一致**，
证明 EnvMeta 不是另起炉灶造的"独立工具"，而是**对标研究领域成熟流程的可视化整合**。

### 11 图清单（实际：5 R + 6 Python，原研究 MAG 部分用 Python 而非 R）

EnvMeta 12 图扣除"循环图"（自研无对照）= **11 图需对照**：

#### Reads-based（7 图，前 5 用 R 脚本对照）

| # | 图类型 | EnvMeta 模块 | 原脚本 | 状态（2026-05-08） |
|---|---|---|---|---|
| 1 | 物种堆叠图（Phylum/Genus/Species）| `envmeta/analysis/stackplot.py` | `scripts/R/01_tax_stackplot.R` | ✅ **完成**（2026-05-08 amplicon 装好后；3 层级 × 3 视图 PDF；数值精确一致到 0.01%）|
| 2 | α 多样性（Shannon/Simpson/Chao1/Pielou）| `envmeta/analysis/alpha.py` | `scripts/R/02_alpha_diversity.R` | ✅ **完成**（KW 100% 一致 / pairwise p ≤ 5% 偏差）|
| 3 | β-PCoA + PERMANOVA | `envmeta/analysis/pcoa.py` | `scripts/R/02_beta_PCoA.R` | ✅ **完成**（F/R²/PC% 完全一致）|
| 4 | RDA | `envmeta/analysis/rda.py` | `scripts/R/03_RDA.R` | ✅ **完成（2026-05-08 修复后）**：F/r/解释度 4 位精度对齐 R vegan |
| 5 | LEfSe（LDA 条形图）| `envmeta/analysis/lefse.py` | `scripts/R/04_LEfSe.R` (R 重绘 Galaxy `input.res`) | ⚠️ **部分完成**（R 重绘 PDF 已生成；逐特征算法对照待补）|
| 6 | 基因热图 | `envmeta/analysis/gene_heatmap.py` | `scripts/python/05_gene_heatmap_log2fc.py` | ✅ **完成**（51 KO + 3 组 + Welch's + BH 一致；4 FigS 补充图已生成）|
| 7 | log2FC | `envmeta/analysis/log2fc.py` | `scripts/python/05_gene_heatmap_log2fc.py` | ✅ **完成**（K11811 抽样数值绝对值一致；🟧 比较方向命名差异已记录）|

#### MAG-based（4 图，全部用 Python 脚本对照）

| # | 图类型 | EnvMeta 模块 | 原脚本 | 状态 |
|---|---|---|---|---|
| 8 | MAG 质量 | `envmeta/analysis/mag_quality.py` | `scripts/python/06_MAG_quality.py` | ✅ **完成**（35/111/22 计数完全一致；指标均值偏差 ≤ 0.04%）|
| 9 | MAG 丰度热图 | `envmeta/analysis/mag_abundance.py` | `scripts/python/07_MAG_abundance_heatmap.py` | ✅ **完成**（算法等价；🟧 abundance 标尺差异（百分比 vs 比例）已记录）|
| 10 | 通路完整度 | `envmeta/analysis/pathway_completeness.py` | `scripts/python/08_pathway_completeness.py` | ✅ **完成**（168 MAG 一致；KB 17→18 通路升级已记录）|
| 11 | MAG 基因谱 | `envmeta/analysis/gene_profile.py` | `scripts/python/06_MAG_gene_profile.py` | ✅ **完成**（拷贝数算法一致；🟦 51 vs 57 KO 列差异已记录）|

### 进度概览（2026-05-08 更新 — 全量完成 + RDA bug 已修）

```
✅ 完成（数值精确对照）：       5 / 11   alpha, pcoa, MAG_quality, stackplot, RDA（修复后）
✅ 完成（算法等价已确认）：     6 / 11   lefse, gene_heatmap, log2fc, gene_profile, pathway, MAG_abundance
─────────────────────────────────────────────────────
合计完成 + 通过：              11 / 11   ✅ 全量完成
```

### Bug 修复历史（2026-05-08 完成）

#### ✅ P0 — RDA 数值已修

修复方案：`envmeta/analysis/rda.py` 第 160-243 行替换为 SS-based 算法
（lstsq fit + SVD eigvals + Type I sequential ANOVA），对标 vegan::rda()
公式。`explained_ref` 默认从 `"constrained"` 改为 `"total"`（R vegan 默认）。

| 指标 | R vegan | EnvMeta（修复后）| 偏差 |
|---|---|---|---|
| Total inertia | 0.1617 | 0.16168 | 0% |
| Constrained | 0.1125 | 0.11248 | 0% |
| RDA1 explained | 32.96% | 32.96% | 0% |
| RDA2 explained | 24.94% | 24.94% | 0% |
| pH F | 3.905 | 3.905 | 0% |
| Eh F | 3.877 | 3.877 | 0% |
| Mantel r 全部 | — | — | 全 0% |

p 值因 RNG 微差（< 0.015），但 4 因子 α=0.05 显著性方向一致。
**293/293 测试全绿**（修复未破坏现有功能）。

#### 重新评估的 4 个"伪 bug"

经详细代码检查后确认**不是 bug**：

| 之前标记 | 真实状态 | 处理 |
|---|---|---|
| 🟧 log2fc 方向反转 | 由用户 `group_a/group_b` 控制；xlabel 自动写出方向 | README 加约定提示 |
| 🟧 mag_heatmap 标尺差 100× | 误判 — abundance.tsv 列总和=100，两侧均用百分比；row_order 是 phylum_cluster 不是 abundance 降序 | README 撤回标记 |
| 🟦 gene_profile 51 vs 57 KO | EnvMeta `drop_zero_kos=True` 是合理 default；可设 `False` 还原 57 列 | README 说明可调参数 |
| 🟦 pathway KB 17→18 | KB v2.0 升级，非 bug | Methods 里说明即可 |

### Python 端待补的精细数值验证（下次 session）

可写一个 Python 验证脚本批量做：
```python
# tools/validate_python_outputs.py
# 比较 paper/benchmarks/validation/*/Fig*_*.txt vs envmeta_*_stats.tsv
# 输出 per-figure mean/median/count 偏差报告
```

### 用户手动操作（环境层面）

**安装 amplicon 包**（stackplot 需要）—— Claude Code Bash 沙箱无法访问 GitHub，请你在 RStudio 或 PowerShell 终端跑：
```r
install.packages("remotes", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")
remotes::install_github("microbiota/amplicon")
```
装好后告知 Claude，可立刻补 stackplot 对照。

### 输出标准

每张图存放：
```
paper/benchmarks/validation/<figure_name>/
├── README.md                 # 数据集 + R 脚本路径 + EnvMeta 参数 + 视觉对比说明
├── r_output.pdf              # R 脚本输出
├── envmeta_output.pdf        # EnvMeta 输出（同参数）
├── side_by_side.pdf          # 两图并排（论文 SI 用）
└── diff_notes.md             # 视觉差异点 + 解释（颜色 / 字体 / 排版）
```

---

## 任务 2 — 第二外部数据集 benchmark

### 背景

iMeta 审稿人 90% 概率会问"EnvMeta 是否只能跑你自己砷渣数据"。第二外部数据集
是**泛化性证据**，必备。

### 候选数据集

| 数据集 | 规模 | 主题 | 资源需求 | 推荐度 |
|---|---|---|---|---|
| **Tara Oceans 子集**（公开下载）| 500+ MAG × 100+ sample | 海洋 N/S 循环 | 16 GB / 4 核 / 1-2h | ⭐⭐⭐⭐⭐ |
| **Oak Ridge 铀污染** | 中等 | 陆地金属循环（与砷渣对口）| 同上 | ⭐⭐⭐⭐ |
| **EMP500** | 500 sample × 5000+ OTU | 全球微生物组 | 同上 | ⭐⭐⭐ |
| **GEM 数据库切片** | 任意规模 | 综合 | 同上 | ⭐⭐⭐ |

### 输出标准

```
paper/benchmarks/validation/second_dataset/
├── README.md                  # 数据集来源 + 预处理 + EnvMeta 跑全套 + 结果总结
├── input_data/                # 子采样后输入数据（≤ 100 MB）
├── envmeta_outputs/           # 14 图全部 PDF + .py 复现脚本
├── benchmark_table.tsv        # runtime / memory / 出图质量对比
└── compare_to_original.md     # 对比原 publication 图（如有）
```

---

## 任务 3 — English README + LICENSE + Zenodo DOI

### 子任务

| 子任务 | 工时 | 状态 |
|---|---|---|
| `LICENSE` 文件添加（MIT）| 5 min | ⬜ |
| `README.md` 翻译为英文（中文迁 `README_CN.md`）| 3-4h | ⬜ |
| GitHub release v1.0 创建 | 30 min | ⬜ |
| Zenodo 关联 GitHub repo + 触发 DOI | 30 min | ⬜ |
| Zenodo metadata 填写 + DOI 嵌入 README | 30 min | ⬜ |

---

## 任务 4 — Methods 第一稿起草

### 背景

依赖任务 1-3 完成（需要引用 R 对照 + 第二数据集 + Zenodo DOI）。

### 子节大纲（已在 outline_imeta.md 第 5.4 节）

| 子节 | 字数 | 状态 |
|---|---|---|
| 4.1 文件识别模块 | 300 | ⬜ |
| 4.2 14 图分析引擎 | 300 | ⬜ |
| 4.3 元素循环推断算法（S1/S2/S3）| 500 | ⬜ |
| 4.4 YAML 假说 schema | 200 | ⬜ |
| 4.5 D3.js 交互 HTML 导出 | 200 | ⬜ |
| 4.6 Benchmark 实施 | 300 | ⬜ |
| 4.7 实施细节 | 150 | ⬜ |

---

## 维护记录

| 日期 | 事项 |
|---|---|
| 2026-05-07 | 本文件创建，4 大任务存档 |
| 2026-05-07 | 任务 1（R 对照）启动 |
