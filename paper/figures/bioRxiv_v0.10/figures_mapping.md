# Paper 3 v0.10 bioRxiv — Figure Asset Mapping

> 2026-05-14 — manuscript v0.10 bioRxiv 投稿前 figure 资产盘点 + 占位生成清单

---

## 图号 ↔ 资产映射表

| Figure | 内容 | 资产路径 | 状态 |
|---|---|---|---|
| **F1** | Five-tier framework architecture | `paper/figures/bioRxiv_v0.10/F1_architecture.png` | 🟡 placeholder 待生成（本目录） |
| **F2** | 14 figures + GUI parameter tuning | `paper/figures/mockups/03_hierarchy_demo.png` (proxy) + `paper/figures/screenshot_*.png` | 🟢 可用 proxy（mockup + screenshots） |
| **F3** | Cycle inference algorithm flowchart | `paper/figures/mockups/10_cascade_with_coupling.png` | 🟢 可用 proxy（Mockup 10）|
| **F4** | Hypothesis scorer + 4-arm calibration | `paper/figures/paper3_hypothesis_scoring/figure_x_calibration_vs_stress.png` | ✅ **已就绪**（real data figure）|
| **F5** | Standalone HTML 4-panel | `paper/figures/bioRxiv_v0.10/F5_HTML_4panel.png` | 🟡 placeholder 待生成 |
| **F6** | Fork Bundle structure | `paper/figures/bioRxiv_v0.10/F6_Bundle.png` | 🟡 placeholder 待生成 |
| **F7** | Case study (arsenic-slag) | `paper/benchmarks/external/liu_2023_coldseep/envmeta_outputs/fig1_cycle_diagram.png` (proxy) | 🟢 可用 proxy（Liu cold-seep cycle） |
| **F8** | Performance scaling envelope | `paper/benchmarks/performance/scaling_curve.png` | ✅ **已就绪**（58-cell benchmark）|
| **F10** | Vs Tools head-to-head | `paper/figures/bioRxiv_v0.10/F10_vs_tools.png` | 🟡 placeholder 待生成 |
| **S_pert** | Perturbation curve (SI) | `paper/benchmarks/external/perturbation/perturbation_curve.png` | ✅ **已就绪** |

**总结**：3/9 已就绪（F4 / F8 / S_pert）+ 4/9 可用 proxy（F2 / F3 / F7 + mockups 资产）+ **4/9 需生成 placeholder**（F1 / F5 / F6 / F10）。

---

## 4 张需生成 placeholder 的内容设计

每张 placeholder = simple matplotlib + text layout，bioRxiv 预印阶段足够；
后续 GPB submission 前用户配真实截图 / 用 Inkscape 配色完善。

### F1 — Five-tier framework architecture

布局：5 层堆叠水平条 + 顶部 Streamlit GUI 层 + 数据流箭头

文字内容：
- L1: General inference engine (S1+S2+S3+null calibration)
- L2: 6-claim YAML schema (pathway_active / inactive / coupling / env_correlation / keystone_in_pathway / group_contrast)
- L3: Plugin framework (post-acceptance)
- L4: Fork Bundle distribution (zip = KB + YAML + config + KEGG snapshot + data)
- L5: KEGG-driven biogeochemical KB (4 elements × 18 pathways × 57 KOs)
- Top: Streamlit GUI orchestration (upload → recognize → analyze → export)

### F5 — Standalone HTML 4-panel SI

布局：2x2 panel grid，每 panel 标签 + 内容描述

- A: Biogeochemical cycle (D3 force simulation, 4 quadrants, draggable nodes)
- B: Hypothesis scoring (sortable claim table, null_p distribution)
- C: Cross-group comparison (CK / A / B side-by-side)
- D: Parameter audit (all S1-S3 intermediate quantities)

中部标注："~400 KB / fully offline / inline D3.js v7"

### F6 — Fork Bundle structure

布局：左 .zip 容器图标 + 右 5 个组件块（KB / YAML / config / KEGG snapshot / sample data）+ 下方 "Load → 5-min reproduction" 流程箭头

文字：
- Bundle.zip ≈ 5-15 MB
- Avoids "version drift" + 100% reproducibility of published figures
- Distribution model: fork-rather-than-community

### F10 — Vs Tools head-to-head comparison

布局：3x3 matrix — 行 = analysis task (PCoA / KO heatmap / Hierarchical HTML) × 列 = tool (EnvMeta / Krona / MicrobiomeAnalyst)

文字：
- 每 cell 标注：tool 名 + 该 task 的输出格式 / 操作步数 (TBD pending Task 3b)
- 第 4 行（unique to EnvMeta）：cycle inference / hypothesis scoring / standalone HTML SI

---

## Placeholder 生成脚本

见 [`generate_placeholders.py`](generate_placeholders.py)（同目录）— matplotlib 简单 box + text 布局，跑一遍生成 F1/F5/F6/F10 四张 PNG。
