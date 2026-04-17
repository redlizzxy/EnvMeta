# 操作效率对比记录

| 分析类型 | 传统方式（代码行数） | 传统方式（步骤数） | 传统方式（耗时） | 一键复现 | EnvMeta（步骤数） | EnvMeta（耗时） | 完成日期 |
|---------|-------------------|------------------|----------------|---------|-----------------|----------------|---------|
| 堆叠图 | 177 (R, 01_tax_stackplot.R) | 5+（写脚本/装 R 包/调参/导出） | ~30 min（含调试） | ✅ pip install + sample_data | 3（上传 2 文件 + 选样式 + 点生成） | ~10 s | 2026-04-13 |
| PCoA | 195 (R, 02_beta_PCoA.R) | 4+（装 vegan/ggrepel/cowplot/写脚本） | ~40 min | ✅ pip install | 3（上传 2 文件 + 点生成） | ~12 s（999 置换） | 2026-04-13 |
| α多样性 | 260 (R, 02_alpha_diversity.R) | 5+（vegan/ggpubr + 2×3 子图布局） | ~35 min | ✅ pip install | 3（上传 2 文件 + 勾指数 + 点生成） | ~5 s | 2026-04-14 |
| RDA排序 | 264 (R, 03_RDA.R) | 6+（vegan/ggrepel/Mantel 循环） | ~60 min（含样本 ID 对齐） | ✅ pip install（skbio RDA + Mantel） | 3（上传 3 文件 + 点生成） | ~2 s | 2026-04-15 |
| LEfSe | 206 (R, 04_LEfSe.R) + Galaxy LEfSe 外部依赖 | 7+（装 Galaxy/R/dplyr + 预跑 + 写脚本） | ~75 min（含 Galaxy 环境准备） | ✅ pip install（内置 KW + LDA，无需 Galaxy） | 3（上传 2 文件 + 调阈值 + 点生成） | ~3 s | 2026-04-16 |
| 基因热图 | 514 (py, 05_gene_heatmap_log2fc.py, 含 log2FC 部分) | 6+（装 matplotlib/调整颜色/位置/字体） | ~90 min | ✅ pip install + 知识库 JSON | 3（上传 2 文件 + 点生成） | ~8 s | 2026-04-13 |
| log2FC | 含在 514 行 05_gene_heatmap_log2fc.py（后半段 ~230 行） | 5+（装 scipy/statsmodels + 配色） | ~50 min | ✅ pip install + 知识库 JSON | 4（上传 2 文件 + 选组对 + 点生成） | ~3 s | 2026-04-14 |
| MAG质量 | 287 (py, 06_MAG_quality.py) | 5+（argparse + bac120/ar53 合并 + 配色/图例布局）| ~45 min | ✅ pip install | 3（上传 2-3 文件 + 点生成）| ~2 s | 2026-04-17 |
| MAG丰度热图 | 544 (py, 07_MAG_abundance_heatmap.py) | 6+（argparse + 硬编码样本组 + 三段配色设计 + bac120/ar53 合并）| ~90 min | ✅ pip install（metadata 驱动组色带）| 3（上传 1-4 文件 + 调 Top-N + 点生成）| ~5 s | 2026-04-20 |
| 通路完整度 | 970 (py, 08_pathway_completeness.py) | 7+（argparse + bac120/ar53 合并 + 3 张图硬编码布局）| ~120 min | ✅ pip install（KB 自动读 18 通路）| 3（上传 2-4 文件 + 样式切换 + 点生成）| ~3 s | 2026-04-17 |
| 基因谱 | 1080 (py, 06_MAG_gene_profile.py) | 7+（argparse + iTOL 解析 + 两语言版 + 4 张子图）| ~150 min | ✅ pip install（KB 自动提供 KO 顺序）| 3（上传 2-4 文件 + 点生成）| ~3 s | 2026-04-17 |
| 共现网络（Gephi 辅助） | 415 (py, 09_cooccurrence_network.py) + Gephi 手动调 | 7+（写脚本 + Gephi 布局 + 手动删标签 + 调参 + 导出）| ~120 min（含 Gephi 反复调）| ✅ pip install + gephi_nodes/edges CSV | 3（上传 2 CSV + 点生成散点 + 一键导出 Gephi 就绪 CSV）| ~3 s（散点图） + 推荐参数省 30 min Gephi 调参 | 2026-04-17 |
| 循环图 | **无现成脚本**（首次实现）| N/A — 需手工画 | ~数小时手工概念图 | ✅ pip install（全自动推断）| 3（上传 + 点生成）| ~5 s | 2026-04-17 |

> **一键复现**列：记录"他人拿到代码/项目后能否一键复现该图"（是/否/需配置环境）。
