"""研究问题决策树 — 「我该用哪张图？」向导。

8 大研究问题 × 每类 2-3 细问题 = ~18 节点，全部映射到已实现的 14 个分析。
app.py 的「图表选择向导」页面按此数据渲染 radio 两步 + 推荐卡片。

数据结构：
    NAVIGATOR: list[CategoryNode]
    CategoryNode = {
        "category": str,                  # 大类名（显示在一级 radio）
        "icon": str,                      # emoji
        "description": str,               # 一句话说明大类
        "subquestions": list[SubQuestion],
    }
    SubQuestion = {
        "q": str,                         # 具体研究问题
        "recommended": list[Recommendation],
        "required_files": list[str],      # FileType.value 字符串
        "tip": str,                       # 可选，小样本 / 特殊情况提醒
    }
    Recommendation = {
        "analysis_id": str,               # 对应 ANALYSIS_INPUTS 键
        "reason": str,                    # 为什么推荐
        "priority": int,                  # 1=首选, 2=次选, 3=补充
    }
"""
from __future__ import annotations


NAVIGATOR: list[dict] = [
    # ─── 1. 组间差异 ───────────────────────────────────────────
    {
        "category": "组间差异分析",
        "icon": "📊",
        "description": "比较 CK/处理组或不同样点间哪些特征有显著差异",
        "subquestions": [
            {
                "q": "哪些物种在组间相对丰度差异显著？",
                "recommended": [
                    {"analysis_id": "lefse", "reason": "LDA 效应量 + 组特异性排序，文献经典方法",
                     "priority": 1},
                    {"analysis_id": "stackplot", "reason": "直观看组间优势物种切换（配合统计使用）",
                     "priority": 2},
                ],
                "required_files": ["abundance_wide", "metadata"],
                "tip": "样本量 <5/组 时 Kruskal-Wallis 可能普遍 >0.05，建议看 LDA 排序而非 p 值",
            },
            {
                "q": "哪些功能基因在组间丰度差异显著？",
                "recommended": [
                    {"analysis_id": "log2fc", "reason": "log2FC + Welch's t-test + BH 校正，基因层面经典方法",
                     "priority": 1},
                    {"analysis_id": "gene_heatmap", "reason": "组均值热图直观看整体模式",
                     "priority": 2},
                ],
                "required_files": ["ko_abundance_wide", "metadata"],
                "tip": "log2FC 需选两组对比（A vs B），多组建议两两对比",
            },
            {
                "q": "哪些 MAG 在组间丰度差异显著？",
                "recommended": [
                    {"analysis_id": "mag_heatmap", "reason": "Top-N MAG 丰度热图 + 门色块，一目了然",
                     "priority": 1},
                    {"analysis_id": "log2fc", "reason": "若关注显著性排序（用 MAG 丰度矩阵做输入）",
                     "priority": 2},
                ],
                "required_files": ["abundance_wide", "mag_taxonomy"],
                "tip": "推荐配 Layer 1 filter_mode=variance 看组间分化大的 MAG",
            },
        ],
    },
    # ─── 2. 群落结构相似性 ───────────────────────────────────
    {
        "category": "群落结构相似性",
        "icon": "🌐",
        "description": "判断不同组/样本的整体微生物群落结构是否有差异",
        "subquestions": [
            {
                "q": "组间群落结构是否显著分离？",
                "recommended": [
                    {"analysis_id": "pcoa", "reason": "PCoA + PERMANOVA 是文献标准流程",
                     "priority": 1},
                ],
                "required_files": ["distance_matrix", "metadata"],
                "tip": "若没有距离矩阵，可从物种丰度用 Bray-Curtis 计算；PERMANOVA p<0.05 表示显著分离",
            },
            {
                "q": "组间群落组成整体有何变化？",
                "recommended": [
                    {"analysis_id": "stackplot", "reason": "Top-N 物种堆叠图，看组间优势物种切换",
                     "priority": 1},
                    {"analysis_id": "pcoa", "reason": "降维可视化 + 统计检验双配",
                     "priority": 2},
                ],
                "required_files": ["abundance_wide", "metadata"],
            },
        ],
    },
    # ─── 3. α/β 多样性 ───────────────────────────────────────
    {
        "category": "α/β 多样性",
        "icon": "🌱",
        "description": "比较不同组的物种丰富度 / 均匀度 / 多样性指数",
        "subquestions": [
            {
                "q": "组间 α 多样性（Shannon / Simpson / Chao1）是否有差异？",
                "recommended": [
                    {"analysis_id": "alpha_boxplot", "reason": "箱线图 + Kruskal-Wallis + 两两 Mann-Whitney",
                     "priority": 1},
                ],
                "required_files": ["alpha_diversity", "metadata"],
                "tip": "样本量 <4/组 时建议只看中位数趋势，不强求显著性",
            },
            {
                "q": "组间 β 多样性（距离）是否分离？",
                "recommended": [
                    {"analysis_id": "pcoa", "reason": "PCoA + PERMANOVA 经典流程",
                     "priority": 1},
                ],
                "required_files": ["distance_matrix", "metadata"],
            },
        ],
    },
    # ─── 4. 功能通路活性 ─────────────────────────────────────
    {
        "category": "功能通路活性",
        "icon": "🧬",
        "description": "看样本/MAG 中元素循环通路是否活跃（KO 覆盖 + 丰度）",
        "subquestions": [
            {
                "q": "组间元素循环基因整体丰度差异？",
                "recommended": [
                    {"analysis_id": "gene_heatmap", "reason": "As/N/S/Fe 四元素 KO × Sample 合并热图",
                     "priority": 1},
                    {"analysis_id": "log2fc", "reason": "基因层面显著性 + 效应量排序",
                     "priority": 2},
                ],
                "required_files": ["ko_abundance_wide", "metadata"],
            },
            {
                "q": "哪些 MAG 能完整跑通哪些通路？",
                "recommended": [
                    {"analysis_id": "pathway", "reason": "MAG × 通路完整度热图 / 气泡图，直观看功能承载者",
                     "priority": 1},
                ],
                "required_files": ["ko_annotation_long"],
                "tip": "completeness ≥ 50% 一般视为'通路在这个 MAG 里完整'",
            },
            {
                "q": "每个 MAG 携带哪些元素循环 KO？",
                "recommended": [
                    {"analysis_id": "gene_profile", "reason": "MAG × KO 拷贝数热图，按元素分块排列",
                     "priority": 1},
                ],
                "required_files": ["ko_annotation_long"],
                "tip": "默认 viridis + blank_zeros，0 值留白避免视觉噪声",
            },
        ],
    },
    # ─── 5. 物种共现网络 ─────────────────────────────────────
    {
        "category": "物种共现网络",
        "icon": "🕸️",
        "description": "识别关键物种（keystone）+ 物种间共现模式",
        "subquestions": [
            {
                "q": "谁是 keystone（度中心性 + 中介中心性双高）？",
                "recommended": [
                    {"analysis_id": "network", "reason": "Degree vs Betweenness 散点图 + 阈值线，keystone 可视化",
                     "priority": 1},
                ],
                "required_files": ["gephi_nodes", "gephi_edges"],
                "tip": "本软件只出散点图 + Gephi 预处理 CSV。完整网络图请在 Gephi 里画（更专业）",
            },
            {
                "q": "如何把结果导入 Gephi 画网络图？",
                "recommended": [
                    {"analysis_id": "network", "reason": "提供 Gephi CSV 预处理（非 keystone Label 自动清空）+ 参数指南",
                     "priority": 1},
                ],
                "required_files": ["gephi_nodes", "gephi_edges"],
                "tip": "在「共现网络图」页底部 expander 有 Gephi 推荐参数和操作步骤",
            },
        ],
    },
    # ─── 6. 元素循环 + 机制验证 ─────────────────────────────
    {
        "category": "元素循环 + 机制验证",
        "icon": "⚛️",
        "description": "EnvMeta 差异化卖点：自动推断元素循环图 + 机制假说评分",
        "subquestions": [
            {
                "q": "样品/组的元素循环结构长什么样？",
                "recommended": [
                    {"analysis_id": "cycle_diagram", "reason": "从 KO 注释自动推断活跃通路 + 合并细胞级联可视化",
                     "priority": 1},
                ],
                "required_files": ["ko_annotation_long"],
                "tip": "支持 As/N/S/Fe 四元素 + 跨元素化学物耦合（如 As(III)↔H₂S→As₂S₃）",
            },
            {
                "q": "我的研究假说是否被数据支持？",
                "recommended": [
                    {"analysis_id": "hypothesis_score", "reason": "YAML 写假说 → 5 类 claim 打分 + null_p + 权重敏感度 + required-veto",
                     "priority": 1},
                ],
                "required_files": ["ko_annotation_long"],
                "tip": "需准备 YAML 假说文件（参考 paper/hypotheses/ 示例）；配 Fork Bundle 一键复现",
            },
            {
                "q": "跨组元素循环有什么差异？",
                "recommended": [
                    {"analysis_id": "cycle_diagram", "reason": "「循环图」页内切组下拉 + 跨组对比表 + 最活 ★ 标记",
                     "priority": 1},
                ],
                "required_files": ["ko_annotation_long", "metadata"],
            },
        ],
    },
    # ─── 7. MAG 质量和丰度 ───────────────────────────────────
    {
        "category": "MAG 质量和丰度",
        "icon": "🦠",
        "description": "检查 MAG 组装质量 + 看哪些 MAG 丰度高",
        "subquestions": [
            {
                "q": "我的 MAG 质量如何？（High/Medium/Low）",
                "recommended": [
                    {"analysis_id": "mag_quality", "reason": "Completeness vs Contamination 散点图 + 阈值标注",
                     "priority": 1},
                ],
                "required_files": ["checkm_quality"],
                "tip": "默认阈值 High=Comp≥90&Cont≤5；MIMAG 标准；可在侧边栏调阈值",
            },
            {
                "q": "哪些 MAG 丰度最高？组间怎么变化？",
                "recommended": [
                    {"analysis_id": "mag_heatmap", "reason": "Top-N MAG × Sample 热图 + 三段非线性配色（长尾分布适配）",
                     "priority": 1},
                ],
                "required_files": ["abundance_wide"],
                "tip": "Layer 1 filter_mode=top_plus_keystone 推荐（Top-N 并集 keystone）",
            },
        ],
    },
    # ─── 8. 环境因子-群落关系 ───────────────────────────────
    {
        "category": "环境因子-群落关系",
        "icon": "🌾",
        "description": "识别哪些环境变量驱动群落结构变化",
        "subquestions": [
            {
                "q": "哪些环境因子显著解释群落变化？",
                "recommended": [
                    {"analysis_id": "rda", "reason": "RDA/CCA 约束排序 + Mantel 逐因子检验",
                     "priority": 1},
                ],
                "required_files": ["abundance_wide", "env_factors", "metadata"],
                "tip": "env_factors 文件需包含 SampleID + Group + ≥2 数值因子；env 与群落样本 ID 要完全对应",
            },
            {
                "q": "单一环境因子（如 Eh / pH）和某通路是否相关？",
                "recommended": [
                    {"analysis_id": "cycle_diagram", "reason": "循环图的 env-pathway Spearman 相关面板 + 置换检验可信度",
                     "priority": 1},
                ],
                "required_files": ["ko_annotation_long", "env_factors"],
                "tip": "循环图会自动标注 strong / suggestive / spurious? 三档可信度（已做 999 次置换）",
            },
        ],
    },
]


def all_analysis_ids() -> set[str]:
    """返回 NAVIGATOR 引用的全部 analysis_id（用于测试数据完整性）。"""
    ids: set[str] = set()
    for cat in NAVIGATOR:
        for sq in cat["subquestions"]:
            for rec in sq["recommended"]:
                ids.add(rec["analysis_id"])
    return ids
