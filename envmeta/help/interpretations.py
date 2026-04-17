"""图表解读文案 — 14 个分析各一段"如何解读"引导。

结构参考假说评分器的 9 档解读模式：

    title             标题
    what_it_shows     这张图回答什么研究问题
    how_to_read       list: 关键视觉元素的含义（轴 / 颜色 / 标注）
    good_signal       什么样的结果可作为正面证据
    warning           常见误判 / 陷阱
    caveats           方法学局限 / 样本量建议

app.py 在每个分析页图形下方调用 `_render_interpretation_expander(analysis_id)`
统一渲染。

使用规范：
- 每条 150-300 字为宜，避免超长
- "good_signal" 用文献经典结果表述，不涉及具体数据集结论
- 保持描述性 / 工具性语气，不下因果结论（与 cycle_diagram 设计一致）
"""
from __future__ import annotations


INTERPRETATIONS: dict[str, dict] = {
    # ═══════════════════════════════════════════════════════════
    # Reads-based
    # ═══════════════════════════════════════════════════════════
    "stackplot": {
        "title": "如何解读物种组成堆叠图",
        "what_it_shows": (
            "每个样本（或组均值）中 Top-N 物种的**相对丰度**分布。是微生物"
            "群落最直观的总览图，适合第一张给读者看。"
        ),
        "how_to_read": [
            "色块高度 = 该物种在样本中的相对丰度（总和 = 100%）",
            "同一颜色跨样本指代**同一物种**（颜色图例在右侧）",
            "灰色 `Others` 条 = Top-N 之外的所有物种之和",
            "x 轴可以是单样本（style=sample）或组均值（style=group）",
        ],
        "good_signal": (
            "明显的组间优势物种切换（如 CK 组 A 属占 40% → 处理组 B 属占 40%）"
            "或 Others 占比 <30%（说明 Top-N 能代表主要群落）。"
        ),
        "warning": (
            "若 Others 占比 > 50%，提示 Top-N 设置过低；建议调到 15-20。"
            "若某样本颜色分布与其他同组样本差异巨大，检查是否数据质量问题或"
            "测序深度不均。"
        ),
        "caveats": (
            "堆叠图展示**相对**丰度，不能直接读出绝对量。若要定量对比，"
            "配合 PCoA + PERMANOVA（整体结构）或 LEfSe（具体物种）使用。"
        ),
    },
    "pcoa": {
        "title": "如何解读 β 多样性 PCoA 图",
        "what_it_shows": (
            "把高维群落距离矩阵降到 2 维散点图，**视觉判断组间是否分离**；"
            "配合 PERMANOVA p 值做统计判断。"
        ),
        "how_to_read": [
            "每个点 = 一个样本；颜色 = 分组",
            "两轴标签 `PC1 (xx.x%)` 中的百分比 = 该轴解释的方差比例",
            "同组样本**聚成一簇**且组间分开 → 群落结构差异显著",
            "图右上角 PERMANOVA 结果：R² + p 值",
        ],
        "good_signal": (
            "PERMANOVA p < 0.05 且 R² > 0.2；视觉上同组样本聚团、不同组明显分离。"
            "PC1 + PC2 合计解释率 > 40% 比较理想。"
        ),
        "warning": (
            "p > 0.05 但视觉上分离 → 可能样本量不足；n < 4/组时 PERMANOVA 检验力很弱。"
            "R² < 0.1 即使 p 显著也说明分组对群落结构解释力有限。"
        ),
        "caveats": (
            "PCoA 是**无约束**排序（不考虑环境变量）。若想知道哪些环境因子驱动群落"
            "变化，请用 RDA/CCA 约束排序。距离选择影响结果：Bray-Curtis 重丰度、"
            "UniFrac 重系统发育，按研究问题选。"
        ),
    },
    "gene_heatmap": {
        "title": "如何解读元素循环基因热图",
        "what_it_shows": (
            "As / N / S / Fe 四大元素循环基因（KO）在各样本的 TPM 热图。"
            "一图看元素循环功能基因的整体丰度模式 + 跨元素对比。"
        ),
        "how_to_read": [
            "行 = KEGG KO（按元素分块 + 通路分组）；列 = 样本（按组排列）",
            "颜色 = log(TPM) 或 z-score（侧边栏可切）；左侧色条标注 KO 所属元素",
            "整行/整列横色 → 该 KO 普遍存在或该样本所有基因高表达",
            "显著组间差异 → 查看左右半图色条对比",
        ],
        "good_signal": (
            "特定元素的基因在处理组明显激活（如 S 循环基因在加硫组变深）、"
            "或特定通路（如 arsC-grx/trx）跨样本一致表达 → 说明功能稳定且响应明确。"
        ),
        "warning": (
            "某 KO 全样本都接近 0 → 可能测序深度不够 / KO 注释过于严格。"
            "z-score 模式下容易把低丰度组内变异放大成视觉噪声，建议先看 log(TPM) 模式。"
        ),
        "caveats": (
            "这张图展示**有多少拷贝**，不是**通路能否完整跑通**。后者请用「代谢通路"
            "完整度」图（基于 MAG × KO 长表计算）。"
        ),
    },
    "alpha_boxplot": {
        "title": "如何解读 α 多样性箱线图",
        "what_it_shows": (
            "每组样本的 Shannon / Simpson / Chao1 等指数分布，**判断处理是否改变**"
            "**群落丰富度或均匀度**。"
        ),
        "how_to_read": [
            "箱体 = Q1-Q3，中线 = 中位数；须 = 1.5×IQR；小点 = 单个样本",
            "Shannon / Simpson 越大 → 越均匀；Chao1 越大 → 越丰富",
            "箱之间横线 + `*` 标星 = 两两 Mann-Whitney U 检验显著",
            "Kruskal-Wallis 总体 p 在图顶（若全组一起检验）",
        ],
        "good_signal": (
            "某组中位数明显高/低于其他组，两两比较有 * 或 **，且趋势与研究假说一致"
            "（如处理组 Shannon 降低 = 选择压力 / 生态位缩窄）。"
        ),
        "warning": (
            "样本量 <4/组 时几乎不可能达到 p<0.05（统计功效不足），此时"
            "**以中位数趋势为主**，不强求显著性。若箱体非常大 → 组内异质性高，"
            "先排查采样/测序批次问题。"
        ),
        "caveats": (
            "α 多样性衡量**单样本内部**；跨样本差异用 β 多样性（PCoA）看。同一样本"
            "用不同指数可能结论不同（Chao1 敏感稀有类群，Simpson 敏感优势类群）。"
        ),
    },
    "log2fc": {
        "title": "如何解读基因差异分析（log2FC）图",
        "what_it_shows": (
            "两组对比下每个 KO 的 log2(倍数变化) + 显著性。快速找出**哪些基因在**"
            "**处理组显著上调/下调**。"
        ),
        "how_to_read": [
            "x 轴 = log2FC；正值 = 组 A 高；负值 = 组 B 高",
            "柱长 = 效应量大小；星号 = padj 显著性（`*` padj<0.05，`**` <0.01，`***` <0.001）",
            "按 4 元素（As/N/S/Fe）分 2×2 子图，柱按 |log2FC| 降序",
            "灰色柱 = padj 不显著（视觉淡化，但仍展示供参考）",
        ],
        "good_signal": (
            "多个同通路 KO 同方向变化（如 arsC 和 arsB 都上调）→ 说明通路级响应"
            "而非单基因噪声。|log2FC| > 2 且 padj < 0.05 通常作为主要差异证据。"
        ),
        "warning": (
            "小样本（n<5/组）下 padj 普遍 >0.05 是正常的（多重检验严厉）。"
            "此时**以 log2FC 效应量排序为主**，并在论文里诚实汇报样本量限制。"
            "pseudocount 避零会把 0-vs-非零变化膨胀成大 log2FC，看到极端值要检查原始丰度。"
        ),
        "caveats": (
            "Welch's t-test 假设近似正态；高度零膨胀的 KO（大部分样本为 0）"
            "结果可能不稳。此时 LEfSe（KW + LDA）更稳健。"
        ),
    },
    "rda": {
        "title": "如何解读 RDA/CCA 排序图",
        "what_it_shows": (
            "有环境变量**约束**的排序，回答「哪些环境因子显著解释群落变化」。"
            "PCoA 是无约束的；RDA/CCA 是有约束的。"
        ),
        "how_to_read": [
            "每点 = 样本；每箭头 = 一个环境因子",
            "**箭头方向** = 该因子梯度方向；**箭头长度** = 解释力度（越长越重要）",
            "箭头间夹角 <30° → 因子高度相关（多重共线）；~90° → 独立；~180° → 负相关",
            "样本点投影到箭头上的位置反映该样本的环境特征值",
        ],
        "good_signal": (
            "前两轴累计解释率 > 30%；Mantel 逐因子检验中有 2-3 个因子 p<0.05 且"
            "箭头长 → 说明该因子是群落变化的关键驱动。同组样本沿某箭头方向聚集 → 该因子"
            "响应一致。"
        ),
        "warning": (
            "env_factors 文件的样本 ID 必须与群落矩阵完全对应（不一致会报错或错误对齐）。"
            "环境因子数量 ≥ 样本数时会过拟合，RDA 约束变得不可解释；建议 n_env ≤ n_sample/2。"
        ),
        "caveats": (
            "RDA 假设线性响应（适合梯度不长的场景）；梯度长（DCA 轴长 > 4 SD）应换 CCA。"
            "与 Mantel 搭配使用：RDA 可视化 + Mantel 给逐因子 p 值。"
        ),
    },
    "lefse": {
        "title": "如何解读 LEfSe 差异分析图",
        "what_it_shows": (
            "找出在组间丰度**显著且效应量大**的生物标志物（biomarker）。LDA 值"
            "综合考虑显著性和效应量，是文献最常用的差异分析方法之一。"
        ),
        "how_to_read": [
            "每条横柱 = 一个分类单元；柱长 = log10(LDA 效应量)",
            "柱颜色 = 该特征在哪一组富集（组颜色图例在图例区）",
            "按组分区、组内按 LDA 降序排列",
            "默认阈值：KW p < 0.05 AND LDA > 2；侧边栏可调",
        ],
        "good_signal": (
            "每组有 5-20 个显著 biomarker；跨分类层级（phylum / class / ... / species）"
            "呈现一致模式（如多个同科物种都在 A 组富集）→ 说明响应有系统性。"
        ),
        "warning": (
            "样本量 <5/组 时 KW 几乎不显著，可放宽到 α=0.1 看趋势。"
            "若某组 biomarker 数量压倒其他组 → 检查是否测序深度不均匀导致偏差。"
        ),
        "caveats": (
            "本软件内置简化 LEfSe（KW + 组均值最大组 + log10(1+1e6×差值) 近似 LDA）。"
            "若论文需要严格复现 Segata 2011 方法，请装 Galaxy LEfSe 做二次验证。"
        ),
    },
    # ═══════════════════════════════════════════════════════════
    # MAG-based
    # ═══════════════════════════════════════════════════════════
    "mag_quality": {
        "title": "如何解读 MAG 质量评估散点图",
        "what_it_shows": (
            "所有 MAG 的完整度 × 污染度散点图，**按 MIMAG 标准**分 High / Medium /"
            " Low 三级，判断可用 MAG 数量。"
        ),
        "how_to_read": [
            "每点 = 一个 MAG；x=Completeness, y=Contamination",
            "点颜色 = 所属门（Phylum）；keystone 以菱形高亮",
            "虚线 = 质量等级阈值（默认 High: Comp≥90 ∧ Cont≤5；Medium: Comp≥50 ∧ Cont≤10）",
            "右上角文字 = 各级 MAG 数量统计",
        ],
        "good_signal": (
            "High 级 MAG 占比 > 30% 说明组装质量好。门分布多样化（10+ 门）说明群落"
            "覆盖广。High + keystone 重合点（高质量 + 功能关键）值得重点分析。"
        ),
        "warning": (
            "High 级 <10 个 → 组装深度不够，下游分析慎用 MAG 水平结论。"
            "某门只有一个 MAG 且质量低 → 分类结论要谨慎（可能是错误合并产物）。"
        ),
        "caveats": (
            "MIMAG 标准用 CheckM / CheckM2 计算；不同工具给出的 Completeness 可能相差"
            "5-10%。Contamination 只反映同一 MAG 内多基因组污染，不捕捉跨 MAG 冗余。"
        ),
    },
    "mag_heatmap": {
        "title": "如何解读 MAG 丰度热图",
        "what_it_shows": (
            "Top-N MAG × Sample 丰度热图，**看哪些 MAG 在哪些样本中丰度高**。"
            "配三段非线性配色适应长尾分布（少数极高 + 多数低值）。"
        ),
        "how_to_read": [
            "行 = MAG（默认按门聚类，门内层次聚类排序）；列 = 样本",
            "颜色 = 相对丰度 %（Blues→YlGn→YlOrRd 三段色，对应低/中/高）",
            "左侧色块 = 门；右上小星号 = keystone 物种",
            "列顶部（若有 metadata）= 组色块（CK/A/B 配色）",
        ],
        "good_signal": (
            "同门内 MAG 在某组集体富集 → 该门响应处理。Top-N 中 keystone 占 20-40%"
            "说明关键物种与高丰度物种显著重合。"
        ),
        "warning": (
            "某 MAG 只在一个样本爆发式丰度（其他接近 0）→ 可能是非典型测序批次问题。"
            "若大部分 MAG 颜色都在蓝色段（低丰度）→ Top-N 过小，建议调到 30-50。"
        ),
        "caveats": (
            "MAG 丰度通常是 read recruit 覆盖度（如 CoverM）；不同工具的绝对值不可直接"
            "跨研究比较，需同样本内归一化。"
        ),
    },
    "pathway": {
        "title": "如何解读代谢通路完整度图",
        "what_it_shows": (
            "每个 MAG 在每条通路的**完整度**（KO 覆盖比例）。回答「谁能跑完整条通路」"
            "这个 MAG 水平的功能问题。"
        ),
        "how_to_read": [
            "行 = MAG；列 = 通路（按元素 As/N/S/Fe 分组）",
            "heatmap 样式：颜色 = 完整度 %",
            "bubble 样式：圆大小 = 完整度，颜色 = 贡献（completeness × abundance）",
            "左色条 = 门；右上小星号 = keystone",
        ],
        "good_signal": (
            "某 MAG 同时在多条元素循环通路 ≥ 50% → 多功能菌。某通路被 5+ MAG 覆盖"
            "(>50%) → 该功能冗余性好，群落鲁棒。"
        ),
        "warning": (
            "全图大部分 < 50% → KO 注释严格或基因组不完整（low-quality MAG 导致）。"
            "单 MAG 某通路 100% 但别的都 0% → 可能是参考基因组高度相似，值得单独验证。"
        ),
        "caveats": (
            "completeness = 通路内**观察到的 KO** / **KB 定义的 KO 全集**。KB 定义来自"
            "KEGG module 种子 + 人工策划。完整度 ≥ 50% 是常用阈值，不是绝对标准。"
        ),
    },
    "gene_profile": {
        "title": "如何解读 MAG 元素循环基因谱",
        "what_it_shows": (
            "每个 MAG 携带**哪些** As/N/S/Fe 元素循环 KO + 拷贝数。比 pathway 更细"
            "（KO 层面）而比 gene_heatmap 更聚焦（MAG 层面）。"
        ),
        "how_to_read": [
            "行 = MAG（按门聚类）；列 = KO（按元素 As→N→S→Fe 分块）",
            "颜色 = log(1 + 拷贝数)；**0 值留白**（blank_zeros 开）清楚区分缺失 vs 低表达",
            "顶部元素色带 + 左侧门色带 + 右上 keystone 星号",
            "列宽自适应；sort_ko_by_coverage=True 时 KO 按覆盖率降序排",
        ],
        "good_signal": (
            "某 MAG 在某元素的 KO 整行高值 → 该元素代谢专家（如 Sulfuricaulis 硫氧化"
            "KO 整行深绿）。同门 MAG 呈现一致 KO 模式 → 系统发育与功能一致。"
        ),
        "warning": (
            "0 值占 > 70% 的 MAG → 基因组不完整或不参与这些元素循环。"
            "某 KO 全为 0（blank 整列）→ 可以考虑从 element_filter 暂时剔除，减少噪声。"
        ),
        "caveats": (
            "拷贝数来自 KO 注释计数，不代表表达量（要看 RNA-seq）。默认 viridis 配色"
            "适合线性呈现；0 值必须留白才能视觉区分「缺失」vs「低表达」。"
        ),
    },
    "network": {
        "title": "如何解读共现网络散点图（Gephi 辅助）",
        "what_it_shows": (
            "Degree vs Betweenness 散点图**可视化 keystone 筛选标准**。注意：完整网络图"
            "请在 Gephi 里画（此处只提供 keystone 筛选可视化 + Gephi 预处理 CSV）。"
        ),
        "how_to_read": [
            "每点 = 一个物种节点；x = Degree（连接数），y = Betweenness（中介中心性）",
            "深蓝大点 = keystone（通常 Degree ≥ 10 且/或 Betweenness ≥ 200）",
            "浅蓝小点 = 非 keystone；keystone 旁标注 Genus 物种名",
            "红色虚线 = 筛选阈值（侧边栏可调）",
        ],
        "good_signal": (
            "Keystone 数量在 10-30% 比较合理；keystone 物种跨门分布（不全是同一门）"
            "→ 网络由多个生态位的关键物种共同维持。"
        ),
        "warning": (
            "Keystone 只有 1-2 个 → 阈值过严，建议降低。Keystone > 50% → 阈值过宽，"
            "失去「关键物种」区分意义。"
        ),
        "caveats": (
            "完整网络可视化请导出 Gephi CSV（侧边栏「下载 Gephi 就绪 CSV」）→ 在 Gephi"
            "用 Fruchterman Reingold 布局（区=10000，重力=10，速度=1）。本图只做"
            "keystone 筛选辅助。"
        ),
    },
    # ═══════════════════════════════════════════════════════════
    # 循环图
    # ═══════════════════════════════════════════════════════════
    "cycle_diagram": {
        "title": "如何解读生物地球化学循环图",
        "what_it_shows": (
            "EnvMeta 差异化卖点。从 KO 注释**自动推断**元素循环活跃通路 + 跨元素耦合 +"
            "承载物种。论文 hero figure 首选。"
        ),
        "how_to_read": [
            "4 个象限 = 4 种元素（As / N / S / Fe）；每个合并细胞 = 一个 MAG 的承载通路",
            "细胞内级联：底物 → 基因椭圆 → 产物；箭头方向 = 反应方向",
            "跨元素虚线（紫/棕/绿）= 化学物耦合（如 As(III)↔H₂S → As₂S₃）",
            "底部面板 = env-pathway Spearman 相关（绿=strong, 黄=suggestive, 红=spurious?）",
            "★ = 跨组最活承载者；✦ = keystone",
        ],
        "good_signal": (
            "主要元素循环都有 MAG 承载（不空）、跨元素耦合连线出现（说明群落协同）、"
            "底部 env 面板 strong 信号与研究假说方向一致。cycle_compare 显示处理组"
            "通路贡献明显 > 对照。"
        ),
        "warning": (
            "某元素象限空 → 该元素无任何 MAG 通路 completeness ≥ 50%（阈值可调）。"
            "env 面板全是 spurious? → 样本量太少或环境梯度差。若 N<6 建议关闭 env 面板。"
        ),
        "caveats": (
            "**描述性输出，非因果断言**。图只显示「谁承载什么 + 什么与什么相关」。"
            "若需「是否支持我的假说」，配合下方的假说评分器（YAML）。相关性经过 999 次"
            "置换检验标定可信度，但不等于因果。"
        ),
    },
    "hypothesis_score": {
        "title": "如何解读机制假说评分（YAML）",
        "what_it_shows": (
            "用户上传 YAML 假说 → EnvMeta 按 5 类 claim 给数据评分。**不是假设检验**，"
            "是 MCDA 加权证据聚合。配三个稳健性指标读更可靠。"
        ),
        "how_to_read": [
            "每条 claim 一行：satisfied / partial / unsatisfied / skipped + 得分 + 权重",
            "Overall score = Σ(w·score) / Σ(w)（仅非 skipped），4 档标签：strong ≥ 0.75 / suggestive ≥ 0.40 / weak > 0 / insufficient",
            "**三指标**：null_p（999 次置换 p）/ weight_robust（±20% OAT 扰动标签是否翻转）/ veto_reasons（required=true 的 claim 若 unsatisfied 则硬否决）",
            "顶部中文一句话解读综合了 label × null_p × robust × veto 9 档判读",
        ],
        "good_signal": (
            "label = strong **且** null_p < 0.05 **且** robust = True → 最强支持"
            "（权重设计与数据贯穿一致且特异）。跨组对比表里处理组显著高于对照组。"
        ),
        "warning": (
            "label = strong 但 null_p > 0.20 → **「strong 但不特异」**：通过率幸运主导，"
            "权重设计没被数据特异支持。label 为 degenerate N/A（全 satisfied）不是坏事；"
            "但 veto_reasons 非空就是硬否决，再高 overall 也不算支持。"
        ),
        "caveats": (
            "**不是假设检验工具**（H0 无法严格定义）。是 MCDA（Multi-Criteria Decision"
            "Analysis）框架下的证据加权聚合，引用 Keeney & Raiffa 1993 / Bradford Hill 1965 /"
            "Fisher 1935。结论语气应为「证据支持等级」而非「假说被证实」。"
        ),
    },
}


def all_analysis_ids() -> set[str]:
    """返回 INTERPRETATIONS 覆盖的全部 analysis_id（测试用）。"""
    return set(INTERPRETATIONS.keys())
