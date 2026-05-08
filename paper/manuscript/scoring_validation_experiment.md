# 假说评分对照实验 — 跨砷数据集复现性测试

> **创建日期**：2026-05-08
> **关联**：[hypothesis_scoring_analysis.md](hypothesis_scoring_analysis.md)（4 条叙事路径分析）
> **设计目的**：用对照实验排除"工具机制本身有问题"vs"Wei 数据广度限制"两种解释，
> 从而决定后续叙事路径（X / W / 重新审视机制）。
> **方法论纯净性**：所有决策**预先 declare**，避免 p-hacking。阈值不调，YAML 不为
> 复现而调，仅看 EnvMeta 默认设定在不同数据集上的输出分布。

---

## 1. 研究问题

```
EnvMeta 假说评分在 Wei 2024 数据上输出 INSUFFICIENT 的根因是什么？

H1 (数据广度假说)：Wei 用 ROCker-only 14 基因 ≠ KEGG 全注释，注释广度不足
                    导致 pathway-completeness 不达阈值。如果在 KEGG 全注释
                    数据集上跑，应输出 STRONG/SUGGESTIVE。

H2 (机制问题假说)：EnvMeta 假说评分机制本身偏严或设计有缺陷，即使在 KEGG
                    全注释数据集上仍倾向输出 INSUFFICIENT。

H3 (混合)：两者都有部分贡献，需细化。
```

## 2. 实验设计 — 三臂对照（adaptive sequential）

### Arm A — Positive Control（已有，无需重跑）

- 数据：作者自己 168 MAG × 10 sample × 57 KEGG KO（KEGG 全注释 + per-sample abundance）
- YAML：`paper/hypotheses/arsenic_steel_slag.yaml`
- 已知结果：**STRONG** label（从 v0.7.0 release 开始就稳定）
- 作用：证明 EnvMeta 评分机制在"高质量 KEGG 全注释"输入下能输出 STRONG，
  排除"机制永远输出 INSUFFICIENT"的极端情况

### Arm B — Treatment 1: ROCker-only（已有，无需重跑）

- 数据：Wei et al. 2024 *Microbiome* (PRJNA1068274), 14 ROCker 基因
- YAML：`paper/benchmarks/external/wei_2024_paddy/wei2024_hypothesis.yaml`
- 已得结果：**INSUFFICIENT** (overall=0.63, null_p=0.9, weight_robust=True)
- 作用：被检验对象 — 是否因 ROCker-only annotation 导致

### Arm C — Treatment 2: 独立第三方 KEGG 全注释（**待跑**）

候选（按推荐度排序，详见 §3）：

| Slot | 候选 | 注释方法 | 主题 |
|---|---|---|---|
| C1 | Liu et al. 2023 *npj Biofilms* | DRAM v1.3.5（含 KEGG）| 深海冷泉砷循环 |
| C2 | Yang et al. 2024 *ISME J* | KofamScan + KEGGdecoder | 矿区贫营养砷酸盐还原 |
| C3 (备选) | Wang et al. 2025 *Soil Ecol Lett* | KofamScan | As/Sb 共污染土 Chloroflexota |

**Sequential testing 策略**：先做 C1，按结果决定是否做 C2。

### 实验矩阵

| Arm | 数据集 | 注释方法 | 预期 label (H1 假说) | 预期 label (H2 假说) | 实测 |
|---|---|---|---|---|---|
| A | 作者 168 MAG | KEGG 全 | STRONG | STRONG | STRONG ✅ |
| B | Wei 2024 | ROCker 14 基因 | INSUFFICIENT | INSUFFICIENT | INSUFFICIENT ✅ |
| C1 | Liu 2023 | DRAM (KEGG) | STRONG/SUGGESTIVE | INSUFFICIENT | ⬜ 待跑 |
| C2 | Yang 2024 | KofamScan | STRONG/SUGGESTIVE | INSUFFICIENT | ⬜ 备选 |

## 3. 候选数据集详情

### Slot C1（首选）：Liu et al. 2023 *npj Biofilms & Microbiomes*

| 项 | 内容 |
|---|---|
| 论文 | "Unexpected genetic and microbial diversity for arsenic cycling in deep sea cold seep sediments" |
| DOI | [10.1038/s41522-023-00382-8](https://doi.org/10.1038/s41522-023-00382-8) |
| 数据 | NCBI BioProject `PRJNA831433`（22 arsenotrophic MAGs）+ 论文 Suppl |
| 注释方法 | **DRAM v1.3.5**（KEGG + Pfam + MEROPS + dbCAN）+ METABOLIC pipeline + GTDB-Tk v2.1.1 r207 |
| 规模 | 87 metagenomes + 33 metatranscriptomes + 22 完整发布的 arsenotrophic MAGs |
| 核心假说 | Asgardarchaeota + 4484-113 + AABM5-125-24 + RBG-13-66-14 等多门微生物在贫氧深海冷泉驱动砷转化 |
| reshape 工时估计 | **6-10 小时**（DRAM `metabolism_summary.xlsx` 直出 KO matrix）|
| 主题 | 深海冷泉砷循环（与 Wei 稻田 + 作者矿渣完全异生境，跨主题 + 同元素双重对照）|

**为什么首选 C1**：
- DRAM 输出格式工整，reshape 最快（< 1 天）
- 87 metagenomes 规模充足
- 22 MAGs 完整发布，无需重跑流程
- 主题异于 Wei 和作者数据 → 跨生境泛化性证据

### Slot C2（次选）：Yang et al. 2024 *ISME J*

| 项 | 内容 |
|---|---|
| 论文 | "Microbially mediated sulfur oxidation coupled with arsenate reduction within oligotrophic mining-impacted habitats" |
| DOI | [10.1093/ismejo/wrae110](https://doi.org/10.1093/ismejo/wrae110) |
| 数据 | NCBI BioProject `PRJNA989741` + Figshare `10.6084/m9.figshare.25189016` |
| 注释方法 | **KofamScan + KEGGdecoder** （KEGG Release 97.0 全注释）|
| 核心假说 | 硫氧化菌（chemolithotrophic SOAsR）在贫营养矿区驱动砷酸盐还原 + 砷释放 |
| reshape 工时估计 | **8-12 小时** |
| 主题 | 矿区贫营养，与作者矿渣场最接近 |

### Slot C3（备选）：Wang et al. 2025 *Soil Ecology Letters*

| 项 | 内容 |
|---|---|
| 论文 | "The hidden diversity and functional potential of *Chloroflexota* genomes in arsenic and antimony co-contaminated soils" |
| DOI | [10.1007/s42832-024-0266-y](https://doi.org/10.1007/s42832-024-0266-y) |
| 数据 | Figshare `10.6084/m9.figshare.26207045.v1`（170 Chloroflexota MAGs）|
| 注释方法 | KofamScan / KofamKOALA（Aramaki 2020）|
| 风险 | Figshare 是否含 KO matrix 待打开确认 — 若仅含 MAG fasta，需重跑 KofamScan（本机 4-8h）|
| reshape 工时估计 | **6-8h（如有 KO 表）/ 12-20h（如需重跑）** |
| 主题 | As/Sb 共污染土 Chloroflexota |

## 4. 预先 declare 的决策树（**禁止事后修改**）

### 第一步：跑 Arm C1 (Liu 2023)

| C1 结果 | 解读 | 下一步动作 |
|---|---|---|
| **STRONG** | H1 强支持。EnvMeta 评分机制在 KEGG 全注释上工作正常；Wei INSUFFICIENT 是数据广度问题。| **停止**。叙事路径锁定 X（重新框定）。可选 W（KB v1.2 加 DNRA + reduced_coverage status）。 |
| **SUGGESTIVE** | H1 中度支持。机制合理，可能边界 case。| **停止**。叙事路径 X，同时论文 Discussion 加一段"评分阈值的边界讨论"。 |
| **WEAK** | H3 混合。机制可能偏严，但不是失效。| **跑 C2** (Yang 2024) 验证是否系统性。|
| **INSUFFICIENT** | H2 警告信号。Liu 数据集 KEGG 全注释仍 INSUFFICIENT。| **必须跑 C2**。如果 C2 也 INSUFFICIENT → 评分机制需重新审视。|

### 第二步（仅在 C1 = WEAK / INSUFFICIENT 时）：跑 Arm C2 (Yang 2024)

| C1 + C2 联合结果 | 解读 | 后续动作 |
|---|---|---|
| C1=INSUFFICIENT + C2=STRONG | C1 是 case-specific（深海冷泉与作者 KB 砷代谢 reaction 略异）| 叙事路径 X 站住，但论文里写 "C1 case study" 段说明特定数据集挑战 |
| C1=INSUFFICIENT + C2=INSUFFICIENT | **评分机制有系统性问题**。可能是 required veto 太硬 / pathway_active 阈值定义不当 / 化学物 coupling 表不全。| **不再走路径 X**。重写假说评分机制 paper section 4.5 — 这段从"工具优势"变成"工具的当前限制 + 未来工作"。可能影响投稿期刊（iMeta 可能仍接受，但叙事大改）|
| C1=WEAK + C2=STRONG | 边界 case，C2 验证机制可工作 | 叙事路径 X，强调 "EnvMeta 评分对数据集质量敏感"，不是缺陷 |
| C1=WEAK + C2=WEAK | 评分机制偏严但不失效 | 路径 W（X + 工具加 reduced_coverage status，可能松动 default 阈值的 motivation 由数据驱动产生）|

### 第三步（仅极端情况）：跑 Arm C3 (Wang 2025) 或扩展

只在 C2 仍无法消歧时考虑。预计不需要走到这一步。

## 5. 投入估算

| 路径 | 数据集 | 工时 | 投入累计 |
|---|---|---|---|
| 最佳（C1=STRONG）| 仅 C1 | 6-10 h | 6-10 h |
| 中等（C1=WEAK/INSUFFICIENT，C2=STRONG）| C1 + C2 | 6-10 + 8-12 = 14-22 h | ≤ 1 周（兼职） |
| 最差（C1+C2 都 INSUFFICIENT）| C1 + C2 + 重审机制 | 14-22 + 5-10 (重写 §4.5) = 19-32 h | ≤ 1.5 周 |

**最大投入 1.5 周**，远低于"瞎调阈值 + 投稿被拒重写" 的潜在成本。

## 6. 防止 p-hacking 的关键约束

⚠️ **以下规则锁定，不可事后修改**：

1. **YAML 模板**：每个数据集用一份独立 YAML，但 claim 结构、阈值、weight、required
   完全沿袭 [arsenic_steel_slag.yaml](../hypotheses/arsenic_steel_slag.yaml) 的默认值。
   仅替换 pathway 名 / coupling species 以适应该数据集的元素覆盖。
2. **不调 min_completeness**：保持 30 默认。
3. **不调 strong/suggestive thresholds**：保持 0.75 / 0.40 默认。
4. **不调 required veto 配置**：每个数据集的 required claim 由作者论文核心论断决定，
   不由 EnvMeta 输出决定。
5. **每个数据集 YAML 必须在跑 EnvMeta **之前**写完并提交 git**（git timestamp 作为
   pre-registration 证据）。
6. **跑出 INSUFFICIENT 后不允许回去改 YAML 让结果变好**。如发现 YAML 写错（claim
   语义不符合作者原文），可改但需 commit 单独标注"修正"，且只允许"原文核实修正"，
   不允许"为复现而调"。

## 7. 最小工作量启动方案

**今天/明天的具体动作**（用户决定后启动）：

```bash
# 1. 用户下载 Liu 2023 Suppl + DRAM 输出
#    https://www.nature.com/articles/s41522-023-00382-8#Sec19
#    保存到 D:\download\liu_2023_coldseep\
#
# 2. 我来：
#    - 写 tools/external_benchmarks/liu2023_reshape.py
#    - 写 paper/benchmarks/external/liu_2023_coldseep/liu2023_hypothesis.yaml（在跑 EnvMeta 之前 commit）
#    - 跑 EnvMeta + 假说评分
#    - 按决策树看 label，决定是否跑 C2
```

预计 24-36 小时内出 C1 结果。

## 8. 输出标准

```
paper/benchmarks/external/liu_2023_coldseep/
├── README.md                  # 数据集 + 复现路径 + EnvMeta 输出 + label
├── input_data_local/          # reshape 输出（gitignore）
├── envmeta_outputs/           # 5 PDF + 1 hypothesis md + stats tsv
├── liu2023_hypothesis.yaml    # 在跑 EnvMeta 之前 commit
└── compare_to_original.md     # 与原文对比 + 与 Wei C1 对比
```

## 9. 维护记录

| 日期 | 事项 |
|---|---|
| 2026-05-08 | 实验设计存档；Arm A + Arm B 已完成；候选 C1/C2/C3 调研完毕 |
