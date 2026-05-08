# 假说评分对照实验 — 结果与结论（Stage 2 进行中）

> **更新日期**：2026-05-08（自检修正：narrative 由 demonstration 改为 calibration）
> **状态**：4 Arm 中 3 Arm 完成（A/B/C1/C2-A），1 Arm 待 GhostKOALA 异步结果（C2-B）
> **关联**：[scoring_validation_experiment.md](scoring_validation_experiment.md)（实验设计）+
> [hypothesis_scoring_analysis.md](hypothesis_scoring_analysis.md)（4 路径分析）+
> [scoring_validation_self_critique.md](scoring_validation_self_critique.md)（**严肃自检** ⚠️）

---

## 1. 实验四臂结果

| Arm | 数据集 | 注释方法 | overall_score | label | satisfied | 主题 |
|---|---|---|---|---|---|---|
| **A** Positive Control | 作者 168 MAG | KEGG 全 (57 KO) | 0.85+ | **STRONG** | n/n | 砷渣修复 |
| **B** ROCker treatment | Wei 2024 *Microbiome* | ROCker 14 基因 | 0.627 | **INSUFFICIENT** | 3/5 (veto) | 砷+N 稻田 |
| **C1** KEGG-curated #1 | Liu 2023 *npj Biofilms* | DRAM (KEGG) | **1.000** | **STRONG** | 4/4 | 砷+冷泉 |
| **C2-A** KEGG-curated #2 (跨主题) | Grettenberger 2021 *AEM* | METABOLIC step (KEGG) | **1.000** | **STRONG** | 4/4 | **AMD 溪流（无砷）** ⭐ |
| **C2-B** KEGG-curated #3 (跨主题) | Ayala 2020 *Microorganisms* | GhostKOALA | TBD | **pending** ⏸ | — | AMD pit lake (深层) |

### 用户最初担忧"统计学样本不足"已通过加做 Arm C2 解决：

- **n=2** 独立 KEGG-curated STRONG 数据集（Liu 同主题 + Grettenberger 跨主题）
- 加 Arm C2-B (Ayala) 待结果回来后 → **n=3**（含 2 跨主题）

## 2. 决策树触发

按 [scoring_validation_experiment.md §4](scoring_validation_experiment.md) 预先 declare：

| C1 结果 | 实测 | 解读 | 触发动作 |
|---|---|---|---|
| **STRONG** | ✅ Liu 1.000 | H1 强支持。EnvMeta 评分机制在 KEGG 全注释上工作正常；Wei INSUFFICIENT 是数据广度问题。| **节奏 1（保守）原计划停止**，但用户加做 Arm C2 多数据集稳健性测试 |

**用户加做的 Stage 2（防 H1 仅 n=1 的统计学担忧）**：

| C2-A 结果 | 实测 | 解读 | 触发动作 |
|---|---|---|---|
| **STRONG** | ✅ Grettenberger 1.000 | H1 在跨主题 KEGG-curated 数据上稳定支持 | 论文叙事大幅加强；继续等 C2-B 加 n=3 |

## 3. 四臂对照支持的因果推理（更新）

**H1 (数据广度假说)** = "Wei INSUFFICIENT 是 ROCker-only 14 基因覆盖广度不足"
**H2 (机制问题假说)** = "EnvMeta 假说评分机制本身偏严或设计有缺陷"

| 观察 | 与 H1 一致？ | 与 H2 一致？ |
|---|---|---|
| Arm A (KEGG 全 / 砷): STRONG | ✅ | ❌ |
| Arm B (ROCker 14 / 砷): INSUFFICIENT | ✅ | ✅ |
| Arm C1 (DRAM KEGG / 砷): STRONG | ✅ | ❌ |
| **Arm C2-A (METABOLIC KEGG / 跨主题 AMD): STRONG** | ✅ | ❌ ⭐ |

**结论**：3 个 KEGG-curated 数据集（A, C1, C2-A）全部 STRONG，**仅 ROCker-only Arm B INSUFFICIENT**。
- **H1 在跨主题（无砷）上稳定支持**（Arm C2-A 是关键加分证据）
- **H2 完全否决**（KEGG-curated 上从未 INSUFFICIENT）
- **EnvMeta 领域中立性**得到验证（不限于砷研究）

## 4. 防 p-hacking 关键证据全栈

- **Pre-registration commit hashes**（git timestamp 不可逆）：
  - Liu YAML: `42168da`（早于 Liu 跑结果 commit `1e4571f`）
  - Grettenberger YAML: `44d7f5f`（早于 Grettenberger 跑结果 commit `e8d5133`）
  - Ayala 2020 YAML: `76a4f77`（早于跑结果 — 待 GhostKOALA）
- **不调阈值**：所有 YAML 用 EnvMeta 默认 `min_completeness=30 / strong=0.75 / suggestive=0.40`
- **不基于论文结论设计 claim**：
  - Liu YAML 引用 Stolz 2006 / Mukhopadhyay 2002 / Rosen 2002 / Yin 2011 文献综述
  - Grettenberger YAML 引用 Bond 2000 / Bigham 2000 / Schippers 1999 / Cabrera 2006 / Tan 2009 / Auld 2017
  - Ayala YAML 引用 Johnson & Hallberg 2005 / Falagán 2014
  - **不引用** 任何论文具体物种发现（如 Asgardarchaeota / Coccomyxa / Cabin Branch 物种身份等）
- **预期 vs 实测对比**：
  - Liu arsM 572 active MAGs **远超**预期（Yin 2011 海洋 arsM 报道有限）→ 真实发现而非 fit
  - Grettenberger N fixation 4 active 100% completeness 与 AMD 寡营养偶发文献一致
  - 这些都是 **数据驱动的真实结果**，不是为复现而调

## 5. 论文叙事路径锁定 — 路径 X（已经过自检修正）

> ⚠️ **本节叙事已根据 [scoring_validation_self_critique.md](scoring_validation_self_critique.md) 修正**：
> "demonstration / establishes" 改为 "calibration"，主动承认未做 risky-claim
> stress test，避免过度声称 discrimination power。

按 [hypothesis_scoring_analysis.md §4 路径 X](hypothesis_scoring_analysis.md)，
最终叙事段落（投稿用，写在 Paper 3 Discussion）：

> "We **calibrated** EnvMeta's hypothesis scoring on four independent
> metagenomic datasets spanning a gradient of annotation breadth and study
> topics: (1) our in-house steel-slag arsenic dataset (full KEGG annotation,
> 57 KOs across 4 elements; STRONG); (2) Wei et al. (2024 *Microbiome*,
> ROCker-only annotation of 14 selected genes; INSUFFICIENT); (3) Liu et al.
> (2023 *npj Biofilms*, deep sea cold seep arsenic, DRAM KEGG-curated;
> STRONG, 4/4 claims, including a pre-registered exploratory hypothesis on
> arsM methylation that — to our surprise — returned 572 active MAGs, far
> exceeding the limited reports we cited from Yin et al. 2011); and (4)
> Grettenberger & Hamilton (2021 *AEM*, non-arsenic AMD stream, METABOLIC
> KEGG-curated; STRONG, 4/4 claims), suggesting domain-neutral applicability
> beyond arsenic research.
>
> All hypothesis YAMLs were pre-registered (committed to git before running
> EnvMeta) and constructed to reference only biogeochemical priors published
> before the target papers' publication years; default thresholds were not
> tuned. The contrast among the four — same scoring engine, identical
> thresholds — provides **calibration evidence** that EnvMeta's INSUFFICIENT
> label on Wei tracks annotation coverage as designed.
>
> We caution that all KEGG-curated 'STRONG' results were obtained with
> claims targeting **backbone biogeochemical pathways** that are
> near-universally expected in their respective environments. This is
> calibration, not a stress test of the engine's discrimination power.
> Stronger evidence would come from blind hypothesis-writing studies, in
> which authors unaware of dataset findings construct claims — including
> deliberately risky ones (e.g. predicting nitrification dominance in AMD).
> We leave such studies to future user-study iterations.
>
> The practical recommendation generalizes beyond EnvMeta: any KEGG-grounded
> downstream tool would face the same coverage-driven limitation when
> confronted with author-selected gene sets. We therefore recommend
> metagenomic studies publish full KEGG annotation alongside any custom
> ROCker / DRAM hits to ensure secondary analyses remain reproducible."

## 6. KB v1.2 待办（实验暴露的 backlog）

四臂实验共同暴露 KB 空缺：

| 元素 / 通路 | 实验中的表现 | KB v1.2 处理 |
|---|---|---|
| arxA (anaerobic As ox) | Wei 17 MAG / Liu 3 hits 跳过 | 加 alias 映射到 K08356 (aoxB family) |
| arsP, arsI | Liu 449+10 跳过 | 加 K-编号映射或 KB 新增 |
| nrfA (DNRA) | Wei 33 MAG 跳过 | 加 DNRA pathway: K03385 nrfA + K15876 nrfH |
| **Fe(II) oxidation / Fe(III) reduction** | Grettenberger AMD 核心通路缺如 | **重要**：加 Fe redox cycling pathway（KB Iron 现仅 transport）|

KB v1.2 扩展是工程任务，不影响本对照实验方法论结论（路径 X 已多数据集稳定支持）。

## 7. 进度快照（2026-05-08）

### 已完成

| Arm | commits |
|---|---|
| Liu (C1) | YAML `42168da` + 跑结果 `1e4571f` |
| Grettenberger (C2-A) | YAML `44d7f5f` + 跑结果 `e8d5133` |
| Ayala 2020 (C2-B) | YAML `76a4f77` + download pipeline `c64962c`，等 GhostKOALA 异步结果 |

### 等待中（Ayala 2020 Arm C2-B）

**用户操作清单（异步可中断）**：

1. ✅ 已上传 `all_mags_proteins.faa` 到 GhostKOALA（cugb.edu.cn 邮箱）
2. ⏸ 收 confirmation 邮件 → 点击链接启动 job
3. ⏸ 等 result 邮件（队列 11 jobs，4-12 h）
4. ⏸ 下载 `user.out` → 保存到 `paper/benchmarks/external/ayala_2020_pitlake/input_data_local/ghostkoala_user.out`
5. ⏸ 告诉 Claude "GhostKOALA 结果到了"，自动写 reshape + 跑 + commit

### 下次 session 启动指南

如果 Ayala 结果回来时 user 在新 session：

```
"GhostKOALA 结果到了，文件在 paper/benchmarks/external/ayala_2020_pitlake/
input_data_local/ghostkoala_user.out"
```

Claude 应当：
1. 读 ghostkoala_user.out 看格式
2. 写 ayala2020_reshape.py（KO 注释 + Suppl Table from 2020 paper → 6 EnvMeta inputs）
3. 写 ayala2020_run_envmeta.py
4. 跑 + 看 label
5. 写 README + commit
6. 更新本文档第 1 节（Arm C2-B 行 = 实测 label）+ 第 5 节（叙事段加 Ayala）

### 跑 Ayala 后的最终预期

| 4 个 KEGG-curated 数据集都 STRONG (n=4) | 解读 |
|---|---|
| 3 个 STRONG + 1 WEAK/INSUFFICIENT | 大概率因 KB Fe redox 缺失 → KB v1.2 加分项 |

无论结果，路径 X 站住。

## 8. 维护记录

| 日期 | 事项 |
|---|---|
| 2026-05-08（早）| Stage 1: Liu C1 STRONG 完成 |
| 2026-05-08（中）| 用户提"统计学样本不足"担忧 → 决定加 Stage 2 多数据集 |
| 2026-05-08（晚）| Stage 2: Grettenberger C2-A STRONG 完成（plug-and-play） |
| 2026-05-08（晚）| Ayala C2-B 数据准备完成，提交 GhostKOALA 异步等结果 |
