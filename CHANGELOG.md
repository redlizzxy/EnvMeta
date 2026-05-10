# 更新日志

本项目遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/) 格式，
版本号采用 [Semantic Versioning](https://semver.org/lang/zh-CN/) 约定（主.次.补丁）。

> 新手用户：在 README 看到提示有新版本时，本地执行 `git pull origin master` 即可更新。
> 完整流程见 [README "更新到新版本"](README.md#-更新到新版本内测阶段) 一节。

---

## [0.9.1] — 2026-05-09 (later)

> **dominance_score 字段兑现 + Paper 3 写作素材完整**：v0.9 Discussion 中
> 标的 future work 已落地工程；Paper 3 投稿核心证据（Methods + Results +
> Discussion + Tables + Figure X）全部就位。

### ✨ 新增

- **`dominance_score` 字段** + 可选 `min_dominance_fraction` 硬阈值
  （`pathway_active` claim type）：
  ```
  dominance_score = pathway.total_contribution /
                    sum(all pathways in element)
  ```
  evidence 字典总是含 `dominance_score` + `element_total_contribution`
  （`pathway_active` + `pathway_inactive` 两 evaluator 都加，信息透明）。
  解决 v0.9 stress test 暴露的二元 `mean_completeness ≥ 50%` 阈值 limit。
  向后兼容（不指定 min_dominance_fraction 即不加阈值）。
  ([fdfae77](https://github.com/redlizzxy/EnvMeta/commit/fdfae77))

- **Stress test v2 YAMLs**（Liu + Ayala，commit `fdfae77`）：Class A reversed
  claim 加 `min_dominance_fraction: 0.20`，其他 claims 不变。v1（commit `50c4687`）
  保留作为 B-tier 历史证据，v2 作为 A-tier 升级证据，git history 完整 audit trail。
  - **B → A 级 discrimination 升级**：Liu 0.625 → **0.250 (weak)** ⭐；Ayala
    0.455 → **0.182 (weak)** ⭐。Class A claim dominance 实测 0.05% (Liu) / 7.08%
    (Ayala) << 20% 阈值 → 正确 `unsatisfied`。**3/3 stress test 现全 A 级 clean
    discrimination**。

- **Paper 3 投稿核心写作素材**：
  - [`paper/manuscript/methods_external_validation.md`](paper/manuscript/methods_external_validation.md)
    Methods §4.6 完整版（~1450 字，7 子节，19 条 Vancouver 引用 + DOI）
  - [`paper/manuscript/results_stress_test_section.md`](paper/manuscript/results_stress_test_section.md)
    Results §X stress test 章节（~800 字，4 子节）
  - [`paper/manuscript/discussion_calibration_discrimination.md`](paper/manuscript/discussion_calibration_discrimination.md)
    Discussion §Y limitation + future work（~640 字，4 子节）
  - [`paper/figures/paper3_hypothesis_scoring/`](paper/figures/paper3_hypothesis_scoring/)
    Table 1 + Table 2 + Figure X（PDF + PNG 600dpi + SVG + matplotlib 复现脚本）

- **AMD diazotrophy 引用 verify 完成**：用户手动 verify Sánchez-España 2008
  DOI ✅；Agent verify Korehi 2014 (主题不对，与 Auld 同类错引模式) + 真正
  AMD diazotrophy 文献 = **Dai 2014 PLoS One** (10.1371/journal.pone.0087976,
  primary metagenomic 742 nif sequences) + **Méndez-García 2015 Front Microbiol**
  (10.3389/fmicb.2015.00475, review). 2 calibration YAML 引用 metadata 修订
  ([cae2de7](https://github.com/redlizzxy/EnvMeta/commit/cae2de7))

### 🧪 测试

- pytest **301/297 全绿**（+4 dominance_score 测试，无回归）

### 📝 paper-side maintenance (2026-05-09 后续 session)

> 本节记录 v0.9.1 后续的 paper-side 修订（无代码 API 变更，仅 metadata + 文档）。

- **KB metadata 字段对齐**：`elements.json` `"version"` 字段从 `"1.1"` → `"2.0"`
  匹配实际 KEGG-driven 4 元素 × 18 通路 × 57 KO 的当前内容（之前因历史原因停留在 v1.1 字符串）
  ([commit `b9c5699`](https://github.com/redlizzxy/EnvMeta/commit/b9c5699))
- **Paper 3 outline 完整整合**（commit `e1da1d9`）：3 段 Paper 3 草稿
  （Methods §4.6 / Results §X / Discussion §Y）+ performance 段落整合到 outline_imeta.md
- **3-tier 性能 benchmark + scaling figure**（commit `e1c69a6`）：sample real-dense /
  Liu real-sparse / Liu synthetic-dense 共 58 measured cells + scaling_curve.{pdf,svg,png}
  + paper-facing performance.md + 中文用户面 docs/performance_zh.md
- **ImageGP 2 ecosystem-extension reframing**（commit `3c9662b`）：发现 ImageGP 2
  (Chen et al. 2024 *iMeta* 3:e239) 是 Liu YX (EIC) 团队作品，全文从 "vs 竞品"改为
  "complement to iMeta visualization ecosystem"；致敬 14 个 sister tools；Abstract /
  Highlights / §5.1 Intro / §5.2.9 / §5.6 Acknowledgments 全部重写
- **Mock review v0.9.1 → v0.9.2 修订**（commit `b9c5699`）：5 大 Major + 5 个 Minor
  温和解修订（4 Resolved + 1 honestly acknowledged）；其中含本节 KB metadata 对齐 +
  As₂S₃ chemistry refs (Newman 1998 / Rodriguez-Freire 2014 / Hollibaugh 2005) +
  999-perm Anderson 2001 justification + Stolz 2006 编号 #22；Recommendation: Major
  Revision → **Minor Revision**
- **Publication strategy 存档**：[`paper/manuscript/publication_strategy.md`](paper/manuscript/publication_strategy.md)
  记录 Path A 改良版（EnvMeta + bioRxiv preprint + 课题论文并行起草）满足 2026 接收
  + 2027-05 见刊硬约束
- **1-day perturbation analysis 落地**：mock review v0.9.2 Major #1 应对（auxiliary
  evidence on author selection bias）。新增 [`tools/external_benchmarks/perturbation_analysis.py`](tools/external_benchmarks/perturbation_analysis.py)
  双模式 runner：within-element + cross-element，N=20 / mode / dataset，3 个
  external calibrations（Liu / Grettenberger / Ayala）共 123 runs。Headline：
  Liu 2023 cross-element 0/20 STRONG（median 0.000）→ 元素级 target accuracy
  机制刚需；within-element 40-50% 仍 STRONG（与 §Y.1 KEGG-coverage-dependent 框架一致）。
  写作素材：[`paper/manuscript/perturbation_analysis_results.md`](paper/manuscript/perturbation_analysis_results.md)；
  集成 Methods §4.6.7 + Results §X.3 + Discussion §Y.3 limitation #1；
  figure: `paper/benchmarks/external/perturbation/perturbation_curve.{pdf,png,svg}`。
  pytest 301/301 仍全绿（无 API 变更）。
- **Mock review v0.9.5 (handling editor) 4 cheap Major 修订**（2026-05-10 third session）：
  - **OpenTimestamps 区块链锚点已落地**（v0.9.5 Major #3）：新增
    [`paper/manuscript/timestamps/`](paper/manuscript/timestamps/) 目录含
    `stamp_via_https.py` runner + 4 commit hash files + 12 .ots-response.bin
    proof files (alice + bob + finney calendars × 4 commits) + `ANCHOR_SUMMARY.md`；
    §4.6.2 现引用 OpenTimestamps anchoring 而非 "future commit"
  - **null_p 重 frame**（v0.9.5 Major #1）：§4.6.1 明确标为 "shuffle-consistency
    diagnostic, explicitly NOT a frequentist p-value"；解释 4-9 claim 离散 satisfaction
    null distribution 必然粗糙；加 WoE 软件 relationship 段
  - **Arm A "saturation" → "partial-perturbation robustness check"** （v0.9.5 Major #4）：
    Results §X.3 + perturbation_results 4 处全部 reframe；明确 6 unperturbed claims
    占总权重 ≈72%，arithmetic 锚定 score
  - **Distance-to-boundary**（v0.9.5 Minor #4）：Results §X.2 + Table 2 描述加
    stress score 到最近 label boundary 的距离（Grett −0.150 / Liu −0.125+0.225 /
    Ayala +0.055 near-boundary 与 §4.6.8 阈值 transition 一致）
  - 测试 **301/301 全绿**
- **Mock review v0.9.4 (independent reviewer) 5 Major 修订**（2026-05-10 second session）：
  - **Major #1 mock-review metadata 泄漏**：删 manuscript 4 文件（methods §4.6.7 +
    perturbation_analysis_results §1/§4 + outline §5.7.1 ref annotations + §5.3 author notes）
    全部 "Mock Review v0.9.x" 标签
  - **Major #2 Arm A baseline 0.919 vs Table 1 = 1.000 数值不一致**：
    扩展 perturbation runner 支持 keystone_df，sample_data 跑出 1.000 STRONG 9/9
    与 Table 1 一致；Results §X.1 "4/4 satisfied" 改为 "Arm A 9/9, externals 4/4"
  - **Major #3 pre-registration institutional-trust-based**：Methods §4.6.2 加段
    诚实承认 author-controlled repo + 3 mitigations + future OpenTimestamps commit
  - **Major #4 默认阈值 0.75/0.40 缺 sensitivity**：新增
    [`tools/external_benchmarks/threshold_sensitivity.py`](tools/external_benchmarks/threshold_sensitivity.py)
    扫 5 阈值 × 8 dataset；新增 Methods §4.6.8 段；4 KEGG-curated arm 全 STRONG
    跨 0.65-0.85 全稳；Arm B 全 INSUFFICIENT 全稳；3 stress 全 weak/suggestive；
    输出 [`paper/benchmarks/external/threshold_sensitivity/`](paper/benchmarks/external/threshold_sensitivity/)
  - **Major #5 Arm A 循环论证 reframe**：Results §X.1 完全重写——明确把 Arm A 标为
    "in-house positive control (engine self-consistency check, not independent
    calibration evidence)"，3 external arm 标为 "primary calibration evidence"；
    Discussion §Y.3 limitation #1 同步 reframe
  - 测试 **301/301 全绿**；预期 mock review v0.9.5 全部 Major Closed
- **Mock review v0.9.3 Major + Minor 全部修订**（2026-05-10 session）：
  - **v0.9.3 Major #1 (Arm A asymmetry) → Closed**：扩展 perturbation runner 支持 Arm A partial
    perturbation（仅 3 个 pathway_active claim，sample_data 数据加载）。Arm A 100% STRONG
    retention 揭示 saturation regime；跨数据集形成 monotonic gradient
    (Arm A 100% → Grettenberger 30% → Ayala 15% → Liu 0%)，Wilson 95% CI 标注
  - **v0.9.2 Major #2 (3-dataset stress caveat) → Closed**：§X.2 加 "2-of-2 within
    3-dataset stress design" caveat
  - **3 个 v0.9.3 新 Minor → Closed**：§4.6.7 加 N=20 rationale + Wilson CI + run_null=False 辩护
  - **5 个 v0.9.2 carried Minor → Closed**：calibration evidence disambiguation (§4.6.4) /
    Newman 1997→1998 reconcile (§5.7.1 #24 改 Newman 1997 *AEM* 63:2022) /
    DRAM-KEGG 区别 (§4.6.3) / dense annotation (§5.2.8) 加 "dense relative to KB target pathways"
    clarification / consistent with 重复软化 (§X.2 末段)
  - **§5.7.1 numbering**：22-entry flat list（mock review v0.9.2 Minor #3）
  - 测试 **301/301 全绿**；预期 mock review v0.9.4 全部 Major Resolved → bioRxiv ready

---

## [0.9.0] — 2026-05-09

> **假说评分对照实验完成 + Stress test discrimination evidence**：4 Arm KEGG-curated
> 数据集 calibration 全 STRONG + n=3 stress test 暴露 discrimination power。
> Paper 3 投稿核心证据全部就位。

### ✨ 新增

- **`pathway_inactive` claim type**（第 6 类）：Popperian falsifiability 主力。
  评估：n_active_mags == 0 → satisfied；> 0 且 mean_comp ≥ max_completeness →
  unsatisfied。default `max_completeness=50`；不破坏旧 YAML（向后兼容）。
  4 个新测试 case，全套 297/293 全绿（[14bc01b](https://github.com/redlizzxy/EnvMeta/commit/14bc01b)）

- **假说写作教程** 双层结构：
  - [`docs/hypothesis_writing_guide.md`](docs/hypothesis_writing_guide.md) — 用户教程
    （含双层 calibration+stress 模板、pre-registration 纪律、pre-prediction 模板、
    6 类 claim 选择指南、Bradford-Hill 对应）
  - [`paper/hypotheses/HYPOTHESIS_DESIGN_PRINCIPLES.md`](paper/hypotheses/HYPOTHESIS_DESIGN_PRINCIPLES.md) — 论文方法学
    （含设计哲学、4 类局限、引用建议）

- **通用 stress runner**：[`tools/external_benchmarks/run_stress_yaml.py`](tools/external_benchmarks/run_stress_yaml.py)
  接受 `--dataset` + `--yaml` 或 `--all`，跑任一 external dataset 的任一 YAML
  → 输出 `fig6_<yaml_stem>_score.md`

### 📊 假说评分对照实验完成（Paper 3 核心证据）

**Arm C2-B 完成（Ayala 2020 pit lake）**：GhostKOALA 1h 完成 → reshape user_ko.txt
（11243 MAG-KO records × 13 MAGs）→ Calibration STRONG (1.000) + Stress SUGGESTIVE (0.455)。
n=4 KEGG-curated calibration 全 STRONG。([60a5be4](https://github.com/redlizzxy/EnvMeta/commit/60a5be4))

**Stress test n=3（Liu / Grettenberger / Ayala）**：

| Dataset | Calibration | Stress | Discrimination |
|---|---|---|---|
| Liu 2023 (冷泉同主题) | STRONG (1.000) | suggestive (0.625) | B 级 |
| Grettenberger 2021 (AMD 跨主题) | STRONG (1.000) | **weak (0.250)** | **A 级** ⭐ |
| Ayala 2020 (pit lake 跨主题) | STRONG (1.000) | suggestive (0.455) | B 级 |

**最强单条证据**：cross-topic `arsenate_reduction_should_dominate` 在 **2/2 无砷
数据集**（Grettenberger + Ayala）都正确 unsatisfied (n=0 active MAGs) → 双重
推翻 universal arsC 担忧 → EnvMeta 领域中立性铁证。
([50c4687](https://github.com/redlizzxy/EnvMeta/commit/50c4687) +
[ba2055c](https://github.com/redlizzxy/EnvMeta/commit/ba2055c) +
[60a5be4](https://github.com/redlizzxy/EnvMeta/commit/60a5be4))

### 📚 文档

- **stress test 双层文档**：
  - [`stress_test_predictions.md`](paper/manuscript/stress_test_predictions.md)（冻结盲预测）
  - [`stress_test_results.md`](paper/manuscript/stress_test_results.md)（实测对照 + 论文叙事段落）

- **假说评分对照实验主结果**：[`scoring_validation_experiment_results.md`](paper/manuscript/scoring_validation_experiment_results.md)
  §10 stress test 章节加 n=3 实测；§5 narrative 由 demonstration 改为 calibration（[1b358fb](https://github.com/redlizzxy/EnvMeta/commit/1b358fb)）

- **严肃自检**：[`scoring_validation_self_critique.md`](paper/manuscript/scoring_validation_self_critique.md)
  Direct/Inferred/Weak 提炼度评估；calibration claim 设计代表性 C 级 自批

- **引用审计**：[`hypothesis_references_audit.md`](paper/manuscript/hypothesis_references_audit.md)
  16 claim × 13 文献 DOI 表；Agent verify 暴露 4 处引用错误（详见 fix(refs)）

### 🐛 修复（hypothesis YAML metadata only）

- **6 个 hypothesis YAML 引用 metadata 修订**（不动 claim 实体）：
  4 处错引：(1) Yin 2011 期刊 ES&T → Plant Physiol；(2) Bothe "2007 FEMS Rev"
  实际不存在 → Bothe et al. 2000 FEMS Microbiol Rev 24:673-690；(3) Cabrera 2006
  期刊 Process Biochem → J Hazard Mater，且主题是金属毒性 SRB 实验非
  acidophilic SRB 综述 → Sánchez-Andrea 2014 升 primary；(4) Auld 2017 主题是
  seasonal community variation 非 diazotrophy → 该 claim N fixation 文献支持
  降级 weak（Methods 应改引 Korehi 2014 / Mendez-Garcia 2015 verified 文献）。
  每 YAML 末尾追加 § REFERENCES 完整 Vancouver + DOI 表。
  ([ddd3098](https://github.com/redlizzxy/EnvMeta/commit/ddd3098))

### 🧪 测试

- pytest **297/297 全绿**（+4 pathway_inactive 测试，未破坏现有功能）

---

## [0.8.2] — 2026-05-08

> R/Python 11 图侧侧对照工作完成 + RDA 数值对齐修复。

### 🐛 修复

- **RDA 数值与 R vegan 不一致**：`skbio.stats.ordination.rda()` 与
  `vegan::rda()` 的 eigval 归一化方式不同，导致 inertia 数值差 16-20×、
  环境因子 ANOVA F 差 3.5×、p 值方向反转（R: pH p=0.002 ↔ EnvMeta: 0.216）。
  修法：保留 `skbio.rda()` 输出的 ordination axes（site / biplot scores），
  但 inertia + ANOVA F 改用 SS-based 公式（`SS(Y) / (n-1)`，对标 vegan）。
  `explained_ref` 默认从 `"constrained"` 改为 `"total"`（对标 R vegan
  `summary()$cont$importance[2,]` 默认）。修复后所有 F / r / 解释度数值
  与 R vegan 4 位精度对齐。
  （[863e886](https://github.com/redlizzxy/EnvMeta/commit/863e886)）

### 📚 文档 / 验证

- **R/Python 侧侧对照 11 图全部完成**（`paper/benchmarks/validation/`）：
  - 数值精确一致（5 图）：alpha / pcoa / mag_quality / stackplot / RDA（修复后）
  - 算法等价（6 图）：lefse / gene_heatmap / log2fc / gene_profile / pathway / mag_heatmap
  - 5 个 R 脚本副本（强制 sans 字体）+ 11 个 README + 论文引用模板
  - 撤回 4 个误判"伪 bug"标记（log2fc 方向 / mag_heatmap 标尺 / gene_profile KO 过滤 / pathway KB 升级，均确认非 bug）
  （[ebe88b5](https://github.com/redlizzxy/EnvMeta/commit/ebe88b5)）

- **Paper 3 投稿前清单**（`paper/manuscript/paper3_pre_submission_checklist.md`）：
  跟踪 4 大任务进度（R 对照 ✅ / 第二数据集 ⬜ / 英文 README ⬜ / Methods 起草 ⬜）

### 🧪 测试

- pytest **293/293 全绿**（修复未破坏现有功能）；test_rda 4 个测试全部通过

---

## [0.8.1] — 2026-04-21

> Mac 端首批内测反馈修复版。聚焦 macOS 安装 / 假说评分 / 交互 HTML 三处实测 bug。

### 🐛 修复

- **HTML 交互导出**：切换分组（A / B / CK）后化学物-通路连线失效（拖拽不跟随、
  hover 不高亮）。根因为 `renderCycle` 清理清单漏了 `#em-chem-links` 容器，
  D3 data join key 复用导致 `.enter()` 返回空集合，DOM 残留旧 `<line>` 绑着
  过期 chemical 引用。修补一行清理 + 加回归测试断言五个兄弟容器都被清理。
  （[6667ae3](https://github.com/redlizzxy/EnvMeta/commit/6667ae3)）

- **假说评分（macOS）**：上传 YAML 文件评分时报 `[Errno 63] File name too long: '#`。
  根因为 `load_hypothesis()` 对字符串入参先 `Path(src).exists()` 试探是否文件路径，
  YAML 正文常达数 KB → macOS `stat()` 抛 `ENAMETOOLONG`，Python 3.11 的
  `Path._IGNORED_ERRNOS` 白名单不含该 errno → 异常上浮到 UI。修法：启发式
  （含换行 / ≥1024 字符 → 直判 YAML 正文）+ try/except OSError 双重保护。
  （[8085d14](https://github.com/redlizzxy/EnvMeta/commit/8085d14)）

### 📚 文档

- **小白安装指南**新增 3 个 Mac 安装常见问题：
  - Q1.1：`CondaToSNonInteractiveError` —— Anaconda 2024 ToS 接受门槛
  - Q1.2：`ResolutionImpossible: protobuf` —— Apple Silicon pip 旧 resolver
    选到无 arm64 wheel 的 protobuf；解法是先 `pip install --upgrade pip`
  - Q1.3：首次启动 `Welcome to Streamlit! ... Email:` —— 留空回车跳过
  （[4742b6c](https://github.com/redlizzxy/EnvMeta/commit/4742b6c)）

- **README 新增"更新到新版本"章节**：本地用户内测期 `git pull origin master`
  完整工作流 + 3 个常见场景（Already up to date / 本地 stash 冲突 / 新依赖）。
  （[977ef83](https://github.com/redlizzxy/EnvMeta/commit/977ef83)）

### 🧹 维护

- 包版本 `__version__` 由历史遗留的 `0.1.0` 同步到 `0.8.1`，与 README 长期使用
  的 "v0.8" 内测版叙事对齐。HTML 导出顶部的 "EnvMeta v…" 徽章会显示新版本。
- 测试：291 → **293 全绿**（新增 2 个回归测试）

---

## [0.8.0] — 2026-04-19

> v0.8 内测版（"Sunday Sprint"）—— 投稿前的功能完整版。

### ✨ 新增

- **新手落地包（S8-ux）**：数据准备指南（11 上游工具 → EnvMeta 输入映射）+
  图表选择向导（按研究问题反查推荐分析）+ 14 图"如何解读"expander +
  样例数据一键加载
- **导出中心（T1）**：4 tab 统一入口（图表 / Bundle / 脚本 / 文档）+ 批量 ZIP
- **HTML 交互导出（T2）**：400-550 KB 独立 HTML（D3 v7 inline 嵌入）+
  per-group 切换 + 化学物精确锚点 + 化学式上下标格式化（v1.3 精修）
- **在线部署**：Streamlit Cloud 自动部署 + 本地 Windows / Linux 行为一致性修复
- **内测素材**：腾讯问卷 + 海报 + 部署指南 + 小白安装指南

### 📚 文档

- CLAUDE.md 拆分：核心规范保留在 [CLAUDE.md](CLAUDE.md)，历史日志移到
  [DEBUG_NOTES.md](DEBUG_NOTES.md)（精简前 1969 行 → 264 行）

---

## [0.7.0] — 2026-04-13 至 2026-04-18

> Phase 3 核心：循环图 + 假说评分 + Fork Bundle + KEGG-driven KB。

### ✨ 新增

- **生物地球化学循环图 v2**：Mockup 10 合并细胞布局 + 跨元素化学物耦合 +
  keystone ★ 标注 + 4 元素 (As/N/S/Fe) × 18 通路 × 57 KO 自动推断
- **统计可信度三件套（S1/S2/S3.5）**：去偏 + 999 次置换 null_p + Saltelli 权重敏感度
- **机制假说 YAML 评分器（S3）**：5 类 claim（pathway_active / coupling_possible /
  env_correlation / keystone_in_pathway / group_contrast）+ Bradford Hill required veto
- **Fork Bundle（S4）**：打包 KB + YAML + config + KEGG 快照 → zip
- **KEGG-driven KB（S5）**：`envmeta kb-build` CLI 从 KEGG snapshot 重建知识库

---

## [0.5.0] — 2026-04-08 前后

> Phase 1 + Phase 2：14 图分析引擎完整版。

### ✨ 新增

- **Phase 1**：7 张 Reads-based 图（堆叠图 / α / PCoA / RDA / LEfSe / 基因热图 /
  log2FC）+ 代码生成器
- **Phase 2**：5 张 MAG-based 图（MAG 质量 / 丰度热图 / 通路完整度 / 基因谱 /
  共现网络 Gephi-prep）+ 4 图共享统一 4 层参数面板

---

## [0.1.0] — 项目骨架

> Phase 0：环境 + 知识库 + 11 文件类型识别 + Streamlit 基础。

---

## 版本号约定

- **主版本号 (X.0.0)**：架构性变化（如 Phase 4 插件框架上线时升 1.0.0）
- **次版本号 (0.X.0)**：新功能 / 新分析图 / 新模块
- **补丁号 (0.0.X)**：bug 修复 / 文档 / 跨平台兼容修复

发表论文后会标记 v1.0.0 + 申请 Zenodo DOI 永久归档。
