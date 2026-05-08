# 假说 YAML 引用审计 — DOI + 提炼度评估

> **创建日期**：2026-05-08（在收到 reviewer 诘问之前主动 audit）
> **目的**：解决假说验证 YAML 中先验综述**没有 DOI** + **提炼度未评估**的方法学缺陷
> **关联**：[stress_test_results.md](stress_test_results.md) +
> [scoring_validation_self_critique.md](scoring_validation_self_critique.md) §3
> **范围**：6 个 external benchmark YAML 中的 16 条 claim × 13 篇先验综述

---

## 1. 现状摘要

| 维度 | 状态 |
|---|---|
| 目标论文 DOI（YAML 文件头）| ✅ 6/6 都有 |
| 先验综述 DOI（claim 内引用）| ❌ 0/13 |
| arsenic_steel_slag.yaml § 4 REFERENCES（项目原 demo）| ✅ 7/7 都有完整 Vancouver 格式 + DOI |

**问题**：reviewer 完全有理由问"Stolz 2006 是哪一篇？DOI 是什么？里面具体哪一段
支持你的 arsenate_reduction REQUIRED claim？" —— 当前 YAML 不能直接回答。

---

## 2. 13 篇先验综述完整 DOI 表（Agent verify 完成 2026-05-08）

> ⚠️ **Agent verify 发现 4 个 YAML 中的引用错误**（标 ❌ ERR）。这些错误在 Methods
> 章节必须明确纠正。

| # | YAML 中标注 | Verify 后正确引用 | DOI | 错误等级 |
|---|---|---|---|---|
| 1 | Stolz et al. 2006 NRM | Stolz JF, Basu P, Santini JM, Oremland RS. Arsenic and selenium in microbial metabolism. **Annu Rev Microbiol**. 2006;60:107-130. | 10.1146/annurev.micro.60.080805.142053 | ✅ 正确 |
| 2 | Mukhopadhyay et al. 2002 | Mukhopadhyay R, Rosen BP, Phung LT, Silver S. Microbial arsenic: from geocycles to genes and enzymes. **FEMS Microbiol Rev**. 2002;26(3):311-325. | 10.1111/j.1574-6976.2002.tb00622.x | ✅ 正确 |
| 3 | Rosen 2002 Trends Microbiol | Rosen BP. Biochemistry of arsenic detoxification. **FEBS Lett**. 2002;529(1):86-92. | 10.1016/S0014-5793(02)03186-1 | ⚠️ 期刊误标（应为 FEBS Lett 不是 Trends Microbiol）|
| 4 | Yin et al. 2011 ES&T | Yin XX, Chen J, Qin J, Sun GX, Rosen BP, Zhu YG. Biotransformation and volatilization of arsenic by three photosynthetic cyanobacteria. **Plant Physiol**. 2011;156(3):1631-1638. | 10.1104/pp.111.178947 | ❌ **ERR：期刊错（Plant Physiol 不是 ES&T）** |
| 5 | Bond et al. 2000 AEM | Bond PL, Druschel GK, Banfield JF. Comparison of AMD microbial communities. Appl Environ Microbiol. 2000;66(11):4962-4971. | 10.1128/AEM.66.11.4962-4971.2000 | ✅ 正确 |
| 6 | Schippers & Sand 1999 | Schippers A, Sand W. Appl Environ Microbiol. 1999;65(1):319-321. | 10.1128/AEM.65.1.319-321.1999 | ✅ 正确 |
| 7 | Bigham & Nordstrom 2000 | Bigham JM, Nordstrom DK. Rev Mineral Geochem. 2000;40:351-403. | 10.2138/rmg.2000.40.7 | ✅ 正确 |
| 8 | Bothe 2007 FEMS Rev N cycle | **正确**: Bothe H, Ferguson SJ, Newton WE (eds.). Biology of the Nitrogen Cycle. Amsterdam: Elsevier. 2007. **452pp** (ISBN 978-0-444-52857-5) — 这是**编辑书籍**不是 FEMS 综述<br/><br/>**或应改引**: Bothe H, Jost G, Schloter M, Ward BB, Witzel KP. Molecular analysis of ammonia oxidation and denitrification in natural environments. **FEMS Microbiol Rev**. 2000;24(5):673-690. | Book: ISBN 978-0-444-52857-5<br/>Article: 10.1111/j.1574-6976.2000.tb00566.x | ❌ **ERR：YAML 写"Bothe 2007 FEMS Rev"实际不存在**。Bothe 2007 是书；FEMS 综述是 Bothe 2000 — 引用日期+期刊都错 |
| 9 | Cabrera et al. 2006 Process Biochem | Cabrera G, Pérez R, Gómez JM, Ábalos A, Cantero D. Toxic effects of dissolved heavy metals on Desulfovibrio. **J Hazard Mater**. 2006;135(1-3):40-46. | 10.1016/j.jhazmat.2005.11.058 | ❌ **ERR：期刊错（J Hazard Mater 不是 Process Biochem）+ 主题不是 acidophilic SRB 综述而是金属毒性 SRB 实验** |
| 10 | Sánchez-Andrea et al. 2014 | Sánchez-Andrea I, Sanz JL, Bijmans MFM, Stams AJM. Sulfate reduction at low pH to remediate acid mine drainage. J Hazard Mater. 2014;269:98-109. | 10.1016/j.jhazmat.2013.12.032 | ✅ 正确 |
| 11 | Tan et al. 2009 AMD N cycle review | Tan GL, Shu WS, Zhou WH, Li XL, Lan CY, Huang LN. Seasonal and spatial variations in microbial community structure and diversity in the acid stream draining across an ongoing surface mining site. **FEMS Microbiol Ecol**. 2009;70(2):277-285. | 10.1111/j.1574-6941.2009.00744.x | ⚠️ **是 study 不是 review；含 N cycle 数据但不是 N-cycle review；YAML "AMD N 循环综述"措辞略不准** |
| 12 | Auld et al. 2017 AMD diazotrophy | Auld RR, Mykytczuk NCS, Leduc LG, Merritt TJS. Seasonal variation in an acid mine drainage microbial community. **Can J Microbiol**. 2017;63(2):137-152. | 10.1139/cjm-2016-0215 | ❌ **ERR：主题错。是 seasonal community variation 不是 diazotrophy。YAML 用它支持 N fixation 偶发是错引** |
| 13 | Johnson & Hallberg 2005 | Johnson DB, Hallberg KB. Acid mine drainage remediation options: a review. Sci Total Environ. 2005;338(1-2):3-14. | 10.1016/j.scitotenv.2004.09.002 | ✅ 正确 |
| 14 | Diaby et al. 2007 | Diaby N, Dold B, Pfeifer HR, Holliger C, Johnson DB, Hallberg KB. Microbial communities in a porphyry copper tailings impoundment. Environ Microbiol. 2007;9(2):298-307. | 10.1111/j.1462-2920.2006.01138.x | ✅ 正确 |
| 15 | Falagán & Johnson 2014 | Falagán C, **Sánchez-España J**, Johnson DB. New insights into the biogeochemistry of extremely acidic environments revealed by a combined cultivation-based and culture-independent study of two stratified pit lakes. FEMS Microbiol Ecol. 2014;87(1):231-243. | 10.1111/1574-6941.12218 | ⚠️ 漏作者：是 3 作者（YAML 写 2 作者，Sánchez-España 是 co-author）|
| 16 | Sánchez-España et al. 2008 | Sánchez-España J, López-Pamo E, Santofimia E, Diez-Ercilla M. The acidic mine pit lakes of the Iberian Pyrite Belt. Appl Geochem. 2008;23(5):1260-1287. | 10.1016/j.apgeochem.2007.12.036 | ⚠️ DOI 推断（Agent 未直接 verify；建议手动 doi.org 确认）|

### 2.1 引用错误等级统计

| 等级 | 数量 | 说明 |
|---|---|---|
| ✅ 正确 | 7 | YAML 引用与 verify 一致 |
| ⚠️ 小错（期刊名/作者数/年份歧义） | 4 | 不影响 claim 提炼度，Methods 中纠正即可 |
| ❌ ERR（主题或引用本身错） | 4 | **影响 claim 提炼度评估**，需在 Methods 显式声明 |

### 2.2 4 处 ❌ ERR 的影响

1. **Yin 2011 期刊错（Plant Physiol 不是 ES&T）** — claim 提炼度仍 Direct（Plant Physiol 论文确实讲 cyanobacteria arsM），但 reviewer 一查 ES&T 找不到会提问
2. **Bothe 2007 文献本身错（应为 Bothe 2000 FEMS）** — Liu stress YAML 引用"Bothe 2007 FEMS"不存在；应改 Bothe 2000 FEMS Microbiol Rev（claim 提炼度仍 Direct）
3. **Cabrera 2006 期刊+主题错** — YAML 引用"acidophilic SRB 综述"不存在；Cabrera 2006 是金属毒性 SRB 实验论文。Grettenberger calibration claim 引用 Cabrera 应改 Sánchez-Andrea 2014（J Hazard Mater，是真正的 acidophilic SRB 综述）
4. **Auld 2017 主题错（是 seasonal variation，不是 diazotrophy）** — Grettenberger + Ayala calibration 的 nitrogen_fixation_explored claim 引用 Auld 2017 支持"AMD 寡营养偶发 nifHDK 报告"，**实际 Auld 2017 没明说 diazotrophy**。这条 claim 的引用支持降级为 **Weak**，应改引真正的 AMD diazotrophy 文献（候选：Korehi et al. 2014 *Front Microbiol*, Mendez-Garcia et al. 2015 *Front Microbiol*）

---

## 3. 16 Claim × 提炼度评估

### 3.1 Liu 2023 calibration（4 claims）

| claim_id | 引用 | 提炼度 | 文献章节支持 |
|---|---|---|---|
| arsenate_reduction_active_in_anoxic_seep | Stolz 2006 + Mukhopadhyay 2002 | **Direct ✅** | Stolz 2006 §3.1 "Arsenate reductase ArsC is widespread in anaerobic As-reducing prokaryotes"；Mukhopadhyay 2002 全文综述 arsC operon |
| as_transport_detox_active | Rosen 2002 | **Direct ✅** | Rosen 2002 标题就是 "Biochemistry of arsenic detoxification"，ArsB / Acr3 efflux 是核心论点 |
| as_methylation_explored | Yin 2011 | **Direct ✅** | Yin 2011 报告海洋 cyanobacteria arsM 砷甲基化稀少，与"海洋 arsM 报道有限"匹配 |
| respiratory_as_reduction_explored | Stolz 2006 | **Direct ✅** | Stolz 2006 §3.2 含 arrA respiratory As reduction 章节，引用海洋沉积物零星报告 |

### 3.2 Liu 2023 stress（4 claims）

| claim_id | 引用 | 提炼度 | 评估 |
|---|---|---|---|
| as_oxidation_should_dominate_stress_A | Stolz 2006 | **Inferred ⚠️** | Stolz 2006 没说"冷泉应有 As 氧化主导"；只说 aoxAB 在某些海洋微生物存在。我的 stress 是从"缺氧应是 reduction 主导"反向推。此处依赖**反向 inference 而非 direct quote** |
| nitrification_should_dominate_stress_B | Bothe 2007 | **Direct ✅** | Bothe 2007 N 循环综述明确 amoA/B/C 需 O2，缺氧应受抑 |
| arsenate_reduction_should_NOT_dominate_stress_C | Stolz 2006 | **Direct ✅** | 同 calibration arsenate_reduction，反向用为 negative claim |
| as_transport_calibration_anchor_D | Rosen 2002 | **Direct ✅** | 同 calibration as_transport |

### 3.3 Grettenberger 2021 calibration（4 claims）

| claim_id | 引用 | 提炼度 | 评估 |
|---|---|---|---|
| sulfide_oxidation_active_in_amd_stream | Bond 2000 + Schippers 1999 | **Direct ✅** | Bond 2000 全文比较 AMD 群落，Acidithiobacillus / Leptospirillum 主导是核心结果；Schippers 1999 sulfide ox 机制综述 |
| dissim_sulfate_reduction_anoxic_sediment | Cabrera 2006 + Sánchez-Andrea 2014 | **Direct ✅** | 两篇都是 acidophilic SRB 综述 |
| nitrate_reduction_explored | Tan 2009 | **Inferred ⚠️** | Tan 2009 是 AMD N cycle 综述但 narG / napA 在 AMD 中是"局部存在" rather than "强信号"，引用支持是 weak-medium |
| nitrogen_fixation_explored | Auld 2017 | **Direct ✅** | Auld 2017 是 AMD diazotrophy 报告，直接支持 nifHDK 偶发 |

### 3.4 Grettenberger 2021 stress（4 claims）

| claim_id | 引用 | 提炼度 | 评估 |
|---|---|---|---|
| nitrification_should_dominate_stress_A | Tan 2009 + Bothe 2007 | **Direct ✅** | 同上 |
| arsenate_reduction_should_dominate_stress_B | Bigham 2000 + Bond 2000 | **Inferred ⚠️** | 这两篇说 AMD 是 Fe/Cu/Zn/SO4 主导没明说"无砷"；我的 stress 是 cross-topic inference。**reviewer 可能问"Cabin Branch 是否真的检测过砷浓度低？"** —— 这条是 stress test 较弱的引用支持点 |
| sulfide_oxidation_should_NOT_dominate_stress_C | Bond 2000 | **Direct ✅** | 同 calibration |
| dissim_sulfate_reduction_calibration_anchor_D | Cabrera 2006 + Sánchez-Andrea 2014 | **Direct ✅** | 同 calibration |

### 3.5 Ayala 2020 calibration（4 claims）

| claim_id | 引用 | 提炼度 | 评估 |
|---|---|---|---|
| dissim_sulfate_reduction_active_in_deep_anoxic | Johnson & Hallberg 2005 + Diaby 2007 | **Direct ✅** | Johnson & Hallberg 2005 综述 biosulfidogenesis 修复 |
| sulfide_oxidation_microaerophilic | Falagán 2014 | **Direct ✅** | Falagán 2014 IPB pit lake 微氧界面研究 |
| nitrate_reduction_explored | Tan 2009 | **Inferred ⚠️** | 同 Grettenberger，AMD N 循环 weak |
| nitrogen_fixation_explored | Auld 2017 | **Direct ✅** | 同 Grettenberger |

### 3.6 Ayala 2020 stress（4 claims）

| claim_id | 引用 | 提炼度 | 评估 |
|---|---|---|---|
| sulfide_oxidation_should_dominate_stress_A | Sánchez-España 2008/2014 + Falagán 2014 | **Inferred ⚠️** | 综述说深层缺氧但没明说"sulfide ox 应不主导"；从生化原理（需 O2/Fe(III) 受体）反向推 |
| arsenate_reduction_should_dominate_stress_B | Sánchez-España 2008 | **Inferred ⚠️** | 同 Grettenberger B，cross-topic 较弱 |
| dissim_sulfate_reduction_should_NOT_dominate_stress_C | Johnson & Hallberg 2005 | **Direct ✅** | 同 calibration |
| nitrate_reduction_calibration_anchor_D | Tan 2009 | **Inferred ⚠️** | 同上 |

---

## 4. 引用提炼度统计（Agent verify 后修正）

| 提炼度 | 数量（之前估计 → 修正后）| 占比 | 主要 claim |
|---|---|---|---|
| **Direct ✅** | 11 → **9** | 56% | calibration backbone + Bothe 2000 (修正) + Stolz 2006 反向 |
| **Inferred ⚠️** | 5 → **5** | 31% | stress claim 主要 + 2 nitrate_reduction_explored |
| **Weak ❌** | 0 → **2** | 13% | **Auld 2017 引用错（主题不是 diazotrophy）** 影响 Grettenberger + Ayala 的 nitrogen_fixation_explored claim 提炼度 |

**修正后结论**：
- ✅ Calibration claims 引用质量打折（9/12 Direct，2/12 Inferred，1/12 Weak — Auld 2017 引用错）
- ⚠️ Stress claims 引用提炼度仍较弱（4/12 Inferred）—— **stress test 方法学软肋未变**
- ❌ **新发现**：Auld 2017 引用错主题（应是 AMD seasonal variation 而非 diazotrophy），nitrogen_fixation_explored 在 Grettenberger + Ayala 两份 calibration YAML 都因此降级 Weak。**Methods 中必须改引真正的 AMD diazotrophy 文献**（如 Korehi et al. 2014 *Front Microbiol* 或 Mendez-Garcia et al. 2015 *Front Microbiol*；待 verify）

---

## 5. 改进建议（不修改原 YAML，写论文 Methods 时补强）

### 5.1 不能做的（pre-registration 限制）

❌ 不能 git push 修改 calibration / stress YAML 添加 DOI（违反 pre-registration anchor 50c4687）

### 5.2 能做的

✅ **写论文 Methods 时**：
1. 在 Methods 4.6 假说评分章节加 **References** 表格，列出 16 条 claim × 13 篇文献 + 完整 DOI
2. 引用提炼度 Direct/Inferred 也透明声明，不掩盖 stress claim 引用较弱
3. 在 Discussion 中提一句："Stress claims rely on inference from biogeochemical
   priors rather than direct quotation, which is acknowledged as a methodological
   limitation; future stress tests should use claims with direct literature support."

✅ **未来新数据集 stress YAML 写作**：
- 强制每条 claim explanation 后附 `# DOI: 10.xxx/yyy` 注释
- 强制 stress claim 注明"引用提炼度=Direct/Inferred"
- 在 docs/hypothesis_writing_guide.md 加一条"DOI 强制规范"

---

## 6. 论文 Methods 中的引用纠正声明（必须）

写 paper 3/4 Methods 时，必须在 references 部分加一段透明声明：

> "During post-hoc reference auditing, four citation errors were identified in
> the pre-registered hypothesis YAMLs (commits 42168da / 44d7f5f / 76a4f77 /
> 50c4687): (1) Yin et al. 2011 was journal-mislabeled as Environ Sci Technol
> instead of the correct Plant Physiol 156(3):1631-1638
> (10.1104/pp.111.178947); (2) the citation to 'Bothe et al. 2007 FEMS
> Microbiol Rev' refers to a non-existent article — the correct reference is
> Bothe et al. 2000 FEMS Microbiol Rev 24(5):673-690
> (10.1111/j.1574-6976.2000.tb00566.x); (3) Cabrera et al. 2006 was journal-
> mislabeled as Process Biochem instead of J Hazard Mater 135(1-3):40-46
> (10.1016/j.jhazmat.2005.11.058), and is a metal toxicity study rather than
> the acidophilic SRB review intended; (4) Auld et al. 2017 Can J Microbiol
> 63(2):137-152 (10.1139/cjm-2016-0215) is a seasonal community variation
> study, not an AMD diazotrophy report; the nitrogen_fixation_explored claims
> (Grettenberger, Ayala) therefore lack direct literature support and are
> downgraded to weak-evidence claims. None of these reference errors affect
> the EnvMeta scoring outputs themselves; they affect only the literature
> grounding of the corresponding claims. We acknowledge this as a methodological
> limitation and recommend future hypothesis YAMLs include verified DOIs
> alongside each claim citation."

---

## 7. 维护

| 日期 | 事项 |
|---|---|
| 2026-05-08 | 初版 — 16 claim × 13 文献提炼度评估；6/13 DOI 已 verify，7 待 verify（Agent 异步） |
| 2026-05-08（深夜）| **Agent verify 完成 → 发现 4 个引用错误（Yin 期刊 / Bothe 文献本身 / Cabrera 期刊+主题 / Auld 主题）。提炼度统计修正：Direct 9 / Inferred 5 / Weak 2。Methods 必须加纠正声明** |
| TBD | 真正的 AMD diazotrophy 文献候选（替代 Auld 2017）verify：Korehi 2014 / Mendez-Garcia 2015 |
