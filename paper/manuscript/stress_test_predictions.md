# Stress Test 盲预测（Pre-registration before run）

> **创建日期**：2026-05-08（在跑任一 stress YAML **之前**写）
> **目的**：把所有 stress YAML 的预测和理由**预先固化**，跑完后只能 diff，不能改 YAML
> **关联**：[docs/hypothesis_writing_guide.md §3 Pre-prediction](../../docs/hypothesis_writing_guide.md)
> **配套 YAMLs**：
> - [liu2023_hypothesis_stress.yaml](../benchmarks/external/liu_2023_coldseep/liu2023_hypothesis_stress.yaml)
> - [grettenberger2021_hypothesis_stress.yaml](../benchmarks/external/grettenberger_2021_amd_stream/grettenberger2021_hypothesis_stress.yaml)
> - [ayala2020_hypothesis_stress.yaml](../benchmarks/external/ayala_2020_pitlake/ayala2020_hypothesis_stress.yaml)

---

## 0. 全局预测

| Stress YAML | 期望 overall_score | 期望 label |
|---|---|---|
| Liu 2023 stress | ~0.20 | weak |
| Grettenberger 2021 stress | ~0.20 | weak |
| Ayala 2020 stress | ~0.20 | weak |

**计算依据**：每 YAML 含 3 stress (期望 unsatisfied=0) + 1 calibration anchor
(期望 satisfied=1.0)，weights = [1.5, 1.0, 1.5, 1.0] 总和 5.0：
overall = (1.5×0 + 1.0×0 + 1.5×0 + 1.0×1.0) / 5.0 = **0.20**

**discrimination 失败信号**：overall ≥ 0.40 (suggestive 或更高)。

---

## 1. Liu 2023 冷泉 stress（4 claims）

| Claim ID | type | expected_status | risk level | reasoning（精简） |
|---|---|---|---|---|
| as_oxidation_should_dominate_stress_A | pathway_active | **unsatisfied** | low | 冷泉缺氧应 As 还原主导，氧化反向 |
| nitrification_should_dominate_stress_B | pathway_active | **unsatisfied** | low | 冷泉缺氧抑制 amoA/B/C |
| arsenate_reduction_should_NOT_dominate_stress_C | pathway_inactive | **unsatisfied** | low-med | backbone As 还原 pathway_inactive 否定，应被 reject |
| as_transport_calibration_anchor_D | pathway_active | **satisfied** | n/a | universal As efflux backbone |

**意外结果分析预案**：
- **A satisfied**：冷泉数据意外含 aoxAB carrier → 不修改 YAML，在 results 中报告并查 KEGG 注释
- **B satisfied**：上层水柱混入或 amoA lateral transfer → 同上
- **C satisfied (= inactive)**：与 calibration YAML 期望矛盾，需要重读 calibration 结果。calibration YAML Liu 已得 STRONG 4/4，所以 As 还原应活跃；如 stress YAML 这条 satisfied (inactive)，说明 stress YAML 跑的是不同 cycle_data 实例 → 排查跑分脚本
- **D unsatisfied**：calibration anchor fail → 整个 stress YAML 跑分系统问题，停下排查

---

## 2. Grettenberger 2021 AMD stress（4 claims）

| Claim ID | type | expected_status | risk level | reasoning |
|---|---|---|---|---|
| nitrification_should_dominate_stress_A | pathway_active | **unsatisfied** | low | AMD pH 2-4 抑制 amoA/B/C |
| arsenate_reduction_should_dominate_stress_B | pathway_active | **unsatisfied** | **med-high** | AMD 无砷主题，但 arsC universal 解毒可能带 |
| sulfide_oxidation_should_NOT_dominate_stress_C | pathway_inactive | **unsatisfied** | low | backbone sulfide ox 否定，应被 reject |
| dissim_sulfate_reduction_calibration_anchor_D | pathway_active | **satisfied** | n/a | 沉积物缺氧 SRB |

**意外结果分析预案**：
- **B satisfied**（最值得关注）：confirm 假说 "arsC 是 universal 解毒型，跨主题数据集都可能 active"
  - 这种"意外" satisfied 不是 EnvMeta 鉴别力问题，而是 **领域错配的 false positive 来自基因功能 promiscuity**
  - 论文 Discussion 中应专门讨论：用户写 cross-topic stress claim 时要避开 universal 解毒型基因
- **A/C satisfied**：数据问题或 KB 注释 leakage，需查
- **D unsatisfied**：anchor fail，停排查

---

## 3. Ayala 2020 pit lake stress（4 claims，待 GhostKOALA）

| Claim ID | type | expected_status | risk level | reasoning |
|---|---|---|---|---|
| sulfide_oxidation_should_dominate_stress_A | pathway_active | **unsatisfied** | med | 深层缺氧但 chemocline 边缘可能微氧 |
| arsenate_reduction_should_dominate_stress_B | pathway_active | **unsatisfied** | **med-high** | 同 Grettenberger，arsC universal 风险 |
| dissim_sulfate_reduction_should_NOT_dominate_stress_C | pathway_inactive | **unsatisfied** | low | backbone dsrAB 否定 |
| nitrate_reduction_calibration_anchor_D | pathway_active | **satisfied** | n/a | 寡营养深层 narG/napA |

---

## 4. 预测 vs 实测 diff 表（跑完后填）

跑完 stress test 后，把每条 claim 的 observed_status 填进**新文件**
`stress_test_results.md`（不要改本预测文件，pre-registration 必须冻结）。

模板：

| dataset | claim_id | expected | observed | match | note |
|---|---|---|---|---|---|
| Liu | as_oxidation_should_dominate_stress_A | unsatisfied | ? | ? | |
| Liu | nitrification_should_dominate_stress_B | unsatisfied | ? | ? | |
| Liu | arsenate_reduction_should_NOT_dominate_stress_C | unsatisfied | ? | ? | |
| Liu | as_transport_calibration_anchor_D | satisfied | ? | ? | |
| Grettenberger | nitrification_should_dominate_stress_A | unsatisfied | ? | ? | |
| Grettenberger | arsenate_reduction_should_dominate_stress_B | unsatisfied | ? | ? | ⚠️ universal arsC 风险 |
| Grettenberger | sulfide_oxidation_should_NOT_dominate_stress_C | unsatisfied | ? | ? | |
| Grettenberger | dissim_sulfate_reduction_calibration_anchor_D | satisfied | ? | ? | |
| Ayala | sulfide_oxidation_should_dominate_stress_A | unsatisfied | ? | ? | |
| Ayala | arsenate_reduction_should_dominate_stress_B | unsatisfied | ? | ? | ⚠️ universal arsC 风险 |
| Ayala | dissim_sulfate_reduction_should_NOT_dominate_stress_C | unsatisfied | ? | ? | |
| Ayala | nitrate_reduction_calibration_anchor_D | satisfied | ? | ? | |

---

## 5. 关键纪律：跑完后 NOT 做的事

1. ❌ 不要修改 stress YAML 的 weight / required / threshold "解决"意外结果
2. ❌ 不要新增"hindsight" claim 解读结果
3. ❌ 不要在 stress_test_results.md 里改 expected（这是 pre-prediction 文件冻结）
4. ✅ 意外结果如实分析，写进 stress_test_results.md 的 "Unexpected outcomes" 章节
5. ✅ KB v1.2 backlog 候选项可以记下（如 "AMD 数据集 arsC 通用解毒应 KB 区分 detox vs respiratory"）
6. ✅ 如果 stress 整体 calibration anchor (D) 都 fail → 停下排查跑分脚本，不是 stress YAML 问题

---

## 6. 维护

| 日期 | 事项 |
|---|---|
| 2026-05-08 | 初版 — 12 条 claim 盲预测固化（3 dataset × 4 claim）|
