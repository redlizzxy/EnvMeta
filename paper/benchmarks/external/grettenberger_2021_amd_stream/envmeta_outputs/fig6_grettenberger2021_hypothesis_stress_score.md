# Grettenberger2021_amd_stream_stress_hypothesis — Stress Hypothesis Score

- overall_score: **0.250**
- label: **WEAK**
- n_satisfied / n_total / n_skipped: 1 / 3 / 1
- null_p: 1.0000 (n=999 permutations)
- weight_robust: True (OAT ±20%)
- veto_reasons: none

## Per-claim results

| ID | Type | Status | Score | Weight |
|---|---|---|---|---|
| nitrification_should_dominate_stress_A | pathway_active | skipped | 0.00 | 1.5 |
| arsenate_reduction_should_dominate_stress_B | pathway_active | unsatisfied | 0.00 | 1.5 |
| sulfide_oxidation_should_NOT_dominate_stress_C | pathway_inactive | unsatisfied | 0.00 | 1.5 |
| dissim_sulfate_reduction_calibration_anchor_D | pathway_active | satisfied | 1.00 | 1.0 |

## Evidence per claim

### nitrification_should_dominate_stress_A

- status=skipped score=0.00 weight=1.5
- evidence: `{'pathway_query': 'Nitrification'}`
- evaluator note: 数据里找不到通路 'Nitrification'（可能被 element_filter 过滤）

### arsenate_reduction_should_dominate_stress_B

- status=unsatisfied score=0.00 weight=1.5
- evidence: `{'pathway_id': 'Arsenate reduction', 'element': 'arsenic', 'n_active_mags': 0, 'mean_completeness': 0.0, 'total_contribution': 0}`
- evaluator note: Arsenate reduction: 无活跃 MAG

### sulfide_oxidation_should_NOT_dominate_stress_C

- status=unsatisfied score=0.00 weight=1.5
- evidence: `{'pathway_id': 'Sulfide oxidation', 'element': 'sulfur', 'n_active_mags': 4, 'mean_completeness': 83.33, 'total_contribution': 333.33, 'max_completeness_threshold': 50.0}`
- evaluator note: Sulfide oxidation: 4 活跃 MAG，平均完整度 83%（≥50%） — 明显违反 negative 预期：通路实际活跃

### dissim_sulfate_reduction_calibration_anchor_D

- status=satisfied score=1.00 weight=1.0
- evidence: `{'pathway_id': 'Dissim. sulfate red.', 'element': 'sulfur', 'n_active_mags': 2, 'mean_completeness': 75.0, 'total_contribution': 150.0}`
- evaluator note: Dissim. sulfate red.: 2 个 MAG 活跃，平均完整度 75%，总贡献 150.0
