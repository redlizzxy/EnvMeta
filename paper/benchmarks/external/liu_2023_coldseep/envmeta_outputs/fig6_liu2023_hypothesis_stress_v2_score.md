# Liu2023_coldseep_stress_hypothesis_v2 — Stress Hypothesis Score

- overall_score: **0.250**
- label: **WEAK**
- n_satisfied / n_total / n_skipped: 1 / 3 / 1
- null_p: 1.0000 (n=999 permutations)
- weight_robust: True (OAT ±20%)
- veto_reasons: none

## Per-claim results

| ID | Type | Status | Score | Weight |
|---|---|---|---|---|
| as_oxidation_should_dominate_stress_A | pathway_active | unsatisfied | 0.00 | 1.5 |
| nitrification_should_dominate_stress_B | pathway_active | skipped | 0.00 | 1.0 |
| arsenate_reduction_should_NOT_dominate_stress_C | pathway_inactive | unsatisfied | 0.00 | 1.5 |
| as_transport_calibration_anchor_D | pathway_active | satisfied | 1.00 | 1.0 |

## Evidence per claim

### as_oxidation_should_dominate_stress_A

- status=unsatisfied score=0.00 weight=1.5
- evidence: `{'pathway_id': 'Arsenite oxidation', 'element': 'arsenic', 'n_active_mags': 2, 'mean_completeness': 50.0, 'total_contribution': 0.3, 'dominance_score': 0.0005, 'element_total_contribution': 618.67}`
- evaluator note: Arsenite oxidation: 活跃但**不主导** (dominance 0.05%/20%, contribution 0.3 vs element total 618.7)

### nitrification_should_dominate_stress_B

- status=skipped score=0.00 weight=1.0
- evidence: `{'pathway_query': 'Nitrification'}`
- evaluator note: 数据里找不到通路 'Nitrification'（可能被 element_filter 过滤）

### arsenate_reduction_should_NOT_dominate_stress_C

- status=unsatisfied score=0.00 weight=1.5
- evidence: `{'pathway_id': 'Arsenate reduction', 'element': 'arsenic', 'n_active_mags': 18, 'mean_completeness': 50.0, 'total_contribution': 6.43, 'max_completeness_threshold': 50.0, 'dominance_score': 0.0104, 'element_total_contribution': 618.67}`
- evaluator note: Arsenate reduction: 18 活跃 MAG，平均完整度 50%（≥50%） — 明显违反 negative 预期：通路实际活跃

### as_transport_calibration_anchor_D

- status=satisfied score=1.00 weight=1.0
- evidence: `{'pathway_id': 'As transport/detox', 'element': 'arsenic', 'n_active_mags': 1, 'mean_completeness': 50.0, 'total_contribution': 0.26, 'dominance_score': 0.0004, 'element_total_contribution': 618.67}`
- evaluator note: As transport/detox: 1 个 MAG 活跃，平均完整度 50%，总贡献 0.3
