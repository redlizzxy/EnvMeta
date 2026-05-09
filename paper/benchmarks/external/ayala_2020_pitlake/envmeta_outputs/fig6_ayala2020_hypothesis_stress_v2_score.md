# Ayala2020_pitlake_deep_stress_hypothesis_v2 — Stress Hypothesis Score

- overall_score: **0.182**
- label: **WEAK**
- n_satisfied / n_total / n_skipped: 1 / 4 / 0
- null_p: 1.0000 (n=999 permutations)
- weight_robust: True (OAT ±20%)
- veto_reasons: none

## Per-claim results

| ID | Type | Status | Score | Weight |
|---|---|---|---|---|
| sulfide_oxidation_should_dominate_stress_A | pathway_active | unsatisfied | 0.00 | 1.5 |
| arsenate_reduction_should_dominate_stress_B | pathway_active | unsatisfied | 0.00 | 1.5 |
| dissim_sulfate_reduction_should_NOT_dominate_stress_C | pathway_inactive | unsatisfied | 0.00 | 1.5 |
| nitrate_reduction_calibration_anchor_D | pathway_active | satisfied | 1.00 | 1.0 |

## Evidence per claim

### sulfide_oxidation_should_dominate_stress_A

- status=unsatisfied score=0.00 weight=1.5
- evidence: `{'pathway_id': 'Sulfide oxidation', 'element': 'sulfur', 'n_active_mags': 1, 'mean_completeness': 66.67, 'total_contribution': 66.67, 'dominance_score': 0.0708, 'element_total_contribution': 941.67}`
- evaluator note: Sulfide oxidation: 活跃但**不主导** (dominance 7.08%/20%, contribution 66.7 vs element total 941.7)

### arsenate_reduction_should_dominate_stress_B

- status=unsatisfied score=0.00 weight=1.5
- evidence: `{'pathway_id': 'Arsenate reduction', 'element': 'arsenic', 'n_active_mags': 0, 'mean_completeness': 0.0, 'total_contribution': 0, 'dominance_score': 0.0, 'element_total_contribution': 850.0}`
- evaluator note: Arsenate reduction: 无活跃 MAG

### dissim_sulfate_reduction_should_NOT_dominate_stress_C

- status=unsatisfied score=0.00 weight=1.5
- evidence: `{'pathway_id': 'Dissim. sulfate red.', 'element': 'sulfur', 'n_active_mags': 8, 'mean_completeness': 84.38, 'total_contribution': 675.0, 'max_completeness_threshold': 50.0, 'dominance_score': 0.7168, 'element_total_contribution': 941.67}`
- evaluator note: Dissim. sulfate red.: 8 活跃 MAG，平均完整度 84%（≥50%） — 明显违反 negative 预期：通路实际活跃

### nitrate_reduction_calibration_anchor_D

- status=satisfied score=1.00 weight=1.0
- evidence: `{'pathway_id': 'Nitrate reduction', 'element': 'nitrogen', 'n_active_mags': 3, 'mean_completeness': 50.0, 'total_contribution': 150.0, 'dominance_score': 0.4286, 'element_total_contribution': 350.0}`
- evaluator note: Nitrate reduction: 3 个 MAG 活跃，平均完整度 50%，总贡献 150.0
