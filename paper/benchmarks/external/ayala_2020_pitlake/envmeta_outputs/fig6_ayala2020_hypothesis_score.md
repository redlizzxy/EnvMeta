# Ayala2020_pitlake_deep_pre-data_hypothesis — Stress Hypothesis Score

- overall_score: **1.000**
- label: **STRONG**
- n_satisfied / n_total / n_skipped: 4 / 4 / 0
- null_p: n/a
- weight_robust: True (OAT ±20%)
- veto_reasons: none

## Per-claim results

| ID | Type | Status | Score | Weight |
|---|---|---|---|---|
| dissim_sulfate_reduction_active_in_deep_anoxic | pathway_active | satisfied | 1.00 | 1.5 |
| sulfide_oxidation_microaerophilic | pathway_active | satisfied | 1.00 | 1.0 |
| nitrate_reduction_explored | pathway_active | satisfied | 1.00 | 0.5 |
| nitrogen_fixation_explored | pathway_active | satisfied | 1.00 | 0.3 |

## Evidence per claim

### dissim_sulfate_reduction_active_in_deep_anoxic

- status=satisfied score=1.00 weight=1.5
- evidence: `{'pathway_id': 'Dissim. sulfate red.', 'element': 'sulfur', 'n_active_mags': 8, 'mean_completeness': 84.38, 'total_contribution': 675.0}`
- evaluator note: Dissim. sulfate red.: 8 个 MAG 活跃，平均完整度 84%，总贡献 675.0

### sulfide_oxidation_microaerophilic

- status=satisfied score=1.00 weight=1.0
- evidence: `{'pathway_id': 'Sulfide oxidation', 'element': 'sulfur', 'n_active_mags': 1, 'mean_completeness': 66.67, 'total_contribution': 66.67}`
- evaluator note: Sulfide oxidation: 1 个 MAG 活跃，平均完整度 67%，总贡献 66.7

### nitrate_reduction_explored

- status=satisfied score=1.00 weight=0.5
- evidence: `{'pathway_id': 'Nitrate reduction', 'element': 'nitrogen', 'n_active_mags': 3, 'mean_completeness': 50.0, 'total_contribution': 150.0}`
- evaluator note: Nitrate reduction: 3 个 MAG 活跃，平均完整度 50%，总贡献 150.0

### nitrogen_fixation_explored

- status=satisfied score=1.00 weight=0.3
- evidence: `{'pathway_id': 'N fixation', 'element': 'nitrogen', 'n_active_mags': 2, 'mean_completeness': 100.0, 'total_contribution': 200.0}`
- evaluator note: N fixation: 2 个 MAG 活跃，平均完整度 100%，总贡献 200.0
