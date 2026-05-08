# Grettenberger2021_amd_stream_pre-data_hypothesis — Hypothesis Score

- overall_score: **1.000**
- label: **STRONG**
- n_satisfied / n_total / n_skipped: 4 / 4 / 0
- null_p: n/a
- weight_robust: True (OAT ±20%)
- veto_reasons: none

## Per-claim results

| ID | Type | Status | Score | Weight |
|---|---|---|---|---|
| sulfide_oxidation_active_in_amd_stream | pathway_active | satisfied | 1.00 | 1.5 |
| dissim_sulfate_reduction_anoxic_sediment | pathway_active | satisfied | 1.00 | 1.0 |
| nitrate_reduction_explored | pathway_active | satisfied | 1.00 | 0.5 |
| nitrogen_fixation_explored | pathway_active | satisfied | 1.00 | 0.3 |

## Evidence per claim

### sulfide_oxidation_active_in_amd_stream

- status=satisfied score=1.00 weight=1.5
- evidence: `{'pathway_id': 'Sulfide oxidation', 'element': 'sulfur', 'n_active_mags': 4, 'mean_completeness': 83.33, 'total_contribution': 333.33}`
- evaluator note: Sulfide oxidation: 4 个 MAG 活跃，平均完整度 83%，总贡献 333.3

### dissim_sulfate_reduction_anoxic_sediment

- status=satisfied score=1.00 weight=1.0
- evidence: `{'pathway_id': 'Dissim. sulfate red.', 'element': 'sulfur', 'n_active_mags': 2, 'mean_completeness': 75.0, 'total_contribution': 150.0}`
- evaluator note: Dissim. sulfate red.: 2 个 MAG 活跃，平均完整度 75%，总贡献 150.0

### nitrate_reduction_explored

- status=satisfied score=1.00 weight=0.5
- evidence: `{'pathway_id': 'Nitrate reduction', 'element': 'nitrogen', 'n_active_mags': 4, 'mean_completeness': 83.33, 'total_contribution': 333.33}`
- evaluator note: Nitrate reduction: 4 个 MAG 活跃，平均完整度 83%，总贡献 333.3

### nitrogen_fixation_explored

- status=satisfied score=1.00 weight=0.3
- evidence: `{'pathway_id': 'N fixation', 'element': 'nitrogen', 'n_active_mags': 4, 'mean_completeness': 100.0, 'total_contribution': 400.0}`
- evaluator note: N fixation: 4 个 MAG 活跃，平均完整度 100%，总贡献 400.0
