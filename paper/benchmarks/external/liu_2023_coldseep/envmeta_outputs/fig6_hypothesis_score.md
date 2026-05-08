# Liu2023_coldseep_pre-data_hypothesis — Hypothesis Score

- overall_score: **1.000**
- label: **STRONG**
- n_satisfied / n_total / n_skipped: 4 / 4 / 0
- null_p: n/a
- weight_robust: True (OAT ±20%)
- veto_reasons: none

## Per-claim results

| ID | Type | Status | Score | Weight |
|---|---|---|---|---|
| arsenate_reduction_active_in_anoxic_seep | pathway_active | satisfied | 1.00 | 1.5 |
| as_transport_detox_active | pathway_active | satisfied | 1.00 | 1.0 |
| as_methylation_explored | pathway_active | satisfied | 1.00 | 0.7 |
| respiratory_as_reduction_explored | pathway_active | satisfied | 1.00 | 0.5 |

## Evidence per claim

### arsenate_reduction_active_in_anoxic_seep

- status=satisfied score=1.00 weight=1.5
- evidence: `{'pathway_id': 'Arsenate reduction', 'element': 'arsenic', 'n_active_mags': 18, 'mean_completeness': 50.0, 'total_contribution': 6.43}`
- evaluator note: Arsenate reduction: 18 个 MAG 活跃，平均完整度 50%，总贡献 6.4

### as_transport_detox_active

- status=satisfied score=1.00 weight=1.0
- evidence: `{'pathway_id': 'As transport/detox', 'element': 'arsenic', 'n_active_mags': 1, 'mean_completeness': 50.0, 'total_contribution': 0.26}`
- evaluator note: As transport/detox: 1 个 MAG 活跃，平均完整度 50%，总贡献 0.3

### as_methylation_explored

- status=satisfied score=1.00 weight=0.7
- evidence: `{'pathway_id': 'As methylation', 'element': 'arsenic', 'n_active_mags': 572, 'mean_completeness': 50.0, 'total_contribution': 600.68}`
- evaluator note: As methylation: 572 个 MAG 活跃，平均完整度 50%，总贡献 600.7

### respiratory_as_reduction_explored

- status=satisfied score=1.00 weight=0.5
- evidence: `{'pathway_id': 'Resp. arsenate red.', 'element': 'arsenic', 'n_active_mags': 17, 'mean_completeness': 50.0, 'total_contribution': 10.99}`
- evaluator note: Resp. arsenate red.: 17 个 MAG 活跃，平均完整度 50%，总贡献 11.0
