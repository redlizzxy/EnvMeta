# Wei2024_paddy_As-N_coupling — Hypothesis Score

- overall_score: **0.627**
- label: **INSUFFICIENT**
- n_satisfied / n_total / n_skipped: 3 / 5 / 0
- null_p: 0.9000 (n=999 permutations)
- weight_robust: True (OAT ±20% perturbation)
- veto_reasons: ['nitrate_reduction_active: unsatisfied (required=true)', 'as_n_coupling_arsenite_nitrate: partial (required=true)']

## Per-claim results

| ID | Type | Status | Score | Weight |
|---|---|---|---|---|
| arsenite_oxidation_active | pathway_active | satisfied | 1.00 | 1.5 |
| nitrate_reduction_active | pathway_active | unsatisfied | 0.00 | 1.5 |
| nitrite_reduction_active | pathway_active | satisfied | 1.00 | 1.0 |
| nitrous_oxide_reduction_active | pathway_active | satisfied | 1.00 | 0.7 |
| as_n_coupling_arsenite_nitrate | coupling_possible | partial | 0.50 | 2.0 |

## Evidence per claim

### arsenite_oxidation_active

- status=satisfied score=1.00 weight=1.5
- evidence: `{'pathway_id': 'Arsenite oxidation', 'element': 'arsenic', 'n_active_mags': 9, 'mean_completeness': 50.0, 'total_contribution': 0.18}`
- evaluator note: Arsenite oxidation: 9 个 MAG 活跃，平均完整度 50%，总贡献 0.2

### nitrate_reduction_active

- status=unsatisfied score=0.00 weight=1.5
- evidence: `{'pathway_id': 'Nitrate reduction', 'element': 'nitrogen', 'n_active_mags': 0, 'mean_completeness': 0.0, 'total_contribution': 0}`
- evaluator note: Nitrate reduction: 无活跃 MAG

### nitrite_reduction_active

- status=satisfied score=1.00 weight=1.0
- evidence: `{'pathway_id': 'Nitrite reduction', 'element': 'nitrogen', 'n_active_mags': 10, 'mean_completeness': 50.0, 'total_contribution': 0.23}`
- evaluator note: Nitrite reduction: 10 个 MAG 活跃，平均完整度 50%，总贡献 0.2

### nitrous_oxide_reduction_active

- status=satisfied score=1.00 weight=0.7
- evidence: `{'pathway_id': 'N2O reduction', 'element': 'nitrogen', 'n_active_mags': 52, 'mean_completeness': 100.0, 'total_contribution': 2.16}`
- evaluator note: N2O reduction: 52 个 MAG 活跃，平均完整度 100%，总贡献 2.2

### as_n_coupling_arsenite_nitrate

- status=partial score=0.50 weight=2.0
- evidence: `{'kb_product': 'As(V)', 'kb_type': 'redox', 'obs_a': True, 'obs_b': False}`
- evaluator note: As(III)↔NO3- 耦合：数据里缺 NO3-（仅一端观测到）
