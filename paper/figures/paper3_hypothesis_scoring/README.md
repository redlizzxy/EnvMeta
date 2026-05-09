# Paper 3 — Hypothesis Scoring Calibration + Stress Test Tables/Figures

> Paper 3 (EnvMeta methodology paper, target: iMeta) 的 Tables + Figures 实物素材
> for the controlled experiment showing EnvMeta hypothesis scoring works as both
> a calibration tool and a discrimination instrument across 4 KEGG-curated
> metagenomic datasets + 3 stress tests.

---

## 文件清单

### Tables (Markdown + TSV 双格式)

| 文件 | 角色 | 论文位置 | 数据来源 |
|---|---|---|---|
| [`table1_calibration.md`](table1_calibration.md) / `.tsv` | Four-Arm calibration 结果（5 行：作者 / Wei / Liu / Grettenberger / Ayala） | Results §X.1 / Methods §4.6.4 | [`scoring_validation_experiment_results.md`](../../manuscript/scoring_validation_experiment_results.md) §1 |
| [`table2_stress.md`](table2_stress.md) / `.tsv` | Three-Arm stress-test 结果（3 行：Liu / Grettenberger / Ayala） | Results §X.2 / Methods §4.6.5 | [`stress_test_results.md`](../../manuscript/stress_test_results.md) §0 + §3 + §4 |

### Figure X (PDF + PNG 600dpi + SVG 三格式)

| 文件 | 用途 |
|---|---|
| [`figure_x_calibration_vs_stress.pdf`](figure_x_calibration_vs_stress.pdf) | **publication-grade** vector PDF，可直接嵌入 paper / SI |
| [`figure_x_calibration_vs_stress.png`](figure_x_calibration_vs_stress.png) | 600 dpi raster，preview / supplementary slides 用 |
| [`figure_x_calibration_vs_stress.svg`](figure_x_calibration_vs_stress.svg) | editable vector，可在 Inkscape / Illustrator 中调整 |
| [`figure_x_generate.py`](figure_x_generate.py) | matplotlib 复现脚本（可重跑产出三格式）|

### 复现命令

```bash
conda activate envmeta
python paper/figures/paper3_hypothesis_scoring/figure_x_generate.py
```

---

## Figure X 内容速览

横向条形图（horizontal bar plot），三对（dataset × {calibration / stress}）：

| Dataset | Calibration | Stress | Discrimination 标注 |
|---|---|---|---|
| Liu 2023 cold seep (same-topic) | 1.000 ▆▆▆▆ | 0.625 ▆▆ | B-tier (binary-threshold limit) |
| **Grettenberger 2021 AMD stream (cross-topic) ★** | 1.000 ▆▆▆▆ | **0.250 ▆** | **A-tier cross-topic n=0 ★** |
| **Ayala 2020 pit lake (cross-topic) ★** | 1.000 ▆▆▆▆ | 0.455 ▆▆ | **B-tier + cross-topic n=0 ★** |

虚线标注 EnvMeta 默认阈值：strong=0.75 / suggestive=0.40。

★ 标注的两行 = 无砷数据集，cross-topic `arsenate_reduction_should_dominate` 实测 n=0 active MAGs → **EnvMeta 领域中立性铁证**。

---

## 配色

| 颜色 | 含义 |
|---|---|
| `#2E86AB` 深蓝 | Calibration overall_score |
| `#E63946` 红 | Stress overall_score |
| `#888888` 灰虚线 | EnvMeta default thresholds (0.40 / 0.75) |

字体：DejaVu Sans（开源 + Linux/Mac/Windows 通用 + ★ 字符兼容）。

---

## 论文 caption 候选（Figure X，~80 字）

> **Figure X. Calibration vs stress overall_score gap, three KEGG-curated
> datasets.** Each pair shows calibration `STRONG` (1.000) and the
> corresponding stress-test overall_score on three pre-registered hypothesis
> YAML pairs. Vertical dashed/dotted lines = EnvMeta default thresholds
> (`strong=0.75`, `suggestive=0.40`). The cross-topic `arsenate_reduction
> should_dominate` claim was rejected with n=0 active MAGs in 2/2 non-arsenic
> datasets (★ Grettenberger AMD stream and Ayala pit lake), providing direct
> evidence of EnvMeta's domain-neutral scoring engine. B-tier annotations
> indicate Liu and Ayala stress tests where the binary `mean_completeness ≥
> 50%` threshold cannot distinguish "dominant" from "detectable but weak"
> pathway activity (resolved by future `dominance_score` field).

---

## 维护记录

| 日期 | 事项 |
|---|---|
| 2026-05-09 | 初版 — Table 1 / Table 2 / Figure X 三件套，Vancouver 引用对齐 [`methods_external_validation.md`](../../manuscript/methods_external_validation.md) |
