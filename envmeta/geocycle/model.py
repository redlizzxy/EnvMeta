"""循环图的轻量数据模型。

这些 dataclass 是 inference 和 renderer 之间的契约，也是 Phase 4 D3 可视化
的 JSON 导出来源。
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class MAGContribution:
    """一个 MAG 对某条通路的贡献。"""
    mag: str
    phylum: str
    completeness: float   # 0-100
    abundance_mean: float  # MAG 在所有样本的平均丰度
    label: str | None = None  # 展示名（优先 Genus；兜底 MAG id）


@dataclass
class PathwayActivity:
    """一条通路在数据集中的活跃度。"""
    pathway_id: str
    element: str
    display_name: str
    reaction: str | None
    n_active_mags: int           # completeness ≥ threshold 的 MAG 数
    mean_completeness: float     # 活跃 MAG 的平均完整度
    total_contribution: float    # Σ completeness × abundance_mean
    contributors: list[MAGContribution] = field(default_factory=list)  # 按贡献降序


@dataclass
class EnvCorrelation:
    """环境因子与通路活性（KO 丰度均值）的相关性。"""
    pathway_id: str
    env_factor: str
    rho: float            # Spearman
    p_value: float
    n_samples: int


@dataclass
class ElementCycle:
    """一个元素（As/N/S/Fe/...）下的所有通路活性。"""
    element_id: str
    display_name: str
    color: str
    pathways: list[PathwayActivity] = field(default_factory=list)

    @property
    def n_active_pathways(self) -> int:
        return sum(1 for p in self.pathways if p.n_active_mags > 0)

    @property
    def total_contribution(self) -> float:
        return sum(p.total_contribution for p in self.pathways)


@dataclass
class SensitivityRow:
    """一条通路在多档完整度阈值下的 Top-1 contributor 稳定性。"""
    pathway_id: str
    element: str
    display_name: str
    thresholds: list[float]
    top1_by_threshold: list[str | None]
    n_active_by_threshold: list[int]
    robust: bool   # True = 三档 Top-1 一致


@dataclass
class CycleData:
    """一次推断的完整结果（供 renderer 消费 + JSON 导出）。"""
    elements: list[ElementCycle] = field(default_factory=list)
    env_correlations: list[EnvCorrelation] = field(default_factory=list)
    full_corr_matrix: list[EnvCorrelation] = field(default_factory=list)  # 不过滤
    sensitivity: list[SensitivityRow] = field(default_factory=list)
    params: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)   # 样本数、MAG 数等汇总

    def to_flat_stats(self):
        """扁平化为长表 DataFrame（供 AnalysisResult.stats 使用）。"""
        import pandas as pd
        rows = []
        for el in self.elements:
            for pw in el.pathways:
                rows.append({
                    "type": "pathway",
                    "element": el.element_id,
                    "pathway_id": pw.pathway_id,
                    "display_name": pw.display_name,
                    "n_active_mags": pw.n_active_mags,
                    "mean_completeness": pw.mean_completeness,
                    "total_contribution": pw.total_contribution,
                    "top_mag": pw.contributors[0].label if pw.contributors else None,
                })
        for ec in self.env_correlations:
            rows.append({
                "type": "env_correlation",
                "pathway_id": ec.pathway_id,
                "display_name": ec.env_factor,
                "mean_completeness": ec.rho,
                "total_contribution": ec.p_value,
            })
        for ec in self.full_corr_matrix:
            rows.append({
                "type": "full_correlation",
                "pathway_id": ec.pathway_id,
                "display_name": ec.env_factor,
                "mean_completeness": ec.rho,
                "total_contribution": ec.p_value,
            })
        for sr in self.sensitivity:
            rows.append({
                "type": "sensitivity",
                "element": sr.element,
                "pathway_id": sr.pathway_id,
                "display_name": sr.display_name,
                "top_mag": "|".join(str(x) for x in sr.top1_by_threshold),
                "n_active_mags": sum(sr.n_active_by_threshold),
                "mean_completeness": 1.0 if sr.robust else 0.0,  # robust flag
            })
        return pd.DataFrame(rows)
