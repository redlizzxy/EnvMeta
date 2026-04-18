"""分析模块通用接口。

所有 `envmeta/analysis/*.py` 下的分析器都返回 AnalysisResult。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import matplotlib.figure
import pandas as pd


@dataclass
class AnalysisResult:
    figure: matplotlib.figure.Figure
    stats: pd.DataFrame | None = None
    params: dict = field(default_factory=dict)
    # 可选的"原始领域对象"—— 分析器若有比 DataFrame 更丰富的结构产出
    # （如 cycle_diagram 的 CycleData、mag_heatmap 的聚类信息），可塞这里。
    # 让下游（导出中心 / HTML 交互导出）不必重跑 analyze 就能访问 typed 结构。
    data: object | None = None
