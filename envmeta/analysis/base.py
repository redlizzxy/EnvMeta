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
