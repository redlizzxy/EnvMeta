"""新手落地包 — 静态数据字典（解读文案 / 研究问题向导 / 文件反向索引）。

数据在三个纯 dict 模块里集中维护；app.py 引用它们渲染 UI。新增分析时只需
在对应字典加一行。
"""
from __future__ import annotations

from envmeta.help.file_analysis_map import ANALYSIS_INPUTS, FILE_TO_ANALYSIS
from envmeta.help.interpretations import INTERPRETATIONS
from envmeta.help.research_navigator import NAVIGATOR

__all__ = ["NAVIGATOR", "INTERPRETATIONS", "FILE_TO_ANALYSIS", "ANALYSIS_INPUTS"]
