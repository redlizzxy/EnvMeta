"""通用参数控件：所有分析页面共享的宽/高/字号/DPI 选择。

使用示例：
    from envmeta.params.common import render_figure_size, render_font_controls
    size = render_figure_size({"width_mm": 160, "height_mm": 100}, prefix="stackplot")
    fonts = render_font_controls({}, prefix="stackplot")
    params = {**size, **fonts, ...}
"""
from __future__ import annotations

import streamlit as st


def render_figure_size(defaults: dict | None = None, *, prefix: str = "") -> dict:
    """在侧边栏渲染宽/高滑块，返回 {"width_mm": int, "height_mm": int}。"""
    d = defaults or {}
    st.sidebar.markdown("**画布**")
    w = st.sidebar.slider(
        "图宽 (mm)", 80, 300, d.get("width_mm", 160), step=10,
        key=f"{prefix}_width_mm",
    )
    h = st.sidebar.slider(
        "图高 (mm)", 60, 250, d.get("height_mm", 100), step=10,
        key=f"{prefix}_height_mm",
    )
    return {"width_mm": w, "height_mm": h}


def render_font_controls(defaults: dict | None = None, *, prefix: str = "") -> dict:
    """渲染字号滑块，返回 {"title_size": int, "label_size": int, "tick_size": int}。"""
    d = defaults or {}
    with st.sidebar.expander("字号", expanded=False):
        title = st.slider("标题", 6, 20, d.get("title_size", 11), key=f"{prefix}_title_size")
        label = st.slider("轴标签", 6, 16, d.get("label_size", 9), key=f"{prefix}_label_size")
        tick = st.slider("刻度", 5, 14, d.get("tick_size", 8), key=f"{prefix}_tick_size")
    return {"title_size": title, "label_size": label, "tick_size": tick}


def render_dpi_selector(default: int = 300, *, prefix: str = "") -> int:
    """DPI 下拉：150 / 300 / 600。"""
    options = [150, 300, 600]
    return st.sidebar.selectbox(
        "导出 DPI", options,
        index=options.index(default),
        key=f"{prefix}_dpi",
    )
