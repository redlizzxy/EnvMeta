"""图形导出。迭代 1 只支持 PNG / PDF，SVG / TIFF 留到迭代 2。"""
from __future__ import annotations

import io
from pathlib import Path

import matplotlib.figure


SUPPORTED_FORMATS = ("png", "pdf")
DEFAULT_DPI = {"png": 300, "pdf": 300}


def export_figure(
    fig: matplotlib.figure.Figure,
    out_path: str | Path,
    fmt: str = "png",
    dpi: int | None = None,
) -> Path:
    """保存 Figure 到磁盘。"""
    fmt = fmt.lower()
    if fmt not in SUPPORTED_FORMATS:
        raise ValueError(f"不支持的格式 {fmt!r}，支持：{SUPPORTED_FORMATS}")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format=fmt, dpi=dpi or DEFAULT_DPI[fmt],
                bbox_inches="tight", facecolor="white")
    return out_path


def export_to_bytes(
    fig: matplotlib.figure.Figure,
    fmt: str = "png",
    dpi: int | None = None,
) -> bytes:
    """把 Figure 序列化为 bytes（供 Streamlit st.download_button 用）。"""
    fmt = fmt.lower()
    if fmt not in SUPPORTED_FORMATS:
        raise ValueError(f"不支持的格式 {fmt!r}，支持：{SUPPORTED_FORMATS}")
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=dpi or DEFAULT_DPI[fmt],
                bbox_inches="tight", facecolor="white")
    buf.seek(0)
    return buf.getvalue()
