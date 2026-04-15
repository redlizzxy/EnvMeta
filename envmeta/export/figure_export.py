"""图形导出：PNG（位图预览）/ PDF + SVG（矢量，SCI 投稿）/ TIFF（300+ dpi，
毕业论文 / 期刊位图要求）。

TIFF 用 matplotlib 的 Pillow 后端，默认走 LZW 压缩。SVG 是纯矢量，不带 dpi。
"""
from __future__ import annotations

import io
from pathlib import Path

import matplotlib.figure


SUPPORTED_FORMATS = ("png", "pdf", "svg", "tiff")
DEFAULT_DPI = {"png": 300, "pdf": 300, "svg": 300, "tiff": 600}
MIME_TYPES = {
    "png":  "image/png",
    "pdf":  "application/pdf",
    "svg":  "image/svg+xml",
    "tiff": "image/tiff",
}


def _normalize_fmt(fmt: str) -> str:
    fmt = fmt.lower()
    if fmt == "tif":
        fmt = "tiff"
    return fmt


def _savefig_kwargs(fmt: str) -> dict:
    """按格式拼 savefig 的额外 kwargs。SVG 不吃 dpi；TIFF 加 pil_kwargs 压缩。"""
    if fmt == "svg":
        return {}
    if fmt == "tiff":
        # Pillow 支持 LZW；保 alpha 通道会很大，关闭
        return {"pil_kwargs": {"compression": "tiff_lzw"}}
    return {}


def export_figure(
    fig: matplotlib.figure.Figure,
    out_path: str | Path,
    fmt: str = "png",
    dpi: int | None = None,
) -> Path:
    """保存 Figure 到磁盘。"""
    fmt = _normalize_fmt(fmt)
    if fmt not in SUPPORTED_FORMATS:
        raise ValueError(f"不支持的格式 {fmt!r}，支持：{SUPPORTED_FORMATS}")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    save_kw = _savefig_kwargs(fmt)
    if fmt != "svg":
        save_kw["dpi"] = dpi or DEFAULT_DPI[fmt]
    fig.savefig(out_path, format=fmt, bbox_inches="tight",
                facecolor="white", **save_kw)
    return out_path


def export_to_bytes(
    fig: matplotlib.figure.Figure,
    fmt: str = "png",
    dpi: int | None = None,
) -> bytes:
    """把 Figure 序列化为 bytes（供 Streamlit st.download_button 用）。"""
    fmt = _normalize_fmt(fmt)
    if fmt not in SUPPORTED_FORMATS:
        raise ValueError(f"不支持的格式 {fmt!r}，支持：{SUPPORTED_FORMATS}")
    buf = io.BytesIO()
    save_kw = _savefig_kwargs(fmt)
    if fmt != "svg":
        save_kw["dpi"] = dpi or DEFAULT_DPI[fmt]
    fig.savefig(buf, format=fmt, bbox_inches="tight",
                facecolor="white", **save_kw)
    buf.seek(0)
    return buf.getvalue()
