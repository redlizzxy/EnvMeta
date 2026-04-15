"""figure_export 多格式支持测试（S2.5-5）。"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from envmeta.export.figure_export import (
    SUPPORTED_FORMATS, export_figure, export_to_bytes,
)


@pytest.fixture
def small_fig():
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.plot([0, 1, 2], [0, 1, 4])
    ax.set_title("test")
    yield fig
    plt.close(fig)


def test_supported_formats_includes_all_four():
    assert set(SUPPORTED_FORMATS) == {"png", "pdf", "svg", "tiff"}


def test_png_bytes_have_png_header(small_fig):
    data = export_to_bytes(small_fig, "png")
    assert data[:8] == b"\x89PNG\r\n\x1a\n"


def test_pdf_bytes_start_with_pdf_magic(small_fig):
    data = export_to_bytes(small_fig, "pdf")
    assert data.startswith(b"%PDF")


def test_svg_bytes_are_xml_svg(small_fig):
    data = export_to_bytes(small_fig, "svg")
    head = data[:256].decode("utf-8", errors="ignore")
    assert "<svg" in head or "<?xml" in head


def test_tiff_bytes_have_tiff_magic(small_fig):
    data = export_to_bytes(small_fig, "tiff")
    # TIFF magic: II*\x00 (little-endian) or MM\x00* (big-endian)
    assert data[:4] in (b"II*\x00", b"MM\x00*")


def test_tif_alias_accepted(small_fig):
    """tif 应该归一化到 tiff。"""
    data = export_to_bytes(small_fig, "tif")
    assert data[:4] in (b"II*\x00", b"MM\x00*")


def test_unsupported_format_raises(small_fig):
    with pytest.raises(ValueError):
        export_to_bytes(small_fig, "bmp")


def test_export_figure_writes_files(tmp_path, small_fig):
    for fmt in ("png", "pdf", "svg", "tiff"):
        p = export_figure(small_fig, tmp_path / f"x.{fmt}", fmt=fmt)
        assert p.exists()
        assert p.stat().st_size > 0


def test_svg_ignores_dpi(small_fig):
    """SVG 是纯矢量，传 dpi 不应报错（应被忽略）。"""
    data = export_to_bytes(small_fig, "svg", dpi=1200)
    head = data[:256].decode("utf-8", errors="ignore")
    assert "<svg" in head or "<?xml" in head
