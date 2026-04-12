#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
论文绘图统一配置文件 (plot_config.py)
================================================================================
所有Python绘图脚本统一导入此文件，确保配色、字体、尺寸、路径一致。

使用方法：
    import sys; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from config.plot_config import *
================================================================================
"""

import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
from scipy import stats

# ============================================================
# 1. 路径配置
# ============================================================
#
# 项目位置:  /home/lizzxy/thesis_project/      (WSL ext4, 快速IO)
# 数据源:    /home/lizzxy/meta2/               (通过软链接引用, 不复制)
# 输出目标:  /mnt/d/workdata/                  (sync_to_win.sh 同步)
#

def get_project_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PROJECT_ROOT = get_project_root()

# 数据源（meta2 原始分析结果，通过 data/raw/ 下的软链接访问）
META2_ROOT  = os.path.expanduser("~/meta2")
DATA_RAW    = os.path.join(PROJECT_ROOT, "data", "raw")       # 软链接到 meta2 子目录
DATA_PROC   = os.path.join(PROJECT_ROOT, "data", "processed") # 脚本中间产物

# 图表输出
FIG_THESIS  = os.path.join(PROJECT_ROOT, "figures", "thesis_cn")
FIG_PAPER   = os.path.join(PROJECT_ROOT, "figures", "paper_en")
FIG_SUPPL   = os.path.join(PROJECT_ROOT, "figures", "supplementary")
TABLES_DIR  = os.path.join(PROJECT_ROOT, "tables")
LOGS_DIR    = os.path.join(PROJECT_ROOT, "logs")

# Windows 输出路径 (D:\workdata)
WIN_TARGET = os.environ.get("THESIS_WIN_TARGET", "/mnt/d/workdata")

for d in [DATA_RAW, DATA_PROC, FIG_THESIS, FIG_PAPER, FIG_SUPPL, TABLES_DIR, LOGS_DIR]:
    os.makedirs(d, exist_ok=True)


# ============================================================
# 2. 样本与分组
# ============================================================

GROUPS = {
    "CK": ["2_1", "2_5", "2_7"],
    "A":  ["8_1", "8_2", "8_3", "8_4"],
    "B":  ["9_1", "9_2", "9_3"],
}
GROUP_ORDER = ["CK", "A", "B"]
ALL_SAMPLES = [s for g in GROUP_ORDER for s in GROUPS[g]]
SAMPLE_TO_GROUP = {s: g for g, samples in GROUPS.items() for s in samples}


# ============================================================
# 3. 配色方案（最终确认版 2026-03-12）
# ============================================================
#
# 柱状图 / PCoA / 元素循环 → 方案B (Earth-Warm)
# 堆叠图 → 方案A (Ocean-Forest)，占比高用浅色(7-12)、占比低用深色(1-6)
# 热图 → 方案D色阶 (RdYlBu_r)
#

# ---- 三组核心配色 (方案B, 用于柱状图/PCoA/元素循环) ----
GROUP_COLORS = {"CK": "#1c9cbd", "A": "#e3943d", "B": "#92181e"}
GROUP_COLORS_LIGHT = {"CK": "#d4e5ef", "A": "#fbd8b5", "B": "#d0916e"}

# ---- 元素循环配色 (与分组色系呼应) ----
ELEMENT_COLORS = {
    "arsenic": "#92181e",   # As → 深红 (同B组)
    "sulfur":  "#e3943d",   # S  → 橙金 (同A组)
    "iron":    "#a35626",   # Fe → 深棕
    "nitrogen":"#1c9cbd",   # N  → 青蓝 (同CK组)
}

# ---- 堆叠图12色 (方案A, 重排: 占比高→浅色, 占比低→深色) ----
# Top1-6 (高丰度物种) → 柔和浅色，避免大面积深色压迫
# Top7-12 (低丰度物种) → 饱和深色，小面积也清晰可辨
STACK_12 = [
    "#bbc4e4", "#d7e7af", "#bf9d6d", "#a4d6c1", "#bce2e8", "#b5b5b6",  # 浅色 (高丰度)
    "#0084c2", "#3ab483", "#00978c", "#21825e", "#6aa3d2", "#3d62ad",   # 深色 (低丰度)
]

# 兼容旧变量名
PALETTE_12 = STACK_12

# ---- 门水平配色 (基于堆叠图色系) ----
PHYLUM_COLORS = {
    "Pseudomonadota": "#bbc4e4", "Chloroflexota": "#0084c2",
    "Acidobacteriota": "#3ab483", "Actinomycetota": "#bf9d6d",
    "Bacteroidota": "#d7e7af", "Gemmatimonadota": "#a4d6c1",
    "Desulfobacterota": "#00978c", "Planctomycetota": "#6aa3d2",
    "Verrucomicrobiota": "#bce2e8", "Nitrospirota": "#21825e",
    "Firmicutes": "#3d62ad", "Other": "#b5b5b6",
}

# ---- 热图色阶 (方案D) ----
HEATMAP_CMAP = "RdYlBu_r"


# ============================================================
# 4. 字体配置
# ============================================================
#
# SCI小论文：Arial / Helvetica (sans-serif)
# 毕业论文：宋体(中文正文) + Times New Roman(英文/数字)
# WSL替代：Liberation Sans ≡ Arial, Liberation Serif ≡ TNR
#

def setup_fonts():
    # --- SCI: Arial ---
    en_sans = "Arial"
    sans_candidates = [
        "/mnt/c/Windows/Fonts/arial.ttf",
        "/mnt/c/Windows/Fonts/arialbd.ttf",
        "/mnt/c/Windows/Fonts/ariali.ttf",
    ]
    for fp in sans_candidates:
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    # fallback: Liberation Sans ≡ Arial
    en_sans_fallback = "Liberation Sans"

    # --- Thesis digits: Times New Roman ---
    en_serif = "Times New Roman"
    serif_candidates = [
        "/mnt/c/Windows/Fonts/times.ttf",
        "/mnt/c/Windows/Fonts/timesbd.ttf",
        "/mnt/c/Windows/Fonts/timesi.ttf",
    ]
    for fp in serif_candidates:
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    en_serif_fallback = "Liberation Serif"

    # --- Thesis CN: SimSun ---
    cn_font = None
    cn_candidates = [
        "/mnt/c/Windows/Fonts/simsun.ttc",
        "/mnt/c/Windows/Fonts/simsunb.ttf",
        "/usr/share/fonts/truetype/wqy/wqy-microhei.ttc",
        "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
    ]
    for fp in cn_candidates:
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
            cn_font = fm.FontProperties(fname=fp).get_name()
            break

    # 检测实际可用字体
    available = {f.name for f in fm.fontManager.ttflist}
    final_sans  = en_sans if en_sans in available else en_sans_fallback
    final_serif = en_serif if en_serif in available else en_serif_fallback

    return final_sans, final_serif, cn_font

EN_FONT_SANS, EN_FONT_SERIF, CN_FONT = setup_fonts()

# 向后兼容
EN_FONT = EN_FONT_SANS


# ============================================================
# 5. Matplotlib 主题
# ============================================================

_BASE_RC = {
    "axes.linewidth": 0.8, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "lines.linewidth": 1.2, "axes.facecolor": "white", "figure.facecolor": "white",
    "savefig.facecolor": "white", "axes.unicode_minus": False,
    "figure.constrained_layout.use": True,
}

def apply_thesis_style():
    """毕业论文风格：中文宋体 + 英文/数字 Times New Roman (serif)"""
    serif_list = [CN_FONT, EN_FONT_SERIF, "DejaVu Serif"] if CN_FONT else [EN_FONT_SERIF, "DejaVu Serif"]
    rc = {**_BASE_RC,
        "font.family": "serif", "font.serif": serif_list,
        "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10,
        "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8,
        "figure.dpi": 150, "savefig.dpi": 300,
        "mathtext.fontset": "stix",
    }
    plt.rcParams.update(rc)

def apply_paper_style():
    """SCI小论文风格：Arial (sans-serif)"""
    rc = {**_BASE_RC,
        "font.family": "sans-serif",
        "font.sans-serif": [EN_FONT_SANS, "Liberation Sans", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "figure.dpi": 150, "savefig.dpi": 600,
        "mathtext.fontset": "dejavusans",
    }
    plt.rcParams.update(rc)


# ============================================================
# 6. 图片尺寸（英寸）
# ============================================================

# SCI
FIG_SINGLE_COL = (3.54, 2.76)   # 90×70 mm
FIG_1_5_COL    = (5.51, 3.94)   # 140×100 mm
FIG_DOUBLE_COL = (7.09, 4.72)   # 180×120 mm
FIG_FULL_PAGE  = (7.09, 9.45)   # 180×240 mm

# 毕业论文
FIG_THESIS_STD  = (6.0, 4.5)
FIG_THESIS_WIDE = (7.5, 4.5)
FIG_THESIS_TALL = (6.0, 7.0)
FIG_THESIS_FULL = (7.5, 9.0)


# ============================================================
# 7. 双版本保存
# ============================================================

def save_figure(fig, fig_id, desc_cn="", desc_en="",
                thesis_size=None, paper_size=None,
                formats_thesis=("pdf",), formats_paper=("pdf", "tiff")):
    """同时保存毕业论文版（中文）和小论文版（英文）"""
    if thesis_size:
        fig.set_size_inches(thesis_size)
    for fmt in formats_thesis:
        fp = os.path.join(FIG_THESIS, f"{fig_id}_{desc_cn}.{fmt}")
        fig.savefig(fp, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
        print(f"  [论文] {fp}")

    if paper_size:
        fig.set_size_inches(paper_size)
    for fmt in formats_paper:
        fp = os.path.join(FIG_PAPER, f"{fig_id}_{desc_en}.{fmt}")
        dpi = 600 if fmt == "tiff" else 300
        fig.savefig(fp, dpi=dpi, bbox_inches="tight", facecolor="white", edgecolor="none")
        print(f"  [Paper] {fp}")


# ============================================================
# 8. 统计工具
# ============================================================

def pairwise_wilcoxon(data_dict, pairs=None):
    if pairs is None:
        groups = list(data_dict.keys())
        pairs = [(groups[i], groups[j]) for i in range(len(groups)) for j in range(i+1, len(groups))]
    results = {}
    for g1, g2 in pairs:
        stat, p = stats.mannwhitneyu(data_dict[g1], data_dict[g2], alternative="two-sided")
        results[(g1, g2)] = (stat, p)
    return results

def sig_label(p):
    if p < 0.001: return "***"
    elif p < 0.01: return "**"
    elif p < 0.05: return "*"
    else: return "ns"
