#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_physicochemical.py — 理化指标柱状图
==============================================================================
输入: ordata.xlsx
输出:
  单独图: Fig1-2a_pH, Fig1-2b_Eh, Fig1-2c_TOC, Fig1-2d_totalAs, Fig1-2e_As形态
  合并图: Fig1-2_physicochemical (2行3列, 第6格留空)

修改记录:
  v2 - 每个指标单独出图 + 合并图
     - 修复Eh显著性标注溢出
     - As形态图例移至右上角外侧
     - 删除As形态浓度图(原f)
==============================================================================
"""

import sys, os
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

# ---- 加载配置 ----
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
try:
    from config.plot_config import *
except ImportError:
    print("[WARN] config not found, using defaults")
    GROUP_COLORS = {"CK": "#1c9cbd", "A": "#e3943d", "B": "#92181e"}
    FIG_THESIS = FIG_PAPER = "."
    DATA_RAW = "."

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.ticker as mticker
import matplotlib.font_manager as fm

# ---- 字体 ----
F_SANS  = "Liberation Sans"
F_SERIF = "Liberation Serif"
F_CN    = "WenQuanYi Micro Hei"

avail = {f.name for f in fm.fontManager.ttflist}
if "Arial" in avail: F_SANS = "Arial"
if "Times New Roman" in avail: F_SERIF = "Times New Roman"
for cn_try in ["SimSun", "WenQuanYi Micro Hei", "WenQuanYi Zen Hei"]:
    if cn_try in avail: F_CN = cn_try; break


# ============================================================
# 1. 数据读取 (同上一版，不变)
# ============================================================

def load_data(xlsx_path):
    xls = pd.ExcelFile(xlsx_path)
    data = {}

    # pH
    df_ph = pd.read_excel(xls, sheet_name='pH')
    data['pH'] = {}
    for _, row in df_ph.iterrows():
        g = str(row.iloc[0]).strip()
        data['pH'][g] = [float(row.iloc[i]) for i in [1,2,3]]

    # Eh
    df_eh = pd.read_excel(xls, sheet_name='Eh')
    data['Eh'] = {}
    for _, row in df_eh.iterrows():
        g = str(row.iloc[0]).strip()
        data['Eh'][g] = [float(row.iloc[i]) for i in [1,2,3]]

    # TOC
    df_toc = pd.read_excel(xls, sheet_name='TOC')
    data['TOC'] = {}
    for _, row in df_toc.iterrows():
        g = str(row.iloc[0]).strip()
        if g.startswith('0'): g = '砷渣'
        data['TOC'][g] = float(row.iloc[1])

    # 总As
    df_as = pd.read_excel(xls, sheet_name='总As')
    samples = df_as[df_as['类型'] == 'Sample'].copy()
    samples = samples[samples['样品名称'].str.match(r'^\d', na=False)]
    samples['ppb'] = pd.to_numeric(samples['浓度 [ ppb ]'], errors='coerce')
    dil = {'0': 300, '2': 800, '8': 600, '9': 10000}
    gmap = {'0': '砷渣', '2': 'CK', '8': 'A', '9': 'B'}
    data['total_As'] = {}
    for pfx, gname in gmap.items():
        sub = samples[samples['样品名称'].str.startswith(pfx + '-') |
                      samples['样品名称'].str.startswith(pfx + '+')]
        data['total_As'][gname] = (sub['ppb'] * dil[pfx] / 1000).values.tolist()

    # As形态
    df_form = pd.read_excel(xls, sheet_name='As形态')
    std_c, std_a = {'As3':[],'As5':[],'MMA':[]}, {'As3':[],'As5':[],'MMA':[]}
    ci_map = {'As3': 12, 'As5': 13, 'MMA': 14}
    ai_map = {'As3': 2, 'As5': 8, 'MMA': 5}
    for r in range(19, 25):
        for sp, ci in ci_map.items():
            std_c[sp].append(float(df_form.iloc[r, ci]))
    for r in range(1, 7):
        for sp, ai in ai_map.items():
            std_a[sp].append(float(df_form.iloc[r, ai]))
    cal = {}
    for sp in ['As3','As5','MMA']:
        s, i, _, _, _ = sp_stats.linregress(std_a[sp], std_c[sp])
        cal[sp] = (s, i)

    srows = [('0-1',9),('0+2',10),('0-3',11),('2-1',14),('2-2',15),('2-3',16),
             ('8-1',17),('8-2',18),('8-3',19),('9-1',20),('9-2',21),('9-3',22)]
    records = []
    for sn, ri in srows:
        row = df_form.iloc[ri]
        pfx, g, d = sn[0], gmap[sn[0]], dil[sn[0]]
        r_dict = {'sample': sn, 'group': g}
        for sp, ai in ai_map.items():
            r_dict[sp] = max((cal[sp][0] * float(row.iloc[ai]) + cal[sp][1]) * d, 0)
        r_dict['total'] = r_dict['As3'] + r_dict['As5'] + r_dict['MMA']
        records.append(r_dict)
    data['As_spec'] = pd.DataFrame(records)

    return data


# ============================================================
# 2. 工具函数
# ============================================================

def sig_star(p):
    if pd.isna(p): return ""
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

def calc_sig(vals_dict, groups):
    pairs = [(groups[i], groups[j]) for i in range(len(groups)) for j in range(i+1, len(groups))]
    results = {}
    for g1, g2 in pairs:
        v1, v2 = vals_dict.get(g1,[]), vals_dict.get(g2,[])
        if len(v1) >= 2 and len(v2) >= 2:
            try: _, p = sp_stats.mannwhitneyu(v1, v2, alternative='two-sided')
            except: p = 1.0
        else: p = np.nan
        results[(g1, g2)] = p
    return results

def add_sig_bracket(ax, x1, x2, y, text, fontfamily, h_frac=0.025, lw=0.8):
    """绘制显著性标注，h_frac 控制竖线高度占 ylim 的比例"""
    ylo, yhi = ax.get_ylim()
    dh = (yhi - ylo) * h_frac
    ax.plot([x1, x1, x2, x2], [y, y+dh, y+dh, y], 'k-', lw=lw, clip_on=False)
    ax.text((x1+x2)/2, y+dh*1.2, text, ha='center', va='bottom',
            fontsize=7, fontfamily=fontfamily, clip_on=False)

def style_ax(ax, ylabel, title_letter, font_main, font_num):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(f"({title_letter})", fontsize=10, fontfamily=font_main,
                 fontweight='bold', loc='left')
    ax.set_ylabel(ylabel, fontsize=9, fontfamily=font_main)
    for t in ax.get_yticklabels():
        t.set_fontfamily(font_num); t.set_fontsize(8)
    ax.tick_params(axis='both', length=3, width=0.6)


# ============================================================
# 3. 单图绘制函数
# ============================================================

GC  = {"CK": "#1c9cbd", "A": "#e3943d", "B": "#92181e"}
GC4 = {"砷渣": "#888888", **GC}
SPEC_COLORS = {'As(III)': '#92181e', 'As(V)': '#1c9cbd', 'MMA': '#e3943d'}

def get_fonts(lang):
    is_cn = (lang == 'cn')
    fm_ = F_CN if is_cn else F_SANS
    fn_ = F_SERIF if is_cn else F_SANS
    return is_cn, fm_, fn_

def get_labels(lang, n=3):
    is_cn = (lang == 'cn')
    if n == 3:
        return ["CK组","A组","B组"] if is_cn else ["CK","A","B"]
    else:
        return ["砷渣","CK组","A组","B组"] if is_cn else ["Slag","CK","A","B"]


def draw_bar3(data_dict, groups, ylabel, ylim, lang, title_letter,
              figsize=(3.5, 3.2)):
    """绘制3组柱状图（pH/Eh），带误差棒和显著性"""
    is_cn, fm_, fn_ = get_fonts(lang)
    fig, ax = plt.subplots(figsize=figsize)
    x = np.arange(len(groups))
    means = [np.mean(data_dict[g]) for g in groups]
    sds   = [np.std(data_dict[g], ddof=1) for g in groups]

    ax.bar(x, means, yerr=sds, color=[GC[g] for g in groups],
           edgecolor='#333', linewidth=0.6, capsize=4, width=0.55,
           error_kw={'linewidth': 0.8, 'capthick': 0.8})
    ax.set_xticks(x)
    ax.set_xticklabels(get_labels(lang, 3), fontsize=9, fontfamily=fm_)
    ax.set_ylim(ylim)
    style_ax(ax, ylabel, title_letter, fm_, fn_)

    # 显著性：从最内层pair开始画，逐层向上
    sigs = calc_sig(data_dict, groups)
    yrange = ylim[1] - ylim[0]
    y_base = max(m + s for m, s in zip(means, sds)) + yrange * 0.03
    step = yrange * 0.08  # 每层间距

    # 排序：先画窄的(相邻pair)，再画宽的(跨组pair)
    sorted_pairs = sorted(sigs.keys(), key=lambda p: abs(groups.index(p[1]) - groups.index(p[0])))
    for i, pair in enumerate(sorted_pairs):
        star = sig_star(sigs[pair])
        if star:
            g1i, g2i = groups.index(pair[0]), groups.index(pair[1])
            add_sig_bracket(ax, g1i, g2i, y_base + i * step, star, fn_)

    # 确保ylim上方留足空间给显著性标注
    needed_top = y_base + len(sorted_pairs) * step + yrange * 0.06
    if needed_top > ylim[1]:
        ax.set_ylim(ylim[0], needed_top)

    fig.tight_layout()
    return fig


def draw_toc(data, lang, figsize=(3.8, 3.2)):
    """TOC柱状图（4组，无误差棒）"""
    is_cn, fm_, fn_ = get_fonts(lang)
    groups = ['砷渣', 'CK', 'A', 'B']
    fig, ax = plt.subplots(figsize=figsize)
    x = np.arange(4)
    vals = [data['TOC'].get(g, 0) for g in groups]

    ax.bar(x, vals, color=[GC4[g] for g in groups],
           edgecolor='#333', linewidth=0.6, width=0.55)
    ax.set_xticks(x)
    ax.set_xticklabels(get_labels(lang, 4), fontsize=9, fontfamily=fm_)
    ax.set_ylim(0, 13)
    style_ax(ax, "TOC (g/kg)" if not is_cn else "TOC（g/kg）", "c", fm_, fn_)

    for i, v in enumerate(vals):
        ax.text(i, v + 0.3, f"{v:.1f}", ha='center', fontsize=7,
                fontfamily=fn_, color='#444')

    fig.tight_layout()
    return fig


def draw_total_as(data, lang, figsize=(3.8, 3.2)):
    """总As柱状图（4组，log Y轴）"""
    is_cn, fm_, fn_ = get_fonts(lang)
    groups = ['砷渣', 'CK', 'A', 'B']
    fig, ax = plt.subplots(figsize=figsize)
    x = np.arange(4)
    means = [np.mean(data['total_As'][g]) for g in groups]
    sds   = [np.std(data['total_As'][g], ddof=1) for g in groups]

    ax.bar(x, means, yerr=sds, color=[GC4[g] for g in groups],
           edgecolor='#333', linewidth=0.6, capsize=4, width=0.55,
           error_kw={'linewidth': 0.8, 'capthick': 0.8})
    ax.set_xticks(x)
    ax.set_xticklabels(get_labels(lang, 4), fontsize=9, fontfamily=fm_)
    ax.set_yscale('log')
    ax.set_ylim(0.5, 500)
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.yaxis.get_major_formatter().set_scientific(False)
    style_ax(ax, "Total As (ppm)" if not is_cn else "总砷含量（ppm）", "d", fm_, fn_)

    for i, m in enumerate(means):
        ax.text(i, m * 1.35, f"{m:.1f}", ha='center', fontsize=6.5,
                fontfamily=fn_, color='#444')

    fig.tight_layout()
    return fig


def draw_as_speciation(data, lang, figsize=(5.0, 3.2)):
    """As形态百分比堆叠图（4组），图例在图右侧外部，无框"""
    is_cn, fm_, fn_ = get_fonts(lang)
    groups = ['砷渣', 'CK', 'A', 'B']
    spec = data['As_spec']
    fig, ax = plt.subplots(figsize=figsize)
    x = np.arange(4)

    pct = {'As(III)':[], 'As(V)':[], 'MMA':[]}
    for g in groups:
        sub = spec[spec['group'] == g]
        tot = sub['total'].mean()
        pct['As(III)'].append(sub['As3'].mean()/tot*100 if tot > 0 else 0)
        pct['As(V)'].append(sub['As5'].mean()/tot*100 if tot > 0 else 0)
        pct['MMA'].append(sub['MMA'].mean()/tot*100 if tot > 0 else 0)

    bottom = np.zeros(4)
    for sp_name in ['As(V)', 'MMA', 'As(III)']:
        vals = pct[sp_name]
        ax.bar(x, vals, bottom=bottom, color=SPEC_COLORS[sp_name],
               edgecolor='white', linewidth=0.3, width=0.55, label=sp_name)
        for i, v in enumerate(vals):
            if v > 2:
                ax.text(i, bottom[i]+v/2, f"{v:.1f}%", ha='center', va='center',
                        fontsize=6, fontfamily=fn_, color='white', fontweight='bold')
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(get_labels(lang, 4), fontsize=9, fontfamily=fm_)
    ax.set_ylim(0, 108)
    style_ax(ax, "Proportion (%)" if not is_cn else "砷形态占比（%）", "e", fm_, fn_)

    # 图例：放到绘图区右侧外部，无边框
    ax.legend(fontsize=7, frameon=False,
              loc='upper left', bbox_to_anchor=(1.02, 1.0),
              prop={'family': fm_, 'size': 7})

    fig.tight_layout(rect=[0, 0, 0.85, 1])  # 右侧留15%空间给图例
    return fig


# ============================================================
# 4. 合并图
# ============================================================

def draw_combined(data, lang):
    """2行3列合并图，第6格(f)留空"""
    is_cn, fm_, fn_ = get_fonts(lang)
    groups_3 = ['CK', 'A', 'B']
    groups_4 = ['砷渣', 'CK', 'A', 'B']

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5) if is_cn else (11, 6.5))

    # ---- (a) pH ----
    ax = axes[0, 0]
    x = np.arange(3)
    means = [np.mean(data['pH'][g]) for g in groups_3]
    sds   = [np.std(data['pH'][g], ddof=1) for g in groups_3]
    ax.bar(x, means, yerr=sds, color=[GC[g] for g in groups_3],
           edgecolor='#333', linewidth=0.6, capsize=4, width=0.55,
           error_kw={'linewidth':0.8, 'capthick':0.8})
    ax.set_xticks(x)
    ax.set_xticklabels(get_labels(lang, 3), fontsize=9, fontfamily=fm_)
    ax.set_ylim(6.5, 8.2)
    style_ax(ax, "pH", "a", fm_, fn_)
    # 显著性
    sigs = calc_sig(data['pH'], groups_3)
    yrange = 8.2 - 6.5
    y_base = max(m+s for m,s in zip(means,sds)) + yrange*0.03
    sorted_p = sorted(sigs.keys(), key=lambda p: abs(groups_3.index(p[1])-groups_3.index(p[0])))
    for i, pair in enumerate(sorted_p):
        star = sig_star(sigs[pair])
        if star:
            add_sig_bracket(ax, groups_3.index(pair[0]), groups_3.index(pair[1]),
                            y_base + i*yrange*0.08, star, fn_)

    # ---- (b) Eh ----
    ax = axes[0, 1]
    means = [np.mean(data['Eh'][g]) for g in groups_3]
    sds   = [np.std(data['Eh'][g], ddof=1) for g in groups_3]
    ax.bar(x, means, yerr=sds, color=[GC[g] for g in groups_3],
           edgecolor='#333', linewidth=0.6, capsize=4, width=0.55,
           error_kw={'linewidth':0.8, 'capthick':0.8})
    ax.set_xticks(x)
    ax.set_xticklabels(get_labels(lang, 3), fontsize=9, fontfamily=fm_)
    # 关键修复：Eh的ylim上限留够空间给3层显著性标注
    ax.set_ylim(170, 210)
    style_ax(ax, "Eh (mV)" if not is_cn else "Eh（mV）", "b", fm_, fn_)
    sigs = calc_sig(data['Eh'], groups_3)
    yrange_eh = 210 - 170
    y_base_eh = max(m+s for m,s in zip(means,sds)) + yrange_eh*0.03
    sorted_p = sorted(sigs.keys(), key=lambda p: abs(groups_3.index(p[1])-groups_3.index(p[0])))
    for i, pair in enumerate(sorted_p):
        star = sig_star(sigs[pair])
        if star:
            add_sig_bracket(ax, groups_3.index(pair[0]), groups_3.index(pair[1]),
                            y_base_eh + i*yrange_eh*0.08, star, fn_)

    # ---- (c) TOC ----
    ax = axes[0, 2]
    x4 = np.arange(4)
    toc_vals = [data['TOC'].get(g, 0) for g in groups_4]
    ax.bar(x4, toc_vals, color=[GC4[g] for g in groups_4],
           edgecolor='#333', linewidth=0.6, width=0.55)
    ax.set_xticks(x4)
    ax.set_xticklabels(get_labels(lang, 4), fontsize=9, fontfamily=fm_)
    ax.set_ylim(0, 13)
    style_ax(ax, "TOC (g/kg)" if not is_cn else "TOC（g/kg）", "c", fm_, fn_)
    for i, v in enumerate(toc_vals):
        ax.text(i, v+0.3, f"{v:.1f}", ha='center', fontsize=7, fontfamily=fn_, color='#444')

    # ---- (d) 总As ----
    ax = axes[1, 0]
    as_means = [np.mean(data['total_As'][g]) for g in groups_4]
    as_sds   = [np.std(data['total_As'][g], ddof=1) for g in groups_4]
    ax.bar(x4, as_means, yerr=as_sds, color=[GC4[g] for g in groups_4],
           edgecolor='#333', linewidth=0.6, capsize=4, width=0.55,
           error_kw={'linewidth':0.8, 'capthick':0.8})
    ax.set_xticks(x4)
    ax.set_xticklabels(get_labels(lang, 4), fontsize=9, fontfamily=fm_)
    ax.set_yscale('log')
    ax.set_ylim(0.5, 500)
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.yaxis.get_major_formatter().set_scientific(False)
    style_ax(ax, "Total As (ppm)" if not is_cn else "总砷含量（ppm）", "d", fm_, fn_)
    for i, m in enumerate(as_means):
        ax.text(i, m*1.35, f"{m:.1f}", ha='center', fontsize=6.5, fontfamily=fn_, color='#444')

    # ---- (e) As形态占比 ----
    ax = axes[1, 1]
    spec = data['As_spec']
    pct = {'As(III)':[], 'As(V)':[], 'MMA':[]}
    for g in groups_4:
        sub = spec[spec['group']==g]
        tot = sub['total'].mean()
        pct['As(III)'].append(sub['As3'].mean()/tot*100 if tot>0 else 0)
        pct['As(V)'].append(sub['As5'].mean()/tot*100 if tot>0 else 0)
        pct['MMA'].append(sub['MMA'].mean()/tot*100 if tot>0 else 0)

    bottom = np.zeros(4)
    for sp_name in ['As(V)', 'MMA', 'As(III)']:
        vals = pct[sp_name]
        ax.bar(x4, vals, bottom=bottom, color=SPEC_COLORS[sp_name],
               edgecolor='white', linewidth=0.3, width=0.55, label=sp_name)
        for i, v in enumerate(vals):
            if v > 2:
                ax.text(i, bottom[i]+v/2, f"{v:.1f}%", ha='center', va='center',
                        fontsize=6, fontfamily=fn_, color='white', fontweight='bold')
        bottom += vals

    ax.set_xticks(x4)
    ax.set_xticklabels(get_labels(lang, 4), fontsize=9, fontfamily=fm_)
    ax.set_ylim(0, 108)
    style_ax(ax, "Proportion (%)" if not is_cn else "砷形态占比（%）", "e", fm_, fn_)

    # ---- (f) 留空 → 用来放(e)的图例 ----
    axes[1, 2].axis('off')
    # 图例放到右侧留空的(f)格子区域，无边框
    ax.legend(fontsize=8, frameon=False,
              loc='upper left', bbox_to_anchor=(1.15, 1.0),
              prop={'family': fm_, 'size': 8})

    fig.tight_layout(w_pad=2.5, h_pad=2.5)
    return fig


# ============================================================
# 5. 保存工具
# ============================================================

def save_fig(fig, out_dir, basename, formats=('pdf', 'png')):
    os.makedirs(out_dir, exist_ok=True)
    for fmt in formats:
        fp = os.path.join(out_dir, f"{basename}.{fmt}")
        fig.savefig(fp, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"  → {fp}")


# ============================================================
# 6. 主程序
# ============================================================

if __name__ == "__main__":
    # 输入文件
    if len(sys.argv) > 1:
        xlsx_path = sys.argv[1]
    else:
        candidates = [
            os.path.join(DATA_RAW, "elementdata", "ordata.xlsx"),
            os.path.join(DATA_RAW, "ordata.xlsx"),
            "/mnt/user-data/uploads/ordata.xlsx",
        ]
        xlsx_path = next((c for c in candidates if os.path.exists(c)), None)
        if xlsx_path is None:
            print("ERROR: ordata.xlsx not found")
            sys.exit(1)

    print(f"Reading: {xlsx_path}")
    data = load_data(xlsx_path)

    for lang in ['en', 'cn']:
        is_cn = (lang == 'cn')
        out_dir = FIG_THESIS if is_cn else FIG_PAPER
        tag = "论文" if is_cn else "Paper"
        print(f"\n=== {tag} version ===")

        # (a) pH 单独
        print(f"[{tag}] pH...")
        fig = draw_bar3(data['pH'], ['CK','A','B'],
                        "pH", (6.5, 8.2), lang, "a")
        save_fig(fig, out_dir,
                 "Fig1-2a_pH" if not is_cn else "Fig1-2a_pH")
        plt.close(fig)

        # (b) Eh 单独
        print(f"[{tag}] Eh...")
        fig = draw_bar3(data['Eh'], ['CK','A','B'],
                        "Eh (mV)" if not is_cn else "Eh（mV）",
                        (170, 210), lang, "b")
        save_fig(fig, out_dir,
                 "Fig1-2b_Eh" if not is_cn else "Fig1-2b_Eh")
        plt.close(fig)

        # (c) TOC 单独
        print(f"[{tag}] TOC...")
        fig = draw_toc(data, lang)
        save_fig(fig, out_dir,
                 "Fig1-2c_TOC" if not is_cn else "Fig1-2c_TOC")
        plt.close(fig)

        # (d) 总As 单独
        print(f"[{tag}] Total As...")
        fig = draw_total_as(data, lang)
        save_fig(fig, out_dir,
                 "Fig1-2d_totalAs" if not is_cn else "Fig1-2d_总As")
        plt.close(fig)

        # (e) As形态 单独
        print(f"[{tag}] As speciation...")
        fig = draw_as_speciation(data, lang)
        save_fig(fig, out_dir,
                 "Fig1-2e_As_speciation" if not is_cn else "Fig1-2e_As形态")
        plt.close(fig)

        # 合并图
        print(f"[{tag}] Combined...")
        fig = draw_combined(data, lang)
        save_fig(fig, out_dir,
                 "Fig1-2_physicochemical" if not is_cn else "Fig1-2_理化指标")
        plt.close(fig)

    print("\nDone! All figures saved.")
