#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07_MAG_abundance_heatmap.py — Fig3-4 Top30 MAG丰度热图
三段非线性配色：Blue 0~0.2, Green 0.2~0.5, Red 0.5~max

输出:
  figures/thesis_cn/Fig3-4_MAG_abundance_top30_CN.pdf + .png
  figures/paper_en/Fig3-4_MAG_abundance_top30_EN.pdf + .png + .tiff

用法:
  python3 scripts/python/07_MAG_abundance_heatmap.py \
      --abund  data/raw/abundance.tsv \
      --bac    data/raw/tax_bac120_summary.tsv \
      --ar     data/raw/tax_ar53_summary.tsv \
      --ks     data/raw/keystone_species.txt \
      --outdir figures \
      --statdir data/processed
"""
import argparse, os, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.font_manager as fm
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

# ================================================================
# 字体
# ================================================================
def setup_fonts():
    avail = {f.name for f in fm.fontManager.ttflist}
    sans = 'DejaVu Sans'
    for c in ['Arial', 'Helvetica', 'Liberation Sans']:
        if c in avail: sans = c; break
    serif = 'DejaVu Serif'
    for c in ['Times New Roman', 'Liberation Serif']:
        if c in avail: serif = c; break
    cn = serif
    for c in ['SimSun', '宋体', 'WenQuanYi Micro Hei', 'Noto Sans CJK JP']:
        if c in avail: cn = c; break
    return sans, serif, cn

F_SANS, F_SERIF, F_CN = setup_fonts()
plt.rcParams['axes.unicode_minus'] = False

def get_fonts(style):
    if style == 'cn': return F_CN, F_SERIF
    return F_SANS, F_SANS

# ================================================================
# 配色
# ================================================================
GROUP_COLORS = {'CK': '#1c9cbd', 'A': '#e3943d', 'B': '#92181e'}

PHYLUM_COLORS = {
    'Pseudomonadota':'#4393C3','Acidobacteriota':'#8B6914','Chloroflexota':'#7FBC41',
    'Bacteroidota_A':'#F4A582','Patescibacteriota':'#C2A5CF','Desulfobacterota':'#D6604D',
    'Actinomycetota':'#B2ABD2','Gemmatimonadota':'#E08214','Planctomycetota':'#FDB863',
    'Myxococcota':'#542788','Nitrospirota':'#1B7837','Zixibacteria':'#A6DBA0',
    'Verrucomicrobiota':'#5AAE61','Methylomirabilota':'#878787','Desulfobacterota_B':'#E78AC3',
    'Desulfobacterota_E':'#FC8D59','Bacteroidota':'#FEE0B6','Thermoproteota':'#B35806',
    'Halobacteriota':'#8073AC','Eisenbacteria':'#35978F','Desulfobacterota_C':'#BF812D',
    'Cyanobacteriota':'#01665E','Chlamydiota':'#80CDC1',
}

TOP10_BC = '#E67E22'
KS_BC = '#8E44AD'
BOTH_BC = '#1ABC9C'

SAMPLE_ORDER = ['2_1','2_5','2_7','8_1','8_2','8_3','8_4','9_1','9_2','9_3']
SAMPLE_GROUP = {'2_1':'CK','2_5':'CK','2_7':'CK',
                '8_1':'A','8_2':'A','8_3':'A','8_4':'A',
                '9_1':'B','9_2':'B','9_3':'B'}

# ================================================================
# 数据加载
# ================================================================
def load_data(args):
    # 丰度
    ab = pd.read_csv(args.abund, sep='\t')
    ab = ab[ab['Genome'] != 'unmapped']
    ab['mean_all'] = ab[SAMPLE_ORDER].mean(axis=1)
    ab = ab.sort_values('mean_all', ascending=False).reset_index(drop=True)
    top30 = ab.head(30).copy()
    top30['rank'] = range(1, 31)

    # GTDB分类
    bac = pd.read_csv(args.bac, sep='\t')
    ar = pd.read_csv(args.ar, sep='\t')
    gtdb = pd.concat([bac[['user_genome','classification']],
                       ar[['user_genome','classification']]], ignore_index=True)
    gtdb = gtdb.rename(columns={'user_genome':'Genome'})

    def gl(cls, pfx):
        for p in cls.split(';'):
            if p.startswith(pfx): return p.replace(pfx, '')
        return ''

    gtdb['Phylum'] = gtdb['classification'].apply(lambda x: gl(x, 'p__'))
    gtdb['Family'] = gtdb['classification'].apply(lambda x: gl(x, 'f__'))
    gtdb['Genus'] = gtdb['classification'].apply(lambda x: gl(x, 'g__'))
    gtdb['Species'] = gtdb['classification'].apply(lambda x: gl(x, 's__'))

    def mk_label(row):
        n = row['Genome'].replace('Mx_All_', '')
        g, s, f = row['Genus'], row['Species'], row['Family']
        if g:
            if s:
                ss = s.replace(g + ' ', '') if g in s else s
                return f"{g} {ss} [{n}]"
            return f"{g} sp. [{n}]"
        elif f:
            return f"{f} (f) [{n}]"
        return f"Unclassified [{n}]"

    gtdb['label'] = gtdb.apply(mk_label, axis=1)
    top30 = top30.merge(gtdb[['Genome','Phylum','label']], on='Genome', how='left')

    # Top10 / Keystone
    top10_mags = set(ab.head(10)['Genome'].tolist())
    ks = pd.read_csv(args.ks, sep='\t')
    ks_mags = set(ks['MAG'].tolist())
    top30['is_top10'] = top30['Genome'].isin(top10_mags)
    top30['is_keystone'] = top30['Genome'].isin(ks_mags)

    return top30, ks_mags, top10_mags


# ================================================================
# 三段非线性配色
# ================================================================
def build_tricolor_cmap(mat):
    """根据数据构建三段配色: Blue 0~0.2, YlGn 0.2~0.5, OrgRed 0.5~max"""
    vmax = np.ceil(mat.max() * 10) / 10  # 向上取整到0.1
    if vmax < 0.6: vmax = 0.6

    # 分段边界：低段细分，高段粗分
    bounds = list(np.arange(0, 0.2, 0.025)) + \
             list(np.arange(0.2, 0.5, 0.05)) + \
             list(np.arange(0.5, min(vmax, 1.0), 0.1))
    if vmax > 1.0:
        bounds += list(np.arange(1.0, vmax + 0.01, 0.5))
    # 去重并确保vmax在内
    bounds = sorted(set([round(b, 3) for b in bounds]))
    if bounds[-1] < vmax:
        bounds.append(round(vmax, 2))

    n_lo = len([b for b in bounds if b <= 0.2]) - 1   # 0~0.2段数
    n_mid = len([b for b in bounds if 0.2 < b <= 0.5]) # 0.2~0.5段数
    n_hi = len([b for b in bounds if b > 0.5])          # 0.5+段数
    n_total = n_lo + n_mid + n_hi

    # 生成颜色节点
    seg_lo = plt.cm.Blues(np.linspace(0.08, 0.55, max(n_lo + 1, 2)))
    seg_mid = plt.cm.YlGn(np.linspace(0.18, 0.65, max(n_mid + 1, 2)))
    seg_hi = plt.cm.YlOrRd(np.linspace(0.40, 0.95, max(n_hi + 1, 2)))

    all_colors = list(seg_lo[:n_lo+1]) + list(seg_mid[1:n_mid+1]) + list(seg_hi[1:n_hi+1])
    # 确保颜色数 = bounds数 - 1 + 1（BoundaryNorm需要）
    cmap = LinearSegmentedColormap.from_list('tri', all_colors, N=256)
    norm = BoundaryNorm(bounds, ncolors=256)

    return cmap, norm, bounds, vmax


# ================================================================
# 绘图
# ================================================================
def draw_figure(top30, style, outpath, cluster=True):
    fm_, fn_ = get_fonts(style)
    is_cn = (style == 'cn')
    n_samples = len(SAMPLE_ORDER)

    mat = top30[SAMPLE_ORDER].values
    labels_y = top30['label'].tolist()
    phyla = top30['Phylum'].tolist()

    # 行聚类排序（门内聚类）
    if cluster:
        # 先按门分组，门内按聚类排序
        top30_sorted = top30.copy()
        phy_order = top30['Phylum'].value_counts().index.tolist()
        top30_sorted['phy_rank'] = top30_sorted['Phylum'].map(
            {p: i for i, p in enumerate(phy_order)})

        new_order = []
        for phy in phy_order:
            idx = top30_sorted[top30_sorted['Phylum'] == phy].index.tolist()
            if len(idx) > 2:
                sub_mat = mat[idx]
                dist = pdist(sub_mat, metric='euclidean')
                link = linkage(dist, method='ward')
                order = leaves_list(link)
                new_order.extend([idx[o] for o in order])
            else:
                new_order.extend(idx)

        mat = mat[new_order]
        labels_y = [labels_y[i] for i in new_order]
        phyla = [phyla[i] for i in new_order]
        is_top10_list = [top30.iloc[i]['is_top10'] for i in new_order]
        is_ks_list = [top30.iloc[i]['is_keystone'] for i in new_order]
        ranks = [top30.iloc[i]['rank'] for i in new_order]
    else:
        is_top10_list = top30['is_top10'].tolist()
        is_ks_list = top30['is_keystone'].tolist()
        ranks = top30['rank'].tolist()

    n_mags = len(mat)
    cmap, norm, bounds, vmax = build_tricolor_cmap(mat)

    # 样品标签
    if is_cn:
        sample_labels = ['CK1','CK2','CK3','A1','A2','A3','A4','B1','B2','B3']
    else:
        sample_labels = ['CK1','CK2','CK3','A1','A2','A3','A4','B1','B2','B3']

    # 尺寸
    fig_h = n_mags * 0.38 + 4
    fig_w = n_samples * 0.8 + 8

    fig = plt.figure(figsize=(fig_w, fig_h))
    # [phylum | Y labels | heatmap | colorbar+legend]
    gs = gridspec.GridSpec(1, 4,
        width_ratios=[0.3, 3.5, n_samples * 0.8, 3.0],
        wspace=0.01)
    ax_phy = fig.add_subplot(gs[0, 0])
    ax_ylab = fig.add_subplot(gs[0, 1])
    ax_heat = fig.add_subplot(gs[0, 2])
    ax_leg = fig.add_subplot(gs[0, 3])

    # ---- 热图 ----
    im = ax_heat.imshow(mat, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')

    ax_heat.set_xticks(range(n_samples))
    ax_heat.set_xticklabels(sample_labels, rotation=45, ha='right',
                             fontsize=8, fontfamily=fn_)
    ax_heat.set_yticks([])
    ax_heat.set_xlim(-0.5, n_samples - 0.5)
    ax_heat.set_ylim(n_mags - 0.5, -0.5)

    # 分组竖线
    ax_heat.axvline(x=2.5, color='black', lw=1.5)
    ax_heat.axvline(x=6.5, color='black', lw=1.5)

    # 顶部分组色条
    for j, s in enumerate(SAMPLE_ORDER):
        grp = SAMPLE_GROUP[s]
        ax_heat.add_patch(Rectangle((j-0.5, -1.8), 1, 1.2,
                          facecolor=GROUP_COLORS[grp], edgecolor='none',
                          clip_on=False, alpha=0.8))
    # 分组标签
    grp_pos = {'CK': 1, 'A': 4.5, 'B': 8}
    for grp, xpos in grp_pos.items():
        ax_heat.text(xpos, -2.8, grp, ha='center', va='center',
                    fontsize=10, fontweight='bold', fontfamily=fn_,
                    color=GROUP_COLORS[grp])

    # 行边框标注
    for i in range(n_mags):
        if is_top10_list[i] and is_ks_list[i]:
            bc, bw = BOTH_BC, 2.0
        elif is_top10_list[i]:
            bc, bw = TOP10_BC, 1.8
        elif is_ks_list[i]:
            bc, bw = KS_BC, 1.8
        else:
            continue
        ax_heat.add_patch(Rectangle((-0.5, i-0.5), n_samples, 1,
                          facecolor='none', edgecolor=bc, linewidth=bw, zorder=4))

    # 门分组线
    prev_phy = None
    for i, phy in enumerate(phyla):
        if phy != prev_phy:
            if prev_phy is not None:
                ax_heat.axhline(y=i-0.5, color='#666', lw=0.6, zorder=3)
            prev_phy = phy

    # ---- 门颜色条 ----
    for i, phy in enumerate(phyla):
        ax_phy.add_patch(Rectangle((0, i-0.5), 1, 1,
                         facecolor=PHYLUM_COLORS.get(phy, '#BDC3C7'), edgecolor='none'))
    ax_phy.set_xlim(0, 1); ax_phy.set_ylim(n_mags-0.5, -0.5); ax_phy.axis('off')

    # ---- Y轴标签 ----
    ax_ylab.set_xlim(0, 1); ax_ylab.set_ylim(n_mags-0.5, -0.5); ax_ylab.axis('off')
    for i in range(n_mags):
        color = PHYLUM_COLORS.get(phyla[i], '#666')
        is_sp = is_top10_list[i] or is_ks_list[i]
        fw = 'bold' if is_sp else 'normal'
        fs = 7.5 if is_sp else 7
        ax_ylab.text(0.98, i, labels_y[i], ha='right', va='center',
                    fontsize=fs, fontfamily=fn_, fontstyle='italic',
                    fontweight=fw, color=color)
        # 标记
        x_cur = 0.01
        if is_top10_list[i]:
            r = int(ranks[i])
            ax_ylab.plot(x_cur, i, marker='o', markersize=6, color=TOP10_BC,
                        markeredgewidth=0, zorder=5, clip_on=False)
            ax_ylab.text(x_cur + 0.025, i, str(r), ha='left', va='center',
                        fontsize=5, fontfamily=fn_, fontweight='bold', color=TOP10_BC)
            x_cur += 0.06
        if is_ks_list[i]:
            ax_ylab.plot(x_cur, i, marker='D', markersize=4, color=KS_BC,
                        markeredgewidth=0, zorder=5, clip_on=False)

    # ---- 图例 ----
    ax_leg.axis('off')
    y_pos = 0.97

    # 门图例
    ttl = '门水平' if is_cn else 'Phylum'
    ax_leg.text(0.05, y_pos, ttl, fontsize=8, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_pos -= 0.02
    sub_phyla = list(dict.fromkeys(phyla))  # 保持出现顺序去重
    for phy in sub_phyla:
        y_pos -= 0.032
        if y_pos < 0.55: break
        ax_leg.add_patch(Rectangle((0.05, y_pos-0.008), 0.06, 0.02,
                         facecolor=PHYLUM_COLORS.get(phy, '#BDC3C7'),
                         edgecolor='none', transform=ax_leg.transAxes, clip_on=False))
        ax_leg.text(0.13, y_pos, phy, fontsize=6, fontfamily=fn_,
                   transform=ax_leg.transAxes, va='center')

    # 色标（竖向，刻度均匀分布，每个色块内部渐变填充）
    y_pos -= 0.04
    ttl2 = '相对丰度 (%)' if is_cn else 'Relative abundance (%)'
    ax_leg.text(0.05, y_pos, ttl2, fontsize=8, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_pos -= 0.015

    # 离散刻度值（均匀排列）
    tick_vals = [0, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0]
    if vmax > 2.0:
        tick_vals.append(round(vmax, 1))
    tick_vals = [v for v in tick_vals if v <= vmax + 0.01]

    bar_x = 0.05
    bar_w = 0.04
    block_h = 0.025  # 每个刻度区间等高
    n_blocks = len(tick_vals) - 1
    bar_h_total = block_h * n_blocks
    sub_steps = 10  # 每个色块内部细分步数，实现渐变

    for b in range(n_blocks):
        lo_val = tick_vals[b]
        hi_val = tick_vals[b + 1]
        for s in range(sub_steps):
            frac = s / sub_steps
            val = lo_val + (hi_val - lo_val) * (frac + 0.5 / sub_steps)
            c = cmap(norm(val))
            y_sub = y_pos - bar_h_total + b * block_h + frac * block_h
            ax_leg.add_patch(Rectangle((bar_x, y_sub), bar_w, block_h / sub_steps,
                             facecolor=c, edgecolor='none',
                             transform=ax_leg.transAxes, clip_on=False))

    # 刻度标签
    for i, tv in enumerate(tick_vals):
        y_tick = y_pos - bar_h_total + i * block_h
        ax_leg.plot([bar_x + bar_w, bar_x + bar_w + 0.012],
                   [y_tick, y_tick],
                   color='black', lw=0.5, transform=ax_leg.transAxes, clip_on=False)
        ax_leg.text(bar_x + bar_w + 0.018, y_tick, f'{tv:.1f}',
                   fontsize=5.5, fontfamily=fn_, ha='left', va='center',
                   transform=ax_leg.transAxes)

    # 色条外框
    ax_leg.add_patch(Rectangle((bar_x, y_pos - bar_h_total), bar_w, bar_h_total,
                     facecolor='none', edgecolor='black', linewidth=0.5,
                     transform=ax_leg.transAxes, clip_on=False))

    y_pos -= bar_h_total + 0.02

    # 分组色标
    y_pos -= 0.03
    ttl3 = '实验分组' if is_cn else 'Treatment group'
    ax_leg.text(0.05, y_pos, ttl3, fontsize=8, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    grp_labels = [('CK','CK组','Group CK'), ('A','A组','Group A'), ('B','B组','Group B')]
    for grp, gname_cn, gname_en in grp_labels:
        y_pos -= 0.035
        ax_leg.add_patch(Rectangle((0.05, y_pos-0.008), 0.06, 0.02,
                         facecolor=GROUP_COLORS[grp], edgecolor='none',
                         transform=ax_leg.transAxes, clip_on=False))
        gname = gname_cn if is_cn else gname_en
        ax_leg.text(0.13, y_pos, gname, fontsize=6, fontfamily=fm_ if is_cn else fn_,
                   transform=ax_leg.transAxes, va='center')

    # 行标注图例
    y_pos -= 0.05
    ttl4 = '行标注' if is_cn else 'Row annotation'
    ax_leg.text(0.05, y_pos, ttl4, fontsize=8, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')

    y_pos -= 0.03
    ax_leg.add_patch(Rectangle((0.05, y_pos-0.008), 0.08, 0.02,
                     facecolor='none', edgecolor=TOP10_BC, linewidth=2,
                     transform=ax_leg.transAxes, clip_on=False))
    ax_leg.plot(0.15, y_pos, marker='o', markersize=6, color=TOP10_BC,
               markeredgewidth=0, transform=ax_leg.transAxes, clip_on=False)
    t10_txt = 'n  丰度排名' if is_cn else 'n  Abundance rank'
    ax_leg.text(0.18, y_pos, t10_txt, fontsize=6,
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='center')

    y_pos -= 0.035
    ax_leg.add_patch(Rectangle((0.05, y_pos-0.008), 0.08, 0.02,
                     facecolor='none', edgecolor=KS_BC, linewidth=2,
                     transform=ax_leg.transAxes, clip_on=False))
    ax_leg.plot(0.15, y_pos, marker='D', markersize=5, color=KS_BC,
               markeredgewidth=0, transform=ax_leg.transAxes, clip_on=False)
    ks_txt = '关键物种' if is_cn else 'Keystone species'
    ax_leg.text(0.18, y_pos, ks_txt, fontsize=6,
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='center')

    y_pos -= 0.035
    ax_leg.add_patch(Rectangle((0.05, y_pos-0.008), 0.08, 0.02,
                     facecolor='none', edgecolor=BOTH_BC, linewidth=2,
                     transform=ax_leg.transAxes, clip_on=False))
    ax_leg.plot(0.15, y_pos, marker='o', markersize=6, color=TOP10_BC,
               markeredgewidth=0, transform=ax_leg.transAxes, clip_on=False)
    ax_leg.plot(0.17, y_pos, marker='D', markersize=5, color=KS_BC,
               markeredgewidth=0, transform=ax_leg.transAxes, clip_on=False)
    both_txt = '两者兼具' if is_cn else 'Both'
    ax_leg.text(0.20, y_pos, both_txt, fontsize=6,
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='center')

    # 标题
    if is_cn:
        title = 'Top 30 MAGs相对丰度热图'
    else:
        title = 'Relative Abundance of Top 30 MAGs'
    fig.suptitle(title, fontsize=14, fontweight='bold',
                fontfamily=fm_ if is_cn else fn_, y=0.98)

    # 保存 PDF + PNG
    dpi = 300 if is_cn else 600
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight', facecolor='white')
    png_path = outpath.rsplit('.', 1)[0] + '.png'
    plt.savefig(png_path, dpi=min(dpi, 300), bbox_inches='tight', facecolor='white')
    print(f"  Saved: {outpath} + .png ({dpi}dpi)")
    plt.close()


# ================================================================
# 数据导出
# ================================================================
def export_stats(top30, stat_dir):
    os.makedirs(stat_dir, exist_ok=True)

    # 1. Top30汇总
    with open(os.path.join(stat_dir, 'Fig3-4_Top30_abundance_summary.txt'), 'w') as f:
        f.write("=== Top 30 MAGs Abundance Summary ===\n\n")
        f.write(f"{'Rank':<5} {'MAG':<15} {'Phylum':<22} {'Species':<32} "
                f"{'CK_mean':>8} {'A_mean':>8} {'B_mean':>8} {'Overall':>8} "
                f"{'Top10':>5} {'KS':>3}\n")
        f.write("-" * 125 + "\n")
        for _, row in top30.iterrows():
            ck = row[['2_1','2_5','2_7']].mean()
            a = row[['8_1','8_2','8_3','8_4']].mean()
            b = row[['9_1','9_2','9_3']].mean()
            t10 = 'Y' if row['is_top10'] else ''
            ks = 'Y' if row['is_keystone'] else ''
            f.write(f"{int(row['rank']):<5} {row['Genome']:<15} {row['Phylum']:<22} "
                    f"{row['label']:<32} {ck:>8.4f} {a:>8.4f} {b:>8.4f} "
                    f"{row['mean_all']:>8.4f} {t10:>5} {ks:>3}\n")

    # 2. 组间差异统计
    with open(os.path.join(stat_dir, 'Fig3-4_group_enrichment.txt'), 'w') as f:
        f.write("=== Group Enrichment Analysis (Top30) ===\n\n")
        f.write("Enrichment defined as: group_mean > 1.5x other groups mean\n\n")
        f.write(f"{'MAG':<15} {'Species':<32} {'CK':>8} {'A':>8} {'B':>8} {'Enriched_in':>12}\n")
        f.write("-" * 90 + "\n")
        for _, row in top30.iterrows():
            ck = row[['2_1','2_5','2_7']].mean()
            a = row[['8_1','8_2','8_3','8_4']].mean()
            b = row[['9_1','9_2','9_3']].mean()
            vals = {'CK': ck, 'A': a, 'B': b}
            enriched = []
            for grp, val in vals.items():
                others = [v for k, v in vals.items() if k != grp]
                if val > 1.5 * np.mean(others):
                    enriched.append(grp)
            enr_str = ','.join(enriched) if enriched else '-'
            f.write(f"{row['Genome']:<15} {row['label']:<32} "
                    f"{ck:>8.4f} {a:>8.4f} {b:>8.4f} {enr_str:>12}\n")

    print(f"  Stats exported to: {stat_dir}")


# ================================================================
# Main
# ================================================================
def main():
    parser = argparse.ArgumentParser(description='Fig3-4 Top30 MAG Abundance Heatmap')
    parser.add_argument('--abund', required=True, help='CoverM abundance.tsv')
    parser.add_argument('--bac',   required=True, help='GTDB bac120 summary')
    parser.add_argument('--ar',    required=True, help='GTDB ar53 summary')
    parser.add_argument('--ks',    required=True, help='Keystone species list')
    parser.add_argument('--outdir', default='figures', help='Output base directory')
    parser.add_argument('--statdir', default='data/processed', help='Stats output directory')
    args = parser.parse_args()

    print(f"Fonts: sans={F_SANS}, serif={F_SERIF}, cn={F_CN}")

    top30, ks_mags, top10_mags = load_data(args)
    print(f"Data: Top 30 MAGs × {len(SAMPLE_ORDER)} samples")
    print(f"  Top10 in Top30: {top30['is_top10'].sum()}")
    print(f"  Keystone in Top30: {top30['is_keystone'].sum()}")

    en_dir = os.path.join(args.outdir, 'paper_en')
    cn_dir = os.path.join(args.outdir, 'thesis_cn')
    os.makedirs(en_dir, exist_ok=True)
    os.makedirs(cn_dir, exist_ok=True)

    # EN
    print("\n=== EN version ===")
    draw_figure(top30, 'en', os.path.join(en_dir, 'Fig3-4_MAG_abundance_top30_EN.pdf'))
    # EN TIFF
    draw_figure(top30, 'en', os.path.join(en_dir, 'Fig3-4_MAG_abundance_top30_EN.tiff'))

    # CN
    print("\n=== CN version ===")
    draw_figure(top30, 'cn', os.path.join(cn_dir, 'Fig3-4_MAG_abundance_top30_CN.pdf'))

    # Stats
    print("\n=== Exporting stats ===")
    export_stats(top30, args.statdir)

    print("\nAll done!")

if __name__ == '__main__':
    main()
