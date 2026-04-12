#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_pathway_completeness.py — Fig3-5 通路完整度
  Fig3-5a: 168 MAG通路完整度热图（标注Top10+Keystone）
  Fig3-5b: Top30气泡图
  Fig3-5c: Keystone气泡图

输出:
  figures/paper_en/Fig3-5a_pathway_completeness_EN.pdf + .png
  figures/paper_en/Fig3-5b_pathway_top30_bubble_EN.pdf + .png
  figures/paper_en/Fig3-5c_pathway_keystone_bubble_EN.pdf + .png
  figures/thesis_cn/ (同上CN版)
  data/processed/Fig3-5_*.txt

用法:
  python3 scripts/python/08_pathway_completeness.py \
      --itol   data/raw/07_element_KO_heatmap.txt \
      --ko     data/raw/kegg_target_only.tsv \
      --bac    data/raw/tax_bac120_summary.tsv \
      --ar     data/raw/tax_ar53_summary.tsv \
      --abund  data/raw/abundance.tsv \
      --ks     data/raw/keystone_species.txt \
      --outdir figures \
      --statdir data/processed
"""
import argparse, os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.font_manager as fm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec

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
# 通路定义：17个功能分类及其包含的KO
# ================================================================
PATHWAY_KOS = {
    'Arsenate reduction':   ['K00537','K03741','K22547','K18701'],
    'Arsenite oxidation':   ['K08355','K08356'],
    'Resp. arsenate red.':  ['K28466','K28467'],
    'As transport/detox':   ['K03893','K03325','K18064','K01551','K25223','K11811'],
    'As methylation':       ['K07755','K25224'],
    'As regulation':        ['K03892'],
    'Nitrate reduction':    ['K00370','K00371','K00374','K02567','K02568','K00367'],
    'Nitrite reduction':    ['K00368','K15864','K00366','K15877'],
    'NO reduction':         ['K04561'],
    'N2O reduction':        ['K00376'],
    'Ammonia oxidation':    ['K10944','K10945','K10946'],
    'N fixation':           ['K02586','K02591'],
    'Assim. sulfate red.':  ['K00380','K00381','K00390','K00392'],
    'Dissim. sulfate red.': ['K00394','K00395','K11180','K11181'],
    'Sulfide oxidation':    ['K17218','K17222','K17223','K17224','K17226','K17227'],
    'Thiosulfate metab.':   ['K01011'],
    'Fe transport':         ['K02013','K02014','K02015','K02016','K02012','K07243'],
    'Fe uptake regulation': ['K03711','K03832'],
}

PATHWAY_ORDER = list(PATHWAY_KOS.keys())

# 通路→元素（显式映射，避免前缀误判）
PATHWAY_ELEMENT = {
    'Arsenate reduction': 'As', 'Arsenite oxidation': 'As',
    'Resp. arsenate red.': 'As', 'As transport/detox': 'As',
    'As methylation': 'As', 'As regulation': 'As',
    'Nitrate reduction': 'N', 'Nitrite reduction': 'N',
    'NO reduction': 'N', 'N2O reduction': 'N',
    'Ammonia oxidation': 'N', 'N fixation': 'N',
    'Assim. sulfate red.': 'S', 'Dissim. sulfate red.': 'S',
    'Sulfide oxidation': 'S', 'Thiosulfate metab.': 'S',
    'Fe transport': 'Fe', 'Fe uptake regulation': 'Fe',
}

ELEM_COLORS = {'As':'#C0392B', 'N':'#2C3E87', 'S':'#1E8449', 'Fe':'#B7600B'}

FUNC_BG = {
    'Arsenate reduction':'#F5D8D8','Arsenite oxidation':'#E0A8A0',
    'Resp. arsenate red.':'#C87870','As transport/detox':'#F0D0B8',
    'As methylation':'#D8A880','As regulation':'#C08060',
    'Nitrate reduction':'#D8E4F5','Nitrite reduction':'#A0B8D8',
    'NO reduction':'#6888B8','N2O reduction':'#E0D0E8',
    'Ammonia oxidation':'#B8A8D0','N fixation':'#8878A8',
    'Assim. sulfate red.':'#D0F0E0','Dissim. sulfate red.':'#88D0A8',
    'Sulfide oxidation':'#48A878','Thiosulfate metab.':'#E0E8B0',
    'Fe transport':'#F5E0C0','Fe uptake regulation':'#C8A070',
}

FUNC_CN = {
    'Arsenate reduction':'砷酸盐还原','Arsenite oxidation':'亚砷酸盐氧化',
    'Resp. arsenate red.':'呼吸型砷酸盐还原','As transport/detox':'砷转运与解毒',
    'As methylation':'砷甲基化','As regulation':'砷调控',
    'Nitrate reduction':'硝酸盐还原','Nitrite reduction':'亚硝酸盐还原',
    'NO reduction':'NO还原','N2O reduction':'N2O还原',
    'Ammonia oxidation':'氨氧化','N fixation':'固氮',
    'Assim. sulfate red.':'硫酸盐同化还原','Dissim. sulfate red.':'硫酸盐异化还原',
    'Sulfide oxidation':'硫化物氧化','Thiosulfate metab.':'硫代硫酸盐代谢',
    'Fe transport':'铁转运系统','Fe uptake regulation':'铁摄取调控',
}

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


def build_abundance_cmap(data_values):
    """根据实际丰度数据构建三段非线性配色，返回 cmap, norm, tick_vals, vmax"""
    from matplotlib.colors import BoundaryNorm
    vmax = max(np.max(data_values), 0.01)  # 避免全零

    # 判断数据分布，决定是否需要非线性
    q90 = np.percentile(data_values[data_values > 0], 90) if (data_values > 0).any() else 0.1
    q50 = np.percentile(data_values[data_values > 0], 50) if (data_values > 0).any() else 0.05

    # 构建边界
    if vmax <= 0.3:
        # 数据范围小，线性分段即可
        bounds = list(np.arange(0, vmax + 0.01, max(vmax/10, 0.01)))
    else:
        # 三段非线性
        bounds = list(np.arange(0, min(0.2, vmax), 0.025))
        if vmax > 0.2:
            bounds += list(np.arange(0.2, min(0.5, vmax), 0.05))
        if vmax > 0.5:
            bounds += list(np.arange(0.5, min(1.0, vmax), 0.1))
        if vmax > 1.0:
            bounds += list(np.arange(1.0, vmax + 0.01, 0.5))

    bounds = sorted(set([round(b, 3) for b in bounds]))
    if len(bounds) < 2:
        bounds = [0, round(vmax, 3)]
    if bounds[-1] < vmax:
        bounds.append(round(vmax, 3))

    # 颜色段
    n_lo = max(len([b for b in bounds if b <= 0.2]) - 1, 0)
    n_mid = max(len([b for b in bounds if 0.2 < b <= 0.5]), 0)
    n_hi = max(len([b for b in bounds if b > 0.5]), 0)

    seg_lo = plt.cm.Blues(np.linspace(0.08, 0.55, max(n_lo + 1, 2)))
    seg_mid = plt.cm.YlGn(np.linspace(0.18, 0.65, max(n_mid + 1, 2)))
    seg_hi = plt.cm.YlOrRd(np.linspace(0.40, 0.95, max(n_hi + 1, 2)))

    colors = list(seg_lo[:max(n_lo + 1, 1)])
    if n_mid > 0:
        colors += list(seg_mid[1:n_mid + 1])
    if n_hi > 0:
        colors += list(seg_hi[1:n_hi + 1])

    if len(colors) < 2:
        colors = [seg_lo[0], seg_lo[-1]]

    cmap = LinearSegmentedColormap.from_list('ab_auto', colors, N=256)
    norm = BoundaryNorm(bounds, ncolors=256)

    # 色标刻度（选取代表性值）
    candidates = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0]
    tick_vals = [v for v in candidates if v <= vmax]
    tick_vals.append(round(vmax, 2))
    tick_vals = sorted(set(tick_vals))

    return cmap, norm, tick_vals, vmax

# ================================================================
# 数据加载
# ================================================================
def load_data(args):
    # MAG KO注释（从kegg_target_only提取presence/absence）
    raw = pd.read_csv(args.ko, sep='\t')
    rows_exp = []
    for _, r in raw.iterrows():
        for ko in str(r['KEGG_ko']).split(','):
            ko_clean = ko.strip().replace('ko:', '')
            rows_exp.append((r['MAG'], ko_clean))
    df_exp = pd.DataFrame(rows_exp, columns=['MAG', 'KO'])
    # 每个MAG拥有的KO集合
    mag_kos = df_exp.groupby('MAG')['KO'].apply(set).to_dict()

    # GTDB分类
    bac = pd.read_csv(args.bac, sep='\t')
    ar = pd.read_csv(args.ar, sep='\t')
    gtdb = pd.concat([bac[['user_genome','classification']],
                       ar[['user_genome','classification']]], ignore_index=True)
    gtdb = gtdb.rename(columns={'user_genome': 'MAG'})

    def gl(cls, pfx):
        for p in cls.split(';'):
            if p.startswith(pfx): return p.replace(pfx, '')
        return ''
    gtdb['Phylum'] = gtdb['classification'].apply(lambda x: gl(x, 'p__'))
    gtdb['Family'] = gtdb['classification'].apply(lambda x: gl(x, 'f__'))
    gtdb['Genus'] = gtdb['classification'].apply(lambda x: gl(x, 'g__'))
    gtdb['Species'] = gtdb['classification'].apply(lambda x: gl(x, 's__'))

    def mk_label(row):
        n = row['MAG'].replace('Mx_All_', '')
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

    # 丰度 → Top10
    ab = pd.read_csv(args.abund, sep='\t')
    ab = ab[ab['Genome'] != 'unmapped']
    scols = [c for c in ab.columns if c != 'Genome']
    ab['mean_all'] = ab[scols].mean(axis=1)
    ab = ab.sort_values('mean_all', ascending=False).reset_index(drop=True)
    top10_rank = {ab.iloc[i]['Genome']: i+1 for i in range(10)}
    top30_mags = set(ab.head(30)['Genome'].tolist())

    # 组平均丰度
    ck_cols = [c for c in ['2_1','2_5','2_7'] if c in ab.columns]
    a_cols  = [c for c in ['8_1','8_2','8_3','8_4'] if c in ab.columns]
    b_cols  = [c for c in ['9_1','9_2','9_3'] if c in ab.columns]
    ab['CK_mean'] = ab[ck_cols].mean(axis=1) if ck_cols else 0
    ab['A_mean']  = ab[a_cols].mean(axis=1) if a_cols else 0
    ab['B_mean']  = ab[b_cols].mean(axis=1) if b_cols else 0
    ab_map = ab.set_index('Genome')[['CK_mean','A_mean','B_mean']].to_dict('index')

    # Keystone
    ks = pd.read_csv(args.ks, sep='\t')
    ks_mags = set(ks['MAG'].tolist())

    # 计算通路完整度矩阵
    all_mags = sorted(gtdb['MAG'].tolist())
    comp_data = []
    for mag in all_mags:
        kos = mag_kos.get(mag, set())
        row = {'MAG': mag}
        for pw in PATHWAY_ORDER:
            pw_kos = set(PATHWAY_KOS[pw])
            detected = kos & pw_kos
            row[pw] = len(detected) / len(pw_kos) * 100
        comp_data.append(row)
    df_comp = pd.DataFrame(comp_data)
    df_comp = df_comp.merge(gtdb[['MAG', 'Phylum', 'label']], on='MAG', how='left')
    df_comp['is_top10'] = df_comp['MAG'].isin(top10_rank.keys())
    df_comp['is_keystone'] = df_comp['MAG'].isin(ks_mags)
    df_comp['top_rank'] = df_comp['MAG'].map(top10_rank)
    df_comp['in_top30'] = df_comp['MAG'].isin(top30_mags)
    # 组平均丰度
    for col in ['CK_mean','A_mean','B_mean']:
        df_comp[col] = df_comp['MAG'].map(lambda m: ab_map.get(m, {}).get(col, 0))
    # 丰度排名（全部MAG）
    ab_rank = {ab.iloc[i]['Genome']: i+1 for i in range(len(ab))}
    df_comp['abund_rank'] = df_comp['MAG'].map(ab_rank)

    # 排序：按门→通路总完整度降序
    df_comp['total_comp'] = df_comp[PATHWAY_ORDER].sum(axis=1)
    phy_order = df_comp['Phylum'].value_counts().index.tolist()
    df_comp['phy_rank'] = df_comp['Phylum'].map({p: i for i, p in enumerate(phy_order)})
    df_comp = df_comp.sort_values(['phy_rank', 'total_comp'],
                                   ascending=[True, False]).reset_index(drop=True)

    return df_comp, phy_order, top10_rank, ks_mags, top30_mags


# ================================================================
# Fig3-5a: 全168 MAG通路完整度热图
# ================================================================
def draw_heatmap(df_comp, phy_order, style, outpath):
    fm_, fn_ = get_fonts(style)
    is_cn = (style == 'cn')
    n_mags = len(df_comp)
    n_pw = len(PATHWAY_ORDER)
    mat = df_comp[PATHWAY_ORDER].values  # 0~100

    # 配色：白→对应元素色
    pw_cmaps = {}
    for pw in PATHWAY_ORDER:
        elem = PATHWAY_ELEMENT[pw]
        ec = ELEM_COLORS[elem]
        pw_cmaps[pw] = LinearSegmentedColormap.from_list(
            pw, ['#FFFFFF', FUNC_BG[pw], ec])

    fig_h = n_mags * 0.13 + 6
    fig_w = n_pw * 0.55 + 10
    fig = plt.figure(figsize=(fig_w, fig_h))

    gs = gridspec.GridSpec(2, 5,
        width_ratios=[0.4, 4.2, n_pw * 0.55, 2.0, 2.8],
        height_ratios=[n_mags, 3.0], wspace=0.02, hspace=0.0)
    ax_phy = fig.add_subplot(gs[0, 0])
    ax_ylab = fig.add_subplot(gs[0, 1])
    ax_heat = fig.add_subplot(gs[0, 2])
    ax_abund = fig.add_subplot(gs[0, 3])
    ax_leg = fig.add_subplot(gs[0, 4])
    ax_xlab = fig.add_subplot(gs[1, 2])

    # ---- 热图 ----
    for j, pw in enumerate(PATHWAY_ORDER):
        cmap_pw = pw_cmaps[pw]
        for i in range(n_mags):
            val = mat[i, j]
            if val > 0:
                color = cmap_pw(val / 100.0)
                ax_heat.add_patch(Rectangle((j-0.45, i-0.45), 0.9, 0.9,
                                  facecolor=color, edgecolor='none', zorder=1))

    for i in range(n_mags + 1):
        ax_heat.axhline(y=i-0.5, color='white', lw=0.25, zorder=2)
    for j in range(n_pw + 1):
        ax_heat.axvline(x=j-0.5, color='white', lw=0.25, zorder=2)

    # 行边框标注
    for i, (_, row) in enumerate(df_comp.iterrows()):
        if row['is_top10'] and row['is_keystone']:
            bc, bw = BOTH_BC, 2.0
        elif row['is_top10']:
            bc, bw = TOP10_BC, 1.8
        elif row['is_keystone']:
            bc, bw = KS_BC, 1.8
        else:
            continue
        ax_heat.add_patch(Rectangle((-0.5, i-0.5), n_pw, 1,
                          facecolor='none', edgecolor=bc, linewidth=bw, zorder=4))

    ax_heat.set_yticks([]); ax_heat.set_xticks([])
    ax_heat.set_xlim(-0.5, n_pw - 0.5)
    ax_heat.set_ylim(n_mags - 0.5, -0.5)

    # 元素分隔线
    elem_starts = {}
    for j, pw in enumerate(PATHWAY_ORDER):
        e = PATHWAY_ELEMENT[pw]
        if e not in elem_starts: elem_starts[e] = j
    for elem, start in elem_starts.items():
        if start > 0:
            ax_heat.plot([start-0.5, start-0.5], [-0.5, n_mags-0.5],
                        color='#333', lw=1.5, zorder=3, clip_on=True)

    # 顶部元素色条
    for j, pw in enumerate(PATHWAY_ORDER):
        elem = PATHWAY_ELEMENT[pw]
        ax_heat.add_patch(Rectangle((j-0.5, -2.5), 1, 1.8,
                          facecolor=ELEM_COLORS[elem], edgecolor='none',
                          clip_on=False, alpha=0.85))
    elem_positions = {'As': [], 'N': [], 'S': [], 'Fe': []}
    for j, pw in enumerate(PATHWAY_ORDER):
        elem_positions[PATHWAY_ELEMENT[pw]].append(j)
    enames = {'As':'砷','N':'氮','S':'硫','Fe':'铁'} if is_cn else \
             {'As':'Arsenic','N':'Nitrogen','S':'Sulfur','Fe':'Iron'}
    for elem, positions in elem_positions.items():
        if positions:
            mid = (positions[0] + positions[-1]) / 2
            ax_heat.text(mid, -3.8, enames[elem], ha='center', va='center',
                        fontsize=8.5, fontweight='bold',
                        fontfamily=fm_ if is_cn else fn_, color=ELEM_COLORS[elem])

    # 门分组线
    prev_phy = None
    for i, (_, row) in enumerate(df_comp.iterrows()):
        if row['Phylum'] != prev_phy:
            if prev_phy is not None:
                ax_heat.axhline(y=i-0.5, color='#888', lw=0.5, zorder=3)
            prev_phy = row['Phylum']

    # ---- CK/A/B 丰度侧栏 ----
    abund_cols = ['CK_mean', 'A_mean', 'B_mean']
    abund_mat = df_comp[abund_cols].values
    ab_cmap, ab_norm, ab_ticks, ab_vmax = build_abundance_cmap(abund_mat)

    ax_abund.imshow(abund_mat, aspect='auto', cmap=ab_cmap, norm=ab_norm,
                    interpolation='nearest')
    ax_abund.set_xlim(-0.5, 2.5); ax_abund.set_ylim(n_mags-0.5, -0.5)
    ax_abund.set_yticks([])
    grp_labels = ['CK组','A组','B组'] if is_cn else ['CK','A','B']
    ax_abund.set_xticks([0, 1, 2])
    ax_abund.set_xticklabels(grp_labels, fontsize=6, fontfamily=fm_ if is_cn else fn_,
                              rotation=0, ha='center')
    ax_abund.tick_params(axis='x', left=False, labelleft=False, bottom=False, labelbottom=False)
    ax_abund.set_xticks([])
    # 组名放色块上方
    for j, (grp, col) in enumerate([('CK','#1c9cbd'),('A','#e3943d'),('B','#92181e')]):
        ax_abund.add_patch(Rectangle((j-0.5, -2.5), 1, 1.8,
                           facecolor=col, edgecolor='none', clip_on=False, alpha=0.8))
        grp_disp = f'{grp}组' if is_cn else grp
        ax_abund.text(j, -3.2, grp_disp, ha='center', va='center',
                     fontsize=6, fontweight='bold', fontfamily=fm_ if is_cn else fn_,
                     color=col, clip_on=False)
    for i in range(n_mags+1):
        ax_abund.axhline(y=i-0.5, color='white', lw=0.25)
    ax_abund.axvline(x=0.5, color='white', lw=0.3)
    ax_abund.axvline(x=1.5, color='white', lw=0.3)

    # ---- 门颜色条 ----
    for i, (_, row) in enumerate(df_comp.iterrows()):
        ax_phy.add_patch(Rectangle((0, i-0.5), 1, 1,
                         facecolor=PHYLUM_COLORS.get(row['Phylum'], '#BDC3C7'), edgecolor='none'))
    ax_phy.set_xlim(0, 1); ax_phy.set_ylim(n_mags-0.5, -0.5); ax_phy.axis('off')

    # ---- Y轴标签 ----
    ax_ylab.set_xlim(0, 1); ax_ylab.set_ylim(n_mags-0.5, -0.5); ax_ylab.axis('off')
    for i, (_, row) in enumerate(df_comp.iterrows()):
        color = PHYLUM_COLORS.get(row['Phylum'], '#666')
        is_sp = row['is_top10'] or row['is_keystone']
        fw = 'bold' if is_sp else 'normal'
        fs = 4.3 if is_sp else 3.8
        ax_ylab.text(0.98, i, row['label'], ha='right', va='center',
                    fontsize=fs, fontfamily=fn_, fontstyle='italic',
                    fontweight=fw, color=color)
        x_cur = 0.01
        if row['is_top10']:
            r = int(row['top_rank'])
            ax_ylab.plot(x_cur, i, marker='o', markersize=4.5, color=TOP10_BC,
                        markeredgewidth=0, zorder=5, clip_on=False)
            ax_ylab.text(x_cur+0.02, i, str(r), ha='left', va='center',
                        fontsize=3.5, fontfamily=fn_, fontweight='bold', color=TOP10_BC)
            x_cur += 0.055
        if row['is_keystone']:
            ax_ylab.plot(x_cur, i, marker='D', markersize=3, color=KS_BC,
                        markeredgewidth=0, zorder=5, clip_on=False)

    # ---- 底部X标签 + 功能色条 ----
    ax_xlab.set_xlim(-0.5, n_pw-0.5); ax_xlab.set_ylim(0, 3.5); ax_xlab.axis('off')
    for j, pw in enumerate(PATHWAY_ORDER):
        c = FUNC_BG.get(pw, '#DDD')
        ax_xlab.add_patch(Rectangle((j-0.5, 2.0), 1, 1.2,
                          facecolor=c, edgecolor='white', linewidth=0.3))
    for elem, start in elem_starts.items():
        if start > 0:
            ax_xlab.plot([start-0.5, start-0.5], [2.0, 3.2], color='#333', lw=1.5, clip_on=True)
    for j, pw in enumerate(PATHWAY_ORDER):
        disp = FUNC_CN.get(pw, pw) if is_cn else pw
        ax_xlab.text(j, 1.8, disp, ha='right', va='top', rotation=45,
                    fontsize=5.5, fontfamily=fm_ if is_cn else fn_, fontstyle='italic')

    # ---- 图例 ----
    ax_leg.axis('off')
    y_pos = 0.98
    # 门图例
    ttl = '门水平' if is_cn else 'Phylum'
    ax_leg.text(0.02, y_pos, ttl, fontsize=7, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_pos -= 0.012
    for phy in phy_order:
        y_pos -= 0.018
        if y_pos < 0.55:
            remaining = len(phy_order) - phy_order.index(phy)
            ax_leg.text(0.12, y_pos, f'+ {remaining} more', fontsize=3.5,
                       fontfamily=fn_, color='#888', transform=ax_leg.transAxes, va='center')
            break
        ax_leg.add_patch(Rectangle((0.02, y_pos-0.005), 0.05, 0.012,
                         facecolor=PHYLUM_COLORS.get(phy, '#BDC3C7'),
                         edgecolor='none', transform=ax_leg.transAxes, clip_on=False))
        ax_leg.text(0.09, y_pos, phy, fontsize=3.5, fontfamily=fn_,
                   transform=ax_leg.transAxes, va='center')

    # 完整度色标（从上到下：100%→0%）
    y_pos -= 0.03
    ttl2 = '通路完整度' if is_cn else 'Pathway completeness'
    ax_leg.text(0.02, y_pos, ttl2, fontsize=6.5, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_pos -= 0.015
    bar_x = 0.02; bar_w = 0.04; block_h = 0.02
    pct_labels = [(100,'100%'),(75,'75%'),(50,'50%'),(25,'25%'),(0,'0%')]
    for pct, lbl in pct_labels:
        grey = int(255 * (1 - pct/100))
        c = f'#{grey:02x}{grey:02x}{grey:02x}'
        ax_leg.add_patch(Rectangle((bar_x, y_pos-0.008), bar_w, block_h,
                         facecolor=c, edgecolor='none',
                         transform=ax_leg.transAxes, clip_on=False))
        ax_leg.text(bar_x + bar_w + 0.01, y_pos, lbl, fontsize=4.5, fontfamily=fn_,
                   transform=ax_leg.transAxes, va='center')
        y_pos -= block_h + 0.003

    # 组平均丰度色标（竖向，从上到下 = 高→低）
    y_pos -= 0.02
    ttl_ab = '组平均丰度' if is_cn else 'Group mean abundance'
    ax_leg.text(0.02, y_pos, ttl_ab, fontsize=6.5, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_pos -= 0.015
    # 用和热图侧栏同一个cmap/norm/ticks
    tick_vals_leg = ab_ticks
    n_blk = len(tick_vals_leg) - 1
    blk_h = 0.018
    sub_st = 8
    tick_display = list(reversed(tick_vals_leg))
    for b in range(n_blk):
        hi_v = tick_display[b]; lo_v = tick_display[b+1]
        for ss in range(sub_st):
            frac = ss / sub_st
            val = hi_v - (hi_v - lo_v) * frac
            c = ab_cmap(ab_norm(val))
            y_blk = y_pos - 0.005 - b * blk_h - frac * (blk_h)
            ax_leg.add_patch(Rectangle((bar_x, y_blk - blk_h/sub_st), bar_w, blk_h/sub_st,
                             facecolor=c, edgecolor='none',
                             transform=ax_leg.transAxes, clip_on=False))
    for ti, tv in enumerate(tick_display):
        y_tk = y_pos - 0.005 - ti * blk_h
        ax_leg.plot([bar_x+bar_w, bar_x+bar_w+0.01], [y_tk, y_tk],
                   color='black', lw=0.5, transform=ax_leg.transAxes, clip_on=False)
        ax_leg.text(bar_x+bar_w+0.015, y_tk, f'{tv:.2f}', fontsize=4, fontfamily=fn_,
                   transform=ax_leg.transAxes, va='center')
    ax_leg.add_patch(Rectangle((bar_x, y_pos - 0.005 - n_blk*blk_h), bar_w, n_blk*blk_h,
                     facecolor='none', edgecolor='black', linewidth=0.5,
                     transform=ax_leg.transAxes, clip_on=False))
    y_pos -= n_blk * blk_h + 0.025

    # 实验分组图例
    ttl_grp = '实验分组' if is_cn else 'Treatment group'
    ax_leg.text(0.02, y_pos, ttl_grp, fontsize=6.5, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    for grp, gname_cn, gname_en in [('CK','CK组','Group CK'),('A','A组','Group A'),('B','B组','Group B')]:
        y_pos -= 0.022
        ax_leg.add_patch(Rectangle((0.02, y_pos-0.006), 0.05, 0.014,
                         facecolor={'CK':'#1c9cbd','A':'#e3943d','B':'#92181e'}[grp],
                         edgecolor='none', transform=ax_leg.transAxes, clip_on=False))
        gname = gname_cn if is_cn else gname_en
        ax_leg.text(0.09, y_pos, gname, fontsize=4.5, fontfamily=fm_ if is_cn else fn_,
                   transform=ax_leg.transAxes, va='center')

    # 行标注
    y_pos -= 0.02
    ttl4 = '行标注' if is_cn else 'Row annotation'
    ax_leg.text(0.02, y_pos, ttl4, fontsize=6.5, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_pos -= 0.02
    ax_leg.add_patch(Rectangle((0.02, y_pos-0.007), 0.06, 0.016,
                     facecolor='none', edgecolor=TOP10_BC, linewidth=2,
                     transform=ax_leg.transAxes, clip_on=False))
    ax_leg.plot(0.10, y_pos, marker='o', markersize=5, color=TOP10_BC,
               markeredgewidth=0, transform=ax_leg.transAxes, clip_on=False)
    t10 = 'n  丰度排名' if is_cn else 'n  Abundance rank'
    ax_leg.text(0.13, y_pos, t10, fontsize=4.5, fontfamily=fm_ if is_cn else fn_,
               transform=ax_leg.transAxes, va='center')
    y_pos -= 0.025
    ax_leg.add_patch(Rectangle((0.02, y_pos-0.007), 0.06, 0.016,
                     facecolor='none', edgecolor=KS_BC, linewidth=2,
                     transform=ax_leg.transAxes, clip_on=False))
    ax_leg.plot(0.10, y_pos, marker='D', markersize=4, color=KS_BC,
               markeredgewidth=0, transform=ax_leg.transAxes, clip_on=False)
    ks_t = '关键物种' if is_cn else 'Keystone species'
    ax_leg.text(0.13, y_pos, ks_t, fontsize=4.5, fontfamily=fm_ if is_cn else fn_,
               transform=ax_leg.transAxes, va='center')
    y_pos -= 0.025
    ax_leg.add_patch(Rectangle((0.02, y_pos-0.007), 0.06, 0.016,
                     facecolor='none', edgecolor=BOTH_BC, linewidth=2,
                     transform=ax_leg.transAxes, clip_on=False))
    ax_leg.plot(0.10, y_pos, marker='o', markersize=5, color=TOP10_BC,
               markeredgewidth=0, transform=ax_leg.transAxes, clip_on=False)
    ax_leg.plot(0.125, y_pos, marker='D', markersize=4, color=KS_BC,
               markeredgewidth=0, transform=ax_leg.transAxes, clip_on=False)
    both_t = '两者兼具' if is_cn else 'Both'
    ax_leg.text(0.155, y_pos, both_t, fontsize=4.5, fontfamily=fm_ if is_cn else fn_,
               transform=ax_leg.transAxes, va='center')

    title = 'MAG通路完整度' if is_cn else 'MAG Pathway Completeness'
    fig.suptitle(title, fontsize=14, fontweight='bold',
                fontfamily=fm_ if is_cn else fn_, y=0.998)

    dpi = 300 if is_cn else 600
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight', facecolor='white')
    png = outpath.rsplit('.', 1)[0] + '.png'
    plt.savefig(png, dpi=min(dpi, 300), bbox_inches='tight', facecolor='white')
    print(f"  Saved: {outpath} + .png")
    plt.close()


# ================================================================
# Fig3-5b/c: 气泡图
# ================================================================
def draw_bubble(df_sub, fig_label, style, outpath):
    fm_, fn_ = get_fonts(style)
    is_cn = (style == 'cn')
    n_mags = len(df_sub)
    n_pw = len(PATHWAY_ORDER)
    mat = df_sub[PATHWAY_ORDER].values

    # Y轴标签
    labels_y = df_sub['label'].tolist()
    phyla = df_sub['Phylum'].tolist()

    # 行高根据数量调整
    row_h = max(0.35, 8.0 / n_mags)
    fig_h = n_mags * row_h + 3.5
    fig_w = n_pw * 0.6 + 9

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(2, 4,
        width_ratios=[0.3, 4.0, n_pw * 0.6, 2.0],
        height_ratios=[n_mags, 2.5], wspace=0.02, hspace=0.02)
    ax_phy = fig.add_subplot(gs[0, 0])
    ax_ylab = fig.add_subplot(gs[0, 1])
    ax_bub = fig.add_subplot(gs[0, 2])
    ax_abund = fig.add_subplot(gs[0, 3])
    ax_xlab = fig.add_subplot(gs[1, 2])

    # ---- 气泡 ----
    for j, pw in enumerate(PATHWAY_ORDER):
        elem = PATHWAY_ELEMENT[pw]
        color = ELEM_COLORS[elem]
        for i in range(n_mags):
            val = mat[i, j]
            if val > 0:
                size = val / 100.0 * 200 + 5  # 5~205
                ax_bub.scatter(j, i, s=size, c=color, alpha=0.7,
                              edgecolors='white', linewidths=0.3, zorder=2)

    ax_bub.set_xlim(-0.5, n_pw - 0.5)
    ax_bub.set_ylim(n_mags - 0.5, -0.5)
    ax_bub.set_xticks([]); ax_bub.set_yticks([])

    # 网格线
    for i in range(n_mags + 1):
        ax_bub.axhline(y=i-0.5, color='#eee', lw=0.3, zorder=0)
    for j in range(n_pw + 1):
        ax_bub.axvline(x=j-0.5, color='#eee', lw=0.3, zorder=0)

    # 元素分隔线
    elem_starts = {}
    for j, pw in enumerate(PATHWAY_ORDER):
        e = PATHWAY_ELEMENT[pw]
        if e not in elem_starts: elem_starts[e] = j
    for elem, start in elem_starts.items():
        if start > 0:
            ax_bub.axvline(x=start-0.5, color='#999', lw=1.0, zorder=1)

    # 顶部元素色条
    for j, pw in enumerate(PATHWAY_ORDER):
        elem = PATHWAY_ELEMENT[pw]
        ax_bub.add_patch(Rectangle((j-0.5, -1.5), 1, 1.0,
                         facecolor=ELEM_COLORS[elem], edgecolor='none',
                         clip_on=False, alpha=0.85))
    enames = {'As':'砷','N':'氮','S':'硫','Fe':'铁'} if is_cn else \
             {'As':'Arsenic','N':'Nitrogen','S':'Sulfur','Fe':'Iron'}
    elem_positions = {'As':[],'N':[],'S':[],'Fe':[]}
    for j, pw in enumerate(PATHWAY_ORDER):
        elem_positions[PATHWAY_ELEMENT[pw]].append(j)
    for elem, positions in elem_positions.items():
        if positions:
            mid = (positions[0]+positions[-1])/2
            ax_bub.text(mid, -2.3, enames[elem], ha='center', va='center',
                       fontsize=8, fontweight='bold',
                       fontfamily=fm_ if is_cn else fn_, color=ELEM_COLORS[elem])

    # ---- 门颜色条 ----
    for i, phy in enumerate(phyla):
        ax_phy.add_patch(Rectangle((0, i-0.5), 1, 1,
                         facecolor=PHYLUM_COLORS.get(phy, '#BDC3C7'), edgecolor='none'))
    ax_phy.set_xlim(0, 1); ax_phy.set_ylim(n_mags-0.5, -0.5); ax_phy.axis('off')

    # ---- Y轴标签 ----
    ax_ylab.set_xlim(0, 1); ax_ylab.set_ylim(n_mags-0.5, -0.5); ax_ylab.axis('off')
    for i in range(n_mags):
        color = PHYLUM_COLORS.get(phyla[i], '#666')
        ax_ylab.text(0.98, i, labels_y[i], ha='right', va='center',
                    fontsize=7, fontfamily=fn_, fontstyle='italic',
                    fontweight='bold', color=color)

    # ---- CK/A/B 丰度侧栏（气泡图右侧） ----
    abund_cols = ['CK_mean', 'A_mean', 'B_mean']
    abund_mat = df_sub[abund_cols].values
    ab_cmap, ab_norm, ab_ticks, ab_vmax = build_abundance_cmap(abund_mat)

    ax_abund.imshow(abund_mat, aspect='auto', cmap=ab_cmap, norm=ab_norm,
                    interpolation='nearest')
    ax_abund.set_xlim(-0.5, 2.5); ax_abund.set_ylim(n_mags-0.5, -0.5)
    ax_abund.set_yticks([])
    grp_lbl = ['CK组','A组','B组'] if is_cn else ['CK','A','B']
    ax_abund.set_xticks([0, 1, 2])
    ax_abund.set_xticklabels(grp_lbl, fontsize=7, fontfamily=fm_ if is_cn else fn_,
                              rotation=0, ha='center')
    ax_abund.tick_params(axis='x', left=False, labelleft=False, bottom=False, labelbottom=False)
    ax_abund.set_xticks([])
    # 组名放色块上方
    for jj, (grp, col) in enumerate([('CK','#1c9cbd'),('A','#e3943d'),('B','#92181e')]):
        ax_abund.add_patch(Rectangle((jj-0.5, -1.5), 1, 1.0,
                           facecolor=col, edgecolor='none', clip_on=False, alpha=0.8))
        grp_disp = f'{grp}组' if is_cn else grp
        ax_abund.text(jj, -2.0, grp_disp, ha='center', va='center',
                     fontsize=6, fontweight='bold', fontfamily=fm_ if is_cn else fn_,
                     color=col, clip_on=False)
    for ii in range(n_mags+1):
        ax_abund.axhline(y=ii-0.5, color='white', lw=0.3)
    ax_abund.axvline(x=0.5, color='white', lw=0.3)
    ax_abund.axvline(x=1.5, color='white', lw=0.3)

    # ---- 底部X标签 ----
    ax_xlab.set_xlim(-0.5, n_pw-0.5); ax_xlab.set_ylim(0, 2.5); ax_xlab.axis('off')
    for j, pw in enumerate(PATHWAY_ORDER):
        c = FUNC_BG.get(pw, '#DDD')
        ax_xlab.add_patch(Rectangle((j-0.5, 1.5), 1, 0.8,
                          facecolor=c, edgecolor='white', linewidth=0.3))
    for elem, start in elem_starts.items():
        if start > 0:
            ax_xlab.plot([start-0.5, start-0.5], [1.5, 2.3], color='#999', lw=1.0, clip_on=True)
    for j, pw in enumerate(PATHWAY_ORDER):
        disp = FUNC_CN.get(pw, pw) if is_cn else pw
        ax_xlab.text(j, 1.3, disp, ha='right', va='top', rotation=45,
                    fontsize=6, fontfamily=fm_ if is_cn else fn_, fontstyle='italic')

    # ---- 图例（用ax_abund右侧空间） ----
    # 完整度气泡大小
    ax_abund.text(3.5, -0.8, '完整度' if is_cn else 'Completeness',
                 fontsize=7, fontweight='bold', fontfamily=fm_ if is_cn else fn_,
                 clip_on=False, va='center')
    for pct in [25, 50, 75, 100]:
        size = pct / 100.0 * 200 + 5
        y_leg = [25,50,75,100].index(pct) * 1.2 + 0.3
        ax_abund.scatter(3.8, y_leg, s=size, c='#888', alpha=0.6,
                        edgecolors='white', linewidths=0.3, clip_on=False, zorder=10)
        ax_abund.text(4.5, y_leg, f'{pct}%', fontsize=6, fontfamily=fn_,
                     va='center', clip_on=False)

    # 元素颜色
    y_leg_e = 5.5
    ax_abund.text(3.5, y_leg_e, '元素' if is_cn else 'Element',
                 fontsize=7, fontweight='bold', fontfamily=fm_ if is_cn else fn_,
                 clip_on=False, va='center')
    for idx, (elem, ename) in enumerate([('As','Arsenic'),('N','Nitrogen'),
                                          ('S','Sulfur'),('Fe','Iron')]):
        y_e = y_leg_e + (idx + 1) * 1.0
        edisp = {'As':'砷','N':'氮','S':'硫','Fe':'铁'}[elem] if is_cn else ename
        ax_abund.scatter(3.8, y_e, s=80, c=ELEM_COLORS[elem], alpha=0.7,
                        edgecolors='white', linewidths=0.3, clip_on=False, zorder=10)
        ax_abund.text(4.5, y_e, edisp, fontsize=6,
                     fontfamily=fm_ if is_cn else fn_, va='center', clip_on=False)

    # 门水平图例（不限制数量，用紧凑间距）
    sub_phyla = list(dict.fromkeys(phyla))
    y_leg_p = y_leg_e + 5.5
    ax_abund.text(3.5, y_leg_p, '门水平' if is_cn else 'Phylum',
                 fontsize=7, fontweight='bold', fontfamily=fm_ if is_cn else fn_,
                 clip_on=False, va='center')
    for idx_p, phy in enumerate(sub_phyla):
        y_p = y_leg_p + (idx_p + 1) * 0.8
        ax_abund.add_patch(Rectangle((3.5, y_p - 0.25), 0.5, 0.5,
                           facecolor=PHYLUM_COLORS.get(phy, '#BDC3C7'),
                           edgecolor='none', clip_on=False, zorder=10))
        ax_abund.text(4.3, y_p, phy, fontsize=5, fontfamily=fn_,
                     va='center', clip_on=False)

    # 组平均丰度色标（竖向，从上到下 = 高→低）
    y_leg_ab = y_leg_p + len(sub_phyla) * 0.8 + 1.5
    ax_abund.text(3.5, y_leg_ab, '组平均丰度' if is_cn else 'Group mean abundance',
                 fontsize=7, fontweight='bold', fontfamily=fm_ if is_cn else fn_,
                 clip_on=False, va='center')
    tick_vals_leg = ab_ticks
    n_blocks = len(tick_vals_leg) - 1
    blk_h = 0.7
    sub_steps = 8
    tick_display = list(reversed(tick_vals_leg))
    for b in range(n_blocks):
        hi_v = tick_display[b]; lo_v = tick_display[b+1]
        for ss in range(sub_steps):
            frac = ss / sub_steps
            val = hi_v - (hi_v - lo_v) * frac
            c = ab_cmap(ab_norm(val))
            y_blk = y_leg_ab + 0.5 + b * blk_h + frac * blk_h
            ax_abund.add_patch(Rectangle((3.5, y_blk), 0.5, blk_h / sub_steps,
                               facecolor=c, edgecolor='none', clip_on=False, zorder=10))
    for ti, tv in enumerate(tick_display):
        y_tk = y_leg_ab + 0.5 + ti * blk_h
        ax_abund.plot([4.0, 4.15], [y_tk, y_tk], color='black', lw=0.5,
                     clip_on=False, zorder=10)
        ax_abund.text(4.3, y_tk, f'{tv:.2f}', fontsize=5, fontfamily=fn_,
                     va='center', clip_on=False)
    ax_abund.add_patch(Rectangle((3.5, y_leg_ab + 0.5), 0.5, n_blocks * blk_h,
                       facecolor='none', edgecolor='black', linewidth=0.5,
                       clip_on=False, zorder=11))

    # 实验分组图例
    y_leg_grp = y_leg_ab + 0.5 + n_blocks * blk_h + 1.5
    ax_abund.text(3.5, y_leg_grp, '实验分组' if is_cn else 'Treatment group',
                 fontsize=7, fontweight='bold', fontfamily=fm_ if is_cn else fn_,
                 clip_on=False, va='center')
    for ig, (grp, gname_cn, gname_en) in enumerate([('CK','CK组','Group CK'),
                                                      ('A','A组','Group A'),
                                                      ('B','B组','Group B')]):
        y_g = y_leg_grp + (ig + 1) * 0.9
        ax_abund.add_patch(Rectangle((3.5, y_g - 0.25), 0.5, 0.5,
                           facecolor={'CK':'#1c9cbd','A':'#e3943d','B':'#92181e'}[grp],
                           edgecolor='none', clip_on=False, zorder=10))
        gname = gname_cn if is_cn else gname_en
        ax_abund.text(4.3, y_g, gname, fontsize=5, fontfamily=fm_ if is_cn else fn_,
                     va='center', clip_on=False)

    title = fig_label + ('通路完整度' if is_cn else ' Pathway Completeness')
    fig.suptitle(title, fontsize=13, fontweight='bold',
                fontfamily=fm_ if is_cn else fn_, y=0.99)

    dpi = 300 if is_cn else 600
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight', facecolor='white')
    png = outpath.rsplit('.', 1)[0] + '.png'
    plt.savefig(png, dpi=min(dpi, 300), bbox_inches='tight', facecolor='white')
    print(f"  Saved: {outpath} + .png")
    plt.close()


# ================================================================
# 数据导出
# ================================================================
def export_stats(df_comp, stat_dir):
    os.makedirs(stat_dir, exist_ok=True)

    # 1. 完整度矩阵（含组丰度和富集方向）
    out = df_comp[['MAG','Phylum','label','is_top10','top_rank','is_keystone',
                    'CK_mean','A_mean','B_mean'] + PATHWAY_ORDER].copy()
    out['total_completeness'] = df_comp[PATHWAY_ORDER].sum(axis=1)
    # 计算富集方向
    def get_enriched(row):
        vals = {'CK': row['CK_mean'], 'A': row['A_mean'], 'B': row['B_mean']}
        enriched = []
        for grp, val in vals.items():
            others = [v for k, v in vals.items() if k != grp]
            if val > 1.5 * np.mean(others) and val > 0.05:
                enriched.append(grp)
        return ','.join(enriched) if enriched else '-'
    out['enriched_in'] = out.apply(get_enriched, axis=1)
    out.to_csv(os.path.join(stat_dir, 'Fig3-5_pathway_completeness_matrix.tsv'),
               sep='\t', index=False)

    # 2. 通路汇总统计
    with open(os.path.join(stat_dir, 'Fig3-5_pathway_summary.txt'), 'w') as f:
        f.write("=== Pathway Completeness Summary ===\n\n")
        f.write(f"Total MAGs: {len(df_comp)}\n")
        f.write(f"Pathways: {len(PATHWAY_ORDER)}\n\n")
        f.write(f"{'Pathway':<25} {'Element':<5} {'KOs':>4} "
                f"{'Mean%':>7} {'Med%':>6} {'MAGs>0':>7} {'MAGs=100':>9}\n")
        f.write("-" * 70 + "\n")
        for pw in PATHWAY_ORDER:
            vals = df_comp[pw].values
            elem = PATHWAY_ELEMENT[pw]
            n_kos = len(PATHWAY_KOS[pw])
            mean_v = vals.mean()
            med_v = np.median(vals)
            n_pos = (vals > 0).sum()
            n_full = (vals >= 100).sum()
            f.write(f"{pw:<25} {elem:<5} {n_kos:>4} "
                    f"{mean_v:>6.1f}% {med_v:>5.0f}% {n_pos:>7} {n_full:>9}\n")

    # 3. Top10通路完整度（含组丰度）
    top10 = df_comp[df_comp['is_top10']].sort_values('top_rank')
    with open(os.path.join(stat_dir, 'Fig3-5_Top10_completeness.txt'), 'w') as f:
        f.write("=== Top 10 MAGs Pathway Completeness + Group Abundance ===\n\n")
        header = f"{'Rank':<5} {'MAG':<15} {'Species':<30} {'CK':>6} {'A':>6} {'B':>6} {'Enrich':>8}"
        for pw in PATHWAY_ORDER:
            short = pw[:12]
            header += f" {short:>12}"
        f.write(header + "\n")
        f.write("-" * (75 + 12 * len(PATHWAY_ORDER)) + "\n")
        for _, row in top10.iterrows():
            enriched = get_enriched(row)
            line = f"{int(row['top_rank']):<5} {row['MAG']:<15} {row['label']:<30} "
            line += f"{row['CK_mean']:>6.3f} {row['A_mean']:>6.3f} {row['B_mean']:>6.3f} {enriched:>8}"
            for pw in PATHWAY_ORDER:
                line += f" {row[pw]:>11.0f}%"
            f.write(line + "\n")

    # 4. Keystone通路完整度（含组丰度）
    ks_df = df_comp[df_comp['is_keystone']]
    with open(os.path.join(stat_dir, 'Fig3-5_Keystone_completeness.txt'), 'w') as f:
        f.write("=== Keystone Species Pathway Completeness + Group Abundance ===\n\n")
        header = f"{'MAG':<15} {'Species':<30} {'CK':>6} {'A':>6} {'B':>6} {'Enrich':>8}"
        for pw in PATHWAY_ORDER:
            short = pw[:12]
            header += f" {short:>12}"
        f.write(header + "\n")
        f.write("-" * (70 + 12 * len(PATHWAY_ORDER)) + "\n")
        for _, row in ks_df.iterrows():
            enriched = get_enriched(row)
            line = f"{row['MAG']:<15} {row['label']:<30} "
            line += f"{row['CK_mean']:>6.3f} {row['A_mean']:>6.3f} {row['B_mean']:>6.3f} {enriched:>8}"
            for pw in PATHWAY_ORDER:
                line += f" {row[pw]:>11.0f}%"
            f.write(line + "\n")

    print(f"  Stats exported to: {stat_dir}")


# ================================================================
# Main
# ================================================================
def main():
    parser = argparse.ArgumentParser(description='Fig3-5 Pathway Completeness')
    parser.add_argument('--ko',     required=True, help='kegg_target_only.tsv')
    parser.add_argument('--bac',    required=True, help='GTDB bac120 summary')
    parser.add_argument('--ar',     required=True, help='GTDB ar53 summary')
    parser.add_argument('--abund',  required=True, help='CoverM abundance.tsv')
    parser.add_argument('--ks',     required=True, help='Keystone species list')
    parser.add_argument('--outdir', default='figures')
    parser.add_argument('--statdir',default='data/processed')
    args = parser.parse_args()

    print(f"Fonts: sans={F_SANS}, serif={F_SERIF}, cn={F_CN}")

    df_comp, phy_order, top10_rank, ks_mags, top30_mags = load_data(args)
    print(f"Data: {len(df_comp)} MAGs × {len(PATHWAY_ORDER)} pathways")

    en_dir = os.path.join(args.outdir, 'paper_en')
    cn_dir = os.path.join(args.outdir, 'thesis_cn')
    os.makedirs(en_dir, exist_ok=True)
    os.makedirs(cn_dir, exist_ok=True)

    # Fig3-5a: 全168 MAG热图
    print("\n=== Fig3-5a: Full 168 MAGs heatmap ===")
    draw_heatmap(df_comp, phy_order, 'en',
                 os.path.join(en_dir, 'Fig3-5a_pathway_completeness_EN.pdf'))
    draw_heatmap(df_comp, phy_order, 'cn',
                 os.path.join(cn_dir, 'Fig3-5a_pathway_completeness_CN.pdf'))

    # Fig3-5b: Top30气泡图
    print("\n=== Fig3-5b: Top30 bubble ===")
    df_top30 = df_comp[df_comp['in_top30']].sort_values('abund_rank').reset_index(drop=True)
    draw_bubble(df_top30, 'Top 30 MAGs', 'en',
                os.path.join(en_dir, 'Fig3-5b_pathway_top30_bubble_EN.pdf'))
    draw_bubble(df_top30, '丰度Top30 MAGs', 'cn',
                os.path.join(cn_dir, 'Fig3-5b_pathway_top30_bubble_CN.pdf'))

    # Fig3-5c: Keystone气泡图
    print("\n=== Fig3-5c: Keystone bubble ===")
    df_ks = df_comp[df_comp['is_keystone']].reset_index(drop=True)
    draw_bubble(df_ks, 'Keystone Species', 'en',
                os.path.join(en_dir, 'Fig3-5c_pathway_keystone_bubble_EN.pdf'))
    draw_bubble(df_ks, '关键物种', 'cn',
                os.path.join(cn_dir, 'Fig3-5c_pathway_keystone_bubble_CN.pdf'))

    # 数据导出
    print("\n=== Exporting stats ===")
    export_stats(df_comp, args.statdir)

    print("\nAll done!")

if __name__ == '__main__':
    main()
