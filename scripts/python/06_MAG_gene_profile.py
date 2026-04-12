#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_MAG_gene_profile.py — Fig3-3 MAG元素循环基因谱热图
输出:
  figures/thesis_cn/Fig3-3_MAG_gene_profile_CN.pdf  (300dpi, SimSun+TNR)
  figures/paper_en/Fig3-3_MAG_gene_profile_EN.pdf   (600dpi, Arial)
  figures/paper_en/Fig3-3_MAG_gene_profile_EN.tiff  (600dpi)

用法:
  python3 scripts/python/06_MAG_gene_profile.py \\
      --itol  data/raw/07_element_KO_heatmap.txt \\
      --ko    data/raw/kegg_target_only.tsv \\
      --bac   data/raw/tax_bac120_summary.tsv \\
      --ar    data/raw/tax_ar53_summary.tsv \\
      --abund data/raw/abundance.tsv \\
      --ks    data/raw/keystone_species.txt \\
      --outdir figures
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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec

# ================================================================
# 字体设置
# ================================================================
def setup_fonts():
    """检测可用字体，返回 (sans, serif, cn) 名称"""
    avail = {f.name for f in fm.fontManager.ttflist}
    # SCI: Arial > Liberation Sans > DejaVu Sans
    sans = 'DejaVu Sans'
    for c in ['Arial', 'Helvetica', 'Liberation Sans']:
        if c in avail: sans = c; break
    # Thesis数字/英文: Times New Roman > Liberation Serif > DejaVu Serif
    serif = 'DejaVu Serif'
    for c in ['Times New Roman', 'Liberation Serif']:
        if c in avail: serif = c; break
    # Thesis中文: SimSun > 宋体 > WenQuanYi
    cn = serif  # fallback
    for c in ['SimSun', '宋体', 'WenQuanYi Micro Hei', 'Noto Sans CJK JP']:
        if c in avail: cn = c; break
    return sans, serif, cn

F_SANS, F_SERIF, F_CN = setup_fonts()

# mathtext字体配置：确保下标等用sans-serif而非默认serif斜体
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = F_SANS
plt.rcParams['mathtext.it'] = F_SANS
plt.rcParams['axes.unicode_minus'] = False

def get_fonts(style):
    """返回 (main_font, num_font) — main用于标签，num用于数字"""
    if style == 'cn':
        return F_CN, F_SERIF   # 中文正文宋体, 数字TNR
    else:
        return F_SANS, F_SANS  # 英文全用Arial

# ================================================================
# 数据配置
# ================================================================
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
ELEM_FILL = {'As':'#C0392B','N':'#2C3E87','S':'#1E8449','Fe':'#B7600B'}
ELEM_CMAPS = {
    'As': LinearSegmentedColormap.from_list('as',['#FFFFFF','#F5B7B1','#C0392B','#7B241C']),
    'N':  LinearSegmentedColormap.from_list('n', ['#FFFFFF','#AEB6BF','#2C3E87','#1A237E']),
    'S':  LinearSegmentedColormap.from_list('s', ['#FFFFFF','#A9DFBF','#1E8449','#0B5345']),
    'Fe': LinearSegmentedColormap.from_list('fe',['#FFFFFF','#F5CBA7','#B7600B','#6E3507']),
}
# Scheme2 Enhanced 功能分类色
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
KO_FUNC = {
    'K00537':('arsC-grx','Arsenate reduction'),'K03741':('arsC-trx','Arsenate reduction'),
    'K22547':('HAC1','Arsenate reduction'),'K18701':('arsC-myc','Arsenate reduction'),
    'K08355':('aoxA','Arsenite oxidation'),'K08356':('aoxB','Arsenite oxidation'),
    'K28466':('arrA','Resp. arsenate red.'),'K28467':('arrB','Resp. arsenate red.'),
    'K03893':('arsB','As transport/detox'),'K03325':('ACR3','As transport/detox'),
    'K18064':('ACR2','As transport/detox'),'K01551':('arsA','As transport/detox'),
    'K25223':('arsJ','As transport/detox'),'K11811':('arsH','As transport/detox'),
    'K07755':('AS3MT','As methylation'),'K25224':('gapdh','As methylation'),
    'K03892':('arsR','As regulation'),
    'K00370':('narG','Nitrate reduction'),'K00371':('narH','Nitrate reduction'),
    'K00374':('narI','Nitrate reduction'),'K02567':('napA','Nitrate reduction'),
    'K02568':('napB','Nitrate reduction'),'K00367':('narB','Nitrate reduction'),
    'K00368':('nirK','Nitrite reduction'),'K15864':('nirS','Nitrite reduction'),
    'K00366':('nirA','Nitrite reduction'),'K15877':('CYP55','Nitrite reduction'),
    'K04561':('norB','NO reduction'),'K00376':('nosZ','N2O reduction'),
    'K10944':('amoA','Ammonia oxidation'),'K10945':('amoB','Ammonia oxidation'),
    'K10946':('amoC','Ammonia oxidation'),
    'K02586':('nifD','N fixation'),'K02591':('nifK','N fixation'),
    'K00380':('cysJ','Assim. sulfate red.'),'K00381':('cysI','Assim. sulfate red.'),
    'K00390':('cysH','Assim. sulfate red.'),'K00392':('sir','Assim. sulfate red.'),
    'K00394':('aprA','Dissim. sulfate red.'),'K00395':('aprB','Dissim. sulfate red.'),
    'K11180':('dsrA','Dissim. sulfate red.'),'K11181':('dsrB','Dissim. sulfate red.'),
    'K17218':('sqr','Sulfide oxidation'),'K17222':('soxA','Sulfide oxidation'),
    'K17223':('soxX','Sulfide oxidation'),'K17224':('soxB','Sulfide oxidation'),
    'K17226':('soxY','Sulfide oxidation'),'K17227':('soxZ','Sulfide oxidation'),
    'K01011':('TST','Thiosulfate metab.'),
    'K02013':('ABC.FEV.A','Fe transport'),'K02014':('TC.FEV.OM','Fe transport'),
    'K02015':('ABC.FEV.P','Fe transport'),'K02016':('ABC.FEV.S','Fe transport'),
    'K02012':('afuA','Fe transport'),'K07243':('FTR','Fe transport'),
    'K03711':('fur','Fe uptake regulation'),'K03832':('tonB','Fe uptake regulation'),
}
# 中文功能名
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

TOP10_BC = '#E67E22'
KS_BC = '#8E44AD'
BOTH_BC = '#1ABC9C'


def build_abundance_cmap(data_values):
    """根据实际丰度数据构建三段非线性配色"""
    from matplotlib.colors import BoundaryNorm
    vmax = max(np.max(data_values), 0.01)
    if vmax <= 0.3:
        bounds = list(np.arange(0, vmax + 0.01, max(vmax/10, 0.01)))
    else:
        bounds = list(np.arange(0, min(0.2, vmax), 0.025))
        if vmax > 0.2:
            bounds += list(np.arange(0.2, min(0.5, vmax), 0.05))
        if vmax > 0.5:
            bounds += list(np.arange(0.5, min(1.0, vmax), 0.1))
        if vmax > 1.0:
            bounds += list(np.arange(1.0, vmax + 0.01, 0.5))
    bounds = sorted(set([round(b, 3) for b in bounds]))
    if len(bounds) < 2: bounds = [0, round(vmax, 3)]
    if bounds[-1] < vmax: bounds.append(round(vmax, 3))

    n_lo = max(len([b for b in bounds if b <= 0.2]) - 1, 0)
    n_mid = max(len([b for b in bounds if 0.2 < b <= 0.5]), 0)
    n_hi = max(len([b for b in bounds if b > 0.5]), 0)
    seg_lo = plt.cm.Blues(np.linspace(0.08, 0.55, max(n_lo+1, 2)))
    seg_mid = plt.cm.YlGn(np.linspace(0.18, 0.65, max(n_mid+1, 2)))
    seg_hi = plt.cm.YlOrRd(np.linspace(0.40, 0.95, max(n_hi+1, 2)))
    colors = list(seg_lo[:max(n_lo+1, 1)])
    if n_mid > 0: colors += list(seg_mid[1:n_mid+1])
    if n_hi > 0: colors += list(seg_hi[1:n_hi+1])
    if len(colors) < 2: colors = [seg_lo[0], seg_lo[-1]]

    cmap = LinearSegmentedColormap.from_list('ab_auto', colors, N=256)
    norm = BoundaryNorm(bounds, ncolors=256)
    candidates = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0]
    tick_vals = sorted(set([v for v in candidates if v <= vmax] + [round(vmax, 2)]))
    return cmap, norm, tick_vals, vmax

# ================================================================
# 数据加载
# ================================================================
def load_data(args):
    # iTOL热图 -> presence/absence
    lines = open(args.itol).readlines()
    ko_labels = []
    for l in lines:
        if l.startswith('FIELD_LABELS'):
            ko_labels = l.strip().split(',')[1:]
            break
    data_rows = []
    in_data = False
    for l in lines:
        if l.strip() == 'DATA':
            in_data = True; continue
        if in_data and l.strip():
            parts = l.strip().split(',')
            data_rows.append([parts[0]] + [int(x) for x in parts[1:]])
    df_ko = pd.DataFrame(data_rows, columns=['MAG'] + ko_labels)

    # 拷贝数
    raw = pd.read_csv(args.ko, sep='\t')
    rows_exp = []
    for _, r in raw.iterrows():
        for ko in str(r['KEGG_ko']).split(','):
            rows_exp.append((r['MAG'], ko.strip()))
    df_exp = pd.DataFrame(rows_exp, columns=['MAG','KO'])
    tmap = {}
    for ko in ko_labels:
        k_id = ko.split('(')[1].rstrip(')')
        tmap[f'ko:{k_id}'] = ko
    df_t = df_exp[df_exp['KO'].isin(tmap.keys())]
    cn = df_t.groupby(['MAG','KO']).size().reset_index(name='copies')
    cn['ko_label'] = cn['KO'].map(tmap)
    copy_pivot = cn.pivot_table(index='MAG', columns='ko_label', values='copies', fill_value=0)

    # GTDB分类
    bac = pd.read_csv(args.bac, sep='\t')
    ar = pd.read_csv(args.ar, sep='\t')
    gtdb = pd.concat([bac[['user_genome','classification']],
                       ar[['user_genome','classification']]], ignore_index=True)
    gtdb = gtdb.rename(columns={'user_genome':'MAG'})
    def gl(cls, pfx):
        for p in cls.split(';'):
            if p.startswith(pfx): return p.replace(pfx,'')
        return ''
    gtdb['Phylum'] = gtdb['classification'].apply(lambda x: gl(x,'p__'))
    gtdb['Family'] = gtdb['classification'].apply(lambda x: gl(x,'f__'))
    gtdb['Genus'] = gtdb['classification'].apply(lambda x: gl(x,'g__'))
    gtdb['Species'] = gtdb['classification'].apply(lambda x: gl(x,'s__'))
    def mk_label(row):
        n = row['MAG'].replace('Mx_All_','')
        g, s, f = row['Genus'], row['Species'], row['Family']
        if g:
            if s:
                ss = s.replace(g+' ','') if g in s else s
                return f"{g} {ss} [{n}]"
            return f"{g} sp. [{n}]"
        elif f:
            return f"{f} (f) [{n}]"
        return f"Unclassified [{n}]"
    gtdb['label'] = gtdb.apply(mk_label, axis=1)
    df = df_ko.merge(gtdb[['MAG','Phylum','label']], on='MAG', how='left')

    # Top10 + Keystone
    ab = pd.read_csv(args.abund, sep='\t')
    ab = ab[ab['Genome'] != 'unmapped']
    scols = [c for c in ab.columns if c != 'Genome']
    ab['mean_all'] = ab[scols].mean(axis=1)
    ab = ab.sort_values('mean_all', ascending=False).reset_index(drop=True)
    top10_rank = {ab.iloc[i]['Genome']: i+1 for i in range(10)}
    ks = pd.read_csv(args.ks, sep='\t')
    ks_mags = set(ks['MAG'].tolist())
    df['top_rank'] = df['MAG'].map(top10_rank)
    df['is_top10'] = df['MAG'].isin(top10_rank.keys())
    df['is_keystone'] = df['MAG'].isin(ks_mags)

    # 组平均丰度
    ck_cols = [c for c in ['2_1','2_5','2_7'] if c in ab.columns]
    a_cols  = [c for c in ['8_1','8_2','8_3','8_4'] if c in ab.columns]
    b_cols  = [c for c in ['9_1','9_2','9_3'] if c in ab.columns]
    ab['CK_mean'] = ab[ck_cols].mean(axis=1) if ck_cols else 0
    ab['A_mean']  = ab[a_cols].mean(axis=1) if a_cols else 0
    ab['B_mean']  = ab[b_cols].mean(axis=1) if b_cols else 0
    ab_map = ab.set_index('Genome')[['CK_mean','A_mean','B_mean']].to_dict('index')
    for col in ['CK_mean','A_mean','B_mean']:
        df[col] = df['MAG'].map(lambda m: ab_map.get(m, {}).get(col, 0))

    # 元素映射
    ko_element = {}
    for i, ko in enumerate(ko_labels):
        if i < 17: ko_element[ko] = 'As'
        elif i < 34: ko_element[ko] = 'N'
        elif i < 49: ko_element[ko] = 'S'
        else: ko_element[ko] = 'Fe'

    # 排序
    df['gene_count'] = df[ko_labels].sum(axis=1)
    phy_order = df['Phylum'].value_counts().index.tolist()
    df['phy_rank'] = df['Phylum'].map({p:i for i,p in enumerate(phy_order)})
    df = df.sort_values(['phy_rank','gene_count'], ascending=[True,False]).reset_index(drop=True)

    # 过滤全零KO
    ko_sums = df[ko_labels].sum()
    active_kos = ko_sums[ko_sums > 0].index.tolist()

    # 拷贝数矩阵
    mat = np.zeros((len(df), len(active_kos)))
    for j, ko in enumerate(active_kos):
        for i, (_,row) in enumerate(df.iterrows()):
            mag = row['MAG']
            if mag in copy_pivot.index and ko in copy_pivot.columns:
                mat[i,j] = copy_pivot.loc[mag, ko]
    mat_log = np.log1p(mat)

    return df, active_kos, ko_labels, ko_element, mat, mat_log, phy_order

# ================================================================
# 绘图
# ================================================================
def draw_figure(df, active_kos, ko_labels, ko_element, mat, mat_log, phy_order,
                style, outpath):
    fm_, fn_ = get_fonts(style)
    is_cn = (style == 'cn')
    max_log = mat_log.max()
    n_mags, n_kos = mat.shape

    # X轴标签
    ko_display = []
    ko_funcs = []
    for ko in active_kos:
        k_id = ko.split('(')[1].rstrip(')')
        if k_id in KO_FUNC:
            gn, func = KO_FUNC[k_id]
            ko_display.append(gn); ko_funcs.append(func)
        else:
            ko_display.append(ko.split('(')[0]); ko_funcs.append('Unknown')

    # 尺寸（加高以容纳完整图例：23门+4拷贝+3行标+8丰度+3分组+18功能≈59项）
    n_legend_items = len(phy_order) + 4 + 3 + 8 + 3 + 18
    min_legend_h = n_legend_items * 0.35 + 8
    fig_h = max(n_mags * 0.13 + 6, min_legend_h)
    fig_w = n_kos * 0.30 + 10
    fig = plt.figure(figsize=(fig_w, fig_h))

    # [phylum | Y labels | heatmap | abundance sidebar | legend]
    gs = gridspec.GridSpec(2, 5,
        width_ratios=[0.4, 4.2, n_kos*0.62, 2.0, 2.8],
        height_ratios=[n_mags, 3.5], wspace=0.02, hspace=0.0)
    ax_phy  = fig.add_subplot(gs[0,0])
    ax_ylab = fig.add_subplot(gs[0,1])
    ax_heat = fig.add_subplot(gs[0,2])
    ax_abund = fig.add_subplot(gs[0,3])
    ax_leg  = fig.add_subplot(gs[0,4])
    ax_func = fig.add_subplot(gs[1,2])

    # ---- 热图 ----
    for j, ko in enumerate(active_kos):
        elem = ko_element[ko]
        cmap_e = ELEM_CMAPS[elem]
        for i in range(n_mags):
            val = mat_log[i,j]
            if val > 0:
                color = cmap_e(min(val / max_log * 0.85 + 0.15, 1.0))
                ax_heat.add_patch(Rectangle((j-0.45, i-0.45), 0.9, 0.9,
                                  facecolor=color, edgecolor='none', zorder=1))
    # 网格
    for i in range(n_mags+1):
        ax_heat.axhline(y=i-0.5, color='white', lw=0.25, zorder=2)
    for j in range(n_kos+1):
        ax_heat.axvline(x=j-0.5, color='white', lw=0.25, zorder=2)

    # 行边框标注
    for i, (_,row) in enumerate(df.iterrows()):
        if row['is_top10'] and row['is_keystone']:
            bc, bw = BOTH_BC, 2.0
        elif row['is_top10']:
            bc, bw = TOP10_BC, 1.8
        elif row['is_keystone']:
            bc, bw = KS_BC, 1.8
        else:
            continue
        ax_heat.add_patch(Rectangle((-0.5, i-0.5), n_kos, 1,
                          facecolor='none', edgecolor=bc, linewidth=bw, zorder=4))

    ax_heat.set_yticks([]); ax_heat.set_xticks([])
    ax_heat.set_xlim(-0.5, n_kos-0.5)
    ax_heat.set_ylim(n_mags-0.5, -0.5)

    # 元素分隔线
    elem_starts = {}
    for j, ko in enumerate(active_kos):
        e = ko_element[ko]
        if e not in elem_starts: elem_starts[e] = j
    for elem, start in elem_starts.items():
        if start > 0:
            ax_heat.plot([start-0.5, start-0.5], [-0.5, n_mags-0.5],
                        color='#333', lw=1.5, zorder=3, clip_on=True)

    # 顶部元素色条
    for j, ko in enumerate(active_kos):
        elem = ko_element[ko]
        ax_heat.add_patch(Rectangle((j-0.5, -2.5), 1, 1.8,
                          facecolor=ELEM_FILL[elem], edgecolor='none', clip_on=False, alpha=0.85))
    elem_positions = {'As':[],'N':[],'S':[],'Fe':[]}
    for j, ko in enumerate(active_kos):
        elem_positions[ko_element[ko]].append(j)
    elem_names_en = {'As':'Arsenic','N':'Nitrogen','S':'Sulfur','Fe':'Iron'}
    elem_names_cn = {'As':'砷','N':'氮','S':'硫','Fe':'铁'}
    enames = elem_names_cn if is_cn else elem_names_en
    for elem, positions in elem_positions.items():
        if positions:
            mid = (positions[0]+positions[-1])/2
            ax_heat.text(mid, -3.8, enames[elem], ha='center', va='center',
                        fontsize=8.5, fontweight='bold',
                        fontfamily=fm_ if is_cn else fn_, color=ELEM_FILL[elem])

    # 门分组线
    prev_phy = None
    for i, (_,row) in enumerate(df.iterrows()):
        if row['Phylum'] != prev_phy:
            if prev_phy is not None:
                ax_heat.axhline(y=i-0.5, color='#888', lw=0.5, zorder=3)
            prev_phy = row['Phylum']

    # ---- CK/A/B 丰度侧栏 ----
    abund_cols = ['CK_mean', 'A_mean', 'B_mean']
    abund_mat = df[abund_cols].values
    ab_cmap, ab_norm, ab_ticks, ab_vmax = build_abundance_cmap(abund_mat)

    ax_abund.imshow(abund_mat, aspect='auto', cmap=ab_cmap, norm=ab_norm,
                    interpolation='nearest')
    ax_abund.set_xlim(-0.5, 2.5); ax_abund.set_ylim(n_mags-0.5, -0.5)
    ax_abund.set_yticks([]); ax_abund.set_xticks([])
    ax_abund.tick_params(bottom=False, labelbottom=False)
    # 顶部分组色条 + 组名在色块上方
    for j, (grp, col) in enumerate([('CK','#1c9cbd'),('A','#e3943d'),('B','#92181e')]):
        ax_abund.add_patch(Rectangle((j-0.5, -2.5), 1, 1.8,
                           facecolor=col, edgecolor='none', clip_on=False, alpha=0.8))
        grp_disp = f'{grp}组' if is_cn else grp
        ax_abund.text(j, -3.2, grp_disp, ha='center', va='center',
                     fontsize=6, fontweight='bold', fontfamily=fm_ if is_cn else fn_,
                     color=col, clip_on=False)
    # 网格
    for i in range(n_mags+1):
        ax_abund.axhline(y=i-0.5, color='white', lw=0.25)
    ax_abund.axvline(x=0.5, color='white', lw=0.3)
    ax_abund.axvline(x=1.5, color='white', lw=0.3)

    # ---- 门颜色条 ----
    for i, (_,row) in enumerate(df.iterrows()):
        ax_phy.add_patch(Rectangle((0, i-0.5), 1, 1,
                         facecolor=PHYLUM_COLORS.get(row['Phylum'],'#BDC3C7'), edgecolor='none'))
    ax_phy.set_xlim(0,1); ax_phy.set_ylim(n_mags-0.5,-0.5); ax_phy.axis('off')

    # ---- Y轴标签 ----
    ax_ylab.set_xlim(0,1); ax_ylab.set_ylim(n_mags-0.5,-0.5); ax_ylab.axis('off')
    for i, (_,row) in enumerate(df.iterrows()):
        color = PHYLUM_COLORS.get(row['Phylum'],'#666')
        is_sp = row['is_top10'] or row['is_keystone']
        fw = 'bold' if is_sp else 'normal'
        fs = 4.3 if is_sp else 3.8
        ax_ylab.text(0.98, i, row['label'], ha='right', va='center',
                    fontsize=fs, fontfamily=fn_, fontstyle='italic', fontweight=fw, color=color)
        x_cur = 0.01
        if row['is_top10']:
            rank = int(row['top_rank'])
            ax_ylab.plot(x_cur, i, marker='o', markersize=4.5, color=TOP10_BC,
                        markeredgewidth=0, zorder=5, clip_on=False)
            ax_ylab.text(x_cur+0.02, i, str(rank), ha='left', va='center',
                        fontsize=3.5, fontfamily=fn_, fontweight='bold', color=TOP10_BC)
            x_cur += 0.055
        if row['is_keystone']:
            ax_ylab.plot(x_cur, i, marker='D', markersize=3, color=KS_BC,
                        markeredgewidth=0, zorder=5, clip_on=False)

    # ---- 底部功能色条 + X标签 ----
    ax_func.set_xlim(-0.5, n_kos-0.5); ax_func.set_ylim(0, 3.5); ax_func.axis('off')
    for j, func in enumerate(ko_funcs):
        c = FUNC_BG.get(func,'#DDDDDD')
        ax_func.add_patch(Rectangle((j-0.5, 2.0), 1, 1.2,
                          facecolor=c, edgecolor='white', linewidth=0.3))
    for elem, start in elem_starts.items():
        if start > 0:
            ax_func.plot([start-0.5, start-0.5], [2.0, 3.2], color='#333', lw=1.5, clip_on=True)
    for j, name in enumerate(ko_display):
        ax_func.text(j, 1.8, name, ha='right', va='top', rotation=45,
                    fontsize=5.5, fontfamily=fn_, fontstyle='italic')

    # ---- 图例（紧凑布局，全部完整显示） ----
    ax_leg.axis('off')
    y_ = 0.99
    S = 0.013  # 行距

    # 门（全部显示）
    ax_leg.text(0.02, y_, '门水平' if is_cn else 'Phylum', fontsize=5.5, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_ -= 0.007
    for phy in phy_order:
        y_ -= S
        ax_leg.add_patch(Rectangle((0.02, y_-0.003), 0.035, 0.008, facecolor=PHYLUM_COLORS.get(phy,'#BDC3C7'),
                         edgecolor='none', transform=ax_leg.transAxes, clip_on=False))
        ax_leg.text(0.06, y_, phy, fontsize=2.8, fontfamily=fn_, transform=ax_leg.transAxes, va='center')

    # 基因拷贝数
    y_ -= 0.018
    ax_leg.text(0.02, y_, '基因拷贝数' if is_cn else 'Gene copy number', fontsize=5.5, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_ -= 0.008
    for elem, ename_en in [('As','Arsenic'),('N','Nitrogen'),('S','Sulfur'),('Fe','Iron')]:
        y_ -= 0.015
        ns = 18; bw = 0.18/ns
        cmap_e = ELEM_CMAPS[elem]
        for s in range(ns):
            frac = s/ns; c = cmap_e(frac*0.85+0.15 if frac>0 else 0)
            ax_leg.add_patch(Rectangle((0.02+s*bw, y_-0.003), bw, 0.008, facecolor=c, edgecolor='none',
                             transform=ax_leg.transAxes, clip_on=False))
        ename = {'As':'砷','N':'氮','S':'硫','Fe':'铁'}[elem] if is_cn else ename_en
        ax_leg.text(0.22, y_, ename, fontsize=3, fontfamily=fm_ if is_cn else fn_,
                   transform=ax_leg.transAxes, va='center')

    # 行标注
    y_ -= 0.018
    ax_leg.text(0.02, y_, '行标注' if is_cn else 'Row annotation', fontsize=5.5, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_ -= 0.013
    ax_leg.add_patch(Rectangle((0.02, y_-0.004), 0.04, 0.010, facecolor='none', edgecolor=TOP10_BC, linewidth=1.2,
                     transform=ax_leg.transAxes, clip_on=False))
    ax_leg.plot(0.07, y_, marker='o', markersize=3.5, color=TOP10_BC, markeredgewidth=0,
               transform=ax_leg.transAxes, clip_on=False)
    ax_leg.text(0.09, y_, 'n 丰度排名' if is_cn else 'n Rank', fontsize=3, fontfamily=fm_ if is_cn else fn_,
               transform=ax_leg.transAxes, va='center')
    y_ -= 0.015
    ax_leg.add_patch(Rectangle((0.02, y_-0.004), 0.04, 0.010, facecolor='none', edgecolor=KS_BC, linewidth=1.2,
                     transform=ax_leg.transAxes, clip_on=False))
    ax_leg.plot(0.07, y_, marker='D', markersize=2.5, color=KS_BC, markeredgewidth=0,
               transform=ax_leg.transAxes, clip_on=False)
    ax_leg.text(0.09, y_, '关键物种' if is_cn else 'Keystone', fontsize=3, fontfamily=fm_ if is_cn else fn_,
               transform=ax_leg.transAxes, va='center')
    y_ -= 0.015
    ax_leg.add_patch(Rectangle((0.02, y_-0.004), 0.04, 0.010, facecolor='none', edgecolor=BOTH_BC, linewidth=1.2,
                     transform=ax_leg.transAxes, clip_on=False))
    ax_leg.text(0.09, y_, '两者兼具' if is_cn else 'Both', fontsize=3, fontfamily=fm_ if is_cn else fn_,
               transform=ax_leg.transAxes, va='center')

    # 组平均丰度色标（竖向，高→低）
    y_ -= 0.02
    ax_leg.text(0.02, y_, '组平均丰度' if is_cn else 'Group abundance', fontsize=5.5, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_ -= 0.012
    bx, bwl = 0.02, 0.03
    tick_rev = list(reversed(ab_ticks))
    n_blk = len(ab_ticks) - 1
    bh = 0.011
    for b in range(n_blk):
        hi_v = tick_rev[b]; lo_v = tick_rev[b+1]
        for ss in range(5):
            frac = ss / 5; val = hi_v - (hi_v - lo_v) * frac; c = ab_cmap(ab_norm(val))
            yb = y_ - b * bh - frac * bh
            ax_leg.add_patch(Rectangle((bx, yb - bh/5), bwl, bh/5, facecolor=c, edgecolor='none',
                             transform=ax_leg.transAxes, clip_on=False))
    for ti, tv in enumerate(tick_rev):
        yt = y_ - ti * bh
        ax_leg.plot([bx+bwl, bx+bwl+0.006], [yt, yt], color='black', lw=0.3,
                   transform=ax_leg.transAxes, clip_on=False)
        ax_leg.text(bx+bwl+0.009, yt, f'{tv:.2f}', fontsize=2.5, fontfamily=fn_,
                   transform=ax_leg.transAxes, va='center')
    ax_leg.add_patch(Rectangle((bx, y_ - n_blk*bh), bwl, n_blk*bh, facecolor='none', edgecolor='black',
                     linewidth=0.3, transform=ax_leg.transAxes, clip_on=False))
    y_ -= n_blk * bh + 0.008

    # 实验分组
    y_ -= 0.012
    ax_leg.text(0.02, y_, '实验分组' if is_cn else 'Treatment group', fontsize=5.5, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    for grp, gcn, gen in [('CK','CK组','Group CK'),('A','A组','Group A'),('B','B组','Group B')]:
        y_ -= S
        ax_leg.add_patch(Rectangle((0.02, y_-0.003), 0.035, 0.008,
                         facecolor={'CK':'#1c9cbd','A':'#e3943d','B':'#92181e'}[grp],
                         edgecolor='none', transform=ax_leg.transAxes, clip_on=False))
        ax_leg.text(0.06, y_, gcn if is_cn else gen, fontsize=3, fontfamily=fm_ if is_cn else fn_,
                   transform=ax_leg.transAxes, va='center')

    # 功能分类
    y_ -= 0.018
    ax_leg.text(0.02, y_, '功能分类' if is_cn else 'Functional category', fontsize=5.5, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_ -= 0.007
    func_groups = {
        'As':['Arsenate reduction','Arsenite oxidation','Resp. arsenate red.',
              'As transport/detox','As methylation','As regulation'],
        'N':['Nitrate reduction','Nitrite reduction','NO reduction',
             'N2O reduction','Ammonia oxidation','N fixation'],
        'S':['Assim. sulfate red.','Dissim. sulfate red.','Sulfide oxidation','Thiosulfate metab.'],
        'Fe':['Fe transport','Fe uptake regulation'],
    }
    for funcs in func_groups.values():
        for func in funcs:
            y_ -= S
            if y_ < 0.003: break
            c = FUNC_BG.get(func, '#DDD')
            ax_leg.add_patch(Rectangle((0.02, y_-0.003), 0.035, 0.008, facecolor=c, edgecolor='#ccc', linewidth=0.2,
                             transform=ax_leg.transAxes, clip_on=False))
            disp = FUNC_CN.get(func, func) if is_cn else func.replace('N2O', 'N2O')
            f_font = fm_ if (is_cn and func in FUNC_CN) else fn_
            ax_leg.text(0.06, y_, disp, fontsize=2.8, fontfamily=f_font, transform=ax_leg.transAxes, va='center')

    # 标题
    if is_cn:
        title = 'MAG元素循环基因谱'
    else:
        title = 'MAG Element Cycle Gene Profiles'
    fig.suptitle(title, fontsize=14, fontweight='bold',
                fontfamily=fm_ if is_cn else fn_, y=0.998)

    # 保存 PDF + PNG
    dpi = 300 if is_cn else 600
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight', facecolor='white')
    # 同步生成PNG
    png_path = outpath.rsplit('.', 1)[0] + '.png'
    plt.savefig(png_path, dpi=min(dpi, 300), bbox_inches='tight', facecolor='white')
    print(f"  Saved: {outpath} + .png ({dpi}dpi)")
    plt.close()

# ================================================================
# 数据导出
# ================================================================
def export_stats(df, active_kos, ko_element, mat, top10_rank, ks_mags, stat_dir):
    """导出供论文描述的统计数据"""
    os.makedirs(stat_dir, exist_ok=True)

    # 1. 全局统计
    with open(os.path.join(stat_dir, 'Fig3-3_gene_profile_summary.txt'), 'w') as f:
        f.write("=== MAG Element Cycle Gene Profile Summary ===\n\n")
        f.write(f"Total MAGs: {len(df)}\n")
        f.write(f"Active KOs (present in >=1 MAG): {len(active_kos)}\n")
        removed = 57 - len(active_kos)
        f.write(f"Removed KOs (absent in all MAGs): {removed}\n\n")

        # 各元素基因数
        for elem in ['As', 'N', 'S', 'Fe']:
            elem_kos = [ko for ko in active_kos if ko_element[ko] == elem]
            f.write(f"{elem}: {len(elem_kos)} active KOs\n")

        # 每个MAG的基因数统计
        gene_counts = df[active_kos].sum(axis=1)
        f.write(f"\nGenes per MAG: mean={gene_counts.mean():.1f}, "
                f"median={gene_counts.median():.0f}, "
                f"min={gene_counts.min()}, max={gene_counts.max()}\n")

        # 各门统计
        f.write(f"\n=== Phylum distribution ===\n")
        for phy, cnt in df['Phylum'].value_counts().items():
            avg_genes = df[df['Phylum']==phy][active_kos].sum(axis=1).mean()
            f.write(f"  {phy}: {cnt} MAGs, avg {avg_genes:.1f} genes\n")

    # 2. 完整MAG×KO presence矩阵（含分类+组丰度+富集方向）
    out_df = df[['MAG', 'Phylum', 'label'] + active_kos].copy()
    out_df['total_genes'] = df[active_kos].sum(axis=1)
    out_df['is_top10'] = df['MAG'].isin(top10_rank.keys())
    out_df['top_rank'] = df['MAG'].map(top10_rank)
    out_df['is_keystone'] = df['MAG'].isin(ks_mags)
    for col in ['CK_mean', 'A_mean', 'B_mean']:
        if col in df.columns:
            out_df[col] = df[col]
    def get_enriched(row):
        ck = row.get('CK_mean', 0); a = row.get('A_mean', 0); b = row.get('B_mean', 0)
        vals = {'CK': ck, 'A': a, 'B': b}
        enriched = []
        for grp, val in vals.items():
            others = [v for k, v in vals.items() if k != grp]
            if val > 1.5 * np.mean(others) and val > 0.05:
                enriched.append(grp)
        return ','.join(enriched) if enriched else '-'
    out_df['enriched_in'] = out_df.apply(get_enriched, axis=1)
    out_df.to_csv(os.path.join(stat_dir, 'Fig3-3_MAG_KO_matrix.tsv'),
                  sep='\t', index=False)

    # 3. Top10专项统计
    top10_df = df[df['MAG'].isin(top10_rank.keys())].copy()
    top10_df['rank'] = top10_df['MAG'].map(top10_rank)
    top10_df = top10_df.sort_values('rank')
    with open(os.path.join(stat_dir, 'Fig3-3_Top10_summary.txt'), 'w') as f:
        f.write("=== Top 10 Abundant MAGs - Functional Summary ===\n\n")
        f.write(f"{'Rank':<5} {'MAG':<15} {'Phylum':<20} {'Species':<30} "
                f"{'As':>3} {'N':>3} {'S':>3} {'Fe':>3} {'Total':>5}\n")
        f.write("-" * 90 + "\n")
        for _, row in top10_df.iterrows():
            as_n = sum(1 for ko in active_kos if ko_element[ko]=='As' and row[ko]>0)
            n_n  = sum(1 for ko in active_kos if ko_element[ko]=='N'  and row[ko]>0)
            s_n  = sum(1 for ko in active_kos if ko_element[ko]=='S'  and row[ko]>0)
            fe_n = sum(1 for ko in active_kos if ko_element[ko]=='Fe' and row[ko]>0)
            total = as_n + n_n + s_n + fe_n
            f.write(f"{int(row['rank']):<5} {row['MAG']:<15} {row['Phylum']:<20} "
                    f"{row['label']:<30} {as_n:>3} {n_n:>3} {s_n:>3} {fe_n:>3} {total:>5}\n")

    # 4. Keystone专项统计
    ks_df = df[df['MAG'].isin(ks_mags)].copy()
    with open(os.path.join(stat_dir, 'Fig3-3_Keystone_summary.txt'), 'w') as f:
        f.write("=== Keystone Species - Functional Summary ===\n\n")
        f.write(f"{'MAG':<15} {'Phylum':<20} {'Species':<30} "
                f"{'As':>3} {'N':>3} {'S':>3} {'Fe':>3} {'Total':>5}\n")
        f.write("-" * 85 + "\n")
        for _, row in ks_df.iterrows():
            as_n = sum(1 for ko in active_kos if ko_element[ko]=='As' and row[ko]>0)
            n_n  = sum(1 for ko in active_kos if ko_element[ko]=='N'  and row[ko]>0)
            s_n  = sum(1 for ko in active_kos if ko_element[ko]=='S'  and row[ko]>0)
            fe_n = sum(1 for ko in active_kos if ko_element[ko]=='Fe' and row[ko]>0)
            total = as_n + n_n + s_n + fe_n
            f.write(f"{row['MAG']:<15} {row['Phylum']:<20} {row['label']:<30} "
                    f"{as_n:>3} {n_n:>3} {s_n:>3} {fe_n:>3} {total:>5}\n")

    # 5. 各KO的MAG覆盖率
    with open(os.path.join(stat_dir, 'Fig3-3_KO_coverage.txt'), 'w') as f:
        f.write("=== KO Gene Coverage across MAGs ===\n\n")
        f.write(f"{'KO_label':<25} {'Gene':<12} {'Element':<5} {'Function':<25} "
                f"{'MAGs_with':>9} {'Pct':>6}\n")
        f.write("-" * 90 + "\n")
        for ko in active_kos:
            k_id = ko.split('(')[1].rstrip(')')
            gn, func = KO_FUNC.get(k_id, (ko.split('(')[0], 'Unknown'))
            elem = ko_element[ko]
            n_with = (df[ko] > 0).sum()
            pct = n_with / len(df) * 100
            f.write(f"{ko:<25} {gn:<12} {elem:<5} {func:<25} {n_with:>9} {pct:>5.1f}%\n")

    print(f"  Stats exported to: {stat_dir}")


def make_subset_figure(df, active_kos, ko_labels, ko_element, copy_pivot,
                       subset_mags, fig_label, style, outpath):
    """为Top10或Keystone子集生成独立热图"""
    fm_, fn_ = get_fonts(style)
    is_cn = (style == 'cn')

    # 筛选子集
    sub = df[df['MAG'].isin(subset_mags)].copy()
    if 'top_rank' in sub.columns and sub['top_rank'].notna().any():
        sub = sub.sort_values('top_rank')
    else:
        sub = sub.sort_values(['phy_rank', 'gene_count'], ascending=[True, False])
    sub = sub.reset_index(drop=True)

    n_mags = len(sub)
    n_kos = len(active_kos)

    # 拷贝数矩阵
    mat_copy = np.zeros((n_mags, n_kos))
    for j, ko in enumerate(active_kos):
        for i, (_, row) in enumerate(sub.iterrows()):
            mag = row['MAG']
            if mag in copy_pivot.index and ko in copy_pivot.columns:
                mat_copy[i, j] = copy_pivot.loc[mag, ko]
    mat_log = np.log1p(mat_copy)
    max_log = mat_log.max() if mat_log.max() > 0 else 1

    # X轴标签
    ko_display = []
    ko_funcs = []
    for ko in active_kos:
        k_id = ko.split('(')[1].rstrip(')')
        if k_id in KO_FUNC:
            gn, func = KO_FUNC[k_id]
            ko_display.append(gn); ko_funcs.append(func)
        else:
            ko_display.append(ko.split('(')[0]); ko_funcs.append('Unknown')

    # 尺寸（子集行少，加高保证图例完整）
    row_h = 0.35 if n_mags <= 15 else 0.25
    n_leg_items = len(sub['Phylum'].unique()) + 4 + 3 + 8 + 3 + 18
    min_h = n_leg_items * 0.28 + 6
    fig_h = max(n_mags * row_h + 5, min_h)
    fig_w = n_kos * 0.30 + 10

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(2, 5,
        width_ratios=[0.4, 4.2, n_kos * 0.62, 2.0, 2.8],
        height_ratios=[n_mags, 3.5], wspace=0.02, hspace=0.0)
    ax_phy  = fig.add_subplot(gs[0, 0])
    ax_ylab = fig.add_subplot(gs[0, 1])
    ax_heat = fig.add_subplot(gs[0, 2])
    ax_abund = fig.add_subplot(gs[0, 3])
    ax_leg  = fig.add_subplot(gs[0, 4])
    ax_func = fig.add_subplot(gs[1, 2])

    # ---- 热图 ----
    for j, ko in enumerate(active_kos):
        elem = ko_element[ko]
        cmap_e = ELEM_CMAPS[elem]
        for i in range(n_mags):
            val = mat_log[i, j]
            if val > 0:
                color = cmap_e(min(val / max_log * 0.85 + 0.15, 1.0))
                ax_heat.add_patch(Rectangle((j-0.45, i-0.45), 0.9, 0.9,
                                  facecolor=color, edgecolor='none', zorder=1))

    for i in range(n_mags + 1):
        ax_heat.axhline(y=i-0.5, color='white', lw=0.4, zorder=2)
    for j in range(n_kos + 1):
        ax_heat.axvline(x=j-0.5, color='white', lw=0.4, zorder=2)

    ax_heat.set_yticks([]); ax_heat.set_xticks([])
    ax_heat.set_xlim(-0.5, n_kos - 0.5)
    ax_heat.set_ylim(n_mags - 0.5, -0.5)

    # 元素分隔线
    elem_starts = {}
    for j, ko in enumerate(active_kos):
        e = ko_element[ko]
        if e not in elem_starts: elem_starts[e] = j
    for elem, start in elem_starts.items():
        if start > 0:
            ax_heat.plot([start-0.5, start-0.5], [-0.5, n_mags-0.5],
                        color='#333', lw=1.5, zorder=3, clip_on=True)

    # 顶部元素色条 — 根据行数动态调整位置
    bar_y = -1.2  # 色条底部y
    bar_h = 0.8   # 色条高度
    label_y = bar_y - 0.8  # 元素名称y
    for j, ko in enumerate(active_kos):
        elem = ko_element[ko]
        ax_heat.add_patch(Rectangle((j-0.5, bar_y), 1, bar_h,
                          facecolor=ELEM_FILL[elem], edgecolor='none', clip_on=False, alpha=0.85))
    elem_positions = {'As':[],'N':[],'S':[],'Fe':[]}
    for j, ko in enumerate(active_kos):
        elem_positions[ko_element[ko]].append(j)
    enames = {'As':'砷','N':'氮','S':'硫','Fe':'铁'} if is_cn else \
             {'As':'Arsenic','N':'Nitrogen','S':'Sulfur','Fe':'Iron'}
    for elem, positions in elem_positions.items():
        if positions:
            mid = (positions[0]+positions[-1])/2
            ax_heat.text(mid, label_y, enames[elem], ha='center', va='center',
                        fontsize=9, fontweight='bold',
                        fontfamily=fm_ if is_cn else fn_, color=ELEM_FILL[elem])

    # ---- CK/A/B 丰度侧栏（子图） ----
    abund_cols2 = ['CK_mean', 'A_mean', 'B_mean']
    abund_mat2 = sub[abund_cols2].values
    ab_cmap2, ab_norm2, ab_ticks2, ab_vmax2 = build_abundance_cmap(abund_mat2)

    ax_abund.imshow(abund_mat2, aspect='auto', cmap=ab_cmap2, norm=ab_norm2,
                    interpolation='nearest')
    ax_abund.set_xlim(-0.5, 2.5); ax_abund.set_ylim(n_mags-0.5, -0.5)
    ax_abund.set_yticks([]); ax_abund.set_xticks([])
    ax_abund.tick_params(bottom=False, labelbottom=False)
    for j2, (grp2, col2) in enumerate([('CK','#1c9cbd'),('A','#e3943d'),('B','#92181e')]):
        ax_abund.add_patch(Rectangle((j2-0.5, -1.5), 1, 1.0,
                           facecolor=col2, edgecolor='none', clip_on=False, alpha=0.8))
        grp_d = f'{grp2}组' if is_cn else grp2
        ax_abund.text(j2, -2.0, grp_d, ha='center', va='center',
                     fontsize=6, fontweight='bold', fontfamily=fm_ if is_cn else fn_,
                     color=col2, clip_on=False)
    for i2 in range(n_mags+1):
        ax_abund.axhline(y=i2-0.5, color='white', lw=0.3)
    ax_abund.axvline(x=0.5, color='white', lw=0.3)
    ax_abund.axvline(x=1.5, color='white', lw=0.3)

    # ---- 门颜色条 ----
    for i, (_, row) in enumerate(sub.iterrows()):
        ax_phy.add_patch(Rectangle((0, i-0.5), 1, 1,
                         facecolor=PHYLUM_COLORS.get(row['Phylum'], '#BDC3C7'), edgecolor='none'))
    ax_phy.set_xlim(0, 1); ax_phy.set_ylim(n_mags-0.5, -0.5); ax_phy.axis('off')

    # ---- Y轴标签（子集行少，字体可大些） ----
    ax_ylab.set_xlim(0, 1); ax_ylab.set_ylim(n_mags-0.5, -0.5); ax_ylab.axis('off')
    for i, (_, row) in enumerate(sub.iterrows()):
        color = PHYLUM_COLORS.get(row['Phylum'], '#666')
        ax_ylab.text(0.98, i, row['label'], ha='right', va='center',
                    fontsize=7, fontfamily=fn_, fontstyle='italic',
                    fontweight='bold', color=color)

    # ---- 底部功能色条 ----
    ax_func.set_xlim(-0.5, n_kos-0.5); ax_func.set_ylim(0, 3.5); ax_func.axis('off')
    for j, func in enumerate(ko_funcs):
        c = FUNC_BG.get(func, '#DDD')
        ax_func.add_patch(Rectangle((j-0.5, 2.0), 1, 1.2,
                          facecolor=c, edgecolor='white', linewidth=0.3))
    for elem, start in elem_starts.items():
        if start > 0:
            ax_func.plot([start-0.5, start-0.5], [2.0, 3.2], color='#333', lw=1.5, clip_on=True)
    for j, name in enumerate(ko_display):
        ax_func.text(j, 1.8, name, ha='right', va='top', rotation=45,
                    fontsize=6.5, fontfamily=fn_, fontstyle='italic')

    # ---- 图例（子图，含所有图例项） ----
    ax_leg.axis('off')
    y_ = 0.98
    S = 0.028

    # 门
    sub_phyla = sub['Phylum'].unique().tolist()
    ax_leg.text(0.02, y_, '门水平' if is_cn else 'Phylum', fontsize=7, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_ -= 0.012
    for phy in sub_phyla:
        y_ -= S
        ax_leg.add_patch(Rectangle((0.02, y_-0.008), 0.05, 0.018,
                         facecolor=PHYLUM_COLORS.get(phy, '#BDC3C7'),
                         edgecolor='none', transform=ax_leg.transAxes, clip_on=False))
        ax_leg.text(0.08, y_, phy, fontsize=5, fontfamily=fn_,
                   transform=ax_leg.transAxes, va='center')

    # 拷贝数渐变
    y_ -= 0.035
    ax_leg.text(0.02, y_, '基因拷贝数' if is_cn else 'Gene copy number', fontsize=7, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_ -= 0.015
    for elem, ename_en in [('As','Arsenic'),('N','Nitrogen'),('S','Sulfur'),('Fe','Iron')]:
        y_ -= 0.032
        ns = 20; bw = 0.22/ns
        cmap_e = ELEM_CMAPS[elem]
        for s in range(ns):
            frac = s/ns; c = cmap_e(frac*0.85+0.15 if frac>0 else 0)
            ax_leg.add_patch(Rectangle((0.02+s*bw, y_-0.006), bw, 0.015,
                             facecolor=c, edgecolor='none', transform=ax_leg.transAxes, clip_on=False))
        ename = {'As':'砷','N':'氮','S':'硫','Fe':'铁'}[elem] if is_cn else ename_en
        ax_leg.text(0.26, y_, ename, fontsize=5, fontfamily=fm_ if is_cn else fn_,
                   transform=ax_leg.transAxes, va='center')

    # 组平均丰度色标（竖向，高→低）
    y_ -= 0.04
    ax_leg.text(0.02, y_, '组平均丰度' if is_cn else 'Group abundance', fontsize=7, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_ -= 0.012
    bx, bwl = 0.02, 0.04
    tick_rev2 = list(reversed(ab_ticks2))
    n_blk2 = len(ab_ticks2) - 1
    bh2 = 0.022
    for b in range(n_blk2):
        hi_v = tick_rev2[b]; lo_v = tick_rev2[b+1]
        for ss in range(6):
            frac = ss / 6; val = hi_v - (hi_v - lo_v) * frac; c = ab_cmap2(ab_norm2(val))
            yb = y_ - b * bh2 - frac * bh2
            ax_leg.add_patch(Rectangle((bx, yb - bh2/6), bwl, bh2/6, facecolor=c, edgecolor='none',
                             transform=ax_leg.transAxes, clip_on=False))
    for ti, tv in enumerate(tick_rev2):
        yt = y_ - ti * bh2
        ax_leg.plot([bx+bwl, bx+bwl+0.01], [yt, yt], color='black', lw=0.4,
                   transform=ax_leg.transAxes, clip_on=False)
        ax_leg.text(bx+bwl+0.015, yt, f'{tv:.2f}', fontsize=4.5, fontfamily=fn_,
                   transform=ax_leg.transAxes, va='center')
    ax_leg.add_patch(Rectangle((bx, y_ - n_blk2*bh2), bwl, n_blk2*bh2, facecolor='none',
                     edgecolor='black', linewidth=0.4, transform=ax_leg.transAxes, clip_on=False))
    y_ -= n_blk2 * bh2 + 0.012

    # 实验分组
    y_ -= 0.015
    ax_leg.text(0.02, y_, '实验分组' if is_cn else 'Treatment group', fontsize=7, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    for grp, gcn, gen in [('CK','CK组','Group CK'),('A','A组','Group A'),('B','B组','Group B')]:
        y_ -= S
        ax_leg.add_patch(Rectangle((0.02, y_-0.008), 0.05, 0.018,
                         facecolor={'CK':'#1c9cbd','A':'#e3943d','B':'#92181e'}[grp],
                         edgecolor='none', transform=ax_leg.transAxes, clip_on=False))
        ax_leg.text(0.08, y_, gcn if is_cn else gen, fontsize=5, fontfamily=fm_ if is_cn else fn_,
                   transform=ax_leg.transAxes, va='center')

    # 功能分类
    y_ -= 0.035
    ax_leg.text(0.02, y_, '功能分类' if is_cn else 'Functional category', fontsize=7, fontweight='bold',
               fontfamily=fm_ if is_cn else fn_, transform=ax_leg.transAxes, va='top')
    y_ -= 0.012
    func_groups = {
        'As':['Arsenate reduction','Arsenite oxidation','Resp. arsenate red.',
              'As transport/detox','As methylation','As regulation'],
        'N':['Nitrate reduction','Nitrite reduction','NO reduction',
             'N2O reduction','Ammonia oxidation','N fixation'],
        'S':['Assim. sulfate red.','Dissim. sulfate red.','Sulfide oxidation','Thiosulfate metab.'],
        'Fe':['Fe transport','Fe uptake regulation'],
    }
    for funcs_list in func_groups.values():
        for func in funcs_list:
            y_ -= S
            if y_ < 0.005: break
            c = FUNC_BG.get(func, '#DDD')
            ax_leg.add_patch(Rectangle((0.02, y_-0.008), 0.05, 0.018,
                             facecolor=c, edgecolor='#ccc', linewidth=0.3,
                             transform=ax_leg.transAxes, clip_on=False))
            disp = FUNC_CN.get(func, func) if is_cn else func.replace('N2O', 'N2O')
            f_font = fm_ if (is_cn and func in FUNC_CN) else fn_
            ax_leg.text(0.08, y_, disp, fontsize=5, fontfamily=f_font,
                       transform=ax_leg.transAxes, va='center')

    # 标题
    if is_cn:
        title = fig_label + '元素循环基因谱'
    else:
        title = fig_label + ' Element Cycle Gene Profiles'
    fig.suptitle(title, fontsize=14, fontweight='bold',
                fontfamily=fm_ if is_cn else fn_, y=0.998)

    dpi = 300 if is_cn else 600
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight', facecolor='white')
    png_path = outpath.rsplit('.', 1)[0] + '.png'
    plt.savefig(png_path, dpi=min(dpi, 300), bbox_inches='tight', facecolor='white')
    print(f"  Saved: {outpath} + .png")
    plt.close()


# ================================================================
# Main
# ================================================================
def main():
    parser = argparse.ArgumentParser(description='Fig3-3/3-6 MAG Gene Profile Heatmaps')
    parser.add_argument('--itol',  required=True, help='iTOL heatmap file (07_element_KO_heatmap.txt)')
    parser.add_argument('--ko',    required=True, help='KEGG target KO file (kegg_target_only.tsv)')
    parser.add_argument('--bac',   required=True, help='GTDB bac120 summary')
    parser.add_argument('--ar',    required=True, help='GTDB ar53 summary')
    parser.add_argument('--abund', required=True, help='CoverM abundance.tsv')
    parser.add_argument('--ks',    required=True, help='Keystone species list')
    parser.add_argument('--outdir',default='figures', help='Output base directory')
    parser.add_argument('--statdir',default='data/processed', help='Stats output directory')
    args = parser.parse_args()

    print(f"Fonts: sans={F_SANS}, serif={F_SERIF}, cn={F_CN}")

    df, active_kos, ko_labels, ko_element, mat, mat_log, phy_order = load_data(args)
    print(f"Data: {len(df)} MAGs × {len(active_kos)} KOs")

    # 获取top10和keystone集合（供stats导出用）
    ab = pd.read_csv(args.abund, sep='\t')
    ab = ab[ab['Genome'] != 'unmapped']
    scols = [c for c in ab.columns if c != 'Genome']
    ab['mean_all'] = ab[scols].mean(axis=1)
    ab = ab.sort_values('mean_all', ascending=False).reset_index(drop=True)
    top10_rank = {ab.iloc[i]['Genome']: i+1 for i in range(10)}
    ks = pd.read_csv(args.ks, sep='\t')
    ks_mags = set(ks['MAG'].tolist())

    # 拷贝数pivot（子图需要）
    raw_ko = pd.read_csv(args.ko, sep='\t')
    rows_exp = []
    for _, r in raw_ko.iterrows():
        for ko in str(r['KEGG_ko']).split(','):
            rows_exp.append((r['MAG'], ko.strip()))
    df_exp = pd.DataFrame(rows_exp, columns=['MAG','KO'])
    tmap = {}
    for ko in ko_labels:
        k_id = ko.split('(')[1].rstrip(')')
        tmap[f'ko:{k_id}'] = ko
    df_t = df_exp[df_exp['KO'].isin(tmap.keys())]
    cn = df_t.groupby(['MAG','KO']).size().reset_index(name='copies')
    cn['ko_label'] = cn['KO'].map(tmap)
    copy_pivot = cn.pivot_table(index='MAG', columns='ko_label', values='copies', fill_value=0)

    en_dir = os.path.join(args.outdir, 'paper_en')
    cn_dir = os.path.join(args.outdir, 'thesis_cn')
    os.makedirs(en_dir, exist_ok=True)
    os.makedirs(cn_dir, exist_ok=True)

    # ==========================================
    # Fig3-3: 全168 MAG基因谱
    # ==========================================
    print("\n=== Fig3-3: Full 168 MAGs ===")
    draw_figure(df, active_kos, ko_labels, ko_element, mat, mat_log, phy_order,
                'en', os.path.join(en_dir, 'Fig3-3_MAG_gene_profile_EN.pdf'))
    draw_figure(df, active_kos, ko_labels, ko_element, mat, mat_log, phy_order,
                'cn', os.path.join(cn_dir, 'Fig3-3_MAG_gene_profile_CN.pdf'))
    # EN TIFF
    draw_figure(df, active_kos, ko_labels, ko_element, mat, mat_log, phy_order,
                'en', os.path.join(en_dir, 'Fig3-3_MAG_gene_profile_EN.tiff'))

    # ==========================================
    # Fig3-6a: Top10 MAGs 基因谱
    # ==========================================
    print("\n=== Fig3-6a: Top 10 Abundant MAGs ===")
    make_subset_figure(df, active_kos, ko_labels, ko_element, copy_pivot,
                       set(top10_rank.keys()), 'Top 10 Abundant MAGs',
                       'en', os.path.join(en_dir, 'Fig3-6a_Top10_gene_profile_EN.pdf'))
    make_subset_figure(df, active_kos, ko_labels, ko_element, copy_pivot,
                       set(top10_rank.keys()), '丰度Top10 MAGs',
                       'cn', os.path.join(cn_dir, 'Fig3-6a_Top10_gene_profile_CN.pdf'))

    # ==========================================
    # Fig3-6b: Keystone MAGs 基因谱
    # ==========================================
    print("\n=== Fig3-6b: Keystone Species ===")
    make_subset_figure(df, active_kos, ko_labels, ko_element, copy_pivot,
                       ks_mags, 'Keystone Species',
                       'en', os.path.join(en_dir, 'Fig3-6b_Keystone_gene_profile_EN.pdf'))
    make_subset_figure(df, active_kos, ko_labels, ko_element, copy_pivot,
                       ks_mags, '关键物种',
                       'cn', os.path.join(cn_dir, 'Fig3-6b_Keystone_gene_profile_CN.pdf'))

    # ==========================================
    # 数据导出
    # ==========================================
    print("\n=== Exporting stats ===")
    export_stats(df, active_kos, ko_element, mat, top10_rank, ks_mags, args.statdir)

    print("\nAll done!")

if __name__ == '__main__':
    main()
