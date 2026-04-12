#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_MAG_quality.py — Fig3-1 MAG质量评估散点图
X=Completeness, Y=Contamination, 按门着色, 大小=基因组大小
高亮Keystone物种, 标注高/中质量区域分界线
"""

import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.patches as mpatches
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')

# 字体
_avail = {f.name for f in fm.fontManager.ttflist}
def _pick(cs, fb='DejaVu Sans'):
    for c in cs:
        if c in _avail: return c
    return fb
F_SANS = _pick(['Arial','Helvetica','Liberation Sans'],'DejaVu Sans')
F_SERIF = _pick(['Times New Roman','Liberation Serif'],'DejaVu Serif')
F_CN = _pick(['SimSun','WenQuanYi Micro Hei'], F_SERIF)

def get_fonts(style):
    if style=='thesis': return F_CN, F_SERIF
    return F_SANS, F_SANS

# 门颜色（Top 12 + Other）
PHYLUM_COLORS = {
    'Pseudomonadota':      '#E74C3C',
    'Acidobacteriota':     '#3498DB',
    'Chloroflexota':       '#2ECC71',
    'Bacteroidota_A':      '#F39C12',
    'Patescibacteriota':   '#9B59B6',
    'Desulfobacterota':    '#1ABC9C',
    'Actinomycetota':      '#E67E22',
    'Gemmatimonadota':     '#34495E',
    'Planctomycetota':     '#D35400',
    'Myxococcota':         '#8E44AD',
    'Nitrospirota':        '#27AE60',
    'Zixibacteria':        '#C0392B',
    'Other':               '#BDC3C7',
}

def load_data(checkm2_file, bac_file, ar_file, keystone_file):
    """加载并合并所有数据"""
    # CheckM2
    ck = pd.read_csv(checkm2_file, sep='\t')
    ck = ck.rename(columns={'Name':'MAG'})
    
    # GTDB分类
    bac = pd.read_csv(bac_file, sep='\t')
    ar = pd.read_csv(ar_file, sep='\t')
    gtdb = pd.concat([bac[['user_genome','classification']], 
                       ar[['user_genome','classification']]], ignore_index=True)
    gtdb = gtdb.rename(columns={'user_genome':'MAG'})
    
    def get_phylum(c):
        for p in c.split(';'):
            if p.startswith('p__'): return p.replace('p__','')
        return 'Unknown'
    gtdb['Phylum'] = gtdb['classification'].apply(get_phylum)
    
    # 合并
    df = ck.merge(gtdb[['MAG','Phylum']], on='MAG', how='left')
    df['Phylum'] = df['Phylum'].fillna('Unknown')
    
    # 门归类（少于3个的归为Other）
    phy_counts = df['Phylum'].value_counts()
    minor_phyla = phy_counts[phy_counts < 3].index.tolist()
    df['Phylum_group'] = df['Phylum'].apply(lambda x: 'Other' if x in minor_phyla else x)
    # 不在预定义颜色中的也归Other
    df['Phylum_group'] = df['Phylum_group'].apply(lambda x: x if x in PHYLUM_COLORS else 'Other')
    
    # Keystone标记
    ks = pd.read_csv(keystone_file, sep='\t')
    df['is_keystone'] = df['MAG'].isin(ks['MAG'])
    df['ks_genus'] = df['MAG'].map(dict(zip(ks['MAG'], ks['Genus'])))
    
    # 质量分类
    def quality(row):
        if row['Completeness']>=90 and row['Contamination']<5: return 'High'
        elif row['Completeness']>=50 and row['Contamination']<10: return 'Medium'
        else: return 'Low'
    df['Quality'] = df.apply(quality, axis=1)
    
    return df

def draw_fig3_1(df, output_prefix, fm_, fn_, style):
    is_cn = (style=='thesis')
    
    fig = plt.figure(figsize=(13, 7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.3], wspace=0.05)
    ax = fig.add_subplot(gs[0])
    ax_leg = fig.add_subplot(gs[1])
    ax_leg.axis('off')
    
    # ---- 质量区域背景 ----
    # 高质量: Completeness≥90, Contamination<5
    ax.axhspan(0, 5, xmin=0, xmax=1, alpha=0.06, color='#27AE60', zorder=0)
    ax.axhspan(5, 10, xmin=0, xmax=1, alpha=0.04, color='#F39C12', zorder=0)
    
    # 分界线
    ax.axhline(y=5, color='#27AE60', ls='--', lw=0.8, alpha=0.6, zorder=1)
    ax.axhline(y=10, color='#E74C3C', ls='--', lw=0.8, alpha=0.6, zorder=1)
    ax.axvline(x=90, color='#27AE60', ls='--', lw=0.8, alpha=0.6, zorder=1)
    ax.axvline(x=50, color='#E74C3C', ls='--', lw=0.8, alpha=0.6, zorder=1)
    
    # 区域标注
    ax.text(95, 0.3, 'High Quality' if not is_cn else '高质量', fontsize=7,
           fontfamily=fm_ if is_cn else fn_, color='#27AE60', fontstyle='italic', alpha=0.7)
    ax.text(70, 6, 'Medium Quality' if not is_cn else '中质量', fontsize=7,
           fontfamily=fm_ if is_cn else fn_, color='#F39C12', fontstyle='italic' if not is_cn else 'normal', alpha=0.7)
    
    # ---- 散点 ----
    # 先画非keystone（小点，半透明）
    non_ks = df[~df['is_keystone']]
    for phy, grp in non_ks.groupby('Phylum_group'):
        color = PHYLUM_COLORS.get(phy, '#BDC3C7')
        sizes = grp['Genome_Size'] / 1e6 * 12  # Mb -> 点大小
        ax.scatter(grp['Completeness'], grp['Contamination'],
                  s=sizes, c=color, alpha=0.55, edgecolors='white',
                  linewidths=0.3, zorder=2, label=phy)
    
    # 再画keystone（大点，黑边，高亮）
    ks = df[df['is_keystone']]
    for _, row in ks.iterrows():
        color = PHYLUM_COLORS.get(row['Phylum_group'], '#BDC3C7')
        size = row['Genome_Size'] / 1e6 * 12
        ax.scatter(row['Completeness'], row['Contamination'],
                  s=size*1.8, c=color, alpha=0.9, edgecolors='black',
                  linewidths=1.2, zorder=4, marker='D')
        
        genus = row['ks_genus']
        if pd.notna(genus):
            label = row['MAG'] if genus == 'Unknown' else genus
            dx, dy = 1.5, 0.4
            if row['Contamination'] < 2: dy = 0.5
            if row['Completeness'] > 95: dx = -3.0
            if row['Completeness'] > 88 and row['Contamination'] > 5: dx = 2.0; dy = -0.6
            ax.annotate(label, (row['Completeness'], row['Contamination']),
                       xytext=(row['Completeness']+dx, row['Contamination']+dy),
                       fontsize=5.5, fontfamily=fn_, fontstyle='italic',
                       fontweight='bold', color='#222',
                       arrowprops=dict(arrowstyle='-', color='#888', lw=0.5),
                       zorder=5)
    
    # ---- 轴 ----
    ax.set_xlim(35, 102)
    ax.set_ylim(-0.5, max(df['Contamination'].max()+1, 12))
    
    if is_cn:
        ax.set_xlabel('完整度 (%)', fontsize=11, fontfamily=fm_)
        ax.set_ylabel('污染度 (%)', fontsize=11, fontfamily=fm_)
    else:
        ax.set_xlabel('Completeness (%)', fontsize=11, fontfamily=fn_)
        ax.set_ylabel('Contamination (%)', fontsize=11, fontfamily=fn_)
    
    # 标题
    hq = (df['Quality']=='High').sum()
    mq = (df['Quality']=='Medium').sum()
    lq = (df['Quality']=='Low').sum()
    if is_cn:
        ax.set_title('MAG质量评估', fontsize=13, fontweight='bold', fontfamily=fm_, pad=18)
        sub = f'共{len(df)}个MAGs（高质量: {hq}, 中质量: {mq}, 低质量: {lq}）'
    else:
        ax.set_title('MAG Quality Assessment', fontsize=13, fontweight='bold', fontfamily=fn_, pad=18)
        sub = f'{len(df)} MAGs (High: {hq}, Medium: {mq}, Low: {lq})'
    ax.text(0.5, 1.015, sub, transform=ax.transAxes, fontsize=7.5,
           fontfamily=fm_ if is_cn else fn_, va='bottom', ha='center', color='#666')
    
    ax.tick_params(labelsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # ---- 图例（画在右侧ax_leg中）----
    from matplotlib.font_manager import FontProperties
    leg_font = FontProperties(family=fn_, size=6.5)
    leg_title_font = FontProperties(family=fn_, size=8, weight='bold')
    
    phy_order = df['Phylum_group'].value_counts().index.tolist()
    phy_handles = []
    for phy in phy_order:
        color = PHYLUM_COLORS.get(phy, '#BDC3C7')
        phy_handles.append(mpatches.Patch(facecolor=color, edgecolor='none',
                                          label=f"{phy} ({(df['Phylum_group']==phy).sum()})"))
    
    leg1 = ax_leg.legend(handles=phy_handles, loc='upper left', bbox_to_anchor=(0.0, 1.0),
                    prop=leg_font, title='Phylum',
                    frameon=False, handlelength=1.2, handleheight=0.8)
    leg1.get_title().set_fontproperties(leg_title_font)
    ax_leg.add_artist(leg1)
    
    # 大小图例
    size_vals = [1, 3, 6, 9]  # Mb
    size_handles = [Line2D([0],[0], marker='o', color='none', markerfacecolor='#999',
                          markersize=np.sqrt(s*12*0.8), label=f'{s} Mb',
                          markeredgecolor='white', markeredgewidth=0.3)
                   for s in size_vals]
    
    leg2 = ax_leg.legend(handles=size_handles, loc='lower left', bbox_to_anchor=(0.0, 0.12),
                    prop=FontProperties(family=fn_, size=7), title='Genome size',
                    frameon=False, handlelength=1.5)
    leg2.get_title().set_fontproperties(leg_title_font)
    ax_leg.add_artist(leg2)
    
    # Keystone标记图例
    ks_handle = Line2D([0],[0], marker='D', color='none', markerfacecolor='#999',
                       markersize=8, markeredgecolor='black', markeredgewidth=1.0,
                       label='Keystone species')
    leg3 = ax_leg.legend(handles=[ks_handle], loc='lower left', bbox_to_anchor=(0.0, 0.03),
                    prop=FontProperties(family=fn_, size=7), frameon=False)
    ax_leg.add_artist(leg3)
    
    for ext in ['pdf','png']:
        fp = f"{output_prefix}.{ext}"
        fig.savefig(fp, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  -> {fp}")
    plt.close(fig)

def save_stats(df, stat_dir, fn_):
    """保存MAG质量统计表"""
    os.makedirs(stat_dir, exist_ok=True)
    
    # 汇总统计
    summary = []
    for q in ['High','Medium','Low']:
        sub = df[df['Quality']==q]
        summary.append({
            'Quality': q,
            'Count': len(sub),
            'Completeness_mean': f"{sub['Completeness'].mean():.1f}",
            'Completeness_range': f"{sub['Completeness'].min():.1f}-{sub['Completeness'].max():.1f}",
            'Contamination_mean': f"{sub['Contamination'].mean():.2f}",
            'Contamination_range': f"{sub['Contamination'].min():.2f}-{sub['Contamination'].max():.2f}",
            'Genome_Size_mean_Mb': f"{sub['Genome_Size'].mean()/1e6:.1f}",
        })
    pd.DataFrame(summary).to_csv(os.path.join(stat_dir, 'MAG_quality_summary.txt'),
                                  sep='\t', index=False)
    
    # 全部MAG详细表
    cols = ['MAG','Completeness','Contamination','Genome_Size','GC_Content',
            'Contig_N50','Total_Contigs','Phylum','Quality','is_keystone']
    df[cols].sort_values(['Quality','Completeness'], ascending=[True, False]).to_csv(
        os.path.join(stat_dir, 'MAG_quality_all.txt'), sep='\t', index=False)
    
    print(f"  -> MAG_quality_summary.txt, MAG_quality_all.txt")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--checkm2', required=True)
    parser.add_argument('--bac120', required=True)
    parser.add_argument('--ar53', required=True)
    parser.add_argument('--keystone', required=True)
    parser.add_argument('-o','--output', default='.')
    parser.add_argument('--stat_dir', default='.')
    parser.add_argument('--style', default='paper')
    args = parser.parse_args()
    
    fm_, fn_ = get_fonts(args.style)
    print(f"Style: {args.style}, main={fm_}, num={fn_}")
    plt.rcParams['axes.unicode_minus'] = False
    os.makedirs(args.output, exist_ok=True)
    
    df = load_data(args.checkm2, args.bac120, args.ar53, args.keystone)
    
    prefix = os.path.join(args.output, 'Fig3-1_MAG_quality')
    print("\n--- Fig3-1 ---")
    draw_fig3_1(df, prefix, fm_, fn_, args.style)
    
    print("\n--- Stats ---")
    save_stats(df, args.stat_dir, fn_)
    
    print("\nDone!")

if __name__=="__main__":
    main()
