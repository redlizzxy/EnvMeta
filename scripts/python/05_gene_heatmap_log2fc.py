#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_gene_heatmap_log2fc.py v6 — Fig2-8 + Fig2-9
"""

import argparse, os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.patches as mpatches
import matplotlib.font_manager as fm
from matplotlib.colors import TwoSlopeNorm
from matplotlib.backends.backend_agg import FigureCanvasAgg
import warnings
warnings.filterwarnings('ignore')

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

GROUP_ORDER = ['CK','A','B']
SAMPLE_MAP = {'2_1':'CK_1','2_5':'CK_2','2_7':'CK_3',
              '8_1':'A_1','8_2':'A_2','8_3':'A_3','8_4':'A_4',
              '9_1':'B_1','9_2':'B_2','9_3':'B_3'}

# Scheme 2 Enhanced
FUNC_BG = {
    'arsenic': {
        'Arsenate reduction':'#F5D8D8','Arsenite oxidation':'#E0A8A0',
        'Resp. arsenate red.':'#C87870','As transport/detox':'#F0D0B8',
        'As methylation':'#D8A880','As regulation':'#C08060',
    },
    'sulfur': {
        'Assim. sulfate red.':'#D0F0E0','Dissim. sulfate red.':'#88D0A8',
        'Sulfide oxidation':'#48A878','Thiosulfate metab.':'#E0E8B0',
    },
    'nitrogen': {
        'Nitrate reduction':'#D8E4F5','Nitrite reduction':'#A0B8D8',
        'NO reduction':'#6888B8','N2O reduction':'#E0D0E8',
        'Ammonia oxidation':'#B8A8D0','N fixation':'#8878A8',
    },
    'iron': {
        'Fe transport':'#F5E0C0','Fe uptake regulation':'#C8A070',
    },
}

# 显示名（N₂O用mathtext）
FUNC_DISP_EN = {
    'Arsenate reduction':'Arsenate reduction','Arsenite oxidation':'Arsenite oxidation',
    'Resp. arsenate red.':'Resp. arsenate red.','As transport/detox':'As transport/detox',
    'As methylation':'As methylation','As regulation':'As regulation',
    'Nitrate reduction':'Nitrate reduction','Nitrite reduction':'Nitrite reduction',
    'NO reduction':'NO reduction','N2O reduction':r'N$_2$O reduction',
    'Ammonia oxidation':'Ammonia oxidation','N fixation':'N fixation',
    'Assim. sulfate red.':'Assim. sulfate red.','Dissim. sulfate red.':'Dissim. sulfate red.',
    'Sulfide oxidation':'Sulfide oxidation','Thiosulfate metab.':'Thiosulfate metab.',
    'Fe transport':'Fe transport','Fe uptake regulation':'Fe uptake regulation',
}
FUNC_DISP_CN = {
    'Arsenate reduction':'砷酸盐还原','Arsenite oxidation':'亚砷酸盐氧化',
    'Resp. arsenate red.':'呼吸型砷酸盐还原','As transport/detox':'砷转运与解毒',
    'As methylation':'砷甲基化','As regulation':'砷调控',
    'Nitrate reduction':'硝酸盐还原','Nitrite reduction':'亚硝酸盐还原',
    'NO reduction':'NO还原','N2O reduction':r'N$_2$O还原',
    'Ammonia oxidation':'氨氧化','N fixation':'固氮',
    'Assim. sulfate red.':'硫酸盐同化还原','Dissim. sulfate red.':'硫酸盐异化还原',
    'Sulfide oxidation':'硫化物氧化','Thiosulfate metab.':'硫代硫酸盐代谢',
    'Fe transport':'铁转运系统','Fe uptake regulation':'铁摄取调控',
}

KO_DB = {
    'arsenic': {
        'K00537':('arsC-grx','Arsenate reduction'),'K03741':('arsC-trx','Arsenate reduction'),
        'K22547':('HAC1','Arsenate reduction'),'K18701':('arsC-myc','Arsenate reduction'),
        'K08355':('aoxA','Arsenite oxidation'),'K08356':('aoxB','Arsenite oxidation'),
        'K28466':('arrA','Resp. arsenate red.'),'K28467':('arrB','Resp. arsenate red.'),
        'K03893':('arsB','As transport/detox'),'K03325':('ACR3','As transport/detox'),
        'K18064':('ACR2','As transport/detox'),'K01551':('arsA','As transport/detox'),
        'K25223':('arsJ','As transport/detox'),'K11811':('arsH','As transport/detox'),
        'K07755':('AS3MT','As methylation'),'K25224':('gapdh','As methylation'),
        'K03892':('arsR','As regulation'),
    },
    'nitrogen': {
        'K00370':('narG','Nitrate reduction'),'K00371':('narH','Nitrate reduction'),
        'K00374':('narI','Nitrate reduction'),'K02567':('napA','Nitrate reduction'),
        'K02568':('napB','Nitrate reduction'),'K00367':('narB','Nitrate reduction'),
        'K00368':('nirK','Nitrite reduction'),'K15864':('nirS','Nitrite reduction'),
        'K00366':('nirA','Nitrite reduction'),'K15877':('CYP55','Nitrite reduction'),
        'K04561':('norB','NO reduction'),
        'K00376':('nosZ','N2O reduction'),
        'K10944':('pmoA-amoA','Ammonia oxidation'),'K10945':('pmoB-amoB','Ammonia oxidation'),
        'K10946':('pmoC-amoC','Ammonia oxidation'),
        'K02586':('nifD','N fixation'),'K02591':('nifK','N fixation'),
    },
    'sulfur': {
        'K00380':('cysJ','Assim. sulfate red.'),'K00381':('cysI','Assim. sulfate red.'),
        'K00390':('cysH','Assim. sulfate red.'),'K00392':('sir','Assim. sulfate red.'),
        'K00394':('aprA','Dissim. sulfate red.'),'K00395':('aprB','Dissim. sulfate red.'),
        'K11180':('dsrA','Dissim. sulfate red.'),'K11181':('dsrB','Dissim. sulfate red.'),
        'K17218':('sqr','Sulfide oxidation'),'K17222':('soxA','Sulfide oxidation'),
        'K17223':('soxX','Sulfide oxidation'),'K17224':('soxB','Sulfide oxidation'),
        'K17226':('soxY','Sulfide oxidation'),'K17227':('soxZ','Sulfide oxidation'),
        'K01011':('TST','Thiosulfate metab.'),
    },
    'iron': {
        'K02013':('ABC.FEV.A','Fe transport'),'K02014':('TC.FEV.OM','Fe transport'),
        'K02015':('ABC.FEV.P','Fe transport'),'K02016':('ABC.FEV.S','Fe transport'),
        'K02012':('afuA','Fe transport'),'K07243':('FTR','Fe transport'),
        'K03711':('fur','Fe uptake regulation'),'K03832':('tonB','Fe uptake regulation'),
    }
}

METAB_EN = {'arsenic':'Arsenic','nitrogen':'Nitrogen','sulfur':'Sulfur','iron':'Iron'}
METAB_CN = {'arsenic':'砷代谢','nitrogen':'氮代谢','sulfur':'硫代谢','iron':'铁代谢'}

# ==============================================================================
# 数据IO
# ==============================================================================
def fix_spf_load(f):
    with open(f,'r') as fh: lines=fh.readlines()
    hdr=lines[0].strip().split('\t')
    sn=hdr[2:] if ('Unannotated' in hdr[0] or 'KEGG' in hdr[0]) else hdr[1:]
    kd={}
    for l in lines[1:]:
        if not l.strip(): continue
        p=l.strip().split('\t'); ki=None; vs=-1
        for i,x in enumerate(p):
            if x.startswith('K') and len(x)==6 and x[1:].isdigit(): ki=x; vs=i+1; break
        if not ki: continue
        v=[]
        for j in range(vs,len(p)):
            try: v.append(float(p[j]))
            except: v.append(0.0)
        while len(v)<len(sn): v.append(0.0)
        v=v[:len(sn)]
        kd[ki]=[(a+b)/2 for a,b in zip(kd[ki],v)] if ki in kd else v
    df=pd.DataFrame(kd,index=sn).T
    print(f"Loaded {len(df)} KOs x {len(df.columns)} samples")
    return df

def load_meta(f):
    m=pd.read_csv(f,sep='\t'); m.columns=m.columns.str.strip()
    m['SampleID']=m['SampleID'].str.strip()
    s2g=dict(zip(m['SampleID'],m['Group'].str.strip()))
    g2s={}
    for s,g in s2g.items(): g2s.setdefault(g,[]).append(s)
    return m,s2g,g2s

def group_mean(df,g2s):
    r=pd.DataFrame()
    for g in GROUP_ORDER:
        c=[x for x in g2s.get(g,[]) if x in df.columns]
        if c: r[g]=df[c].mean(axis=1)
    return r

def calc_stats(df,g2s,g1,g2):
    s1=[c for c in g2s.get(g1,[]) if c in df.columns]
    s2=[c for c in g2s.get(g2,[]) if c in df.columns]
    res=[]
    for ko in df.index:
        v1=df.loc[ko,s1].values.astype(float); v2=df.loc[ko,s2].values.astype(float)
        v1p=v1[v1>0]; v2p=v2[v2>0]
        if len(v1p)==0 or len(v2p)==0:
            res.append({'KO':ko,'log2FC':np.nan,'pvalue':1.0}); continue
        m1,m2=np.mean(v1p),np.mean(v2p)
        lfc=np.log2(m2/m1) if m1>0 and m2>0 else np.nan
        p=stats.ttest_ind(v1p,v2p,equal_var=False)[1] if len(v1p)>1 and len(v2p)>1 else 1.0
        res.append({'KO':ko,'log2FC':lfc,'pvalue':p})
    rdf=pd.DataFrame(res).set_index('KO')
    n=(~rdf['pvalue'].isna()).sum()
    rdf['padj']=(rdf['pvalue']*n).clip(upper=1.0) if n>0 else 1.0
    return rdf

# ==============================================================================
# 测量文本宽度（数据坐标单位）
# ==============================================================================
def measure_label_width(ax, fig, labels, fontsize, fontfamily):
    """测量labels中最宽文本在数据坐标中的宽度"""
    canvas = FigureCanvasAgg(fig)
    renderer = canvas.get_renderer()
    max_w_data = 0
    for lb in labels:
        t = ax.text(0, 0, lb, fontsize=fontsize, fontfamily=fontfamily, fontweight='bold')
        bb = t.get_window_extent(renderer=renderer)
        # 转为数据坐标宽度
        inv = ax.transData.inverted()
        p0 = inv.transform((0, 0))
        p1 = inv.transform((bb.width, 0))
        w_data = abs(p1[0] - p0[0])
        max_w_data = max(max_w_data, w_data)
        t.remove()
    return max_w_data

# ==============================================================================
# Fig2-8: 热图面板
# ==============================================================================
def draw_panel(ax, fig, df_z, ko_info, metab_key, fm_, fn_, style, letter, fixed_bw=None):
    is_cn = (style=='thesis')
    func_bg = FUNC_BG.get(metab_key, {})

    # 按功能排序
    func_seen=[]
    for ko in df_z.index:
        f=ko_info.get(ko,('?','?'))[1]
        if f not in func_seen: func_seen.append(f)
    sorted_kos=[]
    for f in func_seen:
        sorted_kos.extend([ko for ko in df_z.index if ko_info.get(ko,('?','?'))[1]==f])
    df=df_z.loc[sorted_kos]

    ng=len(df); nc=len(df.columns)
    vmax=max(abs(df.values.min()),abs(df.values.max()),0.01)
    norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax)

    # 准备标签
    labels=[]
    for ko in df.index:
        gene=ko_info.get(ko,('?','?'))[0]
        labels.append(f"{ko} ({gene})")

    # 热图只占axes右侧部分，左侧留给色块
    # 色块占axes宽度的45%，热图占55%
    label_frac = 0.45
    heat_left = label_frac
    heat_right = 1.0

    # 在axes坐标中画热图（用inset_axes方式或直接用extent+xlim控制）
    # 简单方案：热图extent设为[0, nc, ng, 0]，然后通过xlim把左侧留出来
    # 色块宽度 = nc * label_frac / (1 - label_frac)
    bw = nc * label_frac / (1 - label_frac)

    im=ax.imshow(df.values,cmap='RdBu_r',norm=norm,aspect='auto',extent=[0,nc,ng,0])

    ax.set_xticks([i+0.5 for i in range(nc)])
    ax.set_xticklabels(df.columns, fontsize=8, fontfamily=fn_)
    ax.set_yticks([])

    # 色块+文字
    for i, ko in enumerate(df.index):
        func=ko_info.get(ko,('?','?'))[1]
        bg=func_bg.get(func,'#F0F0F0')
        ax.add_patch(plt.Rectangle((-bw,i),bw,1,
            facecolor=bg,edgecolor='none',clip_on=False,zorder=3))
        ax.text(-bw/2, i+0.5, labels[i], fontsize=6.5,
               ha='center',va='center',fontfamily=fn_,
               fontweight='bold',color='#222222',
               clip_on=False,zorder=4)

    ax.axvline(x=0,color='#999',linewidth=0.5,zorder=5)
    ax.set_xlim(-bw,nc); ax.set_ylim(ng,0)

    if is_cn:
        ax.set_title(f"({letter}) {METAB_CN.get(metab_key,'')}",
                    fontsize=10,fontweight='bold',fontfamily=fm_,pad=6)
    else:
        ax.set_title(f"({letter}) {METAB_EN.get(metab_key,'')} Metabolism",
                    fontsize=10,fontweight='bold',fontfamily=fn_,pad=6)

    for sp in ['top','right','left']: ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.4)
    ax.tick_params(left=False,bottom=True,length=2)

    return im, func_seen, func_bg

def draw_fig2_8(df_all, g2s, prefix, fm_, fn_, style):
    metabs=['arsenic','sulfur','nitrogen','iron']
    letters=['a','b','c','d']
    is_cn=(style=='thesis')

    fig = plt.figure(figsize=(17, 16))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.55],
                          hspace=0.28, wspace=0.5)

    axes_heat=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),
               fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1])]
    ax_leg=fig.add_subplot(gs[:,2])
    ax_leg.axis('off')

    last_im=None
    all_legends={}

    for idx,metab in enumerate(metabs):
        ax=axes_heat[idx]
        ko_info=KO_DB[metab]
        present=[ko for ko in ko_info if ko in df_all.index]
        if not present:
            ax.text(0.5,0.5,'No data',transform=ax.transAxes,ha='center'); continue

        df_sub=df_all.loc[present]
        df_gm=group_mean(df_sub,g2s)
        df_z=df_gm.apply(lambda x:(x-x.mean())/x.std() if x.std()>0 else x*0, axis=1)

        im,func_seen,func_bg=draw_panel(ax,fig,df_z,ko_info,metab,fm_,fn_,style,letters[idx])
        last_im=im

        items=[]
        for f in func_seen:
            disp = FUNC_DISP_CN.get(f,f) if is_cn else FUNC_DISP_EN.get(f,f)
            items.append((disp,func_bg.get(f,'#F0F0F0')))
        all_legends[metab]=items

    # 图例：右侧列，靠右上角
    ly=0.98
    for metab in metabs:
        if metab not in all_legends: continue
        entries=all_legends[metab]
        mt=METAB_CN.get(metab,metab) if is_cn else METAB_EN.get(metab,metab)

        ax_leg.text(0.15, ly, mt, fontsize=9, fontweight='bold',
                   fontfamily=fm_ if is_cn else fn_,
                   transform=ax_leg.transAxes, va='top')
        ly -= 0.04

        for disp, clr in entries:
            ax_leg.add_patch(mpatches.FancyBboxPatch(
                (0.15, ly-0.016), 0.06, 0.025,
                boxstyle="square,pad=0", facecolor=clr, edgecolor='none',
                transform=ax_leg.transAxes, zorder=10, clip_on=False))
            ax_leg.text(0.23, ly-0.003, disp, fontsize=7,
                       fontfamily=fm_ if is_cn else fn_,
                       transform=ax_leg.transAxes, va='top')
            ly -= 0.032

        ly -= 0.025

    # Colorbar：自动定位在4个热图右侧
    if last_im:
        cbar = fig.colorbar(last_im, ax=axes_heat, shrink=0.35, pad=0.06, aspect=25)
        cbar.set_label('Z-score (TPM)', fontsize=8, fontfamily=fn_)
        cbar.ax.tick_params(labelsize=7)

    for ext in ['pdf','png']:
        fp=f"{prefix}_heatmap_combined.{ext}"
        fig.savefig(fp,dpi=300,bbox_inches='tight',facecolor='white')
        print(f"  -> {fp}")
    plt.close(fig)

# ==============================================================================
# Fig2-9: log2FC
# ==============================================================================
def draw_fig2_9(df_all, g2s, prefix, stat_dir, fm_, fn_, style):
    is_cn=(style=='thesis')
    metabs=['arsenic','sulfur','nitrogen','iron']
    comps=[('CK','A'),('CK','B'),('A','B')]

    C_UP_SIG='#C0392B'; C_DOWN_SIG='#2980B9'
    C_UP_NS='#F0C0C0'; C_DOWN_NS='#C0D8F0'

    fig, axes=plt.subplots(len(metabs), len(comps), figsize=(15, 16))
    all_stats=[]

    for ri,metab in enumerate(metabs):
        ko_info=KO_DB[metab]
        present=[ko for ko in ko_info if ko in df_all.index]
        if not present: continue
        df_sub=df_all.loc[present]

        for ci,(g1,g2) in enumerate(comps):
            ax=axes[ri,ci]
            sdf=calc_stats(df_sub,g2s,g1,g2).dropna(subset=['log2FC']).sort_values('log2FC')

            ylabels=[ko_info.get(ko,('?','?'))[0] for ko in sdf.index]
            colors=[]
            for _,row in sdf.iterrows():
                if row['padj']<0.05:
                    colors.append(C_UP_SIG if row['log2FC']>0 else C_DOWN_SIG)
                else:
                    colors.append(C_UP_NS if row['log2FC']>0 else C_DOWN_NS)

            y_pos=range(len(sdf))
            ax.barh(y_pos,sdf['log2FC'],color=colors,height=0.7,edgecolor='none')
            ax.set_yticks(y_pos)
            ax.set_yticklabels(ylabels,fontsize=6,fontfamily=fn_,fontstyle='italic')

            ax.axvline(x=1,color='#888',ls='--',lw=0.6,alpha=0.7)
            ax.axvline(x=-1,color='#888',ls='--',lw=0.6,alpha=0.7)
            ax.axvline(x=0,color='black',lw=0.4)

            for i,(ko,row) in enumerate(sdf.iterrows()):
                sig='***' if row['padj']<0.001 else '**' if row['padj']<0.01 else '*' if row['padj']<0.05 else ''
                if sig:
                    xp=row['log2FC']+(0.04 if row['log2FC']>0 else -0.04)
                    ax.text(xp,i,sig,fontsize=5,va='center',
                           ha='left' if row['log2FC']>0 else 'right',
                           fontfamily=fn_,color='#333')

            if ri==0:
                ax.set_title(f"{g1} vs {g2}",fontsize=9,fontweight='bold',fontfamily=fn_)
            if ci==0:
                ml=METAB_CN[metab] if is_cn else METAB_EN[metab]
                ax.set_ylabel(ml,fontsize=9,fontweight='bold',fontfamily=fm_ if is_cn else fn_)
            if ri==len(metabs)-1:
                ax.set_xlabel(r'log$_2$FC',fontsize=8,fontfamily=fn_)

            ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
            ax.tick_params(labelsize=6)

            for ko,row in sdf.iterrows():
                gene=ko_info.get(ko,('?','?'))[0]; func=ko_info.get(ko,('?','?'))[1]
                all_stats.append({
                    'Metabolism':metab,'Comparison':f"{g1}_vs_{g2}",
                    'KO':ko,'Gene':gene,'Function':func,
                    'log2FC':row['log2FC'],'pvalue':row['pvalue'],'padj':row['padj']
                })

    fig.tight_layout(h_pad=2.0,w_pad=1.5,rect=[0,0,0.88,1])

    # 图例右上角：醒目色块+清晰标签
    lx=0.90; ly=0.97
    fig.text(lx, ly, 'Significance', fontsize=9, fontweight='bold',
            fontfamily=fn_, transform=fig.transFigure, va='top')
    ly -= 0.03

    legend_items=[
        (C_UP_SIG,'Up-regulated\n(P < 0.05)'),
        (C_UP_NS,'Up-regulated\n(P ≥ 0.05)'),
        (C_DOWN_SIG,'Down-regulated\n(P < 0.05)'),
        (C_DOWN_NS,'Down-regulated\n(P ≥ 0.05)'),
    ]
    for clr,lbl in legend_items:
        fig.patches.append(mpatches.FancyBboxPatch(
            (lx, ly-0.012), 0.018, 0.018,
            boxstyle="square,pad=0", facecolor=clr, edgecolor='none',
            transform=fig.transFigure, zorder=10))
        fig.text(lx+0.025, ly, lbl, fontsize=6.5,
                fontfamily=fn_, transform=fig.transFigure, va='top',
                linespacing=1.3)
        ly -= 0.045

    ly -= 0.01
    fig.text(lx, ly, r'Dashed line: |log$_2$FC| = 1',
            fontsize=6.5, fontfamily=fn_, transform=fig.transFigure,
            va='top', color='#666')

    for ext in ['pdf','png']:
        fp=f"{prefix}_log2fc_combined.{ext}"
        fig.savefig(fp,dpi=300,bbox_inches='tight',facecolor='white')
        print(f"  -> {fp}")
    plt.close(fig)

    os.makedirs(stat_dir,exist_ok=True)
    sf=os.path.join(stat_dir,'gene_log2fc_all_comparisons.txt')
    pd.DataFrame(all_stats).to_csv(sf,sep='\t',index=False)
    print(f"  -> {sf}")

# ==============================================================================
# 单通路overview
# ==============================================================================
def draw_overviews(df_all, g2s, prefix, fm_, fn_, style):
    for metab, ko_info in KO_DB.items():
        present=[ko for ko in ko_info if ko in df_all.index]
        if not present: continue
        df_sub=df_all.loc[present]
        df_z=df_sub.apply(lambda x:(x-x.mean())/x.std() if x.std()>0 else x*0, axis=1)
        df_z.columns=[SAMPLE_MAP.get(c,c) for c in df_z.columns]
        ng=len(df_z)
        fig,ax=plt.subplots(figsize=(9,max(3.5,ng*0.4+1)))
        im,_,_=draw_panel(ax,fig,df_z,ko_info,metab,fm_,fn_,style,'',fixed_bw=None)
        t=ax.get_title().replace('() ','')
        ax.set_title(t,fontsize=11,fontweight='bold',fontfamily=fm_ if style=='thesis' else fn_)
        plt.colorbar(im,ax=ax,shrink=0.6,label='Z-score',pad=0.02)
        fig.tight_layout()
        for ext in ['pdf','png']:
            fig.savefig(f"{prefix}_{metab}_overview.{ext}",dpi=300,bbox_inches='tight',facecolor='white')
        plt.close(fig)
        print(f"  -> {metab} overview")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i','--input',required=True)
    parser.add_argument('--metadata',required=True)
    parser.add_argument('-o','--output',default='.')
    parser.add_argument('--stat_dir',default='.')
    parser.add_argument('--style',default='paper')
    args=parser.parse_args()

    fm_,fn_=get_fonts(args.style)
    print(f"Style: {args.style}, main={fm_}, num={fn_}")
    plt.rcParams['axes.unicode_minus']=False

    os.makedirs(args.output,exist_ok=True)
    df_all=fix_spf_load(args.input)
    _,_,g2s=load_meta(args.metadata)

    p8=os.path.join(args.output,'Fig2-8')
    p9=os.path.join(args.output,'Fig2-9')
    ps=os.path.join(args.output,'FigS')

    print("\n--- Fig2-8 ---")
    draw_fig2_8(df_all,g2s,p8,fm_,fn_,args.style)
    print("\n--- Fig2-9 ---")
    draw_fig2_9(df_all,g2s,p9,args.stat_dir,fm_,fn_,args.style)
    print("\n--- Overview ---")
    draw_overviews(df_all,g2s,ps,fm_,fn_,args.style)
    print("\nDone!")

if __name__=="__main__":
    main()
