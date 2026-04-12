#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09_cooccurrence_network.py — Fig3-7 共现网络图
  Fig3-7a: 网络图（节点大小=degree，Keystone深色+标注，边=正/负相关）
  Fig3-7b: Degree vs Betweenness散点图（Keystone筛选标准）

输出:
  figures/paper_en/Fig3-7_network_EN.pdf + .png
  figures/thesis_cn/Fig3-7_network_CN.pdf + .png
  data/processed/Fig3-7_edge_list.tsv
  data/processed/Fig3-7_network_summary.txt

用法:
  python3 scripts/python/09_cooccurrence_network.py \
      --abund  data/raw/abundance.tsv \
      --bac    data/raw/tax_bac120_summary.tsv \
      --ar     data/raw/tax_ar53_summary.tsv \
      --ks     data/raw/keystone_species.txt \
      --nodes  data/raw/network_node_stats.txt \
      --modules data/raw/network_module_info.txt \
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
from matplotlib.patches import FancyArrowPatch
from scipy.stats import spearmanr
import networkx as nx

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
COLOR_NON_KS = '#B0D4E8'  # 浅蓝，非Keystone
COLOR_KS = '#1B3A5C'       # 深蓝，Keystone
COLOR_POS_EDGE = '#CCCCCC'  # 浅灰，正相关
COLOR_NEG_EDGE = '#444444'  # 深灰，负相关

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

# ================================================================
# 数据加载与网络构建
# ================================================================
def load_and_build(args):
    # 丰度
    ab = pd.read_csv(args.abund, sep='\t')
    ab = ab[ab['Genome'] != 'unmapped'].set_index('Genome')
    sample_cols = [c for c in ab.columns if c != 'Genome']
    mat = ab[sample_cols].values
    mags = ab.index.tolist()
    n = len(mags)

    # Spearman |r| > 0.9, p < 0.05
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            r, p = spearmanr(mat[i], mat[j])
            if abs(r) > 0.9 and p < 0.05:
                edges.append((mags[i], mags[j], r, p))

    # 构建networkx图
    G = nx.Graph()
    for m1, m2, r, p in edges:
        G.add_edge(m1, m2, weight=abs(r), r=r, p=p,
                   edge_type='positive' if r > 0 else 'negative')

    # GTDB分类
    bac = pd.read_csv(args.bac, sep='\t')
    ar = pd.read_csv(args.ar, sep='\t')
    gtdb = pd.concat([bac[['user_genome','classification']],
                       ar[['user_genome','classification']]], ignore_index=True)
    gtdb = gtdb.rename(columns={'user_genome':'MAG'})

    def gl(cls, pfx):
        for p in cls.split(';'):
            if p.startswith(pfx): return p.replace(pfx, '')
        return ''
    gtdb['Phylum'] = gtdb['classification'].apply(lambda x: gl(x, 'p__'))
    gtdb['Family'] = gtdb['classification'].apply(lambda x: gl(x, 'f__'))
    gtdb['Genus'] = gtdb['classification'].apply(lambda x: gl(x, 'g__'))
    tax_map = gtdb.set_index('MAG')[['Phylum','Family','Genus']].to_dict('index')

    # Keystone
    ks = pd.read_csv(args.ks, sep='\t')
    ks_mags = set(ks['MAG'].tolist())

    # 模块信息
    mod = pd.read_csv(args.modules, sep='\t')
    mod_map = mod.set_index('MAG')['Module'].to_dict()

    # 节点属性
    for node in G.nodes():
        info = tax_map.get(node, {})
        G.nodes[node]['phylum'] = info.get('Phylum', '')
        G.nodes[node]['family'] = info.get('Family', '')
        G.nodes[node]['genus'] = info.get('Genus', '')
        G.nodes[node]['is_keystone'] = node in ks_mags
        G.nodes[node]['module'] = mod_map.get(node, 0)
        G.nodes[node]['degree'] = G.degree(node)
        G.nodes[node]['betweenness'] = nx.betweenness_centrality(G)[node]

    # 边表DataFrame
    df_edges = pd.DataFrame(edges, columns=['MAG1','MAG2','Spearman_r','p_value'])

    return G, df_edges, ks_mags, ks


# ================================================================
# Fig3-7: 网络图 + 散点图
# ================================================================
def draw_network(G, ks_mags, ks_df, style, outpath):
    fm_, fn_ = get_fonts(style)
    is_cn = (style == 'cn')

    fig, axes = plt.subplots(1, 2, figsize=(20, 10),
                              gridspec_kw={'width_ratios': [1.4, 1]})

    # ======== Fig3-7a: 网络图 ========
    ax = axes[0]
    ax.set_facecolor('white')

    # 布局：先用Kamada-Kawai获取初始位置，再用spring优化
    try:
        pos_init = nx.kamada_kawai_layout(G)
    except Exception:
        pos_init = None
    pos = nx.spring_layout(G, k=6.0/np.sqrt(len(G.nodes())),
                           iterations=300, seed=42, pos=pos_init)

    # 后处理：消除节点重叠（迭代推开太近的节点）
    nodes_list = list(G.nodes())
    positions = np.array([pos[n] for n in nodes_list])
    min_dist = 0.035  # 最小允许距离
    for _ in range(50):
        moved = False
        for i in range(len(nodes_list)):
            for j in range(i+1, len(nodes_list)):
                dx = positions[i][0] - positions[j][0]
                dy = positions[i][1] - positions[j][1]
                dist = np.sqrt(dx*dx + dy*dy)
                if dist < min_dist and dist > 0:
                    # 推开
                    force = (min_dist - dist) / 2
                    nx_ = dx / dist * force
                    ny_ = dy / dist * force
                    positions[i] += [nx_, ny_]
                    positions[j] -= [nx_, ny_]
                    moved = True
        if not moved:
            break
    for i, n in enumerate(nodes_list):
        pos[n] = tuple(positions[i])

    # 绘制边（先画，在节点下面）
    for u, v, d in G.edges(data=True):
        x0, y0 = pos[u]; x1, y1 = pos[v]
        r_val = d['r']
        alpha = min(abs(r_val) * 0.6, 0.5)
        lw = abs(r_val) * 1.5
        color = COLOR_POS_EDGE if r_val > 0 else COLOR_NEG_EDGE
        ax.plot([x0, x1], [y0, y1], color=color, alpha=alpha, lw=lw, zorder=1)

    # 绘制节点
    nodes_sorted = sorted(G.nodes(), key=lambda n: G.nodes[n]['is_keystone'])  # Keystone后画
    for node in nodes_sorted:
        x, y = pos[node]
        d = G.nodes[node]
        deg = G.degree(node)
        is_ks = d['is_keystone']

        # 大小：按degree缩放
        size = max(deg * 25, 40) if is_ks else max(deg * 15, 20)

        # 颜色
        color = COLOR_KS if is_ks else COLOR_NON_KS
        edge_color = 'white'
        edge_width = 0.8

        ax.scatter(x, y, s=size, c=color, edgecolors=edge_color,
                  linewidths=edge_width, zorder=3, alpha=0.9)

        # Keystone标注名称
        if is_ks:
            genus = d['genus'] or d['family'] or node.replace('Mx_All_','')
            # 截短名称
            if len(genus) > 18: genus = genus[:16] + '..'
            ax.annotate(genus, (x, y), fontsize=5.5, fontfamily=fn_,
                       fontstyle='italic', fontweight='bold',
                       xytext=(5, 5), textcoords='offset points',
                       color=COLOR_KS, zorder=4,
                       bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                                edgecolor='none', alpha=0.7))

    ax.set_xlim(ax.get_xlim()[0]-0.05, ax.get_xlim()[1]+0.05)
    ax.set_ylim(ax.get_ylim()[0]-0.05, ax.get_ylim()[1]+0.05)
    ax.axis('off')

    # 网络统计标注
    n_pos = sum(1 for _, _, d in G.edges(data=True) if d['r'] > 0)
    n_neg = sum(1 for _, _, d in G.edges(data=True) if d['r'] < 0)
    stats_text = f"Nodes: {len(G.nodes())} | Edges: {len(G.edges())}\n"
    stats_text += f"|r| > 0.9, p < 0.05"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
           fontsize=8, fontfamily=fn_, va='top', color='#555',
           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#ddd', alpha=0.8))

    # 图例
    leg_y = 0.02
    ax.scatter([], [], s=80, c=COLOR_NON_KS, edgecolors='white', linewidths=0.5,
              label='Non-Keystone' if not is_cn else '非关键物种')
    ax.scatter([], [], s=150, c=COLOR_KS, edgecolors='white', linewidths=0.5,
              label='Keystone' if not is_cn else '关键物种')
    ax.plot([], [], color=COLOR_POS_EDGE, lw=2, alpha=0.5,
           label='Co-Presence' if not is_cn else '正相关')
    ax.plot([], [], color=COLOR_NEG_EDGE, lw=2, alpha=0.5,
           label='Mutual Exclusive' if not is_cn else '负相关')
    leg = ax.legend(loc='lower left', fontsize=7, frameon=True,
                    facecolor='white', edgecolor='#ddd', framealpha=0.9,
                    prop={'family': fm_ if is_cn else fn_})

    title_a = 'a. MAG共现网络' if is_cn else 'a. MAG Co-occurrence Network'
    ax.set_title(title_a, fontsize=12, fontweight='bold',
                fontfamily=fm_ if is_cn else fn_, pad=10)

    # ======== Fig3-7b: Degree vs Betweenness 散点图 ========
    ax2 = axes[1]
    ax2.set_facecolor('white')

    degrees = [G.degree(n) for n in G.nodes()]
    betweenness_raw = nx.betweenness_centrality(G)
    # 转为非标准化：乘以 (n-1)(n-2)/2
    n_nodes = len(G.nodes())
    scale_factor = (n_nodes - 1) * (n_nodes - 2) / 2
    betweenness = {n: betweenness_raw[n] * scale_factor for n in G.nodes()}

    for node in G.nodes():
        d = G.degree(node)
        b = betweenness[node]
        is_ks = G.nodes[node]['is_keystone']
        color = COLOR_KS if is_ks else COLOR_NON_KS
        size = 80 if is_ks else 30
        alpha = 0.9 if is_ks else 0.5
        zorder = 3 if is_ks else 2
        ax2.scatter(d, b, s=size, c=color, alpha=alpha, edgecolors='white',
                   linewidths=0.5, zorder=zorder)

        # Keystone标注
        if is_ks:
            genus = G.nodes[node]['genus'] or G.nodes[node]['family'] or node.replace('Mx_All_','')
            if len(genus) > 15: genus = genus[:13] + '..'
            ax2.annotate(genus, (d, b), fontsize=5, fontfamily=fn_,
                        fontstyle='italic', fontweight='bold',
                        xytext=(4, 4), textcoords='offset points',
                        color=COLOR_KS, zorder=4)

    # Keystone筛选阈值线
    ax2.axvline(x=10, color='#E74C3C', linestyle='--', alpha=0.4, lw=1)
    ax2.axhline(y=200, color='#E74C3C', linestyle='--', alpha=0.4, lw=1)
    thresh_txt = 'Degree ≥ 10 or\nBetweenness ≥ 200' if not is_cn else '度 ≥ 10 或\n介数中心性 ≥ 200'
    ax2.text(0.97, 0.97, thresh_txt, transform=ax2.transAxes,
            fontsize=7, fontfamily=fm_ if is_cn else fn_, va='top', ha='right',
            color='#E74C3C', alpha=0.7)

    ax2.set_xlabel('Degree' if not is_cn else '度 (Degree)',
                   fontsize=10, fontfamily=fm_ if is_cn else fn_)
    ax2.set_ylabel('Betweenness centrality' if not is_cn else '介数中心性 (Betweenness)',
                   fontsize=10, fontfamily=fm_ if is_cn else fn_)
    ax2.tick_params(labelsize=8)
    for spine in ax2.spines.values():
        spine.set_color('#ccc')

    # 图例
    ax2.scatter([], [], s=80, c=COLOR_KS, edgecolors='white', linewidths=0.5,
               label='Keystone' if not is_cn else '关键物种')
    ax2.scatter([], [], s=30, c=COLOR_NON_KS, edgecolors='white', linewidths=0.5,
               label='Non-Keystone' if not is_cn else '非关键物种')
    ax2.legend(loc='upper left', fontsize=7, frameon=True,
              facecolor='white', edgecolor='#ddd',
              prop={'family': fm_ if is_cn else fn_})

    title_b = 'b. 关键物种筛选' if is_cn else 'b. Keystone Species Selection'
    ax2.set_title(title_b, fontsize=12, fontweight='bold',
                 fontfamily=fm_ if is_cn else fn_, pad=10)

    plt.tight_layout(w_pad=3)

    dpi = 300 if is_cn else 600
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight', facecolor='white')
    png = outpath.rsplit('.', 1)[0] + '.png'
    plt.savefig(png, dpi=min(dpi, 300), bbox_inches='tight', facecolor='white')
    print(f"  Saved: {outpath} + .png")
    plt.close()


# ================================================================
# 数据导出
# ================================================================
def export_stats(G, df_edges, ks_mags, stat_dir):
    os.makedirs(stat_dir, exist_ok=True)

    # 边表
    df_edges.to_csv(os.path.join(stat_dir, 'Fig3-7_edge_list.tsv'),
                    sep='\t', index=False)

    # 网络汇总
    n_pos = sum(1 for _, _, d in G.edges(data=True) if d['r'] > 0)
    n_neg = sum(1 for _, _, d in G.edges(data=True) if d['r'] < 0)
    betweenness_raw = nx.betweenness_centrality(G)
    n_nodes = len(G.nodes())
    scale = (n_nodes - 1) * (n_nodes - 2) / 2

    with open(os.path.join(stat_dir, 'Fig3-7_network_summary.txt'), 'w') as f:
        f.write("=== Co-occurrence Network Summary ===\n\n")
        f.write(f"Threshold: |Spearman r| > 0.9, p < 0.05\n")
        f.write(f"Nodes: {n_nodes}\n")
        f.write(f"Edges: {len(G.edges())} (Positive: {n_pos}, Negative: {n_neg})\n")
        f.write(f"Density: {nx.density(G):.4f}\n")
        f.write(f"Components: {nx.number_connected_components(G)}\n")
        avg_deg = np.mean([G.degree(n) for n in G.nodes()])
        f.write(f"Average degree: {avg_deg:.1f}\n")
        f.write(f"Keystone species in network: {sum(1 for n in G.nodes() if G.nodes[n]['is_keystone'])}\n\n")

        f.write(f"{'MAG':<15} {'Degree':>7} {'Between':>10} {'Phylum':<22} {'Genus':<20} {'KS':>3}\n")
        f.write("-" * 80 + "\n")
        for node in sorted(G.nodes(), key=lambda n: G.degree(n), reverse=True):
            d = G.nodes[node]
            bt = betweenness_raw[node] * scale
            ks_flag = 'Y' if d['is_keystone'] else ''
            f.write(f"{node:<15} {G.degree(node):>7} {bt:>10.1f} "
                    f"{d['phylum']:<22} {d['genus']:<20} {ks_flag:>3}\n")

    print(f"  Stats exported to: {stat_dir}")


# ================================================================
# Main
# ================================================================
def main():
    parser = argparse.ArgumentParser(description='Fig3-7 Co-occurrence Network')
    parser.add_argument('--abund',  required=True)
    parser.add_argument('--bac',    required=True)
    parser.add_argument('--ar',     required=True)
    parser.add_argument('--ks',     required=True)
    parser.add_argument('--nodes',  required=True)
    parser.add_argument('--modules',required=True)
    parser.add_argument('--outdir', default='figures')
    parser.add_argument('--statdir',default='data/processed')
    args = parser.parse_args()

    print(f"Fonts: sans={F_SANS}, serif={F_SERIF}, cn={F_CN}")

    G, df_edges, ks_mags, ks_df = load_and_build(args)
    print(f"Network: {len(G.nodes())} nodes, {len(G.edges())} edges")
    print(f"Keystone in network: {sum(1 for n in G.nodes() if G.nodes[n]['is_keystone'])}")

    en_dir = os.path.join(args.outdir, 'paper_en')
    cn_dir = os.path.join(args.outdir, 'thesis_cn')
    os.makedirs(en_dir, exist_ok=True)
    os.makedirs(cn_dir, exist_ok=True)

    print("\n=== EN version ===")
    draw_network(G, ks_mags, ks_df, 'en',
                 os.path.join(en_dir, 'Fig3-7_network_EN.pdf'))

    print("\n=== CN version ===")
    draw_network(G, ks_mags, ks_df, 'cn',
                 os.path.join(cn_dir, 'Fig3-7_network_CN.pdf'))

    print("\n=== Exporting stats ===")
    export_stats(G, df_edges, ks_mags, args.statdir)

    print("\nAll done!")

if __name__ == '__main__':
    main()
