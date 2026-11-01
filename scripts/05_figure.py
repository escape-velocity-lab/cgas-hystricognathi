"""05: 物种树（NCBI Taxonomy 拓扑，无枝长）+ cGAS 4 位点状态 + LQ。
树只是分类学拓扑，不是带时间的系统发育树；仅用于把近缘物种排在一起看。
用法: python 05_figure.py <data/processed> <results>
"""
import sys, os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

PROC, RES = sys.argv[1:3]
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']  # 坑 18
plt.rcParams['axes.unicode_minus'] = False

d = pd.read_csv(os.path.join(RES, 'site_table_v1.csv'), dtype={'taxid': str, '物种级taxid': str})
m = pd.read_csv(os.path.join(PROC, 'orthologs_x_anage.csv'), dtype=str).set_index('gene_id').lineage
d = d[d.物种代表 == '是'].copy()
SITES = ['D431', 'K479', 'L495', 'K498']
NMR = dict(zip(SITES, 'SEYT')); HUM = dict(zip(SITES, 'DKLK'))
MAMMALIA = '40674'

# ---- 由 lineage 建树 ----
children, leaf_of = {}, {}
for _, r in d.iterrows():
    lin = m[str(r.GeneID)].split(',')
    lin = lin[lin.index(MAMMALIA):] + ['leaf:' + r.taxid]
    for p, c in zip(lin, lin[1:]):
        children.setdefault(p, [])
        if c not in children[p]:
            children[p].append(c)
    leaf_of['leaf:' + r.taxid] = r

def collapse(n):  # 去掉只有一个孩子的中间节点
    while n in children and len(children[n]) == 1:
        n = children[n][0]
    return n

def nleaf(n):
    n = collapse(n)
    return 1 if n not in children else sum(nleaf(c) for c in children[n])

order, segs, leafx = [], [], []
def layout(n, depth):
    n = collapse(n)
    if n not in children:
        order.append(n); leafx.append(depth); return len(order) - 1, depth
    ys = []
    for c in sorted(children[n], key=nleaf):
        y, x = layout(c, depth + 1)
        ys.append((y, x))
    for y, x in ys:
        segs.append(((depth, y), (x, y)))          # 横线
    segs.append(((depth, min(y for y, _ in ys)), (depth, max(y for y, _ in ys))))  # 竖线
    return sum(y for y, _ in ys) / len(ys), depth
layout(MAMMALIA, 0)
maxd = max(s[1][0] for s in segs)
N = len(order)

# ---- 画 ----
SURF, INK, INK2 = '#fcfcfb', '#0b0b0b', '#52514e'
C_NMR, C_HUM, C_OTH, C_BAR = '#eb6834', '#c9c8c2', '#ecebe7', '#2a78d6'
fig = plt.figure(figsize=(13, 0.13 * N + 2.2), facecolor=SURF)
gs = fig.add_gridspec(1, 4, width_ratios=[3.2, 2.1, 1.6, 2.6], wspace=0.02, left=0.01, right=0.99,
                      top=1 - 1.5 / (0.13 * N + 2.2), bottom=0.5 / (0.13 * N + 2.2))
axT, axL, axS, axB = (fig.add_subplot(gs[0, i], facecolor=SURF) for i in range(4))
for (x0, y0), (x1, y1) in segs:
    axT.plot([x0, x1], [y0, y1], color=INK2, lw=0.6, solid_capstyle='butt')
for i, n in enumerate(order):
    axT.plot([leafx[i], maxd + 0.5], [i, i], color='#dcdbd6', lw=0.4, ls=':')
axT.set_xlim(-0.5, maxd + 0.6)

hyst = {'Bathyergidae', 'Caviidae', 'Chinchillidae', 'Octodontidae', 'Hystricidae', 'Erethizontidae', 'Thryonomyidae'}
for i, n in enumerate(order):
    r = leaf_of[n]
    lab = f'{r.物种}'
    col = INK if r.裸鼹鼠型命中数 >= 2 else INK2
    w = 'bold' if r.科 in hyst else 'normal'
    axL.text(0.01, i, f'{lab}  ·  {r.目}/{r.科}', va='center', ha='left', fontsize=5.6, color=col, fontweight=w)
    for j, s in enumerate(SITES):
        aa = r[s]
        fc = C_NMR if aa == NMR[s] else (C_HUM if aa == HUM[s] else C_OTH)
        if r['4位点可用'] != '是':
            fc = SURF
        axS.add_patch(Rectangle((j + 0.04, i - 0.42), 0.92, 0.84, facecolor=fc, edgecolor=SURF if r['4位点可用'] == '是' else '#9a9994', lw=0.6,
                                hatch='////' if r['4位点可用'] != '是' else None))
        axS.text(j + 0.5, i, aa, ha='center', va='center', fontsize=5.4,
                 color='#ffffff' if fc == C_NMR else INK)
    if r.LQ == r.LQ:
        axB.barh(i, r.LQ, height=0.62, color=C_BAR, lw=0)
        if r.寿命来源.startswith('借用'):
            axB.text(r.LQ + 0.05, i, '借用近缘种', va='center', fontsize=5, color=INK2)
    else:
        axB.text(0.05, i, '无 AnAge 记录', va='center', fontsize=5, color=INK2)
for ax in (axT, axL, axS, axB):
    ax.set_ylim(N - 0.5, -0.5)
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
for ax in (axT, axL):
    ax.set_xticks([])
axL.set_xlim(0, 1)
axS.set_xlim(0, 4); axS.set_xticks([0.5, 1.5, 2.5, 3.5]); axS.set_xticklabels(SITES, fontsize=7, color=INK2)
axS.xaxis.tick_top(); axS.tick_params(length=0)
xmax = max(4.0, d.LQ.max() * 1.1)
axB.set_xlim(0, xmax)
for v, ls in ((1, '-'), (1.5, '--'), (0.9, ':')):
    axB.axvline(v, color=INK2, lw=0.6, ls=ls, zorder=0)
axB.xaxis.tick_top(); axB.tick_params(labelsize=7, colors=INK2, length=0)
axB.set_xlabel('寿命商 LQ（实线 1 = 同体型平均；虚线 1.5 = 长寿；点线 0.9 = 短寿）', fontsize=7, color=INK2)
axB.xaxis.set_label_position('top')

fig.text(0.01, 1 - 0.35 / (0.13 * N + 2.2),
         f'裸鼹鼠 cGAS 4 个位点在 {N} 种哺乳动物中的状态（人编号；橙 = 裸鼹鼠型 S/E/Y/T，深灰 = 人型 D/K/L/K，浅灰 = 其他，斜线 = 读不出 / 质量不合格）',
         fontsize=10, color=INK, va='top')
fig.text(0.01, 1 - 0.75 / (0.13 * N + 2.2),
         '树 = NCBI Taxonomy 分类拓扑（无枝长，不是时间树）；粗体 = 豪猪亚目（Hystricognathi）各科；LQ = 最大寿命 ÷ 同体重哺乳动物预期寿命（AnAge Build 15 拟合）',
         fontsize=7.5, color=INK2, va='top')
out = os.path.join(RES, 'fig_site_tree_v1.png')
fig.savefig(out, dpi=170, facecolor=SURF)
print('saved', out, 'species', N)
