"""19b（任务 #19）终版：由 13_t15_figures.py 派生，只改输出位置 / 文件名、去掉图内标题（标题放图注）。计算完全相同。
原说明：13（任务 #15b）: 短报告的两张图（英文标注，准备投 bioRxiv）。
图 1：豪猪亚目 + 啮齿目外群的 Upham 2019 时间树，每个物种 4 个位点色块 + LQ 条；
      标出祖先重建（Mk-ER，在全部 277 种哺乳动物树上做，同 #13）中裸鼹鼠型状态出现的那段枝（P 从 <0.5 变成 ≥0.5）。
图 2：4 个位点命中数 × LQ（全部有 LQ 且 4 位点可读的物种），豪猪亚目高亮。
输入: results/t13_species_sites.csv、data/raw/t13/tree/
输出: results/figures/fig1_hystricognathi_tree.{png,pdf}、fig2_hits_vs_LQ.{png,pdf}、fig1_gain_branches.csv
用法: python 13_t15_figures.py <proj_root>
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t15_lib as T

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
OUT = R('manuscript', 'figures')
os.makedirs(OUT, exist_ok=True)
SITES = ['431', '479', '495', '498']
NMR = dict(zip(SITES, 'SEYT'))
HUM = dict(zip(SITES, 'DKLK'))
HYST_FAM = {'Bathyergidae', 'Caviidae', 'Chinchillidae', 'Octodontidae', 'Erethizontidae', 'Thryonomyidae', 'Petromuridae',
            'Echimyidae', 'Hystricidae', 'Ctenomyidae', 'Dinomyidae', 'Dasyproctidae', 'Cuniculidae'}
OUTGROUPS = ['Ctenodactylus gundi', 'Castor canadensis', 'Dipodomys ordii', 'Jaculus jaculus', 'Nannospalax galili', 'Mus musculus',
             'Rattus norvegicus', 'Peromyscus leucopus', 'Muscardinus avellanarius', 'Sciurus vulgaris', 'Marmota marmota']
COMMON = {'Heterocephalus glaber': 'naked mole-rat', 'Fukomys damarensis': 'Damaraland mole-rat', 'Fukomys anselli': "Ansell's mole-rat",
          'Cavia porcellus': 'guinea pig', 'Cavia tschudii': 'montane guinea pig', 'Cavia aperea': 'Brazilian guinea pig',
          'Chinchilla lanigera': 'chinchilla', 'Octodon degus': 'degu', 'Hydrochoerus hydrochaeris': 'capybara',
          'Ctenodactylus gundi': 'gundi', 'Mus musculus': 'house mouse', 'Homo sapiens': 'human'}

# 配色（dataviz 参考色板）：裸鼹鼠型 = 橙，人型 = 深灰，其他 = 浅灰；豪猪亚目高亮 = 蓝
C_NMR, C_HUM, C_OTH, C_HYS, C_BG = '#eb6834', '#52514e', '#d9d8d4', '#2a78d6', '#9a9890'
INK, INK2 = '#0b0b0b', '#52514e'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.edgecolor': INK2, 'axes.linewidth': 0.6,
                     'xtick.color': INK2, 'ytick.color': INK2, 'pdf.fonttype': 42, 'mathtext.fontset': 'custom', 'mathtext.it': 'Arial:italic', 'mathtext.rm': 'Arial'})

df = pd.read_csv(R('results', 't13_species_sites.csv'))
st0 = df.set_index('物种')

# ---------- 祖先重建（全部物种，同 #13 主树） ----------
full, miss = T.build(R('data', 'raw', 't13', 'tree', T.TREE_MAIN), list(df.物种))
pnmr = {}   # (site, full_node) → P(裸鼹鼠型)
for s in SITES:
    states, post, q, lnL = T.marginal_all(full, st0[s].to_dict())
    i = states.index(NMR[s])
    for n, v in post.items():
        pnmr[(s, n)] = v[i]

# ---------- 图 1 ----------
fig_sp = [sp for sp in df.物种 if st0.loc[sp, '科'] in HYST_FAM] + OUTGROUPS
tree, miss1 = T.build(R('data', 'raw', 't13', 'tree', T.TREE_MAIN), fig_sp)
assert not miss1, miss1

def ladderize(n):
    for c in n.ch:
        ladderize(c)
    n.ch.sort(key=lambda c: len(T.tips_of(c)))
ladderize(tree)
tips = T.tips_of(tree)
ntip = len(tips)
y = {}
for i, sp in enumerate(tips):
    y[sp] = ntip - 1 - i
H = max(T.depth(n) for n in T.postorder(tree) if not n.ch)
X = {n: T.depth(n) - H for n in T.postorder(tree)}     # 以「距今百万年」为横轴（负数）
for n in T.postorder(tree):
    if n.ch:
        y[n] = np.mean([y[c] if c.ch else y[c.name] for c in n.ch])
    else:
        y[n] = y[n.name]

# 每条图上枝（child c）对应全树里的节点 = c 的叶集合在全树上的 MRCA；父节点同理
fmap = {n: T.mrca(full, T.tips_of(n)) for n in T.postorder(tree)}
gains = []
for n in T.postorder(tree):
    if n.up is None:
        continue
    for s in SITES:
        pc, pp = pnmr[(s, fmap[n])], pnmr[(s, fmap[n.up])]
        if pc >= 0.5 and pp < 0.5:
            gains.append(dict(site=f'{HUM[s]}{s}{NMR[s]}', branch_to=', '.join(sorted(T.tips_of(n))) if len(T.tips_of(n)) <= 4 else f'{len(T.tips_of(n))} spp. clade',
                              P_parent=round(pp, 3), P_child=round(pc, 3), age_from=round(-X[n.up], 1), age_to=round(-X[n], 1), node=n))
gdf = pd.DataFrame([{k: v for k, v in g.items() if k != 'node'} for g in gains])
gdf.to_csv(os.path.join(OUT, 'Figure1_gain_branches.csv'), index=False)
print(gdf.to_string())

fig = plt.figure(figsize=(7.2, 0.155 * ntip + 1.3))
gs = fig.add_gridspec(1, 3, width_ratios=[5.2, 1.25, 1.35], wspace=0.04, left=0.02, right=0.98, top=1 - 0.75 / (0.155 * ntip + 1.3),
                      bottom=0.55 / (0.155 * ntip + 1.3))
axT, axS, axL = (fig.add_subplot(gs[0, i]) for i in range(3))

hyst = T.mrca(tree, ['Cavia porcellus', 'Heterocephalus glaber', 'Hystrix cristata'])
hyst_tips = set(T.tips_of(hyst))
# 4/4 物种行加底色
for sp in tips:
    if st0.loc[sp, '命中'] == 4:
        for ax in (axT, axS, axL):
            ax.axhspan(y[sp] - 0.5, y[sp] + 0.5, color='#fbe3d8', lw=0, zorder=0)
for n in T.postorder(tree):
    if n.up is not None:
        col = INK if (set(T.tips_of(n)) <= hyst_tips) else C_BG
        axT.plot([X[n.up], X[n]], [y[n], y[n]], color=col, lw=0.9, solid_capstyle='butt', zorder=2)
    if n.ch:
        ys = [y[c] for c in n.ch]
        col = INK if (set(T.tips_of(n)) <= hyst_tips) else C_BG
        axT.plot([X[n], X[n]], [min(ys), max(ys)], color=col, lw=0.9, zorder=2)
# 裸鼹鼠型出现的枝：按枝合并位点，画橙色粗段 + 标签
by_branch = {}
for g in gains:
    by_branch.setdefault(g['node'], []).append(g['site'][-4:])
for n, ss in by_branch.items():
    x0, x1 = X[n.up], X[n]
    axT.plot([x0, x1], [y[n], y[n]], color=C_NMR, lw=3.2, solid_capstyle='butt', zorder=3)
    lab = '\n'.join(s[:3] + s[3] for s in sorted(ss))
    axT.text((x0 + x1) / 2, y[n] + 0.38, lab, color=INK, fontsize=6.3, ha='center', va='bottom', zorder=5,
             bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.85))
for sp in tips:
    it = r'$\it{' + sp.replace(' ', r'\ ') + '}$' + (f' ({COMMON[sp]})' if sp in COMMON else '')
    axT.text(X[[n for n in T.postorder(tree) if n.name == sp][0]] + 0.8, y[sp], it, va='center', fontsize=6.6,
             color=INK if sp in hyst_tips else INK2)
# 括号：豪猪亚目
yl = [y[s] for s in hyst_tips]
axT.set_xlim(-H * 1.02, 57)
axT.set_ylim(-0.7, ntip - 0.3)
axT.set_yticks([])
for sp_ in ('top', 'right', 'left'):
    axT.spines[sp_].set_visible(False)
axT.set_xticks([-80, -60, -40, -20, 0])
axT.set_xticklabels(['80', '60', '40', '20', '0'])
axT.set_xlabel('Million years ago (Upham et al. 2019)', color=INK2)
axT.plot([55.0, 55.0], [min(yl) - 0.3, max(yl) + 0.3], color=C_HYS, lw=1.6, clip_on=False)
axT.text(55.8, np.mean(yl), 'Hystricognathi', rotation=90, va='center', ha='left', color=C_HYS, fontsize=7.5, clip_on=False)

# 位点色块
for j, s in enumerate(SITES):
    for sp in tips:
        a = st0.loc[sp, s]
        if a not in T.AA:
            axS.add_patch(Rectangle((j + 0.06, y[sp] - 0.42), 0.88, 0.84, fc='white', ec=C_OTH, hatch='////', lw=0.4))
            continue
        c = C_NMR if a == NMR[s] else (C_HUM if a == HUM[s] else C_OTH)
        axS.add_patch(Rectangle((j + 0.06, y[sp] - 0.42), 0.88, 0.84, fc=c, ec='white', lw=0.6))
        axS.text(j + 0.5, y[sp], a, ha='center', va='center', fontsize=6, color='white' if c != C_OTH else INK, fontweight='bold')
axS.set_xlim(0, 4)
axS.set_ylim(-0.7, ntip - 0.3)
axS.set_xticks([0.5, 1.5, 2.5, 3.5])
axS.set_xticklabels([f'{HUM[s]}{s}' for s in SITES], fontsize=6.5, rotation=0)
axS.xaxis.tick_top()
axS.tick_params(length=0)
axS.set_yticks([])
for sp_ in axS.spines.values():
    sp_.set_visible(False)
axS.set_xlabel('Human cGAS position', color=INK2)

# LQ 条
for sp in tips:
    lq = st0.loc[sp, 'LQ']
    if lq == lq:
        axL.barh(y[sp], lq, height=0.62, color=C_HYS if sp in hyst_tips else C_BG, lw=0)
        axL.text(lq + 0.06, y[sp], f'{lq:.2f}', va='center', fontsize=5.6, color=INK2, zorder=4,
                 bbox=dict(boxstyle='square,pad=0.05', fc='white', ec='none', alpha=0.9))
    else:
        axL.text(0.05, y[sp], 'n/a', va='center', fontsize=5.6, color=INK2)
axL.axvline(1, color=INK2, lw=0.6, ls=(0, (2, 2)), zorder=1)
axL.set_xlim(0, 4.6)
axL.set_ylim(-0.7, ntip - 0.3)
axL.set_yticks([])
axL.set_xticks([0, 1, 2, 3, 4])
for sp_ in ('top', 'right', 'left'):
    axL.spines[sp_].set_visible(False)
axL.set_xlabel('Longevity quotient (LQ)', color=INK2)

# 图例
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
leg = [Patch(fc=C_NMR, label='Naked mole-rat state (S/E/Y/T)'), Patch(fc=C_HUM, label='Human state (D/K/L/K)'),
       Patch(fc=C_OTH, label='Other residue'), Patch(fc='white', ec=C_OTH, hatch='////', label='Not resolved'),
       Line2D([], [], color=C_NMR, lw=3.2, label='Branch where the naked mole-rat state is reconstructed to arise\n(ML Mk-ER marginal P: <0.5 → ≥0.5)'),
       Patch(fc='#fbe3d8', label='Species carrying all four naked mole-rat residues')]
fig.legend(handles=leg, loc='lower left', ncol=3, fontsize=6.2, frameon=False, bbox_to_anchor=(0.02, -0.06))
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT, f'Figure1.{ext}'), dpi=300, bbox_inches='tight')
plt.close(fig)

# ---------- 图 2 ----------
d2 = df[df.LQ.notna() & df[SITES].isin(list(T.AA)).all(axis=1)].copy()
d2['hyst'] = d2.科.isin(HYST_FAM)
rng = np.random.default_rng(1)
d2['xj'] = d2.命中 + rng.uniform(-0.22, 0.22, len(d2))
fig, ax = plt.subplots(figsize=(4.6, 3.6))
o = d2[~d2.hyst]
ax.scatter(o.xj, o.LQ, s=16, color=C_BG, alpha=0.55, lw=0, label=f'Other mammals (n = {len(o)})', zorder=2)
h = d2[d2.hyst]
ax.scatter(h.xj, h.LQ, s=26, color=C_HYS, edgecolor='white', lw=0.8, label=f'Hystricognathi (n = {len(h)})', zorder=3)
ax.set_yscale('log')
ax.set_yticks([0.25, 0.5, 1, 2, 4])
ax.set_yticklabels(['0.25', '0.5', '1', '2', '4'])
from matplotlib.ticker import NullFormatter
ax.yaxis.set_minor_formatter(NullFormatter())
ax.axhline(1, color=INK2, lw=0.6, ls=(0, (2, 2)), zorder=1)
LAB = {'Heterocephalus glaber': ('naked mole-rat', (8, 2)), 'Fukomys anselli': ("Ansell's mole-rat", (8, 0)),
       'Fukomys damarensis': ('Damaraland mole-rat', (8, -2)), 'Cavia porcellus': ('guinea pig', (8, 0)),
       'Cavia tschudii': ('montane guinea pig', (8, -2)), 'Homo sapiens': ('human', (6, 0)),
       'Thryonomys swinderianus': ('greater cane rat', (8, -2))}
for sp, (t, off) in LAB.items():
    r = d2[d2.物种 == sp]
    if len(r):
        ax.annotate(t, (r.xj.iloc[0], r.LQ.iloc[0]), xytext=off, textcoords='offset points', fontsize=6.5, color=INK, va='center')
for k in range(5):
    n = (d2.命中 == k).sum()
    ax.text(k, 0.17, f'n={n}', ha='center', fontsize=6, color=INK2)
ax.set_ylim(0.15, 9)
ax.set_xlim(-0.5, 4.9)
ax.set_xticks(range(5))
ax.set_xlabel('Number of naked mole-rat residues at cGAS 431/479/495/498 (human numbering)', color=INK2, fontsize=7)
ax.set_ylabel('Longevity quotient (log scale)', color=INK2)
for sp_ in ('top', 'right'):
    ax.spines[sp_].set_visible(False)
ax.legend(frameon=False, fontsize=6.5, loc='upper left')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT, f'FigureS1.{ext}'), dpi=300, bbox_inches='tight')
print('fig2 n =', len(d2), 'hyst', len(h), d2.groupby('命中').size().to_dict())
