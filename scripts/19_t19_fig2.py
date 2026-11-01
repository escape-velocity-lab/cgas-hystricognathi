"""19b（任务 #19）: 稿件 v3 的 Figure 2「三联证据图」（英文标注，300 dpi，PDF + PNG）。
(a) FEL 全树（275 种）每位点 ω：全部有信息位点的分布 vs 催化三联体 / 锌指配位（阳性对照）vs 4 个位点；FDR（BH）q 值标显著性
(b) 6CT9 人 cGAS–DNA 结构：Cα 二维投影 + 4 个位点与 4 个结构域内谱系候选到 DNA（A 位点）/ 催化三联体的最短距离
(c) cGAS 结构域线性图：4 个位点（上）与 6 个非洲鼹鼠科谱系特异候选（下，A/B/C 可信度）
    任务 #33：R62K 剔除（第四种非洲鼹鼠科 Bathyergus suillus GCA_988425685.1 在 62 位是祖先态 R，见任务 #30）
输入: results/t17_hyphy/fel_all.json、data/raw/t17/pdb/6CT9.cif + 6CTA.cif、results/t15_lineage_candidates.csv
输出: manuscript/figures/Figure2.{png,pdf}、manuscript/figures/fig2_numbers.txt（图上每个数的出处）
回归检查: 4 个位点到 DNA / 催化残基的距离必须与 results/t17_structure_sites.csv 一致（同一算法）
用法: python 19_t19_fig2.py <proj_root>
"""
import sys, os, json
import numpy as np
import pandas as pd
import gemmi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
OUT = R('manuscript', 'figures')
os.makedirs(OUT, exist_ok=True)
log = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); log.append(s)

C_NMR, C_CAND, C_CTRL, C_ZN, C_BG, C_DNA = '#eb6834', '#2a78d6', '#1f8a4c', '#7a5cc2', '#c9c8c3', '#e3b23c'
INK, INK2 = '#0b0b0b', '#52514e'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.edgecolor': INK2, 'axes.linewidth': 0.6,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.labelcolor': INK, 'pdf.fonttype': 42})

SITES = {431: 'D431S', 479: 'K479E', 495: 'L495Y', 498: 'K498T'}
CAT = [225, 227, 319]
ZNL = [390, 396, 397, 404]

def fdr_bh(p):
    p = np.asarray(p, float); o = np.argsort(p); n = len(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

# ---------- (a) FEL 全树 ----------
d = json.load(open(R('results', 't17_hyphy', 'fel_all.json')))
hdr = [h[0] for h in d['MLE']['headers']]
t = pd.DataFrame(d['MLE']['content']['0'], columns=hdr[:len(d['MLE']['content']['0'][0])])
t.insert(0, 'pos', np.arange(1, len(t) + 1))
t['p'] = t['p-value']; t['q'] = fdr_bh(t.p.values)
inf = t[t['Total branch length'] > 0].copy()
inf['omega'] = np.where(inf.alpha > 0, inf.beta / inf.alpha.where(inf.alpha > 0), np.nan)
inf.loc[(inf.alpha == 0) & (inf.beta == 0), 'omega'] = np.nan
FLOOR, CEIL = 1e-3, 1e2
inf['om_plot'] = inf.omega.clip(FLOOR, CEIL)
om = inf.omega.dropna()
P(f'(a) FEL 全树：有信息位点 {len(inf)}；ω 中位数 {om.median():.2f}；FDR q≤0.05 且 β<α（纯化）{int(((inf.q <= 0.05) & (inf.beta < inf.alpha)).sum())}；'
  f'q≤0.1 纯化 {int(((inf.q <= 0.1) & (inf.beta < inf.alpha)).sum())}')
groups = [('Catalytic\ntriad', CAT, C_CTRL), ('Zinc\nthumb', ZNL, C_ZN),
          ('NMR\nfour sites', list(SITES), C_NMR)]
gstat = {}
for name, ps, _ in groups:
    g = inf[inf.pos.isin(ps)]
    nsig = int(((g.q <= 0.05) & (g.beta < g.alpha)).sum())
    gstat[name] = (g.omega.mean(), nsig, len(g))
    P(f'    {name.replace(chr(10), " ")}：' + '，'.join(f'{r.pos} ω={r.omega:.3f} p={r.p:.2g} q={r.q:.2g}' for r in g.itertuples())
      + f' → 平均 ω {g.omega.mean():.2f}，q≤0.05 纯化 {nsig}/{len(g)}')
for s in SITES:
    r = inf[inf.pos == s].iloc[0]
    P(f'    {s} ω 分位 = {np.mean(om < r.omega) * 100:.0f}%')

# ---------- (b) 结构 ----------
PDB = R('data', 'raw', 't17', 'pdb')
def load(code):
    st = gemmi.read_structure(os.path.join(PDB, f'{code}.cif')); st.setup_entities(); st.remove_hydrogens(); return st
is_dna = lambda r: r.name in ('DA', 'DC', 'DG', 'DT')
st = load('6CT9'); stA = load('6CTA')
pc = next([r for r in ch if r.het_flag == 'A'] for ch in st[0] if any(r.name == 'LYS' for r in ch))
res = {r.seqid.num: r for r in pc}
dna = np.array([a.pos.tolist() for ch in st[0] for r in ch if is_dna(r) for a in r])
dnaP = np.array([a.pos.tolist() for ch in st[0] for r in ch if is_dna(r) for a in r if a.name == 'P'])
cat_xyz = np.array([a.pos.tolist() for k in CAT for a in res[k] if a.name not in ('N', 'C', 'O', 'CA')])   # 同 #17：催化残基侧链
zn = np.array([a.pos.tolist() for ch in st[0] for r in ch if r.name == 'ZN' for a in r])
atp = np.array([a.pos.tolist() for ch in stA[0] for r in ch if r.name == 'ATP' for a in r])
def mind(k, B):
    A = np.array([a.pos.tolist() for a in res[k] if a.name not in ('N', 'C', 'O')])   # 同 #17：Cα + 侧链（Gly 只有 Cα）
    return float(np.sqrt(((A[:, None] - B[None]) ** 2).sum(-1)).min())

cand = pd.read_csv(R('results', 't15_lineage_candidates.csv'), encoding='utf-8-sig')
GRADE = {201: 'A', 326: 'A', 258: 'B', 388: 'B', 48: 'C', 19: 'C'}
EXCLUDED = {62: 'R62K：B. suillus 为祖先态 R，非非洲鼹鼠科共有（任务 #30/#33）'}
cand = cand[~cand.人位置.isin(EXCLUDED)].copy()
assert len(cand) == 6, cand.人位置.tolist()
P(f'(c) 谱系候选 {len(cand)} 个：{sorted(int(x) for x in cand.人位置)}；剔除：' + '；'.join(f'{k} {v}' for k, v in EXCLUDED.items()))
cand['anc'] = cand['豪猪亚目祖先ML'].str[0]
cand['lab'] = [f'{a}{int(p)}{x}' for a, p, x in zip(cand.anc, cand.人位置, cand.非洲鼹鼠科状态)]
cand['grade'] = cand.人位置.map(GRADE)
cand = cand.sort_values('人位置')

rows = []
for k, lab in SITES.items():
    rows.append((lab, k, 'site', mind(k, dna), mind(k, cat_xyz)))
for r in cand.itertuples():
    if r.人位置 in res:
        rows.append((f'{r.lab} ({r.grade})', r.人位置, 'cand', mind(r.人位置, dna), mind(r.人位置, cat_xyz)))
dist = pd.DataFrame(rows, columns=['label', 'pos', 'kind', 'dDNA', 'dCAT'])
ref = pd.read_csv(R('results', 't17_structure_sites.csv'), encoding='utf-8-sig')
for (k, _), r in zip(SITES.items(), ref.itertuples()):
    mine = dist[dist.pos == k].iloc[0]
    assert abs(mine.dDNA - r[4]) < 0.06 and abs(mine.dCAT - r[11]) < 0.06, (k, mine.dDNA, r[4], mine.dCAT, r[11])
P('(b) 回归检查通过：4 个位点距离与 t17_structure_sites.csv 一致')
for r in dist.itertuples():
    P(f'    {r.label:<14} 到 DNA(A 位点) {r.dDNA:5.1f} Å   到催化三联体 {r.dCAT:5.1f} Å')
P(f'    结构内候选：{[int(p) for p in dist[dist.kind == "cand"].pos]}；不在 6CT9（157–522）里的：'
  f'{[int(p) for p in cand.人位置 if p not in res]}')

# ---------- 画图 ----------
fig = plt.figure(figsize=(7.2, 6.1))
gs = fig.add_gridspec(2, 3, height_ratios=[3.1, 1.45], width_ratios=[1.0, 1.15, 0.95], hspace=0.36, wspace=0.42,
                      left=0.09, right=0.985, top=0.95, bottom=0.07)
axA = fig.add_subplot(gs[0, 0]); axS = fig.add_subplot(gs[0, 1]); axD = fig.add_subplot(gs[0, 2]); axC = fig.add_subplot(gs[1, :])
def tag(ax, s, x=-0.02):
    ax.text(x, 1.04, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='right', color=INK)

# (a)
ly = np.log10(inf.om_plot.dropna())
vp = axA.violinplot([ly], positions=[0], widths=0.8, showextrema=False)
for b in vp['bodies']:
    b.set_facecolor(C_BG); b.set_edgecolor('none'); b.set_alpha(0.9)
axA.plot([-0.25, 0.25], [np.log10(om.median())] * 2, color=INK2, lw=1)
rng = np.random.default_rng(1)
for i, (name, ps, col) in enumerate(groups, start=1):
    g = inf[inf.pos.isin(ps)]
    xs = i + np.linspace(-0.2, 0.2, len(g))
    for x, r in zip(xs, g.itertuples()):
        sig = (r.q <= 0.05) and (r.beta < r.alpha)
        axA.scatter(x, np.log10(r.om_plot), s=24, marker='o', color=col if sig else 'white', edgecolor=col, lw=1.1, zorder=3)
        if i == 3:
            right = r.pos in (495, 498)
            yo = {431: 0.3, 479: 0.75, 495: 0.25, 498: -0.3}[r.pos]
            axA.annotate(str(r.pos), (x, np.log10(r.om_plot)), xytext=(x + (0.3 if right else -0.3), np.log10(r.om_plot) + yo),
                         fontsize=5.6, color=col, va='center', ha='left' if right else 'right',
                         arrowprops=dict(arrowstyle='-', color=col, lw=0.4, shrinkA=0, shrinkB=2))
    m, ns, n = gstat[name]
    axA.text(i, -4.55, f'mean ω\n{m:.2f}\n{ns}/{n} sig.', ha='center', va='bottom', fontsize=5.4, color=col)
axA.axhline(0, color=INK2, lw=0.6, ls=':')
axA.text(-0.55, 0.06, 'ω = 1', fontsize=5.4, color=INK2, ha='left', va='bottom')
axA.set_xticks([0, 1, 2, 3])
axA.set_xticklabels(['All\nsites'] + [g[0] for g in groups], fontsize=5.8)
axA.set_yticks([-3, -2, -1, 0, 1, 2]); axA.set_yticklabels(['≤0.001', '0.01', '0.1', '1', '10', '≥100'])
axA.set_ylim(-4.6, 2.6); axA.set_xlim(-0.6, 3.6)
axA.set_ylabel('Site ω = dN/dS (FEL, 275 mammals)')
for s_ in ('top', 'right'):
    axA.spines[s_].set_visible(False)
tag(axA, 'a', -0.22)

# (b) 结构投影
ks = sorted(res); C = np.array([res[k]['CA'][0].pos.tolist() for k in ks])
cen = C.mean(0)
_, _, Vt = np.linalg.svd(np.vstack([C, dnaP]) - cen)
proj = lambda X: (np.atleast_2d(X) - cen) @ Vt[:2].T
p = proj(C)
for i in range(len(ks) - 1):
    if ks[i + 1] - ks[i] == 1:
        axS.plot(p[i:i + 2, 0], p[i:i + 2, 1], color='#a9a8a3' if ks[i] < 331 else '#7d7c77', lw=0.55, zorder=1)
d2 = proj(dnaP)
axS.scatter(d2[:, 0], d2[:, 1], s=7, color=C_DNA, edgecolor='none', zorder=2)
c2 = proj(np.array([res[k]['CA'][0].pos.tolist() for k in CAT]))
axS.scatter(c2[:, 0], c2[:, 1], s=22, marker='^', color=C_CTRL, edgecolor='white', lw=0.4, zorder=3)
if len(atp):
    axS.scatter(*proj(atp).mean(0), s=55, marker='*', color=C_CTRL, edgecolor='white', lw=0.4, zorder=3)
axS.scatter(*proj(zn)[0], s=22, marker='D', color=C_ZN, edgecolor='white', lw=0.4, zorder=3)
OFF = {431: (6, 4), 479: (-6, -9), 495: (6, -2), 498: (5, -9), 201: (5, 3), 258: (5, 3), 326: (5, -6), 388: (-5, 4)}
for r in dist.itertuples():
    q = proj(res[r.pos]['CA'][0].pos.tolist())[0]
    col = C_NMR if r.kind == 'site' else C_CAND
    axS.scatter(*q, s=30, color=col, marker='o' if r.kind == 'site' else 's', edgecolor='white', lw=0.6, zorder=4)
    dx, dy = OFF.get(r.pos, (5, 3))
    axS.annotate(str(r.pos), q, xytext=(dx, dy), textcoords='offset points', fontsize=6, color=col, fontweight='bold',
                 ha='right' if dx < 0 else 'left', zorder=5)
axS.set_aspect('equal'); axS.set_xticks([]); axS.set_yticks([])
for s_ in axS.spines.values():
    s_.set_visible(False)
axS.set_title('Human cGAS–DNA (PDB 6CT9), Cα trace', fontsize=6.6, color=INK, pad=2)
leg = [Line2D([], [], marker='o', ls='', color=C_NMR, ms=4.5, label='NM-R "longevity" site'),
       Line2D([], [], marker='s', ls='', color=C_CAND, ms=4.5, label='African mole-rat candidate'),
       Line2D([], [], marker='^', ls='', color=C_CTRL, ms=4.5, label='Catalytic E225/D227/D319'),
       Line2D([], [], marker='*', ls='', color=C_CTRL, ms=6, label='ATP (6CTA)'),
       Line2D([], [], marker='D', ls='', color=C_ZN, ms=3.8, label='Zn$^{2+}$ (zinc thumb)'),
       Line2D([], [], marker='o', ls='', color=C_DNA, ms=3, label='DNA phosphate')]
axS.legend(handles=leg, loc='upper center', bbox_to_anchor=(0.5, -0.0), ncol=2, fontsize=5.4, frameon=False,
           handletextpad=0.2, columnspacing=0.6)
tag(axS, 'b', 0.02)

# (b) 距离
dd = dist.iloc[::-1].reset_index(drop=True)
for i, r in dd.iterrows():
    col = C_NMR if r.kind == 'site' else C_CAND
    axD.plot([r.dDNA, r.dCAT], [i, i], color='#d9d8d4', lw=1.2, zorder=1)
    axD.scatter(r.dDNA, i, s=22, marker='o', color=C_DNA, edgecolor=INK2, lw=0.4, zorder=3)
    axD.scatter(r.dCAT, i, s=22, marker='^', color=C_CTRL, edgecolor=INK2, lw=0.4, zorder=3)
axD.set_yticks(range(len(dd)))
axD.set_yticklabels(dd.label, fontsize=6)
for tl, r in zip(axD.get_yticklabels(), dd.itertuples()):
    tl.set_color(C_NMR if r.kind == 'site' else C_CAND)
axD.axvspan(0, 5, color='#f1f0ec', zorder=0)
axD.text(2.5, len(dd) + 0.05, 'contact\n(≤ 5 Å)', fontsize=5.2, color=INK2, ha='center', va='top')
axD.set_xlim(0, 38); axD.set_ylim(-0.6, len(dd) + 0.1)
axD.set_xlabel('Minimum distance, Cα + side chain (Å)')
axD.axhline(len(dd) - len(SITES) - 0.5, color=INK2, lw=0.5, ls='--')
for s_ in ('top', 'right'):
    axD.spines[s_].set_visible(False)
axD.legend(handles=[Line2D([], [], marker='o', ls='', color=C_DNA, mec=INK2, mew=0.4, ms=4.5, label='to DNA (site A)'),
                    Line2D([], [], marker='^', ls='', color=C_CTRL, mec=INK2, mew=0.4, ms=4.5, label='to catalytic triad')],
           loc='lower right', fontsize=5.4, frameon=False, handletextpad=0.2)

# (c) 结构域图
L = 522
DOM = [(1, 160, 'N-terminal (disordered)', '#efeeea'), (161, 330, 'NTase core', '#dcdcd6'), (331, 522, 'C-terminal domain', '#c9c8c2')]
y0, h = 0, 1
for a, b, nm, col in DOM:
    axC.add_patch(Rectangle((a, y0), b - a + 1, h, fc=col, ec=INK2, lw=0.5))
    axC.text((a + b) / 2, y0 + h / 2, nm, ha='center', va='center', fontsize=6, color=INK)
for a, b in ((173, 215), (384, 407)):
    axC.add_patch(Rectangle((a, y0 - 0.22), b - a + 1, 0.18, fc=C_DNA, ec='none'))
axC.text(194, y0 - 0.32, 'DNA-binding', fontsize=5.2, color='#9a7419', ha='center', va='top')
axC.text(395, y0 - 0.32, 'DNA-binding / zinc thumb', fontsize=5.2, color='#9a7419', ha='center', va='top')
for k in CAT:
    axC.plot([k, k], [y0, y0 + h], color=C_CTRL, lw=1.0)
axC.text(272, y0 + h + 0.08, 'E225/D227/D319', fontsize=5.2, color=C_CTRL, ha='center', va='bottom')
# 4 个位点（上）
for i, (k, lab) in enumerate(SITES.items()):
    xt = 432 + i * 24
    axC.plot([k, k, xt], [y0 + h, y0 + h + 0.35, y0 + h + 0.75], color=C_NMR, lw=0.7)
    axC.text(xt, y0 + h + 0.8, lab, fontsize=6, color=C_NMR, ha='center', va='bottom', fontweight='bold')
axC.text(522, y0 + h + 1.42, 'Chen et al. (2025) four residues: hystricognath-wide (Fig. 1)', fontsize=5.6, color=C_NMR, ha='right')
# 6 个候选（下）
GC = {'A': C_CAND, 'B': '#6aa2e3', 'C': '#a9c6ec'}
xs_lab = {19: 25, 48: 70, 201: 201, 258: 258, 326: 326, 388: 388}
for r in cand.itertuples():
    k = int(r.人位置); col = GC[r.grade]
    axC.plot([k, k, xs_lab[k]], [y0, y0 - 0.75, y0 - 1.1], color=col, lw=0.7)
    axC.text(xs_lab[k], y0 - 1.15, f'{r.lab}\n({r.grade})', fontsize=5.8, color=C_CAND if r.grade != 'C' else '#5f8fcf',
             ha='center', va='top', fontweight='bold' if r.grade == 'A' else 'normal')
axC.text(0, y0 - 2.15, 'African mole-rat lineage-specific substitutions (ancestral hystricognath → Bathyergidae; confidence A/B/C)',
         fontsize=5.6, color=C_CAND, ha='left', va='top')
axC.set_xlim(-5, L + 5); axC.set_ylim(-2.6, 2.6)
axC.set_xticks([1, 100, 200, 300, 400, 522]); axC.set_yticks([])
axC.tick_params(axis='x', labelsize=5.8, pad=1)
for s_ in ('top', 'right', 'left'):
    axC.spines[s_].set_visible(False)
axC.spines['bottom'].set_position(('data', -2.6))
axC.set_xlabel('Human cGAS residue (Q8N884)', fontsize=6.2, labelpad=1)
tag(axC, 'c', -0.01)

for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT, f'Figure2.{ext}'), dpi=300)
open(os.path.join(OUT, 'fig2_numbers.txt'), 'w', encoding='utf-8').write('\n'.join(log) + '\n')
P('图 → manuscript/figures/Figure2.{png,pdf}')
