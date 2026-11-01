"""17a-3（任务 #17）: 汇总 HyPhy FEL / MEME / RELAX 结果，4 个位点 vs cGAS 全部位点。

比对列 = 人 cGAS 编号（12_t17_codon_aln.py 已按人坐标对齐），所以 HyPhy 的第 k 个密码子 = 人 k 位。
FEL：α（同义速率）、β（非同义速率）、ω = β/α、p（β ≠ α 的 LRT）。β < α 显著 = 纯化选择（受约束），β > α = 正选择。
MEME：只检验「有没有一部分分支上 ω > 1」（片段性正选择），不检验约束。
RELAX：K < 1 = 前景（豪猪亚目）选择强度放松，K > 1 = 加强。
输出：results/t17_hyphy_sites.csv、t17_hyphy_summary.txt、t17_hyphy_omega.png
用法: python 16_t17_hyphy_summary.py <proj_root>
"""
import sys, os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
H = lambda f: R('results', 't17_hyphy', f)
lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

SITES = [431, 479, 495, 498]
# 对照组：cGAS 已知功能必需的残基。它们是**阳性对照**——如果 FEL 在这个数据集上能把它们检出强约束，
# 那么 4 个位点的「不显著」就不能用「功效不足」解释。
CTRL = {'催化三联体 E225/D227/D319': [225, 227, 319],
        'site A DNA 结合 K173/R176/K407': [173, 176, 407],
        '锌指配位 H390/C396/C397/C404': [390, 396, 397, 404]}
def fdr_bh(p):
    """Benjamini–Hochberg FDR。522 个位点逐个检验，必须做多重检验校正。"""
    p = np.asarray(p, float); o = np.argsort(p); n = len(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

desc = pd.read_csv(R('results', 't17_descriptive.csv'))

def mle(fn):
    """MEME 的表头用 HTML 实体（&alpha; / &beta;<sup>+</sup>），统一成纯文本列名"""
    d = json.load(open(fn))
    ren = {'&alpha;': 'alpha', '&beta;<sup>1</sup>': 'beta1', 'p<sup>1</sup>': 'p1',
           '&beta;<sup>+</sup>': 'beta+', 'p<sup>+</sup>': 'p+', 'FEL &alpha;': 'FEL alpha', 'FEL &beta;': 'FEL beta'}
    hdr = [ren.get(h[0], h[0]) for h in d['MLE']['headers']]
    rows = d['MLE']['content']['0']
    return d, pd.DataFrame(rows, columns=hdr[:len(rows[0])])

def fel(tag, label, off=0):
    fn = H(f'fel_{tag}.json')
    if not os.path.exists(fn) or os.path.getsize(fn) == 0:
        P(f'（{label}：结果缺失）'); return None
    d, t = mle(fn)
    t.insert(0, '人位置', np.arange(1, len(t) + 1) + off)
    t['omega'] = np.where(t['alpha'] > 0, t['beta'] / t['alpha'].where(t['alpha'] > 0), np.nan)
    t['q'] = fdr_bh(t['p-value'].values)
    t['类别'] = np.where(t['p-value'] <= 0.1, np.where(t['beta'] < t['alpha'], '纯化', '正选择'), '不显著')
    inf = t[t['Total branch length'] > 0] if 'Total branch length' in t else t
    P(f'\n===== FEL（{label}）')
    P(f'  全局 MG94×REV ω = {d["fits"]["Global MG94xREV"]["Rate Distributions"]}')
    P(f'  有信息的位点 {len(inf)} 个：纯化 {np.sum(inf.类别 == "纯化")}（{np.mean(inf.类别 == "纯化"):.0%}），'
      f'正选择 {np.sum(inf.类别 == "正选择")}，不显著 {np.sum(inf.类别 == "不显著")}（原始 p ≤ 0.1）；'
      f'**FDR 校正后 q ≤ 0.1 的有 {int((t.q <= 0.1).sum())} 个**')
    om = inf.omega.replace([np.inf], np.nan).dropna()
    P(f'  每位点 ω（β/α，α>0 的位点）：中位数 {om.median():.3f}，四分位 {om.quantile(.25):.3f}–{om.quantile(.75):.3f}')
    # β 不受 α=0 影响，直接比排名
    for s in SITES:
        r = t[t.人位置 == s].iloc[0]
        rk_b = np.mean(inf.beta < r.beta) * 100
        rk_o = np.mean(om < r.omega) * 100 if np.isfinite(r.omega) else np.nan
        P(f'  {s}: α={r.alpha:.3f} β={r.beta:.3f} ω={r.omega:.3f}  p={r["p-value"]:.4f} q={r.q:.3f} → '
          f'{r.类别 if r.q <= 0.1 else "校正后不显著"}；β 高于 {rk_b:.0f}% 的位点，ω 高于 {rk_o:.0f}% 的位点')
    P('  —— 对照组（判断「不显著」是不是功效不足）——')
    tt = t.rename(columns={'p-value': 'p'})
    def grp(name, ps):
        g = tt[tt.人位置.isin(ps)]
        P(f'    {name}：' + '，'.join(f'{int(r.人位置)} ω={r.omega:.2f} p={r.p:.3f} q={r.q:.3f}' for r in g.itertuples())
          + f'  → 平均 ω {g.omega.mean():.2f}，原始 p 显著纯化 {int(((g.p <= 0.1) & (g.beta < g.alpha)).sum())}/{len(g)}，'
            f'FDR 后 {int(((g.q <= 0.1) & (g.beta < g.alpha)).sum())}/{len(g)}')
    grp('【这 4 个位点】', SITES)
    for n, ps in CTRL.items():
        grp(n, ps)
    g = tt[tt.人位置.isin(SITES)]
    P(f'    全部位点：ω 中位数 {tt.omega.median():.2f}，显著纯化占 {np.mean((tt.p <= 0.1) & (tt.beta < tt.alpha)):.0%}')
    return t

fa = fel('all', '全部分支')
fh = fel('hy', '只检验豪猪亚目 45 条分支（其余分支的 β 作为干扰参数单独估计）')
fe = fel('exon', '敏感性分析：只用人 406–522 的 117 列（275 种全部有数据），前景 = 豪猪亚目', off=405)

def meme(tag, label):
    fn = H(f'meme_{tag}.json')
    if not os.path.exists(fn) or os.path.getsize(fn) == 0:
        P(f'\n（MEME {label}：结果缺失）'); return None
    d, t = mle(fn)
    t.insert(0, '人位置', np.arange(1, len(t) + 1))
    P(f'\n===== MEME（{label}）')
    sig = t[t['p-value'] <= 0.1]
    P(f'  p ≤ 0.1 的位点 {len(sig)} 个：{sig.人位置.tolist()}')
    tt = t.rename(columns={'p-value': 'p'})
    tt['q'] = fdr_bh(tt.p.values)
    P(f'  FDR 校正后 q ≤ 0.1 的位点 {int((tt.q <= 0.1).sum())} 个')
    for s in SITES:
        r = tt[tt.人位置 == s].iloc[0]
        P(f'  {s}: α={r.alpha:.3f} β+={r["beta+"]:.3f} p+={r["p+"]:.3f} LRT={r.LRT:.3f} p={r.p:.4f} q={r.q:.3f}'
          + ('  ← FDR 后仍显著' if r.q <= 0.1 else ''))
    for n, ps in CTRL.items():
        g = tt[tt.人位置.isin(ps)]
        if len(g):
            P(f'    对照 {n}：' + '，'.join(f'{int(r.人位置)} p={r.p:.3f} q={r.q:.3f}' for r in g.itertuples()))
    return t

ma, mh = meme('all', '全部分支'), meme('hy', '只检验豪猪亚目分支')

for tag, label in (('relax', '全长 522 密码子'), ('relax_exon', '敏感性分析：117 列子集（275 种全有数据）')):
    if not os.path.exists(H(f'{tag}.json')) or os.path.getsize(H(f'{tag}.json')) == 0:
        P(f'\n（RELAX {label}：结果缺失）'); continue
    d = json.load(open(H(f'{tag}.json')))
    tr = d['test results']
    P(f'\n===== RELAX（{label}；豪猪亚目 45 条分支 = Test，其余 502 条 = Reference）')
    P(f'  K = {tr["relaxation or intensification parameter"]:.3f}，LRT = {tr["LRT"]:.3f}，p = {tr["p-value"]:.4g}'
      '   （K<1 选择放松 / K>1 加强 / K≈1 两组无差别）')
    for m in ('RELAX null', 'RELAX alternative', 'RELAX partitioned descriptive'):
        if m in d['fits']:
            P(f'  {m}: lnL = {d["fits"][m]["Log Likelihood"]:.2f}')

# 汇总表
out = desc[desc.人位置.isin(range(1, 523))][['人位置', '人残基', '豪猪亚目状态', '豪猪亚目变化次数', '全树变化次数']].copy()
for tag, t in (('全树', fa), ('豪猪亚目', fh), ('豪猪亚目117列', fe)):
    if t is not None:
        out = out.merge(t[['人位置', 'alpha', 'beta', 'omega', 'p-value', 'q', '类别']].rename(columns=lambda c: c if c == '人位置' else f'FEL{tag}_{c}'), on='人位置')
for tag, t in (('全树', ma), ('豪猪亚目', mh)):
    if t is not None:
        out = out.merge(t[['人位置', 'p-value']].rename(columns={'p-value': f'MEME{tag}_p'}), on='人位置')
out.to_csv(R('results', 't17_hyphy_sites.csv'), index=False, encoding='utf-8-sig')
P('\n4 个位点汇总：\n' + out[out.人位置.isin(SITES)].T.to_string())
open(R('results', 't17_hyphy_summary.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')

# 图：ω 分布 + 4 个位点 vs 阳性对照（这张图是整个 17a 的核心论据）
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
LO, HI = 1e-3, 1e2
clipv = lambda x: np.log10(np.clip(x, LO, HI))
panels = [(t, l) for t, l in ((fa, '全部 275 种 · 全部分支'), (fh, '只检验豪猪亚目 45 条分支')) if t is not None]
fig, axes = plt.subplots(len(panels), 1, figsize=(10.5, 4.6 * len(panels)), squeeze=False)
for ax, (t, l) in zip(axes[:, 0], panels):
    tt = t.rename(columns={'p-value': 'p'})
    om = clipv(tt.omega.replace([np.inf], np.nan).dropna())
    n, bins, _ = ax.hist(om, bins=45, color='#D4D8DE', label=f'全部 {len(om)} 个位点')
    top = n.max()
    ax.axvline(0, color='#333', ls=':', lw=1.2, zorder=1)
    ax.text(0.02, top * 0.97, 'ω = 1\n中性', fontsize=8, color='#333', va='top')
    # 阳性对照：已知功能必需的残基
    for ps, c, nm in ((CTRL['催化三联体 E225/D227/D319'], '#1F77B4', '催化三联体'),
                      (CTRL['锌指配位 H390/C396/C397/C404'], '#2CA02C', '锌指配位')):
        g = tt[tt.人位置.isin(ps)]
        ax.scatter(clipv(g.omega.values), [top * 0.70] * len(g), s=95, marker='v', color=c,
                   edgecolor='white', lw=1, zorder=5, label=f'{nm}（功能必需，阳性对照）')
    # 目标：4 个位点
    g = tt[tt.人位置.isin(SITES)]
    ax.scatter(clipv(g.omega.values), [top * 0.84] * len(g), s=150, marker='v', color='#D62728',
               edgecolor='white', lw=1.2, zorder=6, label='这 4 个位点')
    # 标签：按 ω 排序后沿 x 轴强制拉开间距，用引线连回真实位置
    gs = g.sort_values('omega')
    xs = clipv(gs.omega.values)
    lab_x, MIN_GAP = [], 0.45
    for x in xs:
        lab_x.append(x if not lab_x else max(x, lab_x[-1] + MIN_GAP))
    shift = (max(lab_x) - max(xs)) / 2 if lab_x else 0
    lo_lim = np.log10(LO)
    if lab_x and lab_x[0] - shift < lo_lim:          # 别把标签推出左边界
        shift = lab_x[0] - lo_lim
    for (r, x, lx) in zip(gs.itertuples(), xs, lab_x):
        star = ' *' if r.q <= 0.1 else (' (p)' if r.p <= 0.1 else '')
        ax.annotate(f'{int(r.人位置)}{star}', (x, top * 0.84), xytext=(lx - shift, top * 0.99),
                    ha='center', fontsize=9.5, color='#D62728', weight='bold',
                    arrowprops=dict(arrowstyle='-', color='#D62728', lw=0.7, shrinkA=1, shrinkB=4))
    med = np.median(tt.omega.dropna())
    ax.axvline(np.log10(med), color='#888', lw=1.2, ls='--', label=f'全部位点中位数 ω = {med:.2f}')
    ax.set_xlim(np.log10(LO) - 0.15, np.log10(HI) + 0.15)
    ax.set_xticks([-3, -2, -1, 0, 1, 2])
    ax.set_xticklabels(['≤0.001', '0.01', '0.1', '1', '10', '≥100'])
    ax.set_xlabel('每个位点的 ω = dN/dS（对数轴）；越靠左越受约束')
    ax.set_ylabel('位点数')
    ax.set_title(f'FEL：{l}　—　4 个位点落在分布右侧，功能必需残基落在最左端', fontsize=11.5)
    ax.legend(fontsize=8.5, frameon=False, loc='center left', bbox_to_anchor=(0.02, 0.62))
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
axes[-1, 0].text(0, -0.22, '* = FDR 校正后 q ≤ 0.1；(p) = 仅原始 p ≤ 0.1、校正后不显著。ω<1 纯化选择（受约束），ω≈1 中性，ω>1 正选择。'
                 'ω=0 的位点画在最左端（≤0.001）。\n'
                 '阳性对照的意义：同一套数据里功能必需残基被干净地检出强约束，'
                 '说明方法有效——所以 4 个位点的「不显著」不是功效不足。',
                 transform=axes[-1, 0].transAxes, fontsize=8.5, color='#555', va='top')
plt.tight_layout(rect=(0, 0.06, 1, 1))
plt.savefig(R('results', 't17_hyphy_omega.png'), dpi=165)
P('图 → results/t17_hyphy_omega.png')
