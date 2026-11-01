"""28（任务 #23d-b）: 汇总「剔除 17 个 tblastn 重建物种后、258 个物种全树 FEL」的结果，
与 275 物种全树 FEL（fel_all.json）逐位点对照。

要回答的问题（来自 #21 自审）：117 个密码子（人 406-522）是低质量 tblastn 序列唯一贡献数据的区段，
比对误差会抬高 dN、把 ω 推向 1，**正好是我们阴性结论的方向**。只用有真实 CDS 的 258 个物种重跑，
4 个位点的 ω 和显著性是否不变？阳性对照（催化三联体、锌指）是否照样找回？
用法: python 28_t23_fel_cds_only.py <proj_root>
"""
import sys, os, json
import numpy as np
import pandas as pd

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
H = lambda f: R('results', 't17_hyphy', f)
SITES = [431, 479, 495, 498]
CTRL = {'催化三联体 E225/D227/D319': [225, 227, 319], '锌指配位 H390/C396/C397/C404': [390, 396, 397, 404]}
L = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); L.append(s)

def fdr_bh(p):
    p = np.asarray(p, float); o = np.argsort(p); n = len(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

def load(fn, label):
    d = json.load(open(fn))
    hdr = [h[0] for h in d['MLE']['headers']]
    t = pd.DataFrame(d['MLE']['content']['0'], columns=hdr[:len(d['MLE']['content']['0'][0])])
    t.insert(0, 'pos', np.arange(1, len(t) + 1))
    t['omega'] = np.where(t['alpha'] > 0, t['beta'] / t['alpha'].where(t['alpha'] > 0), np.nan)
    # q 值只在「有信息位点」上做 BH（和 16_t17_hyphy_summary.py 口径一致）
    inf = t['Total branch length'] > 0 if 'Total branch length' in t else pd.Series(True, index=t.index)
    t['q'] = np.nan
    t.loc[inf, 'q'] = fdr_bh(t.loc[inf, 'p-value'].values)
    t['sig'] = (t['q'] <= 0.1) & inf
    t['purifying'] = t['sig'] & (t['beta'] < t['alpha'])
    P('%s: %d 个密码子，%d 个有信息位点，显著纯化 %d（%.0f%%），显著正选择 %d；全局 ω %s' % (
        label, len(t), inf.sum(), t['purifying'].sum(), 100 * t['purifying'].sum() / max(inf.sum(), 1),
        (t['sig'] & ~t['purifying']).sum(), d['fits']['Global MG94xREV']['Rate Distributions']))
    return t

f275 = load(H('fel_all.json'), '275 物种全树 FEL（原分析）')
f258 = load(H('fel_cds_only.json'), '258 物种全树 FEL（剔除 17 个 tblastn 重建物种）')

P('')
P('===== 4 个位点 =====')
P('%-6s %-26s %-26s' % ('位点', '275 物种 ω (q)', '258 物种 ω (q)'))
rows = []
for p in SITES:
    a = f275[f275.pos == p].iloc[0]; b = f258[f258.pos == p].iloc[0]
    P('%-6d %-26s %-26s' % (p, '%.2f (%.2f)%s' % (a.omega, a.q, ' *' if a.sig else ''),
                            '%.2f (%.2f)%s' % (b.omega, b.q, ' *' if b.sig else '')))
    rows.append(dict(site=p, omega_275=round(a.omega, 3), q_275=round(a.q, 4), sig_275=bool(a.sig),
                     omega_258=round(b.omega, 3), q_258=round(b.q, 4), sig_258=bool(b.sig)))
P('')
P('===== 阳性对照 =====')
for name, ps in CTRL.items():
    for p in ps:
        a = f275[f275.pos == p].iloc[0]; b = f258[f258.pos == p].iloc[0]
        P('%-34s %-4d  275: ω=%.2f q=%.4f%s   258: ω=%.2f q=%.4f%s' % (
            name, p, a.omega, a.q, ' *' if a.sig else '', b.omega, b.q, ' *' if b.sig else ''))
        rows.append(dict(site=p, omega_275=round(a.omega, 3), q_275=round(a.q, 4), sig_275=bool(a.sig),
                         omega_258=round(b.omega, 3), q_258=round(b.q, 4), sig_258=bool(b.sig), control=name))
P('')
m = f275[['pos', 'omega', 'q']].merge(f258[['pos', 'omega', 'q']], on='pos', suffixes=('_275', '_258'))
both = m.dropna(subset=['omega_275', 'omega_258'])
both = both[np.isfinite(both.omega_275) & np.isfinite(both.omega_258)]
P('全部位点 ω 的相关性（%d 个两边都有估计的位点）: Pearson r = %.3f, Spearman = %.3f' % (
    len(both), both.omega_275.corr(both.omega_258), both.omega_275.corr(both.omega_258, method='spearman')))
P('两组各自的显著纯化位点数: %d vs %d；交集 %d' % (
    f275.purifying.sum(), f258.purifying.sum(),
    len(set(f275[f275.purifying].pos) & set(f258[f258.purifying].pos))))
P('')
P('结论：4 个位点在剔除低质量序列后是否仍然全部不显著 -> %s' % ('是' if not f258[f258.pos.isin(SITES)].sig.any() else '否（需改正文）'))

pd.DataFrame(rows).to_csv(R('results', 't23_fel_cds_only.csv'), index=False, encoding='utf-8-sig')
open(R('results', 't23_fel_cds_only.txt'), 'w', encoding='utf-8').write('\n'.join(L))
P('写出 results/t23_fel_cds_only.csv / .txt')
