"""任务 #35 C：调控性（非催化）阳性对照——有实验验证的 DNA 结合 site B / site C 残基和核小体界面残基。
数据全部来自 Table S5（FEL 全树，BH 在 518 个有信息位点内，q ≤ 0.05 且 β < α 判纯化），不重跑任何模型。
用法: python 35_t35_regctrl.py <proj_root>
"""
import sys, os, itertools
import numpy as np, pandas as pd

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
t = pd.read_csv(R('manuscript', 'supplementary', 'supp_table_S5_selection.csv')).set_index('human_position')
SETS = {
    '四个位点 (Chen et al.)': [431, 479, 495, 498],
    '催化三联体': [225, 227, 319],
    '锌指配体': [390, 396, 397, 404],
    'DNA site B (Zhou 2018)': [347, 350],
    'DNA site C (Xie 2019)': [275, 279, 282, 285, 427, 428, 432],
    '核小体界面 (Cao 2020)': [236, 255],
}
REG = SETS['DNA site B (Zhou 2018)'] + SETS['DNA site C (Xie 2019)'] + SETS['核小体界面 (Cao 2020)']
lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

P('FEL 全树（275 种），判据同主结果：BH q ≤ 0.05 且 β < α\n')
P(f'{"组":<24} {"位点":>5} {"残基":>3} {"Fitch 变化":>9} {"ω":>7} {"p":>8} {"q":>8}  判定')
rows = []
for g, sites in SETS.items():
    for s in sites:
        r = t.loc[s]
        om = float(r.FELall_omega)
        rows.append(dict(group=g, pos=s, aa=r.human_residue, fitch=int(r.n_changes_all_mammals_parsimony),
                         omega=om, p=r.FELall_p, q=r.FELall_q, verdict=r.FELall_verdict))
        P(f'{g:<24} {s:5d} {r.human_residue:>3} {rows[-1]["fitch"]:9d} {om:7.3f} {r.FELall_p:8.4f} {r.FELall_q:8.4f}  {r.FELall_verdict}')
df = pd.DataFrame(rows)
df.to_csv(R('results', 't35', 'regctrl.csv'), index=False, encoding='utf-8')

P('\n组汇总：')
for g in list(SETS) + ['调控性对照合计（site B + site C + 核小体，11 个）']:
    d = df[df.pos.isin(REG)] if g.startswith('调控性') else df[df.group == g]
    P(f'  {g:<40} n={len(d):2d}  ω 均值 {d.omega.mean():.2f}  中位数 {d.omega.median():.2f}  '
      f'显著纯化 {(d.verdict == "purifying").sum()}/{len(d)}  ω ≥ 0.75 的 {(d.omega >= 0.75).sum()}/{len(d)}  '
      f'Fitch 变化中位数 {d.fitch.median():.0f}')

reg = df[df.pos.isin(REG)]
focal = df[df.group.str.startswith('四个')]
# 置换：从「11 个调控残基 + 4 个位点」15 个里随机抽 4 个，均值 ω ≥ 观测值的比例（精确枚举）
pool = pd.concat([reg, focal]).omega.values
obs = focal.omega.mean()
combs = list(itertools.combinations(range(len(pool)), 4))
pe = np.mean([pool[list(c)].mean() >= obs - 1e-12 for c in combs])
# 从 11 个调控残基里随机挑 4 个、全部不显著的概率
ns = (reg.verdict != 'purifying').sum()
from math import comb
p_all_ns = comb(int(ns), 4) / comb(len(reg), 4)
from scipy.stats import mannwhitneyu
u = mannwhitneyu(focal.omega, reg.omega, alternative='greater')
P(f'\n四个位点 vs 11 个调控残基：')
P(f'  ω 均值 {obs:.2f} vs {reg.omega.mean():.2f}；中位数 {focal.omega.median():.2f} vs {reg.omega.median():.2f}')
P(f'  精确置换（15 选 4，共 {len(combs)} 种）：均值 ω ≥ {obs:.2f} 的比例 p = {pe:.3f}')
P(f'  Mann–Whitney 单侧（四位点 ω 更高）p = {u.pvalue:.3f}')
P(f'  从 11 个调控残基里随机挑 4 个、全都不显著的概率 = C({ns},4)/C(11,4) = {p_all_ns:.3f}')
P(f'  调控残基里 ω ≥ 1 的：{", ".join(f"{a}{p} (ω={o:.2f}, Fitch {f})" for a, p, o, f in reg[reg.omega >= 1][["aa", "pos", "omega", "fitch"]].values)}')
open(R('results', 't35', 'regctrl.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
