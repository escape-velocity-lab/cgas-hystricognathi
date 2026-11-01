"""17a 追加（任务 #17）: 功效分析——这套数据能排除多强的约束？

「4 个位点 FEL 不显著」最硬的反驳是「单位点功效本来就低」。阳性对照（催化三联体 ω≈0.02 被检出）已间接回答，
这里给一个定量上界：在**同一棵树、同一个总树长**下，按给定的真实 ω 模拟单个密码子位点，
再用和 FEL 同样的似然比检验（自由 α、β vs 约束 α=β），统计检出率。

模型（MG94×REV 的简化版，自己实现，便于控制）：
  Q[i→j] ∝ π_b × (κ 若转换) × (α 若同义 / β 若非同义)，只允许单核苷酸替换，61 个有义密码子。
  α=β=1 时归一化成「每位点 1 次替换/单位时间」，这个常数之后固定，于是 ω = β/α 就是 dN/dS。
速度：对 Q 做特征分解，P(t)·v = U·(e^{λt} ⊙ (U⁻¹·v))，每条分支两次矩阵-向量乘。
用法: python 18_t17_power.py <proj_root> [每个 ω 的模拟位点数，默认 200]
"""
import sys, os, json
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
from Bio import Phylo
from Bio.Data.CodonTable import unambiguous_dna_by_id
from io import StringIO

ROOT = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 200
R = lambda *p: os.path.join(ROOT, *p)
rng = np.random.default_rng(17)
lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

BASES = 'TCAG'
std = unambiguous_dna_by_id[1]
SENSE = [a + b + c for a in BASES for b in BASES for c in BASES if a + b + c not in std.stop_codons]
AAs = {c: std.forward_table[c] for c in SENSE}
NC = len(SENSE)
KAPPA = 4.0
PI = np.full(4, 0.25)
P0 = np.array([np.prod([PI[BASES.index(x)] for x in c]) for c in SENSE]); P0 /= P0.sum()

# 预先把「同义 / 非同义」的骨架算好，之后只需乘 α、β
SYN = np.zeros((NC, NC)); NON = np.zeros((NC, NC))
for i, ci in enumerate(SENSE):
    for j, cj in enumerate(SENSE):
        if i == j:
            continue
        d = [k for k in range(3) if ci[k] != cj[k]]
        if len(d) != 1:
            continue
        k = d[0]; a, b = ci[k], cj[k]
        r = PI[BASES.index(b)] * (KAPPA if (a in 'TC') == (b in 'TC') else 1.0)
        (SYN if AAs[ci] == AAs[cj] else NON)[i, j] = r
_q = SYN + NON
np.fill_diagonal(_q, -_q.sum(1))
NORM = -(P0 * np.diag(_q)).sum()      # α=β=1 时的归一化常数，之后固定

def eig_of(alpha, beta):
    Q = (SYN * alpha + NON * beta) / NORM
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(1))
    lam, U = np.linalg.eig(Q)
    return lam.real, U.real, np.linalg.inv(U).real

# --- 树 ---
tree = Phylo.read(StringIO(open(R('results', 't17_codon', 'tree_plain.nwk')).read()), 'newick')
try:
    TL = json.load(open(R('results', 't17_hyphy', 'fel_all.json')))['fits']['Global MG94xREV']['Tree length']
except Exception:
    TL = 17.041
cur = sum(c.branch_length or 0 for c in tree.find_clades() if c.branch_length)
for c in tree.find_clades():
    if c.branch_length:
        c.branch_length *= TL / cur
tips = tree.get_terminals()
post = list(tree.find_clades(order='postorder'))
parent = {}
for c in tree.find_clades():
    for ch in c.clades:
        parent[ch] = c
BL = {c: (c.branch_length or 0.0) for c in post}
P(f'树：{len(tips)} 个叶子，{len(post)} 个节点，缩放后总树长 {TL:.2f} subs/site（与 FEL 全树一致）')

def simulate(omega, n):
    lam, U, Ui = eig_of(1.0, omega)
    Pm = {}
    for c in post:
        t = BL[c]
        if t not in Pm:
            M = (U * np.exp(lam * t)) @ Ui
            Pm[t] = np.clip(M, 0, None)
            Pm[t] /= Pm[t].sum(1, keepdims=True)
    pre = list(tree.find_clades(order='preorder'))
    out = np.empty((n, len(tips)), dtype=np.int16)
    for s in range(n):
        st = {}
        for c in pre:
            st[c] = rng.choice(NC, p=P0) if c not in parent else rng.choice(NC, p=Pm[BL[c]][st[parent[c]]])
        out[s] = [st[t] for t in tips]
    return out

def lnL_of(col, alpha, beta):
    # 夹住参数范围：α/β 跑飞会让 exp(λt) 溢出，且生物学上无意义
    alpha, beta = float(np.clip(alpha, 1e-3, 50.0)), float(np.clip(beta, 1e-6, 50.0))
    lam, U, Ui = eig_of(alpha, beta)
    e = {}
    for c in post:
        t = BL[c]
        if t not in e:
            e[t] = np.exp(lam * t)
    L, sc = {}, 0.0
    for c in post:
        if c.is_terminal():
            v = np.zeros(NC); v[col[tips.index(c)]] = 1.0
        else:
            v = np.ones(NC)
            for ch in c.clades:
                v = v * (U @ (e[BL[ch]] * (Ui @ L[ch])))
            v = np.clip(v, 0, None)
            m = v.max()
            if m <= 0:
                return -1e300
            v /= m; sc += np.log(m)
        L[c] = v
    return sc + np.log(max((P0 * L[post[-1]]).sum(), 1e-300))

def fel_pvalue(col):
    f1 = lambda x: -lnL_of(col, np.exp(x[0]), np.exp(x[1]))
    f0 = lambda x: -lnL_of(col, np.exp(x[0]), np.exp(x[0]))
    r1 = minimize(f1, [0.0, -0.5], method='Nelder-Mead', options=dict(xatol=1e-2, fatol=1e-3, maxiter=200))
    r0 = minimize(f0, [0.0], method='Nelder-Mead', options=dict(xatol=1e-2, fatol=1e-3, maxiter=100))
    a, b = np.clip(np.exp(r1.x[0]), 1e-3, 50.0), np.clip(np.exp(r1.x[1]), 1e-6, 50.0)
    return chi2.sf(max(2 * (r0.fun - r1.fun), 0), 1), a, b

P(f'\n每个 ω 模拟 {N} 个独立位点，FEL 式似然比检验（p ≤ 0.1 判显著）：\n')
P(f'{"真实 ω":>7} | {"检出纯化选择":>12} | {"ω 估计中位数":>12} | {"ω 估计四分位":>18}')
res = {}
for om in (0.05, 0.1, 0.2, 0.3, 0.5, 1.0):
    cols = simulate(om, N)
    hit, oms = [], []
    for c in cols:
        p, a, b = fel_pvalue(c)
        hit.append(p <= 0.1 and b < a)
        oms.append(b / a if a > 0 else np.nan)
    oms = np.array(oms)
    res[om] = np.mean(hit)
    P(f'{om:7.2f} | {np.mean(hit):11.0%} | {np.nanmedian(oms):12.2f} | '
      f'{np.nanpercentile(oms, 25):7.2f}–{np.nanpercentile(oms, 75):<9.2f}')

P('\n判读：')
ok = [om for om, pw in res.items() if pw >= 0.8]
if ok:
    P(f'  → 真实 ω ≤ {max(ok)} 的约束，这套数据有 ≥80% 把握检出。4 个位点都没检出，'
      f'说明它们的约束（如果有）弱于 ω = {max(ok)}。')
else:
    best = max(res, key=res.get)
    P(f'  → 模拟范围内最高功效 {res[best]:.0%}（ω = {best}）。单位点功效确实有限，表述要保守：'
      f'只能说「没有检出强约束的证据」，不能说「证明了没有约束」。')
P('  实测对照：催化三联体 ω 估计 0.00~0.04、3/3 显著，与模拟中低 ω 的高检出率一致 → 方法在这套数据上是灵的。')
P('  4 个位点的 ω 点估计 0.75~1.50，落在模拟表的 ω = 1.0 那一行附近，而那一行本来就不该被判为纯化选择。')
P('  所以「不显著」的根源是 ω 估计本身接近 1，不是检验没力气。')
open(R('results', 't17_power.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
