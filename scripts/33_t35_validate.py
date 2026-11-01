"""任务 #35 第 0 步：核对 Python 复刻的 FEL 似然和 HyPhy 一致，并实测支长单位（每密码子 / 每核苷酸）。
用法: python 33_t35_validate.py <proj_root>
"""
import sys, os, json, time
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
sys.path.insert(0, os.path.dirname(__file__))
from t35_lib import Model, read_fasta, tip_matrix

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
fj = R('results', 't17_hyphy', 'fel_all.json')
d = json.load(open(fj))
mle = d['MLE']['content']['0']
seqs = read_fasta(R('results', 't17_codon', 'aln.fasta'))
SITES = [431, 479, 495, 498, 225, 227, 319, 390, 173, 12, 100, 300]
out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)

def fit(m, X):
    f1 = lambda x: -m.loglik(X, np.exp(x[0]), np.exp(x[1]))[0]
    f0 = lambda x: -m.loglik(X, np.exp(x[0]), np.exp(x[0]))[0]
    best1 = min((minimize(f1, x0, method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-5, maxiter=400))
                 for x0 in ([0, -0.5], [0, 0], [-1, -3])), key=lambda r: r.fun)
    r0 = minimize(f0, [0.0], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-5))
    return np.exp(best1.x), -best1.fun, -r0.fun, np.exp(r0.x[0])

for label, conv in (('每密码子', 1.0), ('每核苷酸', 3.0)):
    m = Model(fj)
    # t_b = L_b × conv / E，E = 全局模型（α=1, β=ω_g）下每单位 t 的每密码子替换数
    Q = m.SYN + m.NON * m.global_omega
    E = (m.pi[:, None] * Q).sum()
    m.bl_scale = conv / E
    P(f'\n== 支长约定：{label}（bl_scale = {m.bl_scale:.4f}）')
    P(f'{"位点":>5} | {"HyPhy α":>8} {"复刻 α":>8} | {"HyPhy β":>8} {"复刻 β":>8} | {"HyPhy LRT":>9} {"复刻 LRT":>9} | p(HyPhy) p(复刻)')
    t0 = time.time()
    for s in SITES:
        X = tip_matrix(m, seqs, [s - 1])
        (a, b), l1, l0, a0 = fit(m, X)
        h = mle[s - 1]
        lrt = max(2 * (l1 - l0), 0)
        P(f'{s:5d} | {h[0]:8.3f} {a:8.3f} | {h[1]:8.3f} {b:8.3f} | {h[3]:9.3f} {lrt:9.3f} | {h[4]:.4f} {chi2.sf(lrt, 1):.4f}')
    P(f'  用时 {time.time() - t0:.0f} s')
open(R('results', 't35', 'validate.txt'), 'w', encoding='utf-8').write('\n'.join(out) + '\n')
