"""任务 #35 B：四个位点共享 ω 的 profile-likelihood 95% 置信区间（审稿人 3 的 M3-2b）。

模型：FEL 框架（支长、核苷酸速率、密码子频率固定在 HyPhy 全局拟合值，t35_lib 复刻、已逐位点核对），
组内每个位点有自己的同义速率 α_s，非同义速率 β_s = ω·α_s，ω 全组共享。
  profile ℓ(ω) = Σ_s max_{α_s} lnL_s(α_s, ω·α_s)；95% CI = {ω : 2[ℓ(ω̂) − ℓ(ω)] ≤ 3.841}
同时报告：每个位点单独的 profile CI；「4 个 ω 各自自由 vs 共享」的异质性 LRT（df = n−1）。
对照组：催化三联体、锌指配体、11 个调控性残基（任务 C）。
用法: python 36_t35_profile.py <proj_root> <n_proc>
"""
import sys, os, json, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize_scalar, brentq
from scipy.stats import chi2
sys.path.insert(0, os.path.dirname(__file__))
from t35_lib import Model, read_fasta, tip_matrix

ROOT, NPROC = sys.argv[1], int(sys.argv[2])
R = lambda *p: os.path.join(ROOT, *p)
FJ = R('results', 't17_hyphy', 'fel_all.json')
GROUPS = {
    '四个位点': [431, 479, 495, 498],
    '催化三联体': [225, 227, 319],
    '锌指配体': [390, 396, 397, 404],
    '调控性残基（11）': [347, 350, 275, 279, 282, 285, 427, 428, 432, 236, 255],
}
CRIT = chi2.ppf(0.95, 1)
G = {}


def init():
    os.environ['OMP_NUM_THREADS'] = '1'
    m = Model(FJ)
    Q = m.SYN + m.NON * m.global_omega
    m.bl_scale = 3.0 / (m.pi[:, None] * Q).sum()
    seqs = read_fasta(R('results', 't17_codon', 'aln.fasta'))
    G['m'] = m
    G['X'] = {s: tip_matrix(m, seqs, [s - 1]) for g in GROUPS.values() for s in g}


def site_prof(args):
    """给定位点和 ω，对 log α 求最大似然"""
    s, om = args
    m, X = G['m'], G['X'][s]
    f = lambda la: -m.loglik(X, np.exp(la), om * np.exp(la))[0]
    r = minimize_scalar(f, bounds=(np.log(1e-3), np.log(50)), method='bounded', options=dict(xatol=1e-4))
    return s, om, -r.fun, float(np.exp(r.x))


if __name__ == '__main__':
    t0 = time.time()
    grid = np.concatenate([[0.0], np.exp(np.linspace(np.log(0.003), np.log(10), 48))])
    sites = sorted({s for g in GROUPS.values() for s in g})
    lines = []
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

    with Pool(NPROC, initializer=init) as pool:
        tab = {}
        for s, om, l, a in pool.imap_unordered(site_prof, [(s, om) for s in sites for om in grid], chunksize=2):
            tab[(s, om)] = (l, a)
        P(f'网格 {len(grid)} 点 × {len(sites)} 个位点完成，{time.time() - t0:.0f} s')

        def prof(group, om):
            if all((s, om) in tab for s in group):
                return sum(tab[(s, om)][0] for s in group)
            res = pool.map(site_prof, [(s, om) for s in group])
            for s, o, l, a in res:
                tab[(s, o)] = (l, a)
            return sum(r[2] for r in res)

        def profile_ci(group):
            ys = np.array([prof(group, om) for om in grid])
            k = int(np.argmax(ys))
            # 在网格最大值附近用黄金分割细化 ω̂
            lo_b = grid[max(k - 1, 0)]; hi_b = grid[min(k + 1, len(grid) - 1)]
            if hi_b > lo_b:
                r = minimize_scalar(lambda o: -prof(group, o), bounds=(lo_b, hi_b), method='bounded',
                                    options=dict(xatol=1e-4))
                om_hat, lmax = (r.x, -r.fun) if -r.fun > ys[k] else (grid[k], ys[k])
            else:
                om_hat, lmax = grid[k], ys[k]
            g = lambda o: 2 * (lmax - prof(group, o)) - CRIT
            # 下限
            if om_hat <= 0 or g(0.0) <= 0:
                lo = 0.0
            else:
                below = [o for o in grid if o < om_hat and 2 * (lmax - prof(group, o)) > CRIT]
                lo = brentq(g, max(below), om_hat, xtol=1e-4)
            above = [o for o in grid if o > om_hat and 2 * (lmax - prof(group, o)) > CRIT]
            hi = brentq(g, om_hat, min(above), xtol=1e-4) if above else np.inf
            # 关键参考点上的似然比
            ref = {o: 2 * (lmax - prof(group, o)) for o in (0.02, 0.1, 0.2, 0.3, 0.5, 1.0)}
            return om_hat, lmax, lo, hi, ref

        out = {}
        P(f'\n共享 ω 的 profile likelihood（每个位点自己的 α；FEL 全局参数固定）')
        P(f'{"组":<16} {"n":>2} | {"ω̂":>6} | {"95% CI":>16} | 2ΔlnL 在 ω = 0.2 / 0.3 / 0.5 / 1.0（> 3.84 即被排除）| 异质性 LRT')
        for name, grp in GROUPS.items():
            om_hat, lmax, lo, hi, ref = profile_ci(grp)
            # 异质性：每个位点自由 ω（就是 FEL 的 H1）
            free = 0.0; per = {}
            for s in grp:
                o_s, l_s, lo_s, hi_s, _ = profile_ci([s])
                free += l_s; per[s] = (o_s, lo_s, hi_s)
            het = 2 * (free - lmax); df = len(grp) - 1
            ph = chi2.sf(max(het, 0), df) if df > 0 else np.nan
            P(f'{name:<16} {len(grp):2d} | {om_hat:6.3f} | [{lo:6.3f}, {hi:6.3f}] | '
              f'{ref[0.2]:6.1f} / {ref[0.3]:6.1f} / {ref[0.5]:6.1f} / {ref[1.0]:6.1f} | '
              f'{het:.1f} (df {df}, p = {ph:.3g})')
            for s in grp:
                o_s, lo_s, hi_s = per[s]
                P(f'      位点 {s}: ω̂ = {o_s:.3f}，单独 95% CI [{lo_s:.3f}, {hi_s:.3f}]')
            out[name] = dict(sites=grp, omega_hat=om_hat, ci=[lo, hi], lr_at=ref, het_lrt=het, het_df=df, het_p=ph,
                             per_site={str(s): v for s, v in per.items()})
    P(f'\n总耗时 {time.time() - t0:.0f} s')
    json.dump(out, open(R('results', 't35', 'profile.json'), 'w'), ensure_ascii=False, indent=1, default=float)
    open(R('results', 't35', 'profile.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
