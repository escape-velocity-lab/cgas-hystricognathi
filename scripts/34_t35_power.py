"""任务 #35 A：在与主结果完全相同的判定流程（HyPhy FEL 全对齐 + 有信息位点内 BH，q ≤ 0.05 且 β < α）下重估功效。

旧版（18_t17_power.py）的两个问题：
  1. 判定阈值是未校正 p ≤ 0.1，主结果用的是 BH q ≤ 0.05（名义阈值约 0.02）
  2. 把 HyPhy 的支长（每核苷酸替换数）当成每密码子用，模拟的进化量只有真实的 1/3（33_t35_validate.py 实测）

做法：
  make   —— 用 HyPhy 全局拟合的参数（GTR 速率、CF3x4 频率、MG94 支长，t35_lib 复刻、已与 HyPhy 逐位点核对）
            模拟位点，α 取四个位点 HyPhy 估计的均值 0.65，β = ω·α；把模拟列**接在真实 522 列后面**写成新对齐，交给 HyPhy FEL
  python —— 同样模拟，但用复刻的 FEL 似然比检验（全局参数固定在真实拟合值），样本量更大，作交叉核对
  analyze —— 每个模拟位点：真实有信息位点的 p 值（同一次 HyPhy 运行）+ 该模拟位点 → BH → q ≤ 0.05 且 β < α 判纯化
用法: python 34_t35_power.py <proj_root> make <run_id> <seed>
      python 34_t35_power.py <proj_root> python <n_per_omega> <seed> <n_proc>
      python 34_t35_power.py <proj_root> analyze
"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from t35_lib import Model, read_fasta, tip_matrix_from_idx, SENSE

ROOT, MODE = sys.argv[1], sys.argv[2]
R = lambda *p: os.path.join(ROOT, *p)
OUT = R('results', 't35')
FJ = R('results', 't17_hyphy', 'fel_all.json')
ALPHA = float(os.environ.get('T35_ALPHA', 0.65))   # 四个位点 HyPhy α：0.537 / 0.647 / 0.795 / 0.631，均值 0.65；敏感性分析用环境变量改
SUFFIX = os.environ.get('T35_SUFFIX', '')
LEVELS = [(0.1, 50), (0.2, 50), (0.3, 50), (0.5, 50), (1.0, 25)]   # 每次 HyPhy 运行嵌入的数量
OMEGAS = [0.1, 0.2, 0.3, 0.5, 1.0]


def model():
    m = Model(FJ)
    Q = m.SYN + m.NON * m.global_omega
    m.bl_scale = 3.0 / (m.pi[:, None] * Q).sum()     # HyPhy 支长 = 每核苷酸替换数（validate 实测）
    return m


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    ph = k / n; den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den; h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return (max(0, c - h), min(1, c + h))



from scipy.optimize import minimize
from scipy.stats import chi2
_W = {}


def _init():
    os.environ['OMP_NUM_THREADS'] = '1'
    _W['m'] = model()


def test_col(args):
    om, j, col = args
    m = _W['m']
    X = tip_matrix_from_idx(col[:, None])
    f1 = lambda x: -m.loglik(X, np.exp(x[0]), np.exp(x[1]))[0]
    f0 = lambda x: -m.loglik(X, np.exp(x[0]), np.exp(x[0]))[0]
    r1 = min((minimize(f1, x0, method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-4, maxiter=300))
              for x0 in ([np.log(ALPHA), np.log(ALPHA * 0.5)], [np.log(ALPHA), np.log(ALPHA * 0.05)])),
             key=lambda r: r.fun)
    r0 = minimize(f0, [np.log(ALPHA)], method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-4))
    lrt = max(2 * (r0.fun - r1.fun), 0)
    a, b = np.exp(r1.x)
    return om, j, a, b, chi2.sf(lrt, 1)


if __name__ == '__main__':
    if MODE == 'make':
        run_id, seed = sys.argv[3], int(sys.argv[4])
        rng = np.random.default_rng(seed)
        m = model()
        seqs = read_fasta(R('results', 't17_codon', 'aln.fasta'))
        cols, labels = [], []
        for om, n in LEVELS:
            idx = m.simulate(ALPHA, ALPHA * om, n, rng)       # (n_tips, n)
            cols.append(idx); labels += [om] * n
        idx = np.concatenate(cols, axis=1)
        order = rng.permutation(idx.shape[1])                 # 打乱顺序，避免同一 ω 扎堆
        idx = idx[:, order]; labels = [labels[i] for i in order]
        d = os.path.join(OUT, f'power_{run_id}'); os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, 'aln.fasta'), 'w') as f:
            for k, name in enumerate(m.tipname):
                f.write(f'>{name}\n{seqs[name]}{"".join(SENSE[c] for c in idx[k])}\n')
        n_real = len(seqs[m.tipname[0]]) // 3
        json.dump(dict(seed=seed, alpha=ALPHA, n_real=n_real, sim_omega=labels), open(os.path.join(d, 'meta.json'), 'w'))
        print(f'{d}: 真实 {n_real} 列 + 模拟 {len(labels)} 列，种子 {seed}')

    elif MODE == 'python':
        from multiprocessing import Pool
        n_per, seed, nproc = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        rng = np.random.default_rng(seed)
        m = model()
        jobs = []
        for om in OMEGAS:
            idx = m.simulate(ALPHA, ALPHA * om, n_per, rng)
            jobs += [(om, j, idx[:, j]) for j in range(n_per)]
        t0 = time.time()
        rows = []
        with Pool(nproc, initializer=_init) as pool:
            for k, r in enumerate(pool.imap_unordered(test_col, jobs, chunksize=4)):
                rows.append(r)
                if (k + 1) % 100 == 0:
                    print(f'  {k + 1}/{len(jobs)}  {time.time() - t0:.0f} s', flush=True)
        import csv
        with open(os.path.join(OUT, f'power_python{SUFFIX}.csv'), 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['omega', 'rep', 'alpha', 'beta', 'p'])
            w.writerows(sorted(rows))
        print(f'完成 {len(rows)} 个位点，{time.time() - t0:.0f} s，种子 {seed}')

    elif MODE == 'analyze':
        import csv
        lines = []
        def P(*a):
            s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

        def informative(content):
            return [i for i, r in enumerate(content) if r[0] > 0 or r[1] > 0]

        orig = json.load(open(FJ))['MLE']['content']['0']
        p_orig = np.array([orig[i][4] for i in informative(orig)])
        P(f'主结果：有信息位点 {len(p_orig)} 个，BH q ≤ 0.05 被拒 {(bh(p_orig) <= 0.05).sum()} 个，'
          f'被拒位点中最大的原始 p = {p_orig[bh(p_orig) <= 0.05].max():.4f}（这就是实际的名义阈值）')

        def decide(p_real, p_sim, a, b):
            """把一个模拟位点放进真实位点里做 BH：返回 (q, 判纯化与否)"""
            q = bh(np.append(p_real, p_sim))[-1]
            return q, (q <= 0.05) and (b < a)

        # ---- HyPhy 嵌入
        hy = {om: [] for om in OMEGAS}
        for d in sorted(x for x in os.listdir(OUT) if x.startswith('power_run') and not SUFFIX):   # HyPhy 嵌入只做了 α = 0.65
            jf = os.path.join(OUT, d, 'fel.json')
            if not os.path.exists(jf):
                P(f'  [缺] {d}/fel.json'); continue
            meta = json.load(open(os.path.join(OUT, d, 'meta.json')))
            try:
                c = json.load(open(jf))['MLE']['content']['0']
            except Exception:
                P(f'  [未完成] {d}/fel.json'); continue
            nr = meta['n_real']
            real = c[:nr]
            p_real = np.array([real[i][4] for i in informative(real)])
            dr = np.corrcoef(p_real, p_orig)[0, 1] if len(p_real) == len(p_orig) else np.nan
            P(f'  {d}：种子 {meta["seed"]}，真实部分有信息位点 {len(p_real)} 个（与主结果 p 值相关 r = {dr:.4f}），'
              f'BH 被拒 {(bh(p_real) <= 0.05).sum()} 个')
            for k, om in enumerate(meta['sim_omega']):
                a, b, _, _, p = c[nr + k][:5]
                q, hit = decide(p_real, p, a, b)
                hy[om].append((a, b, p, q, hit, p <= 0.1 and b < a))
        # ---- Python 复刻
        py = {om: [] for om in OMEGAS}
        pf = os.path.join(OUT, f'power_python{SUFFIX}.csv')
        if os.path.exists(pf):
            for r in csv.DictReader(open(pf)):
                a, b, p = float(r['alpha']), float(r['beta']), float(r['p'])
                q, hit = decide(p_orig, p, a, b)
                py[float(r['omega'])].append((a, b, p, q, hit, p <= 0.1 and b < a))

        res = {}
        for name, dat in (('HyPhy FEL（模拟位点嵌入真实对齐，完整流程）', hy), ('Python 复刻 FEL（全局参数固定在真实拟合值）', py)):
            if not any(dat.values()):
                continue
            P(f'\n== {name}；模拟 α = {ALPHA}')
            P(f'{"真实ω":>6} | {"n":>4} | {"BH q≤0.05 且 β<α（主结果判据）":>28} | {"95% CI":>13} | {"旧判据 p≤0.1":>10} | ω 估计中位数（四分位）| 4 位全漏检 (1−π)^4')
            res[name] = {}
            for om in OMEGAS:
                v = dat[om]
                if not v:
                    continue
                n = len(v); k = sum(x[4] for x in v); k01 = sum(x[5] for x in v)
                lo, hi = wilson(k, n)
                oms = np.array([x[1] / x[0] if x[0] > 0 else np.nan for x in v])
                pw = k / n
                res[name][om] = dict(n=n, power=pw, lo=lo, hi=hi, power_p01=k01 / n)
                P(f'{om:6.2f} | {n:4d} | {pw:27.0%} | {lo:5.0%}–{hi:<6.0%} | {k01 / n:10.0%} | '
                  f'{np.nanmedian(oms):.2f}（{np.nanpercentile(oms, 25):.2f}–{np.nanpercentile(oms, 75):.2f}）| '
                  f'{(1 - pw) ** 4:.4f}（CI {(1 - hi) ** 4:.4f}–{(1 - lo) ** 4:.4f}）')
        json.dump({k: {str(o): v for o, v in d.items()} for k, d in res.items()},
                  open(os.path.join(OUT, f'power_summary{SUFFIX}.json'), 'w'), ensure_ascii=False, indent=1)
        open(os.path.join(OUT, f'power{SUFFIX}.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
