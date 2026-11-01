"""10（任务 #13 · 13b/13c）: 合并扩采样后的位点表，在 Upham 2019 哺乳动物时间树上做祖先状态重建。
输入:
  results/site_table_v1.csv            任务 #11 的 RefSeq 位点表（只取「4位点可用=是」且「物种代表=是」）
  results/t13_secondsource.csv         13a：GenBank 注释的美洲豪猪 / 梳趾鼠
  results/t13_tblastn_sites.csv        13b：tblastn 读出的无注释组装（只取用途=13b 的）
  data/raw/t13/tree/*.tre              Upham et al. 2019 MCC 时间树（Completed 5911 种为主，DNA-only 4098 种做敏感性）
方法（每个位点独立）:
  1. Fitch 简约（无序、20 态）：两遍法给出节点的最简约状态集合 + 最少变化次数
  2. 最大似然 Mk-ER（等速率，状态 = 该位点所有物种里出现过的氨基酸），优化速率 q；
     节点边际概率 = 把该节点固定为状态 s 的似然 / 总似然（等根频率）
  3. 二态敏感性：裸鼹鼠型 vs 非裸鼹鼠型，Mk-ER 与 Mk-ARD（两个速率）
输出: results/t13_species_sites.csv、results/t13_asr_nodes.csv、results/t13_asr.txt
用法: python 10_t13_asr.py <proj_root>
"""
import sys, os, re, math
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, minimize
from Bio import Phylo

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
SITES = ['431', '479', '495', '498']
NMR = dict(zip(SITES, 'SEYT'))
HUM = dict(zip(SITES, 'DKLK'))
AA = 'ACDEFGHIKLMNPQRSTVWY'

# ---------- 1. 合并物种表 ----------
v1 = pd.read_csv(R('results', 'site_table_v1.csv'))
v1 = v1[(v1['4位点可用'] == '是') & (v1['物种代表'] == '是')]
rows = []
for _, r in v1.iterrows():
    sp = ' '.join(r.物种.split()[:2])
    rows.append(dict(物种=sp, 科=r.科, 目=r.目, 来源='RefSeq ' + r.RefSeq蛋白, **{s: r[f'{HUM[s]}{s}'] for s in SITES},
                     最大寿命_年=r.最大寿命_年, 成年体重_g=r.成年体重_g, LQ=r.LQ, 寿命来源=r.寿命来源 if isinstance(r.寿命来源, str) else ''))
ss = pd.read_csv(R('results', 't13_secondsource.csv'), dtype=str)
for acc, fam in (('KAM6166271.1', 'Erethizontidae'), ('KAM5256511.1', 'Ctenodactylidae')):
    r = ss[ss.accession == acc].iloc[0]
    rows.append(dict(物种=r.物种, 科=fam, 目='Rodentia', 来源='GenBank ' + acc, **{s: r[s] for s in SITES}))
TB_FAM = {'Hystrix cristata': 'Hystricidae', 'Hystrix brachyura': 'Hystricidae', 'Thryonomys swinderianus': 'Thryonomyidae',
          'Petromus typicus': 'Petromuridae', 'Fukomys anselli': 'Bathyergidae', 'Cavia tschudii': 'Caviidae', 'Cavia aperea': 'Caviidae',
          'Hydrochoerus hydrochaeris': 'Caviidae', 'Dolichotis patagonum': 'Caviidae', 'Dasyprocta punctata': 'Dasyproctidae',
          'Cuniculus paca': 'Cuniculidae', 'Dinomys branickii': 'Dinomyidae', 'Ctenomys sociabilis': 'Ctenomyidae',
          'Myocastor coypus': 'Echimyidae', 'Capromys pilorides': 'Echimyidae', 'Octomys mimax': 'Octodontidae',
          'Tympanoctomys barrerae': 'Octodontidae', 'Coendou prehensilis': 'Erethizontidae'}
tb = pd.read_csv(R('results', 't13_tblastn_sites.csv'), dtype=str)
tb = tb[tb.用途 == '13b']
for _, r in tb.iterrows():
    if r.状态 != 'ok':
        continue
    rows.append(dict(物种=r.物种, 科=TB_FAM[r.物种], 目='Rodentia', 来源=f'tblastn {r.WGS}（{r.组装}）', **{s: r[s] for s in SITES}))
df = pd.DataFrame(rows)
# Hystricognathi 新物种的 AnAge 寿命 + LQ（同 #11 的异速生长公式）
an = pd.read_csv(R('data', 'raw', 'anage', 'anage_data.txt'), sep='\t')
an = an[an.Class == 'Mammalia']
an['nm'] = an.Genus + ' ' + an.Species
fit = an[an['Maximum longevity (yrs)'].notna() & an['Adult weight (g)'].notna() & an['Data quality'].isin(['acceptable', 'high'])]
b, a0 = np.polyfit(np.log10(fit['Adult weight (g)']), np.log10(fit['Maximum longevity (yrs)']), 1)
SYN = {'Erethizon dorsatum': 'Erethizon dorsatus', 'Fukomys anselli': 'Cryptomys anselli'}
for i, r in df.iterrows():
    if r.来源.startswith('RefSeq'):
        continue
    x = an[an.nm == SYN.get(r.物种, r.物种)]
    if len(x) and x['Data quality'].iloc[0] != 'low':
        L, M = x['Maximum longevity (yrs)'].iloc[0], x['Adult weight (g)'].iloc[0]
        df.loc[i, ['最大寿命_年', '成年体重_g', '寿命来源']] = [L, M, 'AnAge ' + x['Data quality'].iloc[0] + ('（异名 ' + SYN[r.物种] + '）' if r.物种 in SYN else '')]
        if L == L and M == M:
            df.loc[i, 'LQ'] = round(L / 10 ** (a0 + b * math.log10(M)), 3)
df['命中'] = [sum(r[s] == NMR[s] for s in SITES) for _, r in df.iterrows()]
df = df.drop_duplicates('物种')

# ---------- 2. 树 ----------
def load_tree(fn):
    t = Phylo.read(fn, 'nexus')
    tips = {}
    for c in t.get_terminals():
        m = re.match(r'([A-Z][a-z]+)_([a-z]+)', c.name or '')
        if m:
            tips.setdefault(f'{m.group(1)} {m.group(2)}', c)
    return t, tips

TREE_SYN = {'Cryptomys damarensis': 'Fukomys damarensis', 'Neogale vison': 'Neovison vison', 'Cnephaeus nilssonii': 'Eptesicus nilssonii',
            'Notamacropus eugenii': 'Macropus eugenii', 'Herpailurus yagouaroundi': 'Puma yagouaroundi',
            'Nannospalax galili': 'Nannospalax ehrenbergi', 'Peromyscus californicus': 'Peromyscus californicus',
            'Alexandromys fortis': 'Microtus fortis', 'Megascapheus bottae': 'Thomomys bottae',
            # Upham 2019 用的是旧分类 / 近缘种代替（仅为上树；物种名映射不改变位点数据）
            'Nannospalax galili': 'Spalax ehrenbergi', 'Cricetulus griseus': 'Cricetulus barabensis', 'Grammomys surdaster': 'Grammomys dolichurus',
            'Pteropus medius': 'Pteropus giganteus', 'Bubalus bubalis': 'Bubalus arnee', 'Equus asinus': 'Equus africanus',
            'Equus przewalskii': 'Equus ferus', 'Cebus imitator': 'Cebus capucinus', 'Carlito syrichta': 'Tarsius syrichta',
            'Piliocolobus tephrosceles': 'Procolobus badius', 'Leucopleurus acutus': 'Lagenorhynchus acutus',
            'Aethalodelphis obliquidens': 'Lagenorhynchus obliquidens', 'Giraffa tippelskirchi': 'Giraffa camelopardalis'}

class N:
    __slots__ = ('ch', 'bl', 'name')
    def __init__(s, ch=(), bl=0.0, name=None):
        s.ch, s.bl, s.name = list(ch), bl, name

def prune(clade, keep):
    """返回只含 keep 里物种的子树（单子节点并入，枝长相加）"""
    if clade.is_terminal():
        m = re.match(r'([A-Z][a-z]+)_([a-z]+)', clade.name or '')
        nm = f'{m.group(1)} {m.group(2)}' if m else None
        return N(bl=clade.branch_length or 0.0, name=keep[nm]) if nm in keep else None
    kids = [k for k in (prune(c, keep) for c in clade.clades) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        kids[0].bl += clade.branch_length or 0.0
        return kids[0]
    return N(kids, clade.branch_length or 0.0)

def build(fn):
    t, tips = load_tree(fn)
    keep, miss = {}, []
    for sp in df.物种:
        tn = sp if sp in tips else TREE_SYN.get(sp)
        if tn in tips and tn not in keep:
            keep[tn] = sp
        else:
            miss.append(sp)
    root = prune(t.root, keep)
    root.bl = 0.0
    return root, miss

def postorder(n, out):
    for c in n.ch:
        postorder(c, out)
    out.append(n)
    return out

def tips_of(n):
    return [x.name for x in postorder(n, []) if not x.ch]

def mrca(root, names):
    names = set(names)
    best = None
    for x in postorder(root, []):
        if names <= set(tips_of(x)):
            if best is None or len(tips_of(x)) < len(tips_of(best)):
                best = x
    return best

# ---------- 3. Fitch ----------
def fitch(root, st):
    po = postorder(root, [])
    S, cost = {}, 0
    for n in po:
        if not n.ch:
            v = st[n.name]
            S[n] = set(AA) if v in '?-X' else {v}
        else:
            inter = set.intersection(*[S[c] for c in n.ch])
            if inter:
                S[n] = inter
            else:
                # 多叉：取出现次数最多的状态
                cnt = {}
                for c in n.ch:
                    for a in S[c]:
                        cnt[a] = cnt.get(a, 0) + 1
                mx = max(cnt.values())
                S[n] = {a for a, k in cnt.items() if k == mx}
                cost += len(n.ch) - mx
    # 第二遍（自上而下）：MPR 集合
    F = {root: S[root]}
    for n in reversed(po):
        for c in n.ch:
            inter = S[c] & F[n]
            F[c] = inter if inter else S[c]
    return F, cost

# ---------- 4. Mk 似然 ----------
def mk_P(Q, t):
    w, V = np.linalg.eig(Q)
    return np.real(V @ np.diag(np.exp(w * t)) @ np.linalg.inv(V))

def loglik(root, tipvec, Q, fix=None):
    """fix: (node, state_index) → 把该节点条件似然乘指示向量"""
    k = Q.shape[0]
    L, logsc = {}, 0.0
    cache = {}
    for n in postorder(root, []):
        if not n.ch:
            v = tipvec[n.name].copy()
        else:
            v = np.ones(k)
            for c in n.ch:
                key = round(c.bl, 9)
                if key not in cache:
                    cache[key] = mk_P(Q, c.bl)
                v *= cache[key] @ L[c]
        if fix is not None and n is fix[0]:
            m = np.zeros(k); m[fix[1]] = 1.0; v = v * m
        s = v.max()
        if s <= 0:
            return -np.inf
        L[n] = v / s
        logsc += math.log(s)
    return logsc + math.log(L[root].mean())

def tipvecs(states, st):
    idx = {a: i for i, a in enumerate(states)}
    out = {}
    for sp, v in st.items():
        x = np.zeros(len(states))
        if v in idx:
            x[idx[v]] = 1
        else:
            x[:] = 1
        out[sp] = x
    return out

def Q_er(k, q):
    Q = np.full((k, k), q)
    np.fill_diagonal(Q, -(k - 1) * q)
    return Q

def fit_er(root, tv, k):
    f = lambda lq: -loglik(root, tv, Q_er(k, math.exp(lq)))
    r = minimize_scalar(f, bounds=(math.log(1e-6), math.log(1.0)), method='bounded')
    return math.exp(r.x), -r.fun

def fit_ard2(root, tv):
    def f(p):
        a, b_ = math.exp(p[0]), math.exp(p[1])
        return -loglik(root, tv, np.array([[-a, a], [b_, -b_]]))
    r = minimize(f, x0=[math.log(1e-3), math.log(1e-3)], method='Nelder-Mead', options={'xatol': 1e-3, 'fatol': 1e-4})
    a, b_ = map(math.exp, r.x)
    return np.array([[-a, a], [b_, -b_]]), -r.fun, (a, b_)

def marginal(root, tv, Q, node, lnL):
    k = Q.shape[0]
    return np.array([math.exp(loglik(root, tv, Q, (node, i)) - lnL) for i in range(k)])

# ---------- 5. 跑 ----------
NODES = [
    ('啮齿目祖先 Rodentia', ['Mus musculus', 'Cavia porcellus', 'Sciurus vulgaris']),
    ('梳趾鼠+豪猪亚目 Ctenohystrica', ['Ctenodactylus gundi', 'Cavia porcellus', 'Heterocephalus glaber']),
    ('豪猪亚目祖先 Hystricognathi', ['Cavia porcellus', 'Heterocephalus glaber', 'Hystrix cristata']),
    ('非洲鼹鼠+蔗鼠+岩鼠（除豪猪科）', ['Heterocephalus glaber', 'Thryonomys swinderianus', 'Petromus typicus']),
    ('非洲鼹鼠科 Bathyergidae s.l.', ['Heterocephalus glaber', 'Fukomys damarensis']),
    ('南美支 Caviomorpha', ['Cavia porcellus', 'Erethizon dorsatum', 'Octodon degus', 'Chinchilla lanigera']),
    ('豚鼠总科 Cavioidea', ['Cavia porcellus', 'Dasyprocta punctata', 'Hydrochoerus hydrochaeris', 'Cuniculus paca']),
    ('豚鼠科 Caviidae', ['Cavia porcellus', 'Dolichotis patagonum', 'Hydrochoerus hydrochaeris']),
    ('豚鼠属 Cavia', ['Cavia porcellus', 'Cavia tschudii']),
    ('毛丝鼠总科 Chinchilloidea', ['Chinchilla lanigera', 'Dinomys branickii']),
    ('八齿鼠总科 Octodontoidea', ['Octodon degus', 'Myocastor coypus', 'Ctenomys sociabilis', 'Capromys pilorides']),
]
out_lines, node_rows = [], []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out_lines.append(s)

df.to_csv(R('results', 't13_species_sites.csv'), index=False, encoding='utf-8-sig')
P(f'物种总数 {len(df)}；其中啮齿目 {int((df.目 == "Rodentia").sum())}；豪猪亚目+梳趾鼠新增 {int((~df.来源.str.startswith("RefSeq")).sum())}')

for tag, fn in (('Completed5911', 'MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre'),
                ('DNAonly4098', 'MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre')):
    root, miss = build(R('data', 'raw', 't13', 'tree', fn))
    names = set(tips_of(root))
    P(f'\n===== 树 {tag}：{len(names)} 个物种上树；缺 {len(miss)}：{", ".join(miss)}')
    nodes = []
    for lab, sps in NODES:
        sps = [s for s in sps if s in names]
        if len(sps) >= 2:
            nd = mrca(root, sps)
            nodes.append((lab, nd))
            P(f'  节点 {lab}: 下含 {len(tips_of(nd))} 种（{", ".join(sorted(tips_of(nd))) if len(tips_of(nd)) <= 30 else "…"}）')
    def nwk(n):
        return n.name.replace(' ', '_') if not n.ch else '(' + ','.join(nwk(c) for c in n.ch) + f'):{n.bl:.1f}'
    ct = [x for l, x in nodes if l.startswith('梳趾鼠')]
    if ct:
        P('  梳趾鼠+豪猪亚目子树拓扑（括号后数字 = 该节点上方枝长，百万年）：' + nwk(ct[0]))
    for s in SITES:
        st = {sp: df.set_index('物种').loc[sp, s] for sp in names}
        F, cost = fitch(root, st)
        states = sorted({v for v in st.values() if v in AA})
        tv = tipvecs(states, st)
        q, lnL = fit_er(root, tv, len(states))
        Q = Q_er(len(states), q)
        # 二态
        tvb = {sp: (np.array([0., 1.]) if v == NMR[s] else (np.array([1., 1.]) if v not in AA else np.array([1., 0.]))) for sp, v in st.items()}
        qb, lnLb = fit_er(root, tvb, 2)
        Qa, lnLa, (ra, rb) = fit_ard2(root, tvb)
        P(f'\n  位点 {HUM[s]}{s}（裸鼹鼠型 {NMR[s]}）：Fitch 最少变化 {cost} 次；Mk-ER {len(states)} 态 q={q:.4g}/Myr lnL={lnL:.2f}；'
          f'二态 ER q={qb:.4g} lnL={lnLb:.2f}；二态 ARD 得={ra:.3g} 失={rb:.3g} lnL={lnLa:.2f}')
        for lab, nd in nodes:
            pm = marginal(root, tv, Q, nd, lnL)
            top = sorted(zip(states, pm), key=lambda x: -x[1])[:3]
            pb = marginal(root, tvb, Q_er(2, qb), nd, lnLb)[1]
            pa = marginal(root, tvb, Qa, nd, lnLa)[1]
            P(f'    {lab:40s} Fitch={"".join(sorted(F[nd])):6s} ML: ' + ' '.join(f'{a}={p:.3f}' for a, p in top) +
              f' | P(裸鼹鼠型) 多态ER={dict(zip(states, pm)).get(NMR[s], 0):.3f} 二态ER={pb:.3f} 二态ARD={pa:.3f}')
            node_rows.append(dict(树=tag, 位点=f'{HUM[s]}{s}{NMR[s]}', 节点=lab, Fitch_MPR=''.join(sorted(F[nd])),
                                  ML最可能=top[0][0], ML最可能概率=round(top[0][1], 4),
                                  P裸鼹鼠型_多态ER=round(dict(zip(states, pm)).get(NMR[s], 0), 4),
                                  P裸鼹鼠型_二态ER=round(pb, 4), P裸鼹鼠型_二态ARD=round(pa, 4), Fitch最少变化=cost))
nr = pd.DataFrame(node_rows)
nr.to_csv(R('results', 't13_asr_nodes.csv'), index=False, encoding='utf-8-sig')
# 4 个位点同时为裸鼹鼠型的概率（按位点独立近似相乘）
P('\n===== 4 位点同时为裸鼹鼠型（位点独立近似，P 相乘）')
for (tag, lab), g in nr.groupby(['树', '节点'], sort=False):
    P(f'  {tag:14s} {lab:40s} 多态ER={g.P裸鼹鼠型_多态ER.prod():.3f} 二态ER={g.P裸鼹鼠型_二态ER.prod():.3f} 二态ARD={g.P裸鼹鼠型_二态ARD.prod():.3f}')
open(R('results', 't13_asr.txt'), 'w', encoding='utf-8').write('\n'.join(out_lines) + '\n')
