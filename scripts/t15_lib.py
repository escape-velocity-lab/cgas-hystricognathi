"""任务 #15 共用：Upham 2019 树的读取 / 剪枝、Fitch 简约、Mk-ER 最大似然边际祖先状态。
函数照抄 10_t13_asr.py（那个脚本 import 即运行，不能直接 import），行为保持一致；新增 marginal_all（一次算所有节点）。
"""
import re, math
import numpy as np
from scipy.optimize import minimize_scalar
from Bio import Phylo

AA = 'ACDEFGHIKLMNPQRSTVWY'
TREE_MAIN = 'MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre'
TREE_DNA = 'MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre'

# 与 10_t13_asr.py 的 TREE_SYN 相同（物种名 → Upham 树里的名字）
TREE_SYN = {'Cryptomys damarensis': 'Fukomys damarensis', 'Neogale vison': 'Neovison vison', 'Cnephaeus nilssonii': 'Eptesicus nilssonii',
            'Notamacropus eugenii': 'Macropus eugenii', 'Herpailurus yagouaroundi': 'Puma yagouaroundi',
            'Alexandromys fortis': 'Microtus fortis', 'Megascapheus bottae': 'Thomomys bottae',
            'Nannospalax galili': 'Spalax ehrenbergi', 'Cricetulus griseus': 'Cricetulus barabensis', 'Grammomys surdaster': 'Grammomys dolichurus',
            'Pteropus medius': 'Pteropus giganteus', 'Bubalus bubalis': 'Bubalus arnee', 'Equus asinus': 'Equus africanus',
            'Equus przewalskii': 'Equus ferus', 'Cebus imitator': 'Cebus capucinus', 'Carlito syrichta': 'Tarsius syrichta',
            'Piliocolobus tephrosceles': 'Procolobus badius', 'Leucopleurus acutus': 'Lagenorhynchus acutus',
            'Aethalodelphis obliquidens': 'Lagenorhynchus obliquidens', 'Giraffa tippelskirchi': 'Giraffa camelopardalis'}


class N:
    __slots__ = ('ch', 'bl', 'name', 'up')
    def __init__(s, ch=(), bl=0.0, name=None):
        s.ch, s.bl, s.name, s.up = list(ch), bl, name, None


def _tipname(c):
    m = re.match(r'([A-Z][a-z]+)_([a-z]+)', c.name or '')
    return f'{m.group(1)} {m.group(2)}' if m else None


def _prune(clade, keep):
    if clade.is_terminal():
        nm = _tipname(clade)
        return N(bl=clade.branch_length or 0.0, name=keep[nm]) if nm in keep else None
    kids = [k for k in (_prune(c, keep) for c in clade.clades) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        kids[0].bl += clade.branch_length or 0.0
        return kids[0]
    return N(kids, clade.branch_length or 0.0)


def build(fn, species):
    """species: 物种名列表 → (根节点, 缺失物种)。叶子名 = 我们的物种名。"""
    t = Phylo.read(fn, 'nexus')
    tips = {}
    for c in t.get_terminals():
        nm = _tipname(c)
        if nm:
            tips.setdefault(nm, c)
    keep, miss = {}, []
    for sp in species:
        tn = sp if sp in tips else TREE_SYN.get(sp)
        if tn in tips and tn not in keep:
            keep[tn] = sp
        else:
            miss.append(sp)
    root = _prune(t.root, keep)
    root.bl = 0.0
    for n in postorder(root):
        for c in n.ch:
            c.up = n
    return root, miss


def postorder(n, out=None):
    out = [] if out is None else out
    for c in n.ch:
        postorder(c, out)
    out.append(n)
    return out


def tips_of(n):
    return [x.name for x in postorder(n) if not x.ch]


def mrca(root, names):
    names = set(names)
    best = None
    for x in postorder(root):
        tx = set(tips_of(x))
        if names <= tx and (best is None or len(tx) < len(tips_of(best))):
            best = x
    return best


def depth(n):
    d = 0.0
    while n.up is not None:
        d += n.bl; n = n.up
    return d


def fitch(root, st):
    po = postorder(root)
    S, cost = {}, 0
    for n in po:
        if not n.ch:
            v = st.get(n.name, '?')
            S[n] = set(AA) if v not in AA else {v}
        else:
            inter = set.intersection(*[S[c] for c in n.ch])
            if inter:
                S[n] = inter
            else:
                cnt = {}
                for c in n.ch:
                    for a in S[c]:
                        cnt[a] = cnt.get(a, 0) + 1
                mx = max(cnt.values())
                S[n] = {a for a, k in cnt.items() if k == mx}
                cost += len(n.ch) - mx
    F = {root: S[root]}
    for n in reversed(po):
        for c in n.ch:
            inter = S[c] & F[n]
            F[c] = inter if inter else S[c]
    return F, cost


def _P_er(k, q, t):
    """等速率 Mk 的转移矩阵有解析解：P_ii = 1/k + (k-1)/k·e^{-kqt}，P_ij = 1/k − 1/k·e^{-kqt}"""
    e = math.exp(-k * q * t)
    P = np.full((k, k), (1 - e) / k)
    np.fill_diagonal(P, 1 / k + (k - 1) / k * e)
    return P


def _down(root, tv, k, q):
    """自下而上条件似然（每节点归一化），返回 (L, logscale, Pcache)"""
    L, logsc, P = {}, 0.0, {}
    for n in postorder(root):
        if not n.ch:
            v = tv[n.name].copy()
        else:
            v = np.ones(k)
            for c in n.ch:
                if c not in P:
                    P[c] = _P_er(k, q, c.bl)
                v *= P[c] @ L[c]
        s = v.max()
        L[n] = v / s
        logsc += math.log(s)
    return L, logsc, P


def tipvecs(states, st, names):
    idx = {a: i for i, a in enumerate(states)}
    out = {}
    for sp in names:
        v = st.get(sp, '?')
        x = np.zeros(len(states))
        if v in idx:
            x[idx[v]] = 1
        else:
            x[:] = 1
        out[sp] = x
    return out


def fit_er(root, tv, k):
    def f(lq):
        L, ls, _ = _down(root, tv, k, math.exp(lq))
        return -(ls + math.log(L[root].mean()))
    r = minimize_scalar(f, bounds=(math.log(1e-6), math.log(1.0)), method='bounded')
    return math.exp(r.x), -r.fun


def marginal_all(root, st):
    """Mk-ER 多态（状态 = 该位点出现过的氨基酸，等根频率）。返回 (states, {node: 概率向量}, q, lnL)。
    上下两遍法：节点边际 ∝ 下方似然 × 上方似然。和 10 的「固定节点重算似然」结果相同，但一次算完所有节点。"""
    names = tips_of(root)
    states = sorted({st.get(s) for s in names} & set(AA))
    if len(states) < 2:
        return states, {n: np.ones(1) for n in postorder(root)}, 0.0, 0.0
    k = len(states)
    tv = tipvecs(states, st, names)
    q, lnL = fit_er(root, tv, k)
    L, _, P = _down(root, tv, k, q)
    U = {root: np.full(k, 1.0 / k)}           # 上方（含根频率）消息
    out = {}
    for n in reversed(postorder(root)):
        post = U[n] * L[n]
        out[n] = post / post.sum()
        for c in n.ch:
            sib = np.ones(k)
            for d in n.ch:
                if d is not c:
                    sib *= P[d] @ L[d]
            m = U[n] * sib
            u = P[c].T @ m
            U[c] = u / u.max()
    return states, out, q, lnL
