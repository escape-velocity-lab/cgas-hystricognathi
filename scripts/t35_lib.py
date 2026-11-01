"""任务 #35 公用库：用 Python 复刻 HyPhy FEL 的单位点似然（MG94×REV + CF3x4），参数全部从 HyPhy 全局拟合里读。

为什么要复刻：功效模拟必须在「和 HyPhy 一模一样的模型、支长单位」下生成数据，profile likelihood 也要在同一模型下做。
复刻得对不对，用真实位点核对：对真实的列做 α、β 最大似然，必须和 fel_all.json 里 HyPhy 的 α、β、LRT 对上
（见 33_t35_validate.py）。

模型（HyPhy models/codon/MG_REV.bf 的写法）：
  Q[i→j] = θ_xy × f_p(y) × (α 若同义 / β 若非同义)，i、j 只在密码子第 p 位不同（x→y），θ 为 GTR 核苷酸速率
  f_p 为第 p 位核苷酸频率（CF3x4 修正后），由 HyPhy 报告的 61 个平衡频率反推（平衡频率 ∝ f_1·f_2·f_3）
FEL 在每条支上：α_branch = α_site × t_b，β_branch = β_site × t_b，t_b 是全局拟合中该支的同义速率参数。
HyPhy 报告的「支长」L_b 是 t_b 换算出的期望替换数，换算常数在 validate 里实测确定（BL_PER_T）。
"""
import json, os
import numpy as np
from Bio import Phylo
from Bio.Data.CodonTable import unambiguous_dna_by_id
from io import StringIO

BASES = 'ACGT'
std = unambiguous_dna_by_id[1]
SENSE = sorted(a + b + c for a in 'ACGT' for b in 'ACGT' for c in 'ACGT'
               if a + b + c not in std.stop_codons)                        # HyPhy 的字母序：AAA, AAC, …
NC = len(SENSE)
CIDX = {c: i for i, c in enumerate(SENSE)}
AA = [std.forward_table[c] for c in SENSE]


class Model:
    def __init__(self, fel_json, bl_scale=None):
        d = json.load(open(fel_json))
        g = d['fits']['Global MG94xREV']
        self.pi = np.array([x[0] for x in g['Equilibrium frequencies']])
        self.pi /= self.pi.sum()
        rd = d['fits']['Nucleotide GTR']['Rate Distributions']
        th = np.zeros((4, 4))
        for k, v in rd.items():
            x, y = k.split('nucleotide ')[1][0], k.split('nucleotide ')[2][0]
            th[BASES.index(x), BASES.index(y)] = th[BASES.index(y), BASES.index(x)] = v
        # 注意：MG94xREV 的核苷酸速率是全局拟合重新估的，JSON 里只给了 GTR 阶段的值；validate 会检查这一点带来的偏差
        self.theta = th
        self.f = self._positional_freqs()
        self._skeleton()
        # 树：输入树字符串带内部节点名，支长取全局 MG94xREV 拟合值
        tstr = d['input']['trees']['0']
        self.tree = Phylo.read(StringIO(tstr), 'newick')
        ba = d['branch attributes']['0']
        for c in self.tree.find_clades():
            if c is self.tree.root:
                c.branch_length = 0.0
                continue
            c.branch_length = ba.get(c.name, {}).get('Global MG94xREV', 0.0)   # HyPhy 删掉的 30 条零长内部支 → 0
        self.post = list(self.tree.find_clades(order='postorder'))
        self.tips = self.tree.get_terminals()
        self.tipname = [t.name for t in self.tips]
        self.L = np.array([c.branch_length for c in self.post])
        self.bl_scale = bl_scale          # 支长 → 「α=β=1 时的时间单位」的换算，validate 定
        self.global_omega = list(g['Rate Distributions'].values())[0][0][0]

    def _positional_freqs(self):
        # log π_c = log f1(a) + log f2(b) + log f3(c) + const，最小二乘反解（F3x4 乘积形式下是精确解）
        X = np.zeros((NC, 12))
        for i, c in enumerate(SENSE):
            for p in range(3):
                X[i, 4 * p + BASES.index(c[p])] = 1
        coef, *_ = np.linalg.lstsq(X, np.log(self.pi), rcond=None)
        f = np.exp(coef).reshape(3, 4)
        f /= f.sum(1, keepdims=True)
        return f

    def _skeleton(self):
        SYN = np.zeros((NC, NC)); NON = np.zeros((NC, NC))
        for i, ci in enumerate(SENSE):
            for j, cj in enumerate(SENSE):
                dpos = [k for k in range(3) if ci[k] != cj[k]]
                if len(dpos) != 1:
                    continue
                p = dpos[0]
                r = self.theta[BASES.index(ci[p]), BASES.index(cj[p])] * self.f[p, BASES.index(cj[p])]
                (SYN if AA[i] == AA[j] else NON)[i, j] = r
        self.SYN, self.NON = SYN, NON
        # α=β=1 时每单位时间的期望替换数（每密码子）
        self.rate11 = (self.pi[:, None] * (SYN + NON)).sum()

    def eig(self, alpha, beta):
        Q = self.SYN * alpha + self.NON * beta
        np.fill_diagonal(Q, -Q.sum(1))
        # 可逆模型：对称化后用 eigh，数值稳定
        s = np.sqrt(self.pi)
        S = (s[:, None] * Q) / s[None, :]
        S = (S + S.T) / 2
        lam, V = np.linalg.eigh(S)
        U = V / s[:, None]          # 右特征向量
        Ui = V.T * s[None, :]       # 左特征向量
        return lam, U, Ui

    def times(self):
        return self.L * self.bl_scale

    def loglik(self, tipstates, alpha, beta):
        """tipstates: (n_tips, n_cols) 的列表，元素是 61 维的指示向量矩阵 → 这里用 (n_tips, NC, n_cols) 的布尔数组。
        返回每列的 lnL。"""
        lam, U, Ui = self.eig(alpha, beta)
        T = self.times()
        ncol = tipstates.shape[2]
        partial = {}
        scale = np.zeros(ncol)
        ti = {t: k for k, t in enumerate(self.tips)}
        for idx, c in enumerate(self.post):
            if c.is_terminal():
                v = tipstates[ti[c]].astype(float)
            else:
                v = np.ones((NC, ncol))
                for ch in c.clades:
                    t = ch.branch_length * self.bl_scale
                    v = v * (U @ (np.exp(lam * t)[:, None] * (Ui @ partial.pop(ch))))
                v = np.maximum(v, 0)            # 特征分解的舍入误差会给出极小的负数
                m = v.max(0)
                dead = m <= 0                   # 该列在此 (α, β) 下不可能（如 β=0 而列里有多种氨基酸）
                m[dead] = 1.0
                v = v / m
                scale += np.log(m)
                scale[dead] = -1e10
            partial[c] = v
        root = partial[self.post[-1]]
        return np.maximum(scale + np.log(np.maximum(self.pi @ root, 1e-300)), -1e10)

    def simulate(self, alpha, beta, n, rng):
        """在同一棵树上按 (α, β) 模拟 n 列，返回 (n_tips, n) 的密码子下标"""
        lam, U, Ui = self.eig(alpha, beta)
        pre = list(self.tree.find_clades(order='preorder'))
        parent = {}
        for c in pre:
            for ch in c.clades:
                parent[ch] = c
        st = {}
        for c in pre:
            if c not in parent:
                st[c] = rng.choice(NC, size=n, p=self.pi)
                continue
            t = c.branch_length * self.bl_scale
            P = (U * np.exp(lam * t)) @ Ui
            P = np.clip(P, 0, None); P /= P.sum(1, keepdims=True)
            cum = P.cumsum(1)
            u = rng.random(n)
            ps = st[parent[c]]
            st[c] = np.minimum((u[:, None] > cum[ps]).sum(1), NC - 1)
        return np.array([st[t] for t in self.tips])


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]; seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def tip_matrix(model, seqs, cols):
    """真实比对 → (n_tips, NC, len(cols)) 指示数组；缺失/歧义密码子按 HyPhy 处理成「全 1」
    （部分歧义如 A-G、ANN 也按全 1，略宽于 HyPhy 的部分解析，影响极小）"""
    out = np.zeros((len(model.tips), NC, len(cols)), dtype=bool)
    for k, t in enumerate(model.tips):
        s = seqs[t.name]
        for j, col in enumerate(cols):
            cod = s[3 * col: 3 * col + 3]
            if cod in CIDX:
                out[k, CIDX[cod], j] = True
            else:
                out[k, :, j] = True
    return out


def tip_matrix_from_idx(idx):
    n_tips, n = idx.shape
    out = np.zeros((n_tips, NC, n), dtype=bool)
    for k in range(n_tips):
        out[k, idx[k], np.arange(n)] = True
    return out
