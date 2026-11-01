"""17a-4（任务 #17）: 不依赖选择模型的描述量——每个位点在树上的最少氨基酸变化次数（Fitch 简约）。

两个范围：
  (1) 豪猪亚目冠群 23 种（剪出子树）。只有 5 种有全长序列，其余 18 种只有最后一个外显子（人 406–522），
      所以「全蛋白平均」用两个口径：人 406–522 全部列（23 种都有数据，最公平）+ 全长列（缺失按「任意状态」处理，偏低）
  (2) 全部 275 种哺乳动物，全长 522 列。
同时报：同义变化不算（这里只数氨基酸变化）；每列的非缺失物种数。
输出：results/t17_descriptive.csv（每列）、t17_descriptive.txt（汇总）
用法: python 15_t17_descriptive.py <proj_root>
"""
import sys, os, json
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t15_lib as T

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

SITES = [431, 479, 495, 498]
aln = {r.id.replace('_', ' '): str(r.seq) for r in SeqIO.parse(R('results', 't17_codon', 'aln.fasta'), 'fasta')}
L = len(next(iter(aln.values()))) // 3
aa = {sp: [str(Seq(s[3 * i:3 * i + 3]).translate()) if '-' not in s[3 * i:3 * i + 3] else '-' for i in range(L)] for sp, s in aln.items()}
hy = json.load(open(R('results', 't17_codon', 'clades.json'), encoding='utf-8'))['hystricognathi']

root, miss = T.build(R('data', 'raw', 't13', 'tree', T.TREE_MAIN), list(aln))
assert not miss, miss
hy_root = T.mrca(root, hy)
assert sorted(T.tips_of(hy_root)) == sorted(hy)

rows = []
for i in range(L):
    pos = i + 1
    st = {sp: aa[sp][i] for sp in aln}
    _, c_all = T.fitch(root, st)
    _, c_hy = T.fitch(hy_root, {sp: st[sp] for sp in hy})
    n_hy = sum(st[sp] not in '-X' for sp in hy)
    n_all = sum(v not in '-X' for v in st.values())
    rows.append(dict(人位置=pos, 人残基=aa['Homo sapiens'][i], 豪猪亚目变化次数=c_hy, 豪猪亚目有数据物种=n_hy,
                     豪猪亚目状态=''.join(sorted({st[sp] for sp in hy} - {'-'})), 全树变化次数=c_all, 全树有数据物种=n_all))
df = pd.DataFrame(rows)
df.to_csv(R('results', 't17_descriptive.csv'), index=False, encoding='utf-8-sig')

def summ(sub, col, label):
    v = sub[col].values
    P(f'  {label}：{len(v)} 列，平均 {v.mean():.2f}，中位数 {np.median(v):.0f}，0 次变化的列占 {np.mean(v == 0):.0%}')
    for s in SITES:
        x = int(df.loc[df.人位置 == s, col].iloc[0])
        pct = np.mean(v < x) * 100, np.mean(v <= x) * 100
        P(f'      {s}: {x} 次（比它少的列 {pct[0]:.0f}%，≤ 它的列 {pct[1]:.0f}%）')
    four = df.loc[df.人位置.isin(SITES), col].values
    # 置换：随机抽 4 列，平均变化次数 ≥ 观测值的概率
    rng = np.random.default_rng(17)
    sims = np.array([rng.choice(v, 4, replace=False).mean() for _ in range(20000)])
    P(f'      4 个位点平均 {four.mean():.2f} 次；随机抽 4 列平均 ≥ 这个值的概率 = {np.mean(sims >= four.mean()):.3f}，≤ 的概率 = {np.mean(sims <= four.mean()):.3f}')

P('===== Fitch 最少氨基酸变化次数')
P('(1) 豪猪亚目冠群 23 种')
ex = df[(df.人位置 >= 406) & (df.人位置 <= 522)]
P(f'  人 406–522 段每列有数据物种数：中位数 {ex.豪猪亚目有数据物种.median():.0f} / 23（最少 {ex.豪猪亚目有数据物种.min()}）')
summ(ex, '豪猪亚目变化次数', '口径 A：人 406–522（23 种都有数据）')
summ(df[df.豪猪亚目有数据物种 >= 5], '豪猪亚目变化次数', '口径 B：全长中 ≥5 种有数据的列（多数只有 5 种，偏低，仅参考）')
P('(2) 全部哺乳动物')
summ(df[df.全树有数据物种 >= 100], '全树变化次数', '全长中 ≥100 种有数据的列')
summ(ex[ex.全树有数据物种 >= 100], '全树变化次数', '只看人 406–522')
P('\n4 个位点明细：')
P(df[df.人位置.isin(SITES)].to_string(index=False))
open(R('results', 't17_descriptive.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
