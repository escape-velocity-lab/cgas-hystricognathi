"""25（任务 #23d-c，来自 #21 自审「需补」第 3 条）: 275 个物种里，哪些只出现在 Upham 的 completed 树里、
在 DNA-only 树里没有——也就是说它们在树上的位置是按分类学插补的，不是分子数据定的。

为什么要做：审稿人会盯最有力的那几个反例。*Cavia tschudii*（山地豚鼠，LQ 0.48）携带 4 个位点里的
裸鼹鼠型残基，是「这些残基和长寿无关」最直接的证据之一。如果它在树上的位置是插补的，审稿人会追问。
（其实残基本身是直接测出来的、与建树无关，但需要把话说清楚。）
输出：results/t23_dnaonly.csv / t23_dnaonly.txt
用法: python 25_t23_dnaonly_check.py <proj_root>
"""
import sys, os, re
import pandas as pd

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
TRE = R('data', 'raw', 't13', 'tree')
COMP = os.path.join(TRE, 'MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')
DNA = os.path.join(TRE, 'MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre')

def tips(path):
    t = open(path, encoding='utf-8').read()
    return {m.rsplit('_', 2)[0]: m for m in re.findall(r'([A-Z][a-z]+_[a-z]+_[A-Z]+_[A-Z]+)', t)} or \
           {m.rsplit('_', 1)[0]: m for m in re.findall(r'([A-Z][a-z]+_[a-z]+_[A-Z]+)', t)}

comp, dna = tips(COMP), tips(DNA)
print('completed 树 %d 个叶子，DNA-only 树 %d 个叶子' % (len(comp), len(dna)))

# Upham 2019 用的是旧分类 / 近缘种代替，直接从 10_t13_asr.py 重用同一张表，不重写
#（坑：不用这张表会把 17 个已上树的物种误判为「不在树上」）
src10 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '10_t13_asr.py'), encoding='utf-8').read()
ns = {}
exec(src10[src10.index('TREE_SYN = {'):src10.index('class N:')], ns)
TREE_SYN = {k.replace(' ', '_'): v.replace(' ', '_') for k, v in ns['TREE_SYN'].items()}
print('同物异名映射 %d 条（复用 10_t13_asr.py）' % len(TREE_SYN))

def lookup(k, d):
    return k in d or TREE_SYN.get(k) in d

src = pd.read_csv(R('results', 't17_codon', 'sources.csv'), encoding='utf-8-sig')
src.columns = ['species', 'order', 'family', 'source', 'cov', 'hystri']
rows = []
for _, r in src.iterrows():
    k = r['species'].replace(' ', '_')
    rows.append(dict(物种=r['species'], 目=r['order'], 科=r['family'], 豪猪亚目=bool(r['hystri']),
                     在completed树=lookup(k, comp), 在DNAonly树=lookup(k, dna),
                     树上用名=TREE_SYN.get(k, k).replace('_', ' '),
                     位置来源='分子数据' if lookup(k, dna) else ('分类学插补' if lookup(k, comp) else '不在树上')))
d = pd.DataFrame(rows)
d.to_csv(R('results', 't23_dnaonly.csv'), index=False, encoding='utf-8-sig')

L = []
P = lambda *a: (print(*a, flush=True), L.append(' '.join(str(x) for x in a)))
P('275 个比对物种里：')
P(d['位置来源'].value_counts().to_string())
P('')
P('【豪猪亚目】逐个列出：')
P(d[d['豪猪亚目']][['物种', '科', '位置来源']].to_string(index=False))
P('')
imp = d[d['位置来源'] == '分类学插补']
P('全部 %d 个位置为分类学插补的物种：' % len(imp))
P(imp[['物种', '目', '科', '豪猪亚目']].to_string(index=False) if len(imp) else '  （无）')
open(R('results', 't23_dnaonly.txt'), 'w', encoding='utf-8').write('\n'.join(L))
