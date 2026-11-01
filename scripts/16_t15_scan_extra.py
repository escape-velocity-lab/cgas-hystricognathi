"""16（任务 #15a 补充）: 在 14 的投影比对上数两种对照：
(1) 放宽档：非洲鼹鼠科 3 种共享 X、其余豪猪亚目恰有 1 种（非豚鼠属）也是 X
(2) 只有裸鼹鼠自己有（两种 Fukomys 和其余豪猪亚目都没有）的状态——裸鼹鼠末端枝上的替换
用法: python 16_t15_scan_extra.py <proj_root>
"""
import sys, os
import pandas as pd
ROOT = sys.argv[1]
A = pd.read_csv(os.path.join(ROOT, 'results', 't15_aln_nmrcoord.tsv'), sep='\t', index_col=0).fillna('-')
c = pd.read_csv(os.path.join(ROOT, 'results', 't15_lineage_candidates.csv'))
fam = pd.read_csv(os.path.join(ROOT, 'results', 't13_species_sites.csv')).set_index('物种')
HF = {'Bathyergidae','Caviidae','Chinchillidae','Octodontidae','Erethizontidae','Thryonomyidae','Petromuridae','Echimyidae','Hystricidae','Ctenomyidae','Dinomyidae','Dasyproctidae','Cuniculidae'}
H = [s for s in A.index if s in fam.index and fam.loc[s,'科'] in HF]
B = ['Heterocephalus glaber','Fukomys damarensis','Fukomys anselli']
O = [s for s in H if s not in B]
nmr = A.loc['Heterocephalus glaber']
relax, nmr_only = [], []
for p in A.columns:
    X = nmr[p]
    if X in '-?': continue
    oth = {s: A.loc[s,p] for s in O if A.loc[s,p] not in '-?'}
    if A.loc['Cavia porcellus',p] in '-?': continue
    k = [s for s,v in oth.items() if v == X]
    fd, fa = A.loc['Fukomys damarensis',p], A.loc['Fukomys anselli',p]
    if fd == X and fa in (X,'?') and len(k) == 1 and not k[0].startswith('Cavia'):
        relax.append((int(p), X, k[0]))
    if fd not in '-?' and fd != X and fa != X and not k:
        nmr_only.append((int(p), X))
print('严格候选', len(c))
print('放宽档（其余豪猪亚目恰 1 种也有）', len(relax), relax)
print('裸鼹鼠独有（Fukomys 与其余豪猪亚目都没有，家豚鼠与达马拉兰鼹鼠可读）', len(nmr_only))
print(nmr_only)
