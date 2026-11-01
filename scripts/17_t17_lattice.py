"""17b 追加（任务 #17）: 479/498 接触相邻 DNA——是真界面还是晶格假象？

6CT9 里 K479、K498 的侧链离「晶体对称生成的相邻 DNA」只有 2.9–3.8 Å，但 PDB 的生物组装（assembly1）
只有 1 蛋白 + 1 DNA 双链，没把这个接触算进去。区分真界面 / 晶格接触的标准做法：
**看多个独立晶型（不同空间群 / 不同晶胞）里同一个接触是否重复出现。**
这里扫 6 套含 DNA 的人 cGAS 结构（6CT9、6CTA、6EDB、6EDC + 两套无 DNA 的 4O68、4KM5 作阴性对照：
无 DNA 的晶体里这些残基接触到什么）。
对每套结构：对每个位点，报 (a) 同一不对称单元内最近的 DNA、(b) 含对称拷贝后最近的 DNA 及所属对称操作。
输出：results/t17_lattice.csv / t17_lattice.txt
用法: python 17_t17_lattice.py <proj_root>
"""
import sys, os
import numpy as np
import pandas as pd
import gemmi

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
PDB = R('data', 'raw', 't17', 'pdb')
lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

SITES = [431, 479, 495, 498]
# 对照残基：site A 的已知 DNA 结合残基（K173/R176/K407）+ 一组随机表面残基
REF_A = [173, 176, 407]
CODES = ['6CT9', '6CTA', '6EDB', '6EDC', '4O68', '4KM5']
is_dna = lambda r: r.name in ('DA', 'DC', 'DG', 'DT')

rows = []
for code in CODES:
    st = gemmi.read_structure(os.path.join(PDB, f'{code}.cif'))
    st.setup_entities(); st.remove_hydrogens()
    m = st[0]
    pchain = next(ch for ch in m if sum(r.het_flag == 'A' for r in ch) > 100)
    res = {r.seqid.num: r for r in pchain if r.het_flag == 'A'}
    has_dna = any(is_dna(r) for ch in m for r in ch)
    ns = gemmi.NeighborSearch(m, st.cell, 6).populate()
    for k in SITES + REF_A:
        if k not in res:
            rows.append(dict(PDB=code, 空间群=st.spacegroup_hm, 位点=k, 残基='缺失')); continue
        r = res[k]
        sc = [a for a in r if a.name not in ('N', 'C', 'O')]
        best = {}
        for a in sc:
            for mk in ns.find_atoms(a.pos, '\0', radius=6.0):
                cra = mk.to_cra(m)
                p = st.cell.find_nearest_pbc_position(a.pos, cra.atom.pos, mk.image_idx)
                d = a.pos.dist(p)
                if d > 6.0:
                    continue
                if is_dna(cra.residue):
                    key = 'DNA同单元' if mk.image_idx == 0 else 'DNA对称'
                elif mk.image_idx and cra.residue.het_flag == 'A':
                    key = '蛋白对称'
                else:
                    continue
                if d < best.get(key, (99, ''))[0]:
                    best[key] = (d, f'{cra.chain.name}/{cra.residue.name}{cra.residue.seqid.num}#{mk.image_idx}')
        rows.append(dict(PDB=code, 空间群=st.spacegroup_hm, 含DNA=has_dna, 位点=k, 残基=r.name,
                         到同单元DNA_Å=round(best.get('DNA同单元', (np.nan,))[0], 1),
                         到对称DNA_Å=round(best.get('DNA对称', (np.nan,))[0], 1), 对称DNA=best.get('DNA对称', (0, ''))[1],
                         到对称蛋白_Å=round(best.get('蛋白对称', (np.nan,))[0], 1)))
df = pd.DataFrame(rows)
df.to_csv(R('results', 't17_lattice.csv'), index=False, encoding='utf-8-sig')

P('===== 4 个位点 + site A 对照残基（173/176/407）在 6 套人 cGAS 晶体里到 DNA 的最短距离（侧链重原子，截断 6 Å）')
P('（「对称」= 晶体对称生成的相邻分子。真界面会在多个独立晶型里重复；晶格接触通常只在一个晶型里出现。）\n')
for code, g in df.groupby('PDB', sort=False):
    sg = g.空间群.iloc[0]
    hd = g.含DNA.iloc[0] if '含DNA' in g else '?'
    P(f'{code}（{sg}，{"含 DNA" if hd else "无 DNA"}）')
    for _, r in g.iterrows():
        if r.残基 == '缺失':
            P(f'    {r.位点}: 缺失'); continue
        tag = ' ← site A 已知 DNA 结合残基' if r.位点 in REF_A else ''
        P(f'    {r.位点}{r.残基}: 同单元 DNA {r.到同单元DNA_Å}，对称拷贝 DNA {r.到对称DNA_Å}（{r.对称DNA}），对称蛋白 {r.到对称蛋白_Å}{tag}')
    P('')

P('===== 判读')
dna_struct = df[(df.含DNA == True) & (df.残基 != '缺失')]
for k in SITES + REF_A:
    g = dna_struct[dna_struct.位点 == k]
    n_same = int((g.到同单元DNA_Å <= 4.5).sum())
    n_sym = int((g.到对称DNA_Å <= 4.5).sum())
    sgs = sorted(set(g.loc[g.到对称DNA_Å <= 4.5, '空间群']))
    P(f'  {k}: {len(g)} 套含 DNA 结构里，同单元 DNA ≤4.5 Å 的 {n_same} 套；对称 DNA ≤4.5 Å 的 {n_sym} 套'
      + (f'（空间群 {sgs}）' if sgs else '') + ('  ← site A 对照' if k in REF_A else ''))
open(R('results', 't17_lattice.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
