"""24（任务 #23d-a，来自 #21 自审「需补」第 2 条）: cGAS 三个已报道 DNA 结合位点（site A/B/C）的残基定义，
以及我们研究的 4 个位点（人 431/479/495/498）相对它们的位置。

为什么要做：审稿人（结构生化方向）会问「你们凭什么说这四个位点不属于任何已报道的 DNA 结合面」。
正文原来只说「不在已知功能位点上」，没有给出 site A/B/C 的定义和出处，等于让审稿人自己去信。
特别要查清一件事：site-C 的 KKH-loop 定义为 427-432，**我们的位点 431 就落在这个区间里**。

位点定义出处：
  site A / site B: Zhou et al. 2018 Cell 175:300（PDB 6CT9），2:2 复合物里每个单体的两个 DNA 结合面；
      代表残基 site A = K173/R176/K407/K411，site B = K347/K350。
  site C:          Xie et al. 2019 PNAS 116:11946（PDB 6EDB），人 cGAS 催化域上的第三个界面，由三段组成：
      alpha-region 261-286（接触 DNA 的碱性残基 Q264/K275/K279/K282/K285）、KRKR-loop 299-302（K299/R300/K301/R302）、
      KKH-loop 427-432（K427/K428/K432）。

site C 在 6EDB 的不对称单元内是接触不到 DNA 的，它由晶体对称生成的相邻 DNA 形成（这正是该文「促进相分离」
的论点所在），所以这里和 17_t17_lattice.py 一样，含对称拷贝一起算。
输出：results/t23_sitec.csv / t23_sitec.txt
用法: python 24_t23_sitec.py <proj_root>
"""
import sys, os
import pandas as pd
import gemmi

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
PDB = R('data', 'raw', 't17', 'pdb')
lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

SITE_DEF = {
    'site A':            [173, 176, 407, 411],
    'site B':            [347, 350],
    'site C alpha(261-286)': [264, 275, 279, 282, 285],
    'site C KRKR(299-302)':  [299, 300, 301, 302],
    'site C KKH(427-432)':   [427, 428, 432],
}
SITEC_SPANS = {'alpha-region': range(261, 287), 'KRKR-loop': range(299, 303), 'KKH-loop': range(427, 433)}
OURS = [431, 479, 495, 498]
CODES = ['6CT9', '6CTA', '6EDB', '6EDC']
is_dna = lambda r: r.name in ('DA', 'DC', 'DG', 'DT')

def min_dna_dist(code):
    """返回 {残基号: (残基名, ASU 内最近 DNA 距离, 含对称拷贝最近 DNA 距离)}"""
    st = gemmi.read_structure(os.path.join(PDB, f'{code}.cif'))
    st.setup_entities(); st.remove_hydrogens()
    m = st[0]
    pchain = next(ch for ch in m if sum(r.het_flag == 'A' for r in ch) > 100)
    res = {r.seqid.num: r for r in pchain if r.het_flag == 'A'}
    asu_dna = [a.pos for ch in m for r in ch if is_dna(r) for a in r]
    ns = gemmi.NeighborSearch(m, st.cell, 6).populate()
    out = {}
    for num, r in res.items():
        sc = [a for a in r if a.name not in ('N', 'C', 'O')]
        d_asu = min((a.pos.dist(x) for a in sc for x in asu_dna), default=None)
        d_sym = None
        for a in sc:
            for mark in ns.find_atoms(a.pos, '\0', radius=6.0):
                cra = mark.to_cra(m)
                if not is_dna(cra.residue): continue
                im = st.cell.find_nearest_pbc_image(a.pos, cra.atom.pos, mark.image_idx)
                d = im.dist()
                if d_sym is None or d < d_sym: d_sym = d
        out[num] = (r.name, d_asu, d_sym)
    return st.spacegroup_hm, out

rows = []
for code in CODES:
    sg, dd = min_dna_dist(code)
    P('=' * 78)
    P(f'{code}  空间群 {sg}')
    P('  已报道位点的代表残基（侧链到 DNA 的最近距离，单位 A；sym = 含晶体对称拷贝）')
    for site, nums in SITE_DEF.items():
        cells = []
        for n in nums:
            if n not in dd: cells.append(f'{n}:缺失'); continue
            nm, da, ds = dd[n]
            cells.append('%s%d asu=%s sym=%s' % (nm, n, f'{da:.1f}' if da else 'na', f'{ds:.1f}' if ds else '>6'))
            rows.append(dict(PDB=code, 空间群=sg, 分类=site, 残基号=n, 残基=nm,
                             ASU内到DNA=round(da, 2) if da else None, 含对称到DNA=round(ds, 2) if ds else None))
        P(f'    {site:<24} ' + ' | '.join(cells))
    P('  我们研究的 4 个位点')
    for n in OURS:
        if n not in dd: P(f'    {n}: 缺失'); continue
        nm, da, ds = dd[n]
        tag = '  <-- 落在 site C 的 KKH-loop(427-432) 区间内' if 427 <= n <= 432 else ''
        P('    %s%d  asu=%s  sym=%s%s' % (nm, n, f'{da:.1f}' if da else 'na', f'{ds:.1f}' if ds else '>6', tag))
        rows.append(dict(PDB=code, 空间群=sg, 分类='本文研究的位点', 残基号=n, 残基=nm,
                         ASU内到DNA=round(da, 2) if da else None, 含对称到DNA=round(ds, 2) if ds else None))
    # KKH-loop 全段逐残基，看 431 在这段里是什么角色
    P('  site C KKH-loop 全段逐残基（看 431 在这段里的角色）')
    seg = []
    for n in SITEC_SPANS['KKH-loop']:
        if n not in dd: seg.append(f'{n}:缺失'); continue
        nm, da, ds = dd[n]
        seg.append('%s%d sym=%s' % (gemmi.find_tabulated_residue(nm).one_letter_code.upper(), n, f'{ds:.1f}' if ds else '>6'))
    P('    ' + '  '.join(seg))

pd.DataFrame(rows).to_csv(R('results', 't23_sitec.csv'), index=False, encoding='utf-8-sig')
open(R('results', 't23_sitec.txt'), 'w', encoding='utf-8').write('\n'.join(lines))
P('')
P('写出 results/t23_sitec.csv 和 t23_sitec.txt')
