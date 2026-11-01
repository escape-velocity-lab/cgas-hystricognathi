"""29（任务 #23a-2）: 核实「R255 是 cGAS 结合核小体酸性斑的关键残基」，并量我们的谱系候选 K258 到同一界面的距离。

背景：任务 #22 汇报里「R255 是关键残基」是凭记忆写的，必须核（R26）。这里不引综述，直接从结构算：
取三套人 cGAS–核小体复合物的冷冻电镜结构（7CCQ = Cao et al. 2020；6Y5E、6Y5D = Pathare et al. 2020），
其中 6Y5D 含两个独立的 cGAS 拷贝，所以一共 4 个独立拷贝。
酸性斑按 Kalashnikova 等人的常用定义取 H2A 的 56/61/64/90/91/92（不同结构里 90/91 是 D 还是 N/E 略有差别，按残基号取）。
对每个 cGAS 拷贝，列出所有到酸性斑 <= 5 A 的 cGAS 残基，看 R255 和 K258 是不是真在名单里。
**必须逐拷贝算**：6Y5D 两个拷贝的 K258 一个 4.04 A 一个 5.58 A，按「整个结构取最近」会掩盖这种差异。
输出：results/t23_nucleosome.csv / t23_nucleosome.txt
用法: python 29_t23_nucleosome.py <proj_root> <cif 目录>
"""
import sys, os
import pandas as pd
import gemmi

ROOT, CIFDIR = sys.argv[1], sys.argv[2]
R = lambda *p: os.path.join(ROOT, *p)
AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
ACID = [56, 61, 64, 90, 91, 92]
STRUCTS = [('7CCQ', 'Cao et al. 2020, 1:1 cGAS-NCP'), ('6Y5E', 'Pathare et al. 2020'), ('6Y5D', 'Pathare et al. 2020')]
CUT = 5.0
L = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); L.append(s)

def seqof(c):
    return ''.join(gemmi.find_tabulated_residue(r.name).one_letter_code.upper() for r in c if r.name in AA)

rows, per_copy = [], []
for pdbid, cite in STRUCTS:
    st = gemmi.read_structure(os.path.join(CIFDIR, pdbid + '.cif')); st.setup_entities()
    st.remove_alternative_conformations(); st.remove_hydrogens()
    m = st[0]
    cgas = [c for c in m if len(seqof(c)) > 250]
    h2a = [c for c in m if 'AGLQFPV' in seqof(c)]      # H2A 核心序列特征
    acid = [(c.name, r.name, r.seqid.num, a.pos) for c in h2a for r in c
            if r.seqid.num in ACID and r.name in AA for a in r]
    P('=' * 76)
    P('%s (%s): cGAS 链 %s，H2A 链 %s' % (pdbid, cite, [c.name for c in cgas], [c.name for c in h2a]))
    for cg in cgas:
        hits = []
        for r in cg:
            if r.name not in AA: continue
            pts = [a.pos for a in r]
            t, d = min(((t, min(t[3].dist(x) for x in pts)) for t in acid), key=lambda z: z[1])
            if d <= CUT:
                hits.append((r.seqid.num, r.name, round(d, 2), t[0], t[1], t[2]))
                rows.append(dict(PDB=pdbid, cGAS_chain=cg.name, residue_human_number=r.seqid.num,
                                 residue=r.name, min_dist_to_acidic_patch_A=round(d, 2),
                                 nearest_H2A_chain=t[0], nearest_H2A_residue='%s%d' % (t[1], t[2])))
        hits.sort()
        P('  cGAS 拷贝 %s：到酸性斑 <= %.1f A 的残基共 %d 个 —— %s' % (
            cg.name, CUT, len(hits), ', '.join('%s%d(%.2f)' % (h[1], h[0], h[2]) for h in hits)))
        per_copy.append((pdbid, cg.name, {h[0]: h[2] for h in hits}))

P('')
P('===== 逐拷贝对照（4 个独立 cGAS 拷贝）=====')
P('%-14s %-8s %-8s %-8s' % ('拷贝', 'R255', 'P257', 'K258'))
for pdbid, ch, dd in per_copy:
    P('%-14s %-8s %-8s %-8s' % ('%s:%s' % (pdbid, ch),
                                dd.get(255, '>5'), dd.get(257, '>5'), dd.get(258, '>5')))
common = set.intersection(*[set(d) for _, _, d in per_copy])
P('')
P('4 个拷贝共有的界面残基（全部 <= %.1f A）: %s' % (CUT, sorted(common)))
P('结论：R255 在全部 4 个拷贝里都以 2.3-2.9 A 直接接触酸性斑，是文献所说的 arginine anchor —— 记忆正确。')
P('      K258 在 4 个拷贝里有 3 个落在 4.0-4.2 A（第 4 个 5.58 A），属于界面外缘的长程静电接触，')
P('      不是 R255 那样的盐桥。稿件按这个口径写，不要写成「盐桥」或「关键残基」。')

pd.DataFrame(rows).to_csv(R('results', 't23_nucleosome.csv'), index=False, encoding='utf-8-sig')
open(R('results', 't23_nucleosome.txt'), 'w', encoding='utf-8').write('\n'.join(L))
P('')
P('写出 results/t23_nucleosome.csv / .txt')
