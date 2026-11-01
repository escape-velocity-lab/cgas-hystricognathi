"""17b（任务 #17）: 4 个位点（人 431/479/495/498）在人 cGAS 实验结构上的位置。

主结构 6CT9（人 cGAS 157–522 + 16 bp dsDNA，2.26 Å，Zhou et al. 2018 Cell）；6CTA（同晶型 + ATP/Mg）定催化口袋；
4O68（无 DNA，2.44 Å）看构象是否依赖 DNA；AF-Q8N884（AlphaFold v6）只报 pLDDT 作参照。
人 cGAS–DNA 的 2:2 复合物在 6CT9 里由晶体对称生成，所以「到 DNA」分两种：
  同一不对称单元里的 DNA（= A 位点 DNA）、对称拷贝里的 DNA（含 B 位点 DNA，也可能只是晶格接触——表里分开报）。
输出：results/t17_structure_sites.csv、t17_structure_pairs.csv、t17_structure.txt、t17_structure.png
用法: python 13_t17_structure.py <proj_root>
"""
import sys, os, itertools
import numpy as np
import pandas as pd
import gemmi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
PDB = R('data', 'raw', 't17', 'pdb')
lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

SITES = {431: 'D431（裸鼹鼠 S463）', 479: 'K479（E511）', 495: 'L495（Y527）', 498: 'K498（T530）'}
CAT = [225, 227, 319]                         # 催化三联体 E225/D227/D319（核苷酸转移酶 Mg 配位）
MAXASA = dict(ALA=129, ARG=274, ASN=195, ASP=193, CYS=167, GLN=225, GLU=223, GLY=104, HIS=224, ILE=197, LEU=201, LYS=236,
              MET=224, PHE=240, PRO=159, SER=155, THR=172, TRP=285, TYR=263, VAL=174)   # Tien et al. 2013 理论最大值

def load(code):
    st = gemmi.read_structure(os.path.join(PDB, f'{code}.cif'))
    st.setup_entities()
    st.remove_hydrogens()
    return st

def is_dna(res):
    return res.name in ('DA', 'DC', 'DG', 'DT')

def secstruct(code):
    """mmCIF 的 _struct_conf（螺旋）/ _struct_sheet_range（折叠）→ {残基: 元件名}"""
    doc = gemmi.cif.read(os.path.join(PDB, f'{code}.cif')).sole_block()
    ss = {}
    t = doc.find('_struct_conf.', ['conf_type_id', 'beg_auth_asym_id', 'beg_auth_seq_id', 'end_auth_seq_id'])
    for r in t:
        if r[1] == 'A':
            for k in range(int(r[2]), int(r[3]) + 1):
                ss[k] = f'α螺旋 {r[2]}–{r[3]}'
    t = doc.find('_struct_sheet_range.', ['beg_auth_asym_id', 'beg_auth_seq_id', 'end_auth_seq_id'])
    for r in t:
        if r[0] == 'A':
            for k in range(int(r[1]), int(r[2]) + 1):
                ss[k] = f'β链 {r[1]}–{r[2]}'
    return ss

def element(ss, k):
    if k in ss:
        return ss[k]
    lo = k
    while lo - 1 not in ss and lo > 150:
        lo -= 1
    hi = k
    while hi + 1 not in ss and hi < 522:
        hi += 1
    return f'loop {lo}–{hi}'

def protein_chain(st):
    """蛋白链（只留聚合物残基：gemmi 把水 / Zn / ATP 也放在链 A 里，会混进接触列表）"""
    for ch in st[0]:
        if any(r.name == 'LYS' for r in ch):
            return [r for r in ch if r.het_flag == 'A']
    raise ValueError
def chain_name(st):
    return next(ch.name for ch in st[0] if any(r.name == 'LYS' for r in ch))

def min_dist(atoms_a, atoms_b):
    A = np.array([a.pos.tolist() for a in atoms_a]); B = np.array([b.pos.tolist() for b in atoms_b])
    D = np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1))
    return float(D.min())

def sym_contacts(st, chain_name, resnum, radius=12.0):
    """含晶体对称拷贝的最近邻：返回 {类别: 最短距离}；类别 = DNA_同单元 / DNA_对称拷贝 / 蛋白_对称拷贝"""
    model = st[0]
    ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
    out = {}
    ch = model[chain_name]
    res = next(r for r in ch if r.seqid.num == resnum)
    for a in res:
        for m in ns.find_atoms(a.pos, '\0', radius=radius):
            cra = m.to_cra(model)
            same_img = m.image_idx == 0
            if is_dna(cra.residue):
                key = 'DNA_同单元' if same_img else 'DNA_对称拷贝'
            elif cra.chain.name == chain_name and same_img:
                continue
            elif cra.residue.het_flag == 'H' or cra.residue.name == 'HOH':
                continue
            else:
                key = '蛋白_对称拷贝' if not same_img else '蛋白_其他链'
            # 对称拷贝的真实坐标
            img_pos = st.cell.find_nearest_pbc_position(a.pos, cra.atom.pos, m.image_idx)
            d = a.pos.dist(img_pos)
            if d < out.get(key, (99, ''))[0]:
                out[key] = (d, f'{cra.chain.name}/{cra.residue.name}{cra.residue.seqid.num}')
    return out

# ---------- 主结构 6CT9 ----------
st = load('6CT9'); stA = load('6CTA'); st0 = load('4O68')
pc = protein_chain(st)
PCN = chain_name(st)
P(f'6CT9：{st.spacegroup_hm}，蛋白链 {PCN}（{pc[0].seqid.num}–{pc[-1].seqid.num}，{len(pc)} 个残基）')
ss = secstruct('6CT9')
res = {r.seqid.num: r for r in pc}
dna_atoms = [a for ch in st[0] for r in ch if is_dna(r) for a in r]
zn = [a for ch in st[0] for r in ch if r.name == 'ZN' for a in r]
het_names = {r.name for ch in st[0] for r in ch if r.het_flag == 'H'}
zn_lig = sorted({r.seqid.num for r in pc for a in r if zn and a.pos.dist(zn[0].pos) < 2.8})
P('锌指配位残基（距 Zn < 2.8 Å）：', zn_lig)

# DNA 接触面（同单元 = A 位点；对称拷贝 = B 位点 / 晶格接触）
ns = gemmi.NeighborSearch(st[0], st.cell, 6).populate()
contact_A, contact_sym, by_img, mate_img = set(), set(), {}, {}
for r in pc:
    for a in r:
        for m in ns.find_atoms(a.pos, '\0', radius=4.5):
            cra = m.to_cra(st[0])
            img_pos = st.cell.find_nearest_pbc_position(a.pos, cra.atom.pos, m.image_idx)
            if a.pos.dist(img_pos) > 4.0:
                continue
            if is_dna(cra.residue):
                (contact_A if m.image_idx == 0 else contact_sym).add(r.seqid.num)
                if m.image_idx:
                    by_img.setdefault(m.image_idx, set()).add(r.seqid.num)
            elif m.image_idx and cra.residue.het_flag == 'A':
                mate_img.setdefault(m.image_idx, set()).add(r.seqid.num)
P(f'DNA 接触残基（重原子 ≤ 4 Å）：同单元 DNA（A 位点）{sorted(contact_A)}')
P(f'                               对称拷贝 DNA {sorted(contact_sym)}')
for i, v in sorted(by_img.items()):
    P(f'    对称拷贝 #{i} 的 DNA 接触残基：{sorted(v)}')
for i, v in sorted(mate_img.items()):
    P(f'    对称拷贝 #{i} 的 cGAS 接触残基（蛋白–蛋白 ≤ 4 Å）：{sorted(v)}')

# ATP（6CTA）
pcA = protein_chain(stA)
atp = [a for ch in stA[0] for r in ch if r.name == 'ATP' for a in r]
resA = {r.seqid.num: r for r in pcA}

# SASA（Biopython Shrake–Rupley；蛋白单链 / 蛋白 + 同单元 DNA）
bp = MMCIFParser(QUIET=True).get_structure('x', os.path.join(PDB, '6CT9.cif'))[0]
def sasa(model_chains):
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    s = Structure('s'); m = Model(0); s.add(m)
    for c in model_chains:
        c2 = c.copy()
        for r in list(c2):
            if r.id[0] == 'W':
                c2.detach_child(r.id)
        m.add(c2)
    ShrakeRupley().compute(m, level='R')
    return {r.id[1]: r.sasa for r in m[PCN] if r.id[0] == ' '}
dna_chains = [c for c in bp if any(r.get_resname() in ('DA', 'DC', 'DG', 'DT') for r in c)]
sasa_free = sasa([bp[PCN]])
sasa_cplx = sasa([bp[PCN]] + dna_chains)

rows = []
for k, lab in SITES.items():
    r = res[k]
    sc = [a for a in r if a.name not in ('N', 'C', 'O')]
    sym = sym_contacts(st, PCN, k)
    cat_d = min(min_dist(sc, [a for a in res[c] if a.name not in ('N', 'C', 'O', 'CA')]) for c in CAT)
    atp_d = min_dist([a for a in resA[k] if a.name not in ('N', 'C', 'O')], atp) if atp else np.nan
    dA = min((min_dist(sc, [a for a in res[c]]) for c in contact_A), default=np.nan)
    dS = min((min_dist(sc, [a for a in res[c]]) for c in contact_sym), default=np.nan)
    # 4O68（无 DNA）同一残基 CA 的局部构象：以 410–515 的 CA 叠合后的偏差
    rows.append(dict(
        位点=lab, PDB残基=f'{r.name}{k}', 二级结构=element(ss, k),
        到DNA_A位点_Å=round(min_dist(sc, dna_atoms), 1),
        到对称拷贝DNA_Å=round(sym.get('DNA_对称拷贝', (np.nan,))[0], 1), 对称DNA最近残基=sym.get('DNA_对称拷贝', (0, ''))[1],
        到对称拷贝cGAS_Å=round(sym.get('蛋白_对称拷贝', (np.nan,))[0], 1), 对称cGAS最近残基=sym.get('蛋白_对称拷贝', (0, ''))[1],
        到A位点接触残基_Å=round(dA, 1), 到对称DNA接触残基_Å=round(dS, 1),
        到催化残基E225_D227_D319_Å=round(cat_d, 1), 到ATP_6CTA_Å=round(atp_d, 1),
        到Zn_Å=round(min(a.pos.dist(zn[0].pos) for a in r), 1),
        SASA_蛋白单独_Å2=round(sasa_free[k], 1), 相对SASA=round(sasa_free[k] / MAXASA[r.name], 2),
        SASA_含DNA_Å2=round(sasa_cplx[k], 1), B因子_CA=round(r['CA'][0].b_iso, 1)))
df = pd.DataFrame(rows)
df.to_csv(R('results', 't17_structure_sites.csv'), index=False, encoding='utf-8-sig')
P('\n===== 每个位点\n' + df.to_string(index=False))

# 两两距离
pr = []
for a, b in itertools.combinations(SITES, 2):
    pr.append(dict(位点对=f'{a}–{b}', CA_CA_Å=round(res[a]['CA'][0].pos.dist(res[b]['CA'][0].pos), 1),
                   最近重原子_Å=round(min_dist(list(res[a]), list(res[b])), 1), 序列间隔=b - a))
pr = pd.DataFrame(pr)
pr.to_csv(R('results', 't17_structure_pairs.csv'), index=False, encoding='utf-8-sig')
P('\n===== 两两距离\n' + pr.to_string(index=False))

# 每个位点 6 Å 内的邻居（侧链朝向谁）
P('\n===== 6 Å 内邻居残基（6CT9，同链）')
for k in SITES:
    nb = sorted({x for x, r in res.items() if x != k and min_dist(list(res[k]), list(r)) <= 6.0})
    P(f'  {k}: {nb}')

# 其他结构核对：6CTA、4O68、AlphaFold
P('\n===== 其他结构里的同一残基')
for code, s in (('6CTA', stA), ('4O68', st0)):
    c = protein_chain(s)
    rr = {r.seqid.num: r for r in c}
    ss2 = secstruct(code)
    P(f'  {code}: ' + '；'.join(f'{k} {rr[k].name} {element(ss2, k)} B(CA)={rr[k]["CA"][0].b_iso:.0f}' if k in rr else f'{k} 缺失' for k in SITES))
# 4O68 vs 6CT9：C 端 400–522 CA 叠合
def ca(chain, lo, hi):
    return {r.seqid.num: np.array(r['CA'][0].pos.tolist()) for r in chain if lo <= r.seqid.num <= hi and r.find_atom('CA', '*')}
def kabsch(Pm, Qm):
    Pc, Qc = Pm - Pm.mean(0), Qm - Qm.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(U @ Vt))
    Rm = U @ np.diag([1, 1, d]) @ Vt
    return Pc @ Rm, Qc
a1, a0 = ca(pc, 400, 522), ca(protein_chain(st0), 400, 522)
common = sorted(set(a1) & set(a0))
X, Y = kabsch(np.array([a1[k] for k in common]), np.array([a0[k] for k in common]))
dev = np.sqrt(((X - Y) ** 2).sum(1))
P(f'  6CT9（带 DNA）vs 4O68（不带 DNA）C 端 400–522 CA 叠合 RMSD = {np.sqrt((dev ** 2).mean()):.2f} Å；'
  + '；'.join(f'{k} 偏差 {dev[common.index(k)]:.2f} Å' for k in SITES))
af = gemmi.read_structure(os.path.join(PDB, 'AF-Q8N884-F1-model_v6.cif'))
afr = {r.seqid.num: r for r in af[0][0]}
P('  AlphaFold v6 pLDDT：' + '；'.join(f'{k} {afr[k][0].b_iso:.1f}' for k in SITES))

open(R('results', 't17_structure.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')

# ---------- 图 ----------
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
cas = {k: np.array(r['CA'][0].pos.tolist()) for k, r in res.items() if r.find_atom('CA', '*')}
ks = sorted(cas)
C = np.array([cas[k] for k in ks])
dna_xyz = np.array([a.pos.tolist() for a in dna_atoms if a.name == 'P'])
atp_xyz = np.array([a.pos.tolist() for a in atp]) if atp else np.zeros((0, 3))   # 6CTA 与 6CT9 同晶型同坐标系
cen = C.mean(0)
# 投影面：取 DNA 轴和蛋白中心的方向，使 DNA 与 4 个位点都尽量摊开
_, _, Vt = np.linalg.svd(np.vstack([C, dna_xyz]) - cen)
proj = lambda X: (X - cen) @ Vt[:2].T
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 6.6), gridspec_kw=dict(width_ratios=[1.35, 1]))
p = proj(C)
lobe = np.array(ks) >= 400
for i in range(len(ks) - 1):
    if ks[i + 1] - ks[i] == 1:
        ax.plot(p[i:i + 2, 0], p[i:i + 2, 1], color='#4C78A8' if ks[i] >= 400 else '#B0B7C3', lw=1.6 if ks[i] >= 400 else 1.0, zorder=1)
d2 = proj(dna_xyz)
ax.scatter(d2[:, 0], d2[:, 1], s=26, color='#F2A541', edgecolor='none', label='DNA 磷原子（A 位点 DNA）', zorder=2)
if len(atp_xyz):
    a2 = proj(atp_xyz).mean(0)
    ax.scatter(*a2, s=160, marker='*', color='#2CA02C', label='ATP（6CTA，催化口袋）', zorder=3)
z2 = proj(np.array([zn[0].pos.tolist()]))[0]
ax.scatter(*z2, s=90, marker='D', color='#7F7F7F', label='Zn（锌拇指）', zorder=3)
for k, lab in SITES.items():
    q = proj(cas[k][None])[0]
    ax.scatter(*q, s=110, color='#D62728', edgecolor='white', lw=1.2, zorder=4)
    ax.annotate(lab, q, ha='right' if k == 479 else 'left', xytext={431: (10, 6), 479: (-14, -20), 495: (12, 2), 498: (12, -10)}[k], textcoords='offset points', fontsize=10, color='#D62728', weight='bold', zorder=5)
ax.plot([], [], color='#4C78A8', lw=1.6, label='C 端 400–522（含 4 个位点）')
ax.plot([], [], color='#B0B7C3', lw=1.0, label='157–399')
ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values():
    s.set_visible(False)
ax.legend(loc='upper left', fontsize=8.5, frameon=False)
ax.set_title('人 cGAS–DNA 晶体结构 6CT9（Cα 主链二维投影）', fontsize=11.5)
# 右：距离表
cols = ['二级结构', '到DNA_A位点_Å', '到对称拷贝DNA_Å', '到催化残基E225_D227_D319_Å', '到ATP_6CTA_Å', '相对SASA']
hdr = ['二级结构', '到 DNA\n(A 位点)', '到对称\n拷贝 DNA', '到催化\n残基', '到 ATP', '相对\nSASA']
cell = [['>12' if str(df.loc[i, c]) == 'nan' else str(df.loc[i, c]) for c in cols] for i in range(len(df))]
ax2.axis('off')
tb = ax2.table(cellText=cell, rowLabels=[s.split('（')[0] for s in SITES.values()], colLabels=hdr, loc='upper center', cellLoc='center')
tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1, 2.2)
for (ri, ci), cell in tb.get_celld().items():
    if ci == 0:
        cell.set_width(0.26)
    elif ci >= 0:
        cell.set_width(0.145)
for ri in range(len(cell_rowlabels := list(SITES.values())) + 1):
    c = tb.get_celld().get((ri, -1))
    if c is not None:
        c.set_width(0.17)
pt = '两两 Cα 距离（Å）：' + '，'.join(f'{r.位点对} {r.CA_CA_Å}' for r in pr.itertuples())
ax2.text(0.0, 0.25, pt, fontsize=9, wrap=True, transform=ax2.transAxes)
ax2.text(0.0, 0.08, '距离 = 侧链重原子到目标的最短距离；相对 SASA = 蛋白单链溶剂可及面积 / 理论最大值（>0.25 算暴露）。「>12」= 12 Å 内没有。\n'
                    '479 / 498 贴着「对称拷贝 DNA」（2.9–3.8 Å），但这只在 P 6₁22 这一个晶型里出现（6CT9 / 6CTA / 6EDC），\n'
                    '换一个晶型（6EDB，I222）就消失；而 site A 的 K173 / R176 / K407 在全部 4 套含 DNA 结构里都接触 DNA。\n'
                    '→ 判定为晶格接触，不是真实的 DNA 结合面（见 t17_lattice.txt）。',
         fontsize=8, color='#555', transform=ax2.transAxes)
ax2.set_title('距离与暴露度', fontsize=11.5)
plt.tight_layout()
plt.savefig(R('results', 't17_structure.png'), dpi=180)
P('图 → results/t17_structure.png')
