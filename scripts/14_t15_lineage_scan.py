"""14（任务 #15a）: 找「非洲鼹鼠科有、豚鼠属和其余豪猪亚目都没有」的 cGAS 氨基酸状态，并做祖先重建确认出现在非洲鼹鼠科主干上。

坐标系：裸鼹鼠 A0AAX6RS70（Chen 2025 用的那条，554 aa；4 个位点 = S463/E511/Y527/T530），另给人 Q8N884 编号。
数据（每个物种一条序列 → 投影到裸鼹鼠坐标，裸鼹鼠没有的插入列丢掉）：
  - #11 的 270 条 RefSeq（只取「4位点可用=是」「物种代表=是」，同一条异构体）：Biopython 全局比对到裸鼹鼠（参数同 #11）
  - GenBank 美洲豪猪 KAM6166271.1、梳趾鼠 KAM5256511.1：同上
  - 18 个无注释豪猪亚目组装：12 号脚本的全长 tblastn（查询 = 裸鼹鼠全长），HSP 的 query 坐标就是裸鼹鼠坐标
    过滤：同一性 ≥ 0.35 且长度 ≥ 25（cGAS N 端本来就分化大）；先取 bitscore 最高的 contig 上全部 HSP，
    其他 contig 的 HSP 只填补尚未覆盖的位置；MB21D2 一类旁系同源（同一性约 0.23）自然被滤掉
扫描条件（严格）：位置 p 上 裸鼹鼠 = 达马拉兰鼹鼠 = 安氏鼹鼠 = X，且其余豪猪亚目（含豚鼠属）凡是能读到的都 ≠ X，且家豚鼠必须能读到。
  放宽档：安氏鼹鼠读不到（?）的也列出，单独标注。
祖先重建：每个候选位置在全部物种（Upham 主树）上做 Mk-ER 多态边际 + Fitch；
  看 P(X) 在 非洲鼹鼠科祖先 / 其父节点（非洲鼹鼠+蔗鼠+岩鼠）/ 豪猪亚目祖先 的值。
输出: results/t15_lineage_candidates.csv、t15_lineage_scan.txt、t15_aln_nmrcoord.tsv（每物种一行的投影比对）
用法: python 14_t15_lineage_scan.py <proj_root>
"""
import sys, os, glob, json
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t15_lib as T

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

NMRSEQ = str(next(SeqIO.parse(R('data', 'raw', 'uniprot', 'A0AAX6RS70.fasta'), 'fasta')).seq)
HUMSEQ = str(next(SeqIO.parse(R('data', 'raw', 'uniprot', 'Q8N884.fasta'), 'fasta')).seq)
L = len(NMRSEQ)
assert (NMRSEQ[462], NMRSEQ[510], NMRSEQ[526], NMRSEQ[529]) == tuple('SEYT')

al = Align.PairwiseAligner()
al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score = -10, -0.5
al.end_insertion_score = al.end_deletion_score = 0

def project(seq):
    """seq 全局比对到裸鼹鼠 → 长度 L 的列表（裸鼹鼠每个位置上的残基，缺失 '-'）"""
    a = al.align(NMRSEQ, seq)[0]
    out = ['-'] * L
    for (t0, t1), (q0, q1) in zip(*a.aligned):
        for k in range(t1 - t0):
            out[t0 + k] = seq[q0 + k]
    return out

# 裸鼹鼠 → 人 编号
_h = project(HUMSEQ)
a = al.align(NMRSEQ, HUMSEQ)[0]
N2H = {}
for (t0, t1), (q0, q1) in zip(*a.aligned):
    for k in range(t1 - t0):
        N2H[t0 + k + 1] = q0 + k + 1
assert [N2H[p] for p in (463, 511, 527, 530)] == [431, 479, 495, 498], [N2H.get(p) for p in (463, 511, 527, 530)]

# ---------- 序列 ----------
aln, meta = {}, {}
v1 = pd.read_csv(R('results', 'site_table_v1.csv'))
v1 = v1[(v1['4位点可用'] == '是') & (v1['物种代表'] == '是')]
prot = {r.id: str(r.seq) for r in SeqIO.parse(R('data', 'raw', 'ncbi', 'pkg', 'ncbi_dataset', 'data', 'protein.faa'), 'fasta')}
for _, r in v1.iterrows():
    sp = ' '.join(r.物种.split()[:2])
    aln[sp] = project(prot[r.RefSeq蛋白].replace('*', 'X'))
    meta[sp] = dict(来源='RefSeq ' + r.RefSeq蛋白, 科=r.科, 目=r.目)
gb = {r.id: str(r.seq) for r in SeqIO.parse(R('data', 'raw', 't13', 'genbank_extra.faa'), 'fasta')}
for acc, sp, fam in (('KAM6166271.1', 'Erethizon dorsatum', 'Erethizontidae'), ('KAM5256511.1', 'Ctenodactylus gundi', 'Ctenodactylidae')):
    aln[sp] = project(gb[acc]); meta[sp] = dict(来源='GenBank ' + acc, 科=fam, 目='Rodentia')

src12 = open(R('src', '12_t15_tblastn_full.py'), encoding='utf-8').read()
import re
JOBS = dict((db, sp) for sp, db in re.findall(r"\('([A-Z][a-z]+ [a-z]+)', '([A-Z0-9]+)'\)", src12))
TB_FAM = {'Hystrix cristata': 'Hystricidae', 'Hystrix brachyura': 'Hystricidae', 'Thryonomys swinderianus': 'Thryonomyidae',
          'Petromus typicus': 'Petromuridae', 'Fukomys anselli': 'Bathyergidae', 'Cavia tschudii': 'Caviidae', 'Cavia aperea': 'Caviidae',
          'Hydrochoerus hydrochaeris': 'Caviidae', 'Dolichotis patagonum': 'Caviidae', 'Dasyprocta punctata': 'Dasyproctidae',
          'Cuniculus paca': 'Cuniculidae', 'Dinomys branickii': 'Dinomyidae', 'Ctenomys sociabilis': 'Ctenomyidae',
          'Myocastor coypus': 'Echimyidae', 'Capromys pilorides': 'Echimyidae', 'Octomys mimax': 'Octodontidae',
          'Tympanoctomys barrerae': 'Octodontidae', 'Coendou prehensilis': 'Erethizontidae'}
cov_rows = []
for f in sorted(glob.glob(R('data', 'raw', 't15', 'tblastn_full', '*.xml'))):
    db = os.path.basename(f)[:-4]
    sp = JOBS[db]
    root = ET.parse(f).getroot()
    hs = []
    for hit in root.iter('Hit'):
        for h in hit.iter('Hsp'):
            idn, alen = int(h.findtext('Hsp_identity')), int(h.findtext('Hsp_align-len'))
            if idn / alen >= 0.35 and alen >= 25:
                hs.append(dict(acc=hit.findtext('Hit_accession'), bits=float(h.findtext('Hsp_bit-score')), qf=int(h.findtext('Hsp_query-from')),
                               qt=int(h.findtext('Hsp_query-to')), q=h.findtext('Hsp_qseq'), s=h.findtext('Hsp_hseq'), id=idn / alen))
    out = ['?'] * L
    if hs:
        best = max(hs, key=lambda x: x['bits'])['acc']
        hs.sort(key=lambda x: (x['acc'] != best, -x['bits']))
        for h in hs:
            qp = h['qf']
            for x, y in zip(h['q'], h['s']):
                if x == '-':
                    continue
                if out[qp - 1] == '?':
                    out[qp - 1] = y if y != '*' else 'X'
                qp += 1
        ncont = len({h['acc'] for h in hs})
    else:
        ncont = 0
    aln[sp] = out
    meta[sp] = dict(来源=f'tblastn 全长 {db}', 科=TB_FAM[sp], 目='Rodentia')
    covd = sum(c not in '?' for c in out)
    cov_rows.append(dict(物种=sp, WGS=db, 覆盖裸鼹鼠位置=covd, 覆盖率=round(covd / L, 3), contig数=ncont,
                         HSP=';'.join(f"{h['acc']}:{h['qf']}-{h['qt']}(id{h['id']:.2f})" for h in hs)))
cov = pd.DataFrame(cov_rows)
P('===== tblastn 全长覆盖（裸鼹鼠 554 位）')
P(cov[['物种', 'WGS', '覆盖裸鼹鼠位置', '覆盖率', 'contig数']].to_string(index=False))
missing_tb = set(TB_FAM) - set(cov.物种)
if missing_tb:
    P('！以下 tblastn 结果缺失（未纳入）：', ', '.join(sorted(missing_tb)))

HYST = {sp for sp, m in meta.items() if m['科'] in set(TB_FAM.values()) | {'Bathyergidae', 'Chinchillidae', 'Octodontidae', 'Erethizontidae', 'Caviidae'}}
BATH = ['Heterocephalus glaber', 'Fukomys damarensis', 'Fukomys anselli']
OTHER_H = sorted(HYST - set(BATH))
P(f'\n豪猪亚目 {len(HYST)} 种；非洲鼹鼠科 3；其余豪猪亚目 {len(OTHER_H)}：{", ".join(OTHER_H)}')

# 投影比对存盘
pd.DataFrame({sp: list(v) for sp, v in aln.items()}, index=range(1, L + 1)).T.to_csv(R('results', 't15_aln_nmrcoord.tsv'), sep='\t')

# ---------- 扫描 ----------
def win_id(sp, p, w=10):
    """物种 sp 在裸鼹鼠位置 p 两侧 ±w 的同一性（只计双方都有残基的位置）"""
    a_ = aln[sp]
    xs = [k for k in range(max(1, p - w), min(L, p + w) + 1) if k != p and a_[k - 1] not in '-?']
    return sum(a_[k - 1] == NMRSEQ[k - 1] for k in xs) / max(1, len(xs)), len(xs)

nonhyst = [sp for sp in aln if sp not in HYST]
cands = []
for p in range(1, L + 1):
    X = NMRSEQ[p - 1]
    fd, fa = aln['Fukomys damarensis'][p - 1], aln['Fukomys anselli'][p - 1]
    if fd != X or fa not in (X, '?'):
        continue
    cp = aln['Cavia porcellus'][p - 1]
    if cp in '-?':
        continue
    oth = {sp: aln[sp][p - 1] for sp in OTHER_H}
    readable = {sp: v for sp, v in oth.items() if v not in '-?'}
    if any(v == X for v in readable.values()):
        continue
    nh = [aln[sp][p - 1] for sp in nonhyst if aln[sp][p - 1] not in '-?']
    cands.append(dict(裸鼹鼠位置=p, 人位置=N2H.get(p, ''), 人残基=HUMSEQ[N2H[p] - 1] if p in N2H else '-', 非洲鼹鼠科状态=X,
                      安氏鼹鼠=fa, 家豚鼠=cp, 其余豪猪亚目可读=f'{len(readable)}/{len(OTHER_H)}',
                      其余豪猪亚目残基=''.join(sorted(set(readable.values()))),
                      其余豪猪亚目明细=' '.join(f'{sp.split()[0][:4]}.{sp.split()[1][:3]}={v}' for sp, v in oth.items()),
                      梳趾鼠=aln['Ctenodactylus gundi'][p - 1],
                      豪猪亚目外哺乳动物有X=f'{sum(v == X for v in nh)}/{len(nh)}',
                      **{f'窗口同一性_{k}': round(win_id(sp, p)[0], 2) for k, sp in (('达马', 'Fukomys damarensis'), ('豚鼠', 'Cavia porcellus'))}))
cd = pd.DataFrame(cands)
P(f'\n===== 候选位置（严格 + 安氏鼹鼠读不到的放宽档）：{len(cd)} 个')

# ---------- 结构域注释 ----------
uj = json.load(open(R('data', 'raw', 't15', 'Q8N884.json'), encoding='utf-8'))
FEATS = [(f['type'], f['location']['start']['value'], f['location']['end']['value'], f.get('description', '') or (f.get('ligand') or {}).get('name', ''))
         for f in uj['features'] if f['type'] in ('Region', 'Motif', 'Binding site', 'Site', 'Modified residue', 'Cross-link', 'Zinc finger')]
def domain(h):
    if not h:
        return '裸鼹鼠特有插入（人无对应）'
    if h <= 160:
        return 'N 端无序区（1–160）'
    if h <= 330:
        return '催化结构域 / NTase 核心（161–330）'
    return 'C 端结构域（331–522，含锌指 390–405）'
def feats(h):
    if not h:
        return ''
    out = []
    for t, s, e, d in FEATS:
        if s <= h <= e and not (t == 'Region' and d == 'Disordered'):
            out.append(f'{t}:{d}' if s == e else f'{t} {s}-{e}:{d}')
        elif t in ('Binding site', 'Site', 'Modified residue', 'Cross-link') and s == e and abs(s - h) <= 2:
            out.append(f'邻近({s}) {t}:{d[:40]}')
    return '; '.join(out)
if len(cd):
    cd['结构域'] = [domain(h if h != '' else None) for h in cd.人位置]
    cd['UniProt注释(人)'] = [feats(h if h != '' else None) for h in cd.人位置]

# ---------- 祖先重建 ----------
full, miss = T.build(R('data', 'raw', 't13', 'tree', T.TREE_MAIN), list(aln))
P(f'树：{len(T.tips_of(full))} 种上树；缺 {miss}')
NB = T.mrca(full, BATH)
NPAR = NB.up
NH = T.mrca(full, ['Cavia porcellus', 'Heterocephalus glaber', 'Hystrix cristata'])
P('非洲鼹鼠科父节点下含：', ', '.join(sorted(T.tips_of(NPAR))))
asr = []
for _, r in cd.iterrows():
    p, X = r.裸鼹鼠位置, r.非洲鼹鼠科状态
    st = {sp: (aln[sp][p - 1] if aln[sp][p - 1] in T.AA else '?') for sp in aln}
    states, post, q, lnL = T.marginal_all(full, st)
    F, cost = T.fitch(full, st)
    g = lambda n: float(dict(zip(states, post[n])).get(X, 0.0))
    top = lambda n: '/'.join(f'{a}{v:.2f}' for a, v in sorted(zip(states, post[n]), key=lambda x: -x[1])[:2])
    asr.append(dict(P_X_非洲鼹鼠科祖先=round(g(NB), 3), P_X_父节点=round(g(NPAR), 3), P_X_豪猪亚目祖先=round(g(NH), 3),
                    父节点ML=top(NPAR), 豪猪亚目祖先ML=top(NH), Fitch_非洲鼹鼠科=''.join(sorted(F[NB])), Fitch_父节点=''.join(sorted(F[NPAR])),
                    主干上出现=('是' if g(NB) >= 0.9 and g(NPAR) < 0.5 else ('倾向是' if g(NB) >= 0.5 and g(NPAR) < 0.5 else '否'))))
if len(cd):
    cd = pd.concat([cd, pd.DataFrame(asr)], axis=1)
    # 裸鼹鼠窗口：±8，豪猪亚目代表
    REP = ['Heterocephalus glaber', 'Fukomys damarensis', 'Fukomys anselli', 'Thryonomys swinderianus', 'Petromus typicus', 'Hystrix cristata',
           'Cavia porcellus', 'Cavia tschudii', 'Chinchilla lanigera', 'Octodon degus', 'Erethizon dorsatum', 'Ctenodactylus gundi', 'Mus musculus', 'Homo sapiens']
    cd.to_csv(R('results', 't15_lineage_candidates.csv'), index=False, encoding='utf-8-sig')
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
    P(cd[['裸鼹鼠位置', '人位置', '人残基', '非洲鼹鼠科状态', '安氏鼹鼠', '家豚鼠', '其余豪猪亚目可读', '其余豪猪亚目残基', '梳趾鼠',
          '豪猪亚目外哺乳动物有X', '窗口同一性_豚鼠', 'P_X_非洲鼹鼠科祖先', 'P_X_父节点', '主干上出现', '结构域']].to_string(index=False))
    P('\n===== 每个候选 ±8 窗口（裸鼹鼠坐标；中间一列为位点）')
    for _, r in cd.iterrows():
        p = r.裸鼹鼠位置
        P(f'\n【裸鼹鼠 {p} / 人 {r.人位置}  {r.非洲鼹鼠科状态}】 {r.结构域}  {r["UniProt注释(人)"]}')
        for sp in REP:
            if sp in aln:
                w = ''.join(aln[sp][k - 1] for k in range(max(1, p - 8), min(L, p + 8) + 1))
                P(f'  {sp:26s} {w}   {aln[sp][p - 1]}')
open(R('results', 't15_lineage_scan.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
cov.to_csv(R('results', 't15_tblastn_full_coverage.csv'), index=False, encoding='utf-8-sig')
