"""04: 读出 cGAS 4 个位点（人编号 D431 / K479 / L495 / K498）在每个哺乳动物里的残基，生成 site_table_v1.csv。
- 每条异构体与人 cGAS NP_612450.2 做全局两两比对（BLOSUM62，开口 -10 / 延伸 -0.5，末端缺口不罚）
- 每个基因选一条异构体：4 个位点能读出的个数最多 → 位点窗口同一性最高 → 比对分最高
- 自检：人 0/4；UniProt A0AAX6RS70 第 463/511/527/530 = S/E/Y/T 且映射到人 431/479/495/498；NCBI 裸鼹鼠 4/4
- 质量标记：位点落在缺口、位点 ±10 aa 窗口同一性低（疑似错误基因模型）、LOW QUALITY / partial 记录、长度异常、CDS 翻译对不上
- LQ = 实测最大寿命 / 10^(a + b·log10 体重)，a、b 在 AnAge 全部哺乳动物（acceptable/high、有体重）上拟合
用法: python 04_site_table.py <data/raw> <data/processed> <results>
"""
import sys, os, re, json, math
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

RAW, PROC, RES = sys.argv[1:4]
PKG = os.path.join(RAW, 'ncbi', 'pkg', 'ncbi_dataset', 'data')
SITES = [431, 479, 495, 498]                      # 人编号
NMR = {431: 'S', 479: 'E', 495: 'Y', 498: 'T'}    # 裸鼹鼠型
HUM = {431: 'D', 479: 'K', 495: 'L', 498: 'K'}
W = 10                                            # 窗口半宽
WIN_MIN = 0.35                                    # 窗口同一性下限：经验分布在 0.20~0.40 间断档（真序列最低 0.40，错误模型 ≤0.30），取 0.35

al = PairwiseAligner()
al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score = -10, -0.5
al.end_gap_score = 0  # 末端缺口不罚

prot = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(PKG, 'protein.faa'), 'fasta')}
gid_of = {}
for r in SeqIO.parse(os.path.join(PKG, 'protein.faa'), 'fasta'):
    g = re.search(r'\[GeneID=(\d+)\]', r.description).group(1)
    gid_of.setdefault(g, []).append(r.id)
H = prot['NP_612450.2']
assert H == ''.join(str(r.seq) for r in SeqIO.parse(os.path.join(RAW, 'uniprot', 'Q8N884.fasta'), 'fasta')), '人 RefSeq 与 UniProt Q8N884 不一致'
assert ''.join(H[s - 1] for s in SITES) == 'DKLK', '人编号自检失败'

def map_to_human(q):
    """返回 dict: 人位置(1-based) → q 上对应残基（缺口为 '-'），以及 q 位置 → 人位置 和比对分"""
    a = al.align(H, q)[0]
    h2q, q2h = {}, {}
    for (hs, he), (qs, qe) in zip(*a.aligned):
        for k in range(he - hs):
            h2q[hs + k + 1] = qs + k + 1
            q2h[qs + k + 1] = hs + k + 1
    return h2q, q2h, a.score

def read_sites(q):
    h2q, q2h, sc = map_to_human(q)
    res, wid = {}, {}
    for s in SITES:
        res[s] = q[h2q[s] - 1] if s in h2q else '-'
        same = n = 0
        for p in range(s - W, s + W + 1):
            if p == s:
                continue
            n += 1
            if p in h2q and q[h2q[p] - 1] == H[p - 1]:
                same += 1
        wid[s] = same / n
    win = ''.join(q[h2q[p] - 1] if p in h2q else '-' for p in range(SITES[0] - W, SITES[-1] + W + 1))
    return res, wid, sc, win

# ---- 自检：论文用的裸鼹鼠 A0AAX6RS70 ----
U = ''.join(str(r.seq) for r in SeqIO.parse(os.path.join(RAW, 'uniprot', 'A0AAX6RS70.fasta'), 'fasta'))
assert len(U) == 554 and ''.join(U[p - 1] for p in (463, 511, 527, 530)) == 'SEYT', 'A0AAX6RS70 编号自检失败'
_, q2h, _ = map_to_human(U)
assert [q2h[p] for p in (463, 511, 527, 530)] == SITES, f'A0AAX6RS70 映射失败 {[q2h.get(p) for p in (463, 511, 527, 530)]}'
print('self-check A0AAX6RS70: 463/511/527/530 = SEYT -> human', SITES, 'OK')

# ---- 元数据 ----
titles = pd.read_csv(os.path.join(RAW, 'ncbi', 'protein_titles.tsv'), sep='\t', dtype=str).set_index('accession').title
orth = {json.loads(l)['gene_id']: json.loads(l) for l in open(os.path.join(RAW, 'ncbi', 'cgas_orthologs_mammalia.jsonl'), encoding='utf-8')}
m = pd.read_csv(os.path.join(PROC, 'orthologs_x_anage.csv'), dtype=str)

# CDS 翻译 → 蛋白序列集合（完整性核对：每条选中的蛋白都应有一条 CDS 翻译得到它）
cds_tr = set()
for r in SeqIO.parse(os.path.join(PKG, 'cds.fna'), 'fasta'):
    s = r.seq
    s = s[:len(s) // 3 * 3]
    cds_tr.add(str(s.translate()).rstrip('*'))

# ---- LQ ----
an = pd.read_csv(os.path.join(RAW, 'anage', 'anage_data.txt'), sep='\t')
an = an[an.Class == 'Mammalia']
an['anage_name'] = an.Genus + ' ' + an.Species
fit = an[an['Maximum longevity (yrs)'].notna() & an['Adult weight (g)'].notna() & an['Data quality'].isin(['acceptable', 'high'])]
b, a0 = np.polyfit(np.log10(fit['Adult weight (g)']), np.log10(fit['Maximum longevity (yrs)']), 1)
print(f'LQ allometry: n={len(fit)}  log10(L) = {a0:.4f} + {b:.4f} * log10(M)')
lq = lambda L, M: L / 10 ** (a0 + b * math.log10(M))

# 借用近缘种寿命（不是 taxid 精确匹配；单独标记，主统计可剔除）
BORROW = {
    'Nannospalax galili': ('Spalax ehrenbergi', 'AnAge 只有 S. ehrenbergi 复合种记录，N. galili 是从该复合种拆出的种'),
    'Pongo abelii': ('Pongo pygmaeus', 'AnAge 只有婆罗洲猩猩记录，苏门答腊猩猩 2001 年才拆分为独立种'),
}
# 同一物种多个亚种的基因记录：只留一个代表
REP_SUBSP = {'9612': 'Canis lupus familiaris'}

rows, wins = [], []
for _, r in m.iterrows():
    gid = r.gene_id
    cands = []
    for acc in gid_of.get(gid, []):
        q = prot[acc]
        res, wid, sc, win = read_sites(q)
        nread = sum(res[s] != '-' for s in SITES)
        cands.append((nread, float(np.mean(list(wid.values()))), sc, acc, res, wid, win))
    cands.sort(key=lambda x: (-x[0], -x[1], -x[2]))
    nread, mw, sc, acc, res, wid, win = cands[0]
    why = '唯一记录' if len(cands) == 1 else f'{len(cands)} 条异构体中位点可读数/窗口同一性最高'
    q = prot[acc]
    t = titles.get(acc, '')
    flags = []
    for s in SITES:
        if res[s] == '-':
            flags.append(f'{s}落在缺口')
        elif wid[s] < WIN_MIN:
            flags.append(f'{s}窗口同一性{wid[s]:.2f}')
    if 'LOW QUALITY' in t:
        flags.append('LOW QUALITY PROTEIN')
    if 'partial' in t:
        flags.append('partial')
    if not 450 <= len(q) <= 650:
        flags.append(f'长度{len(q)}')
    if 'X' in win:
        flags.append('窗口含X')
    cds_ok = q in cds_tr
    if not cds_ok:
        flags.append('CDS翻译对不上')
    ann = (orth[gid].get('annotations') or [{}])[0]
    # 寿命
    src, L, M, dq, an_name = '', math.nan, math.nan, '', ''
    if isinstance(r.anage_name, str):
        src, an_name = 'taxid 精确', r.anage_name
        L, M, dq = float(r['Maximum longevity (yrs)']) if isinstance(r['Maximum longevity (yrs)'], str) else math.nan, \
            float(r['Adult weight (g)']) if isinstance(r['Adult weight (g)'], str) else math.nan, r['Data quality']
    elif r.taxname in BORROW:
        bn, bwhy = BORROW[r.taxname]
        x = an[an.anage_name == bn].iloc[0]
        src, an_name = f'借用近缘种（{bwhy}）', bn
        L, M, dq = x['Maximum longevity (yrs)'], x['Adult weight (g)'], x['Data quality']
    q_lq = lq(L, M) if L == L and M == M else math.nan
    grp = '' if q_lq != q_lq else ('长寿' if q_lq >= 1.5 else ('短寿' if q_lq < 0.9 else '中间'))
    hits = [s for s in SITES if res[s] == NMR[s]]
    rep = '是'
    if r.species_taxid in REP_SUBSP and r.taxname != REP_SUBSP[r.species_taxid]:
        rep = '否（同种亚种重复）'
    usable = (nread == 4) and not any(('缺口' in f) or ('窗口同一性' in f) or ('窗口含X' in f) for f in flags)
    rows.append({
        '物种': r.taxname, '俗名': r.common_name if isinstance(r.common_name, str) else '', '目': r.order, '科': r.family,
        'taxid': r.tax_id, '物种级taxid': r.species_taxid, 'GeneID': gid, 'RefSeq蛋白': acc, '蛋白标题': t,
        '异构体选择': why, '蛋白长度': len(q), '组装accession': ann.get('assembly_accession', ''),
        '注释版本': ann.get('annotation_name', ''), '注释日期': ann.get('annotation_release_date', ''),
        'AnAge学名': an_name, '寿命来源': src, '最大寿命_年': L, '成年体重_g': M, 'AnAge数据质量': dq,
        'LQ': round(q_lq, 3) if q_lq == q_lq else math.nan, '寿命分组': grp,
        'D431': res[431], 'K479': res[479], 'L495': res[495], 'K498': res[498],
        '裸鼹鼠型命中数': len(hits), '命中位点': '/'.join(f'{HUM[s]}{s}{NMR[s]}' for s in hits),
        '窗口同一性_431/479/495/498': '/'.join(f'{wid[s]:.2f}' for s in SITES),
        '质量标记': '; '.join(flags), '4位点可用': '是' if usable else '否', '物种代表': rep,
        'CDS核对': '通过' if cds_ok else '不通过',
    })
    wins.append((r.taxname, acc, win, '/'.join(res[s] for s in SITES), len(hits)))

df = pd.DataFrame(rows).sort_values(['目', '科', '物种']).reset_index(drop=True)
os.makedirs(RES, exist_ok=True)
df.to_csv(os.path.join(RES, 'site_table_v1.csv'), index=False, encoding='utf-8-sig')

# 比对窗口（人 421~508），供人工复查
with open(os.path.join(RES, 'site_windows_v1.txt'), 'w', encoding='utf-8') as f:
    ruler = ''.join('^' if p in SITES else ' ' for p in range(SITES[0] - W, SITES[-1] + W + 1))
    f.write(f'# 人 NP_612450.2 第 {SITES[0]-W}~{SITES[-1]+W} 位；^ = 431/479/495/498；「-」= 该物种此处为缺口\n')
    f.write(f'{"Homo sapiens (ref)":42s} {"":16s} {H[SITES[0]-W-1:SITES[-1]+W]}\n{"":42s} {"":16s} {ruler}\n')
    for sp, acc, win, st, n in sorted(wins, key=lambda x: (-x[4], x[0])):
        f.write(f'{sp:42s} {acc:16s} {win}  {st}  hits={n}\n')

# ---- 自检汇总 ----
def row(sp):
    return df[df.物种 == sp].iloc[0]
for sp, exp in (('Homo sapiens', 'DKLK'), ('Heterocephalus glaber', 'SEYT'), ('Mus musculus', 'DKLR'),
                ('Rattus norvegicus', 'DKLK'), ('Fukomys damarensis', 'SEYT')):
    x = row(sp); got = x.D431 + x.K479 + x.L495 + x.K498
    print(f'control {sp:24s} expect {exp} got {got} {"OK" if got == exp else "FAIL"}  ({x.RefSeq蛋白})')
print('rows', len(df), '| species reps', (df.物种代表 == '是').sum(), '| 4-site usable', (df['4位点可用'] == '是').sum())
