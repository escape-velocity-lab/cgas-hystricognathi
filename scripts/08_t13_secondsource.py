"""08（任务 #13 · 13a）: 5 个关键物种的 cGAS 4 位点第二来源核对 + RefSeq 转录本证据等级。
- 来源 A：RefSeq（任务 #11 选中的那条 + 同基因其他异构体）
- 来源 B：UniProt / Ensembl（独立基因注释流程）、GenBank 提交者注释（KFO37345 等）
- 来源 C：07 的 tblastn（另一套组装，见 08 不处理，由 09 统一读）
- 每条 XP/NP 取对应 XM/NM 的 GenBank 记录，摘出「Supporting evidence includes ... RNAseq ...」那句
输出: results/t13_secondsource.csv、results/t13_windows_13a.txt、data/raw/t13/refseq_mrna_evidence.tsv
用法: python 08_t13_secondsource.py <data/raw> <results>
"""
import sys, os, re, json, time, urllib.request, urllib.parse
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

RAW, RES = sys.argv[1:3]
T13 = os.path.join(RAW, 't13')
op = urllib.request.build_opener(urllib.request.ProxyHandler())
E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

def get(url, data=None):
    for k in range(4):
        try:
            return op.open(url, data=data, timeout=90).read().decode('utf-8', 'replace')
        except Exception as e:
            print('retry', url[:80], e); time.sleep(3)
    raise RuntimeError(url)

def efetch(db, ids, rettype):
    return get(E + 'efetch.fcgi', urllib.parse.urlencode({'db': db, 'id': ','.join(ids), 'rettype': rettype, 'retmode': 'text'}).encode())

H = ''.join(str(r.seq) for r in SeqIO.parse(os.path.join(RAW, 'uniprot', 'Q8N884.fasta'), 'fasta'))
SITES = [431, 479, 495, 498]
NMR = dict(zip(SITES, 'SEYT'))
al = PairwiseAligner()
al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score = -10, -0.5
al.end_gap_score = 0

def sites(q):
    a = al.align(H, q)[0]
    h2q = {}
    for (hs, he), (qs, qe) in zip(*a.aligned):
        for k in range(he - hs):
            h2q[hs + k + 1] = qs + k + 1
    res = {s: (q[h2q[s] - 1] if s in h2q else '-') for s in SITES}
    pos = {s: h2q.get(s) for s in SITES}
    win = {s: ''.join(q[h2q[p] - 1] if p in h2q else '-' for p in range(s - 15, s + 16)) for s in SITES}
    return res, pos, win

# ---- 收集序列 ----
seqs = []  # (物种, 来源, accession, 说明, 序列)
SP = ['Cavia porcellus', 'Chinchilla lanigera', 'Octodon degus', 'Heterocephalus glaber', 'Fukomys damarensis',
      'Erethizon dorsatum', 'Ctenodactylus gundi']
for r in SeqIO.parse(os.path.join(RAW, 'ncbi', 'pkg', 'ncbi_dataset', 'data', 'protein.faa'), 'fasta'):
    m = re.search(r'\[organism=([^\]]+)\]', r.description).group(1)
    if m in SP:
        seqs.append((m, 'RefSeq（NCBI Eukaryotic Annotation Pipeline）', r.id, r.description.split(' ', 1)[1], str(r.seq)))
# GenBank 提交者注释（非 RefSeq）
gb_ids = ['KFO37345.1', 'KAM6166271.1', 'KAM5256511.1']
fa = efetch('protein', gb_ids, 'fasta')
open(os.path.join(T13, 'genbank_extra.faa'), 'w', encoding='utf-8').write(fa)
gbp = efetch('protein', gb_ids, 'gp')
open(os.path.join(T13, 'genbank_extra.gp'), 'w', encoding='utf-8').write(gbp)
for r in SeqIO.parse(os.path.join(T13, 'genbank_extra.faa'), 'fasta'):
    sp = re.search(r'\[([^\]]+)\]\s*$', r.description).group(1)
    seqs.append((sp, 'GenBank 提交者注释', r.id, r.description.split(' ', 1)[1], str(r.seq)))
# UniProt（TrEMBL）
up = get('https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=' + urllib.parse.quote(
    '(gene:CGAS OR gene:MB21D1) AND taxonomy_id:33550'))
open(os.path.join(T13, 'uniprot_hystricognathi.fasta'), 'w', encoding='utf-8').write(up)
for r in SeqIO.parse(os.path.join(T13, 'uniprot_hystricognathi.fasta'), 'fasta'):
    sp = re.search(r'OS=(\S+ \S+)', r.description).group(1)
    seqs.append((sp, 'UniProt TrEMBL', r.id.split('|')[1], r.description.split(' ', 1)[1][:80], str(r.seq)))
# Ensembl compara 直系同源（人 ENSG00000164430）
ens = json.load(open(os.path.join(T13, 'ensembl_hystricognathi_orth.json'), encoding='utf-8'))
for h in ens['data'][0]['homologies']:
    t = h['target']
    sp = t['species'].replace('_female', '').replace('_', ' ').capitalize()
    seqs.append((sp, 'Ensembl（release 当前，compara 直系同源）', t['protein_id'], h['type'], t['align_seq'].replace('-', '')))
# Ensembl 雄性裸鼹鼠组装（HetGla 另一套组装）
g = get('https://rest.ensembl.org/sequence/id/ENSHGLG00100034506?type=protein;multiple_sequences=1;content-type=text/x-fasta')
open(os.path.join(T13, 'ensembl_hetgla_male.faa'), 'w', encoding='utf-8').write(g)
for r in SeqIO.parse(os.path.join(T13, 'ensembl_hetgla_male.faa'), 'fasta'):
    seqs.append(('Heterocephalus glaber', 'Ensembl（雄性个体组装 HetGla_1.0）', r.id, 'ENSHGLG00100034506', str(r.seq)))

rows, wl = [], []
for sp, src, acc, desc, q in seqs:
    res, pos, win = sites(q)
    rows.append({'物种': sp, '来源': src, 'accession': acc, '说明': desc[:90], '长度': len(q),
                 **{f'{s}': res[s] for s in SITES}, '命中': sum(res[s] == NMR[s] for s in SITES),
                 '本序列中的位置': '/'.join(str(pos[s]) for s in SITES)})
    wl.append((sp, src, acc, win, res))
df = pd.DataFrame(rows).sort_values(['物种', '来源', 'accession'])
df.to_csv(os.path.join(RES, 't13_secondsource.csv'), index=False, encoding='utf-8-sig')
print(df.to_string())

with open(os.path.join(RES, 't13_windows_13a.txt'), 'w', encoding='utf-8') as f:
    f.write('# 每个位点 ±15 aa（人编号），第 16 个字符 = 位点本身；「-」= 比对缺口\n')
    for s in SITES:
        f.write(f'\n== 人 {H[s-1]}{s}（裸鼹鼠型 {NMR[s]}）==\n{"Homo sapiens Q8N884":60s} {H[s-16:s+15]}\n{"":60s} {" "*15}^\n')
        for sp, src, acc, win, res in sorted(wl):
            f.write(f'{(sp + " | " + acc + " | " + src)[:60]:60s} {win[s]}  {res[s]}\n')

# ---- RefSeq mRNA 证据 ----
acc_rs = [a for sp, src, a, d, q in seqs if src.startswith('RefSeq')]
gp = efetch('protein', acc_rs, 'gp')
xm = dict(re.findall(r'VERSION\s+(\S+)[\s\S]*?DBSOURCE\s+REFSEQ: accession (\S+)', gp))
nm = efetch('nuccore', list(xm.values()), 'gb')
open(os.path.join(T13, 'refseq_mrna.gb'), 'w', encoding='utf-8').write(nm)
ev = []
for rec in nm.split('//\n'):
    m = re.search(r'VERSION\s+(\S+)', rec)
    if not m:
        continue
    flat = re.sub(r'\s+', ' ', rec)
    s = re.search(r'Supporting evidence includes similarity to: ([^"]*?)(?:"|$)', flat)
    st = re.search(r'(PROVISIONAL|VALIDATED|REVIEWED|INFERRED|MODEL|PREDICTED) REFSEQ', flat)
    org = re.search(r'ORGANISM (\S+ \S+)', flat).group(1)
    prot = [p for p, x in xm.items() if x == m.group(1)]
    ev.append({'物种': org, 'mRNA': m.group(1), '蛋白': prot[0] if prot else '', 'RefSeq状态': st.group(1) if st else '',
               '证据': s.group(1).strip() if s else ''})
evd = pd.DataFrame(ev)
evd.to_csv(os.path.join(T13, 'refseq_mrna_evidence.tsv'), sep='\t', index=False)
print(evd.to_string())
