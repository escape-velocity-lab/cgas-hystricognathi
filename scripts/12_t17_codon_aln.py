"""17a-1（任务 #17）: 全哺乳动物 cGAS 密码子比对（人坐标）+ 剪枝树（豪猪亚目分支打 {Foreground}）。

思路 = PAL2NAL：每条蛋白和人 Q8N884 做两两全局比对，按蛋白比对把对应的 3 个核苷酸搬过来。
  比对坐标 = 人编号 1–522（相对人的插入直接丢掉；这只影响插入位置，不影响 431/479/495/498）。
数据：
  - 物种名单 = #13 的 t13_species_sites.csv（281 种）
  - RefSeq 物种：#11 的 cds.fna 里找翻译后与 #13 所用蛋白完全相同的 CDS
  - 美洲豪猪 / 梳趾鼠：GenBank 蛋白的 coded_by（多外显子）→ efetch 跨度 → 拼出 CDS，核对翻译 == 蛋白
  - 其余 17 个无注释组装：#15 已经 efetch 的最后一个外显子（人 406–522 段），其余列为缺失
  - 巴西豚鼠组装太碎（只读出 431），不纳入
每条都核对：比对上 4 个位点的翻译 == #13 位点表。
用法: python 12_t17_codon_aln.py <proj_root>
"""
import sys, os, re, time, json, urllib.request
import pandas as pd
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t15_lib as T

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
OUT = R('results', 't17_codon')
os.makedirs(OUT, exist_ok=True)
RAW = R('data', 'raw', 't17')
os.makedirs(RAW, exist_ok=True)
lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

HUM = str(next(SeqIO.parse(R('data', 'raw', 'uniprot', 'Q8N884.fasta'), 'fasta')).seq)
L = len(HUM)
SITES = {431: '431', 479: '479', 495: '495', 498: '498'}
al = Align.PairwiseAligner()
al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score = -10, -0.5
al.end_insertion_score = al.end_deletion_score = 0

def to_human(prot, nuc, ref_off=1, ref=HUM):
    """prot（= nuc 的翻译）比对到 ref（人序列片段，起点 ref_off）→ {人位置: 密码子}"""
    a = al.align(ref, prot)[0]
    out = {}
    for (t0, t1), (q0, q1) in zip(*a.aligned):
        for k in range(t1 - t0):
            out[ref_off + t0 + k] = nuc[3 * (q0 + k):3 * (q0 + k) + 3]
    return out

sp_tab = pd.read_csv(R('results', 't13_species_sites.csv'))
cod, src = {}, {}

# --- RefSeq：cds.fna 按翻译匹配 ---
D = R('data', 'raw', 'ncbi', 'pkg', 'ncbi_dataset', 'data')
prot = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(D, 'protein.faa'), 'fasta')}
by_tr = {}
for r in SeqIO.parse(os.path.join(D, 'cds.fna'), 'fasta'):
    s = str(r.seq).upper()
    s = s[:len(s) - len(s) % 3]
    tr = str(Seq(s).translate()).rstrip('*')
    by_tr.setdefault(tr, []).append((r.description, s[:3 * len(tr)]))
for _, r in sp_tab[sp_tab.来源.str.startswith('RefSeq')].iterrows():
    acc = r.来源.split()[1]
    pr = prot[acc]
    hits = by_tr.get(pr, [])
    hit = next((s for d, s in hits if r.物种.split()[0] in d), hits[0][1] if hits else None)
    if hit is None:
        P('！CDS 未找到', r.物种, acc); continue
    cod[r.物种] = to_human(pr, hit)
    src[r.物种] = 'RefSeq CDS ' + acc

# --- GenBank 蛋白（美洲豪猪 / 梳趾鼠）：coded_by → 拼 CDS ---
op = urllib.request.build_opener(urllib.request.ProxyHandler())
def efetch(acc, lo, hi, strand=1):
    f = os.path.join(RAW, f'{acc}_{lo}_{hi}_{strand}.fa')
    if not os.path.exists(f) or os.path.getsize(f) < 50:
        url = (f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&seq_start={lo}&seq_stop={hi}'
               f'&strand={strand}&rettype=fasta&retmode=text')
        for k in range(4):
            try:
                open(f, 'wb').write(op.open(url, timeout=120).read()); break
            except Exception as e:
                P('efetch retry', acc, e); time.sleep(5)
        time.sleep(0.4)
    return str(next(SeqIO.parse(f, 'fasta')).seq).upper()

for rec in SeqIO.parse(R('data', 'raw', 't13', 'genbank_extra.gp'), 'genbank'):
    sp = rec.annotations['organism']
    if sp not in ('Erethizon dorsatum', 'Ctenodactylus gundi'):
        continue
    cdsf = next(f for f in rec.features if f.type == 'CDS')
    cb = cdsf.qualifiers['coded_by'][0].replace('\n', '').replace(' ', '')
    parts = [(a, int(x), int(y)) for a, x, y in re.findall(r'([A-Z0-9]+\.\d+):(\d+)\.\.(\d+)', cb)]
    acc = parts[0][0]
    lo, hi = min(p[1] for p in parts), max(p[2] for p in parts)
    span = efetch(acc, lo, hi, 1)
    nuc = ''.join(span[x - lo:y - lo + 1] for _, x, y in parts)
    if cb.startswith('complement'):
        nuc = str(Seq(nuc).reverse_complement())
    tr = str(Seq(nuc[:len(nuc) - len(nuc) % 3]).translate()).rstrip('*')
    pr = str(rec.seq)
    if tr != pr:
        P('！GenBank CDS 翻译与蛋白不符', sp, len(tr), len(pr)); continue
    cod[sp] = to_human(pr, nuc[:3 * len(pr)])
    src[sp] = f'GenBank {rec.id} coded_by {acc}:{lo}-{hi}（{len(parts)} 外显子）'

# --- 无注释组装：#15 的最后外显子小段 ---
c15 = pd.read_csv(R('results', 't15_paml', 'codon_sources.csv'))
H0, H1 = 406, 522
for _, r in c15[c15.来源.str.startswith('efetch')].iterrows():
    if r.物种 in cod:
        continue
    m = re.match(r'efetch (\S+?):(\d+)-(\d+)\(([+-])\)', r.来源)
    acc, lo, hi, s = m.group(1), m.group(2), m.group(3), '1' if m.group(4) == '+' else '2'
    nuc = str(next(SeqIO.parse(R('data', 'raw', 't15', 'exon_nuc', f'{acc}_{lo}_{hi}_{s}.fa'), 'fasta')).seq).upper()
    nuc = nuc[:len(nuc) - len(nuc) % 3]
    cod[r.物种] = to_human(str(Seq(nuc).translate()), nuc, H0, HUM[H0 - 1:H1])
    src[r.物种] = r.来源 + '（仅最后外显子）'

# --- 组装比对 + 核对 ---
seqs, bad = {}, []
st13 = sp_tab.set_index('物种')
for sp, c in cod.items():
    s = []
    for p in range(1, L + 1):
        x = c.get(p, '---')
        if len(x) < 3 or not set(x) <= set('ACGT'):
            x = '---'
        elif str(Seq(x).translate()) == '*':
            x = '---'
        s.append(x)
    for p, col in SITES.items():
        aa = str(Seq(s[p - 1]).translate()) if s[p - 1] != '---' else '?'
        exp = st13.loc[sp, col]
        if aa != exp and exp != '?':
            bad.append((sp, p, aa, exp))
    seqs[sp] = ''.join(s)
P(f'物种 {len(seqs)}；4 位点核对 vs #13：' + ('全部一致' if not bad else f'不一致 {bad}'))
for b in bad:   # 不一致的位点置缺失（不信任）
    sp, p = b[0], b[1]
    seqs[sp] = seqs[sp][:3 * (p - 1)] + '---' + seqs[sp][3 * p:]
cov = {sp: sum(s[i:i + 3] != '---' for i in range(0, len(s), 3)) for sp, s in seqs.items()}
low = [sp for sp, n in cov.items() if n < 0.3 * L and not src[sp].endswith('（仅最后外显子）')]
P('全长覆盖 <30% 剔除：', low)
for sp in low:
    seqs.pop(sp)

root, miss = T.build(R('data', 'raw', 't13', 'tree', T.TREE_MAIN), list(seqs))
P('上树缺失（剔除）：', miss)
for sp in miss:
    seqs.pop(sp)
# 豪猪亚目 = 裸鼹鼠 + 豚鼠 + 冠豪猪 的 MRCA；前景 = 冠群全部分支 + 冠群的主干（梳趾鼠分出后到冠群祖先这一段）
hy = T.mrca(root, ['Heterocephalus glaber', 'Cavia porcellus', 'Hystrix cristata'])
hy_tips = T.tips_of(hy)
fg = set(T.postorder(hy))
P(f'豪猪亚目冠群 {len(hy_tips)} 种：{sorted(hy_tips)}')
nm = lambda sp: sp.replace(' ', '_')
def nwk(n, mark):
    lab = '{Foreground}' if (mark and n in fg) else ''
    if not n.ch:
        return f'{nm(n.name)}{lab}:{n.bl:.4f}'
    return '(' + ','.join(nwk(c, mark) for c in n.ch) + ')' + lab + f':{n.bl:.4f}'
open(os.path.join(OUT, 'tree_fg.nwk'), 'w').write(nwk(root, True)[:-len(':0.0000')] + ';\n')
open(os.path.join(OUT, 'tree_plain.nwk'), 'w').write(nwk(root, False)[:-len(':0.0000')] + ';\n')
with open(os.path.join(OUT, 'aln.fasta'), 'w') as f:
    for sp in sorted(seqs):
        f.write(f'>{nm(sp)}\n{seqs[sp]}\n')
pd.DataFrame([dict(物种=sp, 目=st13.loc[sp, '目'], 科=st13.loc[sp, '科'], 来源=src[sp], 覆盖密码子=cov[sp],
                   豪猪亚目=sp in hy_tips) for sp in sorted(seqs)]).to_csv(os.path.join(OUT, 'sources.csv'), index=False, encoding='utf-8-sig')
json.dump(dict(hystricognathi=sorted(hy_tips)), open(os.path.join(OUT, 'clades.json'), 'w'), ensure_ascii=False, indent=0)
P(f'输出：{len(seqs)} 条 × {L} 密码子 → results/t17_codon/aln.fasta；树 tree_fg.nwk（{len(fg)} 条前景枝，含冠群主干）')
open(os.path.join(OUT, 'log.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
