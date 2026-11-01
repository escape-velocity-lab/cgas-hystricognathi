"""32（任务 #30A）: 用本地 BLAST+ 在裸鼹鼠 T2T 两个单倍型、RefSeq 裸鼹鼠 cGAS 所在片段、开普沙丘鼹鼠组装里
读人 cGAS 431/479/495/498 的对应残基，三方对照，并核单倍型 / 拷贝数。
- 查询 A = 人 Q8N884 第 400~522 位（与 07 完全相同）；残基读法与 09 相同（HSP 同一性 ≥ 0.45 过滤旁系同源）
- 查询 B = 裸鼹鼠 NP_001431126.1 全长：看全长有几个位点（拷贝）、每个 HSP 与 RefSeq 蛋白的逐位差异
- 再把含四位点的那个外显子的核苷酸序列在各组装间逐碱基比对
用法: python 32_t30_t2t_check.py <t30_genomes_dir> <query_dir> <blast_bin> <out_dir> [标签,标签]
"""
import sys, os, gzip, subprocess, shutil
import xml.etree.ElementTree as ET
from Bio import SeqIO
from Bio.Seq import Seq

G, QD, BIN, OUT = sys.argv[1:5]
ONLY = set(sys.argv[5].split(',')) if len(sys.argv) > 5 else None  # 只跑这些标签（调试 / 分批）
os.makedirs(OUT, exist_ok=True)
Q0 = 400
SITES = [431, 479, 495, 498]
NMR = dict(zip(SITES, 'SEYT'))
# (标签, 组装, 原始文件)
ASM = [
    ('T2T_hap1', 'GCA_988278965.1 mHetGla1.1.hap1', 'hap1.fna.gz'),
    ('T2T_hap2', 'GCA_988279095.1 mHetGla1.1.hap2', 'hap2.fna.gz'),
    ('RefSeq_locus', 'GCF_053883595.1 NW_028258164.1:6796961-6905671', 'refseq_cgas_locus.fa'),
    ('Bsuillus', 'GCA_988425685.1 mBatSui1.1.pri (CGCGHY01)', 'bsui.fasta.gz'),
]
H = ''.join(str(r.seq) for r in SeqIO.parse(os.path.join(QD, 'uniprot', 'Q8N884.fasta'), 'fasta'))
qa = os.path.join(OUT, 'qA_human400_522.fa'); open(qa, 'w').write('>Q8N884_400_522\n' + H[Q0 - 1:522] + '\n')
qb = os.path.join(QD, 't30_query_NP_001431126.fa')
NP = str(next(SeqIO.parse(qb, 'fasta')).seq)
qc = os.path.join(QD, 'uniprot', 'A0AAX6RS70.fasta')
# 稿件 Table 1 的 7 个位置 + 四位点（A0AAX6RS70 坐标）
T1 = {43: 'C', 72: 'S', 86: 'K', 233: 'G', 290: 'Q', 358: 'P', 420: 'S', 463: 'S', 511: 'E', 527: 'Y', 530: 'T'}

def run(*a):
    r = subprocess.run(a, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"{a[0]} failed: {r.stderr[-800:]}")
    return r.stdout

def db_for(tag, fn):
    src = os.path.join(G, fn)
    fa = os.path.join(G, tag + '.fna')
    if not os.path.exists(fa):
        op = gzip.open if fn.endswith('.gz') else open
        with op(src, 'rb') as i, open(fa, 'wb') as o:
            # ENA 的 FASTA 头是 >ENA|CGCGHY010000001|CGCGHY010000001.1 ...，makeblastdb -parse_seqids 不认 ENA 前缀（坑 72）
            for line in i:
                if line.startswith(b'>ENA|'):
                    line = b'>' + line[5:].split(b'|', 1)[1]
                o.write(line)
    if not os.path.exists(fa + '.ndb'):
        run(os.path.join(BIN, 'makeblastdb.exe'), '-in', fa, '-dbtype', 'nucl', '-parse_seqids', '-out', fa)
    return fa

def tbn(db, q, out):
    if not os.path.exists(out):
        run(os.path.join(BIN, 'tblastn.exe'), '-query', q, '-db', db, '-outfmt', '5', '-evalue', '1e-3',
            '-max_target_seqs', '20', '-seg', 'no', '-num_threads', '8', '-out', out)
    hs = []
    for hit in ET.parse(out).getroot().iter('Hit'):
        acc = hit.findtext('Hit_accession'); defl = hit.findtext('Hit_def')
        for h in hit.iter('Hsp'):
            hs.append(dict(acc=acc, defl=defl, bits=float(h.findtext('Hsp_bit-score')), ev=float(h.findtext('Hsp_evalue')),
                           qf=int(h.findtext('Hsp_query-from')), qt=int(h.findtext('Hsp_query-to')),
                           hf=int(h.findtext('Hsp_hit-from')), ht=int(h.findtext('Hsp_hit-to')), frame=int(h.findtext('Hsp_hit-frame')),
                           ident=int(h.findtext('Hsp_identity')), alen=int(h.findtext('Hsp_align-len')),
                           q=h.findtext('Hsp_qseq'), s=h.findtext('Hsp_hseq')))
    return hs

def qmap(h, off):
    m, qp = {}, h['qf'] + off
    for a, b in zip(h['q'], h['s']):
        if a != '-':
            m[qp] = b; qp += 1
    return m

def ntseg(db, acc, a, b, minus):
    s = run(os.path.join(BIN, 'blastdbcmd.exe'), '-db', db, '-entry', acc, '-range', f'{min(a,b)}-{max(a,b)}',
            '-strand', 'minus' if minus else 'plus', '-outfmt', '%s').strip()
    return s

rep = []; exon = {}
for tag, asm, fn in ASM:
    if ONLY and tag not in ONLY:
        continue
    if not os.path.exists(os.path.join(G, fn)):
        rep.append(f'## {tag} — {asm}\n文件缺失：{fn}\n'); continue
    db = db_for(tag, fn)
    A = [h for h in tbn(db, qa, os.path.join(OUT, f'{tag}_qA.xml')) if h['ident'] / h['alen'] >= 0.45]
    B = tbn(db, qb, os.path.join(OUT, f'{tag}_qB.xml'))
    rep.append(f'## {tag} — {asm}\n')
    rep.append(f'查询 A（人 400~522）过滤后 HSP {len(A)} 条：')
    for h in sorted(A, key=lambda x: -x['bits']):
        rep.append(f"- {h['acc']}:{h['hf']}-{h['ht']} frame={h['frame']} id={h['ident']}/{h['alen']} bits={h['bits']:.0f} 人{h['qf']+Q0-1}-{h['qt']+Q0-1}  {h['defl'][:70]}")
    if A:
        best = max(A, key=lambda x: x['bits'])
        rec = []
        for s in SITES:
            cov = sorted([h for h in A if h['qf'] + Q0 - 1 <= s <= h['qt'] + Q0 - 1], key=lambda h: (h['acc'] != best['acc'], -h['bits']))
            if not cov:
                rec.append(f'{s}: ?（无 HSP 覆盖）'); continue
            m = qmap(cov[0], Q0 - 1)
            alt = sorted({qmap(x, Q0 - 1).get(s) for x in cov[1:]} - {None})
            win = ''.join(m.get(p, ' ') for p in range(s - 15, s + 16))
            rec.append(f"{s}: **{m.get(s,'?')}**（裸鼹鼠型 {NMR[s]}）窗口 {win}  其他 HSP 残基 {''.join(alt) or '无'}")
        rep.append('\n四位点：\n' + '\n'.join('- ' + r for r in rec))
        # 外显子核苷酸
        minus = best['frame'] < 0
        exon[tag] = ntseg(db, best['acc'], best['hf'], best['ht'], minus)
        rep.append(f"\n含四位点的外显子（{best['acc']}:{min(best['hf'],best['ht'])}-{max(best['hf'],best['ht'])} {'负' if minus else '正'}链）"
                   f"翻译：{Seq(exon[tag][:len(exon[tag])//3*3]).translate()}")
    # 查询 B：全长
    Bf = [h for h in B if h['ident'] / h['alen'] >= 0.45]
    contigs = sorted({h['acc'] for h in Bf})
    rep.append(f'\n查询 B（NP_001431126.1 全长 {len(NP)} aa）同一性 ≥0.45 的 HSP {len(Bf)} 条，落在 {len(contigs)} 条序列：{", ".join(contigs)}')
    cover = set()
    for h in sorted(Bf, key=lambda x: (x['acc'], x['qf'])):
        cover |= set(range(h['qf'], h['qt'] + 1))
        diffs = [f"{NP[p-1]}{p}{b}" for p, b in qmap(h, 0).items() if b not in ('-', NP[p-1])]
        rep.append(f"- {h['acc']}:{h['hf']}-{h['ht']} f={h['frame']} NP{h['qf']}-{h['qt']} id={h['ident']}/{h['alen']} bits={h['bits']:.0f} 与 RefSeq 蛋白不同：{' '.join(diffs) or '无'}")
    rep.append(f'NP_001431126.1 被覆盖 {len(cover)}/{len(NP)} 位\n')
    # 查询 C：A0AAX6RS70（稿件 Table 1 的坐标系），读 7 个非洲鼹鼠科特有位置 + 四位点
    C = [h for h in tbn(db, qc, os.path.join(OUT, f'{tag}_qC.xml')) if h['ident'] / h['alen'] >= 0.45]
    if C:
        bc = max(C, key=lambda x: x['bits'])['acc']
        out = []
        for p, x in T1.items():
            cov = sorted([h for h in C if h['qf'] <= p <= h['qt']], key=lambda h: (h['acc'] != bc, -h['bits']))
            r = qmap(cov[0], 0).get(p, '?') if cov else '?'
            alt = sorted({qmap(h, 0).get(p) for h in cov[1:]} - {None, r}) if cov else []
            out.append(f"{p}{x}:{r}{'(' + ''.join(alt) + ')' if alt else ''}")
        rep.append('A0AAX6RS70 坐标读数（位置+裸鼹鼠残基:本组装残基）：' + '  '.join(out) + '\n')

# 外显子逐碱基对照
if exon:
    rep.append('## 含四位点外显子的核苷酸对照（以 RefSeq_locus 为参照）\n')
    ref = exon.get('RefSeq_locus')
    for t, s in exon.items():
        line = f'- {t}: 长 {len(s)} nt'
        if ref and len(s) == len(ref):
            d = [f'{ref[i]}{i+1}{s[i]}' for i in range(len(s)) if s[i] != ref[i]]
            line += f'；与 RefSeq 差异 {len(d)} 处 {" ".join(d[:30])}'
        elif ref:
            line += f'；长度与 RefSeq（{len(ref)}）不同，见序列'
        rep.append(line)
    open(os.path.join(OUT, 'exon_nt.fa'), 'w').write(''.join(f'>{t}\n{s}\n' for t, s in exon.items()))
txt = '\n'.join(rep)
open(os.path.join(OUT, 't30_t2t_raw.md'), 'w', encoding='utf-8').write(txt)
print(txt)
