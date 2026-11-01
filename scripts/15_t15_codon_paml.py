"""15（任务 #15c）: 495-Y（裸鼹鼠 527）的位点选择信号——PAML codeml。只当线索，功效低。

数据：啮齿目全部可用物种（#13 的 70 种）在 cGAS **最后一个外显子**上的密码子比对。
  这个外显子覆盖裸鼹鼠 438–554 / 人 406–522，包含全部 4 个位点；#13 已确认它在豪猪亚目里是一条 HSP（不跨内含子）。
  - RefSeq 物种：从 #11 的 cds.fna 里找翻译后与所选蛋白相同的 CDS，按蛋白投影取密码子
  - 无注释组装 + 美洲豪猪 / 梳趾鼠：取 tblastn 覆盖该外显子的 HSP 的基因组坐标，efetch 只下这一小段核酸（约 350 bp，不下基因组），
    翻译后比对到裸鼹鼠外显子，取密码子
  每条都核对「翻译 == 投影蛋白」。
树：Upham 2019 主树剪到这些物种，去根（codeml 自己估枝长）。
检验：
  1. 位点模型 M7 vs M8、M8a vs M8（全体分支上有没有 ω>1 的位点），看 495 的 BEB 后验
  2. 分支-位点模型 A vs A-null：前景枝 = 495-Y 出现的枝（在这棵啮齿类树上用 Mk-ER 边际概率 <0.5 → ≥0.5 判定）
用法: python 15_t15_codon_paml.py <proj_root> <codeml.exe>（codeml 参数保留但本脚本只做数据准备；codeml 见 17）
"""
import sys, os, re, glob, json, math, subprocess, time, urllib.request
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from scipy.stats import chi2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t15_lib as T

ROOT, CODEML = sys.argv[1:3]
R = lambda *p: os.path.join(ROOT, *p)
WD = R('results', 't15_paml')
os.makedirs(WD, exist_ok=True)
NUC = R('data', 'raw', 't15', 'exon_nuc')
os.makedirs(NUC, exist_ok=True)
lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

NMRSEQ = str(next(SeqIO.parse(R('data', 'raw', 'uniprot', 'A0AAX6RS70.fasta'), 'fasta')).seq)
E0, E1 = 438, 554                # 裸鼹鼠外显子范围（含）
EX = NMRSEQ[E0 - 1:E1]
SITE = 527                       # = 人 495
al = Align.PairwiseAligner()
al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score = -10, -0.5
al.end_insertion_score = al.end_deletion_score = 0

def codons_on_ref(ref, ref_off, prot, nuc):
    """prot（nuc 的翻译，等长 ×3）比对到 ref → {裸鼹鼠位置: 密码子}"""
    a = al.align(ref, prot)[0]
    out = {}
    for (t0, t1), (q0, q1) in zip(*a.aligned):
        for k in range(t1 - t0):
            out[ref_off + t0 + k] = nuc[3 * (q0 + k):3 * (q0 + k) + 3]
    return out

sp_tab = pd.read_csv(R('results', 't13_species_sites.csv'))
rod = sp_tab[sp_tab.目 == 'Rodentia']
cod, src = {}, {}

# --- RefSeq：cds.fna ---
cds = [(r.description, str(r.seq).upper()) for r in SeqIO.parse(R('data', 'raw', 'ncbi', 'pkg', 'ncbi_dataset', 'data', 'cds.fna'), 'fasta')]
prot = {r.id: str(r.seq) for r in SeqIO.parse(R('data', 'raw', 'ncbi', 'pkg', 'ncbi_dataset', 'data', 'protein.faa'), 'fasta')}
for _, r in rod[rod.来源.str.startswith('RefSeq')].iterrows():
    acc = r.来源.split()[1]
    pr = prot[acc]
    hit = None
    for d, s in cds:
        if f'[organism={r.物种}' not in d and r.物种.split()[0] not in d:
            continue
        s2 = s[:len(s) - len(s) % 3]
        tr = str(Seq(s2).translate()).rstrip('*')
        if tr == pr:
            hit = s2[:3 * len(pr)]; break
    if hit is None:
        P('！CDS 未找到', r.物种, acc); continue
    cod[r.物种] = codons_on_ref(NMRSEQ, 1, pr, hit)
    src[r.物种] = 'RefSeq CDS ' + acc

# --- 其余：tblastn HSP → efetch 小段核酸 ---
op = urllib.request.build_opener(urllib.request.ProxyHandler())
def efetch(acc, lo, hi, strand):
    f = os.path.join(NUC, f'{acc}_{lo}_{hi}_{strand}.fa')
    if not os.path.exists(f) or os.path.getsize(f) < 50:
        url = (f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&seq_start={lo}&seq_stop={hi}'
               f'&strand={strand}&rettype=fasta&retmode=text')
        for k in range(4):
            try:
                open(f, 'wb').write(op.open(url, timeout=60).read()); break
            except Exception as e:
                P('efetch retry', acc, e); time.sleep(5)
        time.sleep(0.4)
    return str(next(SeqIO.parse(f, 'fasta')).seq).upper()

def best_exon_hsp(xmlf, qlo, qhi):
    root = ET.parse(xmlf).getroot()
    best = None
    for hit in root.iter('Hit'):
        for h in hit.iter('Hsp'):
            qf, qt = int(h.findtext('Hsp_query-from')), int(h.findtext('Hsp_query-to'))
            idn, alen = int(h.findtext('Hsp_identity')), int(h.findtext('Hsp_align-len'))
            if qf <= qlo and qt >= qhi and idn / alen >= 0.45:
                b = float(h.findtext('Hsp_bit-score'))
                if best is None or b > best[0]:
                    best = (b, hit.findtext('Hit_accession'), int(h.findtext('Hsp_hit-from')), int(h.findtext('Hsp_hit-to')),
                            int(h.findtext('Hsp_hit-frame')), idn / alen)
    return best

src12 = open(R('src', '12_t15_tblastn_full.py'), encoding='utf-8').read()
JOBS = {sp: db for sp, db in re.findall(r"\('([A-Z][a-z]+ [a-z]+)', '([A-Z0-9]+)'\)", src12)}
XML = {sp: R('data', 'raw', 't15', 'tblastn_full', f'{db}.xml') for sp, db in JOBS.items()}
# 美洲豪猪 / 梳趾鼠：用 #13 的 tblastn（查询 = 人 400–522，覆盖该外显子；位点 527 ↔ 人 495）
XML_H = {'Erethizon dorsatum': R('data', 'raw', 't13', 'tblastn', 'JAQMCP01.xml'), 'Ctenodactylus gundi': R('data', 'raw', 't13', 'tblastn', 'PVKB01.xml')}
for sp in rod.物种:
    if sp in cod:
        continue
    if sp in XML and os.path.exists(XML[sp]):
        b = best_exon_hsp(XML[sp], SITE, SITE)          # 只要求覆盖 527；整体覆盖率在后面按 75% 过滤
    elif sp in XML_H:
        b = best_exon_hsp(XML_H[sp], 495 - 399, 495 - 399)   # 07 的查询坐标 = 人编号 − 399
    else:
        P('！无数据', sp); continue
    if b is None:
        P('！无覆盖外显子的 HSP', sp); continue
    _, acc, hf, ht, fr, idf = b
    lo, hi = min(hf, ht), max(hf, ht)
    strand = 1 if fr > 0 else 2
    nuc = efetch(acc, lo, hi, strand)
    nuc = nuc[:len(nuc) - len(nuc) % 3]
    pr = str(Seq(nuc).translate())
    c = codons_on_ref(EX, E0, pr, nuc)
    cod[sp] = c
    src[sp] = f'efetch {acc}:{lo}-{hi}({"+" if strand == 1 else "-"}) id={idf:.2f}'

# --- 组装密码子比对（裸鼹鼠 438–554） ---
P(f'外显子 裸鼹鼠 {E0}–{E1}（{E1 - E0 + 1} 个密码子），位点 {SITE} 在比对第 {SITE - E0 + 1} 列')
seqs, bad = {}, []
for sp, c in cod.items():
    s = ''
    for p in range(E0, E1 + 1):
        x = c.get(p, '---')
        if len(x) < 3 or '-' in x:
            x = '---'
        aa = str(Seq(x).translate()) if x != '---' else '-'
        if aa == '*':
            x = '---'
        s += x
    aa527 = str(Seq(s[3 * (SITE - E0):3 * (SITE - E0) + 3]).translate()) if s[3 * (SITE - E0)] != '-' else '-'
    exp = sp_tab.set_index('物种').loc[sp, '495']
    if aa527 != exp:
        bad.append((sp, aa527, exp))
    seqs[sp] = s
P('核对：密码子比对第 527 位翻译 vs #13 位点表 495：' + ('全部一致' if not bad else f'不一致 {bad}'))
cov = {sp: sum(s[i:i + 3] != '---' for i in range(0, len(s), 3)) for sp, s in seqs.items()}
low = [sp for sp, n in cov.items() if n < 0.75 * (E1 - E0 + 1)]
P(f'物种 {len(seqs)}；覆盖 <75% 的剔除：{low}')
for sp in low:
    seqs.pop(sp)
pd.DataFrame([dict(物种=sp, 来源=src[sp], 覆盖密码子=cov[sp]) for sp in cod]).to_csv(os.path.join(WD, 'codon_sources.csv'), index=False, encoding='utf-8-sig')

# --- 树 + 前景枝 ---
root, miss = T.build(R('data', 'raw', 't13', 'tree', T.TREE_MAIN), list(seqs))
P('上树缺失：', miss)
for sp in miss:
    seqs.pop(sp)
st = sp_tab.set_index('物种')['495'].to_dict()
states, post, q, lnL = T.marginal_all(root, {sp: st[sp] for sp in seqs})
iy = states.index('Y')
P('参考：只用啮齿类这棵树重建时，P(Y) 豪猪亚目祖先 = '
  f"{post[T.mrca(root, ['Cavia porcellus', 'Heterocephalus glaber', 'Hystrix cristata'])][iy]:.2f}"
  '（#13 在 277 种哺乳动物树上是 0.09；外群取样少会把祖先往 Y 拉——所以前景枝不用这棵树定，见下）')
# 前景枝事先按 #13 的结论定（豪猪亚目祖先 = H，Y 独立出现；松鼠属另一次）：每个 Y 支系的主干
FG_CLADES = [['Heterocephalus glaber', 'Fukomys damarensis', 'Fukomys anselli'], ['Cavia porcellus', 'Cavia tschudii'],
             ['Chinchilla lanigera'], ['Capromys pilorides'], ['Octodon degus', 'Octomys mimax', 'Tympanoctomys barrerae', 'Ctenomys sociabilis'],
             ['Sciurus vulgaris', 'Sciurus carolinensis']]
fg = set()
for cl in FG_CLADES:
    cl = [x for x in cl if x in seqs]
    n = T.mrca(root, cl) if len(cl) > 1 else [x for x in T.postorder(root) if x.name == cl[0]][0]
    assert set(T.tips_of(n)) == set(cl), (cl, T.tips_of(n))
    fg.add(n)
P('前景枝（事先指定 = #13 的「H 祖先、Y 多次独立出现」情景）：')
for n in fg:
    P(f'  → {", ".join(sorted(T.tips_of(n)))}  （啮齿类树 Mk-ER：P(Y) 父 {post[n.up][iy]:.2f} → 子 {post[n][iy]:.2f}）')
code = {sp: f'S{i:02d}' for i, sp in enumerate(sorted(seqs))}
def nwk(n, mark):
    lab = ' #1' if (mark and n in fg) else ''
    if not n.ch:
        return code[n.name] + lab
    return '(' + ','.join(nwk(c, mark) for c in n.ch) + ')' + lab
def unroot(n, mark):
    # 根是二叉：把一个子节点展开成三叉；被展开的子节点上的前景标记转到另一侧（同一条去根后的枝）
    a, b = n.ch
    if not a.ch:
        a, b = b, a
    m = ' #1' if (mark and (a in fg or b in fg)) else ''
    return '(' + ','.join(nwk(c, mark) for c in a.ch) + ',' + nwk(b, False) + m + ');'
open(os.path.join(WD, 'aln.phy'), 'w').write(f' {len(seqs)} {len(next(iter(seqs.values())))}\n' + ''.join(f'{code[sp]}  {s}\n' for sp, s in sorted(seqs.items())))
open(os.path.join(WD, 'tree_plain.nwk'), 'w').write(f' {len(seqs)} 1\n' + unroot(root, False) + '\n')
open(os.path.join(WD, 'tree_fg.nwk'), 'w').write(f' {len(seqs)} 1\n' + unroot(root, True) + '\n')
json.dump(code, open(os.path.join(WD, 'codes.json'), 'w'), ensure_ascii=False, indent=0)

# codeml 部分移到 17_t15_paml_run.py：68 个物种自由估 135 个枝长时，M7/M8/分支-位点每个要跑数小时（实测 45 分钟无一完成），
# 改为先 M0 估枝长、其余模型固定枝长（fix_blength = 2）。
P('密码子比对与树已写出 → 接着跑 17_t15_paml_run.py')
open(os.path.join(WD, 'prep_log.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
