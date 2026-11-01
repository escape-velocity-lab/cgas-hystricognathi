"""22（任务 #23b）: 用 NCBI 远程 tblastn 在狐蝠科（Old World fruit bats）WGS 组装里读 cGAS 4 个位点。
目的：补上稿件里「我们只查了 4 种有注释的果蝠」这个缺口，特别是 Tian et al. 2023 分析过的 6 种果蝠。
覆盖 NCBI 上有组装的全部 13 个狐蝠科物种（15 套组装），其中 4 个有 RefSeq 注释，用作阳性对照
（已知人 479 位 = H，tblastn 读出来必须一致，否则方法有问题）。
查询序列 = 人 cGAS Q8N884 第 400~522 位（覆盖 431/479/495/498），与任务 #13 完全相同。
不下载基因组，只拿 BLAST 结果 XML。代理用 10090（17890 在 2026-10-21 不通）。
用法: python 22_t23_bat_tblastn.py <data/raw> <out_dir>
"""
import sys, os, time, re, urllib.request, urllib.parse
from Bio import SeqIO

RAW, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
PROXY = os.environ.get('T23_PROXY') or None
op = urllib.request.build_opener(urllib.request.ProxyHandler({'http': PROXY, 'https': PROXY} if PROXY else None))
URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
H = ''.join(str(r.seq) for r in SeqIO.parse(os.path.join(RAW, 'uniprot', 'Q8N884.fasta'), 'fasta'))
Q0 = 400
QUERY = H[Q0 - 1:522]

# (物种, WGS 项目, 组装, 用途)  tian = Tian et al. 2023 分析过的 6 种；ctrl = 有 RefSeq 注释，阳性对照
JOBS = [
    ('Cynopterus sphinx',          'JAKZJJ01', 'GCA_030015415.1', 'tian'),   # Tian 2023 新测的参考基因组
    ('Cynopterus brachyotis',      'SSHV01',   'GCA_009793145.1', 'tian'),
    ('Macroglossus sobrinus',      'PVKZ01',   'GCA_004027375.1', 'tian'),
    ('Pteropus alecto',            'ALWS01',   'GCA_000325575.1', 'tian+ctrl'),
    ('Pteropus vampyrus',          'ABRP02',   'GCA_000151845.2', 'tian+ctrl'),
    ('Rousettus aegyptiacus',      'JACASE02', 'GCA_014176215.2', 'tian+ctrl'),
    ('Pteropus medius',            'CACVBW01', 'GCA_902729225.1', 'ctrl'),
    ('Eidolon helvum',             'AWHC01',   'GCA_000465285.1', 'gap'),    # 稿件里点名的缺口：无 cGAS 注释
    ('Eidolon dupreanum',          'JCBABP01', 'GCA_059749105.1', 'new'),
    ('Eonycteris spelaea',         'JCELFR01', 'GCA_062009626.1', 'new'),
    ('Hypsignathus monstrosus',    'CFOBSM01', 'GCA_985612075.1', 'new'),
    ('Pteropus pselaphon',         'BAAJCM01', 'GCA_056242155.1', 'new'),
    ('Pteropus rufus',             'JBWDIR01', 'GCA_056363655.1', 'new'),
    ('Rousettus leschenaultii',    'JBRKBD02', 'GCA_053538355.2', 'new'),
    ('Rousettus madagascariensis', 'JCASON01', 'GCA_059802175.1', 'new'),
]

def put(db):
    data = urllib.parse.urlencode({'CMD': 'Put', 'PROGRAM': 'tblastn', 'DATABASE': f'WGS_VDB://{db}', 'QUERY': QUERY,
                                   'HITLIST_SIZE': '10', 'EXPECT': '1e-5', 'FILTER': 'F'}).encode()
    t = op.open(URL, data=data, timeout=120).read().decode('utf-8', 'replace')
    return re.search(r'RID = (\S+)', t).group(1)

def poll(rid):
    t = op.open(f'{URL}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}', timeout=60).read().decode('utf-8', 'replace')
    m = re.search(r'Status=(\w+)', t)
    return m.group(1) if m else 'UNKNOWN'

pending = {}
for sp, db, asm, why in JOBS:
    f = os.path.join(OUT, f'{db}.xml')
    if os.path.exists(f) and os.path.getsize(f) > 500:
        print('skip (已有)', sp, db, flush=True); continue
    rid = None
    for k in range(3):
        try:
            rid = put(db); break
        except Exception as e:
            print('put retry', db, repr(e), flush=True); time.sleep(15)
    if rid is None:
        print('PUT FAILED', sp, db, flush=True); continue
    pending[db] = rid
    print(time.strftime('%H:%M:%S'), 'submitted', sp, db, rid, flush=True)
    time.sleep(11)

deadline = time.time() + 3600
while pending and time.time() < deadline:
    time.sleep(30)
    for db, rid in list(pending.items()):
        try:
            st = poll(rid)
        except Exception as e:
            print('poll err', db, repr(e), flush=True); continue
        if st == 'READY':
            x = op.open(f'{URL}?CMD=Get&FORMAT_TYPE=XML&RID={rid}', timeout=180).read()
            open(os.path.join(OUT, f'{db}.xml'), 'wb').write(x)
            print(time.strftime('%H:%M:%S'), 'ready', db, len(x), flush=True)
            del pending[db]
        elif st in ('FAILED', 'UNKNOWN'):
            print(time.strftime('%H:%M:%S'), 'FAILED', db, rid, st, flush=True)
            del pending[db]
        time.sleep(2)
print('done; still pending:', pending)
