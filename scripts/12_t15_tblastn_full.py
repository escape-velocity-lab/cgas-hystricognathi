"""12（任务 #15a）: 用 NCBI 远程 tblastn 在无注释的豪猪亚目 WGS 组装里找 cGAS **全长**（各外显子），
给 15a 的「非洲鼹鼠科有、其余豪猪亚目没有」扫描补齐物种（#13 的 07 只查了人 400~522 一段）。
查询 = 裸鼹鼠 A0AAX6RS70（Chen 2025 用的那条）全长；只拿 BLAST 结果 XML，不下基因组。
EXPECT 放宽到 1e-3、HITLIST 20，以免漏掉短外显子；解析时再按同一性 / contig 过滤。
用法: python 12_t15_tblastn_full.py <data/raw> <out_dir>
"""
import sys, os, time, re, urllib.request, urllib.parse
from Bio import SeqIO

RAW, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
op = urllib.request.build_opener(urllib.request.ProxyHandler())
URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
QUERY = str(next(SeqIO.parse(os.path.join(RAW, 'uniprot', 'A0AAX6RS70.fasta'), 'fasta')).seq)

# 与 07 相同的 13b 组装（物种, WGS）
JOBS = [
    ('Hystrix cristata', 'PVJO01'), ('Hystrix brachyura', 'QZML01'), ('Thryonomys swinderianus', 'BAAFGX01'),
    ('Petromus typicus', 'PVIR01'), ('Fukomys anselli', 'JBPETH01'), ('Cavia tschudii', 'PVKK01'),
    ('Hydrochoerus hydrochaeris', 'PVLA01'), ('Dolichotis patagonum', 'PVJX01'), ('Dasyprocta punctata', 'RJWM01'),
    ('Cuniculus paca', 'RJWT01'), ('Dinomys branickii', 'PVLD01'), ('Ctenomys sociabilis', 'PVKA01'),
    ('Myocastor coypus', 'PVJA01'), ('Capromys pilorides', 'PVKN01'), ('Octomys mimax', 'NDGM01'),
    ('Tympanoctomys barrerae', 'NDGN01'), ('Coendou prehensilis', 'JAPYYN01'), ('Cavia aperea', 'AVPZ01'),
]

def put(db):
    data = urllib.parse.urlencode({'CMD': 'Put', 'PROGRAM': 'tblastn', 'DATABASE': f'WGS_VDB://{db}', 'QUERY': QUERY,
                                   'HITLIST_SIZE': '20', 'EXPECT': '1e-3', 'FILTER': 'F'}).encode()
    t = op.open(URL, data=data, timeout=120).read().decode('utf-8', 'replace')
    return re.search(r'RID = (\S+)', t).group(1)

def poll(rid):
    t = op.open(f'{URL}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}', timeout=60).read().decode('utf-8', 'replace')
    m = re.search(r'Status=(\w+)', t)
    return m.group(1) if m else 'UNKNOWN'

pending = {}
for sp, db in JOBS:
    f = os.path.join(OUT, f'{db}.xml')
    if os.path.exists(f) and os.path.getsize(f) > 500:
        continue
    rid = None
    for k in range(3):
        try:
            rid = put(db); break
        except Exception as e:
            print('put retry', db, e, flush=True); time.sleep(15)
    if rid:
        pending[db] = rid
        print(time.strftime('%H:%M:%S'), 'submitted', sp, db, rid, flush=True)
    time.sleep(11)

while pending:
    time.sleep(30)
    for db, rid in list(pending.items()):
        try:
            st = poll(rid)
        except Exception as e:
            print('poll err', db, e, flush=True); continue
        if st == 'READY':
            x = op.open(f'{URL}?CMD=Get&FORMAT_TYPE=XML&RID={rid}', timeout=180).read()
            open(os.path.join(OUT, f'{db}.xml'), 'wb').write(x)
            print(time.strftime('%H:%M:%S'), 'ready', db, len(x), flush=True)
            del pending[db]
        elif st in ('FAILED', 'UNKNOWN'):
            print(time.strftime('%H:%M:%S'), 'FAILED', db, rid, st, flush=True)
            del pending[db]
        time.sleep(2)
print('done', flush=True)
