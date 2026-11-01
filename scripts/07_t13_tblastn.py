"""07（任务 #13）: 用 NCBI 远程 tblastn 在没有基因注释的豪猪亚目（及外群）WGS 组装里读 cGAS 4 个位点，
同时给有注释的物种找「第二个独立组装」做 13a 的第二来源核对。不下载基因组，只拿 BLAST 结果 XML。
查询序列 = 人 cGAS Q8N884 第 400~522 位（覆盖 431/479/495/498）。
用法: python 07_t13_tblastn.py <data/raw> <out_dir>
"""
import sys, os, time, re, urllib.request, urllib.parse
from Bio import SeqIO

RAW, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
op = urllib.request.build_opener(urllib.request.ProxyHandler())
URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
H = ''.join(str(r.seq) for r in SeqIO.parse(os.path.join(RAW, 'uniprot', 'Q8N884.fasta'), 'fasta'))
Q0 = 400
QUERY = H[Q0 - 1:522]

# (物种, WGS 项目, 组装, 用途)
JOBS = [
    ('Hystrix cristata', 'PVJO01', 'GCA_004026905.1', '13b'),
    ('Hystrix brachyura', 'QZML01', 'GCA_016801275.1', '13b'),
    ('Thryonomys swinderianus', 'BAAFGX01', 'GCA_041154305.1', '13b'),
    ('Petromus typicus', 'PVIR01', 'GCA_004026965.1', '13b'),
    ('Fukomys anselli', 'JBPETH01', 'GCA_053478225.1', '13b'),
    ('Cavia tschudii', 'PVKK01', 'GCA_004027695.1', '13b'),
    ('Cavia aperea', 'AVPZ01', 'GCA_000688575.1', '13b'),
    ('Hydrochoerus hydrochaeris', 'PVLA01', 'GCA_004027455.1', '13b'),
    ('Dolichotis patagonum', 'PVJX01', 'GCA_004027295.1', '13b'),
    ('Dasyprocta punctata', 'RJWM01', 'GCA_004363535.1', '13b'),
    ('Cuniculus paca', 'RJWT01', 'GCA_004365215.1', '13b'),
    ('Dinomys branickii', 'PVLD01', 'GCA_004027595.1', '13b'),
    ('Ctenomys sociabilis', 'PVKA01', 'GCA_004027165.1', '13b'),
    ('Myocastor coypus', 'PVJA01', 'GCA_004027025.1', '13b'),
    ('Capromys pilorides', 'PVKN01', 'GCA_004027915.1', '13b'),
    ('Octomys mimax', 'NDGM01', 'GCA_002564305.1', '13b'),
    ('Tympanoctomys barrerae', 'NDGN01', 'GCA_002564285.1', '13b'),
    ('Coendou prehensilis', 'JAPYYN01', 'GCA_028534165.1', '13b'),
    ('Erethizon dorsatum', 'JAQMCP01', 'GCA_028500835.1', '13a'),   # 注释蛋白 KAM6166271 来自 JAQMCO01，这里用另一套组装
    ('Ctenodactylus gundi', 'PVKB01', 'GCA_004027205.1', '13a'),    # 注释蛋白 KAM5256511 来自 JBMERT01
    ('Cavia porcellus', 'JBPSDR01', 'GCA_052324795.1', '13a'),      # 2025 新组装（RefSeq 用的是 GCF_034190915.1）
    ('Cavia porcellus', 'AAKN02', 'GCA_000151735.1', '13a'),        # 2008 旧组装（Ensembl / UniProt H0VH90 用的）
    ('Heterocephalus glaber', 'CGBMJZ01', 'GCA_988279095.1', '13a'), # 2026 染色体级新组装（BLAST 库里只有 120 条/4.3 Mb，无命中，见汇报新坑）
    ('Heterocephalus glaber', 'JBPUHJ01', 'GCA_052151175.1', '13a'), # 2025 另一个体的组装（代替上一条）
    ('Fukomys damarensis', 'JAAHWF01', 'GCA_012274545.1', '13a'),   # 2020 组装（RefSeq XP_010639479 来自 2014 的 AYUG01）
    ('Fukomys damarensis', 'AYUG01', 'GCA_000743615.1', '13a'),
    ('Chinchilla lanigera', 'AGCD01', 'GCA_000276665.1', '13a'),    # 唯一组装：查 partial 记录的 498 位
    ('Octodon degus', 'AJSA01', 'GCA_000260255.1', '13a'),          # 唯一组装
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
        continue
    for k in range(3):
        try:
            rid = put(db); break
        except Exception as e:
            print('put retry', db, e); time.sleep(15)
    pending[db] = rid
    print(time.strftime('%H:%M:%S'), 'submitted', sp, db, rid, flush=True)
    time.sleep(11)

while pending:
    time.sleep(30)
    for db, rid in list(pending.items()):
        try:
            st = poll(rid)
        except Exception as e:
            print('poll err', db, e); continue
        if st == 'READY':
            x = op.open(f'{URL}?CMD=Get&FORMAT_TYPE=XML&RID={rid}', timeout=120).read()
            open(os.path.join(OUT, f'{db}.xml'), 'wb').write(x)
            print(time.strftime('%H:%M:%S'), 'ready', db, len(x), flush=True)
            del pending[db]
        elif st in ('FAILED', 'UNKNOWN'):
            print(time.strftime('%H:%M:%S'), 'FAILED', db, rid, st, flush=True)
            del pending[db]
        time.sleep(2)
print('done')
