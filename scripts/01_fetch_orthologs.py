"""01: 拉人 CGAS（Gene ID 115004）在哺乳纲（taxid 40674）里的全部 NCBI 直系同源：
   - 直系同源基因报告（jsonl，含 taxid、注释版本、组装号）
   - gene 数据包（蛋白 FASTA + CDS FASTA + product report），只按基因下，不下基因组
用法: python 01_fetch_orthologs.py <输出目录 data/raw/ncbi>
"""
import sys, os, json, time, zipfile, hashlib, urllib.request

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
PROXY = urllib.request.ProxyHandler()
OPEN = urllib.request.build_opener(PROXY)
API = 'https://api.ncbi.nlm.nih.gov/datasets/v2'

def req(url, data=None, raw=False):
    err = None
    for k in range(5):  # 坑 13：一律重试
        try:
            r = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'} if data else {})
            b = OPEN.open(r, timeout=300).read()
            return b if raw else json.loads(b)
        except Exception as e:
            err = e; time.sleep(3 + 5 * k)
    raise err

# 1) 直系同源报告（分页）
reps, tok = [], None
while True:
    u = f'{API}/gene/id/115004/orthologs?taxon_filter=40674&page_size=100' + (f'&page_token={tok}' if tok else '')
    r = req(u)
    reps += [x['gene'] for x in r.get('reports', [])]
    tok = r.get('next_page_token')
    print('orthologs so far', len(reps), 'total_count', r.get('total_count'))
    if not tok:
        break
    time.sleep(0.4)
with open(os.path.join(OUT, 'cgas_orthologs_mammalia.jsonl'), 'w', encoding='utf-8') as f:
    for x in reps:
        f.write(json.dumps(x, ensure_ascii=False) + '\n')

# 2) 按 gene id 下数据包：蛋白、CDS、product report
ids = [x['gene_id'] for x in reps]
body = json.dumps({'gene_ids': [int(i) for i in ids],
                   'include_annotation_type': ['FASTA_PROTEIN', 'FASTA_CDS'],
                   'return_report_types': ['DATASET_REPORT', 'PRODUCT_REPORT']}).encode()
z = req(f'{API}/gene/download', data=body, raw=True)
zp = os.path.join(OUT, 'cgas_gene_package.zip')
open(zp + '.part', 'wb').write(z)
os.replace(zp + '.part', zp)  # 坑 38：先写 .part 再改名
with zipfile.ZipFile(zp) as zz:
    bad = zz.testzip()  # 坑 12：CRC 校验
    print('zip members:', zz.namelist(), 'testzip:', bad)
    zz.extractall(os.path.join(OUT, 'pkg'))
print('n_gene_ids', len(ids), 'zip bytes', len(z), 'md5', hashlib.md5(z).hexdigest())
