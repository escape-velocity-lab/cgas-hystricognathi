"""06（任务 #13）: 拉 RefSeq 直系同源集合以外的 cGAS 蛋白——啮齿目（txid9989）里的 GenBank / 非 RefSeq 记录，
外加豪猪亚目的 RefSeq 全部异构体，给 13a 第二来源核对和 13b 扩采样用。只下蛋白，不下基因组。
用法: python 06_t13_fetch_extra.py <data/raw/t13>
"""
import sys, os, json, time, urllib.request, urllib.parse

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
op = urllib.request.build_opener(urllib.request.ProxyHandler())

def call(ep, **kw):
    data = urllib.parse.urlencode(kw).encode()
    for k in range(4):
        try:
            return op.open(E + ep, data=data, timeout=60).read()
        except Exception as e:
            print('retry', ep, e); time.sleep(3)
    raise RuntimeError(ep)

Q = '(cGAS[All Fields] OR "cyclic GMP-AMP synthase"[All Fields] OR MB21D1[All Fields] OR "Mab-21 domain containing 1"[All Fields])'
for name, term in (('rodentia_all', f'txid9989[Organism:exp] AND {Q}'),):
    r = json.loads(call('esearch.fcgi', db='protein', term=term, retmax=2000, retmode='json'))
    ids = r['esearchresult']['idlist']
    print(name, 'count', r['esearchresult']['count'], 'got', len(ids))
    summ = {}
    for i in range(0, len(ids), 200):
        s = json.loads(call('esummary.fcgi', db='protein', id=','.join(ids[i:i + 200]), retmode='json'))['result']
        summ.update({u: s[u] for u in s['uids']})
    json.dump(summ, open(os.path.join(OUT, f'ncbi_{name}_summary.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=0)
    fa = b''
    for i in range(0, len(ids), 200):
        fa += call('efetch.fcgi', db='protein', id=','.join(ids[i:i + 200]), rettype='fasta', retmode='text')
    open(os.path.join(OUT, f'ncbi_{name}.faa'), 'wb').write(fa)
    print('fasta records', fa.count(b'>'))
