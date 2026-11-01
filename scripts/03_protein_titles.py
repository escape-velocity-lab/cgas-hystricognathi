"""03: 取每条蛋白记录的完整标题（NCBI esummary），用来识别 LOW QUALITY PROTEIN / partial 记录。
datasets 的 protein.faa 表头不带这些字样，必须另取。
用法: python 03_protein_titles.py <protein.faa> <输出 tsv>
"""
import sys, json, time, urllib.request, urllib.parse

FAA, OUT = sys.argv[1], sys.argv[2]
PROXY = urllib.request.ProxyHandler()
OPEN = urllib.request.build_opener(PROXY)
accs = [l[1:].split()[0] for l in open(FAA) if l.startswith('>')]
rows = []
for i in range(0, len(accs), 150):
    chunk = accs[i:i + 150]
    body = urllib.parse.urlencode({'db': 'protein', 'id': ','.join(chunk), 'retmode': 'json'}).encode()
    for k in range(5):
        try:
            r = json.load(OPEN.open('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi', data=body, timeout=120))
            break
        except Exception as e:
            time.sleep(3 + 5 * k)
    else:
        raise RuntimeError('esummary failed')
    res = r['result']
    for uid in res['uids']:
        x = res[uid]
        rows.append((x['accessionversion'], x['title'], str(x.get('slen', '')), x.get('createdate', ''), x.get('updatedate', '')))
    time.sleep(0.5)
with open(OUT, 'w', encoding='utf-8') as f:
    f.write('accession\ttitle\tslen\tcreatedate\tupdatedate\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')
got = {r[0] for r in rows}
print('requested', len(accs), 'got', len(rows), 'missing', [a for a in accs if a not in got][:10])
