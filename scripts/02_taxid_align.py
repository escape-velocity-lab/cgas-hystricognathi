"""02: AnAge 与 NCBI 直系同源按 taxid 对齐（坑 42：不按名字 join）。
1) AnAge 全部哺乳动物学名 → NCBI Taxonomy 解析成 taxid（会把同物异名解析到现名）
2) 直系同源物种 taxid → 物种级 taxid（亚种上卷到种）+ 目 / 科
3) 按物种级 taxid join；名字被 NCBI 改写的、多对一的、对不上的全部列进手工核对清单
4) 手工覆盖表 MANUAL 里每一条都写了理由
用法: python 02_taxid_align.py <data/raw> <data/processed>
"""
import sys, os, json, time, urllib.request, urllib.parse
import pandas as pd

RAW, PROC = sys.argv[1], sys.argv[2]
os.makedirs(PROC, exist_ok=True)
PROXY = urllib.request.ProxyHandler()
OPEN = urllib.request.build_opener(PROXY)
API = 'https://api.ncbi.nlm.nih.gov/datasets/v2'
CACHE = os.path.join(RAW, 'ncbi', 'taxonomy_cache.json')
cache = json.load(open(CACHE, encoding='utf-8')) if os.path.exists(CACHE) else {}

def get(url):
    err = None
    for k in range(5):
        try:
            return json.load(OPEN.open(url, timeout=120))
        except Exception as e:
            err = e; time.sleep(3 + 5 * k)
    raise err

def tax_reports(keys):
    """keys: 名字或 taxid（字符串）。结果按 query 缓存。"""
    todo = [k for k in keys if k not in cache]
    for i in range(0, len(todo), 40):
        chunk = todo[i:i + 40]
        r = get(f'{API}/taxonomy/taxon/' + urllib.parse.quote(','.join(chunk), safe=',') + '/dataset_report')
        for x in r.get('reports', []):
            for q in x.get('query', []):
                cache[q] = x
        for q in chunk:
            cache.setdefault(q, {'errors': 'no report'})
        time.sleep(0.4)
        print('taxonomy', i + len(chunk), '/', len(todo))
    json.dump(cache, open(CACHE, 'w', encoding='utf-8'), ensure_ascii=False)
    return {k: cache[k] for k in keys}

def cls(t, rank):
    return t.get('classification', {}).get(rank, {}).get('name', '')

# ---- AnAge ----
a = pd.read_csv(os.path.join(RAW, 'anage', 'anage_data.txt'), sep='\t')
a = a[a.Class == 'Mammalia'].copy()
a['anage_name'] = a.Genus + ' ' + a.Species
R = tax_reports(list(a.anage_name))
a['anage_taxid'] = [str(R[n].get('taxonomy', {}).get('tax_id', '')) for n in a.anage_name]
a['ncbi_name_of_anage'] = [R[n].get('taxonomy', {}).get('current_scientific_name', {}).get('name', '') for n in a.anage_name]

# 手工覆盖（AnAge 学名 → 正确的 NCBI 物种 taxid）。每条都要有理由。
MANUAL = {
    # 坑 42：NCBI 把 P. giganteus 当成 P. vampyrus（马来狐蝠）的异名；AnAge 的 P. giganteus（44 年）是印度狐蝠，NCBI 现名 P. medius
    'Pteropus giganteus': ('Pteropus medius', 'AnAge 的印度狐蝠；NCBI 把 P. giganteus 解析为 P. vampyrus，是错的'),
}
for an, (nc, why) in MANUAL.items():
    t = tax_reports([nc])[nc]['taxonomy']
    a.loc[a.anage_name == an, 'anage_taxid'] = str(t['tax_id'])
    a.loc[a.anage_name == an, 'ncbi_name_of_anage'] = t['current_scientific_name']['name'] + '（手工）'

# ---- 直系同源物种 ----
orth = [json.loads(l) for l in open(os.path.join(RAW, 'ncbi', 'cgas_orthologs_mammalia.jsonl'), encoding='utf-8')]
o = pd.DataFrame({'gene_id': [x['gene_id'] for x in orth], 'symbol': [x['symbol'] for x in orth],
                  'tax_id': [x['tax_id'] for x in orth], 'taxname': [x['taxname'] for x in orth],
                  'common_name': [x.get('common_name', '') for x in orth]})
T = tax_reports(list(o.tax_id))
def species_tid(tid):
    t = T[tid]['taxonomy']
    if t.get('rank') == 'SPECIES':
        return str(t['tax_id'])
    sp = t.get('classification', {}).get('species', {})
    return str(sp.get('id', t['tax_id']))
o['species_taxid'] = [species_tid(t) for t in o.tax_id]
o['rank'] = [T[t]['taxonomy'].get('rank', '') for t in o.tax_id]
o['order'] = [cls(T[t]['taxonomy'], 'order') for t in o.tax_id]
o['family'] = [cls(T[t]['taxonomy'], 'family') for t in o.tax_id]
o['genus'] = [cls(T[t]['taxonomy'], 'genus') for t in o.tax_id]
o['lineage'] = [','.join(str(p) for p in T[t]['taxonomy'].get('parents', [])) for t in o.tax_id]

# ---- join ----
keep = ['anage_name', 'anage_taxid', 'ncbi_name_of_anage', 'Order', 'Family', 'Common name',
        'Maximum longevity (yrs)', 'Adult weight (g)', 'Data quality', 'Sample size', 'Specimen origin']
aa = a[keep].copy()
dup = aa[aa.anage_taxid != ''].groupby('anage_taxid').anage_name.apply(list)
dup = dup[dup.map(len) > 1]
m = o.merge(aa[aa.anage_taxid != ''], left_on='species_taxid', right_on='anage_taxid', how='left')

# 手工核对清单
rows = []
for _, r in m.iterrows():
    sp2 = ' '.join(r.taxname.split()[:2])
    if pd.isna(r.anage_name):
        # 同名但 taxid 没对上的？（名字 join 能对上而 taxid 不能 → 必须人工看）
        hit = a[a.anage_name == sp2]
        rows.append(dict(taxname=r.taxname, tax_id=r.tax_id, species_taxid=r.species_taxid, 问题='AnAge 无对应 taxid',
                         AnAge同名记录=('; '.join(hit.anage_name + ' → taxid ' + hit.anage_taxid)) if len(hit) else '',
                         处理='不进寿命分析（位点表仍保留）'))
    elif r.anage_name != sp2:
        rows.append(dict(taxname=r.taxname, tax_id=r.tax_id, species_taxid=r.species_taxid,
                         问题=f'学名不同：AnAge「{r.anage_name}」→ NCBI「{r.ncbi_name_of_anage}」',
                         AnAge同名记录='', 处理='taxid 一致，采用（异名 / 亚种上卷）'))
for tid, names in dup.items():
    rows.append(dict(taxname='', tax_id='', species_taxid=tid, 问题='AnAge 多条记录解析到同一 taxid：' + ' / '.join(names),
                     AnAge同名记录='', 处理='见 site_table 备注'))
chk = pd.DataFrame(rows)
m.to_csv(os.path.join(PROC, 'orthologs_x_anage.csv'), index=False, encoding='utf-8-sig')
chk.to_csv(os.path.join(PROC, 'taxid_manual_check.csv'), index=False, encoding='utf-8-sig')
print('orthologs', len(o), 'unique species_taxid', o.species_taxid.nunique())
print('matched to AnAge', m.anage_name.notna().sum(), '| with longevity+weight',
      (m['Maximum longevity (yrs)'].notna() & m['Adult weight (g)'].notna()).sum())
print('AnAge mammals', len(a), 'unresolved names', (a.anage_taxid == '').sum(), 'dup taxids', len(dup))
print(chk.问题.str[:12].value_counts().to_string())
