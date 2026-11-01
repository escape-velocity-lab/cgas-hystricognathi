"""26（任务 #23b）: 解析 22_t23_bat_tblastn.py 的 XML，读每个狐蝠科组装在人 cGAS 431/479/495/498 上的残基。
口径与任务 #13 的 09_t13_parse_tblastn.py 完全一致：
- 查询 = 人 Q8N884 第 400~522 位，HSP 的 query 坐标 + 399 = 人编号
- 旁系同源过滤：同一性 >= 0.45（cGAS 真基因 ~0.6-0.9，MB21D2 一类远缘旁系同源只有 ~0.27）
- 每个位点取覆盖它的 HSP 里 bitscore 最高的那条（优先与全局最佳 contig 同一条）
- 记录 contig、HSP 同一性、位点 ±15 aa 窗口
有 RefSeq 注释的 4 个物种是阳性对照：tblastn 读出来必须和注释蛋白一致（已知 479 位 = H）。
输出：results/t23_bat_sites.csv / t23_bat_sites.txt
用法: python 26_t23_parse_bats.py <proj_root>
"""
import sys, os, glob, re
import xml.etree.ElementTree as ET
import pandas as pd

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
D = R('data', 'raw', 't23', 'tblastn_bats')
Q0 = 400
SITES = [431, 479, 495, 498]
NMR = dict(zip(SITES, 'SEYT'))        # 裸鼹鼠型残基
src22 = open(R('src', '22_t23_bat_tblastn.py'), encoding='utf-8').read()
JOBS = {db: (sp, asm, why) for sp, db, asm, why in
        re.findall(r"\('([^']+)',\s*'([A-Z0-9]+)',\s*'(GCA_[0-9.]+)',\s*'([a-z+]+)'\)", src22)}

rows = []
for db, (sp, asm, why) in sorted(JOBS.items(), key=lambda x: x[1][0]):
    f = os.path.join(D, db + '.xml')
    base = {'物种': sp, '组装': asm, 'WGS': db, '用途': why}
    if not os.path.exists(f) or os.path.getsize(f) < 500:
        rows.append({**base, '状态': '未完成/空结果'}); continue
    root = ET.parse(f).getroot()
    hsps = []
    for hit in root.iter('Hit'):
        acc = hit.findtext('Hit_accession')
        for h in hit.iter('Hsp'):
            hsps.append(dict(acc=acc, bits=float(h.findtext('Hsp_bit-score')), ev=float(h.findtext('Hsp_evalue')),
                             qf=int(h.findtext('Hsp_query-from')), qt=int(h.findtext('Hsp_query-to')),
                             hf=int(h.findtext('Hsp_hit-from')), ht=int(h.findtext('Hsp_hit-to')),
                             ident=int(h.findtext('Hsp_identity')), alen=int(h.findtext('Hsp_align-len')),
                             q=h.findtext('Hsp_qseq'), s=h.findtext('Hsp_hseq')))
    hsps = [h for h in hsps if h['ident'] / h['alen'] >= 0.45]
    if not hsps:
        rows.append({**base, '状态': '无命中（过滤后）'}); continue
    best = max(hsps, key=lambda x: x['bits'])
    def hmap(h):
        m, qp = {}, h['qf'] + Q0 - 1
        for a, b in zip(h['q'], h['s']):
            if a != '-':
                m[qp] = b; qp += 1
        return m
    maps = [(h, hmap(h)) for h in hsps]
    row = {**base, '状态': 'ok', '最佳contig': best['acc'], 'HSP数': len(hsps),
           '最佳HSP同一性': '%d/%d' % (best['ident'], best['alen']), 'bits': best['bits']}
    for p in SITES:
        cand = [(h, m) for h, m in maps if p in m and m[p] != '-']
        if not cand:
            row[str(p)] = '-'; row[f'{p}_窗口'] = ''; row[f'{p}_来源'] = '无覆盖'; continue
        cand.sort(key=lambda x: (x[0]['acc'] == best['acc'], x[0]['bits']), reverse=True)
        h, m = cand[0]
        row[str(p)] = m[p]
        row[f'{p}_窗口'] = ''.join(m.get(q, '.') for q in range(p - 15, p + 16))
        row[f'{p}_来源'] = '%s id=%d/%d bits=%.0f' % (h['acc'], h['ident'], h['alen'], h['bits'])
        # 其他 HSP 给出的不同读数（旁系同源 / 组装歧义的警报）
        other = sorted({m2[p] for h2, m2 in cand[1:] if m2.get(p, '-') != '-'} - {m[p]})
        row[f'{p}_其他HSP残基'] = ','.join(other)
    row['裸鼹鼠型残基数'] = sum(1 for p in SITES if row.get(str(p)) == NMR[p])
    rows.append(row)

d = pd.DataFrame(rows)
d.to_csv(R('results', 't23_bat_sites.csv'), index=False, encoding='utf-8-sig')
L = []
P = lambda *a: (print(*a, flush=True), L.append(' '.join(str(x) for x in a)))
P('狐蝠科（Old World fruit bats）cGAS 四位点 tblastn 读数')
P('裸鼹鼠型残基 = 431:S  479:E  495:Y  498:T；人 = D/K/L/K')
P('')
cols = ['物种', 'WGS', '用途', '状态', '431', '479', '495', '498', '裸鼹鼠型残基数', '最佳HSP同一性']
P(d.reindex(columns=cols).to_string(index=False))
P('')
ok = d[d['状态'] == 'ok']
if len(ok):
    P('读出 479 位的物种数: %d；各残基计数: %s' % (ok['479'].notna().sum(), ok['479'].value_counts().to_dict()))
    P('携带裸鼹鼠型 E479 的物种: %s' % (list(ok.loc[ok['479'] == 'E', '物种']) or '无'))
    P('四个位点全是裸鼹鼠型的物种: %s' % (list(ok.loc[ok['裸鼹鼠型残基数'] == 4, '物种']) or '无'))
    P('')
    P('阳性对照（有 RefSeq 注释，已知 479 = H）:')
    P(ok[ok['用途'].str.contains('ctrl')][['物种', '431', '479', '495', '498']].to_string(index=False))
    P('')
    P('479 位 ±15 窗口:')
    for _, r in ok.iterrows():
        P('  %-30s %s' % (r['物种'], r['479_窗口']))
open(R('results', 't23_bat_sites.txt'), 'w', encoding='utf-8').write('\n'.join(L))
P('')
P('写出 results/t23_bat_sites.csv / .txt')
