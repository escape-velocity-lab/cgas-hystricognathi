"""09（任务 #13）: 解析 07 的 tblastn XML，读每个组装在人 cGAS 431/479/495/498 上的残基。
- 查询 = 人 Q8N884 第 400~522 位，HSP 的 query 坐标 + 399 = 人编号
- 每个位点：取覆盖它的 HSP 里 bitscore 最高的那条（优先与全局最佳 contig 同一条 contig）
- 记录 contig、HSP 同一性、位点 ±15 aa 窗口（来自 HSP 的 hit 串，缺口 '-'，跨 HSP 断开处为空）
用法: python 09_t13_parse_tblastn.py <tblastn_dir> <out_csv>
"""
import sys, os, glob, re
import xml.etree.ElementTree as ET
import pandas as pd

D, OUT = sys.argv[1:3]
Q0 = 400
SITES = [431, 479, 495, 498]
NMR = dict(zip(SITES, 'SEYT'))
# 物种 / 组装 / 用途 直接从 07 的 JOBS 表读（不 import，07 一 import 就会提交 BLAST）
src7 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '07_t13_tblastn.py'), encoding='utf-8').read()
JOBS = {db: (sp, asm, why) for sp, db, asm, why in re.findall(r"\('([^']+)', '([A-Z0-9]+)', '(GCA_[0-9.]+)', '(13[ab])'\)", src7)}

rows = []
for f in sorted(glob.glob(os.path.join(D, '*.xml'))):
    db = os.path.basename(f)[:-4]
    if os.path.getsize(f) < 500:
        rows.append({"物种": JOBS[db][0], "组装": JOBS[db][1], "用途": JOBS[db][2], "WGS": db, "状态": "BLAST 返回空结果"}); continue
    root = ET.parse(f).getroot()
    hsps = []
    for hit in root.iter('Hit'):
        acc = hit.findtext('Hit_accession')
        defl = hit.findtext('Hit_def')
        for h in hit.iter('Hsp'):
            hsps.append(dict(acc=acc, defl=defl, bits=float(h.findtext('Hsp_bit-score')), ev=float(h.findtext('Hsp_evalue')),
                             qf=int(h.findtext('Hsp_query-from')), qt=int(h.findtext('Hsp_query-to')),
                             hf=int(h.findtext('Hsp_hit-from')), ht=int(h.findtext('Hsp_hit-to')), frame=h.findtext('Hsp_hit-frame'),
                             ident=int(h.findtext('Hsp_identity')), alen=int(h.findtext('Hsp_align-len')),
                             q=h.findtext('Hsp_qseq'), s=h.findtext('Hsp_hseq')))
    # 旁系同源过滤：cGAS 真基因的 HSP 同一性 ≥ 0.58（本批 54/93 起），MB21D2 一类远缘旁系同源只有 ~0.27；阈值取 0.45
    hsps = [h for h in hsps if h['ident'] / h['alen'] >= 0.45]
    if not hsps:
        rows.append({"物种": JOBS[db][0], "组装": JOBS[db][1], "用途": JOBS[db][2], "WGS": db, "状态": "无命中"}); continue
    best_contig = max(hsps, key=lambda x: x['bits'])['acc']
    def hmap(h):
        m, qp = {}, h['qf'] + Q0 - 1
        for a, b in zip(h['q'], h['s']):
            if a != '-':
                m[qp] = b; qp += 1
        return m
    rec = {"物种": JOBS[db][0], "组装": JOBS[db][1], "用途": JOBS[db][2], "WGS": db, "状态": "ok", '最佳contig': best_contig, 'HSP数': len(hsps)}
    for s in SITES:
        cov = [h for h in hsps if h['qf'] + Q0 - 1 <= s <= h['qt'] + Q0 - 1]
        if not cov:
            rec[str(s)] = '?'; rec[f'{s}_来源'] = '无HSP覆盖'; continue
        cov.sort(key=lambda h: (h['acc'] != best_contig, -h['bits']))
        h = cov[0]; m = hmap(h)
        rec[str(s)] = m.get(s, '?')
        rec[f'{s}_来源'] = f"{h['acc']}:{h['hf']}-{h['ht']}({h['frame']}) id={h['ident']}/{h['alen']} bits={h['bits']:.0f}"
        rec[f'{s}_窗口'] = ''.join(m.get(p, ' ') for p in range(s - 15, s + 16))
        alt = {hmap(x).get(s) for x in cov[1:]} - {None, rec[str(s)]}
        rec[f'{s}_其他HSP残基'] = ''.join(sorted(alt))
    rec['命中'] = sum(rec.get(str(s)) == NMR[s] for s in SITES)
    rows.append(rec)
df = pd.DataFrame(rows)
df.to_csv(OUT, index=False, encoding='utf-8-sig')
pd.set_option('display.width', 250)
print(df[['物种', '用途', 'WGS', '状态', '最佳contig', 'HSP数'] + [str(s) for s in SITES] + ['命中']].to_string())
for _, r in df.iterrows():
    if r['状态'] == 'ok':
        print(r.WGS, ' | '.join(f"{s}:{r.get(f'{s}_窗口','')} [{r.get(f'{s}_其他HSP残基','')}] {r.get(f'{s}_来源','')}" for s in SITES))
