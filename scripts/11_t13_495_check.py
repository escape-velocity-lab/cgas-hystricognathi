"""11（任务 #13 · 13c 补充）: 495 位祖先状态对模型敏感（多态 ER 说 H，二态 ARD 说 Y），这里用两种不依赖速率模型的办法定：
1. 约束简约（Sankoff，单位代价）：把豪猪亚目祖先固定为 H / Y / N，数豪猪亚目内部最少要几次变化、Y 要独立出现几次
2. 节点年龄（Upham 2019 Completed 树）：豪猪亚目冠群、非洲鼹鼠科冠群、梳趾鼠分歧
用法: python 11_t13_495_check.py <proj_root>
"""
import sys, os, importlib.util
ROOT = sys.argv[1]
# 复用 10 的建树和数据（10 跑一遍约 1 分钟；这里只借它的函数，不重跑重建）
spec = importlib.util.spec_from_file_location('asr', os.path.join(ROOT, 'src', '10_t13_asr.py'))
src = open(os.path.join(ROOT, 'src', '10_t13_asr.py'), encoding='utf-8').read()
src = src.split('# ---------- 5. 跑 ----------')[0]
g = {'__name__': 'asr', 'sys': sys}
sys.argv = [sys.argv[0], ROOT]
exec(compile(src, '10_t13_asr.py', 'exec'), g)
df, build, mrca, postorder, tips_of = g['df'], g['build'], g['mrca'], g['postorder'], g['tips_of']
R = g['R']

root, miss = build(R('data', 'raw', 't13', 'tree', 'MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre'))
hy = mrca(root, ['Cavia porcellus', 'Heterocephalus glaber', 'Hystrix cristata'])
st = df.set_index('物种')['495'].to_dict()
states = sorted({st[t] for t in tips_of(hy) if st[t] not in '?-X'})

def sankoff(n):
    """返回 dict 状态 → 最小代价；同时记录回溯"""
    if not n.ch:
        v = st[n.name]
        return {a: (0 if (v == a or v in '?-X') else float('inf')) for a in states}
    cs = [sankoff(c) for c in n.ch]
    return {a: sum(min(c[b] + (a != b) for b in states) for c in cs) for a in states}

cost = sankoff(hy)
print('495 位，豪猪亚目内部（%d 种）约束简约最少变化次数：' % len(tips_of(hy)))
for a in states:
    print(f'  祖先固定为 {a}: {cost[a]}')

# 在「祖先 = H」和「祖先 = Y」下，用 Fitch 风格的回溯数 Y 的独立获得次数（取一个最简约重建）
def assign(n, a, out):
    out[n] = a
    for c in n.ch:
        cc = sankoff(c)
        b = min(states, key=lambda b: (cc[b] + (a != b), b != a))
        assign(c, b, out)
    return out
for a in ('H', 'Y'):
    if a not in states:
        continue
    asg = assign(hy, a, {})
    gains = [n for n in postorder(hy, []) for c in n.ch if asg[c] == 'Y' and asg[n] != 'Y' for n in [c]]
    losses = [c for n in postorder(hy, []) for c in n.ch if asg[n] == 'Y' and asg[c] != 'Y']
    print(f'  祖先={a} 的一个最简约重建：Y 获得 {len(gains)} 次 → ' + '; '.join('{' + ', '.join(sorted(tips_of(x))) + '}' for x in gains)
          + f'；Y 丢失 {len(losses)} 次 → ' + '; '.join('{' + ', '.join(sorted(tips_of(x))) + '}' for x in losses))

def height(n):
    h, x = 0.0, n
    while x.ch:
        x = x.ch[0]; h += x.bl
    return h
for lab, sps in (('梳趾鼠 vs 豪猪亚目分歧', ['Ctenodactylus gundi', 'Cavia porcellus']),
                 ('豪猪亚目冠群', ['Cavia porcellus', 'Heterocephalus glaber', 'Hystrix cristata']),
                 ('非洲鼹鼠科 s.l. 冠群（裸鼹鼠 vs Fukomys）', ['Heterocephalus glaber', 'Fukomys damarensis']),
                 ('豚鼠属 vs 豚鼠科其他（Cavia vs 水豚/长耳豚鼠）', ['Cavia porcellus', 'Hydrochoerus hydrochaeris']),
                 ('啮齿目冠群', ['Mus musculus', 'Cavia porcellus', 'Sciurus vulgaris'])):
    print(f'  {lab}: 约 {height(mrca(root, sps)):.1f} 百万年前')
