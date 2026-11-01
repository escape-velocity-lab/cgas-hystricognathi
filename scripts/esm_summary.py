"""E03 描述统计：读 esm_score.py 的输出，汇总成 esm_v1_summary.txt（只做描述；物种间不独立，不报 p 值）
用法: python esm_summary.py <项目02根目录> [tag]
"""
import sys, os
import numpy as np
import pandas as pd

ROOT = sys.argv[1]; TAG = sys.argv[2] if len(sys.argv) > 2 else 'v1'
RES = os.path.join(ROOT, 'results')
SITES = ['D431', 'K479', 'L495', 'K498']
NMR = dict(zip(SITES, 'SEYT')); HUM = dict(zip(SITES, 'DKLK'))
out = []
def p(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

def rho(x, y):
    d = pd.DataFrame({'x': x, 'y': y}).dropna()
    if len(d) < 5:
        return np.nan, len(d)
    return round(np.corrcoef(d.x.rank(), d.y.rank())[0, 1], 3), len(d)

# ---- 人背景：裸鼹鼠型在「自然出现过的替换」里的位置（E3）----
aa = pd.read_csv(os.path.join(RES, f'esm_{TAG}_human_site_aa.csv'))
p('## 人背景：4 个位点上裸鼹鼠型 LLR 在该位点「哺乳动物里出现过的非人型残基」中的位置')
for s in SITES:
    d = aa[aa.人位点 == s]
    nat = d[(d.哺乳动物中出现的物种数 > 0) & ~d.是人型]
    allx = d[~d.是人型]
    v = float(d[d.是裸鼹鼠型].LLR_相对人型.iloc[0])
    p(f'{s}{NMR[s]}: LLR={v:+.2f}；人型秩={int(d[d.是人型].秩.iloc[0])}，裸鼹鼠型秩={int(d[d.是裸鼹鼠型].秩.iloc[0])}；'
      f'自然出现过的非人型残基 {len(nat)} 种，LLR 中位数 {nat.LLR_相对人型.median():+.2f}，范围 [{nat.LLR_相对人型.min():+.2f}, {nat.LLR_相对人型.max():+.2f}]，'
      f'裸鼹鼠型分位 {(nat.LLR_相对人型 < v).mean():.2f}；19 种非人型全体中分位 {(allx.LLR_相对人型 < v).mean():.2f}')
    top = d.sort_values('秩').head(5)
    p('   模型最偏好的前 5：' + ', '.join(f'{r.残基}({r.LLR_相对人型:+.2f}, 自然 {r.哺乳动物中出现的物种数} 种)' for r in top.itertuples()))

# ---- 物种 ----
sp = pd.read_csv(os.path.join(RES, f'esm_{TAG}_species.csv'))
u = sp[(sp['4位点可用'] == '是') & sp.PLL_per_res.notna()].copy()
u['log10寿命'] = np.log10(u.最大寿命_年)
p(f'\n## 物种：可用且已打分 {len(u)} 种，其中有 LQ 的 {u.LQ.notna().sum()} 种')
# 三种口径：整条序列；只取比对到人 161~522（催化域）/ 421~508（4 位点所在窗口）的位置。整条序列口径受 N 端长度 / 重复片段干扰
PLLS = {'整条': 'PLL_per_res', '催化域161_522': 'PLL_per_res_催化域161_522', '位点窗421_508': 'PLL_per_res_位点窗421_508'}
for nm, c in PLLS.items():
    p(f'\n### PLL/L 口径：{nm}（Spearman ρ, n）')
    p('分布：', u[c].describe().round(4).to_dict())
    p('与蛋白长度 ρ =', rho(u[c], u.蛋白长度))
    g = u.groupby('目')[c].agg(['count', 'median']).query('count >= 5').sort_values('median', ascending=False)
    p('分目中位数（n≥5）：' + '; '.join(f'{o} {r["median"]:.3f}(n={int(r["count"])})' for o, r in g.iterrows()))
    p('全体：LQ', rho(u[c], u.LQ), '| log10 最大寿命', rho(u[c], u.log10寿命), '| log10 体重', rho(u[c], np.log10(u.成年体重_g)))
    for o in ['Primates', 'Rodentia', 'Chiroptera', 'Artiodactyla', 'Carnivora']:
        d = u[u.目 == o]
        p(f'  {o:14s} 目内：LQ', rho(d[c], d.LQ), '| log10 寿命', rho(d[c], d.log10寿命), '| 长度', rho(d[c], d.蛋白长度))
    p('按寿命分组中位数：', u.groupby('寿命分组')[c].median().round(4).to_dict(), ' n=', u.groupby('寿命分组').size().to_dict())
    p('最高 8 种：' + '; '.join(f'{r.物种}({r.目}, L={r.蛋白长度}, {getattr(r, c):.3f})' for r in u.sort_values(c, ascending=False).head(8).itertuples()))

p('\n### S_nmr = Σ[log p(裸鼹鼠型) − log p(人型)]（在各物种自己的序列上，逐位遮住）')
p('S_nmr 分布：', u.S_nmr.describe().round(3).to_dict())
p('S_nmr 与命中数 ρ =', rho(u.S_nmr, u.裸鼹鼠型命中数))
p('按命中数 S_nmr 中位数：', u.groupby('裸鼹鼠型命中数').S_nmr.median().round(3).to_dict(), ' n=', u.groupby('裸鼹鼠型命中数').size().to_dict())
nh = u[~u.豪猪亚目]
p('S_nmr 与 LQ：全体', rho(u.S_nmr, u.LQ), '| 豪猪亚目以外', rho(nh.S_nmr, nh.LQ), '| 豪猪亚目以外 vs log10 寿命', rho(nh.S_nmr, nh.log10寿命))
for s in SITES:
    c = f'{s}_nmr减人'
    p(f'  {s} 单位点 nmr−人 与 LQ（豪猪亚目以外）', rho(nh[c], nh.LQ), f'；全体中位数 {u[c].median():+.2f}，>0 的物种比例 {(u[c] > 0).mean():.2f}')
p('按寿命分组 S_nmr 中位数：', u.groupby('寿命分组').S_nmr.median().round(3).to_dict())

p('\n### 模型对各物种「自己现有残基」的评价（在自己序列上遮住该位）')
for s in SITES:
    c = f'{s}_自身秩'
    p(f'{s}: 自身残基秩 中位数 {u[c].median():.0f}，秩=1 的比例 {(u[c] == 1).mean():.2f}，秩 ≤3 的比例 {(u[c] <= 3).mean():.2f}')

cols = ['物种', '目', '科', 'LQ', '裸鼹鼠型命中数', 'D431_残基', 'K479_残基', 'L495_残基', 'K498_残基', 'S_nmr', 'PLL_per_res', 'PLL_per_res_催化域161_522', 'PLL_per_res_位点窗421_508']
p('\n### S_nmr 最高 12 种')
p(u.sort_values('S_nmr', ascending=False)[cols].head(12).to_string(index=False))
p('\n### S_nmr 最低 8 种')
p(u.sort_values('S_nmr')[cols].head(8).to_string(index=False))
p('\n### 豪猪亚目 + 对照')
p(u[u.豪猪亚目 | u.物种.isin(['Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Nannospalax galili', 'Myotis brandtii', 'Balaena mysticetus'])][cols].to_string(index=False))
p('\n### LQ 最高 12 种')
p(u.sort_values('LQ', ascending=False)[cols].head(12).to_string(index=False))
skipped = sp[sp.PLL_per_res.isna()]
p(f'\n未打分 / 不可用：{len(sp) - len(u)} 种；其中未打分 {len(skipped)}：', skipped[['物种', '蛋白长度']].to_dict('records'))

open(os.path.join(RES, f'esm_{TAG}_summary.txt'), 'w', encoding='utf-8').write('\n'.join(out) + '\n')
