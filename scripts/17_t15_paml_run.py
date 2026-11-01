"""17（任务 #15c）: codeml 运行与解析（数据由 15 准备在 results/t15_paml/）。
1. M0（单一 ω）自由估枝长
2. 用 M0 的枝长（fix_blength = 2 固定）跑 M7 / M8 / M8a / 分支-位点 A / A-null，5 个并行
   ——固定枝长是大规模扫描常用的近似（参数从 ~140 个降到 ~5 个），会让 LRT 略偏乐观，报告里注明
3. 解析 lnL、LRT、rst 里第 90 列（裸鼹鼠 527 / 人 495）的 BEB 后验
用法: python 17_t15_paml_run.py <proj_root> <codeml.exe>
"""
import sys, os, re, json, time, subprocess
from io import StringIO
from Bio import Phylo
from scipy.stats import chi2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT, CODEML = sys.argv[1:3]
WD = os.path.join(ROOT, 'results', 't15_paml')
COL = 527 - 438 + 1
lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

CTL = """seqfile = ../aln.phy
treefile = ../{tree}
outfile = mlc
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
model = {model}
NSsites = {ns}
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = {fo}
omega = {om}
cleandata = 0
fix_blength = {fb}
"""
def run(name, tree, model, ns, fo, om, fb):
    d = os.path.join(WD, name)
    os.makedirs(d, exist_ok=True)
    f = os.path.join(d, 'mlc')
    if os.path.exists(f) and 'Time used' in open(f, errors='ignore').read():
        return None
    open(os.path.join(d, 'codeml.ctl'), 'w').write(CTL.format(tree=tree, model=model, ns=ns, fo=fo, om=om, fb=fb))
    return subprocess.Popen([CODEML, 'codeml.ctl'], cwd=d, stdout=open(os.path.join(d, 'stdout.txt'), 'w'), stderr=subprocess.STDOUT)

t0 = time.time()
p = run('M0', 'tree_plain.nwk', 0, 0, 0, 0.5, 0)
if p:
    p.wait()
P(f'M0 完成（{time.time() - t0:.0f} s）')
m0 = open(os.path.join(WD, 'M0', 'mlc'), errors='ignore').read()
nw = [l for l in m0.splitlines() if l.strip().startswith('(') and re.search(r'\bS\d\d: ', l)][0].strip()   # 带物种代号和枝长的那一行
tr = Phylo.read(StringIO(nw.replace(': ', ':')), 'newick')
allt = frozenset(c.name for c in tr.get_terminals())
BL = {}
for c in tr.find_clades():
    if c is tr.root:
        continue
    s = frozenset(x.name for x in c.get_terminals())
    key = min(s, allt - s, key=lambda z: (len(z), sorted(z)))
    BL[key] = BL.get(key, 0.0) + (c.branch_length or 0.0)

# 把 tree_plain / tree_fg 的拓扑加上 M0 枝长（按「分割」对应），前景标记照搬
def parse_labeled(txt):
    """解析我们自己写的 newick（无枝长，可能带 ' #1'），返回嵌套结构 (children|name, label)"""
    txt = txt.strip().rstrip(';')
    pos = 0
    def node():
        nonlocal pos
        if txt[pos] == '(':
            pos += 1
            ch = [node()]
            while txt[pos] == ',':
                pos += 1
                ch.append(node())
            assert txt[pos] == ')'
            pos += 1
            nm = None
        else:
            m = re.match(r'[A-Za-z0-9_]+', txt[pos:])
            nm = m.group(0); pos += len(nm)
            ch = None
        lab = ''
        m = re.match(r'\s*#1', txt[pos:])
        if m:
            lab = ' #1'; pos += len(m.group(0))
        return (ch, nm, lab)
    return node()
def tips(n):
    return [n[1]] if n[0] is None else [t for c in n[0] for t in tips(c)]
def emit(n, top=False):
    if n[0] is None:
        s = n[1]
    else:
        s = '(' + ','.join(emit(c) for c in n[0]) + ')'
    if top:
        return s + ';'
    st = frozenset(tips(n))
    key = min(st, allt - st, key=lambda z: (len(z), sorted(z)))
    return f'{s}{n[2]}:{BL[key]:.6f}'
for src, dst in (('tree_plain.nwk', 'tree_plain_bl.nwk'), ('tree_fg.nwk', 'tree_fg_bl.nwk')):
    t = open(os.path.join(WD, src)).read().splitlines()
    open(os.path.join(WD, dst), 'w').write(t[0] + '\n' + emit(parse_labeled(t[1]), top=True) + '\n')
P('M0 树长 =', round(sum(BL.values()), 3), '；带枝长的树写入 tree_plain_bl.nwk / tree_fg_bl.nwk')

RUNS = {'M7': ('tree_plain_bl.nwk', 0, 7, 0, 0.5), 'M8': ('tree_plain_bl.nwk', 0, 8, 0, 1.5), 'M8a': ('tree_plain_bl.nwk', 0, 8, 1, 1),
        'A': ('tree_fg_bl.nwk', 2, 2, 0, 1.5), 'Anull': ('tree_fg_bl.nwk', 2, 2, 1, 1)}
t0 = time.time()
procs = [run(k, *v, 2) for k, v in RUNS.items()]
for p in procs:
    if p:
        p.wait()
P(f'5 个模型（固定枝长）并行完成，{time.time() - t0:.0f} s')

def get(k):
    t = open(os.path.join(WD, k, 'mlc'), errors='ignore').read()
    m = re.search(r'lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*(-?[\d.]+)', t)
    rst = open(os.path.join(WD, k, 'rst'), errors='ignore').read() if os.path.exists(os.path.join(WD, k, 'rst')) else ''
    return dict(lnL=float(m.group(2)), np=int(m.group(1)), mlc=t, rst=rst)
res = {k: get(k) for k in ['M0'] + list(RUNS)}
P('\n===== lnL')
for k, v in res.items():
    P(f'  {k:6s} lnL = {v["lnL"]:.3f}  np = {v["np"]}')
w0 = re.search(r'omega \(dN/dS\) =\s*([\d.]+)', res['M0']['mlc'])
P('  M0 ω =', w0.group(1) if w0 else '?')
def lrt(a, b, df, mix=False):
    D = 2 * (res[a]['lnL'] - res[b]['lnL'])
    p = chi2.sf(max(D, 0), df)
    if mix:
        p = 0.5 * p if D > 0 else 1.0
    return D, p
out = {}
for a, b, df, mix, lab in (('M8', 'M7', 2, False, 'M8 vs M7（χ² df=2）'), ('M8', 'M8a', 1, True, 'M8 vs M8a（50:50 混合 χ²）'),
                           ('A', 'Anull', 1, True, '分支-位点 A vs A-null（50:50 混合 χ²）')):
    D, p = lrt(a, b, df, mix)
    out[lab] = (D, p)
    P(f'  {lab}：2ΔlnL = {D:.3f}，p = {p:.3g}')

def section(t, start):
    i = t.find(start)
    return t[i:] if i >= 0 else ''
# 参数
for k in ('M8', 'A', 'Anull'):
    t = res[k]['mlc']
    m = re.search(r'(Parameters in M8.*?\n\n|MLEs of dN/dS \(w\) for site classes.*?\n\n.*?\n\n)', t, re.S)
    m2 = re.search(r'(site class\s+0.*?background w.*?\nforeground w.*?\n)', t, re.S)
    if m2:
        P(f'\n  {k} 位点类别：\n' + m2.group(1))
    elif m:
        P(f'\n  {k} 参数：\n' + m.group(1))
    if k == 'M8':
        m3 = re.search(r'p0 =\s*([\d.]+)\s+p =\s*([\d.]+)\s+q =\s*([\d.]+)\s*\n\s*\(p1 =\s*([\d.]+)\)\s+w =\s*([\d.]+)', t)
        if m3:
            P(f'  M8：p0={m3.group(1)} β(p={m3.group(2)}, q={m3.group(3)})；p1={m3.group(4)} ω={m3.group(5)}')
# 第 90 列的 BEB（rst）
P(f'\n===== 第 {COL} 列（裸鼹鼠 527 / 人 495）的 BEB')
for k in ('M8', 'A'):
    r = section(res[k]['rst'], 'Bayes Empirical Bayes')
    m = re.search(rf'^\s*{COL}\s+(\S)\s+(.*)$', r, re.M)
    P(f'  {k}: ' + (m.group(0).strip() if m else '（rst 中未找到）'))
    sig = re.findall(r'^\s+(\d+)\s+(\S)\s+([\d.]+)(\**)', section(res[k]['mlc'], 'Bayes Empirical Bayes').split('The grid')[0], re.M)
    # 列号 c → 裸鼹鼠 437+c → 人 405+c（该外显子内裸鼹鼠与人无插入缺失，偏移恒为 32）
    P(f'  {k} mlc 列出的 BEB 位点（列/人编号/后验）：' + ('; '.join(f'{a}/人{int(a) + 405}/{c}{d}' for a, b, c, d in sig) if sig else '无'))
json.dump({k: {'lnL': v['lnL'], 'np': v['np']} for k, v in res.items()} | {'LRT': {k: list(v) for k, v in out.items()}},
          open(os.path.join(WD, 'summary.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
open(os.path.join(WD, 'summary.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
