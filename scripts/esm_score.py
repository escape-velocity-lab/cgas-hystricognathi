"""E03: 用蛋白语言模型 ESM-2 给各物种 cGAS 在 4 个位点（人编号 D431 / K479 / L495 / K498）上打分。
- 人背景：4 个位点上 20 种氨基酸的掩码打分（遮住该位，取 log p），「人型 → 裸鼹鼠型」单点 LLR、秩，
  以及该位点在哺乳动物里自然出现过的其他残基的 LLR（作对照分布）
- 组合：4 个位点全部 15 种子集（单点 / 两点 / 三点 / 四点），比「单点 LLR 加和」与「伪似然差 ΔPLL」，二者之差 = 模型眼里的上位效应
- 反方向：裸鼹鼠背景上「裸鼹鼠型 → 人型」
- 每个物种：整条序列伪似然 PLL、PLL/L，以及只在比对到人 161~522（催化域）/ 421~508（位点窗口）的位置上的 PLL/L（逐位值存 npz）；在自己序列的 4 个对应位置上遮住打分，S_nmr = Σ[log p(裸鼹鼠型) − log p(人型)]
- 记录 GPU 显存峰值、耗时
用法: python esm_score.py <项目02根目录> [--model facebook/esm2_t33_650M_UR50D] [--batch 32] [--cache <HF cache dir>]
"""
import sys, os, re, json, time, math, itertools, argparse
import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from transformers import AutoTokenizer, EsmForMaskedLM

ap = argparse.ArgumentParser()
ap.add_argument('root')
ap.add_argument('--model', default='facebook/esm2_t33_650M_UR50D')
ap.add_argument('--batch', type=int, default=32)
ap.add_argument('--cache', default=None)
ap.add_argument('--limit', type=int, default=0, help='只跑前 N 个物种（冒烟测试）')
ap.add_argument('--tag', default='v1')
ap.add_argument('--dtype', default='fp16', choices=['fp16', 'bf16', 'fp32'], help='权重精度；softmax 一律在 fp32 上算')
args = ap.parse_args()

ROOT = args.root
PKG = os.path.join(ROOT, 'data', 'raw', 'ncbi', 'pkg', 'ncbi_dataset', 'data')
RES = os.path.join(ROOT, 'results')
SITES = [431, 479, 495, 498]
NMR = {431: 'S', 479: 'E', 495: 'Y', 498: 'T'}
HUM = {431: 'D', 479: 'K', 495: 'L', 498: 'K'}
AA = 'ACDEFGHIKLMNPQRSTVWY'
HYSTRICO = {'Bathyergidae', 'Heterocephalidae', 'Caviidae', 'Chinchillidae', 'Octodontidae', 'Hystricidae',
            'Erethizontidae', 'Dasyproctidae', 'Echimyidae', 'Capromyidae', 'Myocastoridae', 'Thryonomyidae',
            'Petromuridae', 'Ctenomyidae', 'Abrocomidae', 'Cuniculidae', 'Dinomyidae', 'Hydrochoeridae'}

t_start = time.time()
# ---- 序列与比对（与 04_site_table.py 同参数）----
al = PairwiseAligner()
al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score = -10, -0.5
al.end_gap_score = 0
prot = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(PKG, 'protein.faa'), 'fasta')}
H = prot['NP_612450.2']
assert ''.join(H[s - 1] for s in SITES) == 'DKLK'

def h2q_map(q):
    a = al.align(H, q)[0]
    m = {}
    for (hs, he), (qs, qe) in zip(*a.aligned):
        for k in range(he - hs):
            m[hs + k + 1] = qs + k + 1
    return m

st = pd.read_csv(os.path.join(RES, 'site_table_v1.csv'))
sp = st[(st['物种代表'] == '是')].copy()
if args.limit:
    keep = sp[sp.物种.isin(['Homo sapiens', 'Heterocephalus glaber', 'Cavia porcellus', 'Mus musculus'])]
    sp = pd.concat([keep, sp.drop(keep.index).head(max(0, args.limit - len(keep)))])

# ---- 模型 ----
os.environ.setdefault('HF_HUB_DISABLE_XET', '1')
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f'torch {torch.__version__} cuda {torch.version.cuda} available={torch.cuda.is_available()}',
      torch.cuda.get_device_name(0) if dev == 'cuda' else '', flush=True)
t0 = time.time()
tok = AutoTokenizer.from_pretrained(args.model, cache_dir=args.cache)
DT = {'fp16': torch.float16, 'bf16': torch.bfloat16, 'fp32': torch.float32}[args.dtype]
model = EsmForMaskedLM.from_pretrained(args.model, cache_dir=args.cache, attn_implementation='sdpa', dtype=DT).to(dev).eval()
print(f'model {args.model} [{args.dtype}, {model.config._attn_implementation}] loaded in {time.time() - t0:.1f}s', flush=True)
if dev == 'cuda':
    torch.cuda.reset_peak_memory_stats()
MASK = tok.mask_token_id
AA_ID = torch.tensor(tok.convert_tokens_to_ids(list(AA)), device=dev)
n_forward = 0

def encode(seq):
    ids = tok(seq)['input_ids']  # [CLS] seq [EOS]
    assert len(ids) == len(seq) + 2, '分词不是一个残基一个 token'
    return torch.tensor(ids, device=dev)

@torch.inference_mode()
def masked_logp(seq, positions):
    """positions: 1-based 序列位置列表；返回 dict pos -> 20 维 log p（按 AA 顺序），每个位置单独遮住"""
    global n_forward
    ids = encode(seq)
    out = {}
    for i in range(0, len(positions), args.batch):
        ps = positions[i:i + args.batch]
        x = ids.repeat(len(ps), 1)
        for b, p in enumerate(ps):
            x[b, p] = MASK  # token 下标 = 序列位置（因为第 0 个是 CLS）
        lg = model(input_ids=x).logits.float()
        n_forward += len(ps)
        lp = torch.log_softmax(lg, -1)
        for b, p in enumerate(ps):
            out[p] = lp[b, p, AA_ID].cpu().numpy()
    return out

@torch.inference_mode()
def pll(seq):
    """伪似然：每个位置单独遮住，取真实残基的 log p（全词表 softmax）。返回 (总和, 标准残基数, 逐位 log p 数组；非标准残基为 nan)"""
    global n_forward
    ids = encode(seq)
    L = len(seq)
    tot, n = 0.0, 0
    per = np.full(L, np.nan, dtype=np.float32)
    for i in range(1, L + 1, args.batch):
        ps = list(range(i, min(L, i + args.batch - 1) + 1))
        x = ids.repeat(len(ps), 1)
        idx = torch.arange(len(ps), device=dev)
        pos = torch.tensor(ps, device=dev)
        true = ids[pos]
        x[idx, pos] = MASK
        lp = torch.log_softmax(model(input_ids=x).logits.float(), -1)
        n_forward += len(ps)
        v = lp[idx, pos, true]
        std = torch.tensor([seq[p - 1] in AA for p in ps], device=dev)
        tot += v[std].sum().item()
        n += int(std.sum().item())
        vv = v.cpu().numpy(); vv[~std.cpu().numpy()] = np.nan
        per[ps[0] - 1:ps[-1]] = vv
    return tot, n, per

def mutate(seq, muts):
    s = list(seq)
    for p, b in muts:
        s[p - 1] = b
    return ''.join(s)

ai = {a: k for k, a in enumerate(AA)}

# ---- 1. 人背景：4 个位点 × 20 种氨基酸 ----
usable = sp[sp['4位点可用'] == '是']
mh = masked_logp(H, SITES)
rows = []
for s in SITES:
    lp = mh[s]
    order = np.argsort(-lp)
    rank = {AA[k]: r + 1 for r, k in enumerate(order)}
    col = f'{HUM[s]}{s}'
    seen = usable[col].value_counts()
    for a in AA:
        rows.append({'人位点': col, '残基': a, 'log_p': round(float(lp[ai[a]]), 4),
                     'LLR_相对人型': round(float(lp[ai[a]] - lp[ai[HUM[s]]]), 4), '秩': rank[a],
                     '是人型': a == HUM[s], '是裸鼹鼠型': a == NMR[s], '哺乳动物中出现的物种数': int(seen.get(a, 0))})
site20 = pd.DataFrame(rows)
site20.to_csv(os.path.join(RES, f'esm_{args.tag}_human_site_aa.csv'), index=False, encoding='utf-8-sig')
print(site20[site20.是裸鼹鼠型 | site20.是人型].to_string(), flush=True)

# ---- 2. 组合：15 个子集，加和 LLR vs ΔPLL；同时遮住 k 个位点的联合打分 ----
pll_h, nh, _ = pll(H)
combo = []
for k in range(1, 5):
    for sub in itertools.combinations(SITES, k):
        muts = [(s, NMR[s]) for s in sub]
        mseq = mutate(H, muts)
        add = sum(float(mh[s][ai[NMR[s]]] - mh[s][ai[HUM[s]]]) for s in sub)
        pm, _, _ = pll(mseq)
        combo.append({'组合': '+'.join(f'{HUM[s]}{s}{NMR[s]}' for s in sub), '位点数': k,
                      '加和LLR': round(add, 4), 'ΔPLL': round(pm - pll_h, 4), '上位效应_ΔPLL减加和': round(pm - pll_h - add, 4)})
combo = pd.DataFrame(combo)
combo.to_csv(os.path.join(RES, f'esm_{args.tag}_human_combos.csv'), index=False, encoding='utf-8-sig')
print(combo.to_string(), flush=True)

# ---- 3. 反方向：裸鼹鼠背景（NCBI 代表序列 + UniProt A0AAX6RS70）----
rev = []
U = ''.join(str(r.seq) for r in SeqIO.parse(os.path.join(ROOT, 'data', 'raw', 'uniprot', 'A0AAX6RS70.fasta'), 'fasta'))
nmr_acc = st[st.物种 == 'Heterocephalus glaber'].RefSeq蛋白.iloc[0]
for name, q in (('NCBI ' + nmr_acc, prot[nmr_acc]), ('UniProt A0AAX6RS70', U)):
    m = h2q_map(q)
    qp = [m[s] for s in SITES]
    assert ''.join(q[p - 1] for p in qp) == 'SEYT', name
    mq = masked_logp(q, qp)
    p0, _, _ = pll(q)
    for s, p in zip(SITES, qp):
        lp = mq[p]
        rev.append({'背景': name, '人位点': f'{HUM[s]}{s}', '该序列位置': p, '方向': f'{NMR[s]}→{HUM[s]}',
                    'LLR_裸鼹鼠型→人型': round(float(lp[ai[HUM[s]]] - lp[ai[NMR[s]]]), 4),
                    '裸鼹鼠型秩': int((lp > lp[ai[NMR[s]]]).sum() + 1), '人型秩': int((lp > lp[ai[HUM[s]]]).sum() + 1)})
    pr, _, _ = pll(mutate(q, [(p, HUM[s]) for s, p in zip(SITES, qp)]))
    rev.append({'背景': name, '人位点': '四点', '该序列位置': '', '方向': 'SEYT→DKLK',
                'LLR_裸鼹鼠型→人型': round(sum(r['LLR_裸鼹鼠型→人型'] for r in rev[-4:]), 4),
                '裸鼹鼠型秩': '', '人型秩': '', 'ΔPLL': round(pr - p0, 4)})
rev = pd.DataFrame(rev)
rev.to_csv(os.path.join(RES, f'esm_{args.tag}_nmr_reverse.csv'), index=False, encoding='utf-8-sig')
print(rev.to_string(), flush=True)

# ---- 4. 每个物种：PLL + 4 个位点上的掩码打分 ----
REGIONS = {'催化域161_522': (161, 522), '位点窗421_508': (421, 508)}  # 人编号；1~160 是无序的 N 端，长度和重复在物种间差很大
perpos = {}
out = []
t_sp = time.time()
for n, (_, r) in enumerate(sp.iterrows()):
    q = prot[r.RefSeq蛋白]
    rec = {k: r[k] for k in ['物种', '俗名', '目', '科', 'RefSeq蛋白', '蛋白长度', '最大寿命_年', '成年体重_g', 'LQ', '寿命分组',
                             '裸鼹鼠型命中数', '4位点可用', '质量标记']}
    rec['豪猪亚目'] = r.科 in HYSTRICO
    if len(q) > 1022:
        rec['备注'] = '超过 ESM-2 最大长度 1022，跳过'
        out.append(rec); continue
    m = h2q_map(q)
    ok_pos = [s for s in SITES if s in m]
    mq = masked_logp(q, [m[s] for s in ok_pos])
    s_nmr = 0.0
    for s in SITES:
        c = f'{HUM[s]}{s}'
        if s not in m or q[m[s] - 1] not in AA:
            rec[f'{c}_残基'] = '-'; continue
        lp = mq[m[s]]
        own = q[m[s] - 1]
        rec[f'{c}_残基'] = own
        rec[f'{c}_logp自身'] = round(float(lp[ai[own]]), 4)
        rec[f'{c}_nmr减人'] = round(float(lp[ai[NMR[s]]] - lp[ai[HUM[s]]]), 4)
        rec[f'{c}_nmr减自身'] = round(float(lp[ai[NMR[s]]] - lp[ai[own]]), 4)
        rec[f'{c}_自身秩'] = int((lp > lp[ai[own]]).sum() + 1)
        s_nmr += float(lp[ai[NMR[s]]] - lp[ai[HUM[s]]])
    rec['S_nmr'] = round(s_nmr, 4) if len(ok_pos) == 4 else math.nan
    p, nstd, per = pll(q)
    rec['PLL'] = round(p, 3); rec['PLL标准残基数'] = nstd; rec['PLL_per_res'] = round(p / nstd, 5)
    # 整条序列的 PLL/L 受 N 端长度 / 重复片段影响（坑见汇报）→ 只在比对到人特定区段的位置上求平均
    for name, (a, b) in REGIONS.items():
        v = np.array([per[m[h] - 1] for h in range(a, b + 1) if h in m])
        v = v[~np.isnan(v)]
        rec[f'PLL_per_res_{name}'] = round(float(v.mean()), 5) if len(v) else math.nan
        rec[f'{name}_覆盖位点数'] = len(v)
    perpos[r.RefSeq蛋白] = per
    out.append(rec)
    if n % 20 == 0:
        el = time.time() - t_sp
        print(f'[{n + 1}/{len(sp)}] {r.物种} L={len(q)} PLL/L={p / nstd:.4f} S_nmr={s_nmr:.2f} '
              f'elapsed {el:.0f}s peak {torch.cuda.max_memory_allocated() / 2**30 if dev == "cuda" else 0:.2f}GiB', flush=True)
spdf = pd.DataFrame(out)
spdf.to_csv(os.path.join(RES, f'esm_{args.tag}_species.csv'), index=False, encoding='utf-8-sig')
np.savez_compressed(os.path.join(RES, f'esm_{args.tag}_species_perpos.npz'), **perpos)  # 每条序列逐位 log p(真实残基)，以后换区段不用重跑

# ---- 5. 资源记录 ----
res = {'model': args.model, 'dtype': args.dtype, 'batch': args.batch, 'device': torch.cuda.get_device_name(0) if dev == 'cuda' else 'cpu',
       'torch': torch.__version__, 'cuda': torch.version.cuda,
       'peak_mem_allocated_GiB': round(torch.cuda.max_memory_allocated() / 2**30, 3) if dev == 'cuda' else None,
       'peak_mem_reserved_GiB': round(torch.cuda.max_memory_reserved() / 2**30, 3) if dev == 'cuda' else None,
       'n_forward_sequences': n_forward, 'species_scored': int(spdf.PLL.notna().sum()) if 'PLL' in spdf else 0,
       'species_loop_seconds': round(time.time() - t_sp, 1), 'total_seconds': round(time.time() - t_start, 1)}
json.dump(res, open(os.path.join(RES, f'esm_{args.tag}_resources.json'), 'w'), indent=1)
print(json.dumps(res, indent=1))
