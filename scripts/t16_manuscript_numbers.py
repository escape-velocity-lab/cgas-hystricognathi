"""t16: 从 results/ 现取稿件里用到的全部数字（坑 43：不凭记忆），并生成补充表 S1~S3。
只读 #11 / #13 的产出；不碰 #14（esm*）和 #15（t15_*、figures/）的文件。
用法: python t16_manuscript_numbers.py <项目根目录>
产出: manuscript/numbers.json、manuscript/numbers.txt、manuscript/supp_table_S1.csv / S2 / S3
"""
import sys, os, json
import numpy as np
import pandas as pd

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
OUT = R('manuscript')
os.makedirs(OUT, exist_ok=True)
SITES = ['431', '479', '495', '498']
NMR = dict(zip(SITES, 'SEYT'))
HUM = dict(zip(SITES, 'DKLK'))
AA = set('ACDEFGHIKLMNPQRSTVWY')
# 与 src/13_t15_figures.py 相同的豪猪亚目科名表（Capromys 在本表里归 Echimyidae）
HYST_FAM = {'Bathyergidae', 'Caviidae', 'Chinchillidae', 'Octodontidae', 'Erethizontidae', 'Thryonomyidae', 'Petromuridae',
            'Echimyidae', 'Hystricidae', 'Ctenomyidae', 'Dinomyidae', 'Dasyproctidae', 'Cuniculidae'}
N = {}

# ---------- #11 位点表（RefSeq 直系同源） ----------
v = pd.read_csv(R('results', 'site_table_v1.csv'), encoding='utf-8-sig')
N['v1_genes'] = len(v)
rep = v[v['物种代表'] == '是']
N['v1_species'] = len(rep)
use = rep[rep['4位点可用'] == '是'].copy()
N['v1_usable_species'] = len(use)
N['v1_unusable'] = rep[rep['4位点可用'] != '是']['物种'].tolist()
N['v1_usable_with_LQ'] = int(use['LQ'].notna().sum())
N['v1_usable_LQ_by_source'] = use[use['LQ'].notna()]['寿命来源'].value_counts().to_dict()
N['v1_LQ_groups'] = use['寿命分组'].value_counts().to_dict()
N['v1_orders'] = int(use['目'].nunique())
use['hyst'] = use['科'].isin(HYST_FAM)
N['v1_hyst_species'] = use[use.hyst]['物种'].tolist()
for s in SITES:
    col = [c for c in use.columns if c.endswith(s) and c[0] in 'DKL'][0]
    vc = use[col].value_counts()
    N[f'v1_site{s}_counts'] = vc.to_dict()
    N[f'v1_site{s}_nmr_all'] = int((use[col] == NMR[s]).sum())
    N[f'v1_site{s}_human'] = int((use[col] == HUM[s]).sum())
    out = use[~use.hyst]
    N[f'v1_site{s}_nmr_outside_hyst'] = int((out[col] == NMR[s]).sum())
    N[f'v1_site{s}_nmr_outside_hyst_species'] = out[out[col] == NMR[s]]['物种'].tolist()
N['v1_n_outside_hyst'] = int((~use.hyst).sum())
hits = use['裸鼹鼠型命中数']
N['v1_hits_dist'] = hits.value_counts().sort_index().to_dict()
N['v1_outside_hyst_hits_ge2'] = int(((~use.hyst) & (hits >= 2)).sum())
N['v1_outside_hyst_hits_1'] = int(((~use.hyst) & (hits == 1)).sum())
N['v1_outside_hyst_hits1_species'] = use[(~use.hyst) & (hits == 1)][['物种', '命中位点', 'LQ']].values.tolist()

# AnAge LQ 拟合（与 04_site_table.py 同一口径重算，核对系数）
an = pd.read_csv(R('data', 'raw', 'anage', 'anage_data.txt'), sep='\t', encoding='latin-1')
an = an[an['Class'] == 'Mammalia']
N['anage_species_total'] = None
fit = an[an['Maximum longevity (yrs)'].notna() & an['Adult weight (g)'].notna() & an['Data quality'].isin(['acceptable', 'high'])]
b, a0 = np.polyfit(np.log10(fit['Adult weight (g)'].astype(float)), np.log10(fit['Maximum longevity (yrs)'].astype(float)), 1)
N['anage_mammals'] = len(an)
N['lq_fit_n'] = len(fit)
N['lq_fit_a'] = round(float(a0), 4)
N['lq_fit_b'] = round(float(b), 4)

# ---------- #13 扩展物种表（RefSeq + GenBank + tblastn） ----------
s = pd.read_csv(R('results', 't13_species_sites.csv'), encoding='utf-8-sig')
N['ext_species'] = len(s)
N['ext_source_counts'] = s['来源'].str.split().str[0].value_counts().to_dict()
N['ext_rodents'] = int((s['目'] == 'Rodentia').sum())
s['hyst'] = s['科'].isin(HYST_FAM)
hy = s[s.hyst]
N['ext_hyst_species'] = len(hy)
N['ext_hyst_families'] = int(hy['科'].nunique())
N['ext_hyst_refseq'] = int(hy['来源'].str.startswith('RefSeq').sum())
N['ext_hyst_genbank'] = int(hy['来源'].str.startswith('GenBank').sum())
N['ext_hyst_tblastn'] = int(hy['来源'].str.startswith('tblastn').sum())
for st in SITES:
    rd = hy[hy[st].isin(AA)]
    N[f'hyst_site{st}_readable'] = len(rd)
    N[f'hyst_site{st}_nmr'] = int((rd[st] == NMR[st]).sum())
    N[f'hyst_site{st}_counts'] = rd[st].value_counts().to_dict()
full = hy[hy[SITES].isin(AA).all(axis=1)]
N['hyst_fully_resolved'] = len(full)
N['hyst_unresolved'] = hy[~hy[SITES].isin(AA).all(axis=1)]['物种'].tolist()
N['hyst_hits_dist'] = full['命中'].value_counts().sort_index().to_dict()
cols = ['物种', '科', '最大寿命_年', '成年体重_g', 'LQ']
N['hyst_4of4'] = full[full['命中'] == 4][cols].values.tolist()
N['hyst_3of4'] = full[full['命中'] == 3][cols + ['431', '479', '495', '498']].values.tolist()
ext_out = s[~s.hyst & s[SITES].isin(AA).all(axis=1)]
N['ext_outside_hyst_resolved'] = len(ext_out)
N['ext_outside_hyst_hits_ge2'] = int((ext_out['命中'] >= 2).sum())
gundi = s[s['物种'] == 'Ctenodactylus gundi'].iloc[0]
N['gundi'] = {k: (None if pd.isna(gundi[k]) else gundi[k]) for k in SITES + ['最大寿命_年', 'LQ', '来源']}

# LQ × 命中数（描述性，与图 2 同口径：LQ 有值且 4 位点全可读）
d2 = s[s['LQ'].notna() & s[SITES].isin(AA).all(axis=1)]
N['fig2_n_total'] = len(d2)
N['fig2_n_hyst'] = int(d2.hyst.sum())
N['fig2_n_other'] = int((~d2.hyst).sum())
desc = {}
for k, g in d2.groupby('命中'):
    desc[int(k)] = dict(n=len(g), median=round(float(g.LQ.median()), 2), min=round(float(g.LQ.min()), 2), max=round(float(g.LQ.max()), 2))
N['fig2_LQ_by_hits'] = desc
hd = d2[d2.hyst]
N['hyst_LQ_median'] = round(float(hd.LQ.median()), 2)
N['hyst_LQ_range'] = [round(float(hd.LQ.min()), 2), round(float(hd.LQ.max()), 2)]
N['nonhyst_LQ_median'] = round(float(d2[~d2.hyst].LQ.median()), 2)
# Spearman 仅作描述（物种非独立，不报 p）
N['fig2_spearman_rho_all_descriptive'] = round(float(d2['命中'].rank().corr(d2.LQ.rank())), 2)
N['hyst_spearman_rho_descriptive'] = round(float(hd['命中'].rank().corr(hd.LQ.rank())), 2)

# ---------- #13 祖先重建 ----------
nodes = pd.read_csv(R('results', 't13_asr_nodes.csv'), encoding='utf-8-sig')
N['asr_rows'] = len(nodes)
key_nodes = ['啮齿目祖先 Rodentia', '梳趾鼠+豪猪亚目 Ctenohystrica', '豪猪亚目祖先 Hystricognathi', '非洲鼹鼠科 Bathyergidae s.l.',
             '南美支 Caviomorpha', '豚鼠属 Cavia']
asr = {}
for _, r in nodes.iterrows():
    if r['节点'] in key_nodes:
        asr[f"{r['树']}|{r['位点']}|{r['节点']}"] = [r['Fitch_MPR'], r['ML最可能'], r['ML最可能概率'], r['P裸鼹鼠型_多态ER'],
                                                   r['P裸鼹鼠型_二态ER'], r['P裸鼹鼠型_二态ARD'], int(r['Fitch最少变化'])]
N['asr'] = asr
N['asr_495_check'] = open(R('results', 't13_495_check.txt'), encoding='utf-8').read()

# ---------- 序列核对 ----------
ss = pd.read_csv(R('results', 't13_secondsource.csv'), encoding='utf-8-sig')
N['secondsource_rows'] = len(ss)
tb = pd.read_csv(R('results', 't13_tblastn_sites.csv'), encoding='utf-8-sig')
N['tblastn_jobs'] = len(tb)
N['tblastn_status'] = tb['状态'].value_counts().to_dict()

# ---------- 补充表 ----------
EN_FIX = {'（AnAge 只有婆罗洲猩猩记录，苏门答腊猩猩 2001 年才拆分为独立种）': ' (AnAge lists only P. pygmaeus; P. abelii was split as a separate species in 2001)',
          '（AnAge 只有 S. ehrenbergi 复合种记录，N. galili 是从该复合种拆出的种）': ' (AnAge lists only the S. ehrenbergi complex, from which N. galili was split)'}
def EN_SRC(x):
    for k, w in EN_FIX.items():
        x = x.replace(k, w)
    return EN_SRC0(x)
EN_SRC0 = lambda x: x.replace('taxid 精确', 'AnAge (exact NCBI taxid match)').replace('借用近缘种', 'borrowed from sister taxon').replace('异名', 'synonym')
# S1：位点表（#11 RefSeq 全部 270 条 + #13 新增 20 种）
s1a = v.rename(columns={'物种': 'species', '俗名': 'common_name', '目': 'order', '科': 'family', 'taxid': 'ncbi_taxid',
                        '物种级taxid': 'species_taxid', 'RefSeq蛋白': 'protein_accession', '蛋白长度': 'protein_length_aa',
                        '组装accession': 'assembly', '注释版本': 'annotation_release', 'AnAge学名': 'anage_name',
                        '最大寿命_年': 'max_lifespan_yr', '成年体重_g': 'adult_mass_g', 'AnAge数据质量': 'anage_data_quality',
                        'LQ': 'LQ', 'D431': 'res_431', 'K479': 'res_479', 'L495': 'res_495', 'K498': 'res_498',
                        '裸鼹鼠型命中数': 'n_NMR_residues', '窗口同一性_431/479/495/498': 'window_identity_431/479/495/498',
                        '质量标记': 'quality_flags', '4位点可用': 'usable', '物种代表': 'species_representative', 'CDS核对': 'cds_check'})
s1a['lifespan_source'] = v['寿命来源'].fillna('').map(EN_SRC)
s1a['sequence_source'] = 'RefSeq (NCBI orthologs of human CGAS GeneID 115004)'
s1a['dataset'] = 'RefSeq set (task 11)'
for c in ['usable', 'species_representative']:
    s1a[c] = s1a[c].map(lambda x: 'yes' if x == '是' else ('no' if x == '否' else ('no (subspecies duplicate)' if str(x).startswith('否') else x)))
s1a['cds_check'] = s1a['cds_check'].astype(str).replace({'通过': 'pass', '不通过': 'fail (premature stop in CDS, far from sites)'})
s1a['quality_flags'] = s1a['quality_flags'].fillna('').astype(str).str.replace('CDS翻译对不上', 'CDS translation mismatch').str.replace('落在缺口', ' in alignment gap').str.replace('窗口同一性', ' window identity ').str.replace('长度', 'length ').str.replace('长度异常', 'length outlier')
keep = ['dataset', 'species', 'common_name', 'order', 'family', 'ncbi_taxid', 'species_taxid', 'sequence_source', 'protein_accession',
        'protein_length_aa', 'assembly', 'annotation_release', 'res_431', 'res_479', 'res_495', 'res_498', 'n_NMR_residues',
        'window_identity_431/479/495/498', 'quality_flags', 'usable', 'species_representative', 'cds_check', 'anage_name',
        'lifespan_source', 'max_lifespan_yr', 'adult_mass_g', 'anage_data_quality', 'LQ']
s1a = s1a[keep]
new = s[~s['来源'].str.startswith('RefSeq')].copy()
s1b = pd.DataFrame({'dataset': 'Hystricognathi/outgroup supplement (task 13)', 'species': new['物种'], 'family': new['科'],
                    'order': new['目'], 'sequence_source': new['来源'].str.replace('tblastn', 'tblastn vs WGS'),
                    'res_431': new['431'], 'res_479': new['479'], 'res_495': new['495'], 'res_498': new['498'],
                    'n_NMR_residues': new['命中'], 'usable': np.where(new[SITES].isin(AA).all(axis=1), 'yes', 'partial (only 431 resolved)'),
                    'species_representative': 'yes', 'lifespan_source': new['寿命来源'].fillna('').map(EN_SRC),
                    'max_lifespan_yr': new['最大寿命_年'], 'adult_mass_g': new['成年体重_g'], 'LQ': new['LQ']})
s1 = pd.concat([s1a, s1b], ignore_index=True)[keep]
s1['hystricognath'] = s1['family'].isin(HYST_FAM).map({True: 'yes', False: 'no'})
for c in ['ncbi_taxid', 'species_taxid', 'protein_length_aa']:
    s1[c] = s1[c].astype('Int64')
s1['sequence_source'] = s1['sequence_source'].str.replace('（', ' (').str.replace('）', ')')
s1.to_csv(os.path.join(OUT, 'supp_table_S1.csv'), index=False, encoding='utf-8')
N['S1_rows'] = len(s1)

# S2：豪猪亚目序列来源与核对（多来源注释 + tblastn）
s2a = pd.DataFrame({'check': 'independent annotation / second source', 'species': ss['物种'],
                    'source': ss['来源'].str.replace('（release 当前，compara 直系同源）', ' (current release, Compara ortholog)')
                    .str.replace('（雄性个体组装 HetGla_1.0）', ' (male assembly HetGla_1.0)')
                    .str.replace('（NCBI Eukaryotic Annotation Pipeline）', ' (NCBI Eukaryotic Annotation Pipeline)')
                    .str.replace('GenBank 提交者注释', 'GenBank (submitter annotation)'),
                    'accession': ss['accession'], 'length_aa': ss['长度'], 'res_431': ss['431'], 'res_479': ss['479'],
                    'res_495': ss['495'], 'res_498': ss['498'], 'n_NMR_residues': ss['命中'],
                    'positions_in_sequence': ss['本序列中的位置'].astype(str).str.replace('None', 'NA'),
                    'note': np.where(ss[SITES].eq('-').all(axis=1), 'short isoform; region absent (alternative splicing), not a conflict', '')})
s2b = pd.DataFrame({'check': 'tblastn (human cGAS 400-522) vs genome assembly (WGS)', 'species': tb['物种'],
                    'source': tb['组装'] + ' / ' + tb['WGS'] + ' (' + tb['用途'].astype(str) + ')', 'accession': tb['最佳contig'],
                    'length_aa': np.nan, 'res_431': tb['431'], 'res_479': tb['479'], 'res_495': tb['495'], 'res_498': tb['498'],
                    'n_NMR_residues': tb['命中'], 'positions_in_sequence': tb['431_来源'].astype(str),
                    'note': 'status=' + tb['状态'].astype(str) + '; HSPs=' + tb['HSP数'].astype(str)})
ev = pd.read_csv(R('data', 'raw', 't13', 'refseq_mrna_evidence.tsv'), sep='\t', encoding='utf-8')
s2c = pd.DataFrame({'check': 'RefSeq transcript evidence', 'species': ev['物种'], 'source': 'RefSeq ' + ev['RefSeq状态'],
                    'accession': ev['mRNA'] + ' / ' + ev['蛋白'], 'length_aa': np.nan, 'res_431': '', 'res_479': '', 'res_495': '',
                    'res_498': '', 'n_NMR_residues': np.nan, 'positions_in_sequence': '', 'note': ev['证据'].fillna('(no automated evidence line)')})
s2 = pd.concat([s2a, s2b, s2c], ignore_index=True)
s2['length_aa'] = s2['length_aa'].astype('Int64')
s2['n_NMR_residues'] = s2['n_NMR_residues'].astype('Int64')
s2 = s2.replace({'无命中': 'no hit'}, regex=True)
s2.to_csv(os.path.join(OUT, 'supp_table_S2.csv'), index=False, encoding='utf-8')
N['S2_rows'] = len(s2)

# S3：祖先重建节点概率
NODE_EN = {'啮齿目祖先 Rodentia': 'Rodentia (crown)', '梳趾鼠+豪猪亚目 Ctenohystrica': 'Ctenohystrica (gundi + Hystricognathi)',
           '豪猪亚目祖先 Hystricognathi': 'Hystricognathi (crown)', '非洲鼹鼠+蔗鼠+岩鼠（除豪猪科）': 'Phiomorpha excl. Hystricidae (mole-rats + cane rat + dassie rat)',
           '非洲鼹鼠科 Bathyergidae s.l.': 'Bathyergidae s.l. (naked mole-rat + Fukomys)', '南美支 Caviomorpha': 'Caviomorpha',
           '豚鼠总科 Cavioidea': 'Cavioidea', '豚鼠科 Caviidae': 'Caviidae', '豚鼠属 Cavia': 'Cavia',
           '毛丝鼠总科 Chinchilloidea': 'Chinchilloidea', '八齿鼠总科 Octodontoidea': 'Octodontoidea'}
s3 = pd.DataFrame({'tree': nodes['树'].replace({'Completed5911': 'Upham 2019 completed (5911 spp.) MCC',
                                                'DNAonly4098': 'Upham 2019 DNA-only (4098 spp.) MCC'}),
                   'site': nodes['位点'], 'node': nodes['节点'].map(NODE_EN), 'fitch_MPR_states': nodes['Fitch_MPR'],
                   'ML_multistate_ER_best': nodes['ML最可能'], 'ML_multistate_ER_best_P': nodes['ML最可能概率'],
                   'P_NMR_multistate_ER': nodes['P裸鼹鼠型_多态ER'], 'P_NMR_binary_ER': nodes['P裸鼹鼠型_二态ER'],
                   'P_NMR_binary_ARD': nodes['P裸鼹鼠型_二态ARD'], 'fitch_min_changes_whole_tree': nodes['Fitch最少变化']})
assert s3['node'].notna().all(), '有节点名没翻译'
s3.to_csv(os.path.join(OUT, 'supp_table_S3.csv'), index=False, encoding='utf-8')
N['S3_rows'] = len(s3)
for nm, d in (('S1', s1), ('S2', s2), ('S3', s3)):
    cjk = d.astype(str).apply(lambda c: c.str.contains(r'[一-鿿]')).to_numpy().sum()
    N[f'{nm}_cells_with_chinese'] = int(cjk)


def conv(o):
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return float(o)
    return str(o)


json.dump(N, open(os.path.join(OUT, 'numbers.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=conv)
with open(os.path.join(OUT, 'numbers.txt'), 'w', encoding='utf-8') as f:
    for k, x in N.items():
        if k not in ('asr', 'asr_495_check'):
            f.write(f'{k}: {x}\n')
print(open(os.path.join(OUT, 'numbers.txt'), encoding='utf-8').read())
