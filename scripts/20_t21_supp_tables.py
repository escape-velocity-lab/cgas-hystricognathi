"""任务 #21b：生成缺的补充表 S5 / S6 / S7，把 S4 转成英文，并把 S1~S7 合并成一个 XLSX。

输入：
  results/t17_hyphy/{fel_all,fel_hy,fel_exon,fel_exon_all,meme_all,meme_hy,relax,relax_exon}.json
  results/t17_descriptive.csv（每位点简约变化次数）
  results/t17_power.txt（功效模拟）
  results/t17_structure_sites.csv、t17_structure_pairs.csv、t17_lattice.csv
  manuscript/figures/fig2_numbers.txt（谱系候选的结构距离，与 Figure 2 同源）
  results/t15_lineage_candidates.csv、t15_lineage_scan.txt、t15_scan_extra.txt
  manuscript/supp_table_S1..S4*.csv
输出：manuscript/supplementary/
  supp_table_S5_selection.csv        522 密码子 × 所有选择检验（含 BH q 值）
  supp_table_S5b_relax_power.csv     RELAX 结果 + 功效模拟
  supp_table_S6_structure.csv        4 个位点 + 谱系候选 + 对照残基的结构环境
  supp_table_S6b_lattice.csv         6 套晶体里的晶格接触自查
  supp_table_S7_lineage_scan.csv     非洲鼹鼠谱系特异扫描的候选
  supp_table_S7b_windows.csv         每个候选在每个豪猪亚目物种的 ±8 比对窗口
  supp_table_S4_chen_figS3K_en.csv   S4 的英文版（原文件表头是中文，投稿不能用）
  supplementary_tables_S1_S7.xlsx    合并本，第一个工作表是目录

所有输出表头和取值一律英文。q 值按每个分析自己的有信息位点数做 BH 校正。
用法：python src/20_t21_supp_tables.py
"""
import csv
import json
import os
import re
import sys

import openpyxl
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
OUT = R('manuscript', 'supplementary')
os.makedirs(OUT, exist_ok=True)

CONTROL = {225: 'catalytic triad', 227: 'catalytic triad', 319: 'catalytic triad',
           390: 'zinc thumb', 396: 'zinc thumb', 397: 'zinc thumb', 404: 'zinc thumb',
           173: 'site-A DNA contact', 176: 'site-A DNA contact', 407: 'site-A DNA contact',
           431: 'naked mole-rat position', 479: 'naked mole-rat position',
           495: 'naked mole-rat position', 498: 'naked mole-rat position'}


def read_csv(path):
    with open(path, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh))


def write_csv(path, header, rows):
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f'  {os.path.basename(path)}: {len(rows)} 行')


def bh(pvals):
    """Benjamini-Hochberg，返回与输入同序的 q 值。"""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q, prev = [0.0] * n, 1.0
    for rank, i in enumerate(reversed(order)):
        k = n - rank
        prev = min(prev, pvals[i] * n / k)
        q[i] = prev
    return q


def load_fel(tag, first_pos=1):
    """返回 {人位置: dict(alpha, beta, omega, p, q, verdict)}；q 在有信息位点内做 BH。"""
    path = R('results', 't17_hyphy', f'{tag}.json')
    if not os.path.exists(path):
        print(f'  [跳过] {tag}.json 不存在')
        return {}
    d = json.load(open(path, encoding='utf-8'))
    head = [h[0] for h in d['MLE']['headers']]
    ia, ib, ip = head.index('alpha'), head.index('beta'), head.index('p-value')
    content = d['MLE']['content']['0']
    info = [i for i, r in enumerate(content) if r[ia] > 0 or r[ib] > 0]
    qs = bh([content[i][ip] for i in info])
    out = {}
    for j, i in enumerate(info):
        a, b, p = content[i][ia], content[i][ib], content[i][ip]
        om = b / a if a > 0 else float('inf')
        out[first_pos + i] = dict(
            alpha=round(a, 4), beta=round(b, 4),
            omega=('inf' if a == 0 else round(om, 4)), p=round(p, 4), q=round(qs[j], 4),
            verdict=('purifying' if qs[j] <= 0.05 and b < a else
                     'diversifying' if qs[j] <= 0.05 and b > a else 'n.s.'))
    return out


def load_meme(tag, first_pos=1):
    path = R('results', 't17_hyphy', f'{tag}.json')
    if not os.path.exists(path):
        return {}
    d = json.load(open(path, encoding='utf-8'))
    head = [h[0] for h in d['MLE']['headers']]
    ip = [i for i, h in enumerate(head) if 'p-value' in h.lower()][0]
    content = d['MLE']['content']['0']
    ps = [r[ip] for r in content]
    qs = bh(ps)
    return {first_pos + i: dict(p=round(ps[i], 4), q=round(qs[i], 4)) for i in range(len(ps))}


# ---------------------------------------------------------------- S5
print('S5 选择检验逐位点')
fel_all = load_fel('fel_all')
fel_hy = load_fel('fel_hy')
fel_ex_fg = load_fel('fel_exon', first_pos=406)
fel_ex_all = load_fel('fel_exon_all', first_pos=406)
fel_cds = load_fel('fel_cds_only')   # 任务 #23d：剔除 17 个 tblastn 重建物种后的 258 物种全树 FEL
meme_all = load_meme('meme_all')
meme_hy = load_meme('meme_hy')

desc = {int(r['人位置']): r for r in read_csv(R('results', 't17_descriptive.csv'))}

S5_HEAD = ['human_position', 'human_residue', 'residue_set', 'hystricognath_states',
           'n_changes_hystricognathi_parsimony', 'n_changes_all_mammals_parsimony',
           'n_species_with_data_all_mammals',
           'FELall_alpha', 'FELall_beta', 'FELall_omega', 'FELall_p', 'FELall_q', 'FELall_verdict',
           'FELhyst_alpha', 'FELhyst_beta', 'FELhyst_omega', 'FELhyst_p', 'FELhyst_q', 'FELhyst_verdict',
           'FELexon117all_alpha', 'FELexon117all_beta', 'FELexon117all_omega', 'FELexon117all_p',
           'FELexon117all_q', 'FELexon117all_verdict',
           'FELexon117hyst_alpha', 'FELexon117hyst_beta', 'FELexon117hyst_omega', 'FELexon117hyst_p',
           'FELexon117hyst_q', 'FELexon117hyst_verdict',
           'FELcds258_alpha', 'FELcds258_beta', 'FELcds258_omega', 'FELcds258_p',
           'FELcds258_q', 'FELcds258_verdict',
           'MEMEall_p', 'MEMEall_q', 'MEMEhyst_p', 'MEMEhyst_q']
K = ['alpha', 'beta', 'omega', 'p', 'q', 'verdict']

rows5 = []
for pos in range(1, 523):
    d = desc.get(pos, {})
    row = [pos, d.get('人残基', ''), CONTROL.get(pos, ''), d.get('豪猪亚目状态', ''),
           d.get('豪猪亚目变化次数', ''), d.get('全树变化次数', ''), d.get('全树有数据物种', '')]
    for src in (fel_all, fel_hy, fel_ex_all, fel_ex_fg, fel_cds):
        r = src.get(pos)
        row += [r[k] for k in K] if r else [''] * 6
    for src in (meme_all, meme_hy):
        r = src.get(pos)
        row += [r['p'], r['q']] if r else ['', '']
    rows5.append(row)
write_csv(os.path.join(OUT, 'supp_table_S5_selection.csv'), S5_HEAD, rows5)

# S5b: RELAX + 功效
relax_rows = []
for tag, aln in (('relax', 'full 522-codon alignment'), ('relax_exon', '117-codon final-exon alignment')):
    path = R('results', 't17_hyphy', f'{tag}.json')
    if not os.path.exists(path):
        continue
    d = json.load(open(path, encoding='utf-8'))
    tr = d['test results']
    fits = d['fits']
    relax_rows.append(['RELAX', aln, 'hystricognath branches (test) vs. all others (reference)',
                       round(tr['relaxation or intensification parameter'], 4),
                       round(tr['LRT'], 4), round(tr['p-value'], 4),
                       round(fits['RELAX null']['Log Likelihood'], 2),
                       round(fits['RELAX alternative']['Log Likelihood'], 2),
                       'K > 1 intensified, K < 1 relaxed, K = 1 no difference'])
write_csv(os.path.join(OUT, 'supp_table_S5b_relax_power.csv'),
          ['analysis', 'alignment', 'branch_sets', 'K_or_true_omega', 'LRT', 'p_value',
           'lnL_null', 'lnL_alternative', 'note'], relax_rows)

power_rows = []
for line in open(R('results', 't17_power.txt'), encoding='utf-8'):
    m = re.match(r'\s*([\d.]+)\s*\|\s*(\d+)%\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)–([\d.]+)', line)
    if m:
        power_rows.append(['FEL power simulation', '200 simulated single codons on the 275-species tree '
                                                   '(total length 17.04 subs/site), MG94xREV',
                           m.group(1), f'{m.group(2)}%', m.group(3), f'{m.group(4)}-{m.group(5)}'])
with open(os.path.join(OUT, 'supp_table_S5b_relax_power.csv'), 'a', encoding='utf-8', newline='') as fh:
    w = csv.writer(fh)
    w.writerow([])
    w.writerow(['analysis', 'simulation_setting', 'true_omega', 'fraction_detected_as_purifying_p<=0.1',
                'median_estimated_omega', 'IQR_estimated_omega'])
    w.writerows(power_rows)
print(f'  supp_table_S5b_relax_power.csv: {len(relax_rows)} 行 RELAX + {len(power_rows)} 行功效')

# ---------------------------------------------------------------- S6
print('S6 结构环境')
SEC = {'loop 425–433': 'loop 425-433', 'loop 479–496': 'loop 479-496', 'α螺旋 497–514': 'alpha-helix 497-514'}
sites = read_csv(R('results', 't17_structure_sites.csv'))
S6_HEAD = ['residue_human', 'naked_mole_rat_numbering', 'residue_set', 'PDB_residue_6CT9',
           'secondary_structure', 'dist_to_site_A_DNA_A', 'dist_to_site_A_contact_residues_A',
           'dist_to_catalytic_triad_A', 'dist_to_ATP_6CTA_A', 'dist_to_Zn_A',
           'dist_to_symmetry_DNA_A', 'nearest_symmetry_DNA', 'dist_to_symmetry_cGAS_A',
           'SASA_protein_alone_A2', 'relative_SASA', 'SASA_with_DNA_A2', 'CA_B_factor']
rows6 = []
for r in sites:
    lab = r['位点']                                   # 例 "D431（裸鼹鼠 S463）"
    human = lab.split('（')[0]
    nmr = lab.split('（')[1].rstrip('）').replace('裸鼹鼠 ', '')
    rows6.append([human, nmr, 'naked mole-rat position', r['PDB残基'], SEC.get(r['二级结构'], r['二级结构']),
                  r['到DNA_A位点_Å'], r['到A位点接触残基_Å'], r['到催化残基E225_D227_D319_Å'],
                  r['到ATP_6CTA_Å'], r['到Zn_Å'], r['到对称拷贝DNA_Å'] or 'none within 6 A',
                  r['对称DNA最近残基'], r['到对称拷贝cGAS_Å'] or 'none within 6 A',
                  r['SASA_蛋白单独_Å2'], r['相对SASA'], r['SASA_含DNA_Å2'], r['B因子_CA']])

# 谱系候选的距离：取 Figure 2 存档数字（与图同源）
cand_lab = {}
for line in open(R('manuscript', 'figures', 'fig2_numbers.txt'), encoding='utf-8'):
    m = re.match(r'\s*([A-Z]\d+[A-Z])\s*\(([ABC])\)\s+到 DNA\(A 位点\)\s*([\d.]+) Å\s+到催化三联体\s*([\d.]+) Å', line)
    if m:
        cand_lab[int(m.group(1)[1:-1])] = (m.group(1), m.group(2), m.group(3), m.group(4))
EXCLUDED_62 = ('excluded: the fourth bathyergid Bathyergus suillus (GCA_988425685.1) carries the '
               'ancestral arginine, so the lysine is not shared across Bathyergidae (task #30/#33)')
GRADE_S6 = {201: 'A', 326: 'A', 258: 'B', 388: 'B', 62: 'C', 48: 'C', 19: 'C'}
cand = read_csv(R('results', 't15_lineage_candidates.csv'))
for c in cand:
    pos = int(c['人位置'])
    # 标签 = 豪猪亚目祖先状态 + 人位置 + 非洲鼹鼠科状态（与 Table 1 和 Figure 2 一致）
    lab = f"{c['豪猪亚目祖先ML'][0]}{pos}{c['非洲鼹鼠科状态']}"
    grade = GRADE_S6.get(pos, '')
    if pos in cand_lab:
        _, grade, dna, cat = cand_lab[pos]
        rows6.append([lab, c['裸鼹鼠位置'], f'African mole-rat lineage candidate (confidence {grade})',
                      '', '', dna, '', cat, '', '', '', '', '', '', '', '', ''])
    else:
        note = f'; {EXCLUDED_62}' if pos == 62 else ''
        rows6.append([lab, c['裸鼹鼠位置'],
                      f'African mole-rat lineage candidate (confidence {grade}; outside or disordered '
                      f'in 6CT9 157-522){note}', 'not modelled', '', '', '', '', '', '', '', '', '', '', '', '', ''])
rows6.append(['pairwise CA-CA distances between the four positions (A)', '', '', '', '',
              '; '.join(f"{p['位点对'].replace('–', '-')}: {p['CA_CA_Å']} (nearest heavy atom {p['最近重原子_Å']})"
                        for p in read_csv(R('results', 't17_structure_pairs.csv')))] + [''] * 11)
write_csv(os.path.join(OUT, 'supp_table_S6_structure.csv'), S6_HEAD, rows6)

lat = read_csv(R('results', 't17_lattice.csv'))
write_csv(os.path.join(OUT, 'supp_table_S6b_lattice.csv'),
          ['PDB', 'space_group', 'contains_DNA', 'position', 'residue', 'dist_to_DNA_same_assembly_A',
           'dist_to_symmetry_DNA_A', 'nearest_symmetry_DNA', 'dist_to_symmetry_protein_A', 'residue_set'],
          [[r['PDB'], r['空间群'], 'yes' if r['含DNA'] == 'True' else 'no', r['位点'], r['残基'],
            r['到同单元DNA_Å'] or '>6', r['到对称DNA_Å'] or '>6', r['对称DNA'] or '',
            r['到对称蛋白_Å'] or '>6', CONTROL.get(int(r['位点']), '')] for r in lat])

# ---------------------------------------------------------------- S7
print('S7 谱系特异扫描')
DOMAIN = {'N 端无序区（1–160）': 'disordered N-terminus (1-160)',
          '催化结构域 / NTase 核心（161–330）': 'catalytic domain / nucleotidyltransferase core (161-330)',
          'C 端结构域（331–522，含锌指 390–405）': 'C-terminal domain (331-522, zinc thumb 390-405)'}
GRADE = {201: 'A', 326: 'A', 258: 'B', 388: 'B', 62: 'C', 48: 'C', 19: 'C'}
S7_HEAD = ['human_position', 'naked_mole_rat_position', 'human_residue',
           'ancestral_hystricognath_state', 'African_mole_rat_state', 'confidence',
           'Fukomys_anselli', 'Cavia_porcellus', 'other_hystricognaths_readable',
           'other_hystricognath_residues', 'per_species_detail', 'Ctenodactylus_gundi',
           'non_hystricognath_mammals_with_this_residue', 'window_identity_Damaraland',
           'window_identity_guinea_pig', 'domain', 'UniProt_annotation_human',
           'P_state_African_mole_rat_ancestor', 'P_state_parent_node', 'P_state_hystricognath_ancestor',
           'ML_state_parent_node', 'ML_state_hystricognath_ancestor', 'assigned_to_stem_branch']
rows7 = []
for c in cand:
    pos = int(c['人位置'])
    rows7.append([pos, c['裸鼹鼠位置'], c['人残基'], c['豪猪亚目祖先ML'].split('/')[0][0],
                  c['非洲鼹鼠科状态'], GRADE.get(pos, ''), c['安氏鼹鼠'], c['家豚鼠'],
                  c['其余豪猪亚目可读'], c['其余豪猪亚目残基'], c['其余豪猪亚目明细'], c['梳趾鼠'],
                  c['豪猪亚目外哺乳动物有X'], c['窗口同一性_达马'], c['窗口同一性_豚鼠'],
                  DOMAIN.get(c['结构域'], c['结构域']), c['UniProt注释(人)'],
                  c['P_X_非洲鼹鼠科祖先'], c['P_X_父节点'], c['P_X_豪猪亚目祖先'],
                  c['父节点ML'], c['豪猪亚目祖先ML'],
                  EXCLUDED_62 if pos == 62 else ('yes' if c['主干上出现'] == '是' else 'tentative')])
extra = open(R('results', 't15_scan_extra.txt'), encoding='utf-8').read()
m = re.search(r'放宽档.*?(\[.*?\])', extra, re.S)
if m:
    rows7.append(['relaxed list: positions where exactly one other hystricognath also carries the '
                  'African mole-rat residue (naked mole-rat numbering, residue, species)',
                  m.group(1).replace("'", '')] + [''] * (len(S7_HEAD) - 2))
m = re.search(r'裸鼹鼠独有.*?\n(\[.*?\])', extra, re.S)
if m:
    rows7.append(['naked-mole-rat-terminal list: positions unique to the naked mole-rat among all sampled '
                  'hystricognaths, i.e. not shared with Fukomys (naked mole-rat numbering, residue)',
                  m.group(1).replace("'", '').replace('\n', ' ')] + [''] * (len(S7_HEAD) - 2))
write_csv(os.path.join(OUT, 'supp_table_S7_lineage_scan.csv'), S7_HEAD, rows7)

# S7b：每个候选 ±8 窗口
rows7b, cur = [], None
for line in open(R('results', 't15_lineage_scan.txt'), encoding='utf-8'):
    m = re.match(r'【裸鼹鼠 (\d+) / 人 (\d+)\s+(\S)】', line)
    if m:
        cur = (m.group(2), m.group(1))
        continue
    if cur:
        m = re.match(r'\s{2}(\S[^\s].*?)\s{2,}(\S+)\s+(\S)\s*$', line)
        if m:
            rows7b.append([cur[0], cur[1], m.group(1).strip(), m.group(2), m.group(3)])
        elif not line.strip():
            cur = None
write_csv(os.path.join(OUT, 'supp_table_S7b_windows.csv'),
          ['human_position', 'naked_mole_rat_position', 'species',
           'alignment_window_minus8_to_plus8_naked_mole_rat_coordinates', 'residue_at_position'], rows7b)

# ---------------------------------------------------------------- S4 英文版
print('S4 英文化')
S4_MAP = [('图注物种名', 'figure_label_common_name'), ('拉丁名', 'species'),
          ('Chen_寿命标注_年', 'Chen_lifespan_label_yr'), ('Chen_分类群标注', 'Chen_clade_label'),
          ('Chen_463', 'Chen_NMR463'), ('Chen_511', 'Chen_NMR511'), ('Chen_527', 'Chen_NMR527'),
          ('Chen_530', 'Chen_NMR530'), ('Chen_裸鼹鼠型命中数', 'Chen_n_NMR_residues'),
          ('本项目_431', 'this_study_431'), ('本项目_479', 'this_study_479'),
          ('本项目_495', 'this_study_495'), ('本项目_498', 'this_study_498'),
          ('本项目_裸鼹鼠型命中数', 'this_study_n_NMR_residues'),
          ('本项目_最大寿命_年', 'this_study_max_lifespan_yr_AnAge'), ('本项目_LQ', 'this_study_LQ'),
          ('本项目_序列来源', 'this_study_sequence_source'), ('逐位点一致', 'per_site_agreement'),
          ('结论', 'verdict')]
s4 = read_csv(R('manuscript', 'supp_table_S4_chen_figS3K.csv'))
write_csv(os.path.join(OUT, 'supp_table_S4_chen_figS3K_en.csv'), [e for _, e in S4_MAP],
          [[{'一致': 'agree', '不一致': 'disagree'}.get(r[c], r[c]) for c, _ in S4_MAP] for r in s4])

# ---------------------------------------------------------------- 合并 XLSX
print('合并 XLSX')
SHEETS = [
    ('S1', 'supp_table_S1.csv', R('manuscript', 'supp_table_S1.csv'),
     'Site table. Every mammalian cGAS record examined: accession, assembly, annotation release, the four '
     'residues in human numbering (431/479/495/498), per-site window identity, quality flags, whether the '
     'record was usable, AnAge lifespan source, maximum lifespan, adult mass and longevity quotient (LQ). '
     'Rows cover the 270 RefSeq orthologue records (268 species) plus 20 additional Ctenohystrica species '
     'from GenBank and tblastn.'),
    ('S2', 'supp_table_S2.csv', R('manuscript', 'supp_table_S2.csv'),
     'Sequence provenance and verification for hystricognaths. Independent annotations (RefSeq / Ensembl / '
     'UniProt), independent genome assemblies, direct tblastn reads of genomic DNA, and transcript evidence. '
     'Used to show that the guinea pig S/E/Y/T call does not depend on one gene model.'),
    ('S3', 'supp_table_S3.csv', R('manuscript', 'supp_table_S3.csv'),
     'Ancestral-state reconstruction. Marginal probabilities at 11 key nodes x 4 sites x 2 trees (Upham et al. '
     '2019 completed 5,911-species MCC tree and DNA-only 4,098-species MCC tree) under Fitch parsimony, a '
     'multistate Mk equal-rates model, and binary naked-mole-rat-state models (equal rates and all rates '
     'different).'),
    ('S4', 'supp_table_S4_chen_figS3K_en.csv', os.path.join(OUT, 'supp_table_S4_chen_figS3K_en.csv'),
     'Cross-check against Fig. S3K of Chen et al. (2025). The 18 species of the published alignment, the '
     'residues read from that figure, our independent calls at the corresponding human positions, AnAge '
     'lifespan and LQ, and per-species agreement. All 72 residue calls agree.'),
    ('S5', 'supp_table_S5_selection.csv', os.path.join(OUT, 'supp_table_S5_selection.csv'),
     'Site-wise selection tests for all 522 human-numbered codons. FEL on all branches and with the '
     'hystricognath branches as foreground, each on the full alignment and on the 117-codon final-exon '
     'alignment that is complete for all 275 species; FEL on all branches of the 258 species that have a '
     'real coding sequence, with the 17 tblastn-reconstructed species dropped; MEME in both branch '
     'settings. alpha = synonymous rate, '
     'beta = non-synonymous rate, omega = beta/alpha. q = Benjamini-Hochberg across the informative sites of '
     'that analysis. Positive-control residue sets (catalytic triad, zinc thumb, site-A DNA contacts) and the '
     'four naked mole-rat positions are marked in residue_set. Parsimony change counts are from the same '
     'alignment.'),
    ('S5b', 'supp_table_S5b_relax_power.csv', os.path.join(OUT, 'supp_table_S5b_relax_power.csv'),
     'RELAX test of selection intensity on hystricognath branches (both alignments) and the FEL power '
     'simulation: 200 single codons simulated per true omega on the 275-species tree at the observed total '
     'tree length, analysed with the same likelihood-ratio test.'),
    ('S6', 'supp_table_S6_structure.csv', os.path.join(OUT, 'supp_table_S6_structure.csv'),
     'Structural environment in the human cGAS-DNA structure (PDB 6CT9, 2.26 A; ATP distances from 6CTA). '
     'Minimum distance over C-alpha and side-chain atoms. The four naked mole-rat positions, the structured '
     'African mole-rat lineage candidates, relative solvent accessibility, B-factors, and the pairwise '
     'distances between the four positions. "not modelled" = outside or disordered in the 157-522 construct.'),
    ('S6b', 'supp_table_S6b_lattice.csv', os.path.join(OUT, 'supp_table_S6b_lattice.csv'),
     'Crystal-lattice check. Shortest distance from each position to DNA in the same assembly and to '
     'symmetry-generated DNA, across six human cGAS crystal forms (6CT9, 6CTA, 6EDB, 6EDC, 4O68, 4KM5). '
     'A genuine interface recurs across crystal forms; the contact at 479/498 appears only in P6(1)22 '
     'crystals, whereas the site-A residues K173/R176/K407 contact DNA in every DNA-bound form.'),
    ('S7', 'supp_table_S7_lineage_scan.csv', os.path.join(OUT, 'supp_table_S7_lineage_scan.csv'),
     'African mole-rat lineage-specific scan. Positions where all three African mole-rats share a residue '
     'carried by no other readable hystricognath, with per-species residues, window identities, UniProt '
     'annotation of the human position, and ancestral-state probabilities used to assign the change to the '
     'African mole-rat stem branch. Seven positions qualified in the 24-species scan; position 62 was '
     'subsequently excluded because a fourth bathyergid (Bathyergus suillus, GCA_988425685.1) carries the '
     'ancestral residue, leaving the six candidates of Table 1. The last two rows give the relaxed candidate '
     'list and the naked-mole-rat-terminal list.'),
    ('S7b', 'supp_table_S7b_windows.csv', os.path.join(OUT, 'supp_table_S7b_windows.csv'),
     'Alignment windows (+/-8 residues, naked mole-rat coordinates) at every candidate position for every '
     'hystricognath, the gundi, mouse and human. These are the raw alignments behind Table 1 and allow a '
     'reader to judge whether a one-residue shift could change any call.'),
    # 任务 #23 新增的三张（由 src/27_t23_supp_tables.py 生成）
    ('S6c', 'supp_table_S6c_dna_binding_sites.csv', os.path.join(OUT, 'supp_table_S6c_dna_binding_sites.csv'),
     'The three reported DNA-binding sites of human cGAS, checked against four DNA-bound crystal forms. '
     'Site A and site B as defined by Zhou et al. (2018); site C (alpha-region 261-286, KRKR-loop 299-302, '
     'KKH-loop 427-432) as defined by Xie et al. (2019). Distances are minimum side-chain-to-DNA, reported '
     'both within the asymmetric unit and over crystal symmetry mates, because site B and site C are formed '
     'by symmetry-related DNA in these crystals. Position 431 lies inside the KKH-loop but, unlike the three '
     'basic residues of that loop, never approaches DNA within 6 A.'),
    ('S8', 'supp_table_S8_pteropodids.csv', os.path.join(OUT, 'supp_table_S8_pteropodids.csv'),
     'The four positions in every Old World fruit bat (Pteropodidae) genome assembly available from NCBI: '
     '15 species in seven genera, including all six pteropodids analysed by Tian et al. (2023). Four species '
     'have an annotated cGAS and serve as positive controls for the tblastn procedure; the other eleven were '
     'read by remote tblastn of human cGAS residues 400-522 against the single WGS project. No assembly '
     'carries the naked mole-rat glutamate at human position 479.'),
    ('S9', 'supp_table_S9_tree_placement.csv', os.path.join(OUT, 'supp_table_S9_tree_placement.csv'),
     'Placement of all 275 aligned species on the Upham et al. (2019) trees. Species present in the '
     'DNA-only 4,098-species tree are placed by molecular data; species present only in the completed '
     '5,911-species tree are placed by taxonomic imputation. Two species are imputed: Fukomys anselli and '
     'Leopardus geoffroyi. Cavia tschudii is placed by molecular data.'),
    ('S10', 'supp_table_S10_nucleosome_interface.csv', os.path.join(OUT, 'supp_table_S10_nucleosome_interface.csv'),
     'The cGAS-nucleosome acidic-patch interface. Every cGAS residue within 5 A of the histone H2A acidic '
     'patch (H2A 56, 61, 64, 90, 91, 92), listed separately for each of the four independent cGAS copies in '
     'three cryo-EM structures: 7CCQ (Cao et al. 2020) and 6Y5E and 6Y5D (Pathare et al. 2020); 6Y5D '
     'contains two cGAS copies. R255, the arginine anchor, and P257 are the only residues present in all '
     'four copies. K258, the position substituted to glutamine in African mole-rats, is present in three of '
     'the four at 4.0-4.2 A.'),
]

wb = openpyxl.Workbook()
toc = wb.active
toc.title = 'Contents'
toc.append(['Supplementary Tables S1-S10'])
toc.append(['cGAS residues in Hystricognathi - Brief Report (draft_v3)'])
toc.append([])
toc.append(['Sheet', 'Table', 'Rows', 'Description'])
toc['A1'].font = Font(bold=True, size=14)
for c in toc[4]:
    c.font = Font(bold=True)

for sheet, fname, path, descr in SHEETS:
    rows = list(csv.reader(open(path, encoding='utf-8-sig', newline='')))
    ws = wb.create_sheet(sheet)
    for row in rows:
        ws.append(row)
    ws.freeze_panes = 'A2'
    for col in range(1, min(ws.max_column, 40) + 1):
        width = max((len(str(ws.cell(r, col).value or '')) for r in range(1, min(ws.max_row, 200) + 1)),
                    default=10)
        ws.column_dimensions[get_column_letter(col)].width = min(max(width + 2, 10), 55)
    for c in ws[1]:
        c.font = Font(bold=True)
    toc.append([sheet, f'Table {sheet}', max(len(rows) - 1, 0), descr])

toc.column_dimensions['A'].width = 8
toc.column_dimensions['B'].width = 12
toc.column_dimensions['C'].width = 8
toc.column_dimensions['D'].width = 120
for row in toc.iter_rows(min_row=4):
    row[3].alignment = Alignment(wrap_text=True, vertical='top')
toc.append([])
toc.append(['Generated by src/20_t21_supp_tables.py; every number traces to results/ in the project '
            'repository (see README.md in the release package).'])

xlsx = os.path.join(OUT, 'supplementary_tables_S1_S10.xlsx')
wb.save(xlsx)
print(f'  {os.path.basename(xlsx)}: {len(SHEETS)} 个数据工作表 + 目录')
