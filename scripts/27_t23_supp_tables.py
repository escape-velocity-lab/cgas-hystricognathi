"""27（任务 #23）: 生成本次新增的三张补充表，表头和取值一律英文（投稿要求）。
  S8  狐蝠科 15 套组装的四位点读数（23b）
  S9  275 个物种在 Upham completed / DNA-only 树上的位置来源（23d-c）
  S6c cGAS 三个已报道 DNA 结合位点（site A/B/C）的残基定义和到 DNA 的距离（23d-a），作为 S6 的补表
用法: python 27_t23_supp_tables.py <proj_root>
"""
import sys, os
import pandas as pd

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
OUT = R('manuscript', 'supplementary')
os.makedirs(OUT, exist_ok=True)
SITES = ['431', '479', '495', '498']

# ---------- S8 ----------
b = pd.read_csv(R('results', 't23_bat_sites.csv'), encoding='utf-8-sig')
use = {'tian': 'Tian et al. 2023', 'tian+ctrl': 'Tian et al. 2023; RefSeq-annotated (positive control)',
       'ctrl': 'RefSeq-annotated (positive control)', 'gap': 'named as a gap in an earlier draft', 'new': 'additional pteropodid'}
s8 = pd.DataFrame({
    'species': b['物种'], 'genus': b['物种'].str.split().str[0],
    'assembly': b['组装'], 'wgs_project': b['WGS'],
    'role': b['用途'].map(use), 'status': b['状态'].map({'ok': 'resolved'}).fillna(b['状态']),
    'res_431': b['431'], 'res_479': b['479'], 'res_495': b['495'], 'res_498': b['498'],
    'n_naked_mole_rat_residues': b['裸鼹鼠型残基数'],
    'best_hsp_identity': b['最佳HSP同一性'], 'best_hsp_bitscore': b['bits'], 'best_contig': b['最佳contig'],
    'window_479_plus_minus_15': b['479_窗口'], 'window_431_plus_minus_15': b['431_窗口'],
    'source_431': b['431_来源'], 'source_479': b['479_来源'], 'source_495': b['495_来源'], 'source_498': b['498_来源'],
})
s8.to_csv(os.path.join(OUT, 'supp_table_S8_pteropodids.csv'), index=False, encoding='utf-8-sig')
print('S8:', len(s8), 'rows;  res_479 =', s8['res_479'].value_counts().to_dict())

# ---------- S9 ----------
d = pd.read_csv(R('results', 't23_dnaonly.csv'), encoding='utf-8-sig')
s9 = pd.DataFrame({
    'species': d['物种'], 'order': d['目'], 'family': d['科'],
    'hystricognath': d['豪猪亚目'].map({True: 'yes', False: 'no'}),
    'in_completed_tree_5911sp': d['在completed树'].map({True: 'yes', False: 'no'}),
    'in_DNA_only_tree_4098sp': d['在DNAonly树'].map({True: 'yes', False: 'no'}),
    'placement': d['位置来源'].map({'分子数据': 'molecular data', '分类学插补': 'taxonomic imputation', '不在树上': 'absent from both trees'}),
})
s9.to_csv(os.path.join(OUT, 'supp_table_S9_tree_placement.csv'), index=False, encoding='utf-8-sig')
print('S9:', len(s9), 'rows;  placement =', s9['placement'].value_counts().to_dict())

# ---------- S6c ----------
c = pd.read_csv(R('results', 't23_sitec.csv'), encoding='utf-8-sig')
cat = {'site A': 'site A (Zhou et al. 2018)', 'site B': 'site B (Zhou et al. 2018)',
       'site C alpha(261-286)': 'site C, alpha-region 261-286 (Xie et al. 2019)',
       'site C KRKR(299-302)': 'site C, KRKR-loop 299-302 (Xie et al. 2019)',
       'site C KKH(427-432)': 'site C, KKH-loop 427-432 (Xie et al. 2019)',
       '本文研究的位点': 'position studied here'}
s6c = pd.DataFrame({
    'PDB': c['PDB'], 'space_group': c['空间群'], 'site': c['分类'].map(cat),
    'residue': c['残基'], 'residue_human_number': c['残基号'],
    'min_dist_to_DNA_in_asymmetric_unit_A': c['ASU内到DNA'],
    'min_dist_to_DNA_incl_symmetry_mates_A': c['含对称到DNA'],
}).sort_values(['PDB', 'site', 'residue_human_number'])
s6c['contacts_DNA_within_4A'] = (s6c['min_dist_to_DNA_incl_symmetry_mates_A'] <= 4.0).map({True: 'yes', False: 'no'})
s6c.to_csv(os.path.join(OUT, 'supp_table_S6c_dna_binding_sites.csv'), index=False, encoding='utf-8-sig')
print('S6c:', len(s6c), 'rows')
# ---------- S10 ----------
n = pd.read_csv(R('results', 't23_nucleosome.csv'), encoding='utf-8-sig')
n = n.rename(columns={'cGAS_chain': 'cGAS_copy'})
n['structure_reference'] = n['PDB'].map({'7CCQ': 'Cao et al. 2020', '6Y5E': 'Pathare et al. 2020', '6Y5D': 'Pathare et al. 2020'})
n = n[['PDB', 'structure_reference', 'cGAS_copy', 'residue', 'residue_human_number',
       'min_dist_to_acidic_patch_A', 'nearest_H2A_chain', 'nearest_H2A_residue']]
n = n.sort_values(['PDB', 'cGAS_copy', 'residue_human_number'])
n.to_csv(os.path.join(OUT, 'supp_table_S10_nucleosome_interface.csv'), index=False, encoding='utf-8-sig')
print('S10:', len(n), 'rows;  residues within 5 A per copy =',
      n.groupby(['PDB', 'cGAS_copy']).size().to_dict())

print('  6EDB 里各位点是否接触 DNA（<=4 A）:')
sub = s6c[s6c.PDB == '6EDB']
print(sub[['site', 'residue', 'residue_human_number', 'min_dist_to_DNA_incl_symmetry_mates_A', 'contacts_DNA_within_4A']].to_string(index=False))
