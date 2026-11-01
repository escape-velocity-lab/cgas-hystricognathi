"""任务 #21c：把要公开的脚本和数据整理到 release/cgas-hystricognathi/。

按 manuscript/submission_checklist.md 第 2 节的「公开 / 不公开」清单复制文件，
生成 requirements.txt、LICENSE 文件和带 sha256 的 MANIFEST.csv。
README.md 和 UPLOAD_STEPS.md 是手写的，本脚本不覆盖它们。

**不上传任何东西**，只在本地准备目录。
用法：python src/21_t21_release_package.py
"""
import csv
import glob
import hashlib
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
REL = R('release', 'cgas-hystricognathi')

# (目标子目录, 源 glob) —— 源路径相对项目根
COPY = [
    ('scripts', 'src/*.py'),
    ('alignments', 'results/t17_codon/aln.fasta'),
    ('alignments', 'results/t17_codon/aln_exon.fasta'),
    ('alignments', 'results/t17_codon/tree_plain.nwk'),
    ('alignments', 'results/t17_codon/tree_fg.nwk'),
    ('alignments', 'results/t17_codon/clades.json'),
    ('alignments', 'results/t17_codon/sources.csv'),
    ('alignments', 'results/t17_codon/log.txt'),
    ('alignments/paml_input', 'results/t15_paml/aln.phy'),
    ('alignments/paml_input', 'results/t15_paml/tree_*.nwk'),
    ('alignments/paml_input', 'results/t15_paml/codon_sources.csv'),
    ('alignments/paml_input', 'results/t15_paml/codes.json'),
    ('alignments/paml_input', 'results/t15_paml/prep_log.txt'),
    ('tables', 'results/site_table_v1.csv'),
    ('tables', 'results/site_windows_v1.txt'),
    ('tables', 'results/t13_species_sites.csv'),
    ('tables', 'results/t13_secondsource.csv'),
    ('tables', 'results/t13_tblastn_sites.csv'),
    ('tables', 'results/t13_asr_nodes.csv'),
    ('tables', 'results/t13_asr.txt'),
    ('tables', 'results/t13_495_check.txt'),
    ('tables', 'results/t15_aln_nmrcoord.tsv'),
    ('tables', 'results/t15_lineage_candidates.csv'),
    ('tables', 'results/t15_lineage_scan.txt'),
    ('tables', 'results/t15_scan_extra.txt'),
    ('tables', 'results/t15_tblastn_full_coverage.csv'),
    ('tables', 'results/t17_descriptive.csv'),
    ('tables', 'results/t17_descriptive.txt'),
    ('tables', 'results/t17_hyphy_sites.csv'),
    ('tables', 'results/t17_hyphy_summary.txt'),
    ('tables', 'results/t17_power.txt'),
    ('tables', 'results/t17_structure.txt'),
    ('tables', 'results/t17_structure_sites.csv'),
    ('tables', 'results/t17_structure_pairs.csv'),
    ('tables', 'results/t17_lattice.csv'),
    ('tables', 'results/t17_lattice.txt'),
    ('tables', 'results/t18_chen_figS3K.csv'),
    ('tables', 'results/t18_numbers.json'),
    ('tables', 'data/processed/orthologs_x_anage.csv'),
    ('tables', 'data/processed/taxid_manual_check.csv'),
    ('tables', 'manuscript/numbers.json'),
    ('selection_tests/hyphy', 'results/t17_hyphy/*.json'),
    ('selection_tests/hyphy', 'results/t17_hyphy/*.log'),
    ('selection_tests/paml', 'results/t15_paml/summary.json'),
    ('selection_tests/paml', 'results/t15_paml/summary.txt'),
    ('protein_language_model', 'results/esm_v1.md'),
    ('protein_language_model', 'results/esm_v1_prediction.md'),
    ('protein_language_model', 'results/esm_v1_summary.txt'),
    ('protein_language_model', 'results/esm_v1_*.csv'),
    ('protein_language_model', 'results/esm_v1_resources.json'),
    ('supplementary', 'manuscript/supplementary/*.csv'),
    ('supplementary', 'manuscript/supplementary/*.xlsx'),
    ('supplementary', 'manuscript/supp_table_S1.csv'),
    ('supplementary', 'manuscript/supp_table_S2.csv'),
    ('supplementary', 'manuscript/supp_table_S3.csv'),
    ('figures', 'manuscript/figures/Figure*.pdf'),
    ('figures', 'manuscript/figures/Figure*.png'),
    ('figures', 'manuscript/figures/Figure1_gain_branches.csv'),
    ('figures', 'manuscript/figures/fig2_numbers.txt'),
]
# PAML 每个模型的输出目录：只收 mlc（结果）和 codeml.ctl（设置），不收 rst/rub 等中间文件
PAML_MODELS = ['M0', 'M7', 'M8', 'M8a', 'A', 'Anull']

MIT = """MIT License

Copyright (c) 2026 [AUTHOR NAME TBD]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

DATA_LICENSE = """The derived data files in this repository (tables/, alignments/, selection_tests/,
supplementary/, figures/) are released under the Creative Commons Attribution 4.0
International licence (CC BY 4.0): https://creativecommons.org/licenses/by/4.0/

They are derived from public databases that carry their own terms:
  NCBI RefSeq / GenBank / WGS   public domain (see NCBI data usage policies)
  Ensembl                       no restriction on use
  UniProt                       CC BY 4.0
  Protein Data Bank (6CT9 etc.) CC0 1.0
  AnAge Build 15 (HAGR)         free for academic use; cite de Magalhaes et al. 2024
  Upham et al. (2019) trees     CC BY 4.0, via VertLife

No sequence, figure or supplementary file from Chen et al. (2025) is redistributed
here. Table S4 contains only the residue identities we read from the published
figure, which are facts, together with our own independent calls.
"""


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


if os.path.exists(REL):
    # 只清空由本脚本生成的部分，保留手写的 README.md / UPLOAD_STEPS.md
    for name in os.listdir(REL):
        p = os.path.join(REL, name)
        if os.path.isdir(p):
            shutil.rmtree(p)
        elif name not in ('README.md', 'UPLOAD_STEPS.md'):
            os.remove(p)
os.makedirs(REL, exist_ok=True)

n = 0
for sub, pattern in COPY:
    dest = os.path.join(REL, sub)
    os.makedirs(dest, exist_ok=True)
    hits = sorted(glob.glob(R(*pattern.split('/'))))
    if not hits:
        print(f'  [警告] 没有匹配：{pattern}')
    for src in hits:
        if '__pycache__' in src:
            continue
        shutil.copy2(src, os.path.join(dest, os.path.basename(src)))
        n += 1
for m in PAML_MODELS:
    dest = os.path.join(REL, 'selection_tests', 'paml', m)
    os.makedirs(dest, exist_ok=True)
    for f in ('mlc', 'codeml.ctl'):
        src = R('results', 't15_paml', m, f)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(dest, f))
            n += 1
print(f'复制 {n} 个文件')

open(os.path.join(REL, 'LICENSE'), 'w', encoding='utf-8').write(MIT)
open(os.path.join(REL, 'LICENSE-DATA.txt'), 'w', encoding='utf-8').write(DATA_LICENSE)

# requirements：只列项目实际 import 的包，版本取当前环境
USED = ['biopython', 'gemmi', 'matplotlib', 'numpy', 'openpyxl', 'pandas', 'scipy', 'tabulate']
freeze = dict(line.split('==') for line in
              subprocess.run([sys.executable, '-m', 'pip', 'freeze'], capture_output=True,
                             text=True).stdout.splitlines() if '==' in line)
with open(os.path.join(REL, 'requirements.txt'), 'w', encoding='utf-8') as fh:
    fh.write('# Python 3.12.10 on Windows 11 (64-bit). Versions are those used for the analyses.\n')
    for pkg in USED:
        key = next((k for k in freeze if k.lower() == pkg.lower()), None)
        fh.write(f'{key}=={freeze[key]}\n' if key else f'{pkg}\n')
    fh.write('# Optional, only for the protein language model (Supplementary Note S1):\n')
    fh.write('# torch==2.9.0+cu128, fair-esm==2.0.0 (ESM-2 650M, esm2_t33_650M_UR50D)\n')
    fh.write('# External programs: HyPhy 2.5.103 (built from source, see README), PAML 4.10.10,\n')
    fh.write('# NCBI BLAST+ 2.16.0 remote tblastn, NCBI Datasets CLI v16.\n')

rows = []
for dirpath, _, files in os.walk(REL):
    for f in sorted(files):
        p = os.path.join(dirpath, f)
        rel = os.path.relpath(p, REL).replace('\\', '/')
        if rel == 'MANIFEST.csv':
            continue
        rows.append([rel, os.path.getsize(p), sha256(p)])
rows.sort()
with open(os.path.join(REL, 'MANIFEST.csv'), 'w', encoding='utf-8', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['path', 'size_bytes', 'sha256'])
    w.writerows(rows)

total = sum(r[1] for r in rows)
print(f'MANIFEST.csv：{len(rows)} 个文件，合计 {total / 1e6:.1f} MB')
print(f'目录：{REL}')
