# cGAS "longevity" residues in Hystricognathi: code and data

Code, alignments, selection-test outputs and supplementary tables for the preprint

> *Guinea pigs carry the naked mole-rat cGAS "longevity" residues without the longevity*
> (bioRxiv DOI to be added; archived snapshot: Zenodo DOI to be added)

**What the project asks.** Chen et al. (2025, *Science* 390, eadp5056) showed that four C-terminal
residues of naked mole-rat cGAS (S463, E511, Y527, T530; human D431, K479, L495, K498) let cGAS
promote homologous-recombination repair, and read them as a longevity adaptation from an
18-mammal alignment. We extended the comparison to 262 mammals and 24 hystricognath rodents and
asked three questions: when did the residues appear (ancestral-state reconstruction), are they
maintained (codon-level selection tests against positive controls), and where do they sit in the
protein (crystal structures). We also scanned full-length cGAS for substitutions that really are
specific to the African mole-rat lineage. All analyses are computational and use public data only.

## Licence

- **Code** (`scripts/`): MIT, see `LICENSE`.
- **Data** (everything else: `alignments/`, `selection_tests/`, `tables/`, `supplementary/`,
  `figures/`, `protein_language_model/`): CC BY 4.0, see `LICENSE-DATA.txt`.

## Layout

| Directory / file | Contents |
|---|---|
| `scripts/` | Analysis scripts, numbered in run order. Docstrings and console output are partly in Chinese; the tables below give the English summary |
| `alignments/` | Codon alignment in human coordinates (`aln.fasta`, 275 species x 522 codons), the 117-codon final exon complete for all species (`aln_exon.fasta`), the 258-species real-CDS subset (`aln_cds_only.fasta`), pruned trees (`tree_plain.nwk`, `tree_fg.nwk` with hystricognath branches labelled `{Foreground}`, `tree_cds_only.nwk`), clade definitions, per-sequence provenance (`sources.csv`), and the PAML input set (`paml_input/`) |
| `selection_tests/hyphy/` | HyPhy JSON outputs, logs and batch files: FEL x 5, MEME x 2, RELAX x 2 |
| `selection_tests/paml/` | codeml `mlc` and `codeml.ctl` for every model run, plus the parsed summary |
| `selection_tests/profile_power_t35/` | Shared-omega profile likelihoods, the FEL power simulation and the regulatory-residue reference set |
| `tables/` | Intermediate result tables the manuscript cites (site table, tblastn reads, ancestral states, lineage scan, structure, lattice, bat sites, nucleosome interface) |
| `supplementary/` | Supplementary Tables S1-S10 as CSV with their own `README.md` (row counts and source file of every table), and Supplementary Notes S1-S2 |
| `figures/` | Figures 1, 2 and S1 (PDF with embedded fonts, 300-dpi PNG) and the numbers printed while drawing them |
| `protein_language_model/` | ESM-2 scoring behind Supplementary Note S1 |
| `ACCESSIONS.csv` | Every input accession: RefSeq/GenBank proteins, genome assemblies and WGS projects read by tblastn, UniProt and PDB entries |
| `MANIFEST.csv` | Every file with size and SHA-256 |

Absolute paths of the analysis machine have been removed from logs and HyPhy JSON files (replaced by
project-relative paths). Raw downloads (NCBI packages, BLAST XML, PDB coordinates, the AnAge dump,
the Upham et al. trees) are public and are not redistributed.

## Input data and accessions

| Source | What | Identifier |
|---|---|---|
| NCBI Datasets | cGAS orthologues across Mammalia | human *CGAS* GeneID 115004, taxid 40674 (270 genes, 474 isoforms) |
| NCBI RefSeq / GenBank | Protein and CDS records | per species in `ACCESSIONS.csv` and `supplementary/TableS1.csv` |
| NCBI WGS | Unannotated hystricognath and pteropodid assemblies read by remote tblastn | per species in `ACCESSIONS.csv` (from Tables S1 and S8); *Bathyergus suillus* GCA_988425685.1 |
| UniProt | Human cGAS Q8N884 (site numbering, positive-control annotation); naked mole-rat A0AAX6RS70 (lineage-scan query) | |
| PDB | Human cGAS crystal structures 6CT9, 6CTA, 4O68, 4KM5, 6EDB, 6EDC; cGAS-nucleosome 7CCQ, 6Y5D, 6Y5E | |
| AnAge Build 15 | Maximum lifespan and adult mass | de Magalhaes et al. 2024 |
| VertLife | Upham et al. (2019) mammal trees: completed 5,911-species and DNA-only 4,098-species MCC | |

## Software

`requirements.txt` lists the Python packages (Python 3.12.10, Windows 11, 64-bit). External programs:
HyPhy 2.5.103 (built from source; validated against the HyPhy regression values for FEL, MEME and
RELAX, and against PAML), PAML 4.10.10 `codeml`, NCBI BLAST+ (remote `tblastn`; local BLAST+ 2.17.0
for `32_t30_t2t_check.py`), NCBI Datasets CLI v16, R 4.6 with `ape` (ancestral states), pandoc
(only for the Word export).

Paths to external programs are read from environment variables: `HYPHY` and `HYPHY_PATH`
(`14_t17_hyphy_run.py`), `PANDOC` (`30_t25_export_docx.py`, `31_t25_supplementary.py`). Network
scripts use the standard `http_proxy` / `https_proxy` variables if set. The 32-bit Windows HyPhy
build hangs on exit after writing its JSON, so `14_t17_hyphy_run.py` watches the output and stops the
process once the JSON is complete; with an official Linux/macOS HyPhy build, call it directly with
the arguments recorded in `selection_tests/hyphy/*.log` / `*.run`.

## How to re-run

Scripts expect the project layout they were written in: run them from the project root, with
`results/`, `data/` and `manuscript/` next to `scripts/` (rename `scripts/` to `src/`). Map this
package onto that layout as: `alignments/` -> `results/t17_codon/`, `alignments/paml_input/` ->
`results/t15_paml/`, `selection_tests/hyphy/` -> `results/t17_hyphy/`,
`selection_tests/profile_power_t35/` -> `results/t35/`, `tables/` -> `results/`.

```sh
python -m venv .venv && .venv/bin/pip install -r requirements.txt

# 1. Orthologues, residue calls, life history (network: NCBI)
python src/01_fetch_orthologs.py data/raw/ncbi
python src/02_taxid_align.py
python src/03_protein_titles.py
python src/04_site_table.py
# 2. Dense hystricognath sampling and ancestral states (network: NCBI BLAST, UniProt, Ensembl)
python src/06_t13_fetch_extra.py; python src/07_t13_tblastn.py; python src/08_t13_secondsource.py
python src/09_t13_parse_tblastn.py; python src/10_t13_asr.py; python src/11_t13_495_check.py
# 3. Codon alignment and selection tests
python src/12_t17_codon_aln.py
python src/14_t17_hyphy_run.py ...          # FEL / MEME / RELAX, see selection_tests/hyphy/*.run
python src/15_t17_descriptive.py; python src/16_t17_hyphy_summary.py
python src/15_t15_codon_paml.py; python src/17_t15_paml_run.py     # rodent-wide codeml
python src/23_t23_drop_tblastn.py <root>; python src/28_t23_fel_cds_only.py <root>
python src/33_t35_validate.py <root>; python src/34_t35_power.py <root>
python src/35_t35_regctrl.py <root>; python src/36_t35_profile.py <root>
# 4. Structure
python src/13_t17_structure.py; python src/17_t17_lattice.py
python src/24_t23_sitec.py <root>; python src/29_t23_nucleosome.py <root>
# 5. African mole-rat lineage scan
python src/12_t15_tblastn_full.py; python src/14_t15_lineage_scan.py; python src/16_t15_scan_extra.py
python src/32_t30_t2t_check.py ...          # Bathyergus suillus and T2T checks (local BLAST+)
# 6. Bats, tree placement, figures and tables
python src/22_t23_bat_tblastn.py; python src/26_t23_parse_bats.py; python src/25_t23_dnaonly_check.py <root>
python src/19_t19_fig1_figS1.py; python src/19_t19_fig2.py
python src/t16_manuscript_numbers.py; python src/t18_numbers.py; python src/t18_figS3K_table.py
python src/20_t21_supp_tables.py; python src/27_t23_supp_tables.py; python src/31_t25_supplementary.py
# (packaging: src/67_t90_submission_package.py in the project repository, not included here)
```

Each script's docstring states its exact arguments, inputs and outputs. The network steps are the
slow ones and their outputs are included in `tables/` and `alignments/`, so the selection tests,
structure analysis and figures can be re-run from this package without network access. HyPhy on
four cores took about eight hours in total, dominated by MEME on the full alignment.

Which output supports which claim:

| Output | Analysis | Manuscript |
|---|---|---|
| `selection_tests/hyphy/fel_all.json` | FEL, all branches, 522 codons | Table 2, Figure 2a, Table S5a |
| `fel_hy.json`, `fel_exon.json` | FEL, hystricognath foreground, 522 / 117 codons | E479 within Hystricognathi; Table S5a |
| `fel_exon_all.json`, `fel_cds_only.json` | FEL, all branches, 117 codons / 258 real-CDS species | robustness; Table S5a |
| `meme_all.json`, `meme_hy.json` | MEME | episodic selection; Table S5a |
| `relax.json`, `relax_exon.json` | RELAX | selection intensity; Table S5b |
| `profile_power_t35/power*.json` | FEL power simulation | Results 5, Methods; Table S5c |
| `profile_power_t35/profile.json` | shared-omega profile likelihood | Abstract (1.23, 0.84-1.87; 0.016, 0.003-0.054); Table S5d |
| `selection_tests/paml/` | codeml, 68 rodents | Supplementary Note S2 |

## Verifying the package

```sh
python - <<'PY'
import csv, hashlib, pathlib
bad = [r['path'] for r in csv.DictReader(open('MANIFEST.csv', encoding='utf-8-sig'))
       if r['path'] != 'MANIFEST.csv'
       and hashlib.sha256(pathlib.Path(r['path']).read_bytes()).hexdigest() != r['sha256']]
print('mismatched files:', bad or 'none')
PY
```

## Citation

Cite the preprint and, for this snapshot, its Zenodo DOI (both to be added). Please also cite the
primary sources: NCBI RefSeq, Ensembl, UniProt, the PDB entries above, AnAge Build 15
(de Magalhaes et al. 2024) and Upham et al. (2019).
