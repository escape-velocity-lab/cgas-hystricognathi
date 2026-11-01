# Supplementary material

Guinea pigs carry the naked mole-rat cGAS "longevity" residues without the longevity
(manuscript v3.1). All tables are UTF-8 CSV. A table that the manuscript cites as one
item but that spans several files is split as a, b, c, d. Paths under "Source" are
relative to the project root; `release/` mirrors them.

| File | Rows | Content | Source |
|---|---|---|---|
| TableS1.csv | 290 | Site table: every cGAS record examined, residues at 431/479/495/498, quality flags, lifespan, LQ | results/site_table_v1.csv + results/t13_species_sites.csv |
| TableS2.csv | 67 | Sequence provenance and verification for hystricognaths | results/t13_secondsource.csv, t13_tblastn_sites.csv |
| TableS3.csv | 88 | Ancestral-state probabilities: 11 nodes x 4 sites x 2 trees, all methods | results/t13_asr_nodes.csv |
| TableS4.csv | 18 | Cross-check against Fig. S3K of Chen et al. (2025): 18 species x 4 positions | results/t18_chen_figS3K.csv |
| TableS5a.csv | 522 | Per-site selection tests, 522 codons: FEL all branches, FEL hystricognath foreground, FEL on the 117-codon final exon (both branch sets), FEL on 258 species with real CDS, MEME (both settings); alpha, beta, omega, p, BH q; positive controls and the 11 validated regulatory residues marked | results/t17_hyphy/{fel_all,fel_hy,fel_exon,fel_exon_all,fel_cds_only,meme_all,meme_hy}.json + results/t17_descriptive.csv; reference-set column from results/t35/regctrl.csv |
| TableS5b.csv | 2 | RELAX on hystricognath branches, full and 117-codon alignments | results/t17_hyphy/relax.json, relax_exon.json |
| TableS5c.csv | 15 | FEL power simulation (the version reported in Results and Methods) | results/t35/power_summary.json, power_summary_alpha0537.json |
| TableS5d.csv | 26 | Shared-omega profile likelihood per residue set (four positions 1.23, catalytic triad 0.016) and per-site intervals | results/t35/profile.json |
| TableS6a.csv | 12 | Structural environment in PDB 6CT9/6CTA: four positions and lineage candidates; distances, SASA, B-factors, pairwise CA distances | results/t17_structure_sites.csv, t17_structure_pairs.csv, manuscript/figures/fig2_numbers.txt, results/t15_lineage_candidates.csv |
| TableS6b.csv | 42 | Crystal-lattice contact check across six crystal forms | results/t17_lattice.csv |
| TableS6c.csv | 82 | Reported DNA-binding sites A/B/C checked in four DNA-bound crystal forms | results/t23_sitec.csv |
| TableS6d.csv | 14 | cGAS-nucleosome acidic-patch interface: every cGAS residue within range of the H2A-H2B acidic patch in 7CCQ, 6Y5E and 6Y5D | results/t23_nucleosome.csv |
| TableS7a.csv | 7 | African mole-rat lineage-specific scan: 7 strict candidates (position 62 excluded by *Bathyergus suillus*), per-species residues, reconstruction probabilities | results/t15_lineage_candidates.csv; position-62 exclusion from results/t30 |
| TableS7b.csv | 98 | Alignment windows (+/-8) at every candidate for every hystricognath, gundi, mouse and human | results/t15_lineage_scan.txt |
| TableS7c.csv | 22 | Relaxed candidate list (4) and naked-mole-rat-terminal list (18) | results/t15_scan_extra.txt |
| TableS8.csv | 15 | Old World fruit bats: the four positions in 15 pteropodid species | results/t23_bat_sites.csv |
| TableS9.csv | 9 | Nine widely cited naked mole-rat longevity mechanisms screened with the same hystricognath panel and the same ancestral-state reconstruction | the table printed in the Supplementary material section of manuscript_v3.1.md |
| TableS10.csv | 275 | Placement of the 275 aligned species on the Upham et al. (2019) trees | results/t23_dnaonly.csv |
| supplementary_notes_S1_S2.md / .docx | - | Supplementary Note S1 (ESM-2 scoring) and Note S2 (PAML at position 495) | manuscript/supplementary/ |

Row counts exclude the header. TableS7a has 7 rows: the scan found seven positions;
position 62 is kept but marked excluded, leaving the six candidates of Table 3.

Two notes on reading the tables together:

- The per-site omega values in TableS5d (1.19, 1.59, 1.41, 0.78) come from the shared-omega
  profile model and differ from the FEL estimates in TableS5a (1.12, 1.50, 1.32, 0.75).
  Table 2 of the manuscript reports the FEL values. The two are different estimators, not a
  disagreement.
- TableS5a q values were recomputed independently from the raw HyPhy JSON with a separate
  Benjamini-Hochberg implementation across all 2,320 cells of the five FEL runs: no mismatch.
