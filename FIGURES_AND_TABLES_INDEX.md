# Figures & tables — source index

Section-by-section index of every figure and table, its file, and the script that
generates it. Figure numbers match the manuscript float order. All generating scripts
live in `scripts/`; figures in `figures/` (main) and `supp_figures/` (supplementary);
editable table exports (CSV + `illex_tables.xlsx`) in `tables_export/`.

## Main figures

| Fig | File (`figures/`) | Script (`scripts/`) | Subject |
|-----|-------------------|---------------------|---------|
| 1 | `fig1_map.png` | `plot_map.py` | Sampling map + chromosome-2 arrangement frequency |
| 2 | `fig2_overview.png` | `plot_fig1_overview.py` | No structure outside chr2; ZZ/ZO sex assignment |
| 3 | `fig3_inversion.png` | `fig2_inversion.py` | The chromosome-2 inversion (PCA, FST, per-arrangement diversity) |
| 4 | `fig4_demography.png` | `fig3_demography_inversion.py` | Demographic history + inversion frequency history |
| 5 | `fig5_selection.png` | `fig7_selection_manhattan.py` | Genome-wide selective-sweep scan (autosomes + chrZ) |
| 6 | `fig6_diversity.png` | `plot_diversity.py` | Uniform diversity; no geographic structure (FST≈0) |

## Supplementary figures

Supplementary figures are numbered by first citation in main.tex (renumbered 2026-11-13). Image filenames keep their earlier numbers.

| Fig | File (`supp_figures/`) | Script (`scripts/`) | Subject |
|-----|------------------------|---------------------|---------|
| S1 | `figS3_winpca.png` | (upstream) | Windowed PCA across chr2 |
| S2 | `figS1_pca_geography.png` | `plot_pca_geography.py` | PCA spread is not geographic structure |
| S3 | `figS2_fst_expanded.png` | `plot_fst_persnp.py` | Per-SNP FST vs neutral null; no local adaptation |
| S4 | `figS19_fst_genomewide.png` | (upstream) | Genome-wide per-SNP FST finds no outlier |
| S5 | `figS6_cline.png` | (upstream) | The inversion is non-clinal |
| S6 | `figS20_chr2_karyotype_fst.png` | `plot_chr2_karyotype_fst.py` | Between-karyotype (AA vs BB) FST along chr2 |
| S7 | `figS11_breakpoints.png` | (upstream) | Inversion breakpoints |
| S8 | `figS4_ld_triangle.png` | (upstream) | LD triangle across chr2 |
| S9 | `figS26_abc_confusion.png` | see supplement.tex (sfig:abcconf) | Confusion matrix of the ABC-RF model choice for the chromosome-2 inversion. |
| S10 | `figS14b_npe_posterior.png`, `figS14c_abc_trajectory.png` | see supplement.tex (sfig:profiles) | Posterior of the balancing model with gene flux for the chromosome-2 inversion |
| S11 | `figS25_abc_ppc.png` | see supplement.tex (sfig:abcppc) | Posterior predictive check of the balancing model with gene flux. |
| S12 | `figS5_joint_sfs.png` | (upstream) | Joint 2D SFS, AA vs BB arrangement |
| S13 | `figS21_inversion_go_qq.png` | `plot_inversion_go_qq.py` (+ `inversion_go_permutation_null.py`) | Inversion-gene GO enrichment is null (permutation QQ) |
| S14 | `figS12_recomb_boxplots.png` | scripts/plot_recomb_boxplots.py | Per-chromosome male recombination rate (200-rep bootstrap-corrected), with r_Z/r_A |
| S15 | `figS17_recomb_selection.png` | `plot_recomb_selection.py` | Recombination, gene proximity and selection efficacy |
| S16 | `figS22_chr1_recomb_landscape.png` | `plot_chr1_recomb_landscape.py` | chr1 male ReLERNN landscape (200-rep corrected) + CIs; masked block and centromere gap; compressed dynamic range |
| S17 | `figS10_mkado.png` | `plot_supp_mkado.py` | Per-gene McDonald–Kreitman neutrality index (volcano) |
| S18 | `figS23_pi_blocks_by_division.png` | see supplement.tex (sfig:piblocks) | The 22 elevated-diversity blocks are shared by every sampling division. |
| S19 | `figS24_block_ld_triangles.png` | see supplement.tex (sfig:blockld) | The elevated-diversity blocks carry no inversion-like linkage disequilibrium. |
| S20 | `figS7_confusion.png` | (upstream) | Sweep-classifier confusion matrix |
| S21 | `figS8_pipeline.png` | `plot_pipeline_schematic.py` | The selective-sweep scan pipeline |
| S22 | `figS13_markov.png` | `fig_markov.py` | Recombination-aware Markov decoder for Tier-1 calls |
| S23 | `figS18_tier1_manhattan.png` | `plot_tier1_manhattan.py` | Genome-wide Tier-1 sweep calls (autosomes + chrZ) |
| S24 | `figS15_gsea.png` | `plot_gsea.py` | GO enrichment of the sweep genes (Tier-1 vs BGS-robust) |
| S25 | `figS16_go_network.png` | `plot_go_network.py` | GO enrichment map (BGS-robust sweep genes) |
| S26 | `figS9_diversity_boxplots.png` | `plot_diversity.py` | Per-chromosome diversity & Tajima's D |

## Tables

Editable exports in `tables_export/` — one CSV per table plus a multi-sheet
`illex_tables.xlsx`. Regenerate with `scripts/export_tables.py` after `make all`.

Supplementary tables are numbered by first citation in main.tex (renumbered 2026-11-13). Exported CSV and sheet names follow these numbers (from main.aux and supplement.aux).

| Table | Export | Source (`.tex`) | Subject |
|-------|--------|-----------------|---------|
| 1 | `Table1_summary.csv` | `main.tex` (tab:summary) | Sample & whole-genome summary |
| 2 | `Table2_inv.csv` | `main.tex` (tab:inv) | The chromosome-2 inversion |
| 3 | `Table3_demog.csv` | `main.tex` (tab:demog) | Demographic model |
| 4 | `Table4_sweep.csv` | `main.tex` (tab:sweep) | High-confidence sweep candidates (incl. chrZ) |
| S1 | `TableS1_sampling.csv` | `sampling_table.tex` (stab:sampling; longtable) | Sampling by collection event and division |
| S2 | `TableS2_geo.csv` | `supplement.tex` (stab:geo) | chr2 karyotype counts by NAFO division |
| S3 | `TableS3_coindetii.csv` | `supplement.tex` (stab:coindetii) | Outgroup usage |
| S4 | `TableS4_invgenes.csv` | `supplement.tex` (stab:invgenes) | The 45 annotated genes captured by the chr2 inversion |
| S5 | `TableS5_piblocks.csv` | `pi_blocks_table.tex` (stab:piblocks) | Anomalous-diversity blocks |
| S6 | `TableS6_goenrich.csv` | `supplement.tex` (stab:goenrich) | Tier-1 GO enrichment |
| S7 | `TableS7_goclusters.csv` | `supplement.tex` (stab:goclusters) | Functional clusters (BGS-robust) |
| S8 | `TableS8_candidates.csv` | `supplement.tex` (stab:candidates) | Tier-2 cross-method-concordant genes |
| S9 | `TableS9_zsex.csv` | `supplement.tex` (stab:zsex) | Z recombination and sex-biased mutation |
| S10 | `TableS10_software.csv` | `software_table.tex` | Software and versions |
