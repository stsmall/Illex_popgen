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

| Fig | File (`supp_figures/`) | Script (`scripts/`) | Subject |
|-----|------------------------|---------------------|---------|
| S1 | `figS1_pca_geography.png` | `plot_pca_geography.py` | PCA spread is not geographic structure |
| S2 | `figS2_fst_expanded.png` | `plot_fst_persnp.py` | Per-SNP FST vs neutral null; no local adaptation |
| S3 | `figS3_winpca.png` | (upstream) | Windowed PCA across chr2 |
| S4 | `figS4_ld_triangle.png` | (upstream) | LD triangle across chr2 |
| S5 | `figS5_joint_sfs.png` | (upstream) | Joint 2D SFS, AA vs BB arrangement |
| S6 | `figS6_cline.png` | (upstream) | The inversion is non-clinal |
| S7 | `figS7_confusion.png` | (upstream) | Sweep-classifier confusion matrix |
| S8 | `figS8_pipeline.png` | `plot_pipeline_schematic.py` | The selective-sweep scan pipeline |
| S9 | `figS9_diversity_boxplots.png` | `plot_diversity.py` | Per-chromosome diversity & Tajima's D |
| S10 | `figS10_mkado.png` | `plot_supp_mkado.py` | Per-gene McDonald–Kreitman neutrality index (volcano) |
| S11 | `figS11_breakpoints.png` | (upstream) | Inversion breakpoints |
| S12 | `figS12_recomb_boxplots.png` | (upstream) | Per-chromosome recombination-rate distribution |
| S13 | `figS13_markov.png` | `fig_markov.py` | Recombination-aware Markov decoder for Tier-1 calls |
| S14 | `figS14_inversion_profiles.png` | `plot_inversion_profiles.py` | Profile-likelihood bounds on the inversion decline |
| S15 | `figS15_gsea.png` | `plot_gsea.py` | GO enrichment of the sweep genes (Tier-1 vs BGS-robust) |
| S16 | `figS16_go_network.png` | `plot_go_network.py` | GO enrichment map (BGS-robust sweep genes) |
| S17 | `figS17_recomb_selection.png` | `plot_recomb_selection.py` | Recombination, gene proximity and selection efficacy |
| S18 | `figS18_tier1_manhattan.png` | `plot_tier1_manhattan.py` | Genome-wide Tier-1 sweep calls (autosomes + chrZ) |
| S19 | `figS19_fst_genomewide.png` | (upstream) | Genome-wide per-SNP FST finds no outlier |

## Tables

Editable exports in `tables_export/` — one CSV per table plus a multi-sheet
`illex_tables.xlsx`. Regenerate with `scripts/export_tables.py`.

| Table | Export | Source (`.tex`) | Subject |
|-------|--------|-----------------|---------|
| 1 | `Table1_summary.csv` | `main.tex` (tab:summary) | Sample & whole-genome summary |
| 2 | `Table2_inv.csv` | `main.tex` (tab:inv) | The chromosome-2 inversion |
| 3 | `Table3_demog.csv` | `main.tex` (tab:demog) | Demographic model |
| 4 | `Table4_sweep.csv` | `main.tex` (tab:sweep) | High-confidence sweep candidates (incl. chrZ) |
| S1 | `TableS1_software.csv` | `software_table.tex` | Software and versions |
| S2 | `TableS2_geo.csv` | `supplement.tex` (stab:geo) | chr2 karyotype counts by NAFO division |
| S3 | `TableS3_coindetii.csv` | `supplement.tex` (stab:coindetii) | Outgroup usage |
| S4 | `TableS4_goenrich.csv` | `supplement.tex` (stab:goenrich) | Tier-1 GO enrichment |
| S5 | `TableS5_goclusters.csv` | `supplement.tex` (stab:goclusters) | Functional clusters (BGS-robust) |
| S6 | `TableS6_candidates.csv` | `supplement.tex` (stab:candidates) | Tier-2 cross-method-concordant genes |
