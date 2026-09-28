# Analysis pipeline

Scripts for every analysis in the manuscript, copied from the working analysis tree (`analysis/steps/`).
Paths inside the scripts refer to that tree on the analysis server; `config.sh` defines them.

| Directory | Contents |
|---|---|
| `callset/grenepipe_config.yaml` | Read trimming (fastp), mapping (BWA-MEM2), duplicate removal (Picard) and GATK HaplotypeCaller settings |
| `steps/00_callset` | Merging, site filters (`filter_variants_one.sh`), high-depth mask (`build_shared_inputs.sh`) |
| `steps/03_karyotype`, `steps/05_karyo_pca_arms` | Karyotype calling, inversion content, anomalous-diversity block tests |
| `steps/04_angsd_chr2`, `steps/08_demography` | ANGSD SAF/realSFS/thetaStat runs (genome-wide, per division, Z), Stairway Plot 2, moments |
| `steps/09_momentsld`, `steps/10_gone2` | LD-based Ne |
| `steps/11_relernn` | ReLERNN recombination map (`run_one_tag.sh`) |
| `steps/13_diploshic`, `steps/14_sweep_seqmodel` | Sweep simulation (SLiM recipes in `scripts/harness/slim`), diploS/HIC training and scanning, outlier scan (`scripts/harness/outlier_scan.py --mask`), Markov decoder (`calibration/`), background-selection screen (`calibration/tier1_bgs_screen.py`), GO enrichment (`calibration/tier1_go.py`) |
| `steps/17_ensemble_scan` | RAiSD, SweepFinder2 and per-SNP plink2 FST scans |
| `mkado` | McDonald-Kreitman inputs and runs |

The chromosome-2 inversion analysis (msinv fits, drift bounds, load test) is in https://github.com/stsmall/msinv under `illex/` and `results/illex/`.
Figure and table scripts are in `../scripts/`.
