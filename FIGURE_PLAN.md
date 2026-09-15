# Figure & table plan — consolidated

Single source of truth for the figure/table architecture (supersedes scattered notes).
Framing: **Fig 1 recapitulates Baker et al.** (PCA/structure + population-genomic stats) so the reader
sees we reproduce it; **everything else is new** (the inversion's evolutionary characterisation + the
selection scan). **Every analysis we ran is represented somewhere** (main or supplement) — nothing dropped.

## MAIN TEXT

### Figures
- **Fig 1 — Baker recapitulation (sample & structure).** *Naive/confirmatory.*
  - A. **Autosomal PCA, chr2 EXCLUDED** (shows: no structure) — *new computation*
  - B. **Z PCA** — *new computation*
  - C. **Z coverage → sex** (342♂/291♀) — exists (`steps/12_sex/sex_assignment.png`)
  - (whole-genome π/θ/Tajima's D in Table 1; **pairwise FST matrix** among divisions as a small panel or Table)
- **Fig 2 — The inversion.**
  - A. **chr2 PCA, 3 stripes** (AA/AB/BB clusters) — data ready (`chr2_karyo_coords.tsv`)
  - B. **FST across chr2** — split from current inversion-object panel a, text→caption
  - C. **per-arrangement diversity** (π_AA, π_BB, dxy) — split from panel b, text→caption
- **Fig 3 — Demography.**
  - A. **Population-size history composite** (Stairway2 + moments + GONE2) — *build composite*
  - B. **Inversion frequency history** (current fitted-history fig: trajectory + Δχ² panels), text→caption
- **Fig 4 — Recombination (single composite, WITH chr2).**
  - A. genome-wide autosomal landscape, **sex-averaged (M+F combined)**, chr2 inversion shaded — use `mapdf.tsv`
  - B. **chrZ male vs female** (sex-biased; kept separate)
- **Fig 5 — Selection.**
  - A. **pipeline schematic** (sims → diploSHIC 5-state classifier → HMM decode → empirical-outlier + concordance)
  - B. **genome-wide sweep Manhattan** + candidates (from `genome.preds`, `outlier_scan/windows.tsv`)

### Tables
- **Table 1** sample & whole-genome summary — DONE
- **Table 2** inversion model/params — DONE
- **Table 3** sweep candidates, grouped by function, Z as rows — DONE

## SUPPLEMENT (all diagnostics + everything tested)
- Fig S: **per-chromosome PCA panels** (each chrom) — *new computation*
- Fig S: **per-chromosome recombination panels**, M+F combined (autosomes), Z separate; **chr2 per-karyotype** — from ReLERNN maps
- Fig S: **mkado McDonald–Kreitman** α results + **DoS scatter plots** (SF3B4 etc.) — from `mkado_illex/results/{coindetii,argentinus}/*`
- Fig S: joint 2-D SFS (colourbar redraw: no interpolation, floor at 1, mask 0-cells) — vendored, redraw pending
- Fig S: LD triangle — DONE (vendored)
- Fig S: windowed PCA (coords to fix) — DONE (vendored)
- Fig S: cline / non-clinality — DONE (vendored)
- Fig S: breakpoint zooms (inversion-object panels c/d) — split pending
- Fig S: confusion matrix — DONE (vendored)
- Table S1 software/versions — DONE
- Table S2 karyotype × NAFO division (geography) — DONE
- Table S3 full candidate set (typeset the outlier_scan tsvs) — pointer done, typeset pending
- Table S: **pairwise FST matrix** among divisions (if not a main panel)

## Generation status / queue
| Item | Data ready? | Needs |
|---|---|---|
| Fig 4 recomb composite | ✅ `mapdf.tsv` + Z maps | plot script |
| Fig 5B Manhattan | ✅ genome.preds + windows.tsv | plot script |
| Fig 2A chr2 PCA | ✅ chr2_karyo_coords.tsv | plot script |
| Fig 2B/C, Fig 3B, S breakpoints, S jSFS redraw | ✅ (msinv) | edit `figures_inversion.py` |
| Fig 3A demography composite | ✅ stairway_out + moments + GONE2 | build composite |
| Fig 5A pipeline schematic | n/a (diagram) | draw (TikZ or matplotlib) |
| Fig 1A/B autosomal + Z PCA, per-chrom PCA, FST matrix | ❌ **new computation** | run PCA/FST (pg_gpu/plink) |
| S mkado α + DoS scatter | ✅ mkado results | plot script |
