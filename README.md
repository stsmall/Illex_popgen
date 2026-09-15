# Manuscript — *Illex illecebrosus* population genomics

LaTeX source for the paper. Self-contained: all figures are vendored into
`figures/`, so it builds and uploads (GitHub, bioRxiv) without external paths.

## Build
```bash
make          # pdflatex + bibtex + 2× pdflatex  ->  main.pdf
make clean    # remove build artifacts
```
Requires `pdflatex` and `bibtex` (TeX Live). `main.pdf` is ~16 pp.

## Files
| File | Purpose |
|---|---|
| `main.tex` | the manuscript (Abstract → Intro → Results → Discussion → Limitations → Methods) |
| `references.bib` | tool + literature citations |
| `software_table.tex` | Methods software/versions table (`\input` by main.tex) |
| `figures/` | the 7 vendored figures (fig1–fig7) |
| `Makefile` | build rules |

## Provenance & sync
- **Source of record is `../ILLEX_PROJECT_MANUSCRIPT_AND_METHODS.md`** (the working
  document with all results/reasoning). This LaTeX is a rendering of it; keep the
  numbers in sync with that file, not the reverse.
- Figure sources and their original locations are in `../FIGURES_AND_TABLES_INDEX.md`.
- Framing: the paper's contribution is the **evolutionary characterisation of the
  chr2 inversion** and the **genome-wide selective-sweep scan**. The
  PCA/FST/admixture structure result **confirms** Baker et al. and is treated as a
  starting point, not a finding.

## Before submission — required
1. **`references.bib`** — complete every entry (DOI/pages) and resolve the rows
   marked `CONFIRM` (esp. `baker2025`: fill the real Baker et al. citation; GONE2,
   SLiM 5, tskit, Pfam release, SQANTI3 version).
2. **Author list, affiliations, acknowledgements, funding** — placeholders in `main.tex`.
3. **Abstract** — tighten to the target journal's word limit.
4. **Numbers** — every value is transcribed from the markdown/audited results; do a
   final pass against `../ILLEX_PROJECT_MANUSCRIPT_AND_METHODS.md` (the 33-vs-34
   concordant count, etc.).
5. Consider a supplementary methods section for the sweep-scan calibration (§8b/§8d).
