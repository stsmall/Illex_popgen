#!/usr/bin/env python
"""Supplementary Fig. S10: per-gene McDonald-Kreitman volcano plot, I. illecebrosus vs I. coindetii.

The figure is mkado's own volcano plot, written by mkado_illex/scripts/run_mkado.sh
(MK vcf --per-gene --volcano results/coindetii/coindetii.volcano.pdf, mkado 0.5.0). mkado plots
-log10(NI) against -log10(p) for genes with a finite NI (Dn, Ds, Pn and Ps all > 0, 2,354 of the
20,999 tested) and draws the Bonferroni line over all tested genes. This script rasterises that PDF
at 200 dpi to supp_figures/figS10_mkado.png. Requires Ghostscript (gs).
"""
import subprocess
from pathlib import Path

PDF = Path("/sietch_colab/data_share/illex/popgen_data/mkado_illex/results/coindetii/coindetii.volcano.pdf")
OUT = Path(__file__).resolve().parent.parent / "supp_figures" / "figS10_mkado.png"

subprocess.run(["gs", "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r200",
                "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", f"-sOutputFile={OUT}", str(PDF)],
               check=True)
print(f"wrote {OUT}")
