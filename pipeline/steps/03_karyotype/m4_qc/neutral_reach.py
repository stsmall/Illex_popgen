#!/usr/bin/env python
"""Neutral probability that 2a, founded at frequency p0 t_inv generations ago, is present today at 0.602 +- 0.02
(review M1). Frequency paths follow the Wright-Fisher diffusion of the msinv trajectory models (Supplementary Methods):
steps of 10 generations, variance p(1-p) dt / (2 Ne), absorbing at 0 and 1, under the neutral baseline of the
inversion fit (Ne 507,288 before 382,957 generations ago, then exponential growth to 9,112,696 at present).
A path is accepted when it is neither lost nor fixed and ends within 0.02 of 0.602, as in the reference tables.
Writes analysis/steps/03_karyotype/m4_qc/neutral_reach.tsv (p0, t_inv, accepted fraction, paths).
"""
import numpy as np
import pandas as pd

N_ANC, N_NOW, T_GROW = 507_288.0, 9_112_696.0, 382_957.0
P_NOW, TOL, DT = 0.602, 0.02, 10
P0 = np.array([0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5])
T_INV = [100_000, 200_000, 390_000, 700_000, 1_000_000, 1_500_000]
N_PATHS = 50_000
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/m4_qc/neutral_reach.tsv"
R_GROW = np.log(N_NOW / N_ANC) / T_GROW


def ne(tau):
    """Effective size tau generations before the present."""
    return N_ANC if tau > T_GROW else N_NOW * np.exp(-R_GROW * tau)


rng = np.random.default_rng(20261113)
rows = []
for t_inv in T_INV:
    p = np.repeat(P0, N_PATHS).astype(np.float64)
    for tau in range(t_inv, 0, -DT):
        live = (p > 0.0) & (p < 1.0)
        sd = np.sqrt(p[live] * (1.0 - p[live]) * DT / (2.0 * ne(tau)))
        p[live] = np.clip(p[live] + sd * rng.standard_normal(live.sum()), 0.0, 1.0)
    ok = (np.abs(p - P_NOW) <= TOL) & (p > 0.0) & (p < 1.0)
    for i, p0 in enumerate(P0):
        acc = ok[i * N_PATHS:(i + 1) * N_PATHS].mean()
        rows.append((p0, t_inv, acc, N_PATHS))
        print(f"p0 {p0:<6} t_inv {t_inv:>9,}  accepted {acc:.5f}", flush=True)
pd.DataFrame(rows, columns=["p0", "t_inv", "accepted", "paths"]).to_csv(OUT, sep="\t", index=False)
