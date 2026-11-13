"""Supp Fig S14 (bottom row) -- individual posterior trajectories of the chr2 inversion.

Balancing model (B) with gene flux. Left: frequency of the derived 2a arrangement against
generations since the onset of selection (t_sel), for 200 trajectories re-simulated from
NPE posterior draws, with pointwise 50% and 95% bands and median of 400. Aligning on the
onset removes the spread in timing that smears the pointwise bands of main-text Fig. 4b.
Right: posterior draws of the barrier age t_inv and the onset t_sel behind those
trajectories (dashed lines are medians).

Data (written by illex/abc/npe/plot_npe.py in the msinv repository, branch
illex/trajectory-abc): data/inversion_abc/traj_aligned_npe_f1.csv, traj_params_npe_f1.csv.
Output: supp_figures/figS14c_abc_trajectory.png
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from figstyle import apply, C, OKABE, despine, save  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

apply()
COL = OKABE[7]     # same accent as Fig. 4b
P_NOW = 0.602  # 2a frequency over the 619 karyotyped individuals (karyotype_final 223/299/97)

al = pd.read_csv(ROOT / "data" / "inversion_abc" / "traj_aligned_npe_f1.csv")
pr = pd.read_csv(ROOT / "data" / "inversion_abc" / "traj_params_npe_f1.csv")

fig = plt.figure(figsize=(7.8, 3.0))
gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1], wspace=0.3)

ax = fig.add_subplot(gs[0, 0])
u = al["u_gen"] / 1e3
for c in [c for c in al.columns if c.startswith("path")]:
    ax.plot(u, al[c], lw=0.35, color=C["ink"], alpha=0.10, zorder=1)
ax.fill_between(u, al["q025"], al["q975"], color=COL, alpha=0.18, lw=0, zorder=2, label="95%")
ax.fill_between(u, al["q25"], al["q75"], color=COL, alpha=0.32, lw=0, zorder=3, label="50%")
ax.plot(u, al["q50"], lw=2.0, color=COL, zorder=4, label="median")
ax.plot([], [], lw=0.8, color=C["ink"], alpha=0.4, label="200 trajectories")
ax.axhline(P_NOW, lw=0.8, ls="--", color=C["ink"])
ax.axvline(0, lw=0.8, ls=":", color=C["muted"])
ax.set_xlim(-10, 40)
ax.set_ylim(0, 0.8)
ax.set_yticks([0, 0.2, 0.4, P_NOW])
ax.set_yticklabels(["0", "0.2", "0.4", f"{P_NOW:.3f}"])
ax.set_xlabel("Thousand generations since onset of selection")
ax.set_ylabel("Frequency of derived\n(2a) arrangement")
ax.legend(loc="center right", fontsize=8)
despine(ax)

ax2 = fig.add_subplot(gs[0, 1])
bins = np.linspace(0, 1500, 61)
ax2.hist(pr["t_inv"] / 1e3, bins=bins, color=C["muted"], alpha=0.6,
         label="barrier forms ($t_\\mathrm{inv}$)")
ax2.hist(pr["t_sel"] / 1e3, bins=bins, histtype="step", color=COL, lw=1.6,
         label="selection begins ($t_\\mathrm{sel}$)")
for v, c in ((np.median(pr["t_inv"]) / 1e3, C["muted"]), (np.median(pr["t_sel"]) / 1e3, COL)):
    ax2.axvline(v, color=c, lw=1, ls="--")
ax2.set_xlim(1500, 0)
ax2.set_yticks([])
ax2.set_xlabel("Thousand generations ago")
ax2.legend(fontsize=7.5, loc="upper left")
despine(ax2, keep=("bottom",))

save(fig, str(ROOT / "supp_figures" / "figS14c_abc_trajectory.png"))
print("median t_inv, t_sel (kgen):", np.median(pr["t_inv"]) / 1e3, np.median(pr["t_sel"]) / 1e3,
      "| median aligned path reaches 0.55 at u =", float(u[np.argmax(al["q50"].values >= 0.55)]), "kgen",
      "| lower 2.5% bound reaches 0.55 at u =", float(u[np.argmax(al["q025"].values >= 0.55)]), "kgen")
