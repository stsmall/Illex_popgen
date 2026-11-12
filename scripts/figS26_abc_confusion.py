"""Supp Fig S26 -- confusion matrix of the ABC random-forest model choice for the chromosome-2 inversion.

Grouped five-fold cross-validation of the post-hoc-weighted random-forest classifier over the four frequency
histories with gene flux (neutral drift N, directional selection still under way D1 or stopped D2, balancing
selection B), 16,899 simulations. Each row is a true model, and the cells give the share of its importance-weighted
mass assigned to each predicted model, so rows sum to 1.

Data (written by illex/abc/abc_infer.R, model-choice section, on the base2 reference table in the msinv repository,
branch illex/trajectory-abc, results/illex/modelchoice_base2/abc_rf): data/inversion_abc/abcrf_confusion_cv.csv.
Output: supp_figures/figS26_abc_confusion.png
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from figstyle import apply, save  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

apply()
c = pd.read_csv(ROOT / "data" / "inversion_abc" / "abcrf_confusion_cv.csv")
order = ["N_1", "D1_1", "D2_1", "B_1"]
lab = {"N_1": "N", "D1_1": "D1", "D2_1": "D2", "B_1": "B"}
m = c.pivot(index="true", columns="predicted", values="w").loc[order, order]
m = m.div(m.sum(axis=1), axis=0)
acc = np.trace(m.values) / len(order)

fig, ax = plt.subplots(figsize=(3.4, 3.0))
im = ax.imshow(m.values, cmap="Blues", vmin=0, vmax=1)
for i in range(len(order)):
    for j in range(len(order)):
        v = m.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8, color="white" if v > 0.55 else "#222222")
ax.set_xticks(range(len(order)), [lab[k] for k in order])
ax.set_yticks(range(len(order)), [lab[k] for k in order])
ax.set_xlabel("Predicted model")
ax.set_ylabel("True model")
ax.set_title(f"Cross-validated accuracy {acc:.2f}", fontsize=9)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("Share of true-model mass")
save(fig, ROOT / "supp_figures" / "figS26_abc_confusion.png")
