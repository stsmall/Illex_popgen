"""Fig 1 map: NW-Atlantic sampling with per-NAFO-division karyotype pies.

Shows sampling geography AND the non-cline (pie composition ~constant across
latitude). Pies sit at the mean sampling latitude/longitude of each division
(Baker et al. 2025 sample metadata). Karyotype pies (AA/AB/BB) sized by sqrt(N).
"""
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/scripts")
from figstyle import apply, KARYO, C, despine
apply()
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import cartopy.crs as ccrs
import cartopy.feature as cfeature

DIV = pd.read_csv("/sietch_colab/data_share/illex/popgen_data/analysis/steps/04_hwe_cline/karyo_by_division.tsv", sep="\t")
# pie positions = mean sampling coordinates of each division (Baker et al. 2025 metadata), so the
# map matches the actual collection sites rather than nominal NAFO-division centroids
META = pd.read_csv("/sietch_colab/data_share/illex/popgen_data/seq_data/baker_2025/docs/Squid_Meta_Sept2024_Simple.txt", sep="\t")
POS = META.groupby("Zone")[["Lat", "Lon"]].mean()
DIV = DIV[DIV["division"].isin(POS.index)].copy()
DIV["lon"] = DIV["division"].map(POS["Lon"]); DIV["mean_lat"] = DIV["division"].map(POS["Lat"])

fig = plt.figure(figsize=(6.4, 7.2))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-79, -45, 33.5, 53], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="#ececec", edgecolor="none")
ax.add_feature(cfeature.COASTLINE.with_scale("50m"), lw=0.5, edgecolor="#9a9a9a")
ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="#f7fbff")
gl = ax.gridlines(draw_labels=True, lw=0.4, color="#dddddd")
gl.top_labels = gl.right_labels = False

def pie_at(ax, lon, lat, fracs, colors, r):
    """draw a karyotype pie (wedges) at map coords; r in degrees."""
    start = 90.0
    for frac, col in zip(fracs, colors):
        theta = 360.0 * frac
        ax.add_patch(Wedge((lon, lat), r, start - theta, start, facecolor=col,
                           edgecolor="white", lw=0.6, transform=ccrs.PlateCarree(),
                           zorder=5))
        start -= theta

Nmax = DIV["N"].max()
for _, row in DIV.iterrows():
    n = row["N"]; tot = row["AA"] + row["AB"] + row["BB"]
    fr = [row["AA"] / tot, row["AB"] / tot, row["BB"] / tot]
    r = 0.55 + 1.1 * np.sqrt(n / Nmax)   # radius in degrees
    pie_at(ax, row["lon"], row["mean_lat"], fr,
           [KARYO["AA"], KARYO["AB"], KARYO["BB"]], r)
    lab = row["division"] if row["N"] >= 5 else f'{row["division"]} (n={int(row["N"])})'
    ax.text(row["lon"], row["mean_lat"] - r - 0.35, lab, ha="center",
            va="top", fontsize=7.5, color=C["ink"], transform=ccrs.PlateCarree(), zorder=6)

# legend
from matplotlib.patches import Patch
leg = [Patch(fc=KARYO["AA"], ec="white", label="AA (standard)"),
       Patch(fc=KARYO["AB"], ec="white", label="AB"),
       Patch(fc=KARYO["BB"], ec="white", label="BB (inverted)")]
ax.legend(handles=leg, loc="lower right", fontsize=8, frameon=True,
          framealpha=0.9, edgecolor="none", title="chr2 karyotype")
ax.set_title("Sampling and chr2 arrangement frequency", loc="left")

out = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/figures/fig1_map.png"
fig.savefig(out, dpi=220, bbox_inches="tight")
print("wrote", out, "| divisions:", ",".join(DIV["division"]))
