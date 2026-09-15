"""Fig 1 map: NW-Atlantic sampling with per-NAFO-division karyotype pies.

Shows sampling geography AND the non-cline (pie composition ~constant across
latitude). Latitude is the empirical division mean_lat; longitude is the NAFO
standard-zone centroid (approximate). Karyotype pies (AA/AB/BB) sized by sqrt(N).
"""
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, KARYO, C, despine
apply()
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import cartopy.crs as ccrs
import cartopy.feature as cfeature

DIV = pd.read_csv("/sietch_colab/data_share/illex/popgen_data/analysis/steps/04_hwe_cline/karyo_by_division.tsv", sep="\t")
# approximate NAFO division centroid longitudes (deg W); latitude uses data mean_lat
LON = {"3K": -52.0, "3N": -49.5, "3O": -52.5, "4R": -59.5, "4S": -62.0,
       "4T": -63.5, "4W": -60.5, "4X": -66.0, "6A": -73.0, "6B": -73.8, "6C": -75.0}
DIV = DIV[DIV["division"].isin(LON)].copy()
DIV["lon"] = DIV["division"].map(LON)

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
    ax.text(row["lon"], row["mean_lat"] - r - 0.35, row["division"], ha="center",
            va="top", fontsize=7.5, color=C["ink"], transform=ccrs.PlateCarree(), zorder=6)

# legend
from matplotlib.patches import Patch
leg = [Patch(fc=KARYO["AA"], ec="white", label="AA (standard)"),
       Patch(fc=KARYO["AB"], ec="white", label="AB"),
       Patch(fc=KARYO["BB"], ec="white", label="BB (inverted)")]
ax.legend(handles=leg, loc="lower right", fontsize=8, frameon=True,
          framealpha=0.9, edgecolor="none", title="chr2 karyotype")
ax.set_title("Sampling and chr2 arrangement frequency", loc="left")

out = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/figures/fig_map.png"
fig.savefig(out, dpi=220, bbox_inches="tight")
print("wrote", out, "| divisions:", ",".join(DIV["division"]))
