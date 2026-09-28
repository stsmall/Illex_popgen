#!/usr/bin/env python3
"""
Plot ReLERNN recombination maps for Illex illecebrosus.

Two figures:
  1. Recombination LANDSCAPE  -- one panel per chromosome, recombRate (cM/Mb)
     vs genomic position, with a Male line and a Female line on each.
       autosomes  <- run_male_auto / run_female_auto
       chrZ       <- run_Z_male (ZZ diploid) / run_Z_female (ZO hemizygous)
  2. Overlapping HISTOGRAM across chromosomes -- per-chromosome distribution of
     window rates overlaid, to surface outlier chromosomes. Per-chrom means are
     z-scored; |z| > 2 chromosomes are highlighted + labelled. Male & female
     shown in separate panels.

BSCORRECTED columns: chrom start end nSites recombRate CI95LO CI95HI
recombRate is per-bp (crossover/bp/gen); cM/Mb = recombRate * 1e8.

Auto-discovers whichever maps exist, so it runs on autosomes-only now and on
autosomes+Z once the Z maps finish. Run again to refresh when chrZ lands.
"""
import os, re, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn"
OUT = os.path.join(RD, "figures")
os.makedirs(OUT, exist_ok=True)
CM_PER_MB = 1e8   # per-bp rate -> cM/Mb

MALE_C, FEMALE_C = "#1f77b4", "#d62728"

# map source dirs: (sex, projectDir). Z appended only if its BSCORRECTED exists.
SOURCES = {
    "Male":   ["run_male_auto", "run_Z_male"],
    "Female": ["run_female_auto", "run_Z_female"],
}


def _find_bscorr(projdir):
    hits = glob.glob(os.path.join(RD, projdir, "proj", "*.PREDICT.BSCORRECTED.txt"))
    return hits[0] if hits else None


def _load(projdir):
    f = _find_bscorr(projdir)
    if not f:
        return None
    df = pd.read_csv(f, sep="\t")
    # chrom comes in as b'1' byte-string reprs -> strip to '1'
    df["chrom"] = df["chrom"].astype(str).str.replace(r"^b'|'$", "", regex=True)
    df["cM_Mb"] = df["recombRate"].astype(float) * CM_PER_MB
    df["mid_Mb"] = (df["start"] + df["end"]) / 2.0 / 1e6
    return df[["chrom", "start", "end", "mid_Mb", "nSites", "cM_Mb"]]


def load_sex(dirs):
    parts = [_load(d) for d in dirs]
    parts = [p for p in parts if p is not None]
    return pd.concat(parts, ignore_index=True) if parts else None


def natural_key(c):
    # sort 1,2,...,45 then Z/others last
    m = re.fullmatch(r"\d+", c)
    return (0, int(c)) if m else (1, c)


def rolling_median(y, w=25):
    if len(y) < 3:
        return y
    w = max(3, min(w, len(y) | 1))
    return pd.Series(y).rolling(w, center=True, min_periods=1).median().values


def landscape(maps, chroms):
    ncol = 6
    nrow = int(np.ceil(len(chroms) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 3.0, nrow * 1.9),
                             squeeze=False)
    for i, chrom in enumerate(chroms):
        ax = axes[i // ncol][i % ncol]
        for sex, color in (("Male", MALE_C), ("Female", FEMALE_C)):
            df = maps.get(sex)
            if df is None:
                continue
            d = df[df["chrom"] == chrom].sort_values("mid_Mb")
            if d.empty:
                continue
            ax.plot(d["mid_Mb"], d["cM_Mb"], color=color, lw=0.3, alpha=0.20)
            ax.plot(d["mid_Mb"], rolling_median(d["cM_Mb"].values),
                    color=color, lw=1.1)
        ax.set_title(f"chr{chrom}", fontsize=8, pad=2)
        ax.tick_params(labelsize=6)
        ax.margins(x=0.02)
    for j in range(len(chroms), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.supxlabel("Position (Mb)", fontsize=10)
    fig.supylabel("Recombination rate (cM/Mb)", fontsize=10)
    handles = [Line2D([0], [0], color=MALE_C, lw=2, label="Male"),
               Line2D([0], [0], color=FEMALE_C, lw=2, label="Female")]
    fig.legend(handles=handles, loc="upper right", fontsize=9, ncol=2,
               frameon=False)
    fig.suptitle("Illex recombination landscape (ReLERNN, cM/Mb) — "
                 "male vs female, rolling median", fontsize=11)
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.97))
    p = os.path.join(OUT, "recomb_landscape_male_female.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def outlier_hist(maps, chroms):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharex=True, sharey=True)
    summary = []
    xmax = 0
    for df in maps.values():
        if df is not None:
            xmax = max(xmax, np.nanpercentile(df["cM_Mb"], 99))
    bins = np.linspace(0, xmax, 60)
    for ax, (sex, color) in zip(axes, (("Male", MALE_C), ("Female", FEMALE_C))):
        df = maps.get(sex)
        if df is None:
            ax.set_title(f"{sex}: (no map yet)")
            continue
        means = {c: df[df["chrom"] == c]["cM_Mb"].mean() for c in chroms
                 if not df[df["chrom"] == c].empty}
        # baseline (mu, sd) from AUTOSOMES only, so chrZ doesn't inflate the SD
        # and mask autosomal outliers; z-score every chrom (incl. Z) against it.
        auto = [v for c, v in means.items() if re.fullmatch(r"\d+", c)]
        mu, sd = np.mean(auto), np.std(auto)
        outs = []
        for c in chroms:
            d = df[df["chrom"] == c]["cM_Mb"].dropna()
            if d.empty:
                continue
            z = (means[c] - mu) / sd if sd > 0 else 0.0
            is_out = abs(z) > 2
            ax.hist(d, bins=bins, density=True, histtype="step",
                    lw=1.8 if is_out else 0.5,
                    alpha=0.9 if is_out else 0.20,
                    color="k" if is_out else color,
                    zorder=5 if is_out else 1)
            summary.append((sex, c, means[c], z, is_out))
            if is_out:
                outs.append((means[c], c, z))
        # stagger outlier labels left->right by mean so they don't overlap
        ymax = ax.get_ylim()[1]
        for k, (m, c, z) in enumerate(sorted(outs)):
            ax.axvline(m, color="k", lw=0.6, ls=":")
            ax.text(m, ymax * (0.96 - 0.11 * (k % 4)), f"chr{c} (z={z:+.1f})",
                    fontsize=7.5, rotation=90, va="top", ha="right")
        ax.set_title(f"{sex}  (mean {mu:.3f} cM/Mb; outliers = |z|>2, bold black)")
        ax.set_xlabel("Window recombination rate (cM/Mb)")
    axes[0].set_ylabel("Density")
    fig.suptitle("Per-chromosome recombination-rate distributions — "
                 "outlier chromosomes highlighted", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = os.path.join(OUT, "recomb_hist_outliers.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    sdf = pd.DataFrame(summary, columns=["sex", "chrom", "mean_cM_Mb",
                                         "zscore", "outlier"])
    sdf = sdf.sort_values(["sex", "zscore"])
    sdf.to_csv(os.path.join(OUT, "per_chrom_mean_recomb.tsv"),
               sep="\t", index=False)
    return p, sdf


def main():
    maps = {sex: load_sex(dirs) for sex, dirs in SOURCES.items()}
    have = {s: (m is not None) for s, m in maps.items()}
    print("maps loaded:", have)
    if not any(have.values()):
        print("no BSCORRECTED maps found yet -- nothing to plot")
        return
    chroms = sorted(
        set().union(*[set(m["chrom"]) for m in maps.values() if m is not None]),
        key=natural_key)
    has_Z = any((m is not None) and ("Z" in set(m["chrom"])) for m in maps.values())
    print(f"chromosomes ({len(chroms)}): {chroms}")
    print(f"chrZ present: {has_Z}")
    p1 = landscape(maps, chroms)
    p2, sdf = outlier_hist(maps, chroms)
    print("wrote:", p1)
    print("wrote:", p2)
    print("\nper-chrom mean cM/Mb (outliers flagged):")
    print(sdf.to_string(index=False))


if __name__ == "__main__":
    main()
