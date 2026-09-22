"""Supplementary sampling table: every collection event (division x date x depth) with source, depth, the
number of individuals analysed here (all 633 sequenced) and the number retained by Baker et al. 2025 (540).
Sources and Baker Ns transcribed from Baker et al. 2025 Table 1; our Ns from their deposited metadata
restricted to the 633 sequenced individuals (ENA PRJEB101095)."""
import os, pandas as pd
A = "/sietch_colab/data_share/illex/popgen_data"
_HERE = os.path.dirname(os.path.abspath(__file__)); _ROOT = os.path.dirname(_HERE)
m = pd.read_csv(f"{A}/seq_data/baker_2025/docs/Squid_Meta_Sept2024_Simple.txt", sep="\t", dtype=str)
ids = set(open(f"{A}/analysis/steps/00_callset/baker_633.txt").read().split())
m = m[m.ID.isin(ids)].copy()
m["date"] = pd.to_datetime(m.Collection_Date, format="%d-%b-%y", errors="coerce")
m.loc[m.date.isna(), "date"] = pd.to_datetime(m.Collection_Date, format="%b-%y", errors="coerce")
m["depth"] = pd.to_numeric(m.Depth, errors="coerce")
m.loc[(m.Zone == "3K") & (m.depth == 4), "depth"] = 5          # Baker Table 1 publishes 5 m for this jig haul
SRC = {"3N": "Canadian trawl survey", "3O": "Canadian trawl survey", "4R": "Canadian trawl survey",
       "4S": "Canadian trawl survey", "4T": "Canadian trawl survey", "4W": "Canadian trawl survey",
       "4X": "Canadian trawl survey", "6A": "US trawl survey", "6B": "Commercial trawl", "6C": "US trawl survey"}
BAKER = {("3K","2023-08-16",3):48, ("3K","2023-08-30",5):46, ("3N","2023-10-20",160):1, ("3O","2023-10-13",204):7,
 ("3O","2023-10-14",613):1, ("3O","2023-10-14",350):1, ("4R","2023-08-07",113):24, ("4R","2023-08-25",245):7,
 ("4R","2023-08-25",233):4, ("4R","2023-08-25",173):2, ("4R","2023-08-28",221):3, ("4S","2023-08-28",228):2,
 ("4S","2023-08-29",251):1, ("4S","2023-08-30",80):1, ("4S","2023-08-30",124):2, ("4S","2023-08-31",155):3,
 ("4S","2023-08-31",152):1, ("4S","2023-09-01",296):1, ("4S","2023-09-01",273):1, ("4S","2023-09-03",243):2,
 ("4T","2023-09-07",27):1, ("4T","2023-09-10",81):1, ("4T","2023-09-10",101):1, ("4T","2023-09-15",50):1,
 ("4T","2023-09-22",39):1, ("4T","2023-09-22",142):1, ("4T","2023-09-22",158):1, ("4T","2023-09-23",341):2,
 ("4T","2023-09-28",62):1, ("4T","2023-09-28",56):1, ("4T","2023-09-29",49):1, ("4W","2023-07-16",84):43,
 ("4X","2022-07",None):76, ("6A","2022-09-22",120):84, ("6B","2023-06-18",74):46, ("6B","2023-08-16",201):48,
 ("6C","2022-09-14",307):73}
g = m.groupby(["Zone", "date", "depth"], dropna=False).size().reset_index(name="n")
rows = []
for r in g.itertuples():
    ds = r.date.strftime("%Y-%m") if r.Zone == "4X" else r.date.strftime("%Y-%m-%d")
    dp = None if pd.isna(r.depth) else int(r.depth)   # Baker truncates half-metre depths
    src = ("Recreational jig" if ds == "2023-08-16" else "Coastal fishery jig") if r.Zone == "3K" else SRC[r.Zone]
    rows.append(dict(_d=r.date, division=r.Zone, date=(r.date.strftime("%B %Y") if r.Zone == "4X" else r.date.strftime("%-d %B %Y")),
                     source=src, depth_m=("not recorded" if dp is None else dp), n_this_study=r.n,
                     n_baker2025=BAKER.get((r.Zone, ds, dp), 0), _lat=m[m.Zone == r.Zone].Lat.astype(float).mean()))
t = pd.DataFrame(rows).sort_values(["_lat", "_d", "depth_m"], ascending=[False, True, True], key=lambda c: c.astype(str).str.zfill(4) if c.name=="depth_m" else c).drop(columns=["_lat","_d"])
os.makedirs(f"{_ROOT}/tables_export", exist_ok=True)
t.to_csv(f"{_ROOT}/tables_export/TableS8_sampling.csv", index=False)
tex = ["\\begin{longtable}{llllrr}",
 "\\caption{\\textbf{Sampling events.} Every collection event (NAFO division, date, depth) with its source, the "
 "number of individuals analysed in this study (all 633 whole-genome sequenced individuals) and the number retained "
 "by \\citet{baker2025} after their genotype-missingness filter (540). Divisions are ordered north to south by mean "
 "sampling latitude. 4X was bulk-sampled across the division with no per-tow depth or position recorded.}"
 "\\label{stab:sampling}\\\\", "\\toprule", "Division & Date & Source & Depth (m) & $N$ (this study) & $N$ (Baker) \\\\",
 "\\midrule", "\\endfirsthead", "\\toprule", "Division & Date & Source & Depth (m) & $N$ (this study) & $N$ (Baker) \\\\",
 "\\midrule", "\\endhead"]
prev = None
for r in t.itertuples():
    tex.append(f"{r.division if r.division != prev else ''} & {r.date} & {r.source} & {r.depth_m} & {r.n_this_study} & {r.n_baker2025} \\\\")
    prev = r.division
tex += ["\\midrule", f"Total & & & & {t.n_this_study.sum()} & {t.n_baker2025.sum()} \\\\", "\\bottomrule", "\\end{longtable}"]
open(f"{_ROOT}/sampling_table.tex", "w").write("\n".join(tex) + "\n")
print(len(t), "events | totals", t.n_this_study.sum(), t.n_baker2025.sum())
print(t.groupby("division")[["n_this_study", "n_baker2025"]].sum().to_string())
