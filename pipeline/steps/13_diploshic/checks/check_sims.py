import glob, os, sys
WD="/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic"
man=f"{WD}/sims/manifest.tsv"
assert os.path.exists(man)
rows=[l.split() for l in open(man)][1:]
kinds={r[1] for r in rows}; subs={int(r[2]) for r in rows}
assert kinds=={"neutral","sweep"}, kinds
assert subs=={-1}|set(range(11)), subs           # neutral(-1) + 11 sweep positions
for r in rows[:5]:
    txt=open(r[0]).read(); assert txt.count("//")>=1 and "segsites" in txt, r[0]
print(f"sims check PASSED ({len(rows)} reps)")
