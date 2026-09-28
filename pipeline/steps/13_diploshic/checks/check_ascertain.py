import sys, numpy as np, json, glob
sys.path.insert(0,"/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic")
import ascertain as A
WD="/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic"
asc=json.load(open(f"{WD}/ascertainment.json"))
raw=sorted(glob.glob(f"{WD}/sims/raw/neutral/*.ms"))[0]
rng=np.random.default_rng(0); A.ascertain_ms(raw, "/tmp/asc_test.ms", rng)
txt=open("/tmp/asc_test.ms").read()
miss=txt.count("?")/max(1,sum(c in "01?" for c in txt))
assert 0.02 < miss < 0.5, miss                       # realistic missingness injected
s_raw=int([l for l in open(raw) if l.startswith("segsites")][0].split()[1])
s_asc=int([l for l in open("/tmp/asc_test.ms") if l.startswith("segsites")][0].split()[1])
assert s_asc <= s_raw, (s_asc, s_raw)                # filtering only removes sites
print(f"ascertain check PASSED (miss={miss:.3f}, segsites {s_raw}->{s_asc})")
