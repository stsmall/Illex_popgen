import json, sys
d = json.load(open(sys.argv[1]))
assert d["n_samples"] == 633, d["n_samples"]
assert 0.0 < d["miss_rate"] < 0.5, d["miss_rate"]
p = sum(x[1] for x in d["depth_hist"]); assert abs(p-1) < 1e-6, p
assert d["depth_hist"][0][0] >= 0 and max(x[0] for x in d["depth_hist"]) <= 60
assert all(0 <= q[1] <= 1 for q in d["het_dropout_by_depth"])
assert d["filters"]["max_f_missing"] == 0.5 and d["filters"]["min_maf"] == 0.01
print("error-model check PASSED")
