import json,os
WD="/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic"
assert os.path.exists(f"{WD}/model.json") and os.path.exists(f"{WD}/model.weights.h5")
m=json.load(open(f"{WD}/test_metrics.json"))
assert m["accuracy"]>0.7, m["accuracy"]                       # 5-class, chance=0.2
assert min(m["per_class_recall"].values())>0.3, m["per_class_recall"]  # no class collapse
print(f"model check PASSED (acc {m['accuracy']:.3f})")
