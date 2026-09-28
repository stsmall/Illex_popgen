import numpy as np, pandas as pd, pytest
from seqmodel.std_diploshic import standardize, prediction_step

SCHEMA = "config/probability_schema.yaml"

def _write(tmp_path, rows, hdr="chrom\tclassifiedWinStart\tclassifiedWinEnd\tbigWinRange\tpredClass\tprob(neutral)\tprob(likedSoft)\tprob(linkedHard)\tprob(soft)\tprob(hard)"):
    p = tmp_path / "pred.tsv"
    p.write_text(hdr + "\n" + "\n".join(rows) + "\n")
    return str(p)

def test_maps_typo_alias_and_order(tmp_path):
    f = _write(tmp_path, ["1\t0\t4000\t-\tneutral\t0.90\t0.03\t0.03\t0.02\t0.02"])
    df = standardize(f, SCHEMA)
    assert list(df.columns) == ["chrom","start","end","p_neutral","p_linked_hard","p_linked_soft","p_hard","p_soft"]
    assert df.chrom.iloc[0] == "1" and df.start.iloc[0] == 0 and df.end.iloc[0] == 4000
    assert abs(df.p_linked_soft.iloc[0] - 0.03) < 1e-9   # from 'likedSoft'

def test_renormalizes_small_drift(tmp_path):
    f = _write(tmp_path, ["1\t0\t4000\t-\tneutral\t0.90\t0.03\t0.03\t0.02\t0.0200001"])
    df = standardize(f, SCHEMA)
    assert abs(df[["p_neutral","p_linked_hard","p_linked_soft","p_hard","p_soft"]].iloc[0].sum() - 1.0) < 1e-9

def test_hardfail_bad_rowsum(tmp_path):
    f = _write(tmp_path, ["1\t0\t4000\t-\tneutral\t0.90\t0.30\t0.30\t0.20\t0.20"])  # sums to 1.9
    with pytest.raises(ValueError, match="row sum"):
        standardize(f, SCHEMA)

def test_hardfail_missing_column(tmp_path):
    f = _write(tmp_path, ["1\t0\t4000\t-\tneutral\t0.9\t0.1"], hdr="chrom\tclassifiedWinStart\tclassifiedWinEnd\tbigWinRange\tpredClass\tprob(neutral)\tprob(other)")
    with pytest.raises(ValueError, match="p_linked_soft|missing|alias"):
        standardize(f, SCHEMA)

def test_prediction_step(tmp_path):
    rows = [f"1\t{s}\t{s+4000}\t-\tneutral\t0.9\t0.03\t0.03\t0.02\t0.02" for s in range(0, 40000, 4000)]
    df = standardize(_write(tmp_path, rows), SCHEMA)
    assert prediction_step(df) == 4000

def test_prediction_step_multichrom(tmp_path):
    rows = ([f"1\t{s}\t{s+4000}\t-\tneutral\t0.9\t0.03\t0.03\t0.02\t0.02" for s in range(0,20000,4000)]
            + [f"2\t{s}\t{s+4000}\t-\tneutral\t0.9\t0.03\t0.03\t0.02\t0.02" for s in range(0,20000,4000)])
    df = standardize(_write(tmp_path, rows), SCHEMA)
    assert prediction_step(df) == 4000
