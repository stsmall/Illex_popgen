import pytest
from seqmodel.std_relernn import standardize, integrated_length_morgans

HDR = "chrom\tstart\tend\tnSites\trecombRate\tCI95LO\tCI95HI"
def _w(tmp_path, rows):
    p = tmp_path/"m.txt"; p.write_text(HDR+"\n"+"\n".join(rows)+"\n"); return str(p)

def test_strips_bytes_chrom_and_selects_which(tmp_path):
    f=_w(tmp_path, ["b'1'\t0\t1000\t50\t2.0e-9\t1.0e-9\t3.0e-9"])
    assert standardize(f,"corrected").rate.iloc[0]==2.0e-9
    assert standardize(f,"lower").rate.iloc[0]==1.0e-9
    assert standardize(f,"upper").rate.iloc[0]==3.0e-9
    assert standardize(f).chrom.iloc[0]=="1"

def test_rejects_overlap(tmp_path):
    f=_w(tmp_path, ["b'1'\t0\t1000\t5\t2e-9\t1e-9\t3e-9","b'1'\t500\t1500\t5\t2e-9\t1e-9\t3e-9"])
    with pytest.raises(ValueError, match="overlap"):
        standardize(f)

def test_rejects_ci_order(tmp_path):
    f=_w(tmp_path, ["b'1'\t0\t1000\t5\t2e-9\t3e-9\t1e-9"])  # lo>hi
    with pytest.raises(ValueError, match="CI"):
        standardize(f)

def test_integrated_length(tmp_path):
    f=_w(tmp_path, ["b'1'\t0\t1000000\t5\t2e-9\t1e-9\t3e-9"])  # 1Mb at 2e-9 = 2e-3 Morgans
    assert abs(integrated_length_morgans(standardize(f),"1") - 2e-3) < 1e-12
