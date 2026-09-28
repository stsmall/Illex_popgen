# tests/test_mapcoord.py
import numpy as np, pandas as pd, pytest
from seqmodel.mapcoord import genetic_positions, deltas, split_blocks

def _map(intervals):  # list of (start,end,rate)
    return pd.DataFrame(intervals, columns=["start","end","rate"]).assign(chrom="1")[["chrom","start","end","rate"]]
def _wins(mids, w=4000):
    return pd.DataFrame({"chrom":"1","start":[m-w//2 for m in mids],"end":[m+w//2 for m in mids]})

def test_cumulative_positions_constant_rate():
    m=_map([(0,100000,1e-8)]); w=_wins([10000,20000,30000])
    g=genetic_positions(w,m,"1")
    assert np.allclose(np.diff(g), 1e-8*10000)   # 10kb steps at 1e-8 => 1e-4 M each

def test_delta_g_and_x():
    m=_map([(0,100000,1e-8)]); w=_wins([10000,20000])
    dg,dx=deltas(w,m,"1"); assert np.allclose(dg,[1e-4]) and np.allclose(dx,[10000])

def test_uncovered_midpoint_raises():
    m=_map([(0,5000,1e-8)]); w=_wins([10000])
    with pytest.raises(ValueError, match="cover|uncovered"):
        genetic_positions(w,m,"1")

def test_split_on_physical_gap():
    m=_map([(0,10_000_000,1e-8)]); w=_wins([100000,104000,5_000_000])  # big jump before 3rd
    blocks=split_blocks(w,m,"1",max_phys_gap=1_000_000)
    assert blocks==[(0,1),(2,2)]

def test_empty_chrom_map_raises():
    import pytest
    m=_map([(0,100000,1e-8)])   # only chrom "1"
    w=_wins([10000])            # query absent chrom "2"
    with pytest.raises(ValueError):
        genetic_positions(w, m, "2")

def test_split_on_internal_map_gap():
    m=_map([(0,100000,1e-8),(200000,300000,1e-8)])  # 100kb hole [100k,200k)
    w=_wins([50000,250000])     # both midpoints covered, 200kb apart
    assert split_blocks(w, m, "1", max_phys_gap=1_000_000)==[(0,0),(1,1)]

def test_windows_multichrom_raises():
    import pandas as pd, pytest
    m=_map([(0,300000,1e-8)])   # chrom "1" only
    w=pd.DataFrame({"chrom":["1","2"],"start":[10000,10000],"end":[14000,14000]})
    with pytest.raises(ValueError):
        genetic_positions(w, m, "1")
    with pytest.raises(ValueError):
        split_blocks(w, m, "1", max_phys_gap=1_000_000)
