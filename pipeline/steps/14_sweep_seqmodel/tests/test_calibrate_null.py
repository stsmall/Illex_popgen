import numpy as np, pandas as pd
from seqmodel.calibrate_null import genome_max, threshold, empirical_p

def test_genome_max_empty_is_zero():
    assert genome_max(pd.DataFrame(columns=["score"]))==0.0

def test_genome_max_value():
    assert genome_max(pd.DataFrame({"score":[1.0,5.0,2.0]}))==5.0

def test_empirical_p_formula():
    nm=[0.,1.,2.,3.,4.]  # B=5
    assert abs(empirical_p(3.5,nm)-(1+1)/(5+1))<1e-12   # one null >=3.5

def test_threshold_higher_quantile():
    nm=list(range(101))  # 0..100
    assert threshold(nm,0.05)>=95
