import numpy as np

def genome_max(regions_df):
    if regions_df is None or len(regions_df)==0 or "score" not in regions_df: return 0.0
    return float(np.max(regions_df["score"].to_numpy(float)))

def threshold(null_maxima, alpha):
    nm=np.asarray(list(null_maxima),float)
    try: return float(np.quantile(nm,1-alpha,method="higher"))
    except TypeError: return float(np.quantile(nm,1-alpha,interpolation="higher"))

def empirical_p(obs, null_maxima):
    nm=np.asarray(list(null_maxima),float); B=len(nm)
    return (1+int(np.sum(nm>=obs)))/(B+1)
