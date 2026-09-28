import numpy as np, pandas as pd, yaml
PROB = ["p_neutral","p_linked_hard","p_linked_soft","p_hard","p_soft"]

def _resolve(headers, aliases, name):
    hits = [h for h in aliases if h in headers]
    if len(hits) != 1:
        raise ValueError(f"column '{name}': expected exactly 1 of {aliases} in file, found {hits} (missing/ambiguous alias)")
    return hits[0]

def standardize(pred_path, schema_path):
    schema = yaml.safe_load(open(schema_path))
    raw = pd.read_csv(pred_path, sep="\t")
    cols, hdr = {}, list(raw.columns)
    for name, al in {**schema["coordinate_columns"], **schema["probability_columns"]}.items():
        cols[name] = raw[_resolve(hdr, al, name)]
    df = pd.DataFrame(cols)
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int); df["end"] = df["end"].astype(int)
    P = df[PROB].to_numpy(float)
    if not np.isfinite(P).all() or (P < 0).any() or (P > 1.0001).any():
        raise ValueError("probabilities non-finite or out of [0,1]")
    s = P.sum(1)
    if np.any(np.abs(s - 1.0) > 0.02):
        raise ValueError(f"row sum off by >0.02 (max dev {np.max(np.abs(s-1)):.3g}) — parsing/model error")
    df[PROB] = P / s[:, None]
    return df[["chrom","start","end"] + PROB].sort_values(["chrom","start"]).reset_index(drop=True)

def prediction_step(df):
    d = df.groupby("chrom")["start"].diff().dropna().astype(int)
    d = d[d > 0]
    return int(d.mode().iloc[0]) if len(d) else 0
