import numpy as np, json
WD="/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic"
_A=json.load(open(f"{WD}/ascertainment.json"))
_D=np.array([d for d,_ in _A["depth_hist"]]); _P=np.array([p for _,p in _A["depth_hist"]])
_HD=dict((d,q) for d,q in _A["het_dropout_by_depth"])
# precomputed vectorized dropout lookup (faster than dict.get per element)
_max_d=int(_D.max())
_drop_arr=np.array([_HD.get(d,0.0) for d in range(_max_d+1)])
def _read_ms(path):
    lines=open(path).read().split("\n"); i=[k for k,l in enumerate(lines) if l.startswith("segsites")][0]
    s=int(lines[i].split()[1]); pos=lines[i+1].split()[1:]
    haps=[l for l in lines[i+2:i+2+ (0 if s==0 else 10**9)] if l and set(l)<= set("01")]
    return s, pos, haps
def ascertain_ms(in_path, out_path, rng):
    s,pos,haps=_read_ms(in_path)
    if s==0:
        open(out_path,"w").write(f"//\nsegsites: 0\n\n"); return
    G=np.frombuffer("".join(haps).encode(), dtype=np.uint8).reshape(len(haps),-1) - ord("0")  # (nhap, s)
    nhap=G.shape[0]; dip=G[:nhap//2*2].reshape(-1,2,s)        # (ndip,2,s)
    geno=dip.sum(1)                                           # 0/1/2 alt count
    depth=_D[rng.choice(len(_D), size=geno.shape, p=_P)]      # per-genotype depth
    miss=depth<2
    het=geno==1
    drop=_drop_arr[depth.astype(int)]
    het_to_hom=het & (rng.random(geno.shape)<drop)
    geno=np.where(het_to_hom, 2*rng.integers(0,2,geno.shape), geno)  # het->0 or 2
    geno=np.where(miss,-1,geno)
    # site filters: MAF>0.01 over non-missing, and polymorphic
    ac=np.where(geno<0,0,geno).sum(0); an=2*(geno>=0).sum(0)
    maf=np.minimum(ac,an-ac)/np.maximum(an,1)
    keep=(maf>0.01)&(ac>0)&(ac<an)
    geno=geno[:,keep]; pos=[p for p,k in zip(pos,keep) if k]
    # re-emit diploid as haplotype pairs; missing -> '?'
    a_arr=np.where(geno==-1, ord("?"), np.where(geno==2, ord("1"), ord("0"))).astype(np.uint8)
    b_arr=np.where(geno==-1, ord("?"), np.where(geno>=1,  ord("1"), ord("0"))).astype(np.uint8)
    rows_a=[a_arr[d].tobytes().decode() for d in range(a_arr.shape[0])]
    rows_b=[b_arr[d].tobytes().decode() for d in range(b_arr.shape[0])]
    with open(out_path,"w") as f:
        f.write("//\n"); f.write(f"segsites: {geno.shape[1]}\n")
        f.write("positions: "+" ".join(pos)+"\n")
        for a,b in zip(rows_a,rows_b): f.write(a+"\n"+b+"\n")
def _worker(args):
    raw, out = args
    import numpy as np, os
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if os.path.exists(out):
        return
    ascertain_ms(raw, out, np.random.default_rng(hash(raw)%2**32))
def ascertain_all():
    import glob, os, multiprocessing as mp
    pairs=[]
    for raw in glob.glob(f"{WD}/sims/raw/*/*.ms"):
        out=raw.replace("/raw/","/asc/"); os.makedirs(os.path.dirname(out),exist_ok=True)
        pairs.append((raw, out))
    with mp.Pool(processes=30) as pool:
        pool.map(_worker, pairs)
