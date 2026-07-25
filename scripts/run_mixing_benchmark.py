from pathlib import Path
import gzip
import math
import numpy as np
import pandas as pd

ROOT=Path(r"D:\TISTA_JTH_manuscript\cross_lineage_validation"); OUT=ROOT/"results"; OUT.mkdir(exist_ok=True)
RNG=np.random.default_rng(20260724)
RNG_NULL=np.random.default_rng(20260725)
BMARK=["CD19","MS4A1","CD79A","CD79B","CD22","CD37","CD74","HLA-DRA","BANK1","BLK","FCRL1","FCRL2","CD180","TNFRSF13C","HVCN1"]
TMARK=["CD3D","CD3E","CD3G","TRAC","LCK","MAL","IL7R","LTB","CD247","TRBC1"]

def gene_map():
    m={}
    with gzip.open(ROOT/"Homo_sapiens.gene_info.gz","rt",encoding="utf-8",errors="ignore") as f:
        for line in f:
            if line.startswith("#"): continue
            z=line.rstrip().split("\t")
            for x in z[5].split("|"):
                if x.startswith("Ensembl:ENSG"): m[x.split(":",1)[1]]=z[2]
    return m
def bh(p):
    p=np.asarray(p); n=len(p); o=np.argsort(p); q=np.empty(n); q[o]=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]; return np.minimum(q,1)
def score(x,genes):
    g=[z for z in genes if z in x.index]; a=x.loc[g]; z=a.sub(a.mean(1),axis=0).div(a.std(1,ddof=1).replace(0,np.nan),axis=0); return z.mean(0)
def fit_gene_models(Y,condition,bscore,donor):
    n=Y.shape[0]
    donor_levels=np.unique(donor)
    donor_terms=np.column_stack([(donor==z).astype(float) for z in donor_levels[1:]])
    X0=np.column_stack([np.ones(n),condition,donor_terms]); X1=np.column_stack([X0,bscore]); B0=np.linalg.lstsq(X0,Y,rcond=None)[0]; B1=np.linalg.lstsq(X1,Y,rcond=None)[0]
    E0=Y-X0@B0; E1=Y-X1@B1; v0=np.sum(E0**2,0)/(n-X0.shape[1]); v1=np.sum(E1**2,0)/(n-X1.shape[1]); se0=np.sqrt(v0*np.linalg.inv(X0.T@X0)[1,1]); se1=np.sqrt(v1*np.linalg.inv(X1.T@X1)[1,1]); z0=np.abs(B0[1]/se0); z1=np.abs(B1[1]/se1); erfc=np.vectorize(math.erfc); p0=erfc(z0/np.sqrt(2)); p1=erfc(z1/np.sqrt(2)); return B0[1],B1[1],p0,p1,E0,E1,X1

mapping=gene_map(); tpm=pd.read_csv(ROOT/"GSE107011_TPM.txt.gz",sep="\t",index_col=0,compression="gzip"); tpm.index=[mapping.get(str(x).split('.')[0],"") for x in tpm.index]; tpm=tpm[tpm.index!=""].groupby(level=0).sum(); tpm=tpm.loc[tpm.mean(1)>0.2]
donors=["DZQV","925L","9JD4","G4YW"]
profiles={}
for d in donors:
    bcols=[c for c in tpm.columns if c.startswith(d+"_B_") and "Plasmablast" not in c]
    other=[c for c in tpm.columns if c.startswith(d+"_") and c not in bcols and "PBMC" not in c and "Progenitor" not in c and "Plasmablast" not in c]
    b=tpm[bcols].mean(1).clip(lower=0); bg=tpm[other].mean(1).clip(lower=0); profiles[d]=(b/b.sum(),bg/bg.sum())

rows=[]; example=None
for rep in range(100):
    mats=[]; cond=[]; weights=[]; ds=[]
    for d in donors:
        b,bg=profiles[d]
        # Overlapping ranges prevent the lineage score from becoming a nearly
        # deterministic surrogate for the binary group label.
        for group,(lo,hi) in enumerate([(0.03,0.14),(0.08,0.19)]):
            for _ in range(12):
                w=RNG.uniform(lo,hi); p=((1-w)*bg+w*b); p=p/p.sum(); lib=int(RNG.lognormal(np.log(8e6),.15)); mats.append(RNG.multinomial(lib,p.values)); cond.append(group); weights.append(w); ds.append(d)
    counts=pd.DataFrame(np.asarray(mats).T,index=tpm.index); log=np.log2(counts.add(.5).div(counts.sum(0).add(1),axis=1)*1e6+1); bs=score(log,BMARK).values; ts=score(log,TMARK).values; Y=log.T.values; c=np.asarray(cond,float)
    donor_array=np.asarray(ds)
    b0,b1,p0,p1,e0,e1,X1=fit_gene_models(Y,c,bs,donor_array); q0=bh(p0); q1=bh(p1); top=np.argsort(b0)[-200:]; att=np.nanmedian((b0[top]-b1[top])/b0[top]); pr2=1-np.sum(e1[:,top]**2)/np.sum(e0[:,top]**2)
    cnull=c.copy()
    for d in donors:
        ix=np.where(donor_array==d)[0]; cnull[ix]=RNG_NULL.permutation(cnull[ix])
    _,_,pn0,_,_,_,_=fit_gene_models(Y,cnull,bs,donor_array); qn0=bh(pn0)
    row={"replicate":rep+1,"n_samples":len(c),"cor_true_weight_B_score":np.corrcoef(weights,bs)[0,1],"cor_group_B_score":np.corrcoef(c,bs)[0,1],"adjusted_model_condition_number":np.linalg.cond(X1),"cor_true_weight_T_score":np.corrcoef(weights,ts)[0,1],"DE_FDR05_unadjusted":int(np.sum(q0<.05)),"DE_FDR05_B_score_adjusted":int(np.sum(q1<.05)),"DE_FDR05_donor_stratified_null_labels":int(np.sum(qn0<.05)),"top200_median_attenuation":att,"top200_partial_r2":pr2}
    rows.append(row)
    if rep==0: example=pd.DataFrame({"donor":ds,"composition_group":c.astype(int),"true_B_fraction":weights,"B_score":bs,"T_score":ts})
res=pd.DataFrame(rows); res.to_csv(OUT/"B_cell_composition_mixing_benchmark_replicates.tsv",sep="\t",index=False); example.to_csv(OUT/"B_cell_composition_mixing_benchmark_example_samples.tsv",sep="\t",index=False)
summ=[]
for col in res.columns[2:]:
    summ.append({"metric":col,"median":res[col].median(),"p025":res[col].quantile(.025),"p975":res[col].quantile(.975)})
pd.DataFrame(summ).to_csv(OUT/"B_cell_composition_mixing_benchmark_summary.tsv",sep="\t",index=False)
print(pd.DataFrame(summ).to_string(index=False))
