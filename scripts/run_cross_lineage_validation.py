from pathlib import Path
import gzip,re,itertools,math,os
import numpy as np
import pandas as pd

REPO=Path(__file__).resolve().parents[1]
ROOT=Path(os.environ.get("TISTA_CROSS_LINEAGE_DATA", REPO/"data"/"cross_lineage"))
OUT=REPO/"results"/"cross_lineage"; OUT.mkdir(parents=True,exist_ok=True)
RNG=np.random.default_rng(20260724)

ERYTHROID=["ALAS2","AHSP","HBB","HBA1","HBA2","CA1","GYPA","SLC4A1","EPB42","BPGM","FECH","KLF1","ANK1","TMOD1","HEMGN"]
B_CELL=["CD19","MS4A1","CD79A","CD79B","CD22","CD37","CD74","HLA-DRA","BANK1","BLK","FCRL1","FCRL2","CD180","TNFRSF13C","HVCN1"]
T_CELL=["CD3D","CD3E","CD3G","TRAC","LCK","MAL","IL7R","LTB","CD247","TRBC1"]
NEUTROPHIL=["CSF3R","FCGR3B","CEACAM8","S100A8","S100A9","MNDA","FPR1","CXCR2","MMP8","OLFM4"]

def gene_map():
    m={}
    with gzip.open(ROOT/"Homo_sapiens.gene_info.gz","rt",encoding="utf-8",errors="ignore") as f:
        for line in f:
            if line.startswith("#"): continue
            z=line.rstrip("\n").split("\t"); symbol=z[2]
            for x in z[5].split("|"):
                if x.startswith("Ensembl:ENSG"): m[x.split(":",1)[1]]=symbol
    return m

def parse_soft(acc):
    txt=gzip.open(ROOT/f"{acc}_family.soft.gz","rt",encoding="utf-8",errors="ignore").read(); rows=[]
    for s in txt.split("^SAMPLE = ")[1:]:
        d={}
        for key,out in [("!Sample_geo_accession","gsm"),("!Sample_title","title")]:
            q=re.search(rf"{re.escape(key)} = (.+)",s); d[out]=q.group(1) if q else ""
        for q in re.findall(r"!Sample_characteristics_ch1 = (.+)",s):
            if ": " in q: k,v=q.split(": ",1); d[k]=v
        rows.append(d)
    return pd.DataFrame(rows)

def collapse_symbols(counts,mapping):
    counts.index=[mapping.get(str(x).split('.')[0],"") for x in counts.index]
    counts=counts[counts.index!=""]
    return counts.groupby(level=0).sum()

def logcpm(counts): return np.log2(counts.add(.5).div(counts.sum(axis=0).add(1),axis=1)*1e6+1)
def score(x,genes):
    g=[z for z in genes if z in x.index]; z=x.loc[g].sub(x.loc[g].mean(axis=1),axis=0).div(x.loc[g].std(axis=1,ddof=1).replace(0,np.nan),axis=0); return z.mean(axis=0),g
def exact_signflip(x):
    x=np.asarray(x,float); obs=abs(x.mean()); vals=[abs(np.mean(x*np.asarray(s))) for s in itertools.product([-1,1],repeat=len(x))]; return np.mean(np.asarray(vals)>=obs-1e-12)
def exact_between(a,b):
    x=np.r_[a,b]; k=len(a); obs=abs(np.mean(a)-np.mean(b)); vals=[]; total=math.comb(len(x),k)
    if total>200000:
        null=[]
        for _ in range(100000):
            ix=RNG.choice(len(x),k,replace=False); mask=np.zeros(len(x),bool); mask[ix]=1; null.append(abs(x[mask].mean()-x[~mask].mean()))
        return (np.sum(np.asarray(null)>=obs-1e-12)+1)/(len(null)+1)
    for ix in itertools.combinations(range(len(x)),k):
        m=np.zeros(len(x),bool); m[list(ix)]=1; vals.append(abs(x[m].mean()-x[~m].mean()))
    return (np.sum(np.asarray(vals)>=obs-1e-12)+1)/(len(vals)+1)
def partial_r2(Y,X0,X1):
    e0=Y-X0@np.linalg.lstsq(X0,Y,rcond=None)[0]; e1=Y-X1@np.linalg.lstsq(X1,Y,rcond=None)[0]; s0=np.sum(e0**2); s1=np.sum(e1**2); return 1-s1/s0,s0,s1
def medcorr(Y,x):
    yc=Y-Y.mean(0); xc=x-x.mean(); den=np.sqrt(np.sum(yc**2,0)*np.sum(xc**2)); r=np.divide(xc@yc,den,out=np.full(Y.shape[1],np.nan),where=den>0); return np.nanmedian(r),r
def matched_permutation(delta,genes,observed,lineage_delta,n=5000):
    means=delta.mean(axis=1); bins=pd.qcut(means.rank(method="first"),10,labels=False); selected=pd.Index(genes); counts=bins.loc[selected].value_counts().to_dict(); vals=[]
    for _ in range(n):
        pick=[]
        for b,k in counts.items(): pick.extend(RNG.choice(bins.index[bins==b],k,replace=False))
        vals.append(np.nanmedian([np.corrcoef(delta.loc[g].values,lineage_delta)[0,1] for g in pick]))
    return (np.sum(np.asarray(vals)>=observed-1e-12)+1)/(n+1),np.nanmedian(vals)

mapping=gene_map()

# EPO: one platform; Base2 to EPO4 around day 10-16.
epo=pd.read_csv(ROOT/"GSE186294_Illumina_counts.txt.gz",sep=r"\s+",index_col=0,compression="gzip")
epo.columns=[str(x).strip('"') for x in epo.columns]; epo=collapse_symbols(epo,mapping); em=parse_soft("GSE186294"); em=em[em.title.str.contains("Illumina")].copy(); em["sample_num"]=em.title.str.extract(r"Sample (\d+)").astype(int); em["column"]=em.sample_num.astype(str); em["subject"]=(em.sample_num.sub(1)//5+1).astype(str); em["visit"]=em["time (day)"].str.extract(r"^(Base1|Base2|EPO3|EPO4|Post7)")
elog=logcpm(epo); escore,eg=score(elog,ERYTHROID); eneg_t,_=score(elog,T_CELL); eneg_n,_=score(elog,NEUTROPHIL)
ew=em[em.visit.isin(["Base2","EPO4"])].pivot(index="subject",columns="visit",values="column").dropna(); epat=ew.index.tolist(); edelta=pd.DataFrame({p:elog[ew.loc[p,"EPO4"]]-elog[ew.loc[p,"Base2"]] for p in epat}); esd=np.array([escore[ew.loc[p,"EPO4"]]-escore[ew.loc[p,"Base2"]] for p in epat])
etop=edelta.mean(axis=1).nlargest(200).index; egenes=[g for g in etop if g not in eg]; EY=edelta.loc[egenes].T.values; eX0=np.ones((len(epat),1)); eX1=np.column_stack([np.ones(len(epat)),esd]); epr2,es0,es1=partial_r2(EY,eX0,eX1); emc,_=medcorr(EY,esd); eperm,enull=matched_permutation(edelta,egenes,emc,esd)

# Fully nested EPO LOPO.
ecv0=ecv1=0
for i,p in enumerate(epat):
    tr=[q for q in epat if q!=p]; top=edelta[tr].mean(1).nlargest(200).index; genes=[g for g in top if g not in eg]; ytr=edelta.loc[genes,tr].T.values; yte=edelta.loc[genes,p].values; x=np.array([esd[epat.index(q)] for q in tr]); x0=np.ones((len(tr),1)); x1=np.column_stack([np.ones(len(tr)),x]); p0=(np.ones((1,1))@np.linalg.lstsq(x0,ytr,rcond=None)[0]).ravel(); p1=(np.array([[1,esd[i]]])@np.linalg.lstsq(x1,ytr,rcond=None)[0]).ravel(); ecv0+=np.sum((yte-p0)**2); ecv1+=np.sum((yte-p1)**2)
ecvr2=1-ecv1/ecv0

# Rituximab: Baseline to week 26, active vs placebo.
bc=pd.read_csv(ROOT/"GSE112594_counts.txt.gz",sep="\t",index_col=0,compression="gzip"); bc=collapse_symbols(bc,mapping); bm=parse_soft("GSE112594"); bm["column"]=bm.title; blog=logcpm(bc); bscore,bg=score(blog,B_CELL); bt,_=score(blog,T_CELL); bn,_=score(blog,NEUTROPHIL)
bw=bm[bm.visit.isin(["Baseline","Wk 26"])].pivot(index=["subject id","treatment"],columns="visit",values="column").dropna(); bpat=[x[0] for x in bw.index]; bgroup=np.array([x[1] for x in bw.index]); bdelta=pd.DataFrame({x[0]:blog[bw.loc[x,"Wk 26"]]-blog[bw.loc[x,"Baseline"]] for x in bw.index}); bsd=np.array([bscore[bw.loc[x,"Wk 26"]]-bscore[bw.loc[x,"Baseline"]] for x in bw.index])
active=[bpat[i] for i in range(len(bpat)) if bgroup[i]=="Active"]; placebo=[bpat[i] for i in range(len(bpat)) if bgroup[i]=="Placebo"]; effect=bdelta[active].mean(axis=1)-bdelta[placebo].mean(axis=1); btop=effect.nsmallest(200).index; bgenes=[g for g in btop if g not in bg]; BY=bdelta.loc[bgenes,bpat].T.values; treat=(bgroup=="Active").astype(float); bX0=np.column_stack([np.ones(len(bpat)),treat]); bX1=np.column_stack([bX0,bsd]); bpr2,bs0,bs1=partial_r2(BY,bX0,bX1); bmc,_=medcorr(BY,bsd); bperm,bnull=matched_permutation(bdelta,bgenes,bmc,bsd)

# Fully nested rituximab LOPO.
bcv0=bcv1=0
for i,p in enumerate(bpat):
    tr=np.arange(len(bpat))!=i; tp=[bpat[k] for k in np.where(tr)[0] if bgroup[k]=="Active"]; cp=[bpat[k] for k in np.where(tr)[0] if bgroup[k]=="Placebo"]; top=(bdelta[tp].mean(1)-bdelta[cp].mean(1)).nsmallest(200).index; genes=[g for g in top if g not in bg]; ytr=bdelta.loc[genes,[bpat[k] for k in np.where(tr)[0]]].T.values; yte=bdelta.loc[genes,p].values; x0=np.column_stack([np.ones(tr.sum()),treat[tr]]); x1=np.column_stack([x0,bsd[tr]]); p0=(np.array([[1,treat[i]]])@np.linalg.lstsq(x0,ytr,rcond=None)[0]).ravel(); p1=(np.array([[1,treat[i],bsd[i]]])@np.linalg.lstsq(x1,ytr,rcond=None)[0]).ravel(); bcv0+=np.sum((yte-p0)**2); bcv1+=np.sum((yte-p1)**2)
bcvr2=1-bcv1/bcv0

# Bootstrap partial R2 and score changes.
def bootstrap(Y,group,sd,n=10000):
    out=[]; levels=np.unique(group); inds=[np.where(group==z)[0] for z in levels]
    for _ in range(n):
        ix=np.concatenate([RNG.choice(q,len(q),replace=True) for q in inds]); g=group[ix]; x=sd[ix]; y=Y[ix]; X0=np.column_stack([np.ones(len(ix))]+[(g==z).astype(float) for z in levels[1:]]); X1=np.column_stack([X0,x]); r,_,_=partial_r2(y,X0,X1); out.append([r]+[x[g==z].mean() for z in levels])
    return np.asarray(out),levels
eb,el=bootstrap(EY,np.array(["EPO"]*len(epat)),esd); bb,bl=bootstrap(BY,bgroup,bsd)

summary=[]
def add(dataset,lineage,n,score_change,p,pr2,boot,mc,perm,null,cvr2,extra=""):
    lo,hi=np.nanpercentile(boot[:,0],[2.5,97.5]); summary.append({"dataset":dataset,"lineage":lineage,"paired_n":n,"lineage_score_change":score_change,"primary_test_p":p,"partial_r2":pr2,"bootstrap_partial_r2_low":lo,"bootstrap_partial_r2_high":hi,"median_gene_correlation":mc,"matched_set_empirical_p":perm,"matched_set_null_median":null,"nested_lopo_predictive_r2":cvr2,"notes":extra})
add("GSE186294","erythroid",len(epat),esd.mean(),exact_signflip(esd),epr2,eb,emc,eperm,enull,ecvr2,"Base2 to EPO4; Illumina platform")
ad=bsd[bgroup=="Active"]; pl=bsd[bgroup=="Placebo"]
add("GSE112594","B_cell",len(bpat),ad.mean()-pl.mean(),exact_between(ad,pl),bpr2,bb,bmc,bperm,bnull,bcvr2,f"Active n={len(ad)}; placebo n={len(pl)}; baseline to week 26")
pd.DataFrame(summary).to_csv(OUT/"cross_lineage_composition_dominance_summary.tsv",sep="\t",index=False)

pd.DataFrame({"dataset":"GSE186294","subject":epat,"group":"EPO","lineage_score_change":esd}).to_csv(OUT/"GSE186294_subject_changes.tsv",sep="\t",index=False)
pd.DataFrame({"dataset":"GSE112594","subject":bpat,"group":bgroup,"lineage_score_change":bsd}).to_csv(OUT/"GSE112594_subject_changes.tsv",sep="\t",index=False)
specificity=pd.DataFrame([
    {"dataset":"GSE186294","contrast":"Base2_to_EPO4","score":"erythroid","mean_change":esd.mean(),"permutation_p":exact_signflip(esd)},
    {"dataset":"GSE186294","contrast":"Base2_to_EPO4","score":"T_cell","mean_change":np.mean([eneg_t[ew.loc[p,"EPO4"]]-eneg_t[ew.loc[p,"Base2"]] for p in epat]),"permutation_p":exact_signflip([eneg_t[ew.loc[p,"EPO4"]]-eneg_t[ew.loc[p,"Base2"]] for p in epat])},
    {"dataset":"GSE186294","contrast":"Base2_to_EPO4","score":"neutrophil","mean_change":np.mean([eneg_n[ew.loc[p,"EPO4"]]-eneg_n[ew.loc[p,"Base2"]] for p in epat]),"permutation_p":exact_signflip([eneg_n[ew.loc[p,"EPO4"]]-eneg_n[ew.loc[p,"Base2"]] for p in epat])},
])
for label,values in [("B_cell",bsd),("T_cell",np.array([bt[bw.loc[x,"Wk 26"]]-bt[bw.loc[x,"Baseline"]] for x in bw.index])),("neutrophil",np.array([bn[bw.loc[x,"Wk 26"]]-bn[bw.loc[x,"Baseline"]] for x in bw.index]))]:
    aa=values[bgroup=="Active"]; pp=values[bgroup=="Placebo"]
    specificity.loc[len(specificity)]={"dataset":"GSE112594","contrast":"active_minus_placebo_change","score":label,"mean_change":aa.mean()-pp.mean(),"permutation_p":exact_between(aa,pp)}
specificity.to_csv(OUT/"cross_lineage_score_specificity.tsv",sep="\t",index=False)
pd.DataFrame({"signature":"erythroid","gene":pd.Series(eg)}).to_csv(OUT/"locked_erythroid_signature.tsv",sep="\t",index=False)
pd.DataFrame({"signature":"B_cell","gene":pd.Series(bg)}).to_csv(OUT/"locked_B_cell_signature.tsv",sep="\t",index=False)
print(pd.DataFrame(summary).to_string(index=False))
