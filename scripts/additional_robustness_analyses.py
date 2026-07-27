from pathlib import Path
import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("TISTA_PROCESSED_DATA", REPO / "data" / "processed"))
OUT = REPO / "results" / "primary"
OUT.mkdir(parents=True, exist_ok=True)

MARKERS = ["PPBP","PF4","GP1BA","GP9","ITGA2B","ITGB3","TUBB1","TREML1","RGS18","SDPR","SPARC","CLU"]

def exact_between(r, n):
    values=np.r_[r,n]; k=len(r); obs=np.mean(r)-np.mean(n); null=[]
    for ix in itertools.combinations(range(len(values)),k):
        m=np.zeros(len(values),bool); m[list(ix)]=True
        null.append(values[m].mean()-values[~m].mean())
    extreme=np.sum(np.abs(null)>=abs(obs)-1e-12)
    return (extreme+1)/(len(null)+1)

def exact_signflip(x):
    obs=abs(np.mean(x)); null=[]
    for signs in itertools.product([-1,1],repeat=len(x)):
        null.append(abs(np.mean(x*np.asarray(signs))))
    return np.mean(np.asarray(null)>=obs-1e-12)

def score(logcpm, genes):
    genes=[g for g in genes if g in logcpm.index]
    z=logcpm.loc[genes].sub(logcpm.loc[genes].mean(axis=1),axis=0).div(logcpm.loc[genes].std(axis=1,ddof=1),axis=0)
    return z.mean(axis=0)

def paired_delta(meta, values):
    x=meta[["gsm","patient","response","time"]].copy(); x["score"]=x.gsm.map(values)
    w=x.pivot(index=["patient","response"],columns="time",values="score").dropna(subset=["pre","1wk"])
    w["delta"]=w["1wk"]-w["pre"]
    return w.reset_index()

def main():
    meta=pd.read_csv(DATA/"GSE112278_metadata.tsv",sep="\t")
    counts=pd.read_csv(DATA/"GSE112278_counts.tsv",sep="\t",index_col=0)[meta.gsm]
    lib=counts.sum(axis=0); logcpm=np.log2(counts.add(.5).div(lib.add(1),axis=1)*1e6+1)
    ranks=pd.read_csv(DATA/"external_platelet_marker_ranks.tsv",sep="\t")
    ranks["primary_symbol"]=ranks.gene.replace({"CAVIN2":"SDPR"})
    strict=ranks[(ranks.GSE302674_percentile>=.99)&(ranks.GSE262073_percentile>=.99)].primary_symbol.tolist()
    strict=[g for g in strict if g in MARKERS]
    variants={"Locked 12-gene":MARKERS,"Strict dual-reference":strict}
    for g in MARKERS: variants[f"Leave out {g}"]=[x for x in MARKERS if x!=g]
    all_scores={name:score(logcpm,genes) for name,genes in variants.items()}
    base=all_scores["Locked 12-gene"]
    rows=[]
    for name,s in all_scores.items():
        p=paired_delta(meta,s.to_dict()); r=p[p.response=="responder"].delta.values; n=p[p.response=="nonresponder"].delta.values
        rows.append({"variant":name,"n_genes":len(variants[name]),"responder_mean_delta":r.mean(),"nonresponder_mean_delta":n.mean(),
                     "difference":r.mean()-n.mean(),"between_exact_p":exact_between(r,n),"responder_signflip_p":exact_signflip(r),
                     "sample_score_correlation_with_locked":np.corrcoef(s,base)[0,1]})
    out=pd.DataFrame(rows); out.to_csv(OUT/"signature_definition_robustness.tsv",sep="\t",index=False)
    audit=pd.read_csv(DATA/"GSE112278_week1_composition_adjusted_gene_audit.tsv",sep="\t")
    top=audit.nlargest(208,"responder_mean_delta")
    overlap=set(MARKERS)&set(top.gene)
    no_overlap=top[~top.gene.isin(MARKERS)]
    pd.DataFrame([{
        "top_set_n":len(top),"locked_score_genes_in_top_set":len(overlap),"overlap_genes":",".join(sorted(overlap)),
        "overlap_free_n":len(no_overlap),"overlap_free_median_correlation":no_overlap.delta_score_correlation.median(),
        "overlap_free_median_attenuation":no_overlap.attenuation_fraction.median()
    }]).to_csv(OUT/"mathematical_overlap_audit.tsv",sep="\t",index=False)
    ranks.to_csv(OUT/"external_reference_marker_transportability.tsv",sep="\t",index=False)

    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":600})
    fig,axs=plt.subplots(1,3,figsize=(11,4))
    q=ranks[ranks.primary_symbol.isin(MARKERS)].copy(); y=np.arange(len(q))
    axs[0].scatter(q.GSE302674_percentile,y,label="GSE302674",s=28); axs[0].scatter(q.GSE262073_percentile,y,label="GSE262073",s=28,marker="s")
    axs[0].set_yticks(y,q.primary_symbol); axs[0].invert_yaxis(); axs[0].set_xlim(.89,1.003); axs[0].set_xlabel("Within-reference expression percentile"); axs[0].legend(frameon=False,fontsize=7); axs[0].set_title("External reference support")
    loo=out[out.variant.str.startswith("Leave out")]
    axs[1].scatter(loo.responder_mean_delta,loo.variant.str.replace("Leave out ","",regex=False),label="Responder")
    axs[1].scatter(loo.nonresponder_mean_delta,loo.variant.str.replace("Leave out ","",regex=False),label="Nonresponder")
    axs[1].axvline(0,color="#A0AEC0",lw=.8); axs[1].set_xlabel("Paired score change"); axs[1].set_title("Leave-one-marker-out stability"); axs[1].legend(frameon=False,fontsize=7)
    axs[2].scatter(loo.sample_score_correlation_with_locked,loo.variant.str.replace("Leave out ","",regex=False),color="#2F855A")
    axs[2].set_xlim(.95,1.001); axs[2].set_xlabel("Correlation with locked score"); axs[2].set_title("Signature-definition transportability")
    for i,a in enumerate(axs): a.text(-.15,1.04,chr(65+i),transform=a.transAxes,fontweight="bold",fontsize=12)
    fig.tight_layout()
    for ext in ["png","pdf","svg"]: fig.savefig(OUT/f"Supplementary_Figure_signature_robustness.{ext}",bbox_inches="tight")
    plt.close(fig)
    print(out.to_string(index=False))

if __name__=="__main__": main()
