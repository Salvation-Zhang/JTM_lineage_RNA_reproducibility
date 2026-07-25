from pathlib import Path
import csv
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import landscape, A4
from reportlab.lib.colors import HexColor, black, white

ROOT=Path(r"D:\TISTA_JTH_manuscript\upgrade_analysis")
RES=ROOT/"results"; OUT=ROOT/"output"/"pdf"; OUT.mkdir(parents=True,exist_ok=True)
PDF=ROOT.parent/"figures"/"Figure_5_composition_dominance.pdf"
PDF.parent.mkdir(parents=True,exist_ok=True)

def read_tsv(name):
    with open(RES/name,encoding="utf-8") as f: return list(csv.DictReader(f,delimiter="\t"))
boot={r["estimand"]:r for r in read_tsv("patient_bootstrap_and_variance_attribution.tsv")}
var=read_tsv("composition_variance_attribution.tsv")[0]
sens=read_tsv("young_platelet_signature_size_sensitivity.tsv")
month=read_tsv("complete_case_month1.tsv")

BLUE=HexColor("#2B6CB0"); RED=HexColor("#C53030"); GREY=HexColor("#718096"); LIGHT=HexColor("#EDF2F7"); DARK=HexColor("#1A202C")
w,h=landscape(A4); c=canvas.Canvas(str(PDF),pagesize=(w,h)); c.setTitle("Composition dominance and sensitivity analyses")

def panel(x,y,pw,ph,label,title):
    c.setStrokeColor(HexColor("#CBD5E0")); c.roundRect(x,y,pw,ph,6,stroke=1,fill=0)
    c.setFillColor(DARK); c.setFont("Helvetica-Bold",12); c.drawString(x+10,y+ph-18,label)
    c.setFont("Helvetica-Bold",10); c.drawString(x+30,y+ph-18,title)
def axes(x,y,pw,ph,ymin,ymax):
    l=x+48; b=y+38; r=x+pw-16; t=y+ph-38; c.setStrokeColor(GREY); c.line(l,b,l,t); c.line(l,b,r,b)
    for i in range(5):
        yy=b+(t-b)*i/4; val=ymin+(ymax-ymin)*i/4; c.setFont("Helvetica",7); c.setFillColor(GREY); c.drawRightString(l-5,yy-2,f"{val:.1f}"); c.setStrokeColor(HexColor("#E2E8F0")); c.line(l,yy,r,yy)
    return l,b,r,t

margin=28; gap=14; pw=(w-2*margin-gap)/2; ph=(h-2*margin-gap)/2

# A Bootstrap estimates and held-out reconstruction R2
x=margin; y=margin+ph+gap; panel(x,y,pw,ph,"A","Variance attribution and patient-bootstrap uncertainty"); l,b,r,t=axes(x,y,pw,ph,0,1.2)
items=[("overlap_free_multivariate_partial_r2","Partial R^2"),("overlap_free_median_gene_correlation","Median correlation"),("overlap_free_median_attenuation","Median attenuation")]
for i,(k,lab) in enumerate(items):
    q=boot[k]; v=float(q["estimate"]); lo=float(q["bootstrap_95ci_low"]); hi=float(q["bootstrap_95ci_high"]); xx=l+(r-l)*(i+.5)/3
    c.setStrokeColor(BLUE); c.setLineWidth(1.3); c.line(xx,b+(t-b)*lo/1.2,xx,b+(t-b)*hi/1.2); c.line(xx-4,b+(t-b)*lo/1.2,xx+4,b+(t-b)*lo/1.2); c.line(xx-4,b+(t-b)*hi/1.2,xx+4,b+(t-b)*hi/1.2); c.setFillColor(BLUE); c.circle(xx,b+(t-b)*v/1.2,3,fill=1,stroke=0); c.setFillColor(DARK); c.setFont("Helvetica",7); c.drawCentredString(xx,b-13,lab)
c.setFont("Helvetica",7.5); c.drawString(l,b-27,f"Nested LOPO held-out transcriptome R^2 = {float(var['nested_lopo_incremental_predictive_r2']):.3f}; genes reselected within each fold.")

# B young platelet sensitivity
x=margin+pw+gap; y=margin+ph+gap; panel(x,y,pw,ph,"B","Young-platelet residual sensitivity"); l,b,r,t=axes(x,y,pw,ph,-0.5,0.5)
q=[z for z in sens if z["score"]=="age_residual" and z["residualization"]=="baseline_fitted"]
for grp,color,key in [("Responder",BLUE,"responder_mean_change"),("Nonresponder",RED,"nonresponder_mean_change")]:
    pts=[]
    for z in q:
        size=float(z["requested_tail_size"]); v=float(z[key]); xx=l+(r-l)*(size-25)/75; yy=b+(t-b)*(v+.5); pts.append((xx,yy)); c.setFillColor(color); c.circle(xx,yy,3,fill=1,stroke=0)
    c.setStrokeColor(color); c.setLineWidth(1.4)
    for p1,p2 in zip(pts,pts[1:]): c.line(*p1,*p2)
for size in [25,50,100]: c.setFillColor(DARK); c.setFont("Helvetica",7); c.drawCentredString(l+(r-l)*(size-25)/75,b-13,str(size))
c.setFont("Helvetica",7.5); c.setFillColor(BLUE); c.drawString(l,b-27,"Responder"); c.setFillColor(RED); c.drawString(l+62,b-27,"Nonresponder"); c.setFillColor(DARK); c.drawRightString(r,b-27,"Genes requested per tail")

# C month 1 complete cases
x=margin; y=margin; panel(x,y,pw,ph,"C","Complete-case month-1 platelet-score changes"); l,b,r,t=axes(x,y,pw,ph,-0.8,0.8)
contrasts=["1mon_minus_pre","1mon_minus_1wk"]
for i,con in enumerate(contrasts):
    for j,(grp,color) in enumerate([("responder",BLUE),("nonresponder",RED)]):
        z=next(v for v in month if v["contrast"]==con and v["group"]==grp); val=float(z["mean_change"]); center=l+(r-l)*(i+.5)/2+(j-.5)*28; zero=b+(t-b)*.5; yy=b+(t-b)*(val+.8)/1.6; c.setFillColor(color); c.rect(center-11,min(zero,yy),22,abs(yy-zero),fill=1,stroke=0); c.setFillColor(DARK); c.setFont("Helvetica",6.5); c.drawCentredString(center,yy+(5 if val>=0 else -9),f"n={z['n_complete']}")
c.setFont("Helvetica",7); c.drawCentredString(l+(r-l)*.25,b-14,"Month 1 - baseline"); c.drawCentredString(l+(r-l)*.75,b-14,"Month 1 - week 1")

# D framework
x=margin+pw+gap; y=margin; panel(x,y,pw,ph,"D","Operational composition-dominance framework")
steps=["1  External lineage score","2  Within-patient coupling","3  Effect attenuation and partial R^2","4  Lineage and matched-set controls","5  Patient-exclusion stability","6  Multiplicity-corrected residual test"]
bx=x+58; bw=pw-92; top=y+ph-48
for i,s in enumerate(steps):
    yy=top-i*29; c.setFillColor(LIGHT); c.setStrokeColor(GREY); c.roundRect(bx,yy-13,bw,21,5,fill=1,stroke=1); c.setFillColor(DARK); c.setFont("Helvetica",8); c.drawString(bx+10,yy-6,s)
    if i<5: c.setStrokeColor(GREY); c.line(bx+bw/2,yy-13,bx+bw/2,yy-21)
c.setFont("Helvetica-Oblique",7); c.drawString(x+20,y+16,"A bulk response is termed composition-dominant only when all six checks are concordant.")

c.setFont("Helvetica",7); c.setFillColor(GREY); c.drawRightString(w-margin,12,"Exploratory sensitivity analysis; patient is the resampling unit.")
c.save(); print(PDF)
