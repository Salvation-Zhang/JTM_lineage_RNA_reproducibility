from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import landscape, A4
from reportlab.lib.colors import HexColor, white

REPO=Path(__file__).resolve().parents[1]
OUT=REPO/"figures"; OUT.mkdir(parents=True,exist_ok=True)
PDF=OUT/"Figure_1_study_design_JTM.pdf"
W,H=landscape(A4); c=canvas.Canvas(str(PDF),pagesize=(W,H))
navy=HexColor("#183B56"); blue=HexColor("#2D7DD2"); green=HexColor("#2E8B57"); red=HexColor("#C94C4C"); gold=HexColor("#E6A23C"); gray=HexColor("#667085"); pale=HexColor("#F4F7FA")

def txt(x,y,s,size=8,color=navy,bold=False,anchor="left"):
    c.setFillColor(color); c.setFont("Helvetica-Bold" if bold else "Helvetica",size)
    if anchor=="center": c.drawCentredString(x,y,s)
    elif anchor=="right": c.drawRightString(x,y,s)
    else: c.drawString(x,y,s)
def box(x,y,w,h,title,lines,color):
    c.setFillColor(white); c.setStrokeColor(color); c.setLineWidth(1.5); c.roundRect(x,y,w,h,7,fill=1,stroke=1)
    txt(x+w/2,y+h-20,title,10,color,True,"center")
    for i,line in enumerate(lines): txt(x+w/2,y+h-39-i*13,line,7.5,gray,False,"center")
def arrow(x1,y1,x2,y2,color=gray):
    c.setStrokeColor(color); c.setLineWidth(1.2); c.line(x1,y1,x2,y2); c.line(x2,y2,x2-6,y2+3); c.line(x2,y2,x2-6,y2-3)

c.setFillColor(pale); c.roundRect(25,25,W-50,H-50,10,fill=1,stroke=0)
txt(42,H-48,"A",14,navy,True); txt(70,H-46,"Study design and externally specified lineage-aware framework",13,navy,True)

# Evidence sources
box(45,H-155,185,78,"Primary clinical cohort",["GSE112278: eltrombopag-treated ITP","17 patients; 46 whole-blood samples","baseline, week 1, month 1"],blue)
box(45,H-275,185,78,"External score references",["GSE302674 + GSE262073","purified platelet RNA-seq","12-gene platelet RNA score"],green)
box(45,H-395,185,78,"External validation",["GSE186294: EPO / erythroid recovery","GSE112594: rituximab / B-cell depletion","fixed canonical lineage scores"],red)
box(45,H-515,185,78,"Known-truth benchmark",["GSE107011 purified immune profiles","donor-matched composition-only mixtures","100 replicates + null-label calibration"],gold)

# Central analysis
box(315,H-205,210,105,"Primary platelet composition audit",["paired score change and exact tests","top-gene coupling + matched sets","effect attenuation + conditional partial R2","partially nested held-out reconstruction"],blue)
box(315,H-365,210,105,"Transport across lineages",["target-lineage score direction","overlap-free partial R2","matched-set comparison","nested held-out reconstruction"],red)
box(315,H-500,210,80,"Composition-only falsification test",["donor-adjusted DE before vs after score","model condition number","donor-stratified null labels"],gold)

arrow(230,H-116,315,H-155,blue); arrow(230,H-236,315,H-180,green)
arrow(230,H-356,315,H-315,red); arrow(230,H-476,315,H-460,gold)

# Outputs and interpretation
box(610,H-175,185,95,"Primary clinical finding",["platelet score rises in responders","72.1% median effect attenuation","partial R2 = 0.771","held-out R2 = 0.702"],blue)
box(610,H-330,185,95,"Cross-lineage evidence",["EPO / erythroid partial R2 = 0.356","rituximab / B-cell partial R2 = 0.198","positive held-out R2 in both cohorts"],red)
box(610,H-485,185,95,"Known-truth result",["median 3,037 unadjusted DE genes","median 0 after B-cell-score adjustment","null labels: median 0 discoveries"],gold)
arrow(525,H-153,610,H-128,blue); arrow(525,H-313,610,H-283,red); arrow(525,H-460,610,H-438,gold)

txt(W/2,48,"Interpretation: lineage-derived RNA is biologically meaningful, but bulk treatment signatures do not by themselves establish cell-intrinsic reprogramming.",9,navy,True,"center")
txt(W-38,31,"RNA scores are abundance proxies, not measured cell counts.",6.5,gray,False,"right")
c.save(); print(PDF)
