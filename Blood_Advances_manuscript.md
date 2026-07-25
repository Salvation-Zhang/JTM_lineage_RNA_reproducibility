# Therapy-driven lineage shifts can dominate whole-blood treatment transcriptomes

Ye Zhang^1^ and Jun Peng^1,2^

^1^ Department of Hematology, Qilu Hospital of Shandong University, Jinan, Shandong 250012, China  
^2^ Shandong Provincial Key Laboratory of Immunohematology, Qilu Hospital of Shandong University, Jinan, China

Correspondence: Jun Peng, Department of Hematology, Qilu Hospital of Shandong University, Jinan, Shandong 250012, China; email: junpeng@163.com

## Short title

Lineage shifts in blood transcriptomes

## Abstract

Whole-blood transcriptomes combine cell-intrinsic regulation with changes in circulating-cell and cell-derived RNA abundance, creating an identifiability problem when therapy directly alters a blood lineage. We reanalyzed 46 longitudinal whole-blood RNA-sequencing samples from 17 eltrombopag-treated patients with immune thrombocytopenia using an externally defined 12-gene platelet RNA score. The score increased from baseline to week 1 in responders (n=11; mean standardized change, 0.886) but not nonresponders (n=4; −0.371; exact between-group P=0.0381). Among the 208 genes most induced in responders, median correlation with platelet-score change was 0.916, and score adjustment attenuated the median response effect by 72.1%; both exceeded 5,000 expression-matched gene sets (empirical P=0.00020). After direct score-gene overlap was removed, platelet-score change explained 77.1% of residual variation in this selected gene set (bootstrap 95% CI, 52.8%–88.8%) and reduced fully nested leave-one-patient-out reconstruction error by 70.2%. The framework transported to erythropoietin-associated erythroid recovery (n=10; partial R2=0.356) and rituximab-associated B-cell depletion (n=35; partial R2=0.198). In 100 composition-only experiments generated from purified immune-cell profiles, the B-cell score tracked the known B-cell fraction (median r=0.967). Donor-adjusted conventional testing identified a median 3,037 FDR-significant genes despite fixed within-cell profiles; lineage-score adjustment reduced this to 0 (2.5th–97.5th percentile, 0–2). Therapy-driven lineage shifts can therefore generate broad whole-blood treatment signatures without requiring cell-intrinsic reprogramming. Lineage-aware analysis should be routine when interventions alter blood-cell composition.

### Keywords

immune thrombocytopenia; eltrombopag; rituximab; erythropoietin; whole-blood transcriptomics; cell composition

## Key Points

- Platelet RNA recovery accounted for most of the selected responder-associated whole-blood gene set during eltrombopag treatment in ITP.
- Cross-lineage cohorts and composition-only mixtures show why treatment signatures require explicit lineage-aware interpretation.

## Introduction

Whole-blood transcriptomics provides a practical window into systemic disease and treatment response, but its signals arise from two sources unless they are modeled explicitly: altered transcription within cells and altered abundance of circulating cells or cell-derived RNA.[1] The latter can be biologically meaningful while still confounding mechanistic interpretation. This problem is particularly acute when a therapy restores or depletes a blood lineage whose RNA contributes to the assayed specimen.

Immune thrombocytopenia (ITP) is characterized by immune-mediated platelet destruction and impaired platelet production. Eltrombopag, a thrombopoietin-receptor agonist, increases platelet production and can produce marked changes in circulating platelet abundance. Platelets are anucleate but contain a diverse and dynamically regulated RNA repertoire.[2] Consequently, a whole-blood transcriptional signature observed during successful eltrombopag treatment could report platelet recovery itself, immune-cell reprogramming, or both.

We therefore reanalyzed a longitudinal whole-blood transcriptomic cohort of eltrombopag-treated ITP using a lineage-aware compositional framework. To reduce outcome-driven circularity, we defined the platelet RNA proxy from two independent purified-platelet datasets rather than from treatment-associated genes in the primary cohort. We then transported the same analytical logic to erythropoietin-associated erythroid recovery and rituximab-associated B-cell depletion and tested it in composition-only mixtures generated from purified human immune-cell profiles. This design distinguishes a dataset-specific platelet observation from a testable property of bulk blood transcriptomes.

## Methods

### Study design and datasets

The primary cohort was GSE112278, comprising 46 whole-blood 3′ RNA-sequencing samples from 17 patients with chronic ITP sampled before eltrombopag, after 1 week, and after 1 month.[3] Response groups were those assigned in the source study using criteria modified from the International Working Group guidelines; we did not redefine response from the transcriptomic data. In the source definition, complete response required a platelet count >100 × 10^9/L and no bleeding; response required a count of 30–100 × 10^9/L, at least a two-fold increase from baseline, and no bleeding; nonresponse was defined as a count <30 × 10^9/L, failure to double the baseline count within 90 days of treatment, or bleeding. Available samples included 12, 11, and 9 responder samples and 5, 4, and 5 nonresponder samples at the three time points, respectively. Missing visits were retained in longitudinal patient-fixed-effect models; paired week-1 analyses used the 15 patients with both baseline and week-1 samples. The source report noted that three responders received intravenous immunoglobulin 4–8 days before pretreatment RNA collection; sample-linked platelet counts at the RNA draws were unavailable, and no count values were reconstructed from the published figure.

GSE302674 and GSE262073 were used as independent purified-platelet RNA-sequencing references.[9,10] GSE43177 and GSE46922 were used for exploratory analyses in purified peripheral-blood T cells.[11,12] GSE196676, containing 56,312 bone-marrow CD34-positive cells from four newly diagnosed untreated patients with ITP and four healthy donors, was analyzed using donor-level pseudobulk.[13]

For cross-lineage validation, GSE186294 contributed paired Illumina whole-blood profiles from 10 healthy men before erythropoietin and at the expected erythroid response window (Base2 to EPO4, approximately days 10–16).[4] GSE112594 contributed paired baseline and week-26 whole-blood profiles from 25 rituximab-treated and 10 placebo-treated participants with new-onset type 1 diabetes.[5,6] Fixed 15-gene erythroid and B-cell signatures were specified from canonical lineage markers without selecting genes on treatment effects in the validation cohorts. Scores were interpreted as lineage-derived RNA-abundance proxies, not measured cell counts.

GSE107011 supplied RNA-sequencing profiles of FACS-purified immune populations from four donors for an independent composition-only benchmark.[7] For each donor, purified B-cell profiles were mixed with the donor-matched mean non-B immune-cell profile. Low and high composition groups were sampled from overlapping B-cell-fraction ranges (3%–14% and 8%–19%), with 12 samples per group per donor, variable library size, and multinomial sampling. One hundred independent experiments were generated; by construction, no within-cell transcriptional change occurred.

### External platelet RNA score

A platelet marker panel was evaluated in both purified-platelet references, and a 12-gene score detectable in the primary cohort was fixed before testing longitudinal outcomes: PPBP, PF4, GP1BA, GP9, ITGA2B, ITGB3, TUBB1, TREML1, RGS18, SDPR/CAVIN2, SPARC, and CLU. Within GSE112278, normalized expression values were standardized gene-wise and averaged to obtain a sample-level platelet RNA score. The score was interpreted as an RNA-abundance proxy and not as an estimate of platelet count or platelet fraction.

### Primary paired analysis

For each patient with baseline and week-1 samples, we calculated the change in platelet RNA score. Mean changes were summarized separately in responders and nonresponders. The responder change was tested using an exact sign-flip test, and the between-group difference used exact label permutation.

### Gene-level composition audit

For every gene, the patient-level week-1 minus baseline log-counts-per-million change was modeled as a function of response group and change in the external platelet RNA score. The 208 genes with the largest mean induction among responders were used to summarize coupling to platelet-score change and attenuation of the response coefficient. To test whether these summaries merely reflected expression level and gene-set size, 5,000 random 208-gene sets were sampled after matching by mean-expression decile.

### Variance attribution and patient-level uncertainty

To quantify platelet-associated RNA beyond response-group membership, we formed a patient-by-gene matrix of baseline-to-week-1 changes for the 208 genes most induced in responders and removed the 10 genes directly shared with the fixed platelet score. A reduced multivariate model included an intercept and response group; the full model additionally included patient-level platelet-score change. Incremental multivariate partial R2 was the proportional reduction in residual sum of squares after adding platelet-score change. Uncertainty in the responder change, between-group difference, partial R2, median gene-level coupling, and median coefficient attenuation was estimated using 10,000 stratified bootstrap resamples, with patients as the resampling unit.

For fully nested leave-one-patient-out transcriptome reconstruction, responder-induced genes were reselected within each training fold, direct platelet-score overlap was removed, reduced and platelet-adjusted models were fitted, and expression changes were reconstructed for the held-out patient using that patient's observed platelet RNA-score change. Reconstruction errors were aggregated across folds, and incremental held-out R2 was calculated from the reduction in residual sum of squares. This analysis evaluated transport of the composition association to held-out transcriptomes; it was not a clinical-outcome prediction model.

### Specificity and technical sensitivity

Prespecified erythroid, neutrophil, monocyte, T-cell, and B-cell scores served as negative controls. Their paired changes were analyzed with the same exact tests used for the platelet score. A patient-fixed-effect model tested the responder-by-week-1 effect on the platelet score before and after simultaneous adjustment for log library size and all five nonplatelet scores.

### Cross-lineage transport and composition-only benchmark

In GSE186294, paired Base2-to-EPO4 changes were tested by exact sign flipping. In GSE112594, active-minus-placebo differences in paired baseline-to-week-26 change were tested by label permutation; because exhaustive enumeration was infeasible, P values used 100,000 seeded permutations. For each dataset, the 200 most directionally responsive genes were selected, direct signature overlap was removed, and incremental multivariate partial R2 was estimated after adding the lineage score to a reduced treatment model. Patient bootstrap confidence intervals used 10,000 resamples, matched-expression gene-set tests used 5,000 permutations, and fully nested leave-one-patient-out analyses repeated gene selection within each training fold.

In the mixing benchmark, conventional group differential-expression models were compared with models additionally containing the fixed B-cell score; both included donor fixed effects. Across 100 replicates, we recorded correlation with the known B-cell fraction, group–score correlation, adjusted-model condition number, the number of genes at FDR<0.05, top-200 effect attenuation, and partial R2. Overlapping fraction ranges were used to avoid near-deterministic collinearity between group and score. As a null calibration, group labels were independently permuted within each donor in every replicate and the donor-adjusted differential-expression analysis was repeated.

The primary clinical hypothesis concerned platelet-associated RNA recovery in GSE112278. Cross-lineage transport and the composition-only benchmark were external validation analyses. Because the top responsive gene sets were selected from the same cohort-level contrasts subsequently summarized, in-sample coupling, attenuation, and partial R2 are conditional descriptive estimates rather than unbiased genome-wide effect estimates. Fully nested leave-one-patient-out analyses, which repeated gene selection within each training fold, were used to evaluate held-out transcriptome reconstruction. These analyses did not predict treatment response because the response label and lineage-score change were included as model inputs.

### Full longitudinal and stability analyses

Gene-wise models included patient fixed effects, week-1 and month-1 indicators, responder-by-time interactions, and the external platelet RNA score. Benjamini–Hochberg correction was applied separately to week-1 and month-1 interaction tests. Leave-one-patient-out analyses repeated the paired composition audit after excluding each of the 15 paired patients.

### Complete-case month-1 and platelet-maturation sensitivity analyses

Month-1 analyses included only patients observed at both visits being contrasted and did not impute missing samples. Paired changes were tested using exact sign-flip tests, and responder-versus-nonresponder differences used exact label permutation.

An exploratory platelet-maturation analysis used GSE126448, containing paired FACS-sorted immature/high-RNA and mature/low-RNA platelets from four donors.[8] Signatures containing up to 25, 50, or 100 genes per tail were selected independently of treatment outcome by requiring consistent direction across all four donor pairs and adequate expression. Young-minus-mature contrasts were evaluated before and after residualization against the fixed platelet RNA score. Residualization was performed both across all samples and using a regression fitted only in baseline samples. These specifications were treated as sensitivity analyses rather than independent confirmatory tests.

### Exploratory secondary analyses

As a secondary audit of the previously proposed metabolic interpretation, 12 Reactome pathways were fixed before analysis.[14] Whole-blood pathway scores were tested after platelet-score adjustment, and purified T-cell (GSE43177 and GSE46922) and bone-marrow CD34-positive data (GSE196676) were examined for directional consistency. Bone-marrow counts were aggregated by donor, never by cell, and exact tests used donor-level allocations. These analyses were not considered replication of eltrombopag response and are reported primarily in the Supplement.

## Results

### An externally defined platelet RNA score rises selectively in responders

The two purified-platelet references ranked the canonical platelet genes near the top of their respective transcriptomes, supporting transport of the externally defined signature into GSE112278. In the 15 patients with paired baseline and week-1 samples, the platelet RNA score increased by a mean of 0.886 standardized units in responders and decreased by 0.371 units in nonresponders. The responder change was supported by exact sign flipping (P=0.00195), and the responder–nonresponder difference was significant by exact label permutation (P=0.0381). Patient-level trajectories showed that the group result was not attributable to a single extreme responder.

### The responder-induced transcriptome is strongly coupled to platelet RNA recovery

Among the 208 genes with the largest week-1 induction in responders, the median patient-level correlation with platelet-score change was 0.916. Including platelet-score change in the paired gene model reduced the median responder coefficient by 72.1%. These observations were not reproduced by expression-matched random gene sets: the null medians were −0.014 for correlation and 0.030 for attenuation, with empirical upper-tail P=0.00020 for both statistics.

### Platelet RNA change explains most variation in the responder-induced transcriptome

After removal of the 10 genes directly shared with the platelet score, 198 responder-induced genes remained for multivariate variance attribution. Addition of platelet-score change to a model already containing response group reduced residual variation by 77.1% (multivariate partial R2=0.771; patient-bootstrap 95% CI, 0.528–0.888). The overlap-free median gene-level correlation was 0.909 (95% CI, 0.733–0.967), and median response-effect attenuation was 71.3%, although its bootstrap interval was wide (8.8%–106.8%). In fully nested leave-one-patient-out analysis, with responder-induced genes reselected within every training fold, adding the held-out patient's observed platelet-score change reduced transcriptome-reconstruction error by 70.2%. This was not clinical-response prediction.

The responder platelet-score increase was 0.886 standardized units (patient-bootstrap 95% CI, 0.496–1.390), and the responder-versus-nonresponder difference was 1.257 units (95% CI, 0.078–2.413), emphasizing both the stability of the direction and the uncertainty imposed by four paired nonresponders.

### The longitudinal change is platelet-specific and technically robust

None of the five nonplatelet lineage scores showed both a significant responder increase and a responder–nonresponder difference. In contrast, the platelet score remained associated with responder status at week 1 after simultaneous adjustment for library size and erythroid, neutrophil, monocyte, T-cell, and B-cell scores (β=1.334; P=0.00531), similar to the base estimate (β=1.262; P=0.00593).

The composition finding was stable to removing every paired patient. Across 15 leave-one-patient-out folds, the mean responder score change ranged from 0.671 to 0.976, the median correlation of the top 208 genes with score change ranged from 0.863 to 0.935, and median response-effect attenuation ranged from 60.1% to 91.4%.

### The result is insensitive to signature definition and direct gene overlap

To determine whether a single highly abundant platelet transcript drove the result, we recalculated the score after removing each of the 12 markers in turn. Every leave-one-marker-out score remained almost perfectly correlated with the locked score across samples (r=0.997–1.000). Across all variants, the mean week-1 change remained positive in responders (range, 0.857–0.926) and negative in nonresponders (range, −0.400 to −0.346); exact between-group P values ranged from 0.0329 to 0.0483. A stricter 11-gene subset supported at high abundance in both purified-platelet references produced the same pattern (responder change, 0.895; nonresponder change, −0.346; exact P=0.0410; r with locked score=0.999). Ten of the locked score genes occurred among the top 208 responder-induced genes. After removing all 10 overlapping genes, the remaining 198 genes retained a median correlation of 0.916 with platelet-score change and a median attenuation of 77.9%, indicating that direct mathematical overlap did not explain the result.

### No reproducible gene-level response remains after platelet-score adjustment

In the full 46-sample longitudinal model with patient fixed effects and platelet-score adjustment, no individual gene reached FDR<0.05 for either the responder-by-week-1 or responder-by-month-1 interaction. This negative result argues against selecting a residual hub gene from nominal associations and indicates that the dominant reproducible whole-blood response is compositional in this cohort.

### Month-1 data support persistence but not a precisely resolved trajectory

Among complete cases, the platelet RNA score remained above baseline at month 1 in responders (n=9; mean change, 0.496; exact sign-flip P=0.00781), whereas no baseline-to-month-1 change was observed in nonresponders (n=5; mean change, 0.002; P=1.00). The between-group difference did not reach significance (P=0.172). Among patients observed at both week 1 and month 1, responder scores decreased modestly (n=8; mean change, −0.192; P=0.0859), while the between-group difference for this interval was borderline (P=0.0544). Thus, complete cases support persistence of platelet-associated RNA recovery but do not define a statistically resolved week-1 peak.

### Composition dominance generalizes across erythroid restoration and B-cell depletion

In GSE186294, the fixed erythroid RNA score increased from Base2 to EPO4 in the 10 paired participants (mean standardized change, 1.187; exact sign-flip P=0.00195). After removal of signature overlap, erythroid-score change explained 35.6% of residual multivariate variation in the 200 most induced genes (10,000-bootstrap 95% CI, 26.9%–66.9%). The median gene-level correlation with erythroid-score change was 0.463 and exceeded 5,000 expression-matched gene sets (empirical P=0.00020); fully nested leave-one-participant-out held-out R2 was positive but modest (0.066).

In GSE112594, the baseline-to-week-26 B-cell score change was lower after rituximab than after placebo (active n=25; placebo n=10; difference, −1.912; 100,000-label permutation P=0.00001). B-cell-score change explained 19.8% of residual multivariate variation in the 200 most depleted genes (bootstrap 95% CI, 9.1%–38.4%). Median gene-level coupling was 0.664, exceeding matched gene sets (empirical P=0.00020), and fully nested held-out R2 was 0.158.

Other lineage scores sometimes changed in the opposite direction, including a T-cell-score decrease after erythropoietin and relative increase after rituximab. These findings are expected under closed-sum compositional substitution and preclude interpreting non-target lineage scores as universally invariant negative controls. The transport result is instead that the intervention-targeted lineage score changed in the biologically expected direction and captured a significant component of held-out transcriptome variation in both settings.

### A composition-only benchmark recapitulates broad treatment-associated differential expression

Across 100 simulated experiments using purified profiles from four donors, the fixed B-cell score closely tracked the known B-cell fraction (median r=0.967; 2.5th–97.5th percentile, 0.957–0.972). Overlapping fraction distributions limited group–score correlation to a median of 0.598 (0.487–0.703), and the donor-adjusted model remained numerically stable (median condition number, 5.83; 5.53–6.53).

Despite the absence of any within-cell transcriptional change by construction, donor-adjusted conventional group testing produced a median of 3,037 FDR-significant genes (2.5th–97.5th percentile, 1,933–4,113). Adding the fixed B-cell score reduced this to a median of 0 genes (0–2), attenuated the median top-200 group effect by 100.7% (98.8%–102.7%), and explained 89.1% of residual top-200 variation (87.1%–90.8%). In the donor-stratified null-label calibration, the median number of FDR-significant genes was 0 (0–1), arguing against systematic false-positive inflation by the analysis workflow. Thus, composition change alone was sufficient to generate a large apparent treatment transcriptome, and an externally specified lineage score recovered the known source without pathological model collinearity.

### A putative immature-platelet signal was not robustly separable from platelet abundance

Across 25-, 50-, and 100-gene definitions, the unadjusted young-minus-mature contrast increased in responders and differed from nonresponders. After residualization against platelet RNA abundance across all samples, responder changes were small and nonsignificant for every signature size (mean changes, −0.046 to 0.100; P=0.404–0.650), and between-group P values ranged from 0.203 to 0.477. Baseline-fitted residualization produced larger estimates for the 50- and 100-gene signatures, but responder sign-flip tests remained nonsignificant (P=0.0762 and P=0.0703), with nominal between-group P values arising only in selected sensitivity specifications. These data are compatible with an early contribution from RNA-rich platelets but do not establish maturation independent of overall platelet-derived RNA abundance.

### Exploratory datasets do not support a general treatment-linked FAO mechanism

Platelet adjustment markedly attenuated apparent whole-blood lipid-pathway effects. A fatty-acid-oxidation difference was detected in one purified T-cell chronicity comparison, but not in a separate ITP-versus-control T-cell dataset or in donor-level bone-marrow CD34-positive pseudobulk. Because none of these datasets measured longitudinal eltrombopag response, they neither replicate nor support a general treatment-linked FAO mechanism (Supplementary Results and Supplementary Figure 4).

## Discussion

This study indicates that platelet-derived RNA is the dominant measured axis within the selected responder-associated whole-blood gene set during eltrombopag treatment, and that the underlying interpretive problem extends beyond ITP. The same externally specified analytical logic detected erythroid restoration after erythropoietin and B-cell depletion after rituximab, while a composition-only benchmark generated thousands of apparently differential genes in the complete absence of within-cell regulation. Together, these analyses move the conclusion from a single-cohort correction to a cross-lineage principle: when therapy changes a blood lineage, lineage-derived RNA can become a major axis of the bulk treatment transcriptome.

The result is biologically coherent rather than merely technical. Eltrombopag is intended to restore platelet production, and platelet RNA is therefore part of the therapeutic phenotype captured by whole blood. The interpretive problem arises when recovery of lineage-derived RNA is labeled as cell-intrinsic molecular reprogramming. Similar bias may occur whenever treatment changes the abundance of a lineage that contributes RNA to a bulk specimen.

The finding was also robust to how the platelet signature was defined. Neither PF4, PPBP, nor any other individual marker was required to recover the patient-level pattern, and a stricter dual-reference signature yielded nearly identical scores. Moreover, removal of every score gene from the summarized responder-induced outcome set preserved the strong correlation and attenuation. These analyses reduce two common sources of circularity in bulk-expression composition studies: dependence on a single canonical marker and direct reuse of predictor genes in the outcome summary.

The quantitative variance analysis extends this audit beyond coefficient attenuation. Platelet-score change accounted for approximately three quarters of residual multivariate variation in the overlap-free responder-induced gene set after response group was included, and the incremental effect remained large when patients were held out and gene selection was repeated within each training fold. These estimates are conditional on the selected responder-induced transcriptome and must not be interpreted as the proportion of all whole-blood transcriptomic variance attributable to platelets. They nevertheless identify platelet-associated RNA as the dominant measured axis within the treatment-responsive signal.

The platelet-maturation analysis further illustrates an identifiability limit of bulk whole-blood RNA sequencing. Immature platelets contain more RNA, making platelet age and platelet RNA abundance biologically and statistically coupled. Although unadjusted young-platelet signatures rose in responders, the effect was not consistently retained after abundance adjustment. Direct immature platelet fraction, reticulated platelet measurements, or purified platelet profiling would be required to distinguish increased platelet mass from a shift in platelet-age structure.

More generally, we used a six-component audit for composition dominance: an externally defined target-lineage score changed in the expected direction; patient-level changes coupled to that score; adjustment produced substantial attenuation and partial R2; expression-matched gene sets failed to reproduce the coupling; the score improved held-out transcriptome reconstruction; and no multiplicity-corrected residual signal remained in the primary cohort. The external cohorts show that effect magnitude is context-dependent: partial and held-out R2 values were smaller outside ITP, as expected when treatment has broader biological effects or the target lineage contributes less RNA. Composition dominance should therefore be quantified rather than assumed to be all-or-none.

The composition-only benchmark supplies the causal element that observational cohorts cannot. Because cellular profiles were held fixed and only mixture proportions changed, every conventional group difference was necessarily compositional. The near-complete collapse of differential expression after lineage-score adjustment demonstrates that a broad signature and its removal can arise without invoking cell-intrinsic reprogramming. The overlapping mixture design is important: it preserved independent variation between group and score, produced low model condition numbers, and avoided the unstable inference that occurs when a composition proxy nearly determines treatment group. The donor-stratified null-label calibration yielded essentially no discoveries, supporting calibration of the analysis workflow under no group–composition association.

Cross-lineage score changes also clarify the limits of simplistic specificity tests. In a closed mixture, expansion of one lineage necessarily reduces the relative contribution of others, so opposite-direction changes in T-cell or other scores may be genuine compositional consequences rather than failed negative controls. Accordingly, the strongest evidence is not invariance of every non-target score, but concordance between the intervention-targeted lineage, the direction of perturbation, broad gene-level coupling, matched-set specificity, and held-out reconstruction.

Our analysis also defines what cannot be concluded. None of the lineage RNA scores is a calibrated cell count or fraction. The simulated benchmark validates behavior against known mixture fractions, but it does not retrospectively supply measured counts to the clinical cohorts. Sample-linked platelet counts therefore remain the most important missing validation. Moreover, absence of FDR-significant residual genes does not prove that eltrombopag has no immune effects; it shows that this small public cohort does not resolve a reproducible gene-level effect after accounting for the dominant platelet-associated signal.

The secondary metabolic audit reached the same interpretive boundary: apparent whole-blood lipid-pathway changes were attenuated after platelet adjustment, and discordant purified-cell datasets could not establish that eltrombopag restores T-cell FAO.

Although patient bootstrap and fully nested leave-one-patient-out analyses support stability, limitations include the small primary cohort, particularly the paired nonresponder group, incomplete longitudinal sampling, retrospective analysis of public data, and lack of an independent longitudinal TPO-receptor-agonist cohort. Three responders received intravenous immunoglobulin 4–8 days before pretreatment RNA collection, creating potential baseline perturbation that cannot be separated with the public data. The transported cohorts test a general principle but are not clinical replications of eltrombopag response. Partial R2 values are conditional on directionally selected gene sets rather than genome-wide variance estimates. The benchmark uses idealized mixtures of purified immune profiles and cannot represent all technical and biological interactions in whole blood. Sample-linked platelet counts, immature platelet fractions, detailed time-varying rescue-therapy data, and purified platelet validation remain unavailable.

In conclusion, platelet-associated RNA dominated the selected responder-associated whole-blood gene set during eltrombopag treatment, and analogous composition-associated signals arose across lineage restoration, lineage depletion, and controlled cell mixtures. Explicit modeling of lineage-derived RNA should be a routine identifiability step in transcriptomic biomarker studies whenever the intervention itself changes blood-cell abundance.

## Data sharing statement

All source datasets are publicly available through the Gene Expression Omnibus under accessions GSE112278, GSE302674, GSE262073, GSE186294, GSE112594, GSE107011, GSE126448, GSE43177, GSE46922, and GSE196676. Analysis code, fixed signature definitions, processed subject-level changes, benchmark replicate results, and figure source data will be deposited at [PUBLIC REPOSITORY URL; VERSIONED RELEASE REQUIRED BEFORE SUBMISSION].

## Research ethics

This study was a secondary analysis of publicly available deidentified datasets and involved no new participant recruitment or specimen collection. Ethics approvals and informed-consent procedures for the original studies are described in their source publications.

Consent for publication: Not applicable.

Acknowledgments: None.

## Authorship

Contribution: Y.Z. conceptualized the study; developed the methodology and software; curated and formally analyzed the data; performed the investigation and visualization; and wrote the original draft. J.P. supervised and administered the project and reviewed and edited the manuscript. Both authors read and approved the final manuscript.

Conflict-of-interest disclosure: The authors declare no competing financial interests.

Funding: This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## References

1. Newman AM, Liu CL, Green MR, et al. Robust enumeration of cell subsets from tissue expression profiles. Nat Methods. 2015;12(5):453-457. doi:10.1038/nmeth.3337.
2. Rowley JW, Oler AJ, Tolley ND, et al. Genome-wide RNA-seq analysis of human and mouse platelet transcriptomes. Blood. 2011;118(14):e101-e111. doi:10.1182/blood-2011-03-339705.
3. Zhang H, Zhang BM, Guo X, et al. Blood transcriptome and clonal T-cell correlates of response and non-response to eltrombopag therapy in a cohort of patients with chronic immune thrombocytopenia. Haematologica. 2020;105(3):e129-e132. doi:10.3324/haematol.2019.226688.
4. Wang G, Kitaoka T, Crawford A, et al. Cross-platform transcriptomic profiling of the response to recombinant human erythropoietin. Sci Rep. 2021;11(1):21705. doi:10.1038/s41598-021-00608-9.
5. Linsley PS, Greenbaum CJ, Rosasco M, Presnell S, Herold KC, Dufort MJ. Elevated T cell levels in peripheral blood predict poor clinical response following rituximab treatment in new-onset type 1 diabetes. Genes Immun. 2019;20(4):293-307. doi:10.1038/s41435-018-0032-1.
6. Dufort MJ, Greenbaum CJ, Speake C, Linsley PS. Cell type-specific immune phenotypes predict loss of insulin secretion in new-onset type 1 diabetes. JCI Insight. 2019;4(4):125556. doi:10.1172/jci.insight.125556.
7. Monaco G, Lee B, Xu W, et al. RNA-Seq signatures normalized by mRNA abundance allow absolute deconvolution of human immune cell types. Cell Rep. 2019;26(6):1627-1640.e7. doi:10.1016/j.celrep.2019.01.041.
8. Hille L, Lenz M, Vlachos A, et al. Ultrastructural, transcriptional, and functional differences between human reticulated and non-reticulated platelets. J Thromb Haemost. 2020;18(8):2034-2046. doi:10.1111/jth.14895.
9. Garshick MS, Drenkova K, Kazatsker F, et al. Platelet activation and a platelet biosignature are associated with cardiovascular risk in patients with controlled psoriasis. Arterioscler Thromb Vasc Biol. 2025;45(11):2086-2096. doi:10.1161/ATVBAHA.125.322574.
10. Banerjee M, Rowley JW, Stubben CJ, et al. Prospective, international, multisite comparison of platelet isolation techniques for genome-wide transcriptomics: communication from the SSC of the ISTH. J Thromb Haemost. 2024;22(10):2922-2934. doi:10.1016/j.jtha.2024.06.017.
11. Jernås M, Nookaew I, Wadenvik H, Olsson B. MicroRNA regulate immunological pathways in T-cells in immune thrombocytopenia (ITP). Blood. 2013;121(11):2095-2098. doi:10.1182/blood-2012-12-471250.
12. Jernås M, Hou Y, Strömberg Célind F, et al. Differences in gene expression and cytokine levels between newly diagnosed and chronic pediatric ITP. Blood. 2013;122(10):1789-1792. doi:10.1182/blood-2013-05-502807.
13. Liu Y, Zuo X, Chen P, et al. Deciphering transcriptome alterations in bone marrow hematopoiesis at single-cell resolution in immune thrombocytopenia. Signal Transduct Target Ther. 2022;7(1):347. doi:10.1038/s41392-022-01167-9.
14. Gillespie M, Jassal B, Stephan R, et al. The Reactome pathway knowledgebase 2022. Nucleic Acids Res. 2022;50(D1):D687-D692. doi:10.1093/nar/gkab1028.

## Figure legends

### Figure 1. Study design and fixed analytical framework

The primary longitudinal cohort (GSE112278) was analyzed using patient-level pairing and fixed effects, with a platelet RNA score defined from two purified-platelet references before outcome testing. The same framework was transported to erythropoietin-associated erythroid recovery (GSE186294) and rituximab-associated B-cell depletion (GSE112594). Purified immune profiles from GSE107011 were used for a composition-only benchmark in which cell-intrinsic expression was fixed by construction.

### Figure 2. Patient-level platelet RNA score trajectories

Patient trajectories from baseline to week 1 and month 1 are shown separately for responders and nonresponders. The paired baseline-to-week-1 change was positive in responders (n=11; mean 0.886 standardized units) and negative in nonresponders (n=4; mean −0.371). Exact responder sign-flip P=0.00195; exact between-group P=0.0381.

### Figure 3. Coupling and attenuation of the responder-associated transcriptome

Each point represents one of the 208 genes most induced at week 1 in responders. The horizontal axis shows correlation between patient-level gene-expression change and platelet-score change; the vertical axis shows attenuation of the response coefficient after platelet-score adjustment. Marginal summaries compare the observed medians with 5,000 expression-matched random gene sets.

### Figure 4. Platelet specificity against nonplatelet blood-lineage controls

Mean paired baseline-to-week-1 score changes are shown for responders and nonresponders for the external platelet signature and five prespecified lineage controls. Error bars indicate standard errors. Only the platelet signature showed a significant responder increase and responder–nonresponder difference.

### Figure 5. Cross-lineage transport and composition-only benchmark

(A) Three intervention–lineage settings: eltrombopag/platelet restoration, erythropoietin/erythroid restoration, and rituximab/B-cell depletion. (B) Intervention-targeted lineage-score changes with patient-level observations. (C) Incremental multivariate partial R2 and fully nested leave-one-patient-out held-out transcriptome R2 after removal of direct signature overlap. Held-out lineage-score change and treatment group were model inputs; this is not clinical-outcome prediction. (D) Composition-only B-cell benchmark across 100 replicates. Donor-adjusted conventional differential expression produced a median 3,037 FDR-significant genes, whereas B-cell-score adjustment produced a median of 0; points show replicates and bars show medians with 2.5th–97.5th percentiles. Scores are RNA-abundance proxies and not measured cell counts.

### Supplementary Figure 1. Robustness to platelet-signature definition

Panel A shows the expression percentile of each locked marker in the two independent purified-platelet references. Panel B shows paired week-1 score changes after removing each marker in turn. Panel C shows the correlation of each leave-one-marker-out score with the locked 12-gene score across all primary-cohort samples.

### Supplementary Figure 2. Leave-one-patient-out stability

Leave-one-patient-out estimates are shown for responder platelet-score change, median top-208 correlation with platelet-score change, and median attenuation after adjustment. The final panel contrasts observed correlation and attenuation statistics with expression-matched permutation null medians.

### Supplementary Figure 3. Exploratory analysis of a putative young, RNA-rich platelet wave

(A) Independent derivation of young- and mature-platelet signatures in four paired donors from GSE126448. (B) Patient-level trajectories of the young-minus-mature score after residualization against the locked platelet RNA-abundance score; thick lines indicate group means. (C) Mean baseline-to-week-1 changes in platelet RNA abundance, the unadjusted young-minus-mature contrast, and the abundance-adjusted age residual. Exact sign-flip and label-permutation tests were used because of the small cohort. The analysis is exploratory and does not establish an abundance-independent platelet-maturation effect.

### Supplementary Figure 4. Orthogonal immunometabolic audit

Locked pathway effects are summarized in purified T cells and donor-level bone-marrow pseudobulk. The datasets address different disease contrasts and do not test longitudinal eltrombopag response.

### Supplementary Figure 5. Primary-cohort variance attribution and operational audit

(A) Multivariate partial R2, median gene-level correlation, and median response-effect attenuation after removing genes directly overlapping the platelet score; bars show patient-bootstrap 95% confidence intervals and the annotation reports fully nested leave-one-patient-out held-out transcriptome R2. (B) Platelet-abundance-adjusted young-platelet residuals across 25-, 50-, and 100-gene definitions. (C) Complete-case month-1 platelet-score changes. (D) Six-component operational audit for composition dominance. The age and month-1 panels are exploratory.

 
