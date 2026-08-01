# Supplementary Methods

## Reproducibility and software

The original analysis release was run in the environment pinned in `requirements.txt`. The primary-cohort deepening analyses were run with Python 3.13.14, NumPy 2.4.6, pandas 3.0.3, SciPy 1.18.0, statsmodels 0.14.6, and Matplotlib 3.11.1; this environment is recorded separately in `environment_primary_deepening.txt`. Scripts are available at https://github.com/Salvation-Zhang/JTM_lineage_RNA_reproducibility, release v1.2.0, and archived in Zenodo under concept DOI 10.5281/zenodo.21632172 and version-specific DOI 10.5281/zenodo.21722268. Random procedures used fixed seeds recorded in the scripts. Source expression matrices and gene-annotation files were not redistributed when downloaded directly from GEO or NCBI; accession numbers and retrieval metadata are listed in `README_reproducibility.md`.

## Primary cohort preprocessing

The source study's modified International Working Group response definition was retained: complete response, platelet count >100 × 10⁹/L with no bleeding; response, platelet count 30–100 × 10⁹/L with at least a twofold increase from baseline and no bleeding; and nonresponse, platelet count <30 × 10⁹/L, failure to double baseline within 90 days, or bleeding. Response labels were not re-estimated from transcriptomic data.

GSE112278 3′-end RNA-sequencing counts were transformed as log₂[((Cᵍᵢ + 0.5)/(Lᵢ + 1)) × 10⁶ + 1], where Cᵍᵢ denotes the raw count for gene g in sample i and Lᵢ denotes the total raw library size for sample i. The platelet score was the mean of gene-wise standardized expression for PPBP, PF4, GP1BA, GP9, ITGA2B, ITGB3, TUBB1, TREML1, RGS18, SDPR/CAVIN2, SPARC, and CLU. This score is an RNA-abundance proxy and was not calibrated to platelet counts or fractions.

## Exact tests and resampling

Paired score changes used exact sign-flip tests. Responder-versus-nonresponder contrasts used exact label permutations when feasible. The larger rituximab contrast used 100,000 seeded Monte Carlo label permutations. Monte Carlo P values were calculated as (b + 1)/(B + 1). Patient-bootstrap intervals used 10,000 stratified resamples. Matched-expression null sets used 5,000 seeded resamples.

## Gene-level and multivariate analyses

The responder-associated set consisted of the 208 genes with the largest mean week-1 induction in responders. Because this selection and the descriptive summaries use the same cohort-level contrast, in-sample coupling, attenuation, and partial R² are conditional descriptive estimates. Direct overlap with the platelet score was removed before multivariate analyses. Leave-one-patient-out folds repeated gene selection, overlap removal, and model fitting within the training set. Marker-gene standardization used full-cohort parameters and was not repeated within folds; consequently, the reconstruction analysis was only partially nested. The held-out patient's observed lineage score was used as an explanatory input; the resulting R² quantified transcriptome reconstruction and was not a clinical-outcome prediction.

## Cross-lineage validation

Erythroid and B-cell signatures were fixed from canonical markers without selecting genes on validation treatment effects. The same overlap-free partial R² and held-out reconstruction logic was applied to GSE186294 and GSE112594. These cohorts support transport of the framework to other treatment-lineage settings but are not clinical replications of the primary eltrombopag finding. Opposite changes in non-target scores were interpreted as closed-sum compositional substitution rather than universal negative-control failure.

## Platelet RNA-score specificity model

The platelet RNA score was modeled across all 46 primary-cohort samples using patient fixed effects, week-1 and month-1 indicators, and responder-by-time interactions. The adjusted model additionally included standardized log library size and standardized erythroid, neutrophil, monocyte, T-cell, and B-cell RNA scores. Classical ordinary-least-squares inference and a sensitivity analysis using patient-clustered CR1 covariance with finite-sample correction were reported. Clustered tests used a t reference with 16 degrees of freedom, corresponding to 17 patient clusters. Condition numbers were calculated directly from each fitted design matrix. The executable analysis, complete sample-level input score table, numerical output, runtime versions, input path, official GEO URL, and source-archive SHA-256 are provided in `scripts/run_specificity_model.py`, `results/primary/specificity_model_sample_scores.tsv`, `results/primary/specificity_model_results.tsv`, and `environment_specificity_model.txt`.

## Composition-only benchmark

For each of four donors in GSE107011, the mean purified B-cell profile was mixed with the donor-matched mean non-B immune-cell profile. B-cell fractions were sampled from overlapping intervals (3%–14% and 8%–19%), with 12 samples per group per donor, variable library size, and multinomial sampling. Donor fixed effects were included in both group models. A donor-stratified null-label calibration independently permuted group labels within donor in every replicate.

## Exploratory analyses

Month-1 analyses were complete-case only. Young-platelet signatures were derived independently in paired purified platelet samples and residualized against the fixed platelet RNA score. Immunometabolic analyses were secondary, non-replicative audits and are not evidence of a longitudinal eltrombopag mechanism.

## Recent-IVIG exclusion and selection-independent analysis

The paired primary-cohort analysis was repeated after exclusion of R7, R10, and R11, the three responders reported to have received intravenous immunoglobulin 4–8 days before pretreatment RNA collection. Effect direction and magnitude were emphasized because the exclusion analysis retained eight responders and four nonresponders. A separate selection-independent analysis used every gene with mean log2 counts per million of at least 1 and nonzero variation in paired baseline-to-week-1 change. Locked platelet-score genes were excluded. Incremental multivariate partial R2 was estimated across the complete patient-by-gene change matrix without selecting genes on the responder contrast. Because pooled residual sums of squares weight genes according to their change variance, the model was repeated after dividing each gene's paired-change vector by its sample standard deviation.

## Rank-based composition sensitivity

Centered single-sample rank-based marker-abundance scores were calculated as twice the mean within-sample percentile rank of detected markers minus 1. Scores were calculated for platelets, erythroid cells, neutrophils, monocytes, T cells, and B cells. The rank-based platelet score used the same locked markers as the primary z-score and therefore represented an alternative score construction rather than an independent biological signature. Genome-wide models included response group and standardized paired changes in the five nonplatelet scores before addition of standardized rank-based platelet-score change. The condition number of the full design matrix was recorded. These scores were treated as transcriptomic marker-abundance proxies rather than calibrated cell fractions.
