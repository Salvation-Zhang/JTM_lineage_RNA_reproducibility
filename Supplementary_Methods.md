# Supplementary Methods

## Reproducibility and software

All analyses were run with Python 3 and the scripts in `reproducibility/scripts/`. Random procedures used fixed seeds recorded in the scripts. Source expression matrices and gene-annotation files were not redistributed in this release when they were downloaded directly from GEO or NCBI; accession numbers and retrieval metadata are listed in `README_reproducibility.md`.

## Primary cohort preprocessing

GSE112278 3'-end RNA-sequencing counts were converted to log2 counts per million after a 0.5 pseudocount and library-size normalization. The platelet score was the mean of gene-wise standardized expression for PPBP, PF4, GP1BA, GP9, ITGA2B, ITGB3, TUBB1, TREML1, RGS18, SDPR/CAVIN2, SPARC, and CLU. This score is an RNA-abundance proxy and was not calibrated to platelet counts.

## Exact tests and resampling

Paired score changes used exact sign-flip tests. Responder-versus-nonresponder contrasts used exact label permutations when feasible. The larger rituximab contrast used 100,000 seeded label permutations. Patient-bootstrap intervals used 10,000 stratified resamples. Matched-expression null sets used 5,000 seeded resamples.

## Gene-level and multivariate analyses

The responder-associated set consisted of the 208 genes with the largest mean week-1 induction in responders. Because this selection and the descriptive summaries use the same cohort-level contrast, in-sample coupling, attenuation, and partial R2 are conditional descriptive estimates. Direct overlap with the platelet score was removed before multivariate analyses. Fully nested leave-one-patient-out folds repeated gene selection and overlap removal within the training set. The held-out patient's observed lineage score was used as an explanatory input; the resulting R2 is transcriptome reconstruction, not clinical-outcome prediction.

## Cross-lineage validation

Erythroid and B-cell signatures were fixed from canonical markers without selecting genes on validation treatment effects. The same overlap-free partial-R2 and held-out reconstruction logic was applied to GSE186294 and GSE112594. Opposite changes in non-target scores were interpreted as closed-sum compositional substitution rather than universal negative-control failure.

## Composition-only benchmark

For each of four donors in GSE107011, the mean purified B-cell profile was mixed with the donor-matched mean non-B immune-cell profile. B-cell fractions were sampled from overlapping intervals (3%–14% and 8%–19%), with 12 samples per group per donor, variable library size, and multinomial sampling. Donor fixed effects were included in both group models. A donor-stratified null-label calibration independently permuted group labels within donor in every replicate.

## Exploratory analyses

Month-1 analyses were complete-case only. Young-platelet signatures were derived independently in paired purified platelet samples and residualized against the fixed platelet RNA score. Immunometabolic analyses were secondary, non-replicative audits and are not evidence of longitudinal eltrombopag mechanism.
