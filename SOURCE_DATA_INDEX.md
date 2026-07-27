# Source-data index

## Local source-data layout

Downloaded GEO/NCBI source files are not redistributed. By default, scripts look for:

- `data/primary_raw/`: GSE112278 per-sample count files.
- `data/GSE126448/`: paired immature/high-RNA and mature/low-RNA platelet files.
- `data/cross_lineage/`: GSE107011, GSE112594, GSE186294, and the NCBI human gene-information file.
- `data/processed/`: processed primary-cohort matrices used by legacy robustness and figure scripts.

Equivalent external directories can be supplied with `TISTA_PRIMARY_RAW`, `TISTA_GSE126448_DATA`, `TISTA_CROSS_LINEAGE_DATA`, and `TISTA_PROCESSED_DATA`.

## Primary ITP analyses

- `results/primary/complete_case_month1.tsv`: complete-case month-1 sensitivity results.
- `results/primary/composition_variance_attribution.tsv`: conditional variance attribution.
- `results/primary/patient_bootstrap_and_variance_attribution.tsv`: patient bootstrap and held-out reconstruction summary.
- `results/primary/signature_definition_robustness.tsv`: platelet-signature sensitivity analyses.
- `results/primary/external_reference_marker_transportability.tsv`: external marker transportability.
- `results/primary/mathematical_overlap_audit.tsv`: overlap audit.
- `results/primary/young_platelet_signature_size_sensitivity.tsv`: young-platelet signature-size analysis.

## Primary-cohort deepening analyses

- `results/primary_deepening/IVIG_exclusion_platelet_score.tsv`: paired platelet RNA-score sensitivity after exclusion of R7, R10, and R11.
- `results/primary_deepening/IVIG_exclusion_genomewide_all_paired_patients.tsv`: genome-wide summary in all paired patients.
- `results/primary_deepening/IVIG_exclusion_genomewide_exclude_recent_IVIG_R7_R10_R11.tsv`: genome-wide summary after recent-IVIG exclusion.
- `results/primary_deepening/genomewide_selection_independent_summary.tsv`: selection-independent expressed-gene-universe summary.
- `results/primary_deepening/genomewide_selection_independent_gene_audit.tsv.gz`: complete gene-level selection-independent audit.
- `results/primary_deepening/rank_based_sample_composition_scores.tsv`: sample-level centered rank-based marker-abundance scores.
- `results/primary_deepening/rank_based_platelet_patient_changes.tsv`: paired rank-based platelet-score changes.
- `results/primary_deepening/rank_based_composition_sensitivity_summary.tsv`: rank-based platelet analysis with five nonplatelet lineage controls.
- `results/primary_deepening/Supplementary_Figure_primary_deepening.pdf`: three-panel sensitivity figure.

## Cross-lineage validation and computational benchmark

- `results/cross_lineage/GSE112594_subject_changes.tsv`: subject-level erythropoietin cohort changes.
- `results/cross_lineage/GSE186294_subject_changes.tsv`: subject-level rituximab cohort changes.
- `results/cross_lineage/cross_lineage_composition_dominance_summary.tsv`: cross-lineage results.
- `results/cross_lineage/cross_lineage_score_specificity.tsv`: lineage specificity.
- `results/cross_lineage/B_cell_composition_mixing_benchmark_replicates.tsv`: all 100 benchmark replicates, including donor-stratified null-label counts.
- `results/cross_lineage/B_cell_composition_mixing_benchmark_summary.tsv`: benchmark summary.
- `results/cross_lineage/B_cell_composition_mixing_benchmark_example_samples.tsv`: example mixture sample metadata.

## Young/RNA-rich platelet analyses

- `results/young_platelet/GSE112278_platelet_age_scores.tsv`: sample-level scores in the primary cohort.
- `results/young_platelet/GSE126448_paired_platelet_age_gene_statistics.tsv`: paired purified-platelet statistics.
- `results/young_platelet/platelet_age_longitudinal_tests.tsv`: longitudinal tests.
- `results/young_platelet/platelet_age_signature_audit.tsv`: signature audit.

## Fixed signatures

- `signatures/locked_B_cell_signature.tsv`
- `signatures/locked_erythroid_signature.tsv`
- `signatures/locked_platelet_age_signatures.tsv`
