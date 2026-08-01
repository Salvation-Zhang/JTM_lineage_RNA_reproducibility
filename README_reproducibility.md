# JTM lineage-RNA reproducibility release

This repository contains the analysis scripts, supplementary methods and results, source-data tables, and figure-generation resources for the manuscript:

**Therapy-associated shifts in lineage-derived RNA can shape apparent whole-blood treatment signatures**

The manuscript is intended for submission to the *Journal of Translational Medicine*.

## Reproduce

1. Install the packages listed in `requirements.txt`. The primary-cohort extension analyses used the exact versions recorded in `environment_primary_deepening.txt`.
2. Place the source files under `data/` using the structure documented in `SOURCE_DATA_INDEX.md`, or set `TISTA_PRIMARY_RAW`, `TISTA_GSE126448_DATA`, `TISTA_CROSS_LINEAGE_DATA`, and `TISTA_PROCESSED_DATA` to the corresponding source-data directories.
3. Run `python scripts/run_cross_lineage_validation.py` and `python scripts/run_mixing_benchmark.py`.
4. Run `python scripts/run_depth_upgrade.py`, `python scripts/run_young_platelet_wave.py`, and `python scripts/run_primary_deepening_analyses.py`.
5. Run the figure-generation scripts as required.

Randomized procedures use fixed seeds specified in the corresponding scripts. The composition-only benchmark includes donor fixed effects and donor-stratified null-label calibration.

## Public data accessions

The analyses use publicly available data from GSE112278, GSE302674, GSE262073, GSE186294, GSE112594, GSE107011, GSE126448, GSE43177, GSE46922, and GSE196676. Gene symbols were mapped using the NCBI *Homo sapiens* gene-information file where required.

## Release and citation

The public code repository is:

https://github.com/Salvation-Zhang/JTM_lineage_RNA_reproducibility

The current archived release is `v1.2.0`.

- Zenodo concept DOI: https://doi.org/10.5281/zenodo.21632172
- Version-specific DOI for `v1.2.0`: https://doi.org/10.5281/zenodo.21722268

The concept DOI resolves to the complete version history. Use the version-specific DOI when citing the exact `v1.2.0` analysis release.

## Scope and interpretation

Lineage scores are RNA-abundance proxies, not measured cell counts or cell fractions. Sample-linked platelet counts were unavailable for the primary cohort. The controlled composition-only benchmark evaluates the analysis under known mixtures but does not replace clinical calibration using measured cell counts, platelet fractions, platelet age structure, or RNA content per platelet.

The cross-lineage analyses evaluate the behavior of the analytical approach in other treatment-lineage settings. They are not clinical replications of the primary eltrombopag-associated finding. Leave-one-patient-out analyses repeated gene selection and model fitting within each training fold, but marker-gene standardization used full-cohort parameters; the reconstruction analysis was therefore not fully nested and should not be interpreted as treatment-response prediction.
