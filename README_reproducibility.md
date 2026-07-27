# Blood Advances reproducibility release

This directory contains the manuscript, supplementary methods/results, source-data tables, figure PDFs, and scripts for:

**Therapy-driven lineage shifts can shape apparent whole-blood treatment signatures**

## Reproduce

1. Install the packages listed in `requirements.txt`. The added primary-deepening analyses used the exact versions in `environment_primary_deepening.txt`.
2. Place the source files under `data/` using the structure documented in `SOURCE_DATA_INDEX.md`, or set `TISTA_PRIMARY_RAW`, `TISTA_GSE126448_DATA`, `TISTA_CROSS_LINEAGE_DATA`, and `TISTA_PROCESSED_DATA` to the corresponding source-data directories.
3. Run `python scripts/run_cross_lineage_validation.py` and `python scripts/run_mixing_benchmark.py`.
4. Run `python scripts/run_depth_upgrade.py`, `python scripts/run_young_platelet_wave.py`, and `python scripts/run_primary_deepening_analyses.py`.
5. Run the figure-generation scripts as required.

The scripts use fixed random seeds. The benchmark includes donor fixed effects and a donor-stratified null-label calibration.

## Public data accessions

GSE112278, GSE302674, GSE262073, GSE186294, GSE112594, GSE107011, GSE126448, GSE43177, GSE46922, and GSE196676. Gene symbols were mapped using the NCBI Homo sapiens gene-information file.

## Release status

The public code repository is https://github.com/Salvation-Zhang/TISTA_Blood_Advances_reproducibility. This directory prepares release `v1.1.0`. The concept DOI for all versions is https://doi.org/10.5281/zenodo.21557317. The version-specific DOI will be added to the submission files after Zenodo archives the GitHub release.

## Scope

Lineage scores are RNA-abundance proxies, not measured cell counts. The primary cohort has no sample-linked platelet count table; the computational benchmark validates behavior under known mixtures but does not replace clinical count calibration.
