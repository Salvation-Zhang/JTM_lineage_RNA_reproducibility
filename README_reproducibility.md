# Blood Advances reproducibility release

This directory contains the manuscript, supplementary methods/results, source-data tables, figure PDFs, and scripts for:

**Therapy-driven lineage shifts can dominate whole-blood treatment transcriptomes**

## Reproduce

1. Install the packages listed in `reproducibility/requirements.txt`.
2. Place the GEO matrices and NCBI human gene annotation listed below in `reproducibility/data/`.
3. Run `run_cross_lineage_validation.py`.
4. Run `run_mixing_benchmark.py`.
5. Run the figure-generation scripts.

The scripts use fixed random seeds. The benchmark includes donor fixed effects and a donor-stratified null-label calibration.

## Public data accessions

GSE112278, GSE302674, GSE262073, GSE186294, GSE112594, GSE107011, GSE126448, GSE43177, GSE46922, and GSE196676. Gene symbols were mapped using the NCBI Homo sapiens gene-information file.

## Release status

This is a Zenodo/GitHub-ready release. Before submission, publish the `reproducibility/` directory and result tables to a public repository and replace the manuscript data-sharing placeholder with the permanent URL and version/DOI.

## Scope

Lineage scores are RNA-abundance proxies, not measured cell counts. The primary cohort has no sample-linked platelet count table; the computational benchmark validates behavior under known mixtures but does not replace clinical count calibration.
