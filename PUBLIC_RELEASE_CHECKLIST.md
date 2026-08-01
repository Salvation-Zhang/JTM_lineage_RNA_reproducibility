# Public release checklist

## Release status

- Release version: `v1.2.0`
- Repository: `https://github.com/Salvation-Zhang/JTM_lineage_RNA_reproducibility`
- Version-specific DOI: `10.5281/zenodo.21722268`
- Concept DOI: `10.5281/zenodo.21632172`
- Analysis scripts, fixed signatures, subject-level results, 100 benchmark replicates, null-label results, figures, and source-data tables are included.
- `CITATION.cff`, `LICENSE`, `requirements.txt`, and checksums are included.

## Reproducibility and metadata checks

1. Confirm that the GitHub repository name, README, `CITATION.cff`, and release notes all identify the JTM manuscript and `v1.2.0`.
2. Confirm that the Zenodo v1.2.0 record uses the JTM title and the renamed GitHub URL while retaining the existing DOIs.
3. Confirm that manuscript, supplementary material, cover letter, and data-sharing statements use the same repository URL and DOI pair.
4. Confirm that no current analysis file describes the leave-one-patient-out reconstruction as fully nested or as clinical-response prediction.
5. Add the recalculated platelet RNA-score specificity-model script, its result table, and `environment_specificity_model.txt`; verify that they reproduce the reported classical OLS and patient-clustered CR1 estimates.

## Author confirmations before submission

- Confirm that all authors approved the final submitted version.
- Confirm Jun Peng's preferred degree/credential if the journal submission system requires one; the letter currently uses the safe name-only form.
- Verify the exact modified IWG response threshold from the original study Supplementary Appendix. Retain the wording that response grouping follows the source study.
