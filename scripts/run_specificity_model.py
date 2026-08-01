from pathlib import Path
import glob
import hashlib
import os
import platform
import re
import sys

import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm


REPO = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("TISTA_PRIMARY_RAW", REPO / "data" / "primary_raw"))
OUT = REPO / "results" / "primary"
OUT.mkdir(parents=True, exist_ok=True)

PLATELET = ["PPBP", "PF4", "GP1BA", "GP9", "ITGA2B", "ITGB3", "TUBB1", "TREML1", "RGS18", "SDPR", "SPARC", "CLU"]
LINEAGES = {
    "erythroid": ["ALAS2", "AHSP", "HBB", "HBA1", "HBA2", "CA1", "GYPA", "SLC4A1", "EPB42", "BPGM", "FECH", "KLF1", "ANK1", "TMOD1", "HEMGN"],
    "neutrophil": ["CSF3R", "FCGR3B", "CEACAM8", "S100A8", "S100A9", "MNDA", "FPR1", "CXCR2", "MMP8", "OLFM4"],
    "monocyte": ["CTSS", "FCN1", "VCAN", "CTSD", "LILRB1", "LST1", "S100A10", "FCER1G", "TYMP", "LGALS3"],
    "T_cell": ["CD3D", "CD3E", "CD3G", "TRAC", "LCK", "MAL", "IL7R", "LTB", "CD247", "TRBC1"],
    "B_cell": ["CD19", "MS4A1", "CD79A", "CD79B", "CD22", "CD37", "CD74", "HLA-DRA", "BANK1", "BLK", "FCRL1", "FCRL2", "CD180", "TNFRSF13C", "HVCN1"],
}


def load_primary():
    series, rows = [], []
    for path in sorted(glob.glob(str(RAW / "GSM*.gz"))):
        match = re.search(r"(GSM\d+)_((?:R|NR)\d+)_(pre|1wk|1mon)_", Path(path).name)
        if not match:
            continue
        values = pd.read_csv(path, sep="\t", header=None, index_col=0, compression="gzip")
        values = pd.to_numeric(values.iloc[:, 0], errors="coerce").fillna(0)
        values.name = match.group(1)
        series.append(values)
        rows.append({
            "gsm": match.group(1),
            "patient": match.group(2),
            "response": "nonresponder" if match.group(2).startswith("NR") else "responder",
            "time": match.group(3),
            "source_file": Path(path).name,
        })
    if len(series) != 46:
        raise RuntimeError(f"Expected 46 GSE112278 samples, found {len(series)} in {RAW}")
    metadata = pd.DataFrame(rows)
    counts = pd.concat(series, axis=1).fillna(0)[metadata.gsm]
    library_size = counts.sum(axis=0)
    logcpm = np.log2(counts.add(0.5).div(library_size.add(1), axis=1) * 1e6 + 1)
    return metadata, logcpm, library_size


def signature_score(logcpm, genes):
    detected = [gene for gene in genes if gene in logcpm.index]
    expression = logcpm.loc[detected]
    standardized = expression.sub(expression.mean(axis=1), axis=0).div(
        expression.std(axis=1, ddof=1).replace(0, np.nan), axis=0
    )
    return standardized.mean(axis=0), detected


def design_matrix(data, adjusted):
    patient = pd.get_dummies(data["patient"], prefix="patient", drop_first=True, dtype=float)
    week1 = (data["time"] == "1wk").astype(float)
    month1 = (data["time"] == "1mon").astype(float)
    responder = (data["response"] == "responder").astype(float)
    columns = {
        "intercept": np.ones(len(data)),
        **{name: patient[name].to_numpy() for name in patient.columns},
        "week1": week1,
        "month1": month1,
        "responder_x_week1": responder * week1,
        "responder_x_month1": responder * month1,
    }
    if adjusted:
        log_library_size = np.log(data["library_size"])
        columns["log_library_size"] = (log_library_size - log_library_size.mean()) / log_library_size.std(ddof=1)
        for name in LINEAGES:
            score = data[f"{name}_score"]
            columns[f"{name}_score"] = (score - score.mean()) / score.std(ddof=1)
    return pd.DataFrame(columns, index=data.index, dtype=float)


def fit_model(data, adjusted):
    x = design_matrix(data, adjusted)
    fit = sm.OLS(data["platelet_score"], x).fit()
    clustered = fit.get_robustcov_results(
        cov_type="cluster",
        groups=data["patient"],
        use_correction=True,
        df_correction=True,
        use_t=True,
    )
    term = x.columns.get_loc("responder_x_week1")
    return {
        "model": "adjusted" if adjusted else "base",
        "n_samples": len(data),
        "n_patients": data["patient"].nunique(),
        "design_columns": x.shape[1],
        "residual_df_classical": fit.df_resid,
        "responder_x_week1_beta": fit.params.iloc[term],
        "classical_ols_se": fit.bse.iloc[term],
        "classical_ols_t": fit.tvalues.iloc[term],
        "classical_ols_p": fit.pvalues.iloc[term],
        "patient_clustered_cr1_se": clustered.bse[term],
        "patient_clustered_cr1_t": clustered.tvalues[term],
        "patient_clustered_cr1_p": clustered.pvalues[term],
        "patient_clustered_reference_df": data["patient"].nunique() - 1,
        "condition_number": np.linalg.cond(x.to_numpy()),
    }


def main():
    metadata, logcpm, library_size = load_primary()
    scores = {}
    detected = {}
    scores["platelet"], detected["platelet"] = signature_score(logcpm, PLATELET)
    for name, genes in LINEAGES.items():
        scores[name], detected[name] = signature_score(logcpm, genes)

    data = metadata.copy()
    data["library_size"] = data["gsm"].map(library_size)
    data["platelet_score"] = data["gsm"].map(scores["platelet"])
    for name in LINEAGES:
        data[f"{name}_score"] = data["gsm"].map(scores[name])
    data.to_csv(OUT / "specificity_model_sample_scores.tsv", sep="\t", index=False)

    results = pd.DataFrame([fit_model(data, False), fit_model(data, True)])
    results.to_csv(OUT / "specificity_model_results.tsv", sep="\t", index=False)

    env_lines = [
        f"python={platform.python_version()}",
        f"python_executable={sys.executable}",
        f"numpy={np.__version__}",
        f"pandas={pd.__version__}",
        f"scipy={scipy.__version__}",
        f"statsmodels={statsmodels.__version__}",
        f"input_directory={RAW.resolve()}",
        "input_accession=GSE112278",
        "input_url=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112278/suppl/GSE112278_RAW.tar",
        "input_archive_sha256=" + (hashlib.sha256((REPO / "data" / "GSE112278_RAW.tar").read_bytes()).hexdigest() if (REPO / "data" / "GSE112278_RAW.tar").exists() else "archive_not_retained"),
        f"input_sample_files={len(data)}",
        "outcome=12-gene platelet RNA score computed from gene-wise standardized log2 CPM",
        "base_model=platelet_score ~ patient fixed effects + week1 + month1 + responder:week1 + responder:month1",
        "adjusted_model=base model + standardized log library size + standardized erythroid, neutrophil, monocyte, T-cell, and B-cell RNA scores",
        "classical_inference=ordinary least squares covariance and residual degrees of freedom",
        "clustered_inference=patient-clustered CR1 covariance with finite-sample correction and t reference on 16 cluster degrees of freedom",
        "condition_number=2-norm condition number of the fitted design matrix",
        "output_results=results/primary/specificity_model_results.tsv",
        "output_sample_scores=results/primary/specificity_model_sample_scores.tsv",
    ]
    (REPO / "environment_specificity_model.txt").write_text("\n".join(env_lines) + "\n", encoding="utf-8")
    print(results.to_string(index=False))
    print("Detected signature genes:", detected)


if __name__ == "__main__":
    main()
