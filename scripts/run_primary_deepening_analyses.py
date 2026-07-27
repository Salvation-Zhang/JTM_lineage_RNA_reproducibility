from pathlib import Path
import os
import glob
import itertools
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


REPO = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("TISTA_PRIMARY_RAW", REPO / "data" / "primary_raw"))
OUT = REPO / "results" / "primary_deepening"
OUT.mkdir(parents=True, exist_ok=True)

PLATELET = ["PPBP", "PF4", "GP1BA", "GP9", "ITGA2B", "ITGB3", "TUBB1", "TREML1", "RGS18", "SDPR", "SPARC", "CLU"]
LINEAGES = {
    "erythroid": ["ALAS2", "AHSP", "HBB", "HBA1", "HBA2", "CA1", "GYPA", "SLC4A1", "EPB42", "BPGM", "FECH", "KLF1", "ANK1", "TMOD1", "HEMGN"],
    "neutrophil": ["CSF3R", "FCGR3B", "CEACAM8", "S100A8", "S100A9", "MNDA", "FPR1", "CXCR2", "MMP8", "OLFM4"],
    "monocyte": ["CTSS", "FCN1", "VCAN", "CTSD", "LILRB1", "LST1", "S100A10", "FCER1G", "TYMP", "LGALS3"],
    "T_cell": ["CD3D", "CD3E", "CD3G", "TRAC", "LCK", "MAL", "IL7R", "LTB", "CD247", "TRBC1"],
    "B_cell": ["CD19", "MS4A1", "CD79A", "CD79B", "CD22", "CD37", "CD74", "HLA-DRA", "BANK1", "BLK", "FCRL1", "FCRL2", "CD180", "TNFRSF13C", "HVCN1"],
}
IVIG_PATIENTS = {"R7", "R10", "R11"}


def load_two_col_gz(path):
    x = pd.read_csv(path, sep="\t", header=None, index_col=0, compression="gzip")
    x = pd.to_numeric(x.iloc[:, 0], errors="coerce").fillna(0)
    x.name = Path(path).name.split("_")[0]
    return x


def load_primary():
    series, rows = [], []
    for path in sorted(glob.glob(str(RAW / "GSM*.gz"))):
        match = re.search(r"(GSM\d+)_((?:R|NR)\d+)_(pre|1wk|1mon)_", Path(path).name)
        if not match:
            continue
        series.append(load_two_col_gz(path))
        rows.append({
            "gsm": match.group(1),
            "patient": match.group(2),
            "response": "nonresponder" if match.group(2).startswith("NR") else "responder",
            "time": match.group(3),
        })
    meta = pd.DataFrame(rows)
    counts = pd.concat(series, axis=1).fillna(0)[meta.gsm]
    logcpm = np.log2(counts.add(0.5).div(counts.sum(axis=0).add(1), axis=1) * 1e6 + 1)
    return meta, counts, logcpm


def zscore_signature(logcpm, genes):
    genes = [g for g in genes if g in logcpm.index]
    x = logcpm.loc[genes]
    z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), genes


def rank_signature(logcpm, genes):
    """Two-sided centered single-sample rank score, scaled approximately to [-1, 1]."""
    genes = [g for g in genes if g in logcpm.index]
    ranks = logcpm.rank(axis=0, method="average", pct=True)
    return 2 * ranks.loc[genes].mean(axis=0) - 1, genes


def paired_matrix(meta, values):
    x = meta[["gsm", "patient", "response", "time"]].copy()
    x["value"] = x.gsm.map(values)
    wide = x.pivot(index=["patient", "response"], columns="time", values="value").dropna(subset=["pre", "1wk"])
    wide["delta"] = wide["1wk"] - wide["pre"]
    return wide.reset_index()


def expression_delta(meta, logcpm):
    wide = meta[meta.time.isin(["pre", "1wk"])].pivot(index=["patient", "response"], columns="time", values="gsm").dropna()
    patients = [idx[0] for idx in wide.index]
    response = np.asarray([idx[1] for idx in wide.index])
    delta = pd.DataFrame({idx[0]: logcpm[row["1wk"]] - logcpm[row["pre"]] for idx, row in wide.iterrows()})
    return wide, patients, response, delta


def exact_between(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    values = np.r_[a, b]
    observed = abs(a.mean() - b.mean())
    null = []
    for chosen in itertools.combinations(range(len(values)), len(a)):
        mask = np.zeros(len(values), bool)
        mask[list(chosen)] = True
        null.append(abs(values[mask].mean() - values[~mask].mean()))
    extreme = np.sum(np.asarray(null) >= observed - 1e-12)
    return (extreme + 1) / (len(null) + 1)


def fit_gene_models(y, response, platelet_delta, covariates=None):
    response_term = (np.asarray(response) == "responder").astype(float)
    x0 = np.column_stack([np.ones(len(response)), response_term])
    if covariates is not None and covariates.shape[1]:
        x0 = np.column_stack([x0, covariates])
    x1 = np.column_stack([x0, platelet_delta])
    b0 = np.linalg.lstsq(x0, y, rcond=None)[0]
    b1 = np.linalg.lstsq(x1, y, rcond=None)[0]
    e0 = y - x0 @ b0
    e1 = y - x1 @ b1
    response_index = 1
    attenuation = np.divide(
        b0[response_index] - b1[response_index],
        b0[response_index],
        out=np.full(y.shape[1], np.nan),
        where=np.abs(b0[response_index]) > 1e-12,
    )
    partial_r2 = 1 - np.sum(e1 ** 2) / np.sum(e0 ** 2)
    return b0[response_index], b1[response_index], attenuation, partial_r2


def patient_level_correlations(y, platelet_delta):
    yc = y - y.mean(axis=0)
    pc = platelet_delta - platelet_delta.mean()
    denominator = np.sqrt(np.sum(yc ** 2, axis=0) * np.sum(pc ** 2))
    return np.divide(pc @ yc, denominator, out=np.full(y.shape[1], np.nan), where=denominator > 0)


def ivig_exclusion(meta, logcpm):
    zscore, _ = zscore_signature(logcpm, PLATELET)
    paired = paired_matrix(meta, zscore.to_dict())
    rows = []
    for label, keep in [
        ("all_paired_patients", np.ones(len(paired), bool)),
        ("exclude_recent_IVIG_R7_R10_R11", ~paired.patient.isin(IVIG_PATIENTS)),
    ]:
        q = paired.loc[keep]
        responder = q.loc[q.response == "responder", "delta"].to_numpy()
        nonresponder = q.loc[q.response == "nonresponder", "delta"].to_numpy()
        rows.append({
            "analysis": label,
            "responder_n": len(responder),
            "nonresponder_n": len(nonresponder),
            "responder_mean_change": responder.mean(),
            "nonresponder_mean_change": nonresponder.mean(),
            "between_group_difference": responder.mean() - nonresponder.mean(),
            "two_sided_exact_label_permutation_p": exact_between(responder, nonresponder),
        })
    pd.DataFrame(rows).to_csv(OUT / "IVIG_exclusion_platelet_score.tsv", sep="\t", index=False)

    wide, patients, response, delta = expression_delta(meta, logcpm)
    expressed = logcpm.mean(axis=1) >= 1
    genes = logcpm.index[expressed & ~logcpm.index.isin(PLATELET)]
    for label, keep_patients in [
        ("all_paired_patients", patients),
        ("exclude_recent_IVIG_R7_R10_R11", [p for p in patients if p not in IVIG_PATIENTS]),
    ]:
        ix = np.asarray([patients.index(p) for p in keep_patients])
        y = delta.loc[genes, keep_patients].T.to_numpy()
        r = response[ix]
        pdelta = np.asarray([
            zscore[wide.loc[(p, response[patients.index(p)]), "1wk"]] - zscore[wide.loc[(p, response[patients.index(p)]), "pre"]]
            for p in keep_patients
        ])
        b0, b1, attenuation, partial_r2 = fit_gene_models(y, r, pdelta)
        corr = patient_level_correlations(y, pdelta)
        pd.DataFrame([{
            "analysis": label,
            "patient_n": len(keep_patients),
            "gene_n": len(genes),
            "genomewide_multivariate_partial_r2": partial_r2,
            "median_gene_platelet_delta_correlation": np.nanmedian(corr),
            "median_response_coefficient_attenuation": np.nanmedian(attenuation),
            "cor_responder_effect_before_vs_after_adjustment": np.corrcoef(b0, b1)[0, 1],
        }]).to_csv(OUT / f"IVIG_exclusion_genomewide_{label}.tsv", sep="\t", index=False)


def genomewide_audit(meta, logcpm):
    wide, patients, response, delta = expression_delta(meta, logcpm)
    platelet, _ = zscore_signature(logcpm, PLATELET)
    pdelta = np.asarray([platelet[row["1wk"]] - platelet[row["pre"]] for _, row in wide.iterrows()])
    mean_logcpm = logcpm.mean(axis=1)
    delta_sd = delta.std(axis=1, ddof=1)
    eligible = (mean_logcpm >= 1) & (delta_sd > 0) & ~logcpm.index.isin(PLATELET)
    genes = logcpm.index[eligible]
    y = delta.loc[genes, patients].T.to_numpy()
    b0, b1, attenuation, partial_r2 = fit_gene_models(y, response, pdelta)
    y_standardized = y / np.std(y, axis=0, ddof=1, keepdims=True)
    _, _, _, standardized_partial_r2 = fit_gene_models(y_standardized, response, pdelta)
    corr = patient_level_correlations(y, pdelta)
    table = pd.DataFrame({
        "gene": genes,
        "mean_log2_cpm": mean_logcpm.loc[genes].to_numpy(),
        "delta_sd": delta_sd.loc[genes].to_numpy(),
        "response_beta_unadjusted": b0,
        "response_beta_platelet_adjusted": b1,
        "attenuation_fraction": attenuation,
        "correlation_with_platelet_score_change": corr,
    })
    table.to_csv(OUT / "genomewide_selection_independent_gene_audit.tsv.gz", sep="\t", index=False, compression="gzip")
    positive = b0 > 0
    summary = pd.DataFrame([{
        "gene_universe_definition": "mean log2 CPM >=1; nonzero paired-change SD; locked platelet-score genes excluded",
        "gene_n": len(table),
        "positive_responder_effect_gene_n": int(positive.sum()),
        "genomewide_multivariate_partial_r2": partial_r2,
        "gene_standardized_genomewide_partial_r2": standardized_partial_r2,
        "median_correlation_all_genes": np.nanmedian(corr),
        "median_attenuation_all_genes": np.nanmedian(attenuation),
        "median_correlation_positive_effect_genes": np.nanmedian(corr[positive]),
        "median_attenuation_positive_effect_genes": np.nanmedian(attenuation[positive]),
        "fraction_positive_effect_genes_attenuated_gt_50pct": np.nanmean(attenuation[positive] > 0.5),
        "cor_responder_effect_before_vs_after_adjustment": np.corrcoef(b0, b1)[0, 1],
    }])
    summary.to_csv(OUT / "genomewide_selection_independent_summary.tsv", sep="\t", index=False)


def rank_composition_sensitivity(meta, logcpm):
    scores = {"platelet": rank_signature(logcpm, PLATELET)[0]}
    for name, genes in LINEAGES.items():
        scores[name] = rank_signature(logcpm, genes)[0]
    sample_scores = pd.DataFrame(scores)
    sample_scores.index.name = "gsm"
    sample_scores.to_csv(OUT / "rank_based_sample_composition_scores.tsv", sep="\t")

    wide, patients, response, delta = expression_delta(meta, logcpm)
    score_delta = {}
    for name, values in scores.items():
        score_delta[name] = np.asarray([values[row["1wk"]] - values[row["pre"]] for _, row in wide.iterrows()])
    covariates = np.column_stack([score_delta[name] for name in LINEAGES])
    # Standardize paired changes so coefficient scales are comparable and conditioning is explicit.
    covariates = (covariates - covariates.mean(axis=0)) / covariates.std(axis=0, ddof=1)
    platelet_delta = score_delta["platelet"]
    platelet_delta = (platelet_delta - platelet_delta.mean()) / platelet_delta.std(ddof=1)

    mean_logcpm = logcpm.mean(axis=1)
    eligible = (mean_logcpm >= 1) & (delta.std(axis=1, ddof=1) > 0) & ~logcpm.index.isin(PLATELET)
    genes = logcpm.index[eligible]
    y = delta.loc[genes, patients].T.to_numpy()
    b0, b1, attenuation, partial_r2 = fit_gene_models(y, response, platelet_delta, covariates)
    corr = patient_level_correlations(y, platelet_delta)
    paired_platelet = pd.DataFrame({"patient": patients, "response": response, "rank_platelet_delta": score_delta["platelet"]})
    responder = paired_platelet.loc[paired_platelet.response == "responder", "rank_platelet_delta"].to_numpy()
    nonresponder = paired_platelet.loc[paired_platelet.response == "nonresponder", "rank_platelet_delta"].to_numpy()
    pd.DataFrame([{
        "method": "single-sample centered rank marker-abundance score with five prespecified nonplatelet lineage controls",
        "patient_n": len(patients),
        "gene_n": len(genes),
        "responder_rank_platelet_mean_change": responder.mean(),
        "nonresponder_rank_platelet_mean_change": nonresponder.mean(),
        "two_sided_exact_between_group_p": exact_between(responder, nonresponder),
        "genomewide_incremental_partial_r2_after_nonplatelet_composition_adjustment": partial_r2,
        "median_gene_correlation_with_rank_platelet_change": np.nanmedian(corr),
        "median_response_effect_attenuation": np.nanmedian(attenuation),
        "design_condition_number": np.linalg.cond(np.column_stack([np.ones(len(response)), (response == "responder").astype(float), covariates, platelet_delta])),
    }]).to_csv(OUT / "rank_based_composition_sensitivity_summary.tsv", sep="\t", index=False)
    paired_platelet.to_csv(OUT / "rank_based_platelet_patient_changes.tsv", sep="\t", index=False)


def main():
    meta, _, logcpm = load_primary()
    ivig_exclusion(meta, logcpm)
    genomewide_audit(meta, logcpm)
    rank_composition_sensitivity(meta, logcpm)
    make_figure()
    print(f"Results written to {OUT}")


def make_figure():
    ivig = pd.read_csv(OUT / "IVIG_exclusion_platelet_score.tsv", sep="\t")
    genome = pd.read_csv(OUT / "genomewide_selection_independent_summary.tsv", sep="\t").iloc[0]
    rank = pd.read_csv(OUT / "rank_based_composition_sensitivity_summary.tsv", sep="\t").iloc[0]
    changes = pd.read_csv(OUT / "rank_based_platelet_patient_changes.tsv", sep="\t")

    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 600,
    })
    colors = {"responder": "#2166AC", "nonresponder": "#B2182B"}
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))

    ax = axes[0]
    labels = ["All paired\npatients", "Exclude R7,\nR10, R11"]
    x = np.arange(2)
    width = 0.34
    ax.bar(x - width / 2, ivig.responder_mean_change, width, color=colors["responder"], label="Responder")
    ax.bar(x + width / 2, ivig.nonresponder_mean_change, width, color=colors["nonresponder"], label="Nonresponder")
    for i, row in ivig.iterrows():
        height = max(row.responder_mean_change, row.nonresponder_mean_change)
        ax.text(i, height + 0.07, f"P={row.two_sided_exact_label_permutation_p:.3f}", ha="center", fontsize=8)
    ax.axhline(0, color="#666666", lw=0.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Mean platelet RNA-score change")
    ax.set_title("Recent-IVIG exclusion")
    ax.set_ylim(-0.45, 1.28)
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2)

    ax = axes[1]
    metrics = [
        genome.genomewide_multivariate_partial_r2,
        genome.fraction_positive_effect_genes_attenuated_gt_50pct,
        genome.cor_responder_effect_before_vs_after_adjustment,
    ]
    names = ["Genome-wide\npartial $R^2$", "Positive genes\n>50% attenuated", "Effect-vector\ncorrelation"]
    ax.bar(np.arange(3), metrics, color=["#4D9221", "#C2A5CF", "#9970AB"])
    ax.set_ylim(0, 1)
    ax.set_xticks(np.arange(3), names)
    ax.set_ylabel("Estimate")
    ax.set_title("Selection-independent gene universe")
    for i, value in enumerate(metrics):
        ax.text(i, value + 0.035, f"{value:.2f}", ha="center", fontsize=8)

    ax = axes[2]
    rng = np.random.default_rng(20260727)
    for j, group in enumerate(["responder", "nonresponder"]):
        values = changes.loc[changes.response == group, "rank_platelet_delta"].to_numpy()
        jitter = rng.uniform(-0.08, 0.08, len(values))
        ax.scatter(np.full(len(values), j) + jitter, values, s=31, color=colors[group], alpha=0.9)
        ax.hlines(values.mean(), j - 0.22, j + 0.22, color="black", lw=1.5)
    ax.axhline(0, color="#666666", lw=0.8)
    ax.set_xticks([0, 1], ["Responder", "Nonresponder"])
    ax.set_ylabel("Rank-based platelet-score change")
    ax.set_title("Alternative score construction")
    ax.text(
        0.5,
        0.97,
        f"Exact P={rank.two_sided_exact_between_group_p:.3f}\nAdjusted partial $R^2$={rank.genomewide_incremental_partial_r2_after_nonplatelet_composition_adjustment:.2f}",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=8,
    )

    for i, ax in enumerate(axes):
        ax.text(-0.14, 1.05, chr(65 + i), transform=ax.transAxes, fontsize=12, fontweight="bold")
    fig.tight_layout(w_pad=2.1)
    for extension in ["pdf", "png", "svg"]:
        fig.savefig(OUT / f"Supplementary_Figure_primary_deepening.{extension}", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
