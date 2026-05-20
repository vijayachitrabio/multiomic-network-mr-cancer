#!/usr/bin/env python3

"""CPTAC-BRCA proteomics and immune-deconvolution validation."""

from pathlib import Path
import warnings

import cptac
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


PROJECT_DIR = Path("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project")
OUT_DIR = PROJECT_DIR / "results" / "cptac"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TARGETS = ["EFNA1", "TNFRSF6B", "ATRAID", "ITIH3", "IL34", "FGF5", "APOE"]
IMMUNE_PRIORITY = {
    "xcell": [
        "CD8+ T-cells",
        "Macrophages",
        "Macrophages M1",
        "Macrophages M2",
        "Dendritic cells",
        "B-cells",
        "Tregs",
        "ImmuneScore",
        "StromaScore",
        "MicroenvironmentScore",
    ],
    "cibersort": [
        "T cells CD8",
        "Macrophages M0",
        "Macrophages M1",
        "Macrophages M2",
        "Dendritic cells activated",
        "B cells naive",
        "B cells memory",
        "T cells regulatory (Tregs)",
    ],
}


def flatten_proteomics_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if isinstance(out.columns, pd.MultiIndex):
        out.columns = [str(col[0]) for col in out.columns]
    else:
        out.columns = [str(col) for col in out.columns]
    return out


def fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    ok = np.isfinite(pvals)
    out = np.full(len(pvals), np.nan)
    if ok.any():
        out[ok] = multipletests(pvals[ok], method="fdr_bh")[1]
    return out


def correlate(prot: pd.DataFrame, immune: pd.DataFrame, source: str) -> pd.DataFrame:
    common = prot.index.intersection(immune.index)
    rows = []
    prot2 = prot.loc[common]
    immune2 = immune.loc[common]
    priority = [c for c in IMMUNE_PRIORITY[source] if c in immune2.columns]
    if not priority:
        priority = list(immune2.columns)

    for gene in prot2.columns:
        for feature in priority:
            x = prot2[gene]
            y = immune2[feature]
            ok = x.notna() & y.notna()
            if ok.sum() < 20:
                continue
            rho, p = spearmanr(x[ok], y[ok])
            rows.append(
                {
                    "gene": gene,
                    "immune_source": source,
                    "immune_feature": feature,
                    "n": int(ok.sum()),
                    "rho": rho,
                    "p": p,
                }
            )
    res = pd.DataFrame(rows)
    if not res.empty:
        res["fdr"] = fdr(res["p"])
    return res


def main():
    warnings.filterwarnings("ignore")
    brca = cptac.Brca(no_internet=True)

    proteomics = flatten_proteomics_columns(brca.get_proteomics(source="umich"))
    present = [gene for gene in TARGETS if gene in proteomics.columns]
    missing = [gene for gene in TARGETS if gene not in proteomics.columns]
    target_prot = proteomics[present].copy()
    target_prot.to_csv(OUT_DIR / "cptac_brca_target_protein_abundance_umich.csv")

    coverage = pd.DataFrame(
        {
            "gene": TARGETS,
            "available_in_umich_proteomics": [gene in present for gene in TARGETS],
        }
    )
    coverage.to_csv(OUT_DIR / "cptac_brca_target_protein_coverage.csv", index=False)

    summary = target_prot.agg(["count", "mean", "std", "median"]).T.reset_index()
    summary = summary.rename(columns={"index": "gene"})
    summary.to_csv(OUT_DIR / "cptac_brca_target_protein_summary.csv", index=False)

    xcell = brca.get_xcell(source="washu")
    cibersort = brca.get_cibersort(source="washu")
    xcell.to_csv(OUT_DIR / "cptac_brca_xcell_scores.csv")
    cibersort.to_csv(OUT_DIR / "cptac_brca_cibersort_scores.csv")

    cors = pd.concat(
        [
            correlate(target_prot, xcell, "xcell"),
            correlate(target_prot, cibersort, "cibersort"),
        ],
        ignore_index=True,
    )
    if not cors.empty:
        cors = cors.sort_values(["fdr", "p", "gene"])
    cors.to_csv(OUT_DIR / "cptac_brca_target_protein_immune_correlations.csv", index=False)

    print("Present proteins:", ", ".join(present))
    print("Missing proteins:", ", ".join(missing))
    print("Wrote outputs to", OUT_DIR)
    if not cors.empty:
        print(cors.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
