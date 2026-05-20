# Multi-Omic Network Mendelian Randomization for Gynecological Cancers

This repository contains the analysis code and final output figures/tables for the study investigating the genetic and clinical links between common gynecological conditions (primarily breast and endometrial cancer). 

The pipeline leverages Mendelian Randomization (MR), Colocalization, Multi-omics (proteomics, metabolomics, transcriptomics), pathway enrichment (MAGMA), and single-cell RNA sequencing validation to map disease networks.

## Repository Layout

```text
scripts/               # Comprehensive R and Python scripts for MR, colocalization, validation, and figure generation
results/               # Analysis results, pathway annotations, and validation summaries
submission_2026-05-20/ # Final generated manuscript figures (main and supplementary) and tables
```

## Methodology

The codebase covers the following systematic workflow:
1. **Data Preparation**: Extraction and harmonization of pQTL and mQTL instruments alongside cancer GWAS summary statistics (e.g., FinnGen, openGWAS).
2. **Mendelian Randomization**: Comprehensive bi-directional MR analyses linking proteins, metabolites, and gynecological cancers (including ER subtypes).
3. **Colocalization**: Bayesian colocalization to prioritize shared causal variants.
4. **Validation & Replication**: Replication in independent cohorts (ARIC, UKB-PPP) and external validation (TCGA, CPTAC proteomics, TISCH scRNA-seq).
5. **Pathway Analysis**: MAGMA-based gene-set and pathway enrichment.

## Reproducibility

The `scripts/` directory contains numbered scripts (`00_setup_env.R` through `59_supplementary_mr_design_figure.R`) that document the exact execution order required to reproduce the analysis and generate the final manuscript figures. 

## Maintainer

Vijayachitra Modhukur
