# Multi-Omic Network Mendelian Randomization for Female Cancers

This repository contains the analysis code and final output figures/tables for the study investigating the genetic and clinical links between common gynecological conditions (primarily breast and endometrial cancer). 

The pipeline leverages Mendelian Randomization (MR), Colocalization, Multi-omics (proteomics, metabolomics, transcriptomics), pathway enrichment (MAGMA), and single-cell RNA sequencing validation to map disease networks.

## Quick Start

1. Clone the repository: `git clone https://github.com/vijayachitrabio/multiomic-network-mr-cancer.git`
2. Open the R project and restore dependencies: `renv::restore()`
3. Download the required data (see Data Availability below).
4. Run the scripts sequentially starting from `scripts/00_setup_env.R`.

## Dependencies

This project relies on several R and Python packages. 
R package dependencies are managed using `renv`. You can install all required packages by running `renv::restore()`.

## Data Availability

**Note:** Raw data files and intermediate results are *not* included in this repository due to size constraints and data access agreements. You must download the required input data (pQTLs, mQTLs, GWAS summary stats) before running the pipeline. 

For detailed information on data sources and how to acquire them, please see [docs/data_sources.md](docs/data_sources.md).

## Reproducibility

To ensure full reproducibility, the analytical pipeline has been strictly organized into sequential scripts. All random seeds are set within the scripts where applicable. 

## Script Execution Order

The `scripts/` directory contains numbered scripts (`00_` through `59_`) that document the exact execution order. For a detailed breakdown of the workflow and what each script does, please consult [docs/workflow.md](docs/workflow.md).

## Outputs

Generated figures and tables are structured as follows:
* `manuscript_outputs/figures/main`: Main manuscript figures
* `manuscript_outputs/figures/supplementary`: Supplementary figures
* `manuscript_outputs/tables/supplementary`: Supplementary tables

## Maintainer

Vijayachitra Modhukur
