# Data Sources

This pipeline relies on several large, publicly available (or controlled-access) datasets. Ensure you have obtained and properly formatted these files before running the pipeline. 

## 1. pQTL and mQTL Data
* **FinnGen**: Protein and metabolite summary statistics can be downloaded from the FinnGen summary statistic portal.
* **OpenGWAS**: Used for extraction of harmonized instruments where available (e.g., via `ieugwasr`).
* **ARIC / UKB-PPP**: Replication pQTL datasets (requires respective data access approvals).

## 2. Gynecological Cancer GWAS
* **Breast Cancer (BCAC)**: Summary statistics available from the Breast Cancer Association Consortium.
* **Endometrial Cancer**: Summary statistics available from relevant consortia/OpenGWAS.

## 3. Validation Datasets
* **TCGA (The Cancer Genome Atlas)**: Transcriptomics data (e.g., BRCA) for survival and expression validation.
* **CPTAC**: Proteomics data for breast cancer subtypes.
* **TISCH**: Single-cell RNA sequencing datasets for tumor microenvironment expression (e.g., EMTAB8107).

## 4. MAGMA Annotations
* SNP location files and gene annotation data required for MAGMA pathway analysis must be downloaded and placed in the `results/pathway/magma/inputs` folder as described in the `scripts/29_prepare_magma_inputs.R` script.
