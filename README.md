# Multi-omic triangulation of circulating proteins and metabolites in breast cancer susceptibility

**Interactive Data Explorer:** [View the live app](https://vijayachitrabio.github.io/multiomic-network-mr-cancer/)

**Vijayachitra Modhukur** et al. | University of Tartu | 2026

---
<img width="1536" height="1024" alt="image" src="https://github.com/user-attachments/assets/9818a24c-0c54-417f-be71-aa70941815db" />


## Overview

This repository contains analysis code, summary outputs, and manuscript-supporting materials for a proteome-wide Mendelian randomization study integrating circulating protein genetics, cancer GWAS, metabolite colocalization, mediation MR, gene-level triangulation, and tumour-context annotation.

Prior proteome-wide MR studies have established the value of plasma protein genetics for breast cancer target prioritization [4,5]. The present study extends that literature by moving beyond association-level discovery toward multi-layer evidence triangulation, integrating MR with fine-mapping-aware colocalization, MAGMA gene-level support, metabolite colocalization and mediation, cross-platform pQTL sensitivity analyses, and tumour-context annotation. This distinction is important because genetically proxied protein–cancer associations may arise from linkage disequilibrium, assay-specific pQTL architecture, or pleiotropic regional effects. The evidence hierarchy used here was therefore designed to separate protein candidates with stronger shared genetic support from secondary signals supported by gene-level, replication, metabolic, or tumour-context evidence alone.

The study prioritizes genetically supported circulating protein and metabolic pathways associated with breast cancer susceptibility, with endometrial and ovarian cancers included as hormone-related comparator outcomes.

The analysis combines:

- cis-pQTL Mendelian randomization of circulating proteins
- protein–cancer colocalization using coloc.abf and SuSiE
- MAGMA gene-level triangulation
- metabolite–cancer colocalization
- two-step protein → metabolite → breast cancer mediation MR
- cross-platform pQTL/MR replication where instruments were available
- tumour-context annotation using TCGA-BRCA, CPTAC-BRCA, TISCH single-cell RNA-seq, and Human Protein Atlas resources

---

## Key findings

- 701 circulating proteins were screened using cis-pQTL instruments from the FinnGen R10 Olink panel.
- Seventeen protein–cancer associations survived false discovery rate correction: 16 for breast cancer and one for endometrial cancer.
- Colocalization-supported breast cancer candidates included **EFNA1**, **TNFRSF6B**, **ATRAID**, and **FGF5**.
- **UMOD** showed a provisional breast cancer signal supported by coloc.abf only.
- **ABO** represented a distinct endometrial cancer comparator signal.
- **SNX15** and **PM20D1** were supported by MR and Bonferroni-significant MAGMA gene-level evidence but were not classified as colocalized protein–cancer signals.
- Protein–metabolite–breast cancer mediation analyses highlighted pathways involving branched-chain amino acid and glycine-related metabolic traits.
- Tumour-context analyses supported immune, stromal, and metabolic interpretations for selected prioritized proteins.
- A methodological finding was that coloc.abf missed or misclassified colocalization at multi-signal loci such as **EFNA1** and **ATRAID**, whereas SuSiE-based colocalization resolved shared signals.

---

## Study design



The workflow includes:

1. Proteome-wide cis-pQTL MR of circulating proteins against breast, endometrial, and ovarian cancer GWAS.
2. Multiple-testing correction across protein–cancer pairs.
3. Protein–cancer colocalization using coloc.abf and SuSiE.
4. MAGMA gene-level triangulation.
5. Metabolite–cancer colocalization across NMR metabolic traits.
6. Two-step protein → metabolite → breast cancer mediation MR.
7. Cross-platform replication using ARIC SomaScan and OpenGWAS INTERVAL where available.
8. Tumour-context annotation using TCGA, CPTAC, TISCH, and Human Protein Atlas resources.

---

## Repository structure

```text
├── scripts/                    # Numbered R/bash scripts
├── data/                       # Input data, not tracked
│   ├── pqtl/                   # FinnGen Olink pQTL summary statistics
│   ├── cancer_gwas/            # Breast, endometrial, and ovarian cancer GWAS
│   └── metabolomics/           # NMR metabolite GWAS
├── results/
│   ├── figures/                # Main and supplementary figures
│   ├── tables/                 # Supplementary and manuscript tables
│   ├── phase2_protein_cancer/  # Proteome-wide MR screen results
│   ├── validation/             # Colocalization, MAGMA, integrated evidence
│   ├── mediation/              # Two-step mediation MR
│   ├── replication/            # ARIC, OpenGWAS, and external pQTL replication
│   ├── tcga_immune/            # TCGA-BRCA expression and immune correlations
│   ├── cptac/                  # CPTAC-BRCA proteomics
│   ├── scrna/                  # TISCH single-cell RNA-seq annotation
│   ├── bidirectional/          # Reverse-direction MR sensitivity analysis
│   └── mvmr/                   # MVMR feasibility assessment
