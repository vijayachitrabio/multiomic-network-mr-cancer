# Multi-Omic Network Mendelian Randomization for Female Cancers

🚀 **Interactive Data Explorer:** [View the Live Shiny App](https://vijayachitrabio.github.io/multiomic-network-mr-cancer/)

**Vijayachitra Modhukur** et al. | University of Tartu | 2026

> 

---

## Overview

This repository contains all analysis code and output files for a proteome-wide two-sample Mendelian randomisation (MR) study integrating proteomics, metabolomics, colocalization, gene-level triangulation, mediation MR, and multi-layer biological validation to identify causal circulating protein candidates for breast, endometrial, and ovarian cancer.

**Key findings:**
- 701 circulating proteins screened using cis-pQTL instruments from FinnGen Olink (N = 619)
- 17 protein–cancer associations survived FDR correction (16 breast, 1 endometrial)
- **6 Tier 1 colocalization-supported candidates**: EFNA1, TNFRSF6B, ATRAID, FGF5, UMOD, ABO (PPH4 ≥ 0.80)
- **2 MAGMA-supported candidates**: SNX15, PM20D1 (Bonferroni-significant gene-level support)
- 5 protein → metabolite (BCAA/Glycine) → breast cancer mediation paths
- Multi-layer validation: TCGA-BRCA (N=1,097), CPTAC-BRCA (N=121), TISCH scRNA-seq, Human Protein Atlas, ARIC SomaScan + OpenGWAS INTERVAL replication
- Key methodological finding: **coloc.abf missed 2 of 8 colocalisations** (EFNA1, ATRAID); SuSiE essential for multi-signal loci
  <img width="1376" height="768" alt="image" src="https://github.com/user-attachments/assets/bda45f6b-da0e-4f17-bc87-ecc027a710f6" />


### Study Design
<img width="3672" height="2380" alt="image" src="https://github.com/user-attachments/assets/4392a5cc-ea77-4dcf-a6e8-4ea4c0c10ed8" />


---

## Repository Structure

```
├── scripts/                    # Numbered R and bash scripts (00–59)
├── data/                       # Input data (not tracked — see Data Availability)
│   ├── pqtl/                   # FinnGen Olink pQTL summary stats
│   ├── cancer_gwas/            # Breast/EC/OvC GWAS summary stats
│   └── metabolomics/           # Nightingale NMR metabolite GWAS
├── results/
│   ├── figures/                # All figures (fig1–fig15, sfig1–6)
│   ├── tables/                 # Supplementary tables (STable1–12)
│   ├── phase2_protein_cancer/  # MR screen results
│   ├── validation/             # Coloc, MAGMA, integrated evidence
│   ├── mediation/              # Two-step mediation MR
│   ├── replication/            # ARIC + OpenGWAS replication
│   ├── tcga_immune/            # TCGA-BRCA expression + immune correlations
│   ├── cptac/                  # CPTAC-BRCA proteomics
│   ├── scrna/                  # TISCH scRNA-seq cell-type enrichment
│   ├── bidirectional/          # Reverse-direction MR
│   └── mvmr/                   # MVMR feasibility
```

---

## Analysis Pipeline

Scripts are numbered and should be run sequentially:

| Script range | Stage | Description |
|---|---|---|
| `00–09` | Setup | Environment, data download, preprocessing |
| `10–19` | Phase 1 | Metabolite–cancer MR screen |
| `20–29` | Phase 2 | Protein–cancer MR screen (701 proteins × 3 cancers) |
| `30–35` | Phase 3 | Protein→metabolite MR (step 1 mediation) |
| `36–42` | Colocalization | coloc.abf + coloc.susie for priority proteins and metabolites |
| `43–44` | MAGMA | Gene-level triangulation |
| `45–46` | Figures | Coloc method comparison, ARIC replication forest |
| `47–48` | Bidirectional MR | Reverse-direction sensitivity analysis |
| `49–50` | Replication | deCODE download attempt, OpenGWAS replication figure |
| `51–59` | Validation | TCGA, CPTAC, scRNA-seq, HPA, integrated evidence |

---

## Data Availability

Raw input data are **not included** due to size and access restrictions. All sources are publicly available:

| Data | Source | Access |
|---|---|---|
| FinnGen Olink pQTL (N=619) | FinnGen R10 | [finngen.fi](https://www.finngen.fi/en/access_results) |
| Breast cancer GWAS | BCAC GCST90018757 (N=228,951) | [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/) |
| Endometrial cancer GWAS | GCST006464 (N~12,906) | [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/) |
| Ovarian cancer GWAS | GCST90016665 (N~25,509) | [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/) |
| NMR metabolomics GWAS | Nightingale Health / MRC-IEU | [IEU Open GWAS](https://gwas.mrcieu.ac.uk/) |
| ARIC SomaScan pQTL | Atherosclerosis Risk in Communities | [dbGaP](https://www.ncbi.nlm.nih.gov/gap/) |
| OpenGWAS INTERVAL SomaScan | Sun et al. 2018, Nature | [IEU Open GWAS](https://gwas.mrcieu.ac.uk/) |
| TCGA-BRCA RNA-seq | TCGA via Xena Browser | [xenabrowser.net](https://xenabrowser.net) |
| CPTAC-BRCA proteomics | CPTAC-3 (UMich) | [cptac-data-portal.georgetown.edu](https://cptac-data-portal.georgetown.edu) |
| TISCH scRNA-seq | EMTAB-8107 | [tisch.comp-genomics.org](http://tisch.comp-genomics.org) |
| Human Protein Atlas | HPA v24 | [proteinatlas.org](https://www.proteinatlas.org) |
| 1000 Genomes EUR LD | 1KGP Phase 3 | [internationalgenome.org](https://www.internationalgenome.org) |

---

## Dependencies

**R packages** (managed via `renv`):
```r
renv::restore()
```
Key packages: `TwoSampleMR` (v0.5.7), `coloc` (v5.2), `susieR`, `data.table`, `ggplot2`, `ggrepel`, `Rsamtools`, `GenomicRanges`, `magma` (external)

**Python** (for manuscript generation):
```bash
pip install python-docx
```

**External tools:**
- MAGMA v1.10 — [ctglab.nl/software/magma](https://ctglab.nl/software/magma)

---

## Quick Start

```bash
# 1. Clone
git clone https://github.com/vijayachitrabio/multiomic-network-mr-cancer.git
cd multiomic-network-mr-cancer

# 2. Restore R environment
Rscript -e "renv::restore()"

# 3. Download input data (see Data Availability above)
#    Place files in data/ following the structure in scripts/00_setup_env.R

# 4. Run pipeline
Rscript scripts/20_proteome_wide_mr_screen.R   # Phase 2 MR screen
Rscript scripts/37_coloc_snx15_pm20d1.R        # Colocalization
Rscript scripts/41_master_evidence_table.R     # Evidence integration
```

---

## Citation

> Modhukur V et al. Multi-omic triangulation of circulating proteins identifies novel breast cancer causal candidates. *Manuscript in preparation.* 2026.

---

## Contact

**Vijayachitra Modhukur**  
 University of Tartu, Estonia  
GitHub: [@vijayachitrabio](https://github.com/vijayachitrabio)

---

## License

Code: MIT License  
Data: subject to original data source access agreements (see Data Availability)
