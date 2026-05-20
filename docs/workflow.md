# Analytical Workflow

The pipeline is designed to be executed sequentially, from script `00` to `59`. Below is an overview of the key phases in the analysis:

## Phase 1: Environment and Setup
* `00_setup_env.R`: Initializes libraries and folder structures.

## Phase 2: Data Downloading and Extraction
* `01_download_pqtl.R` & `01b_download_finngen_olink_pqtl.R`: Fetches primary pQTL data.
* `02c_extract_gwas_instruments.R`: Identifies and formats independent genetic instruments.
* `03_download_cancer_gwas.R`: Downloads the outcome GWAS (Breast/Endometrial cancer).

## Phase 3: Harmonization
* `04_harmonise_all.R` to `04d_harmonise_protein_metabolite.R`: Harmonizes the exposure instruments with the various outcomes.

## Phase 4: Mendelian Randomization (MR)
* `05_protein_cancer_mr.R`: Primary MR analysis of proteins on cancer risk.
* `07_metabolite_cancer_mr.R`: Primary MR analysis of metabolites on cancer risk.
* `09_protein_metabolite_mr.R`: MR analysis investigating protein-metabolite networks.

## Phase 5: Colocalization and Steiger Filtering
* `06_protein_cancer_coloc.R` & `08_metabolite_cancer_coloc.R`: Bayesian colocalization for prioritization.
* `20_steiger_directionality.R`: Steiger filtering to establish causal direction.

## Phase 6: Pathway Analysis (MAGMA)
* `29_prepare_magma_inputs.R` to `31_summarise_magma_results.R`: Runs gene-set and pathway enrichment using MAGMA.

## Phase 7: Replication and Validation
* `34_ukbppp_replication_mr.R` & `38_aric_pqtl_replication_mr.R`: Replicates findings in UKB-PPP and ARIC.
* `51_tcga_brca_expression_immune_validation.R` to `53_tisch_scrna_brca_validation.R`: In silico clinical validation using transcriptomic/proteomic profiles.

## Phase 8: Evidence Integration & Figures
* `41_master_evidence_table.R`: Compiles all tiers of evidence.
* `54_integrated_evidence_map.R` to `59_supplementary_mr_design_figure.R`: Generates final summary figures.
