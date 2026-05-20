#!/usr/bin/env Rscript

# Script 01: Download and QC UKB-PPP pQTLs
# Follows Phase 1 instructions from README

set.seed(42)

# Ensure required packages are loaded
if (!require("arrow", quietly = TRUE)) install.packages("arrow", repos = "https://cloud.r-project.org")
if (!require("here", quietly = TRUE)) install.packages("here", repos = "https://cloud.r-project.org")
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)
library(arrow)
library(here)

# Define paths
out_dir <- here::here("data", "pqtl")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_file <- file.path(out_dir, "pqtl_instruments.csv")

cat("Connecting to UKB-PPP AWS Open Data bucket...\n")
# UKB-PPP is available at s3://open-data-ukbppp/
# The discovery sumstats are typically partitioned by protein or chunk.
# Open connection to the parquet dataset using explicit bucket config
dataset <- tryCatch({
  s3 <- arrow::S3FileSystem$create(anonymous = TRUE, region = "us-east-1")
  arrow::open_dataset(s3$path("ukbiobank.opendata.sagebase.org/discovery/ukbppp_gwas_discovery_results/"))
}, error = function(e) {
  cat("Error connecting to S3. Please ensure 'arrow' is built with S3 support.\n")
  stop(e)
})

cat("Querying and filtering cis-pQTLs...\n")
# Filter: 
# 1. Cis only (often flagged in the dataset or defined by distance to TSS)
# 2. Exclude MHC region (chr 6: 28Mb - 34Mb)
# Note: Schema depends on exact UKB-PPP parquet structure. Assumed common naming convention.
pqtl_data <- dataset |>
  filter(
    cis_trans == "cis" | (distance_to_gene <= 1000000),  # Assuming distance/flag is available
    !(chromosome == "6" & base_pair_location >= 28000000 & base_pair_location <= 34000000),
    p_value < 5e-8 # Initial strict threshold for instruments
  ) |>
  collect() |>
  as.data.table()

cat(sprintf("Extracted %d initial cis-pQTL associations.\n", nrow(pqtl_data)))

# Compute F-statistic
cat("Computing F-statistics...\n")
pqtl_data[, F_stat := (beta^2) / (standard_error^2)]

# Filter F > 10
pqtl_data <- pqtl_data[F_stat > 10]
cat(sprintf("Retained %d associations with F > 10.\n", nrow(pqtl_data)))

# LD clumping per protein
cat("Performing LD clumping (r2 < 0.001, 10Mb window)...\n")
# Standardize names for TwoSampleMR
pqtl_for_clump <- format_data(
  pqtl_data,
  type = "exposure",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  phenotype_col = "protein_id",
  chr_col = "chromosome",
  pos_col = "base_pair_location"
)

# Clump using ieugwasr (which queries the OpenGWAS API or local PLINK)
clumped_pqtl <- clump_data(
  pqtl_for_clump,
  clump_r2 = 0.001,
  clump_kb = 10000,
  pop = "EUR"
)

cat(sprintf("After LD clumping, %d independent cis-pQTL instruments remain.\n", nrow(clumped_pqtl)))

# Save to file
fwrite(clumped_pqtl, out_file)
saveRDS(clumped_pqtl, file.path(out_dir, "pqtl_instruments.rds"))

cat("Script 01 complete.\n")
sessionInfo()
