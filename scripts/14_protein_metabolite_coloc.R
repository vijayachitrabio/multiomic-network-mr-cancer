#!/usr/bin/env Rscript

# Script 10: Colocalization (Protein -> Metabolite)

set.seed(42)

if (!require("here", quietly = TRUE)) install.packages("here", repos = "https://cloud.r-project.org")
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(coloc)
library(here)

res_dir <- here::here("results", "phase4_protein_metabolite")
out_dir <- here::here("results", "phase4_protein_metabolite", "coloc")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

sig_file <- file.path(res_dir, "protein_metabolite_mr_results_significant.csv")
if (!file.exists(sig_file)) {
  cat("No significant MR results found. Skipping colocalization.\n")
  quit(status = 0)
}

sig_res <- fread(sig_file)
cat(sprintf("Found %d significant pairs. Running colocalization framework...\n", nrow(sig_res)))

coloc_results <- data.table()

for (i in 1:nrow(sig_res)) {
  exp_name <- sig_res$exposure[i]
  out_name <- sig_res$outcome[i]
  
  coloc_results <- rbind(coloc_results, list(
    exposure = exp_name,
    outcome = out_name,
    PP.H3 = NA_real_,
    PP.H4 = NA_real_,
    status = "Pending regional extraction"
  ))
}

fwrite(coloc_results, file.path(out_dir, "protein_metabolite_coloc_summary.csv"))
cat("Colocalization framework setup complete.\n")
