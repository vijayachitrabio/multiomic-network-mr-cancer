#!/usr/bin/env Rscript

# 74_gckr_sensitivity_lolo.R
# Performs a Leave-One-Locus-Out sensitivity analysis on the GCKR locus (rs1260326)
# for metabolite -> breast cancer MR. This tests if the mediation effect is entirely
# driven by this massive pleiotropic locus.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
})

cat("Starting GCKR (rs1260326) sensitivity analysis...\n")

harm_file <- "data/harmonised/harmonised_metabolite_Breast_GCST90018757.rds"
if (!file.exists(harm_file)) stop("Harmonised data not found: ", harm_file)

harm_dat <- readRDS(harm_file)

# Filter for the metabolites identified in the mediation pathways
target_mets <- c("Gly", "Total_BCAA", "TG_by_PG")
harm_dat <- harm_dat[harm_dat$exposure %in% target_mets, ]

if (nrow(harm_dat) == 0) stop("No target metabolites found in harmonised data.")

cat(sprintf("Loaded %d instruments across %d metabolites.\n", nrow(harm_dat), length(unique(harm_dat$exposure))))

# Run original MR
cat("Running original MR...\n")
res_orig <- mr(harm_dat, method_list = c("mr_ivw", "mr_wald_ratio"))

# The GCKR locus lead SNP is rs1260326 (chr 2). 
gckr_snp <- "rs1260326"

# Check if rs1260326 is actually present
present <- gckr_snp %in% harm_dat$SNP
cat(sprintf("Is rs1260326 present in instruments? %s\n", present))

# Run MR without GCKR
cat("Running GCKR-excluded MR...\n")
harm_dat_no_gckr <- harm_dat[harm_dat$SNP != gckr_snp, ]
res_no_gckr <- mr(harm_dat_no_gckr, method_list = c("mr_ivw", "mr_wald_ratio"))

# Format results side-by-side
res_orig <- setDT(res_orig)
res_no_gckr <- setDT(res_no_gckr)

combined <- merge(
  res_orig[, .(exposure, nsnp, b, pval)],
  res_no_gckr[, .(exposure, nsnp_no_gckr = nsnp, b_no_gckr = b, pval_no_gckr = pval)],
  by = "exposure",
  all = TRUE
)

combined[, b_change_pct := round(((b_no_gckr - b) / b) * 100, 1)]
combined[, interpretation := ifelse(pval_no_gckr < 0.05, 
                                    "Signal remains directionally consistent",
                                    "Signal disappears (GCKR-driven)")]

out_file <- "results/tables/STable_GCKR_Sensitivity.csv"
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
fwrite(combined, out_file)

cat("\nGCKR Sensitivity Results:\n")
print(combined)
cat("\nSaved to:", out_file, "\n")
