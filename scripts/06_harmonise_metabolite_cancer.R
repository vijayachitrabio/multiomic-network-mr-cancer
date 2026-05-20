#!/usr/bin/env Rscript

# Script 04b: Targeted harmonisation — Metabolite -> Cancer only
#
# Runs only the Metabolite->Cancer section of script 04, using the new
# GWAS-based instruments from script 02c (mqtl_gwas_instruments.csv).
# Run this after 02c completes; faster than re-running all of script 04.

set.seed(42)

if (!require("data.table",  quietly=TRUE)) install.packages("data.table",  repos="https://cloud.r-project.org")
if (!require("TwoSampleMR", quietly=TRUE)) install.packages("TwoSampleMR", repos="https://cloud.r-project.org")
library(data.table)
library(TwoSampleMR)

project_dir <- "."
out_dir     <- file.path(project_dir, "data", "harmonised")
cancer_dir  <- file.path(project_dir, "data", "cancer_gwas")

clean_name <- function(path) {
  x <- basename(path)
  x <- sub("\\.h\\.tsv\\.gz$", "", x)
  x <- sub("\\.tsv\\.gz$", "", x)
  x
}
label_outcome <- function(dat, nm) { dat$outcome <- nm; dat$id.outcome <- nm; dat }

# ---------------------------------------------------------------
# Load GWAS-based mQTL instruments
# ---------------------------------------------------------------
instr_file <- file.path(project_dir, "data", "mqtl", "mqtl_gwas_instruments.csv")
if (!file.exists(instr_file)) stop("Instrument file not found: run script 02c first.\n  Expected: ", instr_file)

cat("Loading mQTL instruments from:", basename(instr_file), "\n")
instr_raw <- fread(instr_file)
instr_raw <- instr_raw[!is.na(SNP) & SNP != "" & !is.na(beta) & !is.na(se) &
                       !is.na(effect_allele) & !is.na(other_allele) & !is.na(pval)]
cat(sprintf("Instruments: %d rows across %d metabolites\n",
            nrow(instr_raw), length(unique(instr_raw$metabolite))))

setDF(instr_raw)
exp_dat <- format_data(
  instr_raw,
  type              = "exposure",
  snp_col           = "SNP",
  beta_col          = "beta",
  se_col            = "se",
  eaf_col           = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col  = "other_allele",
  pval_col          = "pval",
  log_pval          = FALSE,
  samplesize_col    = "n",
  phenotype_col     = "metabolite"
)
cat(sprintf("Formatted exposure: %d rows\n\n", nrow(exp_dat)))

# ---------------------------------------------------------------
# Cancer GWAS files
# ---------------------------------------------------------------
cancer_files <- list.files(cancer_dir, pattern="\\.tsv\\.gz$", full.names=TRUE)
cat(sprintf("Found %d cancer GWAS files.\n\n", length(cancer_files)))

harm_log <- list()

for (cf in cancer_files) {
  cancer_name <- clean_name(cf)
  cat(sprintf("=== %s ===\n", cancer_name))

  out_dat_raw <- fread(cf)
  setDF(out_dat_raw)

  if ("hm_rsid" %in% colnames(out_dat_raw)) {
    out_dat <- format_data(out_dat_raw, type="outcome", snp_col="hm_rsid",
                           beta_col="hm_beta", se_col="standard_error",
                           eaf_col="hm_effect_allele_frequency",
                           effect_allele_col="hm_effect_allele",
                           other_allele_col="hm_other_allele", pval_col="p_value")
  } else if ("rsid" %in% colnames(out_dat_raw)) {
    out_dat <- format_data(out_dat_raw, type="outcome", snp_col="rsid",
                           beta_col="beta", se_col="standard_error",
                           eaf_col="effect_allele_frequency",
                           effect_allele_col="effect_allele",
                           other_allele_col="other_allele", pval_col="p_value")
  } else {
    out_dat <- format_data(out_dat_raw, type="outcome", snp_col="SNP",
                           beta_col="beta", se_col="standard_error",
                           eaf_col="effect_allele_frequency",
                           effect_allele_col="effect_allele",
                           other_allele_col="other_allele", pval_col="p_value")
  }
  out_dat <- label_outcome(out_dat, cancer_name)
  rm(out_dat_raw); gc()

  cat(sprintf("  Outcome variants loaded: %d\n", nrow(out_dat)))

  harm_dat <- tryCatch(
    harmonise_data(exposure_dat=exp_dat, outcome_dat=out_dat, action=2),
    error=function(e) { cat("  Harmonisation ERROR:", conditionMessage(e), "\n"); data.frame() }
  )

  n_before <- nrow(exp_dat)
  n_after  <- nrow(harm_dat)
  cat(sprintf("  Harmonised: %d / %d instruments retained\n", n_after, n_before))

  rds_path <- file.path(out_dir, paste0("harmonised_metabolite_", cancer_name, ".rds"))
  saveRDS(harm_dat, rds_path)
  cat(sprintf("  Saved: %s\n\n", basename(rds_path)))

  harm_log[[cancer_name]] <- data.table(
    exposure="All_Metabolites", outcome=cancer_name,
    type="Metabolite->Cancer", n_before=n_before, n_after=n_after,
    n_dropped=n_before - n_after
  )
}

# Update harmonisation log
log_file  <- file.path(out_dir, "harmonisation_log.csv")
new_log   <- rbindlist(harm_log, fill=TRUE)
if (file.exists(log_file)) {
  old_log <- fread(log_file)
  old_log <- old_log[type != "Metabolite->Cancer"]
  full_log <- rbind(old_log, new_log, fill=TRUE)
} else {
  full_log <- new_log
}
fwrite(full_log, log_file)
cat("Harmonisation log updated:", log_file, "\n")

cat("\nDone. Next: run script 07_metabolite_cancer_mr.R\n")
sessionInfo()
