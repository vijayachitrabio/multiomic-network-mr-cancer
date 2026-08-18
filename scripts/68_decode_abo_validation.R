#!/usr/bin/env Rscript

# Validate ABO protein causality using deCODE Genetics summary statistics

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(dplyr)
})

project_dir <- "."
decode_dir <- "../SUMMARY DATA-proteins"
decode_file <- file.path(decode_dir, "9253_52_ABO_BGAT.txt.gz")
out_dir <- file.path(project_dir, "results", "decode")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("1. Extracting cis-pQTLs for ABO from deCODE...\n")
# ABO is on Chromosome 9. We extract +/- 1.5MB around the gene (chr9:133.2M)
# and only keep genome-wide significant SNPs (Pval < 5e-8)
# Chrom=1, Pos=2, Pval=8
cmd <- sprintf("gzcat '%s' | awk 'NR==1 || (($1==\"chr9\" || $1==\"9\") && $2>=131500000 && $2<=135000000 && $8 < 5e-8)'", decode_file)
cat("Running command:", cmd, "\n")
pqtl_raw <- fread(cmd)

cat(sprintf("Extracted %d genome-wide significant cis-SNPs from deCODE.\n", nrow(pqtl_raw)))

if (nrow(pqtl_raw) == 0) {
  stop("No significant cis-pQTLs found for ABO in deCODE.")
}

cat("2. Formatting exposure data for TwoSampleMR...\n")
# Columns: Chrom Pos Name rsids effectAllele otherAllele Beta Pval min_log10_pval SE N ImpMAF
# Name is the unique variant name, we can use rsids if available, otherwise Name.
pqtl_raw[, SNP := ifelse(rsids != "" & !is.na(rsids), rsids, Name)]

exp_dat <- format_data(
  as.data.frame(pqtl_raw),
  type = "exposure",
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF",
  pval_col = "Pval",
  samplesize_col = "N",
  chr_col = "Chrom",
  pos_col = "Pos"
)
exp_dat$exposure <- "ABO (deCODE)"

# The original manuscript uses two specific instruments for ABO:
target_snps <- c("rs975381715", "rs1381383189")
exp_clumped <- exp_dat[exp_dat$SNP %in% target_snps, ]
cat(sprintf("Retained %d specific instrument(s) for ABO.\n", nrow(exp_clumped)))

cat("4. Extracting outcome data (FinnGen Endometrial Cancer)...\n")
outcome_file <- file.path(project_dir, "data", "cancer_gwas", "Endometrial_GCST006464.h.tsv.gz")
out_raw <- fread(outcome_file)

# Only keep exposure SNPs that are actually present in the outcome!
shared_snps <- intersect(exp_clumped$SNP, out_raw$hm_rsid)
exp_clumped <- exp_clumped[exp_clumped$SNP %in% shared_snps, ]

out_dat <- format_data(
  as.data.frame(out_raw[hm_rsid %in% exp_clumped$SNP]),
    type = "outcome",
    snp_col = "hm_rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
  )

out_dat$outcome <- "Endometrial Cancer"

cat("5. Harmonizing data...\n")
harm_dat <- harmonise_data(exposure_dat = exp_clumped, outcome_dat = out_dat)

cat("6. Running Mendelian Randomization...\n")
mr_res <- mr(harm_dat, method_list = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median"))

# Calculate F-statistic
harm_dat$F_stat <- (harm_dat$beta.exposure^2) / (harm_dat$se.exposure^2)

# Generate ORs
mr_res <- generate_odds_ratios(mr_res)

print(mr_res)

fwrite(mr_res, file.path(out_dir, "decode_abo_mr_results.csv"))
fwrite(harm_dat, file.path(out_dir, "decode_abo_mr_harmonised.csv"))
cat("Success! Results saved to", out_dir, "\n")
