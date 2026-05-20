#!/usr/bin/env Rscript

# Script 04: Harmonise all exposure-outcome pairs
# Follows Phase 1 instructions from README

set.seed(42)

if (!require("here", quietly = TRUE)) install.packages("here", repos = "https://cloud.r-project.org")
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(here)

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"

clean_name <- function(path) {
  x <- basename(path)
  x <- sub("\\.h\\.tsv\\.gz$", "", x)
  x <- sub("\\.tsv\\.gz$", "", x)
  x <- sub("\\.txt\\.gz$", "", x)
  x <- sub("\\.vcf\\.gz$", "", x)
  x <- sub("_std\\.csv$", "", x)
  x <- sub("\\.csv$", "", x)
  x
}

label_outcome <- function(dat, outcome_name) {
  if (nrow(dat) > 0) {
    dat$outcome <- outcome_name
    dat$id.outcome <- outcome_name
  }
  dat
}

make_variant_key <- function(chr, pos, allele1, allele2) {
  chr <- sub("^chr", "", as.character(chr))
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  a_min <- ifelse(allele1 <= allele2, allele1, allele2)
  a_max <- ifelse(allele1 <= allele2, allele2, allele1)
  paste(chr, pos, a_min, a_max, sep=":")
}

uses_variant_key_exposure <- function(dat) {
  "variant_id.exposure" %in% names(dat) &&
    "SNP" %in% names(dat) &&
    any(grepl("^[0-9XYMT]+:", dat$SNP))
}

out_dir <- file.path(project_dir, "data", "harmonised")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

log_file <- file.path(out_dir, "harmonisation_log.csv")
# Initialize log
harm_log <- data.table(exposure=character(), outcome=character(), type=character(), n_before=integer(), n_after=integer(), n_dropped=integer())

cat("Loading exposure instruments...\n")
pqtl_file <- file.path(project_dir, "data", "pqtl", "pqtl_instruments.csv")
if(file.exists(pqtl_file)) {
  pqtls <- fread(pqtl_file)
  if ("beta.exposure" %in% colnames(pqtls)) {
    exp_dat <- pqtls
  } else {
    exp_dat <- format_data(pqtls, type="exposure")
  }
  exp_uses_variant_key <- uses_variant_key_exposure(exp_dat)
  if (exp_uses_variant_key) {
    cat("Detected variant-key pQTL instruments; outcomes will be matched by chr:pos:allele-pair keys.\n")
  }
} else {
  cat("Warning: pQTL instruments not found. Skipping pQTL harmonisation.\n")
}

# 1. Protein -> Cancer Harmonisation
cat("Starting Protein -> Cancer harmonisation...\n")
cancer_dir <- file.path(project_dir, "data", "cancer_gwas")
cancer_files <- list.files(cancer_dir, pattern = "\\.tsv\\.gz$|\\.csv$|\\.txt\\.gz$|\\.vcf\\.gz$", full.names = TRUE)

if(exists("exp_dat") && length(cancer_files) > 0) {
  for (cf in cancer_files) {
    cancer_name <- clean_name(cf)
    cat(sprintf("Processing outcome: %s\n", cancer_name))
    
    out_dat_raw <- fread(cf)
    if (isTRUE(exp_uses_variant_key) &&
        all(c("chromosome", "base_pair_location", "effect_allele", "other_allele") %in% colnames(out_dat_raw))) {
      out_dat_raw[, SNP_variant_key := make_variant_key(
        chromosome, base_pair_location, effect_allele, other_allele
      )]
      out_dat_raw <- out_dat_raw[SNP_variant_key %in% exp_dat$SNP]
      setDF(out_dat_raw)
      out_dat <- format_data(
        out_dat_raw,
        type="outcome",
        snp_col="SNP_variant_key",
        beta_col="beta",
        se_col="standard_error",
        eaf_col="effect_allele_frequency",
        effect_allele_col="effect_allele",
        other_allele_col="other_allele",
        pval_col="p_value"
      )
      # format_data() lowercases SNP col; restore uppercase to match variant-key exp_dat$SNP
      out_dat$SNP <- toupper(out_dat$SNP)
    } else if ("hm_rsid" %in% colnames(out_dat_raw)) {
      setDF(out_dat_raw)
      out_dat <- format_data(out_dat_raw, type="outcome", snp_col="hm_rsid", beta_col="hm_beta", se_col="standard_error", eaf_col="hm_effect_allele_frequency", effect_allele_col="hm_effect_allele", other_allele_col="hm_other_allele", pval_col="p_value")
    } else if ("rsid" %in% colnames(out_dat_raw)) {
      setDF(out_dat_raw)
      out_dat <- format_data(out_dat_raw, type="outcome", snp_col="rsid", beta_col="beta", se_col="standard_error", eaf_col="effect_allele_frequency", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p_value")
    } else {
      setDF(out_dat_raw)
      out_dat <- format_data(out_dat_raw, type="outcome", snp_col="SNP", beta_col="beta", se_col="standard_error", eaf_col="effect_allele_frequency", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p_value")
    }
    out_dat <- label_outcome(out_dat, cancer_name)

    # Harmonise
    harm_dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2)
    
    # Save RDS per protein (or in one large file if small enough, but README says "per protein-cancer")
    # To avoid 2923 files per cancer, we save one grouped file per cancer which can be subset later, 
    # OR save per protein as requested. Saving grouped is much faster.
    saveRDS(harm_dat, file.path(out_dir, paste0("harmonised_protein_", cancer_name, ".rds")))
    
    # Log
    n_bef <- nrow(exp_dat)
    n_aft <- nrow(harm_dat)
    harm_log <- rbind(harm_log, list("All_Proteins", cancer_name, "Protein->Cancer", n_bef, n_aft, n_bef - n_aft))
  }
}

# 2. Protein -> Metabolite Harmonisation
cat("Starting Protein -> Metabolite harmonisation...\n")
mqtl_gwas_dir <- file.path(project_dir, "data", "mqtl", "mqtl_full_gwas")
mqtl_files <- list.files(mqtl_gwas_dir, pattern = "_full_regenie\\.tsv\\.gz$", full.names = TRUE)

if(exists("exp_dat") && length(mqtl_files) > 0) {
  for (mf in mqtl_files) {
    met_name <- gsub("_full_regenie\\.tsv\\.gz$", "", basename(mf))
    cat(sprintf("Processing outcome: %s\n", met_name))
    
    out_dat_raw <- fread(
      mf,
      select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
                 "N", "BETA", "SE", "LOG10P")
    )
    if (isTRUE(exp_uses_variant_key)) {
      out_dat_raw[, SNP_variant_key := make_variant_key(CHROM, GENPOS, ALLELE0, ALLELE1)]
      out_dat_raw <- out_dat_raw[SNP_variant_key %in% exp_dat$SNP]
      snp_col <- "SNP_variant_key"
    } else {
      out_dat_raw <- out_dat_raw[ID %in% exp_dat$SNP]
      snp_col <- "ID"
    }
    setDF(out_dat_raw)

    if (nrow(out_dat_raw) > 0) {
      # Format according to Alasoo/Rahu regenie columns.
      out_dat <- format_data(
        out_dat_raw,
        type = "outcome",
        snp_col = snp_col,
        beta_col = "BETA",
        se_col = "SE",
        eaf_col = "A1FREQ",
        effect_allele_col = "ALLELE1",
        other_allele_col = "ALLELE0",
        pval_col = "LOG10P",
        log_pval = TRUE,
        samplesize_col = "N",
        chr_col = "CHROM",
        pos_col = "GENPOS"
      )
      out_dat <- label_outcome(out_dat, met_name)
      # format_data() lowercases SNP col; restore uppercase to match variant-key exp_dat$SNP
      if (isTRUE(exp_uses_variant_key)) out_dat$SNP <- toupper(out_dat$SNP)
      harm_dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2)
    } else {
      harm_dat <- data.frame()
    }

    saveRDS(harm_dat, file.path(out_dir, paste0("harmonised_protein_", met_name, ".rds")))
    
    n_bef <- nrow(exp_dat)
    n_aft <- nrow(harm_dat)
    harm_log <- rbind(harm_log, list("All_Proteins", met_name, "Protein->Metabolite", n_bef, n_aft, n_bef - n_aft))
  }
}

# 3. Metabolite -> Cancer Harmonisation
# Uses GWAS-based instruments (p<5e-8 + 500kb clumping) from script 02c.
# Replaces broken SuSiE approach where variant coordinates didn't overlap
# between SuSiE credible set files and full GWAS summary stats.
cat("Starting Metabolite -> Cancer harmonisation...\n")
gwas_instr_file <- file.path(project_dir, "data", "mqtl", "mqtl_gwas_instruments.csv")
susie_file      <- file.path(project_dir, "data", "mqtl", "mqtl_susie_instruments.csv")

mqtl_instr_file <- if (file.exists(gwas_instr_file)) gwas_instr_file else susie_file
cat(sprintf("Using mQTL instrument file: %s\n", basename(mqtl_instr_file)))

if(file.exists(mqtl_instr_file) && length(cancer_files) > 0) {
  susie_dat_raw <- fread(mqtl_instr_file)
  n_susie_raw <- nrow(susie_dat_raw)

  # Detect column name scheme (new GWAS-based vs old SuSiE-based)
  if ("SNP" %in% names(susie_dat_raw)) {
    # New GWAS-based instrument file (from script 02c)
    core_cols <- c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "pval", "metabolite")
    missing_cols <- setdiff(core_cols, names(susie_dat_raw))
    if (length(missing_cols) > 0) stop("mQTL instrument file missing columns: ", paste(missing_cols, collapse=", "))
    susie_dat_raw <- susie_dat_raw[!is.na(SNP) & SNP != "" & !is.na(beta) & !is.na(se) &
                                   !is.na(effect_allele) & !is.na(other_allele) & !is.na(pval)]
    cat(sprintf("Retained %d complete instruments from %d rows.\n", nrow(susie_dat_raw), n_susie_raw))
    setDF(susie_dat_raw)
    susie_exp <- format_data(
      susie_dat_raw,
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
  } else {
    # Old SuSiE-based instrument file (fallback)
    core_cols <- c("ID", "BETA", "SE", "A1FREQ", "ALLELE1", "ALLELE0", "LOG10P", "metabolite")
    missing_cols <- setdiff(core_cols, names(susie_dat_raw))
    if (length(missing_cols) > 0) stop("mQTL instrument file missing columns: ", paste(missing_cols, collapse=", "))
    susie_dat_raw <- susie_dat_raw[!is.na(ID) & ID != "" & !is.na(BETA) & !is.na(SE) &
                                   !is.na(ALLELE1) & !is.na(ALLELE0) & !is.na(LOG10P)]
    cat(sprintf("Retained %d complete SuSiE instruments from %d rows.\n", nrow(susie_dat_raw), n_susie_raw))
    setDF(susie_dat_raw)
    if (!"sample_n" %in% names(susie_dat_raw)) susie_dat_raw$sample_n <- 246683
    susie_exp <- format_data(
      susie_dat_raw,
      type              = "exposure",
      snp_col           = "ID",
      beta_col          = "BETA",
      se_col            = "SE",
      eaf_col           = "A1FREQ",
      effect_allele_col = "ALLELE1",
      other_allele_col  = "ALLELE0",
      pval_col          = "LOG10P",
      log_pval          = TRUE,
      samplesize_col    = "sample_n",
      phenotype_col     = "metabolite"
    )
  }
  
  for (cf in cancer_files) {
    cancer_name <- clean_name(cf)
    cat(sprintf("Processing outcome: %s\n", cancer_name))
    
    out_dat_raw <- fread(cf)
    setDF(out_dat_raw)
    if ("hm_rsid" %in% colnames(out_dat_raw)) {
      out_dat <- format_data(out_dat_raw, type="outcome", snp_col="hm_rsid", beta_col="hm_beta", se_col="standard_error", eaf_col="hm_effect_allele_frequency", effect_allele_col="hm_effect_allele", other_allele_col="hm_other_allele", pval_col="p_value")
    } else if ("rsid" %in% colnames(out_dat_raw)) {
      out_dat <- format_data(out_dat_raw, type="outcome", snp_col="rsid", beta_col="beta", se_col="standard_error", eaf_col="effect_allele_frequency", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p_value")
    } else {
      out_dat <- format_data(out_dat_raw, type="outcome", snp_col="SNP", beta_col="beta", se_col="standard_error", eaf_col="effect_allele_frequency", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p_value")
    }
    out_dat <- label_outcome(out_dat, cancer_name)
    
    harm_dat <- harmonise_data(exposure_dat = susie_exp, outcome_dat = out_dat, action = 2)
    saveRDS(harm_dat, file.path(out_dir, paste0("harmonised_metabolite_", cancer_name, ".rds")))
    
    n_bef <- nrow(susie_exp)
    n_aft <- nrow(harm_dat)
    harm_log <- rbind(harm_log, list("All_Metabolites", cancer_name, "Metabolite->Cancer", n_bef, n_aft, n_bef - n_aft))
  }
}

fwrite(harm_log, log_file)
cat("Script 04 complete.\n")
sessionInfo()
