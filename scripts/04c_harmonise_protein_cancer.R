#!/usr/bin/env Rscript

# Script 04c: Protein -> Cancer harmonisation only (FinnGen R10 Olink 701 proteins)
#
# Fixes: format_data() lowercases the SNP_variant_key column; toupper() restores
# uppercase so exp_dat$SNP matches out_dat$SNP for harmonise_data().
#
# Runs ~2-5 min per cancer GWAS (1062 pQTL instruments, 3 cancer files).
# Use after 04_harmonise_all.R when only the protein->cancer section needs re-running.

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
out_dir     <- file.path(project_dir, "data", "harmonised")
cancer_dir  <- file.path(project_dir, "data", "cancer_gwas")

make_variant_key <- function(chr, pos, allele1, allele2) {
  chr    <- sub("^chr", "", as.character(chr))
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  a_min  <- ifelse(allele1 <= allele2, allele1, allele2)
  a_max  <- ifelse(allele1 <= allele2, allele2, allele1)
  paste(chr, pos, a_min, a_max, sep = ":")
}

clean_name <- function(path) {
  x <- basename(path)
  x <- sub("\\.h\\.tsv\\.gz$", "", x)
  x <- sub("\\.tsv\\.gz$",     "", x)
  x
}

# ---------------------------------------------------------------
# Load FinnGen R10 pQTL instruments
# ---------------------------------------------------------------
pqtl_file <- file.path(project_dir, "data", "pqtl", "pqtl_instruments.csv")
if (!file.exists(pqtl_file)) stop("pQTL instrument file not found: ", pqtl_file)

cat("Loading pQTL instruments...\n")
pqtls   <- fread(pqtl_file)
exp_dat <- pqtls          # already has beta.exposure, se.exposure etc.
cat(sprintf("  %d instruments across %d proteins\n",
            nrow(exp_dat), length(unique(exp_dat$exposure))))

# Confirm variant-key format
stopifnot("SNP" %in% names(exp_dat), any(grepl("^[0-9XYM]+:", exp_dat$SNP)))
cat("  SNP format (sample):", paste(head(exp_dat$SNP, 3), collapse = ", "), "\n\n")

# ---------------------------------------------------------------
# Protein -> Cancer harmonisation
# ---------------------------------------------------------------
cancer_files <- list.files(cancer_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
cat(sprintf("Found %d cancer GWAS files.\n\n", length(cancer_files)))

harm_log <- list()

for (cf in cancer_files) {
  cancer_name <- clean_name(cf)
  cat(sprintf("=== %s ===\n", cancer_name))

  out_dat_raw <- fread(cf)
  cat(sprintf("  Loaded %d variants\n", nrow(out_dat_raw)))

  # Determine which columns are available for variant-key construction
  has_pos_cols <- all(c("chromosome", "base_pair_location", "effect_allele", "other_allele") %in%
                        colnames(out_dat_raw))
  has_hm_cols  <- all(c("hm_chrom", "hm_pos", "hm_effect_allele", "hm_other_allele") %in%
                        colnames(out_dat_raw))

  if (has_hm_cols) {
    # Prefer harmonised GRCh38 columns when available
    out_dat_raw[, SNP_variant_key := make_variant_key(
      hm_chrom, hm_pos, hm_effect_allele, hm_other_allele)]
    ea_col  <- "hm_effect_allele"
    oa_col  <- "hm_other_allele"
    b_col   <- "hm_beta"
    se_col_ <- "standard_error"
    eaf_col <- "hm_effect_allele_frequency"
    pv_col  <- "p_value"
  } else if (has_pos_cols) {
    out_dat_raw[, SNP_variant_key := make_variant_key(
      chromosome, base_pair_location, effect_allele, other_allele)]
    ea_col  <- "effect_allele"
    oa_col  <- "other_allele"
    b_col   <- "beta"
    se_col_ <- "standard_error"
    eaf_col <- "effect_allele_frequency"
    pv_col  <- "p_value"
  } else {
    cat("  WARNING: no positional columns for variant-key matching — skipping.\n\n")
    next
  }

  out_dat_raw <- out_dat_raw[SNP_variant_key %in% exp_dat$SNP]
  cat(sprintf("  After variant-key filter: %d rows\n", nrow(out_dat_raw)))

  if (nrow(out_dat_raw) == 0) {
    cat("  No matching variants — writing empty RDS.\n\n")
    saveRDS(data.frame(), file.path(out_dir, paste0("harmonised_protein_", cancer_name, ".rds")))
    harm_log[[cancer_name]] <- data.table(
      exposure="All_Proteins", outcome=cancer_name,
      type="Protein->Cancer", n_before=nrow(exp_dat), n_after=0, n_dropped=nrow(exp_dat))
    next
  }

  setDF(out_dat_raw)
  out_dat <- format_data(
    out_dat_raw,
    type             = "outcome",
    snp_col          = "SNP_variant_key",
    beta_col         = b_col,
    se_col           = se_col_,
    eaf_col          = eaf_col,
    effect_allele_col = ea_col,
    other_allele_col  = oa_col,
    pval_col         = pv_col
  )
  # format_data() lowercases the SNP column — restore uppercase to match exp_dat$SNP
  out_dat$SNP       <- toupper(out_dat$SNP)
  out_dat$outcome   <- cancer_name
  out_dat$id.outcome <- cancer_name

  cat(sprintf("  format_data outcome nrow: %d\n", nrow(out_dat)))
  cat(sprintf("  SNP overlap exp vs out: %d\n", length(intersect(exp_dat$SNP, out_dat$SNP))))

  harm_dat <- tryCatch(
    harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2),
    error = function(e) { cat("  harmonise_data ERROR:", conditionMessage(e), "\n"); data.frame() }
  )
  cat(sprintf("  Harmonised: %d rows (%d mr_keep)\n\n",
              nrow(harm_dat), if (nrow(harm_dat) > 0) sum(harm_dat$mr_keep) else 0))

  rds_path <- file.path(out_dir, paste0("harmonised_protein_", cancer_name, ".rds"))
  saveRDS(harm_dat, rds_path)

  harm_log[[cancer_name]] <- data.table(
    exposure="All_Proteins", outcome=cancer_name,
    type="Protein->Cancer", n_before=nrow(exp_dat), n_after=nrow(harm_dat),
    n_dropped=nrow(exp_dat) - nrow(harm_dat))
}

# Update harmonisation log
log_file <- file.path(out_dir, "harmonisation_log.csv")
new_log  <- rbindlist(harm_log, fill = TRUE)
if (file.exists(log_file)) {
  old_log <- fread(log_file)
  old_log <- old_log[type != "Protein->Cancer"]
  full_log <- rbind(old_log, new_log, fill = TRUE)
} else {
  full_log <- new_log
}
fwrite(full_log, log_file)
cat("Harmonisation log updated:", log_file, "\n")
cat("\nDone. Next: run scripts/05_protein_cancer_mr.R\n")
sessionInfo()
