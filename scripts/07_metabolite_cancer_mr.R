#!/usr/bin/env Rscript

# Script 07: Batch MR (Metabolite -> Cancer)
#
# QC FILTERS APPLIED (updated 2026-05-05):
#   1. F-statistic > 10 per instrument (F = beta^2 / se^2)
#      Retained for correctness with both legacy SuSiE and GWAS-based mQTLs.
#   2. Steiger filtering is attempted only when TwoSampleMR can compute
#      non-missing directions. If directions are all missing for binary cancer
#      outcomes, instruments are retained and the QC log marks Steiger as not
#      applied instead of silently pretending a filter occurred.
#   3. FDR correction applied within each cancer outcome separately.

set.seed(42)

if (!require("here",       quietly = TRUE)) install.packages("here",       repos = "https://cloud.r-project.org")
if (!require("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")
if (!require("tidyverse",  quietly = TRUE)) install.packages("tidyverse",  repos = "https://cloud.r-project.org")
if (!require("TwoSampleMR",quietly = TRUE)) install.packages("TwoSampleMR",repos = "https://cloud.r-project.org")

library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(here)

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
in_dir  <- file.path(project_dir, "data", "harmonised")
out_dir <- file.path(project_dir, "results", "phase4_metabolite_cancer")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

F_THRESHOLD <- 10   # minimum instrument F-statistic

log_file <- file.path(in_dir, "harmonisation_log.csv")
if (file.exists(log_file)) {
  harm_log   <- fread(log_file)
  target_files <- harm_log[type == "Metabolite->Cancer" & n_after > 0,
                            file.path(in_dir, paste0("harmonised_metabolite_", outcome, ".rds"))]
  harm_files <- target_files[file.exists(target_files)]
} else {
  harm_files <- list.files(in_dir, pattern = "^harmonised_metabolite_.*\\.rds$", full.names = TRUE)
}

mr_results_all  <- data.table()
qc_log          <- data.table()

cat("Starting Batch MR for Metabolite -> Cancer (with F>10 and conditional Steiger filter)...\n")

for (f in harm_files) {
  cat("Processing:", basename(f), "\n")
  dat <- readRDS(f)
  if (nrow(dat) == 0) next
  dat <- dat[dat$mr_keep == TRUE & !is.na(dat$mr_keep), ]
  if (nrow(dat) == 0) {
    cat("  No mr_keep instruments — skipping.\n")
    next
  }

  n_before <- nrow(dat)

  # ------------------------------------------------------------------
  # Filter 1: F-statistic > 10
  # ------------------------------------------------------------------
  dat <- dat[!is.na(dat$beta.exposure) & !is.na(dat$se.exposure) & dat$se.exposure > 0, ]
  dat$F_stat <- dat$beta.exposure^2 / dat$se.exposure^2
  n_after_f  <- sum(dat$F_stat > F_THRESHOLD, na.rm = TRUE)
  dat        <- dat[dat$F_stat > F_THRESHOLD, ]

  if (nrow(dat) == 0) {
    cat(sprintf("  All %d instruments failed F>%d — skipping.\n", n_before, F_THRESHOLD))
    qc_log <- rbind(qc_log, data.table(
      file = basename(f), n_before, n_after_f = 0, n_after_steiger = 0,
      steiger_applied = FALSE, skipped = TRUE
    ), fill = TRUE)
    next
  }

  cat(sprintf("  F-filter: %d → %d instruments (dropped %d with F≤%d)\n",
              n_before, nrow(dat), n_before - nrow(dat), F_THRESHOLD))

  # ------------------------------------------------------------------
  # Filter 2: Steiger filtering
  # ------------------------------------------------------------------
  steiger_applied <- FALSE
  dat_steiger <- tryCatch({
    s <- steiger_filtering(dat)
    if (!"steiger_dir" %in% names(s) || all(is.na(s$steiger_dir))) {
      cat("  Steiger directions unavailable — retaining F-filtered instruments.\n")
      dat
    } else {
      steiger_applied <<- TRUE
      s[s$steiger_dir == TRUE & !is.na(s$steiger_dir), ]
    }
  }, error = function(e) {
    cat("  Steiger filtering failed:", conditionMessage(e), "— skipping Steiger.\n")
    dat
  })

  n_after_steiger <- nrow(dat_steiger)
  if (steiger_applied) {
    cat(sprintf("  Steiger filter: %d → %d instruments (dropped %d with wrong direction)\n",
                nrow(dat), n_after_steiger, nrow(dat) - n_after_steiger))
  } else {
    cat(sprintf("  Steiger filter: not applied; retained %d F-filtered instruments\n",
                n_after_steiger))
  }

  qc_log <- rbind(qc_log, data.table(
    file = basename(f), n_before, n_after_f, n_after_steiger,
    steiger_applied, skipped = FALSE
  ), fill = TRUE)

  if (nrow(dat_steiger) == 0) {
    cat("  No instruments passed Steiger filter — skipping.\n")
    next
  }

  # ------------------------------------------------------------------
  # MR
  # ------------------------------------------------------------------
  res <- tryCatch(
    mr(dat_steiger, method_list = c("mr_wald_ratio", "mr_ivw",
                                    "mr_egger_regression", "mr_weighted_median")),
    error = function(e) { cat("  MR failed:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(res) && nrow(res) > 0) {
    qc_by_exposure <- as.data.table(dat_steiger)[, .(
      F_stat_mean = mean(F_stat, na.rm = TRUE),
      n_instruments = .N
    ), by = .(id.exposure, exposure)]
    res <- merge(as.data.table(res), qc_by_exposure,
                 by = c("id.exposure", "exposure"), all.x = TRUE, sort = FALSE)
    # Flag single-SNP Wald ratio results: pleiotropy cannot be tested
    res$single_snp_only <- res$nsnp == 1
    res <- generate_odds_ratios(res)
    mr_results_all <- rbind(mr_results_all, res, fill = TRUE)
  }
}

# ------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------
fwrite(qc_log, file.path(out_dir, "metabolite_cancer_qc_log.csv"))

if (nrow(mr_results_all) > 0) {
  mr_results_all[, fdr := p.adjust(pval, method = "BH"), by = outcome]
  fwrite(mr_results_all, file.path(out_dir, "metabolite_cancer_mr_results_full.csv"))
  sig <- mr_results_all[fdr < 0.05]
  fwrite(sig, file.path(out_dir, "metabolite_cancer_mr_results_significant.csv"))
  cat(sprintf("\nCompleted. %d total MR rows; %d FDR<0.05\n", nrow(mr_results_all), nrow(sig)))
  cat("Saved to:", out_dir, "\n")
} else {
  cat("No valid MR results after QC filters.\n")
  fwrite(mr_results_all, file.path(out_dir, "metabolite_cancer_mr_results_full.csv"))
  fwrite(mr_results_all, file.path(out_dir, "metabolite_cancer_mr_results_significant.csv"))
}

sessionInfo()
