#!/usr/bin/env Rscript

# Script 05: Batch MR (Protein -> Cancer)
#
# QC FILTERS APPLIED (updated 2026-05-05):
#   1. F-statistic > 10 per instrument (F = beta^2 / se^2)
#      All pilot pQTL instruments already pass (F range: 45–234).
#      Filter retained for correctness when the full UKB-PPP panel is loaded.
#   2. Steiger filtering NOT applied here. Steiger filtering requires proper
#      liability-scale r² for binary cancer outcomes (needs n_cases, n_controls,
#      and prevalence). Without these parameters, TwoSampleMR::steiger_filtering()
#      systematically drops all protein instruments — a false negative.
#      Instead, sensitivity is assessed via MR-Egger and weighted median.
#   3. FDR correction applied within each cancer outcome separately.
#   4. Methods: Wald ratio (1 SNP), IVW + MR-Egger + Weighted median (>1 SNP).

set.seed(42)

if (!require("here",       quietly = TRUE)) install.packages("here",       repos = "https://cloud.r-project.org")
if (!require("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")
if (!require("tidyverse",  quietly = TRUE)) install.packages("tidyverse",  repos = "https://cloud.r-project.org")
if (!require("TwoSampleMR",quietly = TRUE)) install.packages("TwoSampleMR",repos = "https://cloud.r-project.org")

library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(here)

project_dir <- "."
in_dir  <- file.path(project_dir, "data", "harmonised")
out_dir <- file.path(project_dir, "results", "phase2_protein_cancer")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

F_THRESHOLD <- 10

log_file <- file.path(in_dir, "harmonisation_log.csv")
if (file.exists(log_file)) {
  harm_log     <- fread(log_file)
  target_files <- harm_log[type == "Protein->Cancer" & n_after > 0,
                            file.path(in_dir, paste0("harmonised_protein_", outcome, ".rds"))]
  harm_files   <- target_files[file.exists(target_files)]
} else {
  harm_files <- list.files(in_dir, pattern = "^harmonised_protein_.*\\.rds$", full.names = TRUE)
}

mr_results_all <- data.table()
qc_log         <- data.table()

cat("Starting Batch MR for Protein -> Cancer (with F>10 filter; no Steiger for binary outcomes)...\n")

for (f in harm_files) {
  cat("Processing:", basename(f), "\n")
  dat <- readRDS(f)
  if (nrow(dat) == 0) { cat("  Empty — skipping.\n"); next }
  dat <- dat[dat$mr_keep == TRUE & !is.na(dat$mr_keep), ]
  if (nrow(dat) == 0) { cat("  No mr_keep instruments — skipping.\n"); next }

  n_before <- nrow(dat)

  # ------------------------------------------------------------------
  # Filter 1: F-statistic > 10
  # ------------------------------------------------------------------
  dat <- dat[!is.na(dat$beta.exposure) & !is.na(dat$se.exposure) & dat$se.exposure > 0, ]
  dat$F_stat <- dat$beta.exposure^2 / dat$se.exposure^2
  dat        <- dat[dat$F_stat > F_THRESHOLD, ]

  if (nrow(dat) == 0) {
    cat(sprintf("  All instruments failed F>%d — skipping.\n", F_THRESHOLD))
    next
  }

  n_after_f <- nrow(dat)
  cat(sprintf("  F-filter: %d → %d instruments\n", n_before, n_after_f))

  qc_log <- rbind(qc_log, data.table(
    file = basename(f), n_before, n_after_f
  ), fill = TRUE)

  # ------------------------------------------------------------------
  # MR
  # ------------------------------------------------------------------
  res <- tryCatch(
    mr(dat, method_list = c("mr_wald_ratio", "mr_ivw",
                            "mr_egger_regression", "mr_weighted_median")),
    error = function(e) { cat("  MR failed:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(res) && nrow(res) > 0) {
    qc_by_exposure <- as.data.table(dat)[, .(
      F_stat_mean = mean(F_stat, na.rm = TRUE),
      n_instruments = .N
    ), by = .(id.exposure, exposure)]
    res <- merge(as.data.table(res), qc_by_exposure,
                 by = c("id.exposure", "exposure"), all.x = TRUE, sort = FALSE)
    res$single_snp_only <- res$nsnp == 1
    res <- generate_odds_ratios(res)
    mr_results_all <- rbind(mr_results_all, res, fill = TRUE)
  }
}

# ------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------
fwrite(qc_log, file.path(out_dir, "protein_cancer_qc_log.csv"))

if (nrow(mr_results_all) > 0) {
  mr_results_all[, fdr := p.adjust(pval, method = "BH"), by = outcome]
  fwrite(mr_results_all, file.path(out_dir, "protein_cancer_mr_results_full.csv"))
  sig <- mr_results_all[fdr < 0.05]
  fwrite(sig, file.path(out_dir, "protein_cancer_mr_results_significant.csv"))
  cat(sprintf("\nCompleted. %d total MR rows; %d FDR<0.05\n", nrow(mr_results_all), nrow(sig)))
  cat("Saved to:", out_dir, "\n")
} else {
  cat("No valid MR results generated.\n")
  fwrite(mr_results_all, file.path(out_dir, "protein_cancer_mr_results_full.csv"))
  fwrite(mr_results_all, file.path(out_dir, "protein_cancer_mr_results_significant.csv"))
}

sessionInfo()
