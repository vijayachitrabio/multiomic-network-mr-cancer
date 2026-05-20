#!/usr/bin/env Rscript

# Script 09: Batch MR (Protein -> Metabolite)
#
# Primary estimator:
#   - MR-RAPS for multi-SNP protein instruments, to reduce bias from sample
#     overlap between UKB-PPP pQTLs and UKB metabolite GWAS. The pilot pQTL
#     set has some 2-SNP proteins; these use the simple MR-RAPS model because
#     overdispersion estimation is unstable with only two variants.
#   - Wald ratio for single-SNP instruments.
#
# Multiple testing:
#   - BH FDR is applied within each protein across the 56 metabolites.

set.seed(42)

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
local_lib <- file.path(project_dir, "r_libs")
if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE)
.libPaths(c(local_lib, .libPaths()))

required_pkgs <- c("here", "data.table", "tidyverse", "TwoSampleMR", "mr.raps")
repos <- c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org")
for (p in required_pkgs) {
  if (!require(p, quietly = TRUE, character.only = TRUE)) {
    install.packages(p, lib = local_lib, repos = repos)
  }
  library(p, character.only = TRUE)
}

in_dir  <- file.path(project_dir, "data", "harmonised")
out_dir <- file.path(project_dir, "results", "phase3_protein_metabolite")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

pick_numeric <- function(x, candidates) {
  for (nm in candidates) {
    if (!is.null(x[[nm]]) && length(x[[nm]]) > 0 && is.numeric(x[[nm]][1])) {
      return(as.numeric(x[[nm]][1]))
    }
  }
  NA_real_
}

run_primary_mr <- function(dat) {
  dat <- dat[dat$mr_keep == TRUE, ]
  dat <- dat[!is.na(dat$beta.exposure) & !is.na(dat$beta.outcome) &
               !is.na(dat$se.exposure) & !is.na(dat$se.outcome) &
               dat$se.exposure > 0 & dat$se.outcome > 0, ]
  if (nrow(dat) == 0) return(NULL)

  base <- data.table(
    id.exposure = dat$id.exposure[1],
    id.outcome  = dat$id.outcome[1],
    outcome     = dat$outcome[1],
    exposure    = dat$exposure[1],
    nsnp        = nrow(dat),
    F_stat_mean = mean(dat$beta.exposure^2 / dat$se.exposure^2, na.rm = TRUE),
    n_instruments = nrow(dat),
    single_snp_only = nrow(dat) == 1
  )

  if (nrow(dat) == 1) {
    res <- mr(dat, method_list = "mr_wald_ratio")
    if (nrow(res) == 0) return(NULL)
    out <- as.data.table(res)[, .(id.exposure, id.outcome, outcome, exposure,
                                  method, nsnp, b, se, pval)]
    out <- merge(out, base[, .(id.exposure, id.outcome, exposure, F_stat_mean,
                               n_instruments, single_snp_only)],
                 by = c("id.exposure", "id.outcome", "exposure"),
                 all.x = TRUE, sort = FALSE)
    return(out)
  }

  use_overdispersion <- nrow(dat) >= 3
  raps_dat <- data.frame(
    beta.exposure = dat$beta.exposure,
    beta.outcome  = dat$beta.outcome,
    se.exposure   = dat$se.exposure,
    se.outcome    = dat$se.outcome
  )

  fit <- tryCatch(
    suppressWarnings(suppressMessages(mr.raps::mr.raps(
      raps_dat,
      diagnostics = FALSE,
      over.dispersion = use_overdispersion,
      loss.function = "tukey"
    ))),
    error = function(e) {
      message("  MR-RAPS failed for ", dat$exposure[1], " -> ", dat$outcome[1],
              ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) return(NULL)

  b <- pick_numeric(fit, c("beta.hat", "beta_hat", "beta"))
  se <- pick_numeric(fit, c("beta.se", "beta_se", "se"))
  pval <- pick_numeric(fit, c("beta.p.value", "beta_p.value", "p.value", "pval"))
  if (is.na(pval) && !is.na(b) && !is.na(se) && se > 0) {
    pval <- 2 * pnorm(-abs(b / se))
  }
  if (is.na(b) || is.na(se) || is.na(pval)) return(NULL)

  cbind(
    base[, .(id.exposure, id.outcome, outcome, exposure)],
    data.table(
      method = if (use_overdispersion) "MR-RAPS overdispersed" else "MR-RAPS simple",
      nsnp = nrow(dat), b = b, se = se, pval = pval
    ),
    base[, .(F_stat_mean, n_instruments, single_snp_only)]
  )
}

# Find only harmonised Protein -> Metabolite files from the harmonisation log.
log_file <- file.path(in_dir, "harmonisation_log.csv")
if (file.exists(log_file)) {
  harm_log <- fread(log_file)
  target_files <- harm_log[type == "Protein->Metabolite" & n_after > 0,
                           file.path(in_dir, paste0("harmonised_protein_", outcome, ".rds"))]
  harm_files <- target_files[file.exists(target_files)]
} else {
  harm_files <- list.files(in_dir, pattern = "^harmonised_protein_.*\\.rds$", full.names = TRUE)
}

mr_results_all <- data.table()
qc_log <- data.table()

cat("Starting Batch MR for Protein -> Metabolite (MR-RAPS primary)...\n")
for (f in harm_files) {
  cat("Processing:", basename(f), "\n")
  dat <- readRDS(f)
  if (nrow(dat) == 0) next

  dat <- as.data.table(dat)
  split_dat <- split(dat, dat$exposure)

  for (exposure_name in names(split_dat)) {
    sub <- split_dat[[exposure_name]]
    res <- run_primary_mr(sub)

    qc_log <- rbind(qc_log, data.table(
      file = basename(f),
      outcome = sub$outcome[1],
      exposure = exposure_name,
      n_harmonised = nrow(sub),
      n_mr_keep = sum(sub$mr_keep == TRUE, na.rm = TRUE),
      method = if (!is.null(res) && nrow(res) > 0) res$method[1] else NA_character_,
      mr_success = !is.null(res) && nrow(res) > 0
    ), fill = TRUE)

    if (!is.null(res) && nrow(res) > 0) {
      mr_results_all <- rbind(mr_results_all, res, fill = TRUE)
    }
  }
}

fwrite(qc_log, file.path(out_dir, "protein_metabolite_qc_log.csv"))

if (nrow(mr_results_all) > 0) {
  mr_results_all[, fdr := p.adjust(pval, method = "BH"), by = exposure]

  fwrite(mr_results_all, file.path(out_dir, "protein_metabolite_mr_results_full.csv"))
  sig_results <- mr_results_all[fdr < 0.05]
  fwrite(sig_results, file.path(out_dir, "protein_metabolite_mr_results_significant.csv"))

  cat(sprintf("Completed. %d total MR rows; %d FDR<0.05\n",
              nrow(mr_results_all), nrow(sig_results)))
  cat("Saved results to", out_dir, "\n")
} else {
  cat("No valid MR results generated.\n")
  fwrite(mr_results_all, file.path(out_dir, "protein_metabolite_mr_results_full.csv"))
  fwrite(mr_results_all, file.path(out_dir, "protein_metabolite_mr_results_significant.csv"))
}

sessionInfo()
