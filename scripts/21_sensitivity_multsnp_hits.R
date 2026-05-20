#!/usr/bin/env Rscript
# Script 17: Sensitivity analyses for multi-SNP Phase 2 protein→cancer FDR hits
#
# Targets (nsnp >= 2, FDR < 0.05):
#   ABO   → Endometrial (2 SNPs)
#   KLB   → Breast      (2 SNPs)
#   PM20D1 → Breast     (2 SNPs)
#   IL34  → Breast      (2 SNPs)
#
# Methods run per pair:
#   - IVW (re-run, matches Phase 2 primary)
#   - Weighted median (robust to ≤50% invalid instruments)
#   - MR-Egger (intercept = pleiotropy test; slope = directional sensitivity)
#   - Leave-one-out (which SNP drives the signal?)
#   - Single-SNP forest plot data
#
# Outputs:
#   results/sensitivity/sensitivity_multsnp_results.csv
#   results/sensitivity/sensitivity_egger_intercepts.csv
#   results/sensitivity/sensitivity_loo_results.csv
#   results/sensitivity/sensitivity_singlesnp_results.csv

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

project_dir  <- "."
harm_dir     <- file.path(project_dir, "data", "harmonised")
out_dir      <- file.path(project_dir, "results", "sensitivity")
dir.create(out_dir, showWarnings = FALSE)

# Multi-SNP FDR<0.05 pairs
targets <- data.frame(
  protein  = c("ABO",   "KLB",                "PM20D1",             "IL34"),
  cancer   = c("Endometrial_GCST006464", "Breast_GCST90018757",
               "Breast_GCST90018757",   "Breast_GCST90018757"),
  harm_rds = c("harmonised_protein_Endometrial_GCST006464.rds",
               "harmonised_protein_Breast_GCST90018757.rds",
               "harmonised_protein_Breast_GCST90018757.rds",
               "harmonised_protein_Breast_GCST90018757.rds"),
  stringsAsFactors = FALSE
)

all_mr   <- list()
all_egg  <- list()
all_loo  <- list()
all_ss   <- list()

for (i in seq_len(nrow(targets))) {
  prot   <- targets$protein[i]
  cancer <- targets$cancer[i]
  rds    <- file.path(harm_dir, targets$harm_rds[i])

  cat(sprintf("\n=== %s -> %s ===\n", prot, cancer))

  harm_all <- readRDS(rds)
  harm <- harm_all[harm_all$exposure == prot & harm_all$mr_keep, ]
  nsnp <- nrow(harm)
  cat(sprintf("  Instruments: %d\n", nsnp))

  if (nsnp < 2) {
    cat("  Skipping — fewer than 2 SNPs\n")
    next
  }

  # ---- MR methods: IVW always; WM and Egger only if nsnp>=3 ----
  methods_use <- "mr_ivw"
  if (nsnp >= 3) methods_use <- c(methods_use, "mr_weighted_median", "mr_egger_regression")
  mr_res <- mr(harm, method_list = methods_use)
  mr_res$exposure <- prot
  mr_res$outcome  <- cancer
  mr_res$or       <- exp(mr_res$b)
  mr_res$or_lci95 <- exp(mr_res$b - 1.96 * mr_res$se)
  mr_res$or_uci95 <- exp(mr_res$b + 1.96 * mr_res$se)
  all_mr[[paste(prot, cancer, sep = "_")]] <- mr_res
  cat(sprintf("  MR results (%d methods):\n", nrow(mr_res)))
  print(mr_res[, c("method", "nsnp", "b", "se", "pval", "or")])

  # ---- Egger intercept (pleiotropy): only if nsnp>=3 ----
  if (nsnp >= 3) {
    egg <- mr_pleiotropy_test(harm)
    egg$exposure <- prot; egg$outcome <- cancer
    all_egg[[paste(prot, cancer, sep = "_")]] <- egg
    cat(sprintf("  MR-Egger intercept: %.4f (p=%.3f)\n", egg$egger_intercept, egg$pval))
  } else {
    cat("  MR-Egger: not applicable (need >=3 SNPs)\n")
    all_egg[[paste(prot, cancer, sep = "_")]] <- data.frame(
      exposure=prot, outcome=cancer, egger_intercept=NA_real_,
      se=NA_real_, pval=NA_real_, note="nsnp<3")
  }

  # ---- Per-SNP Wald ratios — directional consistency ----
  cat("  Per-SNP Wald ratios (directional consistency check):\n")
  wald_list <- lapply(seq_len(nrow(harm)), function(j) {
    h1 <- harm[j, ]
    wr <- mr(h1, method_list = "mr_wald_ratio")
    data.frame(SNP = h1$SNP, b = wr$b, se = wr$se, pval = wr$pval,
               or = exp(wr$b), exposure = prot, outcome = cancer)
  })
  wald_dt <- do.call(rbind, wald_list)
  wald_dt$same_direction_as_ivw <- sign(wald_dt$b) == sign(mr_res$b[mr_res$method == "Inverse variance weighted"])
  print(wald_dt[, c("SNP", "b", "se", "pval", "or", "same_direction_as_ivw")])
  all_ss[[paste(prot, cancer, sep = "_")]] <- wald_dt

  # ---- Leave-one-out ----
  loo <- mr_leaveoneout(harm)
  loo$exposure <- prot; loo$outcome  <- cancer
  all_loo[[paste(prot, cancer, sep = "_")]] <- loo
  cat(sprintf("  Leave-one-out: %d rows\n", nrow(loo)))
  print(loo[, c("SNP", "b", "se", "p")])

  # ---- Single-SNP ----
  ss <- mr_singlesnp(harm)
  ss$exposure <- prot
  ss$outcome  <- cancer
  all_ss[[paste(prot, cancer, sep = "_")]] <- ss
}

# ---- Save ----
mr_dt  <- rbindlist(all_mr,  fill = TRUE)
egg_dt <- rbindlist(all_egg, fill = TRUE)
loo_dt <- rbindlist(all_loo, fill = TRUE)
ss_dt  <- rbindlist(all_ss,  fill = TRUE)

fwrite(mr_dt,  file.path(out_dir, "sensitivity_multsnp_results.csv"))
fwrite(egg_dt, file.path(out_dir, "sensitivity_egger_intercepts.csv"))
fwrite(loo_dt, file.path(out_dir, "sensitivity_loo_results.csv"))
fwrite(ss_dt,  file.path(out_dir, "sensitivity_singlesnp_results.csv"))

cat("\n=== SENSITIVITY SUMMARY ===\n")
cat("Methods compared (IVW / Weighted Median / MR-Egger):\n")
print(mr_dt[, .(exposure, outcome, method, nsnp, b=round(b,4), se=round(se,4),
                pval=signif(pval,3), or=round(or,4))])
cat("\nEgger intercepts (H0: no directional pleiotropy):\n")
print(egg_dt[, .(exposure, outcome, egger_intercept=round(egger_intercept,5), pval=round(pval,3))])

cat("\nOutputs written to:", out_dir, "\n")
cat("Done. Next: review sensitivity_multsnp_results.csv\n")
sessionInfo()
