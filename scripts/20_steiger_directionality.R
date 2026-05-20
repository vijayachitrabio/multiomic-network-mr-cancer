#!/usr/bin/env Rscript
# Script 20: Steiger directionality tests for all Phase 2 protein-cancer hits
#
# Purpose:
#   Verify that genetic variants used as pQTL instruments explain more variance
#   in the protein (exposure) than in the cancer outcome.
#   Steiger filtering: keep only SNPs where r2_exp > r2_out → protects against
#   reverse-causation bias (cancer liability affecting measured protein levels).
#
# Method:
#   Uses TwoSampleMR::steiger_filtering() on harmonised data
#   Requires prevalence (s) and sample sizes for R2 estimation.
#
# For cancer outcomes (binary, case-control):
#   r2 from GWAS beta is estimated on the liability scale using Burgess formula
#   s = n_cases / n_total
#
# For protein exposure:
#   r2 from pQTL beta estimated with N = 619 (FinnGen R10 Olink)
#
# Outputs:
#   results/sensitivity/steiger_directionality_results.csv
#   results/sensitivity/steiger_mr_filtered_results.csv  — MR after removing reversed SNPs

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
out_dir     <- file.path(project_dir, "results", "sensitivity")
dir.create(out_dir, showWarnings = FALSE)

# ── Phase 2 hits
breast_hits <- c("SNX15","EFNA1","FGF5","UMOD","SWAP70","ATRAID","TNFRSF6B",
                 "ITIH3","KLB","PM20D1","TSPAN8","FGFR4","IL34","APOE",
                 "CGREF1","INHBB")
ec_hits <- c("ABO")

all_hits <- data.frame(
  protein = c(breast_hits, ec_hits),
  cancer  = c(rep("Breast", length(breast_hits)), "Endometrial"),
  stringsAsFactors = FALSE
)

# ── Cancer metadata (sample sizes for r2 estimation)
cancer_meta <- data.frame(
  outcome    = c("Breast_GCST90018757", "Endometrial_GCST006464", "Ovarian_GCST90016665"),
  n_cases    = c(122977, 12906, 25509),
  n_controls = c(105974, 108979, 40138),
  stringsAsFactors = FALSE
)
cancer_meta$n_total <- cancer_meta$n_cases + cancer_meta$n_controls
cancer_meta$s       <- cancer_meta$n_cases / cancer_meta$n_total

outcome_map <- c(
  Breast       = "Breast_GCST90018757",
  Endometrial  = "Endometrial_GCST006464",
  Ovarian      = "Ovarian_GCST90016665"
)

# ── FinnGen Olink pQTL N
N_pqtl <- 619L

# ── Load harmonised data
harm_dir <- file.path(project_dir, "data", "harmonised")
rds_files <- list.files(harm_dir, pattern = "^harmonised_protein_.*\\.rds$", full.names = TRUE)

# Keep only main cancer outcomes (not ER subtypes)
main_outcomes <- c("Breast_GCST90018757", "Endometrial_GCST006464", "Ovarian_GCST90016665")
rds_files <- rds_files[!grepl("BreastER", rds_files)]

cat(sprintf("Found %d harmonised RDS files\n", length(rds_files)))

steiger_all    <- list()
mr_filtered_all <- list()

for (rds_path in rds_files) {
  harm_all <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(harm_all) || nrow(harm_all) == 0) next

  outcome_id <- unique(harm_all$id.outcome)[1]
  if (!outcome_id %in% main_outcomes) next

  meta <- cancer_meta[cancer_meta$outcome == outcome_id, ]
  if (nrow(meta) == 0) next

  # Only proteins in our hit list for this outcome
  cancer_label <- names(outcome_map)[outcome_map == outcome_id]
  hits_here    <- all_hits$protein[all_hits$cancer == cancer_label]
  prots_here   <- intersect(hits_here, unique(harm_all$exposure))

  cat(sprintf("\n=== %s (%s) ===\n", outcome_id, cancer_label))

  for (prot in prots_here) {
    harm <- harm_all[harm_all$exposure == prot & harm_all$mr_keep, ]
    if (nrow(harm) == 0) next

    # Add Steiger inputs: ncase/ncontrol for outcome, N for exposure
    harm$units.outcome    <- "log odds"
    harm$units.exposure   <- "SD"
    harm$ncase.outcome    <- meta$n_cases
    harm$ncontrol.outcome <- meta$n_controls
    harm$samplesize.outcome <- meta$n_total
    harm$samplesize.exposure <- N_pqtl
    harm$prevalence.outcome  <- meta$s

    stei <- tryCatch({
      res <- steiger_filtering(harm)
      res
    }, error = function(e) {
      cat(sprintf("  %s: steiger_filtering error: %s\n", prot, conditionMessage(e)))
      NULL
    })
    if (is.null(stei)) next

    # Summarise per protein
    n_snp       <- nrow(stei)
    n_fwd       <- sum(stei$steiger_dir, na.rm = TRUE)
    n_rev       <- n_snp - n_fwd
    cat(sprintf("  %s: %d SNPs — %d forward, %d reversed\n", prot, n_snp, n_fwd, n_rev))

    stei_dt <- as.data.table(stei)
    stei_dt[, protein  := prot]
    stei_dt[, cancer   := cancer_label]
    stei_dt[, outcome  := outcome_id]
    steiger_all[[paste(prot, cancer_label)]] <- stei_dt

    # Run MR on Steiger-filtered SNPs (forward direction only)
    harm_fwd <- stei[stei$steiger_dir & stei$mr_keep, ]
    if (nrow(harm_fwd) == 0) {
      cat(sprintf("  %s: no SNPs remain after Steiger filter — skipping filtered MR\n", prot))
      next
    }

    method <- if (nrow(harm_fwd) >= 2) "mr_ivw" else "mr_wald_ratio"
    res_fwd <- tryCatch(mr(harm_fwd, method_list = method), error = function(e) NULL)
    if (is.null(res_fwd) || nrow(res_fwd) == 0) next

    res_fwd$protein      <- prot
    res_fwd$cancer       <- cancer_label
    res_fwd$n_snp_total  <- n_snp
    res_fwd$n_snp_fwd    <- n_fwd
    res_fwd$n_snp_rev    <- n_rev
    res_fwd$or           <- exp(res_fwd$b)
    res_fwd$or_lci95     <- exp(res_fwd$b - 1.96 * res_fwd$se)
    res_fwd$or_uci95     <- exp(res_fwd$b + 1.96 * res_fwd$se)
    mr_filtered_all[[paste(prot, cancer_label)]] <- as.data.table(res_fwd)
  }
}

# ── Save
stei_out <- rbindlist(steiger_all, fill = TRUE)
mr_filt_out <- rbindlist(mr_filtered_all, fill = TRUE)

fwrite(stei_out, file.path(out_dir, "steiger_directionality_results.csv"))
fwrite(mr_filt_out, file.path(out_dir, "steiger_mr_filtered_results.csv"))

cat("\n=== STEIGER SUMMARY ===\n")
cat("Comparing full IVW vs Steiger-filtered MR:\n\n")

# Load original Phase 2 results for comparison
phase2 <- fread(file.path(project_dir, "results", "phase2_protein_cancer",
                           "protein_cancer_mr_results_full.csv"))
phase2_hits <- phase2[fdr < 0.05 & exposure %in% c(breast_hits, ec_hits)]
phase2_hits[, cancer := fifelse(grepl("Breast", outcome), "Breast",
                         fifelse(grepl("Endometrial", outcome), "Endometrial", "Ovarian"))]

comp <- merge(
  phase2_hits[, .(protein=exposure, cancer, b_full=round(b,4), se_full=round(se,5),
                   pval_full=signif(pval,3), or_full=round(or,4), nsnp_full=nsnp)],
  mr_filt_out[, .(protein, cancer, b_filt=round(b,4), se_filt=round(se,5),
                   pval_filt=signif(pval,3), or_filt=round(or,4),
                   nsnp_filt=nsnp, n_snp_rev)],
  by = c("protein","cancer"), all.x = TRUE
)
comp[, direction_consistent := fifelse(!is.na(b_full) & !is.na(b_filt),
                                        sign(b_full) == sign(b_filt), NA)]
comp[, or_change := round(abs(or_filt - or_full), 4)]
setorder(comp, cancer, pval_full)

print(comp[, .(protein, cancer, or_full, or_filt, n_snp_rev, direction_consistent, or_change,
               pval_full, pval_filt)])

n_consistent <- sum(comp$direction_consistent, na.rm = TRUE)
n_total_comp <- sum(!is.na(comp$direction_consistent))
cat(sprintf("\nDirection consistent after Steiger filter: %d / %d\n", n_consistent, n_total_comp))

fwrite(comp, file.path(out_dir, "steiger_comparison_table.csv"))

cat("\nOutputs:\n")
cat("  ", file.path(out_dir, "steiger_directionality_results.csv"), "\n")
cat("  ", file.path(out_dir, "steiger_mr_filtered_results.csv"), "\n")
cat("  ", file.path(out_dir, "steiger_comparison_table.csv"), "\n")
cat("Done.\n")
