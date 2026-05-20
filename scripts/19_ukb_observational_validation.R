#!/usr/bin/env Rscript
# Script 19: Observational case-control protein comparison in UKB-PPP
#
# Uses:
#   olink_mat.rds        — 52,995 UKB participants × 2,918 proteins (NPX, log2 scale)
#   cancer_reg_all.rds   — UKB cancer registry (eid, ICD prefix)
#
# Design:
#   Cases   = participants with cancer in registry (C50=Breast, C54/C55=EC, C56=OC)
#   Controls = participants with NO cancer of any type in registry
#   Test:  linear model: protein_level ~ case_status (unadjusted; no covariates in scope)
#
# Key check: does the observational direction match the MR causal direction?
#   MR predicts: higher protein → higher cancer risk → cases should have HIGHER protein
#   MR predicts: higher protein → lower cancer risk  → cases should have LOWER protein
#
# Outputs:
#   results/observational/ukb_observational_results.csv
#   results/observational/ukb_triangulation_summary.csv

set.seed(42)
suppressPackageStartupMessages(library(data.table))

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
out_dir     <- file.path(project_dir, "results", "observational")
dir.create(out_dir, showWarnings = FALSE)

# ── MR results for hit proteins (direction reference)
mr_ref <- data.table(
  protein  = c("SNX15","EFNA1","FGF5","UMOD","SWAP70","ATRAID","TNFRSF6B",
               "ITIH3","KLB","PM20D1","TSPAN8","FGFR4","IL34","APOE",
               "CGREF1","INHBB","ABO"),
  cancer   = c(rep("Breast",16), "Endometrial"),
  mr_beta  = c(-0.0866, 0.1272,-0.0413, 0.0242, 0.0434,-0.0449, 0.0352,
               -0.0416, 0.0144, 0.0126,-0.0348,-0.0122, 0.0238, 0.0173,
                0.0265, 0.0253, 0.0474),
  mr_or    = c(0.917, 1.136, 0.960, 1.024, 1.044, 0.956, 1.036,
               0.959, 1.015, 1.013, 0.966, 0.988, 1.024, 1.018,
               1.027, 1.026, 1.049),
  mr_fdr   = c(8.07e-13,1.29e-11,1.29e-11,6.90e-05,1.21e-03,1.67e-03,
               1.81e-03,3.25e-03,3.25e-03,4.69e-03,5.99e-03,1.62e-02,
               1.63e-02,2.17e-02,2.66e-02,3.06e-02,3.06e-03)
)
# MR-predicted direction: positive beta = higher protein → more cancer → cases have MORE
mr_ref[, mr_pred_direction := fifelse(mr_beta > 0, "cases_higher", "cases_lower")]

# ── Load data
cat("Loading UKB-PPP Olink matrix (52,995 × 2,918)...\n")
olink <- readRDS(file.path(project_dir, "olink_mat.rds"))
olink_dt <- as.data.table(olink, keep.rownames = "eid")
olink_dt[, eid := as.integer(eid)]
cat(sprintf("  %d participants, %d proteins\n", nrow(olink_dt), ncol(olink_dt) - 1))

cat("Loading cancer registry...\n")
reg <- as.data.table(readRDS(file.path(project_dir, "cancer_reg_all.rds")))
cat(sprintf("  %d registry records\n", nrow(reg)))

# ── Define case/control groups
breast_eids <- unique(reg$eid[grepl("^C50", reg$prefix)])
ec_eids     <- unique(reg$eid[grepl("^C54|^C55", reg$prefix)])
oc_eids     <- unique(reg$eid[grepl("^C56", reg$prefix)])
any_cancer  <- unique(reg$eid)
control_eids <- olink_dt$eid[!olink_dt$eid %in% any_cancer]

cat(sprintf("\nBreast cases with Olink:      %d\n", sum(olink_dt$eid %in% breast_eids)))
cat(sprintf("EC cases with Olink:          %d\n", sum(olink_dt$eid %in% ec_eids)))
cat(sprintf("OC cases with Olink:          %d\n", sum(olink_dt$eid %in% oc_eids)))
cat(sprintf("Controls (no cancer, Olink):  %d\n", length(control_eids)))

# ── Run case-control comparison for each protein × cancer
run_comparison <- function(case_eids, cancer_label, proteins) {
  proteins_l <- tolower(proteins)
  cases    <- olink_dt[eid %in% case_eids,    .SD, .SDcols = c("eid", proteins_l)]
  controls <- olink_dt[eid %in% control_eids, .SD, .SDcols = c("eid", proteins_l)]

  n_cases    <- nrow(cases)
  n_controls <- nrow(controls)
  cat(sprintf("  %s: %d cases vs %d controls\n", cancer_label, n_cases, n_controls))

  results <- lapply(proteins, function(prot) {
    prot_l <- tolower(prot)
    if (!prot_l %in% names(cases)) return(NULL)

    y_case <- cases[[prot_l]]
    y_ctrl <- controls[[prot_l]]

    # Remove NAs
    y_case <- y_case[!is.na(y_case)]
    y_ctrl <- y_ctrl[!is.na(y_ctrl)]
    if (length(y_case) < 5 || length(y_ctrl) < 5) return(NULL)

    tt <- t.test(y_case, y_ctrl)
    lm_dat <- data.frame(
      y      = c(y_case, y_ctrl),
      status = c(rep(1L, length(y_case)), rep(0L, length(y_ctrl)))
    )
    lm_fit <- lm(y ~ status, data = lm_dat)
    lm_sum <- summary(lm_fit)$coefficients

    data.table(
      protein          = prot,
      cancer           = cancer_label,
      n_cases          = length(y_case),
      n_controls       = length(y_ctrl),
      mean_cases       = round(mean(y_case), 4),
      mean_controls    = round(mean(y_ctrl), 4),
      delta_npx        = round(mean(y_case) - mean(y_ctrl), 4),
      beta_lm          = round(lm_sum["status","Estimate"], 5),
      se_lm            = round(lm_sum["status","Std. Error"], 5),
      pval_lm          = lm_sum["status","Pr(>|t|)"],
      t_stat           = round(tt$statistic, 3),
      obs_direction    = fifelse(mean(y_case) > mean(y_ctrl), "cases_higher", "cases_lower")
    )
  })
  rbindlist(results, fill = TRUE)
}

# ── Run for each cancer
proteins_breast <- mr_ref$protein[mr_ref$cancer == "Breast"]
proteins_ec     <- mr_ref$protein[mr_ref$cancer == "Endometrial"]

cat("\n=== Breast cancer comparison ===\n")
res_breast <- run_comparison(breast_eids, "Breast", proteins_breast)

cat("\n=== Endometrial cancer comparison ===\n")
res_ec <- run_comparison(ec_eids, "Endometrial", proteins_ec)

cat("\n=== Ovarian cancer comparison (bonus: test breast proteins) ===\n")
res_oc <- run_comparison(oc_eids, "Ovarian", proteins_breast)

all_obs <- rbindlist(list(res_breast, res_ec, res_oc), fill = TRUE)

# ── FDR correction within cancer
all_obs[, fdr_obs := p.adjust(pval_lm, method = "BH"), by = cancer]

# ── Merge with MR reference and check triangulation
all_obs <- merge(all_obs, mr_ref[, .(protein, cancer, mr_beta, mr_or, mr_fdr, mr_pred_direction)],
                 by = c("protein","cancer"), all.x = TRUE)

all_obs[, triangulates := fifelse(
  !is.na(obs_direction) & !is.na(mr_pred_direction),
  obs_direction == mr_pred_direction, NA)]

# ── Save
fwrite(all_obs, file.path(out_dir, "ukb_observational_results.csv"))

# ── Summary table
cat("\n=== TRIANGULATION SUMMARY ===\n")
cat("(Do observational directions match MR predictions?)\n\n")

triang <- all_obs[cancer == "Breast"][order(mr_fdr)]
print(triang[, .(protein, mr_or=round(mr_or,4), mr_pred_direction,
                 delta_npx, obs_direction, pval_lm=signif(pval_lm,3),
                 fdr_obs=signif(fdr_obs,3), triangulates)])

cat("\n--- ABO → Endometrial ---\n")
print(all_obs[cancer == "Endometrial",
              .(protein, mr_or, mr_pred_direction, delta_npx,
                obs_direction, pval_lm=signif(pval_lm,3), triangulates)])

# Overall triangulation score
n_total     <- sum(!is.na(all_obs[cancer!="Ovarian"]$triangulates))
n_agree     <- sum(all_obs[cancer!="Ovarian"]$triangulates, na.rm=TRUE)
cat(sprintf("\nDirectional agreement (MR vs observational): %d / %d (%.0f%%)\n",
            n_agree, n_total, 100*n_agree/n_total))

# ── Compact summary file
summary_dt <- all_obs[cancer != "Ovarian",
  .(protein, cancer, mr_or=round(mr_or,4), mr_pred_direction,
    delta_npx, obs_direction, pval_obs=signif(pval_lm,3),
    fdr_obs=signif(fdr_obs,3), triangulates, n_cases, n_controls)]
setorder(summary_dt, cancer, -triangulates, pval_obs)
fwrite(summary_dt, file.path(out_dir, "ukb_triangulation_summary.csv"))

cat("\nOutputs:\n")
cat("  ", file.path(out_dir, "ukb_observational_results.csv"), "\n")
cat("  ", file.path(out_dir, "ukb_triangulation_summary.csv"), "\n")
cat("Done.\n")
