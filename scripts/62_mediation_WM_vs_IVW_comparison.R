#!/usr/bin/env Rscript
# Script 62: WM vs IVW comparison for all 6 mediation paths (step-2 only)
# Step-2: protein -> metabolite MR using WM (primary) and IVW (sensitivity)
# Date: 2026-05-25

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
})

PROJ <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"

# ── The 6 mediation paths from Table 2 ─────────────────────────────────────
# path: protein -> metabolite -> Breast cancer
# Step 1 (protein -> cancer) betas are in mediation_mr_results.csv
# Step 2 (protein -> metabolite) we re-run IVW from harmonised data

paths <- list(
  list(protein="IL34",     metabolite="Total_BCAA", cancer="Breast"),
  list(protein="EFNA1",    metabolite="Total_BCAA", cancer="Breast"),
  list(protein="ATRAID",   metabolite="TG_by_PG",   cancer="Breast"),
  list(protein="ITIH3",    metabolite="Gly",         cancer="Breast"),
  list(protein="APOE",     metabolite="Gly",         cancer="Breast"),
  list(protein="TNFRSF6B", metabolite="Total_BCAA", cancer="Breast")
)

# Load existing WM mediation results (step-2 betas already there)
med_res <- read.csv(file.path(PROJ, "results/mediation_mr_results.csv"),
                    stringsAsFactors=FALSE)

# Load step-1 (protein->cancer) IVW betas from STable2
stab2 <- read.csv(file.path(PROJ, "results/tables/STable2_17_FDR_hits_complete.csv"),
                  stringsAsFactors=FALSE)

rows <- list()

for (pth in paths) {
  cat("\n===", pth$protein, "->", pth$metabolite, "-> Breast ===\n")

  # Load protein-metabolite harmonised data
  rds_file <- file.path(PROJ, "data/harmonised",
                        paste0("harmonised_protein_", pth$metabolite, ".rds"))
  if (!file.exists(rds_file)) {
    cat("  SKIP: harmonised file not found:", rds_file, "\n"); next
  }
  dat <- readRDS(rds_file)

  # Subset to this protein
  sub <- dat[grepl(pth$protein, dat$exposure, ignore.case=TRUE) &
               dat$mr_keep == TRUE, ]
  cat("  N SNPs for step-2:", nrow(sub), "\n")

  if (nrow(sub) == 0) { cat("  SKIP: no SNPs\n"); next }

  # Step-2: Run both WM and IVW
  res_all <- mr(sub, method_list=c("mr_ivw", "mr_weighted_median"))

  ivw_row <- res_all[res_all$method == "Inverse variance weighted", ]
  wm_row  <- res_all[res_all$method == "Weighted median", ]

  # If only 1 SNP, both collapse to Wald ratio
  if (nrow(sub) == 1) {
    wald_row <- mr(sub, method_list="mr_wald_ratio")
    ivw_row <- wald_row; wm_row <- wald_row
    cat("  Single SNP: using Wald ratio for both\n")
  }

  # Step-1 beta (protein -> cancer IVW from STable2)
  s1 <- stab2[grepl(pth$protein, stab2$protein, ignore.case=TRUE) &
                 grepl("Breast", stab2$cancer, ignore.case=TRUE), ]
  alpha     <- if (nrow(s1)>0) s1$beta[1]  else NA
  se_alpha  <- if (nrow(s1)>0) s1$se[1]    else NA

  # Existing WM step-2 from mediation results
  wm_existing <- med_res[grepl(pth$protein, med_res$protein, ignore.case=TRUE) &
                            grepl(pth$metabolite, med_res$metabolite, ignore.case=TRUE), ]

  # Step-2 estimates
  beta_wm  <- if (nrow(wm_row)>0)  wm_row$b[1]   else NA
  se_wm    <- if (nrow(wm_row)>0)  wm_row$se[1]  else NA
  p_wm     <- if (nrow(wm_row)>0)  wm_row$pval[1] else NA
  nsnp_wm  <- if (nrow(wm_row)>0)  wm_row$nsnp[1] else nrow(sub)

  beta_ivw <- if (nrow(ivw_row)>0) ivw_row$b[1]   else NA
  se_ivw   <- if (nrow(ivw_row)>0) ivw_row$se[1]  else NA
  p_ivw    <- if (nrow(ivw_row)>0) ivw_row$pval[1] else NA

  # Step-3 (metabolite->cancer): grab from mediation_mr_results
  beta_met  <- if (nrow(wm_existing)>0) wm_existing$p3_b[1]  else NA

  # Indirect effect: product of coefficients (alpha * beta_step2)
  # Using WM for step-2 (primary) and IVW for step-2 (sensitivity)
  # Step-1 is always from protein->cancer IVW

  indirect_wm  <- if (!is.na(alpha) && !is.na(beta_wm))  alpha * beta_wm  else NA
  se_ind_wm    <- if (!is.na(se_alpha) && !is.na(se_wm))
    sqrt((alpha^2 * se_wm^2) + (beta_wm^2 * se_alpha^2)) else NA
  p_ind_wm     <- if (!is.na(indirect_wm) && !is.na(se_ind_wm))
    2*pnorm(-abs(indirect_wm/se_ind_wm)) else NA

  indirect_ivw <- if (!is.na(alpha) && !is.na(beta_ivw)) alpha * beta_ivw else NA
  se_ind_ivw   <- if (!is.na(se_alpha) && !is.na(se_ivw))
    sqrt((alpha^2 * se_ivw^2) + (beta_ivw^2 * se_alpha^2)) else NA
  p_ind_ivw    <- if (!is.na(indirect_ivw) && !is.na(se_ind_ivw))
    2*pnorm(-abs(indirect_ivw/se_ind_ivw)) else NA

  # Existing p_indirect from manuscript (from mediation_mr_results, WM primary)
  p_ind_existing <- if (nrow(wm_existing)>0) wm_existing$p_indirect[1] else NA

  cat(sprintf("  Step-2 WM:  beta=%7.4f SE=%6.4f p=%8.4e\n", beta_wm,  se_wm,  p_wm))
  cat(sprintf("  Step-2 IVW: beta=%7.4f SE=%6.4f p=%8.4e\n", beta_ivw, se_ivw, p_ivw))
  cat(sprintf("  Indirect WM:  %7.4f (p=%8.4e)\n", indirect_wm,  p_ind_wm))
  cat(sprintf("  Indirect IVW: %7.4f (p=%8.4e)\n", indirect_ivw, p_ind_ivw))
  cat(sprintf("  Existing p_indirect (manuscript): %8.4e\n", p_ind_existing))

  # Direction consistency
  dir_consistent <- if (!is.na(beta_wm) && !is.na(beta_ivw))
    ifelse(sign(beta_wm) == sign(beta_ivw), "YES", "NO") else NA

  rows[[length(rows)+1]] <- data.frame(
    protein               = pth$protein,
    metabolite            = pth$metabolite,
    cancer                = pth$cancer,
    n_snps_step2          = nrow(sub),

    # Step-1 (protein->cancer, IVW from STable2)
    step1_beta_IVW        = round(alpha, 4),
    step1_se_IVW          = round(se_alpha, 4),

    # Step-2 WM (primary)
    step2_beta_WM         = round(beta_wm, 4),
    step2_se_WM           = round(se_wm, 4),
    step2_p_WM            = signif(p_wm, 3),
    step2_method_WM       = if (nrow(sub)==1) "Wald ratio" else "Weighted median",

    # Step-2 IVW (sensitivity)
    step2_beta_IVW        = round(beta_ivw, 4),
    step2_se_IVW          = round(se_ivw, 4),
    step2_p_IVW           = signif(p_ivw, 3),
    step2_method_IVW      = if (nrow(sub)==1) "Wald ratio" else "Inverse variance weighted",

    # Direction concordance step-2
    step2_direction_WM_IVW = dir_consistent,

    # Indirect effect (product of coefficients)
    indirect_WM           = round(indirect_wm, 6),
    se_indirect_WM        = round(se_ind_wm, 6),
    p_indirect_WM         = signif(p_ind_wm, 3),
    OR_indirect_WM        = round(exp(indirect_wm), 6),

    indirect_IVW          = round(indirect_ivw, 6),
    se_indirect_IVW       = round(se_ind_ivw, 6),
    p_indirect_IVW        = signif(p_ind_ivw, 3),
    OR_indirect_IVW       = round(exp(indirect_ivw), 6),

    # Existing manuscript WM p_indirect (cross-check)
    p_indirect_manuscript = signif(p_ind_existing, 3),

    # Sensitivity flag
    IVW_confirms_WM       = ifelse(!is.na(p_ind_ivw) && p_ind_ivw < 0.05, "YES",
                              ifelse(!is.na(p_ind_ivw), "NO", NA)),
    stringsAsFactors = FALSE
  )
}

# ── Combine and save ─────────────────────────────────────────────────────────
results_df <- do.call(rbind, rows)

out_path <- file.path(PROJ, "results/tables/STable_Mediation_WM_vs_IVW_2026-05-25.csv")
write.csv(results_df, out_path, row.names=FALSE)
cat("\n\nFull results saved to:", out_path, "\n\n")

# ── Print summary comparison table ───────────────────────────────────────────
cat(sprintf("%-10s %-12s  %6s  %8s  %8s  %-5s  %8s  %8s  %-3s  %-3s\n",
    "Protein", "Metabolite", "nSNP",
    "p_WM", "p_IVW", "Dir?",
    "p_ind_WM", "p_ind_IVW", "WM_sig", "IVW_sig"))
cat(paste(rep("-", 100), collapse=""), "\n")

for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  cat(sprintf("%-10s %-12s  %6d  %8.2e  %8.2e  %-5s  %8.2e  %8.2e  %-3s  %-3s\n",
      r$protein, r$metabolite, r$n_snps_step2,
      r$step2_p_WM, r$step2_p_IVW, r$step2_direction_WM_IVW,
      r$p_indirect_WM, r$p_indirect_IVW,
      ifelse(r$p_indirect_WM < 0.05, "YES", "NO"),
      ifelse(r$IVW_confirms_WM == "YES", "YES", "NO")))
}
