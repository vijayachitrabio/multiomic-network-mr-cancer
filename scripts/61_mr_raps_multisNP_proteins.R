#!/usr/bin/env Rscript
# Script 61: MR-RAPS sensitivity analysis for multi-SNP proteins
# Proteins: ABO (endometrial), KLB, PM20D1, IL34 (breast)
# Date: 2026-05-25

suppressPackageStartupMessages(library(mr.raps))

PROJ <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"

dat_breast <- readRDS(file.path(PROJ, "data/harmonised/harmonised_protein_Breast_GCST90018757.rds"))
dat_endo   <- readRDS(file.path(PROJ, "data/harmonised/harmonised_protein_Endometrial_GCST006464.rds"))

targets <- list(
  list(protein = "KLB",    cancer = "Breast",      dat = dat_breast),
  list(protein = "PM20D1", cancer = "Breast",      dat = dat_breast),
  list(protein = "IL34",   cancer = "Breast",      dat = dat_breast),
  list(protein = "ABO",    cancer = "Endometrial", dat = dat_endo)
)

# Load IVW reference
stab <- read.csv(file.path(PROJ, "results/tables/STable2_17_FDR_hits_complete.csv"),
                 stringsAsFactors = FALSE)

rows <- list()

for (tgt in targets) {
  sub <- tgt$dat[grepl(tgt$protein, tgt$dat$exposure, ignore.case = TRUE) &
                   tgt$dat$mr_keep == TRUE, ]

  cat("\n===", tgt$protein, "->", tgt$cancer, "===\n")
  cat("N SNPs:", nrow(sub), "\n")

  if (nrow(sub) < 2) { cat("SKIP: <2 SNPs\n"); next }

  # Build data frame in format mr.raps expects
  input_df <- data.frame(
    beta.exposure = sub$beta.exposure,
    beta.outcome  = sub$beta.outcome,
    se.exposure   = sub$se.exposure,
    se.outcome    = sub$se.outcome,
    stringsAsFactors = FALSE
  )

  # Run three variants (suppress diagnostic plots)
  pdf(NULL)  # suppress plots

  raps_l2 <- tryCatch(
    mr.raps(input_df, diagnostics = FALSE, over.dispersion = FALSE, loss.function = "l2"),
    error = function(e) { cat("  l2 failed:", conditionMessage(e), "\n"); NULL }
  )

  raps_over <- tryCatch(
    mr.raps(input_df, diagnostics = FALSE, over.dispersion = TRUE, loss.function = "l2"),
    error = function(e) { cat("  overdispersed failed:", conditionMessage(e), "\n"); NULL }
  )

  raps_huber <- tryCatch(
    mr.raps(input_df, diagnostics = FALSE, over.dispersion = TRUE, loss.function = "huber"),
    error = function(e) { cat("  huber failed:", conditionMessage(e), "\n"); NULL }
  )

  dev.off()

  # IVW reference from STable2
  ivw_row <- stab[grepl(tgt$protein, stab$protein, ignore.case = TRUE) &
                    grepl(sub("_GCST.*", "", tgt$cancer), stab$cancer, ignore.case = TRUE), ]
  ivw_OR <- if (nrow(ivw_row) > 0) ivw_row$OR[1] else NA
  ivw_p  <- if (nrow(ivw_row) > 0) ivw_row$pvalue[1] else NA

  # Print results
  for (nm in c("l2 (no overdispersion)", "l2 (overdispersed)", "huber (robust)")) {
    obj <- switch(gsub(" .*", "", nm),
                  "l2" = if (grepl("overdispersed", nm)) raps_over else raps_l2,
                  "huber" = raps_huber)
    if (!is.null(obj)) {
      b <- obj$beta.hat; se <- obj$beta.se
      p <- 2 * pnorm(-abs(b / se))
      tau2 <- if (!is.null(obj$tau2)) obj$tau2 else NA
      cat(sprintf("  %-28s OR=%6.4f [%6.4f-%6.4f]  p=%8.4e  tau2=%s\n",
                  nm, exp(b), exp(b-1.96*se), exp(b+1.96*se), p,
                  if (is.na(tau2)) "N/A" else sprintf("%.6f", tau2)))
    }
  }

  # Save primary (robust huber) and simple estimate
  primary <- if (!is.null(raps_huber)) raps_huber else raps_l2
  if (!is.null(primary)) {
    b <- primary$beta.hat; se <- primary$beta.se
    p <- 2 * pnorm(-abs(b / se))
    tau2 <- if (!is.null(primary$tau2)) primary$tau2 else NA

    rows[[length(rows) + 1]] <- data.frame(
      protein        = tgt$protein,
      cancer         = tgt$cancer,
      n_snps         = nrow(sub),
      snps           = paste(sub$SNP, collapse = "; "),
      IVW_OR         = round(ivw_OR, 4),
      IVW_pvalue     = ivw_p,
      RAPS_method    = if (!is.null(raps_huber)) "huber_overdispersed" else "l2_simple",
      RAPS_beta      = round(b, 4),
      RAPS_SE        = round(se, 4),
      RAPS_OR        = round(exp(b), 4),
      RAPS_OR_95lo   = round(exp(b - 1.96*se), 4),
      RAPS_OR_95hi   = round(exp(b + 1.96*se), 4),
      RAPS_pvalue    = round(p, 6),
      tau2           = if (!is.na(tau2)) round(tau2, 6) else NA,
      pleiotropy_flag = if (!is.na(tau2) && tau2 > 0.01) "possible" else "none",
      direction_consistent_with_IVW = if (!is.na(ivw_OR)) ifelse(sign(b) == sign(log(ivw_OR)), "YES", "NO") else NA,
      stringsAsFactors = FALSE
    )
  }
}

# Save results
results_df <- do.call(rbind, rows)
out_path <- file.path(PROJ, "results/tables/STable_MR_RAPS_sensitivity_2026-05-25.csv")
write.csv(results_df, out_path, row.names = FALSE)
cat("\n\nResults saved to:", out_path, "\n\n")

# Summary table
cat(sprintf("%-10s %-12s %5s  %-22s  %8s  %8s  %-8s  %s\n",
            "Protein", "Cancer", "nSNP", "MR-RAPS OR [95%CI]", "p-RAPS", "p-IVW", "tau2", "Consistent"))
cat(paste(rep("-", 100), collapse=""), "\n")
for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  cat(sprintf("%-10s %-12s %5d  %5.4f [%5.4f-%5.4f]  %8.4e  %8.4e  %-8s  %s\n",
              r$protein, r$cancer, r$n_snps,
              r$RAPS_OR, r$RAPS_OR_95lo, r$RAPS_OR_95hi,
              r$RAPS_pvalue, r$IVW_pvalue,
              if (is.na(r$tau2)) "N/A" else as.character(r$tau2),
              r$direction_consistent_with_IVW))
}
