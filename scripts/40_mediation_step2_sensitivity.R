#!/usr/bin/env Rscript
## Script 40: Formal two-step mediation MR — step-2 sensitivity table
## ─────────────────────────────────────────────────────────────────────
## Uses existing MR results (no new data fetch):
##   Step 1: protein → metabolite  (phase3_protein_metabolite/)
##   Step 2: metabolite → cancer   (phase4_metabolite_cancer/)
##   Total:  protein → cancer      (phase2_protein_cancer/)
##
## Indirect effect: b_indirect = b1 × b2
## Delta-method SE: sqrt(b1² × se2² + b2² × se1² + se1² × se2²)
## Proportion mediated: b_indirect / b_total × 100%
## Z-test for indirect: b_indirect / se_indirect
##
## Step-1 sensitivity: Wald ratio (1 SNP) or MR-RAPS (2 SNPs) — as available
## Step-2 sensitivity: IVW + Weighted median (both reported per path)
## Flag: single-SNP step-1 paths are Wald-only; caution warranted

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
})

proj    <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
out_dir <- file.path(proj, "results/mediation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ── Mediation triplets to test ────────────────────────────────────────────
## (protein, metabolite, cancer) — top 6 from mediation_validation_queue
triplets <- tribble(
  ~protein,   ~metabolite,  ~cancer,
  "IL34",     "Total_BCAA", "Breast_GCST90018757",
  "EFNA1",    "Total_BCAA", "Breast_GCST90018757",
  "TNFRSF6B", "Total_BCAA", "Breast_GCST90018757",
  "ATRAID",   "TG_by_PG",   "Breast_GCST90018757",
  "APOE",     "Gly",        "Breast_GCST90018757",
  "ITIH3",    "Gly",        "Breast_GCST90018757"
)

## ── Load MR result tables ──────────────────────────────────────────────────
step1_all <- fread(file.path(proj, "results/phase3_protein_metabolite/protein_metabolite_mr_results_full.csv"))
step2_all <- fread(file.path(proj, "results/phase4_metabolite_cancer/metabolite_cancer_mr_results_full.csv"))
total_all <- fread(file.path(proj, "results/phase2_protein_cancer/protein_cancer_mr_results_full.csv"))

## Step-2 methods to report
step2_methods <- c("Inverse variance weighted", "Weighted median")

## ── Delta-method indirect effect ──────────────────────────────────────────
indirect <- function(b1, se1, b2, se2) {
  b_ind  <- b1 * b2
  se_ind <- sqrt(b1^2 * se2^2 + b2^2 * se1^2 + se1^2 * se2^2)
  z      <- b_ind / se_ind
  pval   <- 2 * pnorm(-abs(z))
  list(b = b_ind, se = se_ind, z = z, pval = pval,
       ci_lo = b_ind - 1.96 * se_ind,
       ci_hi = b_ind + 1.96 * se_ind)
}

## ── Build results ─────────────────────────────────────────────────────────
rows <- list()

for (i in seq_len(nrow(triplets))) {
  prot <- triplets$protein[i]
  met  <- triplets$metabolite[i]
  can  <- triplets$cancer[i]

  ## Step 1: best available method for this protein-metabolite pair
  s1 <- step1_all[exposure == prot & outcome == met & !is.na(b)] |>
    arrange(desc(nsnp)) |>  # prefer more SNPs
    slice(1)

  if (nrow(s1) == 0) {
    message(sprintf("⚠ No step-1 result: %s → %s", prot, met)); next
  }

  ## Total effect
  tot <- total_all[exposure == prot & grepl(sub("_GCST.*","",can), outcome) &
                   !is.na(b)] |>
    arrange(desc(nsnp)) |>
    slice(1)
  b_total <- if (nrow(tot) > 0) tot$b[1] else NA_real_

  ## Step-2: test each method
  for (meth2 in step2_methods) {
    s2 <- step2_all[exposure == met & outcome == can &
                    method == meth2 & !is.na(b)]
    if (nrow(s2) == 0) next

    ind <- indirect(s1$b[1], s1$se[1], s2$b[1], s2$se[1])
    prop_med <- if (!is.na(b_total) && b_total != 0) {
      round(100 * ind$b / b_total, 1)
    } else NA_real_

    rows[[length(rows)+1]] <- tibble(
      protein          = prot,
      metabolite       = met,
      cancer           = "Breast",
      ## Step 1
      step1_method     = s1$method[1],
      step1_nsnp       = s1$nsnp[1],
      step1_b          = round(s1$b[1],   5),
      step1_se         = round(s1$se[1],  5),
      step1_p          = signif(s1$pval[1], 3),
      ## Step 2
      step2_method     = meth2,
      step2_nsnp       = s2$nsnp[1],
      step2_b          = round(s2$b[1],   5),
      step2_se         = round(s2$se[1],  5),
      step2_p          = signif(s2$pval[1], 3),
      ## Indirect
      b_indirect       = round(ind$b,    6),
      se_indirect      = round(ind$se,   6),
      ci95_lo          = round(ind$ci_lo,6),
      ci95_hi          = round(ind$ci_hi,6),
      p_indirect       = signif(ind$pval, 3),
      ## Total + proportion mediated
      b_total          = round(b_total,  5),
      prop_med_pct     = prop_med,
      ## Reliability flags
      step1_single_snp = (s1$nsnp[1] <= 1),
      step2_consistent = (sign(s2$b[1]) == sign(ind$b)),
      note             = case_when(
        s1$nsnp[1] <= 1 & meth2 == "Weighted median"
          ~ "Wald ratio step-1; WM step-2 is more robust for step-2 pleiotropy",
        s1$nsnp[1] <= 1
          ~ "Wald ratio step-1 — single-SNP; interpret with caution",
        TRUE ~ "Multi-SNP step-1 + multi-SNP step-2"
      )
    )
  }
}

results <- bind_rows(rows)

## ── Significance assessment ───────────────────────────────────────────────
results <- results |>
  mutate(
    sig_indirect = p_indirect < 0.05,
    both_steps_sig = step1_p < 0.05 & step2_p < 0.05,
    coloc_supported = protein %in% c("EFNA1","TNFRSF6B","ATRAID"),  # from script 39
    evidence_grade = case_when(
      coloc_supported & both_steps_sig & sig_indirect & !step1_single_snp
        ~ "Strong (coloc + multi-SNP MR + indirect p<0.05)",
      coloc_supported & sig_indirect
        ~ "Moderate (coloc + indirect p<0.05; Wald step-1)",
      sig_indirect & both_steps_sig
        ~ "Suggestive (indirect p<0.05; no coloc support)",
      sig_indirect
        ~ "Weak (indirect p<0.05; one step non-significant)",
      TRUE ~ "Not supported"
    )
  )

## ── Save ─────────────────────────────────────────────────────────────────
out_full    <- file.path(out_dir, "mediation_step2_sensitivity.csv")
write_csv(results, out_full)
message(sprintf("✓ Saved → %s", out_full))

## ── Print summary ─────────────────────────────────────────────────────────
cat("\n╔══════════════════════════════════════════════════════╗\n")
cat("  TWO-STEP MEDIATION MR — SENSITIVITY TABLE\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")

## Main display: IVW step-2 rows
main <- results |>
  filter(step2_method == "Inverse variance weighted") |>
  select(protein, metabolite, step1_method, step1_nsnp,
         step2_nsnp, b_indirect, p_indirect, prop_med_pct,
         b_total, evidence_grade)

print(main, n = Inf)

cat("\n── Step-2 sensitivity (IVW vs Weighted Median) ───────\n")
sens <- results |>
  select(protein, metabolite, step2_method, step2_nsnp,
         step2_b, step2_p, b_indirect, p_indirect, prop_med_pct) |>
  arrange(protein, metabolite, step2_method)
print(sens, n = Inf)

## ── Compact manuscript-ready table ───────────────────────────────────────
cat("\n── Manuscript summary (IVW step-2, main results) ─────\n")
ms <- results |>
  filter(step2_method == "Inverse variance weighted") |>
  mutate(
    indirect_fmt = sprintf("%.4f (95%%CI %.4f, %.4f; p=%s)",
                           b_indirect, ci95_lo, ci95_hi,
                           formatC(p_indirect, format="e", digits=2)),
    prop_fmt     = ifelse(is.na(prop_med_pct), "—",
                          sprintf("%.1f%%", prop_med_pct))
  ) |>
  select(protein, metabolite, step1_method, step1_nsnp,
         indirect_fmt, prop_fmt, evidence_grade)
print(ms, n = Inf)

## Also save WM-step2 comparison CSV
out_wm <- file.path(out_dir, "mediation_step2_wm_comparison.csv")
write_csv(
  results |> select(protein, metabolite, step2_method, step2_b, step2_p,
                    b_indirect, se_indirect, p_indirect, prop_med_pct,
                    evidence_grade),
  out_wm
)
message(sprintf("✓ WM comparison → %s", out_wm))
