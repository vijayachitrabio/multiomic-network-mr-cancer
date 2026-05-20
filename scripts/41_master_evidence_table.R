#!/usr/bin/env Rscript
## Script 41: Master evidence table — all 17 Phase-2 proteins
## ─────────────────────────────────────────────────────────────
## Consolidates: MR (phase2) + Coloc (STable8) + MAGMA (script 29-31)
##   + Mediation (script 40) + ER subtype (script 18) + Steiger (script 20)
##   + Druggability (script 22) → STable_master + tier assignments

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
})

proj    <- "."
out_dir <- file.path(proj, "results/tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ── 1. Phase-2 MR (protein → cancer; breast + endometrial) ───────────────
mr_raw <- fread(file.path(proj, "results/phase2_protein_cancer/protein_cancer_mr_results_significant.csv"))

# Keep one row per protein: prefer breast, fall back to endometrial (ABO)
mr <- mr_raw |>
  mutate(cancer_short = ifelse(grepl("Breast", outcome), "Breast", "Endometrial")) |>
  group_by(exposure) |>
  arrange(ifelse(cancer_short == "Breast", 0L, 1L)) |>   # breast first
  slice(1) |>
  ungroup() |>
  select(protein = exposure, cancer_mr = cancer_short,
         mr_method = method, mr_nsnp = nsnp,
         mr_b = b, mr_se = se, mr_pval = pval,
         mr_or = or, mr_or_lo = or_lci95, mr_or_hi = or_uci95, mr_fdr = fdr)

## ── 2. Coloc — unified table from script 42 ───────────────────────────────
coloc_all <- fread(file.path(proj, "results/tables/STable8_protein_coloc.csv")) |>
  select(protein,
         coloc_PPH3_abf = PPH3_abf,
         coloc_PPH4_abf = PPH4_abf,
         coloc_interpretation_abf = coloc_interp_abf,
         coloc_PPH4_susie = PPH4_susie,
         coloc_n_pairs = n_coloc_pairs,
         coloc_PPH4_best = PPH4_best,
         coloc_method_best = best_method,
         coloc_verdict = final_verdict) |>
  mutate(
    coloc_verdict = case_when(
      grepl("STRONG", coloc_verdict) ~ "STRONG",
      grepl("MODERATE", coloc_verdict) ~ "MODERATE",
      grepl("DISTINCT", coloc_verdict) ~ "DISTINCT",
      TRUE ~ "INSUFFICIENT"
    )
  )

## ── 3. MAGMA ─────────────────────────────────────────────────────────────
magma <- fread(file.path(proj, "results/tables/STable6_magma_mr_hit_lookup.csv")) |>
  select(protein, magma_breast_p = breast_p,
         magma_breast_bonf = breast_bonferroni_sig,
         magma_breast_rank = breast_rank)

## ── 4. Mediation (WM step-2 as primary) ──────────────────────────────────
med_raw <- fread(file.path(proj, "results/mediation/mediation_step2_sensitivity.csv"))
med <- med_raw |>
  filter(step2_method == "Weighted median") |>
  mutate(mediation_supported = p_indirect < 0.05 & step1_p < 0.05) |>
  group_by(protein) |>
  # keep best-supported path per protein
  arrange(p_indirect) |>
  slice(1) |>
  ungroup() |>
  select(protein, med_metabolite = metabolite,
         med_p_indirect_wm = p_indirect,
         med_prop_pct = prop_med_pct,
         med_supported = mediation_supported)

## ── 5. ER subtype ─────────────────────────────────────────────────────────
er <- fread(file.path(proj, "results/er_subtype/er_subtype_comparison.csv")) |>
  select(protein = exposure, er_pattern,
         or_ERpos, pval_ERpos, or_ERneg, pval_ERneg)

## ── 6. Steiger ────────────────────────────────────────────────────────────
steiger <- fread(file.path(proj, "results/sensitivity/steiger_comparison_table.csv")) |>
  filter(grepl("Breast", cancer)) |>
  select(protein, steiger_ok = direction_consistent,
         steiger_n_reversed = n_snp_rev)

## ── 7. Druggability ───────────────────────────────────────────────────────
drug <- fread(file.path(proj, "results/pathway/opentargets_druggability.csv")) |>
  select(protein, n_drugs = n_known_drugs,
         tractability_SM, tractability_AB,
         top_drug, top_drug_phase)

## ── 8. Assemble ───────────────────────────────────────────────────────────
master <- mr |>
  left_join(coloc_all,  by = "protein") |>
  left_join(magma,      by = "protein") |>
  left_join(med,        by = "protein") |>
  left_join(er,         by = "protein") |>
  left_join(steiger,    by = "protein") |>
  left_join(drug,       by = "protein")

## ── 9. Tier assignment ────────────────────────────────────────────────────
master <- master |>
  mutate(
    tier = case_when(
      # Tier 1: MR + strong coloc (SuSiE PPH4 ≥ 0.8)
      !is.na(coloc_PPH4_best) & coloc_PPH4_best >= 0.8
        ~ "Tier 1 — MR + strong coloc",
      # Tier 2a: MR + MAGMA Bonferroni (SNX15, PM20D1)
      !is.na(magma_breast_bonf) & magma_breast_bonf &
        (is.na(coloc_PPH4_best) | coloc_PPH4_best < 0.5)
        ~ "Tier 2a — MR + MAGMA Bonferroni",
      # Tier 2b: MR + moderate coloc
      !is.na(coloc_PPH4_best) & coloc_PPH4_best >= 0.5
        ~ "Tier 2b — MR + moderate coloc",
      # Tier 2c: MR + coloc tested but not supported / not tested
      !is.na(coloc_verdict)
        ~ "Tier 2c — MR only (coloc not supportive)",
      # Tier 2d: not yet tested for coloc
      TRUE ~ "Tier 2d — MR only (coloc pending)"
    ),
    tier_short = case_when(
      grepl("Tier 1", tier)  ~ "T1",
      grepl("Tier 2a", tier) ~ "T2a",
      grepl("Tier 2b", tier) ~ "T2b",
      grepl("Tier 2c", tier) ~ "T2c",
      TRUE                   ~ "T2d"
    )
  )

## ── 10. Save full master table ────────────────────────────────────────────
out_master <- file.path(out_dir, "STable_master_evidence.csv")
write_csv(master, out_master)
message(sprintf("✓ Master table → %s  (%d proteins, %d columns)",
                out_master, nrow(master), ncol(master)))

## ── 11. Print summary ─────────────────────────────────────────────────────
cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("  MASTER EVIDENCE TABLE — ALL 17 PROTEINS\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

display <- master |>
  arrange(tier_short, mr_pval) |>
  mutate(
    OR_fmt      = sprintf("%.3f (%.3f–%.3f)", mr_or, mr_or_lo, mr_or_hi),
    coloc_fmt   = ifelse(is.na(coloc_PPH4_best), "—",
                         sprintf("%.3f (%s)", coloc_PPH4_best, coloc_method_best)),
    magma_fmt   = ifelse(is.na(magma_breast_p), "—",
                         ifelse(magma_breast_bonf,
                                sprintf("%.2e ✓Bonf", magma_breast_p),
                                sprintf("%.2e", magma_breast_p))),
    med_fmt     = ifelse(is.na(med_p_indirect_wm), "—",
                         ifelse(med_supported,
                                sprintf("%s p=%.3f (%.1f%%)",
                                        med_metabolite, med_p_indirect_wm, med_prop_pct),
                                sprintf("%s NS", med_metabolite))),
    drug_fmt    = ifelse(n_drugs == 0, "novel", sprintf("%d drugs", n_drugs))
  ) |>
  select(tier_short, protein, mr_nsnp, OR_fmt, coloc_fmt,
         magma_fmt, med_fmt, er_pattern, drug_fmt)

print(as_tibble(display), n = Inf)

## ── 12. Tier summary ─────────────────────────────────────────────────────
cat("\n── Tier breakdown ────────────────────────────────────────\n")
tier_sum <- master |>
  count(tier) |>
  arrange(tier)
print(tier_sum)

## ── 13. Manuscript-ready compact version (key columns only) ───────────────
compact <- master |>
  arrange(tier_short, mr_pval) |>
  transmute(
    Tier          = tier_short,
    Protein       = protein,
    `MR OR (95% CI)` = sprintf("%.3f (%.3f–%.3f)", mr_or, mr_or_lo, mr_or_hi),
    `MR p`        = formatC(mr_pval, format="e", digits=2),
    `MR FDR`      = formatC(mr_fdr,  format="e", digits=2),
    `MR SNPs`     = mr_nsnp,
    `Coloc PPH4`  = ifelse(is.na(coloc_PPH4_best), "—",
                            sprintf("%.3f (%s)", coloc_PPH4_best, coloc_method_best)),
    `Coloc verdict` = ifelse(is.na(coloc_verdict), "not tested", coloc_verdict),
    `MAGMA p`     = ifelse(is.na(magma_breast_p), "—",
                            sprintf("%.2e%s", magma_breast_p,
                                    ifelse(!is.na(magma_breast_bonf) & magma_breast_bonf," *",""))),
    `Mediation`   = ifelse(is.na(med_metabolite), "—",
                            ifelse(med_supported,
                                   sprintf("→%s ✓", med_metabolite),
                                   sprintf("→%s NS", med_metabolite))),
    `Prop.med%`   = ifelse(is.na(med_prop_pct), "—",
                            sprintf("%.1f%%", med_prop_pct)),
    `ER pattern`  = ifelse(is.na(er_pattern), "—", er_pattern),
    `Steiger OK`  = ifelse(is.na(steiger_ok), "—",
                            ifelse(steiger_ok, "✓", "✗")),
    `Drugs`       = ifelse(n_drugs == 0, "novel", as.character(n_drugs))
  )

out_compact <- file.path(out_dir, "STable7_master_evidence_compact.csv")
write_csv(compact, out_compact)
message(sprintf("✓ Compact table → %s", out_compact))

cat("\n── Compact manuscript table ───────────────────────────────\n")
print(as_tibble(compact), n = Inf)
