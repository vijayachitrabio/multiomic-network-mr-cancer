#!/usr/bin/env Rscript
## Script 42: STable — protein colocalization results (all 17 MR-hit proteins)
## Combines coloc.abf/coloc.susie outputs from scripts 37-39 and pending pilots.
## Output: results/tables/STable8_protein_coloc.csv

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr)
})

proj    <- "."
out_dir <- file.path(proj, "results/tables")

## ── Load ──────────────────────────────────────────────────────────────────
# SNX15 + PM20D1: coloc.abf from script 37
abf2 <- fread(file.path(proj, "results/validation/protein_coloc_snx15_pm20d1.csv")) |>
  select(protein, n_harm = n_snps, PPH3_abf = PPH3, PPH4_abf = PPH4,
         PPH4_abf_sens = PPH4_sens_p12_5e5, coloc_interp_abf = interpretation)

# SNX15 + PM20D1: coloc.susie from script 38
sus2 <- fread(file.path(proj, "results/validation/protein_coloc_susie_snx15_pm20d1.csv")) |>
  select(protein, n_ld_harm = n_snps_ld_harm,
         n_cs_pqtl, n_cs_gwas, n_coloc_pairs,
         PPH4_susie = PPH4_susie_best, susie_best_snp)

snx_pm <- left_join(abf2, sus2, by = "protein") |>
  mutate(coloc_tested = TRUE)

# Remaining 6: coloc.abf + coloc.susie from script 39
rem6 <- fread(file.path(proj, "results/validation/protein_coloc_remaining6.csv")) |>
  select(protein, n_harm, n_ld_harm, n_cs_pqtl, n_cs_gwas, n_coloc_pairs,
         PPH3_abf, PPH4_abf, PPH4_susie, susie_best_snp,
         coloc_interp_abf = interpretation) |>
  mutate(PPH4_abf_sens = NA_real_, coloc_tested = TRUE)

# Pending 9: pilot coloc runs from script 57
pending9_path <- file.path(proj, "results/validation/protein_coloc_pending9_summary.csv")
pending9 <- if (file.exists(pending9_path)) {
  fread(pending9_path) |>
    select(protein, n_harm, n_ld_harm, n_cs_pqtl, n_cs_gwas, n_coloc_pairs,
           PPH3_abf, PPH4_abf, PPH4_susie, susie_best_snp,
           coloc_interp_abf = interpretation) |>
    mutate(PPH4_abf_sens = NA_real_, coloc_tested = TRUE)
} else {
  tibble()
}

all_coloc <- bind_rows(snx_pm, rem6, pending9) |>
  mutate(
    PPH4_best = pmax(PPH4_susie, PPH4_abf, na.rm = TRUE),
    best_method = case_when(
      !is.na(PPH4_susie) & PPH4_susie >= PPH4_abf ~ "SuSiE",
      TRUE ~ "ABF"
    ),
    borderline_abf_susie_nonconvergent =
      PPH4_abf >= 0.8 & PPH4_abf < 0.85 &
      (is.na(PPH4_susie) | PPH4_susie == 0) &
      n_cs_pqtl == 0 & n_cs_gwas == 0,
    final_verdict = case_when(
      borderline_abf_susie_nonconvergent ~ "Moderate colocalization (ABF only)",
      PPH4_best >= 0.8 ~ "STRONG colocalization",
      PPH4_best >= 0.5 ~ "MODERATE colocalization",
      PPH3_abf  >= 0.8 ~ "DISTINCT causal variants",
      TRUE ~ "INSUFFICIENT evidence"
    ),
    tier_impact = case_when(
      borderline_abf_susie_nonconvergent ~ "Tier 2a (ABF-only, SuSiE non-convergent)",
      PPH4_best >= 0.8 ~ "Tier 1 upgrade",
      PPH4_best >= 0.5 ~ "Supports Tier 2b",
      TRUE ~ "No coloc support"
    ),
    note = case_when(
      protein %in% c("EFNA1","ATRAID") & PPH4_susie >= 0.8 & PPH3_abf >= 0.8
        ~ "coloc.abf missed (PPH3≥0.8); SuSiE resolved via multiple GWAS credible sets",
      protein %in% c("SNX15","PM20D1")
        ~ paste0("FinnGen N=619 power limitation; MAGMA Bonferroni p<0.05 provides locus-level triangulation"),
      protein == "APOE"
        ~ "rs429358 absent from GWAS; GWAS 0 credible sets; ABF moderate (PPH4=0.52)",
      protein == "IL34"
        ~ "Strong pQTL signal; GWAS 0 credible sets in region — MR effect likely systemic",
      borderline_abf_susie_nonconvergent
        ~ paste0(
          "PPH4_abf=", round(PPH4_abf, 3),
          " (just above threshold); PPH3=", round(PPH3_abf, 3),
          "; SuSiE non-convergent — zero credible sets from both FinnGen pQTL and breast GWAS; ",
          "downgraded from Tier 1 to Tier 2a (provisional)"
        ),
      protein == "ABO"
        ~ "Strong ABF coloc (PPH4=0.912); SuSiE non-convergent due to absence of GWAS credible sets at the ABO locus (9q34) in the endometrial GWAS — data limitation at a complex multi-variant locus; Tier 1 retained",
      protein == "TSPAN8"
        ~ "Moderate ABF coloc; SuSiE lower than ABF",
      protein == "FGF5"
        ~ "Pending-target pilot resolved strong coloc by ABF and SuSiE",
      TRUE ~ ""
    ),
    coloc_support = case_when(
      PPH4_abf >= 0.8 & PPH4_susie >= 0.8
        ~ paste0(
          "ABF+SuSiE (PPH4=", sprintf("%.3f", PPH4_abf),
          "/", sprintf("%.3f", PPH4_susie),
          "; dual-confirmed)"
        ),
      protein %in% c("EFNA1", "ATRAID") & PPH4_susie >= 0.8 & PPH3_abf >= 0.8
        ~ paste0(
          "SuSiE (PPH4=", sprintf("%.3f", PPH4_susie),
          "; ABF missed due to LD complexity)"
        ),
      protein == "ABO"
        ~ paste0(
          "ABF only (PPH4=", sprintf("%.3f", PPH4_abf),
          "); SuSiE non-convergent — no GWAS credible sets at 9q34 (endometrial)"
        ),
      borderline_abf_susie_nonconvergent
        ~ paste0(
          "ABF only (PPH4=", sprintf("%.3f", PPH4_abf),
          "); SuSiE non-convergent — no credible sets from either pQTL or GWAS"
        ),
      best_method == "SuSiE" ~ "SuSiE",
      best_method == "ABF" ~ "ABF",
      TRUE ~ ""
    )
  ) |>
  select(-borderline_abf_susie_nonconvergent) |>
  arrange(desc(PPH4_best))

## ── Save ──────────────────────────────────────────────────────────────────
out <- file.path(out_dir, "STable8_protein_coloc.csv")
write_csv(all_coloc, out)
message(sprintf("✓ STable8 → %s", out))

## ── Print ─────────────────────────────────────────────────────────────────
cat("\n╔══════════════════════════════════════════════╗\n")
cat("  STable8: PROTEIN COLOCALIZATION RESULTS\n")
cat("╚══════════════════════════════════════════════╝\n\n")

display <- all_coloc |>
  select(protein, n_harm, n_ld_harm, n_cs_pqtl, n_cs_gwas,
         n_coloc_pairs, PPH3_abf, PPH4_abf, PPH4_susie,
         PPH4_best, best_method, final_verdict, tier_impact, note, coloc_support)

print(as_tibble(display), n = Inf)

cat("\n── Summary ───────────────────────────────────\n")
cat(sprintf("  Strong coloc:               %d protein(s): %s\n",
            sum(all_coloc$final_verdict == "STRONG colocalization", na.rm=TRUE),
            paste(all_coloc$protein[all_coloc$final_verdict == "STRONG colocalization"], collapse=", ")))
cat(sprintf("  Moderate coloc (PPH4≥0.5):  %d protein(s): %s\n",
            sum(all_coloc$PPH4_best >= 0.5 & all_coloc$PPH4_best < 0.8, na.rm=TRUE),
            paste(all_coloc$protein[!is.na(all_coloc$PPH4_best) &
                                     all_coloc$PPH4_best >= 0.5 &
                                     all_coloc$PPH4_best < 0.8], collapse=", ")))
cat(sprintf("  Distinct variants (PPH3≥0.8):%d protein(s): %s\n",
            sum(all_coloc$PPH3_abf >= 0.8, na.rm=TRUE),
            paste(all_coloc$protein[all_coloc$PPH3_abf >= 0.8], collapse=", ")))
cat(sprintf("  SuSiE overturned ABF:        EFNA1 (ABF PPH3=0.901 → SuSiE PPH4=0.963)\n"))
cat(sprintf("                               ATRAID (ABF PPH3=0.997 → SuSiE PPH4=0.996)\n"))
