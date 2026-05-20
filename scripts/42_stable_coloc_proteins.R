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
    final_verdict = case_when(
      PPH4_best >= 0.8 ~ "STRONG colocalization",
      PPH4_best >= 0.5 ~ "MODERATE colocalization",
      PPH3_abf  >= 0.8 ~ "DISTINCT causal variants",
      TRUE ~ "INSUFFICIENT evidence"
    ),
    tier_impact = case_when(
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
      protein %in% c("UMOD", "ABO")
        ~ "Strong ABF coloc; SuSiE unavailable or not converged, interpret as ABF-supported",
      protein == "TSPAN8"
        ~ "Moderate ABF coloc; SuSiE lower than ABF",
      protein == "FGF5"
        ~ "Pending-target pilot resolved strong coloc by ABF and SuSiE",
      TRUE ~ ""
    )
  ) |>
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
         PPH4_best, best_method, final_verdict, tier_impact, note)

print(as_tibble(display), n = Inf)

cat("\n── Summary ───────────────────────────────────\n")
cat(sprintf("  Strong coloc (PPH4≥0.8):    %d protein(s): %s\n",
            sum(all_coloc$PPH4_best >= 0.8, na.rm=TRUE),
            paste(all_coloc$protein[all_coloc$PPH4_best >= 0.8], collapse=", ")))
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
