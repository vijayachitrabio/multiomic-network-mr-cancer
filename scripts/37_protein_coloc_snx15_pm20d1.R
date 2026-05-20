#!/usr/bin/env Rscript
## Script 37: Protein-side colocalization for SNX15 and PM20D1 (v2 — GRCh38 direct match)
## ─────────────────────────────────────────────────────────────────────────────────────
## FinnGen pQTL (GRCh38) matched directly to GWAS (GRCh38 harmonised .h file).
## No liftover required.
##
## FinnGen pQTL URL:
##   storage.googleapis.com/finngen-public-data-r10/omics/proteomics/release_2023_03_02/
##   data/Olink/pQTL/Olink_Batch1_{PROTEIN}.txt.gz
## File columns (no header): chr pos variant_id ref alt alt_freq beta se t_stat p log10_p [n]
## Coordinates: GRCh38
##
## GWAS: Breast_GCST90018757.h.tsv.gz — GWAS Catalog harmonised → GRCh38
## Confirmed: base_pair_location = GRCh38 (hm_coordinate_conversion present, hm_code 65=forward)

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(Rsamtools)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
})

set.seed(2026)
proj    <- "."
out_dir <- file.path(proj, "results/validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(proj, "data/pqtl/priority_regions"), showWarnings = FALSE)

## ── Configuration ────────────────────────────────────────────────────────────
WINDOW_BP <- 500000L          # ±500 kb cis window
MIN_SNPS  <- 50L
BASE_URL  <- "https://storage.googleapis.com/finngen-public-data-r10/omics/proteomics/release_2023_03_02/data/Olink/pQTL"

## GWAS (GRCh38 harmonised)
gwas_cfg <- list(
  cancer     = "Breast",
  file       = file.path(proj, "data/cancer_gwas/Breast_GCST90018757.h.tsv.gz"),
  n_cases    = 122977L,
  n_controls = 105974L
)
gwas_cfg$n_total <- gwas_cfg$n_cases + gwas_cfg$n_controls
gwas_cfg$s       <- gwas_cfg$n_cases / gwas_cfg$n_total

## Proteins (GRCh38 instrument positions from FinnGen)
targets <- list(
  SNX15  = list(chr = "11", pos_hg38 = 65025679L),
  PM20D1 = list(chr =  "1", pos_hg38 = 205847765L)
)

## ── Tabix extractor (GRCh38 region) ─────────────────────────────────────────
extract_pqtl_region <- function(protein, chr, pos_hg38, window = WINDOW_BP) {
  url    <- sprintf("%s/Olink_Batch1_%s.txt.gz", BASE_URL, protein)
  start  <- max(1L, pos_hg38 - window)
  end    <- pos_hg38 + window
  region <- GRanges(chr, IRanges(start, end))

  message(sprintf("  Tabix fetch: chr%s:%d-%d from %s", chr, start, end, protein))
  tbx   <- TabixFile(url)
  lines <- scanTabix(tbx, param = region)[[1]]
  message(sprintf("  %d lines returned", length(lines)))
  if (length(lines) == 0) return(NULL)

  x <- fread(text = paste(lines, collapse = "\n"), header = FALSE)
  col_names <- c("chr","pos","variant_id","ref","alt","alt_freq",
                 "beta","se","t_stat","p","log10_p","n")
  if (ncol(x) == 11) col_names <- col_names[-12]
  setnames(x, col_names[seq_len(ncol(x))])
  if (!"n" %in% names(x)) x[, n := 619L]
  x[, protein := protein]
  x[!is.na(beta) & !is.na(se) & se > 0 & !is.na(alt_freq) & alt_freq > 0 & alt_freq < 1]
}

## ── GWAS region loader (GRCh38) ───────────────────────────────────────────────
load_gwas_region <- function(chr, start_hg38, end_hg38) {
  message(sprintf("  Loading breast GWAS chr%s:%d-%d (GRCh38)...", chr, start_hg38, end_hg38))
  gwas <- read_tsv(gwas_cfg$file, show_col_types = FALSE,
                   col_types = cols_only(
                     chromosome              = col_integer(),
                     base_pair_location      = col_integer(),
                     effect_allele           = col_character(),
                     other_allele            = col_character(),
                     beta                    = col_double(),
                     standard_error          = col_double(),
                     effect_allele_frequency = col_double(),
                     p_value                 = col_double()
                   )) |>
    filter(chromosome == as.integer(chr),
           base_pair_location >= start_hg38,
           base_pair_location <= end_hg38) |>
    rename(pos = base_pair_location,
           ea = effect_allele, oa = other_allele,
           beta_gwas = beta, se_gwas = standard_error,
           eaf_gwas = effect_allele_frequency, p_gwas = p_value) |>
    filter(!is.na(beta_gwas), !is.na(se_gwas), se_gwas > 0,
           !is.na(eaf_gwas), eaf_gwas > 0, eaf_gwas < 1)
  message(sprintf("  %d GWAS variants in region", nrow(gwas)))
  gwas
}

## ── Harmonise pQTL ↔ GWAS (both in GRCh38) ──────────────────────────────────
harmonise_datasets <- function(pqtl, gwas) {
  flip <- c(A="T", T="A", C="G", G="C")
  merged <- inner_join(
    pqtl |> mutate(pos = as.integer(pos)),
    gwas |> mutate(pos = as.integer(pos)),
    by = "pos"
  ) |>
    mutate(
      ea_p  = toupper(alt), oa_p  = toupper(ref),
      ea_g  = toupper(ea),  oa_g  = toupper(oa),
      ea_pf = flip[ea_p],   oa_pf = flip[oa_p],
      match_d  = (ea_p == ea_g  & oa_p == oa_g),
      match_s  = (ea_p == oa_g  & oa_p == ea_g),
      match_fl = (!is.na(ea_pf) & ea_pf == ea_g & oa_pf == oa_g),
      match_fs = (!is.na(ea_pf) & ea_pf == oa_g & oa_pf == ea_g),
      palin    = (ea_p == flip[oa_p])
    ) |>
    filter(!palin | (alt_freq > 0.1 & alt_freq < 0.9)) |>
    filter(match_d | match_s | match_fl | match_fs) |>
    mutate(
      beta_gwas_h = if_else(match_d | match_fl, beta_gwas, -beta_gwas),
      eaf_gwas_h  = if_else(match_d | match_fl, eaf_gwas, 1 - eaf_gwas)
    )
  message(sprintf("  %d harmonised SNPs", nrow(merged)))
  merged
}

## ── COLOC ─────────────────────────────────────────────────────────────────────
run_coloc <- function(harm, protein, n_pqtl = 619) {
  if (nrow(harm) < MIN_SNPS) {
    message(sprintf("  ⚠ Only %d SNPs — below MIN_SNPS=%d, skipping", nrow(harm), MIN_SNPS))
    return(NULL)
  }
  snp_ids <- paste0("chr", harm$chr, ":", harm$pos, ":", harm$ea_p, ":", harm$oa_p)

  D1 <- list(   # pQTL (quantitative)
    beta    = harm$beta,
    varbeta = harm$se^2,
    snp     = snp_ids,
    type    = "quant",
    N       = n_pqtl,
    MAF     = pmin(harm$alt_freq, 1 - harm$alt_freq)
  )
  D2 <- list(   # GWAS (case-control)
    beta    = harm$beta_gwas_h,
    varbeta = harm$se_gwas^2,
    snp     = snp_ids,
    type    = "cc",
    N       = gwas_cfg$n_total,
    s       = gwas_cfg$s,
    MAF     = pmin(harm$eaf_gwas_h, 1 - harm$eaf_gwas_h)
  )

  message(sprintf("  Running coloc.abf with %d SNPs ...", nrow(harm)))

  ## Default priors (p12 = 1e-5)
  res <- coloc.abf(D1, D2)
  pp  <- res$summary
  message(sprintf("  [p12=1e-5]  PPH0=%.3f  PPH1=%.3f  PPH2=%.3f  PPH3=%.3f  PPH4=%.3f",
                  pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"],
                  pp["PP.H3.abf"], pp["PP.H4.abf"]))

  ## Sensitive priors (p12 = 5e-5) — for reporting
  res_s <- coloc.abf(D1, D2, p12 = 5e-5)
  pp_s  <- res_s$summary
  message(sprintf("  [p12=5e-5]  PPH0=%.3f  PPH1=%.3f  PPH2=%.3f  PPH3=%.3f  PPH4=%.3f",
                  pp_s["PP.H0.abf"], pp_s["PP.H1.abf"], pp_s["PP.H2.abf"],
                  pp_s["PP.H3.abf"], pp_s["PP.H4.abf"]))

  top_snp_df <- tryCatch({
    df <- as.data.frame(res$results)
    df[order(-df$SNP.PP.H4.abf), ][1, , drop = FALSE]
  }, error = function(e) NULL)

  coloc_lead_snp   <- if (!is.null(top_snp_df)) top_snp_df$snp[1]         else NA_character_
  coloc_lead_PP_H4 <- if (!is.null(top_snp_df)) round(top_snp_df$SNP.PP.H4.abf[1], 4) else NA_real_

  tibble(
    protein              = protein,
    cancer               = "Breast",
    n_snps               = nrow(harm),
    PPH0                 = round(pp["PP.H0.abf"],   4),
    PPH1                 = round(pp["PP.H1.abf"],   4),
    PPH2                 = round(pp["PP.H2.abf"],   4),
    PPH3                 = round(pp["PP.H3.abf"],   4),
    PPH4                 = round(pp["PP.H4.abf"],   4),
    PPH4_sens_p12_5e5    = round(pp_s["PP.H4.abf"], 4),
    PPH3_sens_p12_5e5    = round(pp_s["PP.H3.abf"], 4),
    coloc_lead_snp       = coloc_lead_snp,
    coloc_lead_PP_H4     = coloc_lead_PP_H4,
    interpretation       = case_when(
      pp["PP.H4.abf"] >= 0.8  ~ "STRONG colocalization",
      pp["PP.H4.abf"] >= 0.5  ~ "MODERATE colocalization",
      pp_s["PP.H4.abf"] >= 0.5 ~ "MODERATE colocalization (sensitive prior)",
      pp["PP.H3.abf"] >= 0.5  ~ "DISTINCT causal variants",
      TRUE ~ "INSUFFICIENT evidence"
    )
  )
}

## ═══════════════════════════════════════════════════════════════════════════════
## MAIN
## ═══════════════════════════════════════════════════════════════════════════════
results <- list()

for (protein in names(targets)) {
  cfg <- targets[[protein]]
  message(sprintf("\n═══ %s ═══", protein))

  ## Step 1: Extract pQTL region (GRCh38)
  pqtl_raw <- tryCatch(
    extract_pqtl_region(protein, cfg$chr, cfg$pos_hg38),
    error = function(e) { message("  Tabix ERROR: ", e$message); NULL }
  )
  if (is.null(pqtl_raw) || nrow(pqtl_raw) == 0) next
  message(sprintf("  %d pQTL variants extracted", nrow(pqtl_raw)))

  ## Save regional pQTL
  out_file <- file.path(proj, "data/pqtl/priority_regions",
                        paste0(protein, "_pqtl_regions.tsv.gz"))
  fwrite(pqtl_raw, out_file, sep = "\t")
  message(sprintf("  Saved → %s", out_file))

  ## Step 2: Load GWAS region (GRCh38 coords from pQTL)
  gwas_reg <- tryCatch(
    load_gwas_region(cfg$chr, min(pqtl_raw$pos), max(pqtl_raw$pos)),
    error = function(e) { message("  GWAS load ERROR: ", e$message); NULL }
  )
  if (is.null(gwas_reg) || nrow(gwas_reg) == 0) next

  ## Step 3: Harmonise (both GRCh38 — match on pos directly)
  harm <- tryCatch(
    harmonise_datasets(pqtl_raw, gwas_reg),
    error = function(e) { message("  Harmonise ERROR: ", e$message); NULL }
  )
  if (is.null(harm) || nrow(harm) < MIN_SNPS) {
    message(sprintf("  ✗ Insufficient harmonised SNPs (%d)", if(is.null(harm)) 0L else nrow(harm)))
    next
  }

  ## Step 4: COLOC
  res <- tryCatch(
    run_coloc(harm, protein),
    error = function(e) { message("  COLOC ERROR: ", e$message); NULL }
  )
  if (!is.null(res)) {
    results[[protein]] <- res
    message(sprintf("  ✓ PPH4 = %.4f — %s", res$PPH4, res$interpretation))
  }
}

## ── Save results ──────────────────────────────────────────────────────────────
if (length(results) > 0) {
  final <- bind_rows(results)
  out_path <- file.path(out_dir, "protein_coloc_snx15_pm20d1.csv")
  write_csv(final, out_path)
  message(sprintf("\n✓ Saved → %s", out_path))

  cat("\n═══ PROTEIN-SIDE COLOC RESULTS ═══\n")
  final |>
    select(protein, n_snps, PPH0, PPH1, PPH2, PPH3, PPH4,
           coloc_lead_snp, coloc_lead_PP_H4, interpretation) |>
    print(n = Inf)
} else {
  message("\n✗ No results produced")
}
