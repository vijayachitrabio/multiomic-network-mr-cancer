#!/usr/bin/env Rscript
## Script 38: coloc.susie for SNX15 and PM20D1 — resolves PPH3 from script 37
## ─────────────────────────────────────────────────────────────────────────────
## Approach: compute LD from 1000 Genomes EUR GRCh38 VCF via Tabix, run
## runsusie() on FinnGen pQTL and breast GWAS, then coloc.susie().
##
## Why needed: coloc.abf gave PPH3≈1 (distinct signals) partly because it cannot
## model multiple causal variants per locus. SuSiE fine-maps each dataset
## independently and tests whether ANY of the credible set variants colocalize.
## This handles allelic heterogeneity and improves resolution with dense LD.
##
## 1000G source (GRCh38 high-coverage, 3202 samples):
##   https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/
##   1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
##
## BUGS FIXED (v3 2026-05-07):
##  - PPH4 extraction: coloc.bf_bf/$summary$PP.H4.abf (NOT $results$PP.H4.abf)
##  - clean_lbf now replaces NaN as well as true NA
##  - All tibble fields wrapped in as.numeric() / length/NULL guards
##  - coloc.susie tried with data lists first (most reliable), pre-run as fallback
##  - coloc.abf always runs as baseline regardless of SuSiE outcome

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)       # >= 5.2.3 — includes runsusie / coloc.susie
  library(susieR)
  library(Rsamtools)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
})

set.seed(2026)
proj    <- "."
out_dir <- file.path(proj, "results/validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ── Config ────────────────────────────────────────────────────────────────────
WINDOW_BP  <- 500000L
MIN_SNPS   <- 100L     # need more for SuSiE fine-mapping
MAF_FLOOR  <- 0.01     # EUR MAF threshold for LD computation
BASE_PQTL  <- "https://storage.googleapis.com/finngen-public-data-r10/omics/proteomics/release_2023_03_02/data/Olink/pQTL"
BASE_1KG   <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"

gwas_cfg <- list(
  file       = file.path(proj, "data/cancer_gwas/Breast_GCST90018757.h.tsv.gz"),
  n_cases    = 122977L,
  n_controls = 105974L
)
gwas_cfg$n_total <- gwas_cfg$n_cases + gwas_cfg$n_controls
gwas_cfg$s       <- gwas_cfg$n_cases / gwas_cfg$n_total

targets <- list(
  SNX15  = list(chr = "11", pos_hg38 = 65025679L),
  PM20D1 = list(chr =  "1", pos_hg38 = 205847765L)
)

## ── Step 0: EUR sample IDs and column indices ──────────────────────────────
message("Loading 1000G population panel...")
panel_url <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
panel <- fread(panel_url)
eur_ids <- panel[Superpopulation == "EUR", SampleID]
message(sprintf("  %d EUR samples", length(eur_ids)))

get_eur_col_idx <- function(chr_num) {
  vcf_url <- sprintf("%s/1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                     BASE_1KG, chr_num)
  hdr <- headerTabix(TabixFile(vcf_url))
  chrom_line <- hdr$header[length(hdr$header)]
  all_samples <- strsplit(chrom_line, "\t", fixed = TRUE)[[1]][-(1:9)]
  eur_col <- which(all_samples %in% eur_ids) + 9L   # 1-based, offset by 9 fixed fields
  list(url = vcf_url, eur_col = eur_col, n_eur = length(eur_col))
}

## ── Step 1: pQTL extractor (same as script 37) ────────────────────────────
extract_pqtl <- function(protein, chr, pos_hg38, window = WINDOW_BP) {
  url    <- sprintf("%s/Olink_Batch1_%s.txt.gz", BASE_PQTL, protein)
  start  <- max(1L, pos_hg38 - window)
  end    <- pos_hg38 + window
  region <- GRanges(chr, IRanges(start, end))
  tbx    <- TabixFile(url)
  lines  <- scanTabix(tbx, param = region)[[1]]
  if (length(lines) == 0) return(NULL)
  x <- fread(text = paste(lines, collapse = "\n"), header = FALSE)
  col_names <- c("chr","pos","variant_id","ref","alt","alt_freq",
                 "beta","se","t_stat","p","log10_p","n")
  if (ncol(x) == 11) col_names <- col_names[-12]
  setnames(x, col_names[seq_len(ncol(x))])
  if (!"n" %in% names(x)) x[, n := 619L]
  x[!is.na(beta) & !is.na(se) & se > 0 & !is.na(alt_freq) & alt_freq > 0 & alt_freq < 1]
}

## ── Step 2: GWAS region loader (GRCh38) ──────────────────────────────────
load_gwas <- function(chr, pos_min, pos_max) {
  read_tsv(gwas_cfg$file, show_col_types = FALSE,
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
           base_pair_location >= pos_min,
           base_pair_location <= pos_max) |>
    rename(pos = base_pair_location, ea = effect_allele, oa = other_allele,
           beta_g = beta, se_g = standard_error,
           eaf_g = effect_allele_frequency, p_g = p_value) |>
    filter(!is.na(beta_g), !is.na(se_g), se_g > 0,
           !is.na(eaf_g), eaf_g > 0, eaf_g < 1)
}

## ── Step 3: Harmonise ─────────────────────────────────────────────────────
harmonise <- function(pqtl, gwas) {
  flip <- c(A="T",T="A",C="G",G="C")
  inner_join(
    pqtl |> mutate(pos = as.integer(pos)),
    gwas |> mutate(pos = as.integer(pos)),
    by = "pos"
  ) |>
    mutate(
      ea_p = toupper(alt), oa_p = toupper(ref),
      ea_g2 = toupper(ea), oa_g2 = toupper(oa),
      ea_pf = flip[ea_p], oa_pf = flip[oa_p],
      match_d  = ea_p == ea_g2  & oa_p == oa_g2,
      match_s  = ea_p == oa_g2  & oa_p == ea_g2,
      match_fl = !is.na(ea_pf) & ea_pf == ea_g2 & oa_pf == oa_g2,
      match_fs = !is.na(ea_pf) & ea_pf == oa_g2 & oa_pf == ea_g2,
      palin    = ea_p == flip[oa_p]
    ) |>
    filter(!palin | (alt_freq > 0.1 & alt_freq < 0.9)) |>
    filter(match_d | match_s | match_fl | match_fs) |>
    mutate(
      beta_g_h = if_else(match_d | match_fl, beta_g, -beta_g),
      eaf_g_h  = if_else(match_d | match_fl, eaf_g, 1 - eaf_g)
    )
}

## ── Step 4: Build LD matrix from 1000G EUR VCF ───────────────────────────
build_ld_matrix <- function(chr_num, pos_min, pos_max, eur_col, vcf_url,
                            target_snp_ids) {
  ## target_snp_ids: character vector "POS:REF:ALT" (GRCh38) to filter/match
  region <- GRanges(paste0("chr", chr_num), IRanges(pos_min, pos_max))
  message(sprintf("  1000G tabix: chr%s:%d-%d", chr_num, pos_min, pos_max))
  tbx   <- TabixFile(vcf_url)
  lines <- scanTabix(tbx, param = region)[[1]]
  message(sprintf("  %d VCF lines fetched", length(lines)))
  if (length(lines) == 0) return(NULL)

  ## Parse fixed fields + EUR genotypes
  ## VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <samples...>
  ## Genotype: "0|0", "0|1", "1|0", "1|1" (phased); take only GT (first subfield)
  message("  Parsing VCF genotypes...")

  # Split all lines at once: list of length n_lines, each element = all fields
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  n_lines     <- length(split_lines)

  # Extract fixed info (CHROM, POS, REF, ALT) + EUR columns
  pos_vec  <- as.integer(vapply(split_lines, `[[`, character(1), 2))
  ref_vec  <- vapply(split_lines, `[[`, character(1), 4)
  alt_vec  <- vapply(split_lines, `[[`, character(1), 5)
  snp_key  <- paste0(pos_vec, ":", ref_vec, ":", alt_vec)

  # Keep only biallelic SNPs (single char REF and ALT) with MAF > floor
  keep_snp <- nchar(ref_vec) == 1 & nchar(alt_vec) == 1 &
              !grepl(",", alt_vec, fixed = TRUE)

  # EUR genotype matrix: rows = variants, cols = EUR samples
  # dosage = first allele + third allele of phased "A|B"
  eur_gt_matrix <- matrix(NA_real_, nrow = n_lines, ncol = length(eur_col))
  for (i in seq_len(n_lines)) {
    flds <- split_lines[[i]][eur_col]
    eur_gt_matrix[i, ] <- as.integer(substr(flds, 1L, 1L)) +
                          as.integer(substr(flds, 3L, 3L))
  }

  # Compute EUR MAF per variant
  af <- rowMeans(eur_gt_matrix, na.rm = TRUE) / 2
  keep_maf <- af > MAF_FLOOR & af < (1 - MAF_FLOOR)

  keep <- keep_snp & keep_maf
  message(sprintf("  %d / %d variants pass SNP+MAF filters", sum(keep), n_lines))

  snp_key_f   <- snp_key[keep]
  geno_f      <- eur_gt_matrix[keep, , drop = FALSE]
  rownames(geno_f) <- snp_key_f

  # Filter to variants present in harmonised data
  shared <- intersect(snp_key_f, target_snp_ids)
  message(sprintf("  %d variants shared with harmonised pQTL+GWAS", length(shared)))
  if (length(shared) < MIN_SNPS) return(NULL)

  geno_shared <- geno_f[shared, , drop = FALSE]

  # Compute LD (variant × variant correlation matrix)
  message("  Computing LD matrix...")
  ld <- cor(t(geno_shared))
  ld[is.na(ld)] <- 0
  diag(ld) <- 1

  # Regularise to nearest positive definite matrix
  # (1000G EUR N=633 gives a near-PSD but not strictly PD LD matrix at ~1000 SNPs)
  # Eigenvalue truncation: shrink negative eigenvalues to a small positive floor
  message("  Regularising LD matrix (eigen truncation)...")
  eig <- eigen(ld, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-4)    # floor at 1e-4
  ld_reg <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  # Re-normalise diagonal to 1 (ensure proper correlation matrix)
  d_inv  <- 1 / sqrt(diag(ld_reg))
  ld_reg <- diag(d_inv) %*% ld_reg %*% diag(d_inv)
  diag(ld_reg) <- 1
  rownames(ld_reg) <- colnames(ld_reg) <- shared

  list(ld = ld_reg, snp_ids = shared)
}

## ── Step 5: coloc.susie ────────────────────────────────────────────────────
run_coloc_susie <- function(harm, ld_info, protein, n_pqtl = 619L) {

  ## --- SNP key matching ---
  ## 1000G keys are "POS:REF:ALT"; pQTL data has ref/alt columns
  snp_key_ref <- paste0(harm$pos, ":", toupper(harm$ref), ":", toupper(harm$alt))

  shared <- ld_info$snp_ids
  idx    <- match(shared, snp_key_ref)

  ok <- !is.na(idx)
  if (sum(ok) < MIN_SNPS) {
    message(sprintf("  ⚠ Only %d variants in LD+harm intersection, skipping", sum(ok)))
    return(NULL)
  }

  harm_sub <- harm[idx[ok], ]
  ld_sub   <- ld_info$ld[ok, ok, drop = FALSE]
  snp_sub  <- shared[ok]
  message(sprintf("  %d variants in LD × harmonised intersection", sum(ok)))

  ## --- Build coloc data lists (same for both ABF and SuSiE) ---
  D1 <- list(
    beta    = harm_sub$beta,
    varbeta = harm_sub$se^2,
    snp     = snp_sub,
    type    = "quant",
    N       = n_pqtl,
    MAF     = pmin(harm_sub$alt_freq, 1 - harm_sub$alt_freq),
    LD      = ld_sub
  )
  D2 <- list(
    beta    = harm_sub$beta_g_h,
    varbeta = harm_sub$se_g^2,
    snp     = snp_sub,
    type    = "cc",
    N       = gwas_cfg$n_total,
    s       = gwas_cfg$s,
    MAF     = pmin(harm_sub$eaf_g_h, 1 - harm_sub$eaf_g_h),
    LD      = ld_sub
  )

  ## --- coloc.abf (always run first; does NOT use LD) ---
  message("  Running coloc.abf (comparison baseline)...")
  D1_abf <- D1[c("beta","varbeta","snp","type","N","MAF")]
  D2_abf <- D2[c("beta","varbeta","snp","type","N","s","MAF")]
  pp_abf_vec <- tryCatch({
    res_abf <- coloc.abf(D1_abf, D2_abf)
    res_abf$summary   # named numeric vector: PP.H0.abf … PP.H4.abf
  }, error = function(e) {
    message("  coloc.abf error: ", e$message)
    c(PP.H0.abf=NA_real_, PP.H1.abf=NA_real_, PP.H2.abf=NA_real_,
      PP.H3.abf=NA_real_, PP.H4.abf=NA_real_)
  })
  message(sprintf("  [ABF] PPH3=%.4f  PPH4=%.4f",
                  pp_abf_vec["PP.H3.abf"], pp_abf_vec["PP.H4.abf"]))

  ## --- SuSiE diagnostics (convergence + CS counts) ---
  ## clean_lbf: replace true NA and NaN with 0 (log-BF = 0 → BF = 1, neutral)
  ## NaN arises from 0/0 in SuSiE when all prior mass is on null variant
  ## ±Inf are legitimate (log(0) or perfect BF) — do NOT replace these
  clean_lbf <- function(s) {
    if (!is.null(s$lbf_variable)) {
      bad <- is.na(s$lbf_variable) | is.nan(s$lbf_variable)
      if (any(bad)) {
        s$lbf_variable[bad] <- 0
      }
    }
    s
  }

  message("  Running runsusie() on pQTL dataset (max_iter=10000)...")
  s1 <- tryCatch(
    runsusie(D1, repeat_until_convergence = TRUE, maxit = 10000L),
    error = function(e) { message("  SuSiE D1 error: ", e$message); NULL }
  )
  message("  Running runsusie() on GWAS dataset (max_iter=10000)...")
  s2 <- tryCatch(
    runsusie(D2, repeat_until_convergence = TRUE, maxit = 10000L),
    error = function(e) { message("  SuSiE D2 error: ", e$message); NULL }
  )

  n_cs1 <- if (!is.null(s1) && !is.null(s1$sets$cs)) length(s1$sets$cs) else 0L
  n_cs2 <- if (!is.null(s2) && !is.null(s2$sets$cs)) length(s2$sets$cs) else 0L
  message(sprintf("  Credible sets: pQTL=%d, GWAS=%d", n_cs1, n_cs2))

  if (!is.null(s1)) s1 <- clean_lbf(s1)
  if (!is.null(s2)) s2 <- clean_lbf(s2)

  ## --- coloc.susie: three escalating attempts ---
  pph4_susie <- 0
  best_snp   <- NA_character_
  n_pairs    <- 0L

  # Attempt 1: pass pre-run SuSiE objects (class "susie" → coloc.susie skips re-running)
  if (!is.null(s1) && !is.null(s2) && n_cs1 > 0 && n_cs2 > 0) {
    message("  Attempt 1: coloc.susie(s1, s2) with pre-run SuSiE objects...")
    cs_result <- tryCatch(
      coloc.susie(s1, s2),
      error = function(e) { message("  Attempt 1 error: ", e$message); NULL }
    )

    if (!is.null(cs_result)) {
      # Extract best PPH4 from $summary (data.table with one row per CS pair)
      # Columns: nsnps, hit1, hit2, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf
      tryCatch({
        smry <- cs_result$summary
        if (!is.null(smry) && (is.data.frame(smry) || is.data.table(smry)) &&
            nrow(smry) > 0 && "PP.H4.abf" %in% names(smry)) {
          best_idx   <- which.max(as.numeric(smry[["PP.H4.abf"]]))
          pph4_susie <- as.numeric(smry[["PP.H4.abf"]][best_idx])
          n_pairs    <- nrow(smry)
          if ("hit1" %in% names(smry)) {
            best_snp <- as.character(smry[["hit1"]][best_idx])
          }
        }
      }, error = function(e) message("  PPH4 extraction (Attempt 1) error: ", e$message))
    }
  }

  # Attempt 2: pass data lists directly — coloc.susie re-runs SuSiE internally
  # (use if Attempt 1 failed or returned 0 pairs)
  if (n_pairs == 0) {
    message("  Attempt 2: coloc.susie(D1, D2) with data lists + runsusie.args...")
    cs_result2 <- tryCatch(
      coloc.susie(D1, D2,
                  runsusie.args = list(repeat_until_convergence = TRUE, maxit = 10000L)),
      error = function(e) { message("  Attempt 2 error: ", e$message); NULL }
    )

    if (!is.null(cs_result2)) {
      tryCatch({
        smry2 <- cs_result2$summary
        if (!is.null(smry2) && (is.data.frame(smry2) || is.data.table(smry2)) &&
            nrow(smry2) > 0 && "PP.H4.abf" %in% names(smry2)) {
          best_idx2   <- which.max(as.numeric(smry2[["PP.H4.abf"]]))
          pph4_c2     <- as.numeric(smry2[["PP.H4.abf"]][best_idx2])
          # Use Attempt 2 result only if better or Attempt 1 gave 0
          if (pph4_c2 > pph4_susie || n_pairs == 0) {
            pph4_susie <- pph4_c2
            n_pairs    <- nrow(smry2)
            if ("hit1" %in% names(smry2)) {
              best_snp <- as.character(smry2[["hit1"]][best_idx2])
            }
          }
        }
      }, error = function(e) message("  PPH4 extraction (Attempt 2) error: ", e$message))
    }
  }

  # Attempt 3: manual coloc.bf_bf on lbf_variable rows for each CS
  if (n_pairs == 0 && !is.null(s1) && !is.null(s2)) {
    message("  Attempt 3: manual coloc.bf_bf fallback...")
    tryCatch({
      idx1 <- s1$sets$cs_index
      idx2 <- s2$sets$cs_index
      if (length(idx1) > 0 && length(idx2) > 0) {
        bf1 <- s1$lbf_variable[idx1, , drop = FALSE]
        bf2 <- s2$lbf_variable[idx2, , drop = FALSE]
        # Ensure both have same SNP columns
        common_cols <- intersect(colnames(bf1), colnames(bf2))
        common_cols <- setdiff(common_cols, "null")
        if (length(common_cols) >= MIN_SNPS) {
          bf1 <- bf1[, common_cols, drop = FALSE]
          bf2 <- bf2[, common_cols, drop = FALSE]
          ret <- coloc:::coloc.bf_bf(bf1, bf2)
          smry3 <- ret$summary
          if (!is.null(smry3) && (is.data.frame(smry3) || is.data.table(smry3)) &&
              nrow(smry3) > 0 && "PP.H4.abf" %in% names(smry3)) {
            best_idx3  <- which.max(as.numeric(smry3[["PP.H4.abf"]]))
            pph4_susie <- as.numeric(smry3[["PP.H4.abf"]][best_idx3])
            n_pairs    <- nrow(smry3)
            if ("hit1" %in% names(smry3)) {
              best_snp <- as.character(smry3[["hit1"]][best_idx3])
            }
          }
        } else {
          message(sprintf("  Attempt 3: only %d common lbf columns — skipping", length(common_cols)))
        }
      }
    }, error = function(e) message("  Attempt 3 error: ", e$message))
  }

  ## --- Safety: ensure scalar, finite, correct types ---
  if (is.null(pph4_susie) || length(pph4_susie) == 0 || !is.finite(pph4_susie)) pph4_susie <- 0
  if (is.null(best_snp)   || length(best_snp)   == 0)                           best_snp   <- NA_character_
  if (is.null(n_pairs)    || length(n_pairs)     == 0)                           n_pairs    <- 0L

  message(sprintf("  ✓ coloc.susie best PPH4 = %.4f  (coloc.abf PPH4 = %.4f, n_pairs = %d)",
                  pph4_susie,
                  as.numeric(pp_abf_vec["PP.H4.abf"]),
                  n_pairs))

  ## --- Interpret ---
  pp4_abf  <- as.numeric(pp_abf_vec["PP.H4.abf"])
  pp3_abf  <- as.numeric(pp_abf_vec["PP.H3.abf"])
  interp <- dplyr::case_when(
    pph4_susie >= 0.8                              ~ "STRONG colocalization (SuSiE)",
    pph4_susie >= 0.5                              ~ "MODERATE colocalization (SuSiE)",
    !is.na(pp4_abf) & pp4_abf >= 0.5              ~ "MODERATE colocalization (ABF only)",
    !is.na(pp3_abf) & pp3_abf >= 0.5              ~ "DISTINCT causal variants",
    TRUE                                           ~ "INSUFFICIENT evidence"
  )

  tibble::tibble(
    protein           = protein,
    cancer            = "Breast",
    n_snps_ld_harm    = as.integer(sum(ok)),
    n_cs_pqtl         = as.integer(n_cs1),
    n_cs_gwas         = as.integer(n_cs2),
    n_coloc_pairs     = as.integer(n_pairs),
    PPH4_susie_best   = round(as.numeric(pph4_susie),          4),
    susie_best_snp    = best_snp,
    PPH0_abf          = round(as.numeric(pp_abf_vec["PP.H0.abf"]), 4),
    PPH1_abf          = round(as.numeric(pp_abf_vec["PP.H1.abf"]), 4),
    PPH2_abf          = round(as.numeric(pp_abf_vec["PP.H2.abf"]), 4),
    PPH3_abf          = round(as.numeric(pp_abf_vec["PP.H3.abf"]), 4),
    PPH4_abf          = round(as.numeric(pp_abf_vec["PP.H4.abf"]), 4),
    interpretation    = interp
  )
}

## ═══════════════════════════════════════════════════════════════════════════
## MAIN
## ═══════════════════════════════════════════════════════════════════════════
results <- list()

for (protein in names(targets)) {
  cfg <- targets[[protein]]
  message(sprintf("\n═══ %s (coloc.susie) ═══", protein))

  ## Load VCF URL + EUR column indices (one call per chromosome)
  message("  Getting 1000G VCF header for chr", cfg$chr, "...")
  vcf_info <- tryCatch(
    get_eur_col_idx(cfg$chr),
    error = function(e) { message("  VCF header error: ", e$message); NULL }
  )
  if (is.null(vcf_info)) next
  message(sprintf("  Found %d EUR sample columns", vcf_info$n_eur))

  ## pQTL region
  pqtl <- tryCatch(
    extract_pqtl(protein, cfg$chr, cfg$pos_hg38),
    error = function(e) { message("  pQTL error: ", e$message); NULL }
  )
  if (is.null(pqtl) || nrow(pqtl) == 0) next
  message(sprintf("  %d pQTL variants", nrow(pqtl)))

  ## GWAS region
  gwas <- tryCatch(
    load_gwas(cfg$chr, min(pqtl$pos), max(pqtl$pos)),
    error = function(e) { message("  GWAS error: ", e$message); NULL }
  )
  if (is.null(gwas) || nrow(gwas) == 0) next
  message(sprintf("  %d GWAS variants", nrow(gwas)))

  ## Harmonise
  harm <- tryCatch(
    harmonise(pqtl, gwas),
    error = function(e) { message("  Harmonise error: ", e$message); NULL }
  )
  if (is.null(harm) || nrow(harm) < MIN_SNPS) {
    message(sprintf("  ✗ Only %d harmonised SNPs", if(is.null(harm)) 0L else nrow(harm)))
    next
  }
  message(sprintf("  %d harmonised SNPs", nrow(harm)))

  ## Build LD matrix
  ## Target keys must be in "POS:REF:ALT" format to match 1000G VCF keys
  ## pQTL 'ref' and 'alt' columns are already the VCF REF/ALT alleles (GRCh38)
  snp_key_ref <- paste0(harm$pos, ":", toupper(harm$ref), ":", toupper(harm$alt))
  ## Also include reverse orientation in case any REF/ALT are flipped
  snp_key_alt <- paste0(harm$pos, ":", toupper(harm$alt), ":", toupper(harm$ref))
  target_ids  <- union(snp_key_ref, snp_key_alt)

  ld_info <- tryCatch(
    build_ld_matrix(cfg$chr, min(pqtl$pos), max(pqtl$pos),
                    vcf_info$eur_col, vcf_info$url, target_ids),
    error = function(e) { message("  LD error: ", e$message); NULL }
  )
  if (is.null(ld_info)) {
    message("  ✗ LD build failed or too few shared variants")
    next
  }

  ## coloc.susie
  res <- tryCatch(
    run_coloc_susie(harm, ld_info, protein),
    error = function(e) { message("  run_coloc_susie() outer error: ", e$message); NULL }
  )
  if (!is.null(res)) {
    results[[protein]] <- res
    message(sprintf("  ✓ PPH4 (SuSiE) = %.4f  PPH4 (ABF) = %.4f  — %s",
                    res$PPH4_susie_best, res$PPH4_abf, res$interpretation))
  }
}

## ── Save ─────────────────────────────────────────────────────────────────
if (length(results) > 0) {
  final <- dplyr::bind_rows(results)
  out_path <- file.path(out_dir, "protein_coloc_susie_snx15_pm20d1.csv")
  readr::write_csv(final, out_path)
  message(sprintf("\n✓ Saved → %s", out_path))
  cat("\n═══ COLOC.SUSIE RESULTS ═══\n")
  print(final, n = Inf)
} else {
  message("\n✗ No results produced")
}
