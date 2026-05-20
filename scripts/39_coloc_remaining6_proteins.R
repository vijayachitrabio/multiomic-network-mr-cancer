#!/usr/bin/env Rscript
## Script 39: coloc.abf + coloc.susie for remaining 6 Phase-2 protein hits
## ─────────────────────────────────────────────────────────────────────────
## Proteins: IL34, EFNA1, TNFRSF6B, APOE, ATRAID, ITIH3
## pQTL data already cached in data/pqtl/priority_regions/
## Same pipeline as scripts 37 (abf) + 38 (susie); GRCh38 throughout.
## Run with optional arg: Rscript 39_... PROTEIN  (to do one at a time)

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
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

## ── Config ────────────────────────────────────────────────────────────────
WINDOW_BP <- 500000L
MIN_SNPS  <- 50L
MAF_FLOOR <- 0.01
BASE_1KG  <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"

gwas_cfg <- list(
  file       = file.path(proj, "data/cancer_gwas/Breast_GCST90018757.h.tsv.gz"),
  n_cases    = 122977L,
  n_controls = 105974L
)
gwas_cfg$n_total <- gwas_cfg$n_cases + gwas_cfg$n_controls
gwas_cfg$s       <- gwas_cfg$n_cases / gwas_cfg$n_total

## Lead positions (GRCh38) from cached pQTL files (lowest p-value per region)
targets <- list(
  IL34     = list(chr = "16", pos_hg38 = 70660097L),
  EFNA1    = list(chr =  "1", pos_hg38 = 155158049L),
  TNFRSF6B = list(chr = "20", pos_hg38 = 63706054L),
  APOE     = list(chr = "19", pos_hg38 = 44908684L),
  ATRAID   = list(chr =  "2", pos_hg38 = 27163931L),
  ITIH3    = list(chr =  "3", pos_hg38 = 52775509L)
)

## Optional: run single protein from command line
## Rscript 39_... IL34
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1 && args[1] %in% names(targets)) {
  targets <- targets[args[1]]
  message(sprintf("Running single protein: %s", args[1]))
}

## ── EUR sample index (cached per chromosome) ──────────────────────────────
message("Loading 1000G population panel...")
panel_url <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
panel   <- fread(panel_url)
eur_ids <- panel[Superpopulation == "EUR", SampleID]
message(sprintf("  %d EUR samples", length(eur_ids)))

vcf_cache <- list()   # reuse VCF header across proteins on same chromosome

get_vcf_info <- function(chr_num) {
  key <- as.character(chr_num)
  if (!is.null(vcf_cache[[key]])) return(vcf_cache[[key]])
  vcf_url <- sprintf("%s/1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                     BASE_1KG, chr_num)
  hdr <- headerTabix(TabixFile(vcf_url))
  chrom_line  <- hdr$header[length(hdr$header)]
  all_samples <- strsplit(chrom_line, "\t", fixed = TRUE)[[1]][-(1:9)]
  eur_col     <- which(all_samples %in% eur_ids) + 9L
  info <- list(url = vcf_url, eur_col = eur_col, n_eur = length(eur_col))
  vcf_cache[[key]] <<- info
  info
}

## ── pQTL loader (from disk cache) ────────────────────────────────────────
load_pqtl <- function(protein) {
  f <- file.path(proj, "data/pqtl/priority_regions",
                 paste0(protein, "_pqtl_regions.tsv.gz"))
  x <- fread(f)
  x[!is.na(beta) & !is.na(se) & se > 0 & !is.na(alt_freq) & alt_freq > 0 & alt_freq < 1]
}

## ── GWAS region loader ────────────────────────────────────────────────────
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

## ── Harmonise ─────────────────────────────────────────────────────────────
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
      ea_pf = flip[ea_p],  oa_pf = flip[oa_p],
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

## ── LD matrix builder ─────────────────────────────────────────────────────
build_ld_matrix <- function(chr_num, pos_min, pos_max, eur_col, vcf_url,
                            target_snp_ids) {
  region <- GRanges(paste0("chr", chr_num), IRanges(pos_min, pos_max))
  message(sprintf("  1000G tabix: chr%s:%d-%d", chr_num, pos_min, pos_max))
  lines <- scanTabix(TabixFile(vcf_url), param = region)[[1]]
  message(sprintf("  %d VCF lines fetched", length(lines)))
  if (length(lines) == 0) return(NULL)

  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  n_lines     <- length(split_lines)
  pos_vec     <- as.integer(vapply(split_lines, `[[`, character(1), 2))
  ref_vec     <- vapply(split_lines, `[[`, character(1), 4)
  alt_vec     <- vapply(split_lines, `[[`, character(1), 5)
  snp_key     <- paste0(pos_vec, ":", ref_vec, ":", alt_vec)

  keep_snp <- nchar(ref_vec) == 1 & nchar(alt_vec) == 1 &
              !grepl(",", alt_vec, fixed = TRUE)

  message("  Parsing genotypes...")
  geno <- matrix(NA_real_, nrow = n_lines, ncol = length(eur_col))
  for (i in seq_len(n_lines)) {
    flds   <- split_lines[[i]][eur_col]
    geno[i, ] <- as.integer(substr(flds, 1L, 1L)) + as.integer(substr(flds, 3L, 3L))
  }

  af       <- rowMeans(geno, na.rm = TRUE) / 2
  keep_maf <- af > MAF_FLOOR & af < (1 - MAF_FLOOR)
  keep     <- keep_snp & keep_maf
  message(sprintf("  %d / %d pass SNP+MAF", sum(keep), n_lines))

  snp_key_f        <- snp_key[keep]
  geno_f           <- geno[keep, , drop = FALSE]
  rownames(geno_f) <- snp_key_f

  shared <- intersect(snp_key_f, target_snp_ids)
  message(sprintf("  %d shared with harmonised data", length(shared)))
  if (length(shared) < MIN_SNPS) return(NULL)

  geno_s <- geno_f[shared, , drop = FALSE]

  message("  Computing + regularising LD...")
  ld <- cor(t(geno_s)); ld[is.na(ld)] <- 0; diag(ld) <- 1
  eig <- eigen(ld, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-4)
  ld_reg <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  d_inv  <- 1 / sqrt(diag(ld_reg))
  ld_reg <- diag(d_inv) %*% ld_reg %*% diag(d_inv)
  diag(ld_reg) <- 1
  rownames(ld_reg) <- colnames(ld_reg) <- shared

  list(ld = ld_reg, snp_ids = shared)
}

## ── coloc.abf + coloc.susie ───────────────────────────────────────────────
clean_lbf <- function(s) {
  if (!is.null(s$lbf_variable)) {
    bad <- is.na(s$lbf_variable) | is.nan(s$lbf_variable)
    if (any(bad)) s$lbf_variable[bad] <- 0
  }
  s
}

extract_pph4 <- function(smry) {
  ## smry = $summary from coloc.susie/coloc.bf_bf (data.table with PP.H4.abf col)
  if (is.null(smry)) return(list(pph4 = 0, snp = NA_character_, n = 0L))
  if (is.numeric(smry) && "PP.H4.abf" %in% names(smry))
    return(list(pph4 = as.numeric(smry["PP.H4.abf"]), snp = NA_character_, n = 1L))
  if ((is.data.frame(smry) || is.data.table(smry)) &&
      nrow(smry) > 0 && "PP.H4.abf" %in% names(smry)) {
    best <- which.max(as.numeric(smry[["PP.H4.abf"]]))
    return(list(
      pph4 = as.numeric(smry[["PP.H4.abf"]][best]),
      snp  = if ("hit1" %in% names(smry)) as.character(smry[["hit1"]][best]) else NA_character_,
      n    = nrow(smry)
    ))
  }
  list(pph4 = 0, snp = NA_character_, n = 0L)
}

run_coloc <- function(harm, ld_info, protein, n_pqtl = 619L) {

  ## LD match
  snp_key_ref <- paste0(harm$pos, ":", toupper(harm$ref), ":", toupper(harm$alt))
  shared      <- ld_info$snp_ids
  idx         <- match(shared, snp_key_ref)
  ok          <- !is.na(idx)
  n_ok        <- sum(ok)

  if (n_ok < MIN_SNPS) {
    message(sprintf("  ⚠ Only %d variants in LD ∩ harm — skipping SuSiE (will do ABF only)", n_ok))
  }

  harm_sub <- harm[idx[ok], ]
  ld_sub   <- ld_info$ld[ok, ok, drop = FALSE]
  snp_sub  <- shared[ok]

  ## Build data lists
  make_d <- function(use_ld = TRUE) {
    D1 <- list(beta=harm_sub$beta, varbeta=harm_sub$se^2, snp=snp_sub,
               type="quant", N=n_pqtl,
               MAF=pmin(harm_sub$alt_freq, 1-harm_sub$alt_freq))
    D2 <- list(beta=harm_sub$beta_g_h, varbeta=harm_sub$se_g^2, snp=snp_sub,
               type="cc", N=gwas_cfg$n_total, s=gwas_cfg$s,
               MAF=pmin(harm_sub$eaf_g_h, 1-harm_sub$eaf_g_h))
    if (use_ld) { D1$LD <- ld_sub; D2$LD <- ld_sub }
    list(D1=D1, D2=D2)
  }

  ## --- coloc.abf (no LD needed) ---
  message("  coloc.abf...")
  ds      <- make_d(use_ld = FALSE)
  pp_abf  <- tryCatch({
    coloc.abf(ds$D1, ds$D2)$summary
  }, error = function(e) {
    message("  ABF error: ", e$message)
    c(PP.H0.abf=NA_real_, PP.H1.abf=NA_real_, PP.H2.abf=NA_real_,
      PP.H3.abf=NA_real_, PP.H4.abf=NA_real_)
  })
  message(sprintf("  [ABF] PPH3=%.4f  PPH4=%.4f",
                  pp_abf["PP.H3.abf"], pp_abf["PP.H4.abf"]))

  ## --- coloc.susie (needs LD, skip if too few variants) ---
  pph4_s <- 0; best_snp <- NA_character_; n_pairs <- 0L
  n_cs1  <- 0L; n_cs2 <- 0L

  if (n_ok >= MIN_SNPS) {
    ds_ld <- make_d(use_ld = TRUE)

    message("  runsusie pQTL...")
    s1 <- tryCatch(runsusie(ds_ld$D1, repeat_until_convergence=TRUE, maxit=10000L),
                   error=function(e){message("  D1 err: ",e$message);NULL})
    message("  runsusie GWAS...")
    s2 <- tryCatch(runsusie(ds_ld$D2, repeat_until_convergence=TRUE, maxit=10000L),
                   error=function(e){message("  D2 err: ",e$message);NULL})

    if (!is.null(s1)) { s1 <- clean_lbf(s1); n_cs1 <- length(s1$sets$cs %||% list()) }
    if (!is.null(s2)) { s2 <- clean_lbf(s2); n_cs2 <- length(s2$sets$cs %||% list()) }
    message(sprintf("  CS: pQTL=%d  GWAS=%d", n_cs1, n_cs2))

    if (!is.null(s1) && !is.null(s2) && n_cs1 > 0 && n_cs2 > 0) {
      ## Attempt 1: pre-run objects
      csr <- tryCatch(coloc.susie(s1, s2),
                      error=function(e){message("  susie err: ",e$message);NULL})
      ## Attempt 2: data lists
      if (is.null(csr)) {
        csr <- tryCatch(
          coloc.susie(ds_ld$D1, ds_ld$D2,
                      runsusie.args=list(repeat_until_convergence=TRUE, maxit=10000L)),
          error=function(e){message("  susie(D) err: ",e$message);NULL})
      }
      ## Attempt 3: coloc.bf_bf
      if (is.null(csr)) {
        csr <- tryCatch({
          idx1 <- s1$sets$cs_index; idx2 <- s2$sets$cs_index
          bf1  <- s1$lbf_variable[idx1,,drop=FALSE]
          bf2  <- s2$lbf_variable[idx2,,drop=FALSE]
          bf1[is.nan(bf1)] <- 0; bf2[is.nan(bf2)] <- 0
          coloc:::coloc.bf_bf(bf1, bf2)
        }, error=function(e){message("  bf_bf err: ",e$message);NULL})
      }

      if (!is.null(csr)) {
        ex <- extract_pph4(csr$summary)
        pph4_s   <- ex$pph4;  best_snp <- ex$snp;  n_pairs <- ex$n
      }
    }
  }

  ## Safety
  if (!is.finite(pph4_s) || is.null(pph4_s) || length(pph4_s)==0) pph4_s <- 0
  if (is.null(best_snp) || length(best_snp)==0) best_snp <- NA_character_

  pp4 <- as.numeric(pp_abf["PP.H4.abf"])
  pp3 <- as.numeric(pp_abf["PP.H3.abf"])
  interp <- case_when(
    pph4_s >= 0.8                       ~ "STRONG coloc (SuSiE)",
    pph4_s >= 0.5                       ~ "MODERATE coloc (SuSiE)",
    !is.na(pp4) & pp4 >= 0.8           ~ "STRONG coloc (ABF only)",
    !is.na(pp4) & pp4 >= 0.5           ~ "MODERATE coloc (ABF only)",
    !is.na(pp3) & pp3 >= 0.5           ~ "DISTINCT causal variants",
    TRUE                                ~ "INSUFFICIENT evidence"
  )

  message(sprintf("  → PPH4 SuSiE=%.4f  ABF=%.4f  [%s]", pph4_s, pp4, interp))

  tibble::tibble(
    protein        = protein,
    cancer         = "Breast",
    n_harm         = nrow(harm),
    n_ld_harm      = as.integer(n_ok),
    n_cs_pqtl      = as.integer(n_cs1),
    n_cs_gwas      = as.integer(n_cs2),
    n_coloc_pairs  = as.integer(n_pairs),
    PPH4_susie     = round(as.numeric(pph4_s), 4),
    susie_best_snp = best_snp,
    PPH0_abf       = round(as.numeric(pp_abf["PP.H0.abf"]), 4),
    PPH1_abf       = round(as.numeric(pp_abf["PP.H1.abf"]), 4),
    PPH2_abf       = round(as.numeric(pp_abf["PP.H2.abf"]), 4),
    PPH3_abf       = round(as.numeric(pp_abf["PP.H3.abf"]), 4),
    PPH4_abf       = round(as.numeric(pp_abf["PP.H4.abf"]), 4),
    interpretation = interp
  )
}

## ── null-coalescing helper ─────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

## ═══════════════════════════════════════════════════════════════════════════
## MAIN
## ═══════════════════════════════════════════════════════════════════════════
results <- list()
out_path <- file.path(out_dir, "protein_coloc_remaining6.csv")

for (protein in names(targets)) {
  cfg <- targets[[protein]]
  message(sprintf("\n══════════════════════════════════════\n  %s\n══════════════════════════════════════", protein))

  ## VCF header (cached per chromosome)
  vcf_info <- tryCatch(get_vcf_info(cfg$chr),
                       error=function(e){message("VCF err: ",e$message);NULL})
  if (is.null(vcf_info)) { message("  ✗ skipping"); next }
  message(sprintf("  chr%s: %d EUR columns", cfg$chr, vcf_info$n_eur))

  ## pQTL (from disk)
  pqtl <- tryCatch(load_pqtl(protein),
                   error=function(e){message("pQTL err: ",e$message);NULL})
  if (is.null(pqtl) || nrow(pqtl)==0) { message("  ✗ no pQTL"); next }
  message(sprintf("  %d pQTL variants (cached)", nrow(pqtl)))

  ## GWAS
  gwas <- tryCatch(load_gwas(cfg$chr, min(pqtl$pos), max(pqtl$pos)),
                   error=function(e){message("GWAS err: ",e$message);NULL})
  if (is.null(gwas) || nrow(gwas)==0) { message("  ✗ no GWAS"); next }
  message(sprintf("  %d GWAS variants", nrow(gwas)))

  ## Harmonise
  harm <- tryCatch(harmonise(pqtl, gwas),
                   error=function(e){message("Harm err: ",e$message);NULL})
  if (is.null(harm) || nrow(harm) < MIN_SNPS) {
    message(sprintf("  ✗ only %d harmonised", if(is.null(harm)) 0L else nrow(harm))); next
  }
  message(sprintf("  %d harmonised SNPs", nrow(harm)))

  ## LD
  target_ids <- union(
    paste0(harm$pos, ":", toupper(harm$ref), ":", toupper(harm$alt)),
    paste0(harm$pos, ":", toupper(harm$alt), ":", toupper(harm$ref))
  )
  ld_info <- tryCatch(
    build_ld_matrix(cfg$chr, min(pqtl$pos), max(pqtl$pos),
                    vcf_info$eur_col, vcf_info$url, target_ids),
    error=function(e){message("LD err: ",e$message);NULL}
  )
  if (is.null(ld_info)) { message("  ✗ LD failed"); next }

  ## Coloc
  res <- tryCatch(run_coloc(harm, ld_info, protein),
                  error=function(e){message("Coloc err: ",e$message);NULL})
  if (!is.null(res)) {
    results[[protein]] <- res

    ## Append to CSV immediately (so partial results survive if script crashes)
    if (file.exists(out_path)) {
      readr::write_csv(res, out_path, append = TRUE)
    } else {
      readr::write_csv(res, out_path)
    }
    message(sprintf("  ✓ saved to %s", basename(out_path)))
  }
}

## ── Final summary ─────────────────────────────────────────────────────────
if (length(results) > 0) {
  final <- dplyr::bind_rows(results)
  cat("\n╔══════════════════════════════════════╗\n")
  cat("  COLOC RESULTS — REMAINING 6 PROTEINS\n")
  cat("╚══════════════════════════════════════╝\n")
  print(final |> select(protein, n_harm, n_ld_harm, n_cs_pqtl, n_cs_gwas,
                         n_coloc_pairs, PPH4_susie, PPH3_abf, PPH4_abf, interpretation),
        n = Inf)
} else {
  message("✗ No results produced")
}
