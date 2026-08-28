#!/usr/bin/env Rscript
## Script 84: PAIRED SuSiE control -- as-published vs palindrome-dropped
## ─────────────────────────────────────────────────────────────────────────────
## WHY THIS IS NECESSARY
## Script 82 ran the SuSiE arm only in the palindrome-dropped configuration and
## returned EFNA1 PPH4 = NA (zero GWAS credible sets) and ATRAID PPH4 = 0.0015,
## against published values of 0.963 and 0.996. That difference has TWO possible
## causes and script 82 cannot separate them:
##
##   (a) removing ambiguous palindromic variants genuinely destroys the shared
##       credible-set signal at these two loci, or
##   (b) this reimplementation does not reproduce the original SuSiE arm at all
##       (LD construction, variant ordering, convergence path), in which case the
##       comparison is uninformative rather than alarming.
##
## The ABF arm (script 83) settled this for itself by running both configurations
## and reproducing the published ABF posteriors to 3 decimals. This script does
## the same for SuSiE: identical pipeline, identical LD matrix, run TWICE per
## locus, differing ONLY in the palindromic filter. Anything that shows up in the
## difference is then attributable to the filter and nothing else.
##
## Scope: EFNA1 and ATRAID -- the two proteins whose Tier 1 assignment rests on
## SuSiE rather than ABF. TNFRSF6B is included as a positive control because its
## palindrome-dropped SuSiE value (0.900) already sits close to published (0.885).
##
## Output: results/validation/protein_coloc_susie_palindromic_paired.csv

suppressPackageStartupMessages({
  library(data.table); library(coloc); library(susieR)
  library(Rsamtools); library(GenomicRanges); library(dplyr); library(readr)
})
set.seed(2026)
proj <- "."
out_dir <- file.path(proj, "results/validation")
WINDOW_BP <- 500000L
MIN_SNPS  <- 50L
BASE_1KG <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
`%||%` <- function(a, b) if (is.null(a)) b else a

GW <- list(file = file.path(proj, "data/cancer_gwas/Breast_GCST90018757.h.tsv.gz"),
           n_cases = 122977L, n_controls = 105974L)
GW$n_total <- GW$n_cases + GW$n_controls
GW$s <- GW$n_cases / GW$n_total

TARGETS <- c("EFNA1", "ATRAID", "TNFRSF6B")

message("Loading 1000G EUR panel...")
panel <- fread("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt")
eur_ids <- panel[Superpopulation == "EUR", SampleID]

get_vcf_info <- function(chr_num) {
  url <- sprintf("%s/1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                 BASE_1KG, chr_num)
  hdr <- headerTabix(TabixFile(url))
  samples <- strsplit(hdr$header[length(hdr$header)], "\t", fixed = TRUE)[[1]][-(1:9)]
  list(url = url, eur_col = which(samples %in% eur_ids) + 9L)
}

load_pqtl <- function(p) {
  x <- fread(file.path(proj, "data/pqtl/priority_regions", paste0(p, "_pqtl_regions.tsv.gz")))
  x[!is.na(beta) & !is.na(se) & se > 0 & !is.na(alt_freq) & alt_freq > 0 & alt_freq < 1]
}

gwas_all <- NULL
get_gwas_region <- function(chr, lo, hi) {
  if (is.null(gwas_all)) {
    message("  reading breast GWAS (once)...")
    g <- read_tsv(GW$file, show_col_types = FALSE,
          col_types = cols_only(chromosome = col_integer(),
            base_pair_location = col_integer(), effect_allele = col_character(),
            other_allele = col_character(), beta = col_double(),
            standard_error = col_double(), effect_allele_frequency = col_double(),
            p_value = col_double())) |>
      rename(pos = base_pair_location, ea = effect_allele, oa = other_allele,
             beta_g = beta, se_g = standard_error,
             eaf_g = effect_allele_frequency, p_g = p_value) |>
      filter(!is.na(beta_g), !is.na(se_g), se_g > 0, !is.na(eaf_g), eaf_g > 0, eaf_g < 1)
    gwas_all <<- as.data.table(g)
  }
  gwas_all[chromosome == as.integer(chr) & pos >= lo & pos <= hi]
}

harmonise <- function(pqtl, gwas, drop_palin) {
  flip <- c(A = "T", T = "A", C = "G", G = "C")
  h <- inner_join(pqtl |> mutate(pos = as.integer(pos)),
                  gwas |> mutate(pos = as.integer(pos)), by = "pos") |>
    mutate(ea_p = toupper(alt), oa_p = toupper(ref),
           ea_g2 = toupper(ea), oa_g2 = toupper(oa),
           ea_pf = flip[toupper(alt)], oa_pf = flip[toupper(ref)],
           match_d  = ea_p == ea_g2 & oa_p == oa_g2,
           match_s  = ea_p == oa_g2 & oa_p == ea_g2,
           match_fl = !is.na(ea_pf) & ea_pf == ea_g2 & oa_pf == oa_g2,
           match_fs = !is.na(ea_pf) & ea_pf == oa_g2 & oa_pf == ea_g2,
           palin    = ea_p == flip[oa_p])
  h <- if (drop_palin) filter(h, !palin)
       else filter(h, !palin | (alt_freq > 0.1 & alt_freq < 0.9))
  h |> filter(match_d | match_s | match_fl | match_fs) |>
    mutate(beta_g_h = if_else(match_d | match_fl, beta_g, -beta_g),
           eaf_g_h  = if_else(match_d | match_fl, eaf_g, 1 - eaf_g))
}

build_ld <- function(chr, lo, hi, eur_col, url, want) {
  lines <- tryCatch(scanTabix(TabixFile(url),
                    param = GRanges(paste0("chr", chr), IRanges(lo, hi)))[[1]],
                    error = function(e) character(0))
  if (!length(lines)) return(NULL)
  ids <- character(0); dos <- list()
  for (ln in lines) {
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(f) < 10) next
    id <- paste0(f[2], ":", toupper(f[4]), ":", toupper(f[5]))
    if (!(id %in% want)) next
    gt <- f[eur_col]
    d <- vapply(gt, function(g) {
      a <- substr(g, 1, 1); b <- substr(g, 3, 3)
      if (a == "." || b == ".") return(NA_real_)
      as.numeric(a) + as.numeric(b) }, numeric(1), USE.NAMES = FALSE)
    if (mean(is.na(d)) > 0.05) next
    if (length(unique(d[!is.na(d)])) < 2) next
    ids <- c(ids, id); dos[[length(dos) + 1]] <- d
  }
  if (length(ids) < 2) return(NULL)
  M <- do.call(cbind, dos); colnames(M) <- ids
  M <- M[, !duplicated(colnames(M)), drop = FALSE]
  ld <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))
  ld[is.na(ld)] <- 0; diag(ld) <- 1
  list(ld = ld, ids = colnames(M))
}

clean_lbf <- function(s) {
  if (!is.null(s$lbf_variable)) {
    b <- is.na(s$lbf_variable) | is.nan(s$lbf_variable)
    if (any(b)) s$lbf_variable[b] <- 0
  }; s
}
best_pph4 <- function(sm) {
  if (is.null(sm)) return(list(p = NA_real_, snp = NA_character_, n = 0L))
  if (is.numeric(sm) && "PP.H4.abf" %in% names(sm))
    return(list(p = as.numeric(sm["PP.H4.abf"]), snp = NA_character_, n = 1L))
  if ((is.data.frame(sm) || is.data.table(sm)) && nrow(sm) > 0 && "PP.H4.abf" %in% names(sm)) {
    i <- which.max(as.numeric(sm[["PP.H4.abf"]]))
    return(list(p = as.numeric(sm[["PP.H4.abf"]][i]),
                snp = if ("hit1" %in% names(sm)) as.character(sm[["hit1"]][i]) else NA_character_,
                n = nrow(sm)))
  }
  list(p = NA_real_, snp = NA_character_, n = 0L)
}

## one SuSiE colocalization run on an already-harmonised, already-LD-matched set
susie_run <- function(h, ld, ids) {
  D1 <- list(beta = h$beta, varbeta = h$se^2, snp = ids, type = "quant",
             N = 619L, MAF = pmin(h$alt_freq, 1 - h$alt_freq), LD = ld)
  D2 <- list(beta = h$beta_g_h, varbeta = h$se_g^2, snp = ids, type = "cc",
             N = GW$n_total, s = GW$s,
             MAF = pmin(h$eaf_g_h, 1 - h$eaf_g_h), LD = ld)
  s1 <- tryCatch(runsusie(D1, repeat_until_convergence = TRUE, maxit = 10000L),
                 error = function(e) NULL)
  s2 <- tryCatch(runsusie(D2, repeat_until_convergence = TRUE, maxit = 10000L),
                 error = function(e) NULL)
  n1 <- if (!is.null(s1)) { s1 <- clean_lbf(s1); length(s1$sets$cs %||% list()) } else 0L
  n2 <- if (!is.null(s2)) { s2 <- clean_lbf(s2); length(s2$sets$cs %||% list()) } else 0L
  out <- list(pph4 = NA_real_, snp = NA_character_, pairs = 0L, n_cs1 = n1, n_cs2 = n2)
  if (n1 > 0 && n2 > 0) {
    csr <- tryCatch(coloc.susie(s1, s2), error = function(e) NULL)
    if (is.null(csr)) csr <- tryCatch(coloc.susie(D1, D2,
        runsusie.args = list(repeat_until_convergence = TRUE, maxit = 10000L)),
        error = function(e) NULL)
    if (is.null(csr)) csr <- tryCatch({
      b1 <- s1$lbf_variable[s1$sets$cs_index, , drop = FALSE]
      b2 <- s2$lbf_variable[s2$sets$cs_index, , drop = FALSE]
      b1[is.nan(b1)] <- 0; b2[is.nan(b2)] <- 0
      coloc:::coloc.bf_bf(b1, b2) }, error = function(e) NULL)
    if (!is.null(csr)) {
      e <- best_pph4(csr$summary)
      out$pph4 <- e$p; out$snp <- e$snp; out$pairs <- e$n
    }
  }
  out
}

res <- list()
for (prot in TARGETS) {
  message(sprintf("\n========== %s ==========", prot))
  pq <- load_pqtl(prot)
  lead <- pq[which.min(p)]
  chr <- as.character(lead$chr[1]); lp <- as.integer(lead$pos[1])
  lo <- lp - WINDOW_BP; hi <- lp + WINDOW_BP
  gw <- get_gwas_region(chr, lo, hi)
  pqw <- pq[pos >= lo & pos <= hi]

  vi <- get_vcf_info(chr)

  for (cfg in c("as_published", "no_palindromic")) {
    h <- harmonise(pqw, gw, drop_palin = (cfg == "no_palindromic"))
    key <- paste0(h$pos, ":", toupper(h$ref), ":", toupper(h$alt))
    ldi <- build_ld(chr, lo, hi, vi$eur_col, vi$url, key)
    if (is.null(ldi)) { message(sprintf("  [%s] no LD", cfg)); next }
    idx <- match(ldi$ids, key); ok <- !is.na(idx)
    hs <- h[idx[ok], ]; lds <- ldi$ld[ok, ok, drop = FALSE]; ids <- ldi$ids[ok]
    if (nrow(hs) < MIN_SNPS) { message(sprintf("  [%s] too few", cfg)); next }
    message(sprintf("  [%s] %d LD-matched variants; running SuSiE...", cfg, nrow(hs)))
    r <- susie_run(hs, lds, ids)
    message(sprintf("  [%s] CS pQTL=%d GWAS=%d  PPH4=%s  pairs=%d", cfg,
                    r$n_cs1, r$n_cs2,
                    ifelse(is.na(r$pph4), "NA", sprintf("%.4f", r$pph4)), r$pairs))
    res[[paste(prot, cfg)]] <- data.table(
      protein = prot, config = cfg, n_ld_matched = nrow(hs),
      n_cs_pqtl = r$n_cs1, n_cs_gwas = r$n_cs2, n_coloc_pairs = r$pairs,
      PPH4_susie = r$pph4, susie_best_snp = r$snp)
    fwrite(rbindlist(res, fill = TRUE),
           file.path(out_dir, "protein_coloc_susie_palindromic_paired.csv"))
  }
}
out <- rbindlist(res, fill = TRUE)
fwrite(out, file.path(out_dir, "protein_coloc_susie_palindromic_paired.csv"))
message("\n=== PAIRED SuSiE RESULT ===")
print(out)
