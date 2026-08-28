#!/usr/bin/env Rscript
## Script 86: FAITHFUL SuSiE reproduction + LD-stability sweep
## ─────────────────────────────────────────────────────────────────────────────
## CORRECTION TO SCRIPTS 82/84/85
## Those three scripts rebuilt the LD matrix but omitted three steps that the
## original pipeline (scripts 38/39, build_ld_matrix) actually performs:
##
##   1. biallelic-SNV restriction: nchar(ref)==1 & nchar(alt)==1 & no "," in alt
##   2. a 1000G EUR allele-frequency filter, MAF_FLOOR = 0.01 < af < 0.99
##      -- note this is on the REFERENCE panel frequency, not the summary-stat MAF
##   3. eigenvalue regularisation of the LD matrix: floor eigenvalues at 1e-4,
##      reconstruct, then rescale to a correlation matrix
##
## Step 3 is the important one. runsusie() is very sensitive to the conditioning
## of the supplied LD matrix, and feeding it a near-singular correlation matrix
## produces exactly the pathologies seen in scripts 84/85: zero credible sets on
## the GWAS side (EFNA1) or a credible-set split that collapses pairwise PPH4
## (ATRAID). Those runs therefore understated the published result because of an
## omission in the reimplementation, NOT because of anything in the manuscript.
##
## This script replicates all three steps and then re-asks the stability question
## on top of a faithful baseline.
##
## Configurations:
##   faithful      exact replication of scripts 38/39
##   no_palin      faithful + all palindromic variants dropped
##   eurboot1..3   faithful, LD from three independent 80% EUR subsamples
##   no_regularise faithful but WITHOUT step 3, to demonstrate its effect
##
## Output: results/validation/susie_faithful_ld_stability.csv

suppressPackageStartupMessages({
  library(data.table); library(coloc); library(susieR)
  library(Rsamtools); library(GenomicRanges); library(dplyr); library(readr)
})
set.seed(2026)
proj <- "."
out_dir <- file.path(proj, "results/validation")
OUT <- file.path(out_dir, "susie_faithful_ld_stability.csv")

WINDOW_BP <- 500000L
MIN_SNPS  <- 50L
MAF_FLOOR <- 0.01
BASE_1KG <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
`%||%` <- function(a, b) if (is.null(a)) b else a

GW <- list(file = file.path(proj, "data/cancer_gwas/Breast_GCST90018757.h.tsv.gz"),
           n_cases = 122977L, n_controls = 105974L)
GW$n_total <- GW$n_cases + GW$n_controls
GW$s <- GW$n_cases / GW$n_total

PUBLISHED <- c(EFNA1 = 0.963, ATRAID = 0.996, TNFRSF6B = 0.885)
TARGETS <- names(PUBLISHED)
CONFIGS <- list(
  list(name = "faithful",      palin = FALSE, boot = 0L, reg = TRUE),
  list(name = "no_palin",       palin = TRUE,  boot = 0L, reg = TRUE),
  list(name = "eurboot1",       palin = FALSE, boot = 1L, reg = TRUE),
  list(name = "eurboot2",       palin = FALSE, boot = 2L, reg = TRUE),
  list(name = "eurboot3",       palin = FALSE, boot = 3L, reg = TRUE),
  list(name = "no_regularise",  palin = FALSE, boot = 0L, reg = FALSE)
)

message("Loading 1000G EUR panel...")
panel <- fread("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt")
eur_ids <- panel[Superpopulation == "EUR", SampleID]

vcf_info <- function(chr) {
  url <- sprintf("%s/1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                 BASE_1KG, chr)
  hdr <- headerTabix(TabixFile(url))
  samples <- strsplit(hdr$header[length(hdr$header)], "\t", fixed = TRUE)[[1]][-(1:9)]
  list(url = url, eur_col = which(samples %in% eur_ids) + 9L)
}

## stream once; keep the biallelic-SNV filter from the original at parse time
fetch_dosages <- function(chr, lo, hi, eur_col, url) {
  message(sprintf("  streaming chr%s:%d-%d (once)...", chr, lo, hi))
  lines <- tryCatch(scanTabix(TabixFile(url),
                    param = GRanges(paste0("chr", chr), IRanges(lo, hi)))[[1]],
                    error = function(e) character(0))
  if (!length(lines)) return(NULL)
  ids <- character(0); cols <- list(); nskip <- 0L
  for (ln in lines) {
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(f) < 10) next
    ref <- toupper(f[4]); alt <- toupper(f[5])
    ## original keep_snp: biallelic single-nucleotide only
    if (nchar(ref) != 1L || nchar(alt) != 1L || grepl(",", alt, fixed = TRUE)) {
      nskip <- nskip + 1L; next
    }
    gt <- f[eur_col]
    d <- as.integer(substr(gt, 1L, 1L)) + as.integer(substr(gt, 3L, 3L))
    ids <- c(ids, paste0(f[2], ":", ref, ":", alt)); cols[[length(cols) + 1]] <- as.numeric(d)
  }
  if (!length(ids)) return(NULL)
  M <- do.call(cbind, cols); colnames(M) <- ids
  M <- M[, !duplicated(colnames(M)), drop = FALSE]
  message(sprintf("   %d biallelic SNVs cached (%d non-SNV/multiallelic skipped), %d samples",
                  ncol(M), nskip, nrow(M)))
  M
}

build_ld <- function(M, want, boot_seed, regularise) {
  S <- M
  if (boot_seed > 0L) {
    set.seed(1000L + boot_seed)
    S <- S[sample(seq_len(nrow(S)), floor(0.8 * nrow(S))), , drop = FALSE]
  }
  ## original keep_maf: reference-panel allele frequency filter
  af <- colMeans(S, na.rm = TRUE) / 2
  S <- S[, af > MAF_FLOOR & af < (1 - MAF_FLOOR), drop = FALSE]
  shared <- intersect(colnames(S), want)
  if (length(shared) < MIN_SNPS) return(NULL)
  S <- S[, shared, drop = FALSE]
  ld <- cor(S); ld[is.na(ld)] <- 0; diag(ld) <- 1
  if (regularise) {
    ## original: floor eigenvalues, reconstruct, rescale to correlation
    eig <- eigen(ld, symmetric = TRUE)
    eig$values <- pmax(eig$values, 1e-4)
    ldr <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    dinv <- 1 / sqrt(diag(ldr))
    ldr <- diag(dinv) %*% ldr %*% diag(dinv)
    diag(ldr) <- 1
    rownames(ldr) <- colnames(ldr) <- shared
    ld <- ldr
  }
  list(ld = ld, ids = shared)
}

load_pqtl <- function(p) {
  x <- fread(file.path(proj, "data/pqtl/priority_regions", paste0(p, "_pqtl_regions.tsv.gz")))
  x[!is.na(beta) & !is.na(se) & se > 0 & !is.na(alt_freq) & alt_freq > 0 & alt_freq < 1]
}
gwas_all <- NULL
gwas_region <- function(chr, lo, hi) {
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
    mutate(ea_p = toupper(alt), oa_p = toupper(ref), ea_g2 = toupper(ea), oa_g2 = toupper(oa),
           ea_pf = flip[toupper(alt)], oa_pf = flip[toupper(ref)],
           match_d  = ea_p == ea_g2 & oa_p == oa_g2,
           match_s  = ea_p == oa_g2 & oa_p == ea_g2,
           match_fl = !is.na(ea_pf) & ea_pf == ea_g2 & oa_pf == oa_g2,
           match_fs = !is.na(ea_pf) & ea_pf == oa_g2 & oa_pf == ea_g2,
           palin    = ea_p == flip[oa_p])
  h <- if (drop_palin) filter(h, !palin) else filter(h, !palin | (alt_freq > 0.1 & alt_freq < 0.9))
  h |> filter(match_d | match_s | match_fl | match_fs) |>
    mutate(beta_g_h = if_else(match_d | match_fl, beta_g, -beta_g),
           eaf_g_h  = if_else(match_d | match_fl, eaf_g, 1 - eaf_g))
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

res <- list()
for (prot in TARGETS) {
  message(sprintf("\n############ %s (published %.3f) ############", prot, PUBLISHED[[prot]]))
  pq <- load_pqtl(prot)
  lead <- pq[which.min(p)]
  chr <- as.character(lead$chr[1]); lp <- as.integer(lead$pos[1])
  lo <- lp - WINDOW_BP; hi <- lp + WINDOW_BP
  gw <- gwas_region(chr, lo, hi); pqw <- pq[pos >= lo & pos <= hi]
  vi <- vcf_info(chr)
  M <- fetch_dosages(chr, lo, hi, vi$eur_col, vi$url)
  if (is.null(M)) next

  for (cf in CONFIGS) {
    h <- harmonise(pqw, gw, cf$palin)
    key <- paste0(h$pos, ":", toupper(h$ref), ":", toupper(h$alt))
    li <- build_ld(M, key, cf$boot, cf$reg)
    if (is.null(li)) { message(sprintf("  [%-14s] no LD", cf$name)); next }
    idx <- match(li$ids, key); ok <- !is.na(idx)
    hs <- h[idx[ok], ]; ld <- li$ld[ok, ok, drop = FALSE]; ids <- li$ids[ok]
    if (nrow(hs) < MIN_SNPS) { message(sprintf("  [%-14s] too few", cf$name)); next }

    D1 <- list(beta = hs$beta, varbeta = hs$se^2, snp = ids, type = "quant",
               N = 619L, MAF = pmin(hs$alt_freq, 1 - hs$alt_freq), LD = ld)
    D2 <- list(beta = hs$beta_g_h, varbeta = hs$se_g^2, snp = ids, type = "cc",
               N = GW$n_total, s = GW$s, MAF = pmin(hs$eaf_g_h, 1 - hs$eaf_g_h), LD = ld)
    pp <- tryCatch(coloc.abf(D1[setdiff(names(D1), "LD")], D2[setdiff(names(D2), "LD")])$summary,
                   error = function(e) c(PP.H3.abf = NA, PP.H4.abf = NA))
    s1 <- tryCatch(runsusie(D1, repeat_until_convergence = TRUE, maxit = 10000L), error = function(e) NULL)
    s2 <- tryCatch(runsusie(D2, repeat_until_convergence = TRUE, maxit = 10000L), error = function(e) NULL)
    n1 <- if (!is.null(s1)) { s1 <- clean_lbf(s1); length(s1$sets$cs %||% list()) } else 0L
    n2 <- if (!is.null(s2)) { s2 <- clean_lbf(s2); length(s2$sets$cs %||% list()) } else 0L
    pph4 <- NA_real_; snp <- NA_character_; pairs <- 0L
    if (n1 > 0 && n2 > 0) {
      csr <- tryCatch(coloc.susie(s1, s2), error = function(e) NULL)
      if (is.null(csr)) csr <- tryCatch({
        b1 <- s1$lbf_variable[s1$sets$cs_index, , drop = FALSE]
        b2 <- s2$lbf_variable[s2$sets$cs_index, , drop = FALSE]
        b1[is.nan(b1)] <- 0; b2[is.nan(b2)] <- 0
        coloc:::coloc.bf_bf(b1, b2) }, error = function(e) NULL)
      if (!is.null(csr)) { e <- best_pph4(csr$summary); pph4 <- e$p; snp <- e$snp; pairs <- e$n }
    }
    message(sprintf("  [%-14s] n=%4d  CS pQTL=%d GWAS=%d  SuSiE PPH4=%s   (published %.3f)",
                    cf$name, nrow(hs), n1, n2,
                    ifelse(is.na(pph4), "NA", sprintf("%.4f", pph4)), PUBLISHED[[prot]]))
    res[[paste(prot, cf$name)]] <- data.table(
      protein = prot, published_PPH4_susie = PUBLISHED[[prot]], config = cf$name,
      ld_regularised = cf$reg, drop_palindromic = cf$palin, eur_subsample = cf$boot,
      n_ld_matched = nrow(hs), n_cs_pqtl = n1, n_cs_gwas = n2, n_coloc_pairs = pairs,
      PPH4_susie = pph4, susie_best_snp = snp,
      PPH3_abf = as.numeric(pp["PP.H3.abf"]), PPH4_abf = as.numeric(pp["PP.H4.abf"]))
    fwrite(rbindlist(res, fill = TRUE), OUT)
  }
  rm(M); gc(verbose = FALSE)
}
out <- rbindlist(res, fill = TRUE)
fwrite(out, OUT)
message("\n============ FAITHFUL REPRODUCTION + STABILITY ============")
print(out[, .(protein, config, n_ld_matched, n_cs_pqtl, n_cs_gwas, PPH4_susie, published_PPH4_susie)])
