#!/usr/bin/env Rscript
## Script 90: SuSiE LD sensitivity to a Finnish (FIN) LD reference
## ─────────────────────────────────────────────────────────────────────────────
## The pQTL data are from Finnish participants. The expanded 1000 Genomes panel
## contains 99 Finnish (FIN) samples, so LD can be rebuilt from a
## population-matched reference. Pipeline identical to script 89; only the LD
## sample set changes: (a) all FIN (99), (b) FIN founders, and (c) three random
## draws of 99 non-FIN EUR samples as a size-matched control, to separate an
## ancestry effect from small-panel instability.
## Loci: EFNA1, ATRAID, TNFRSF6B (strong coloc.susie support).
## Output: results/validation/susie_finnish_ld.csv

suppressPackageStartupMessages({
  library(data.table); library(coloc); library(susieR)
  library(Rsamtools); library(GenomicRanges); library(dplyr); library(readr)
})
set.seed(2026)
proj <- "."
out_dir <- file.path(proj, "results/validation")
OUT <- file.path(out_dir, "susie_finnish_ld.csv")

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

message("Loading 1000G sample manifest...")
panel <- fread("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt")
eur_all <- panel[Superpopulation == "EUR", SampleID]
eur_unrelated <- panel[Superpopulation == "EUR" & FatherID == "0" & MotherID == "0", SampleID]
fin_all <- panel[Population == "FIN", SampleID]
fin_unrelated <- panel[Population == "FIN" & FatherID == "0" & MotherID == "0", SampleID]
nonfin_eur <- setdiff(eur_all, fin_all)
rand_ctrl <- lapply(1:3, function(i) sample(nonfin_eur, length(fin_all)))
message(sprintf("  FIN total: %d | FIN founders: %d | non-FIN EUR pool: %d", length(fin_all), length(fin_unrelated), length(nonfin_eur)))
message(sprintf("  EUR total: %d | EUR unrelated founders: %d | excluded as related: %d",
                length(eur_all), length(eur_unrelated), length(eur_all) - length(eur_unrelated)))

vcf_info <- function(chr) {
  url <- sprintf("%s/1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                 BASE_1KG, chr)
  hdr <- headerTabix(TabixFile(url))
  samples <- strsplit(hdr$header[length(hdr$header)], "\t", fixed = TRUE)[[1]][-(1:9)]
  eur_col <- which(samples %in% eur_all) + 9L
  eur_sample_order <- samples[eur_col - 9L]  # sample IDs in the same order as eur_col / M's rows
  list(url = url, eur_col = eur_col, eur_sample_order = eur_sample_order)
}

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

## M's rows are in eur_sample_order; row_subset restricts to those sample IDs
build_ld <- function(M, want, row_ids, row_subset = NULL) {
  S <- M
  if (!is.null(row_subset)) {
    keep_rows <- row_ids %in% row_subset
    S <- S[keep_rows, , drop = FALSE]
  }
  af <- colMeans(S, na.rm = TRUE) / 2
  S <- S[, af > MAF_FLOOR & af < (1 - MAF_FLOOR), drop = FALSE]
  shared <- intersect(colnames(S), want)
  if (length(shared) < MIN_SNPS) return(NULL)
  S <- S[, shared, drop = FALSE]
  ld <- cor(S); ld[is.na(ld)] <- 0; diag(ld) <- 1
  eig <- eigen(ld, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-4)
  ldr <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  dinv <- 1 / sqrt(diag(ldr))
  ldr <- diag(dinv) %*% ldr %*% diag(dinv)
  diag(ldr) <- 1
  rownames(ldr) <- colnames(ldr) <- shared
  list(ld = ldr, ids = shared, n_samples = nrow(S))
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
harmonise <- function(pqtl, gwas) {
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
  h <- filter(h, !palin | (alt_freq > 0.1 & alt_freq < 0.9))
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

## Each SuSiE fit is capped at 10 min so non-convergent specifications cannot stall the run;
## a capped fit is recorded as timed_out = TRUE, not as "no credible set".
TIMED_OUT <- FALSE
fit_susie <- function(D) {
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE))
  tryCatch({
    setTimeLimit(elapsed = 600, transient = TRUE)
    runsusie(D, repeat_until_convergence = TRUE, maxit = 10000L)
  }, error = function(e) {
    if (grepl("time limit", conditionMessage(e), ignore.case = TRUE)) TIMED_OUT <<- TRUE
    NULL
  })
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
  h <- harmonise(pqw, gw)
  key <- paste0(h$pos, ":", toupper(h$ref), ":", toupper(h$alt))

  for (cfg in list(list(name = "fin_all_99", subset = fin_all),
                   list(name = "fin_founders", subset = fin_unrelated),
                   list(name = "nonfin_eur_rand1_99", subset = rand_ctrl[[1]]),
                   list(name = "nonfin_eur_rand2_99", subset = rand_ctrl[[2]]),
                   list(name = "nonfin_eur_rand3_99", subset = rand_ctrl[[3]]))) {
    li <- build_ld(M, key, vi$eur_sample_order, cfg$subset)
    if (is.null(li)) { message(sprintf("  [%-18s] no LD", cfg$name)); next }
    idx <- match(li$ids, key); ok <- !is.na(idx)
    hs <- h[idx[ok], ]; ld <- li$ld[ok, ok, drop = FALSE]; ids <- li$ids[ok]
    if (nrow(hs) < MIN_SNPS) { message(sprintf("  [%-18s] too few", cfg$name)); next }

    D1 <- list(beta = hs$beta, varbeta = hs$se^2, snp = ids, type = "quant",
               N = 619L, MAF = pmin(hs$alt_freq, 1 - hs$alt_freq), LD = ld)
    D2 <- list(beta = hs$beta_g_h, varbeta = hs$se_g^2, snp = ids, type = "cc",
               N = GW$n_total, s = GW$s, MAF = pmin(hs$eaf_g_h, 1 - hs$eaf_g_h), LD = ld)
    pp <- tryCatch(coloc.abf(D1[setdiff(names(D1), "LD")], D2[setdiff(names(D2), "LD")])$summary,
                   error = function(e) c(PP.H3.abf = NA, PP.H4.abf = NA))
    TIMED_OUT <<- FALSE
    s1 <- fit_susie(D1)
    s2 <- fit_susie(D2)
    timed_out <- TIMED_OUT
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
    message(sprintf("  [%-18s] n_samples=%d n_snp=%4d CS pQTL=%d GWAS=%d  SuSiE PPH4=%s  (published %.3f)",
                    cfg$name, li$n_samples, nrow(hs), n1, n2,
                    ifelse(is.na(pph4), "NA", sprintf("%.4f", pph4)), PUBLISHED[[prot]]))
    res[[paste(prot, cfg$name)]] <- data.table(
      protein = prot, published_PPH4_susie = PUBLISHED[[prot]], config = cfg$name,
      n_samples = li$n_samples, n_ld_matched = nrow(hs), n_cs_pqtl = n1, n_cs_gwas = n2,
      n_coloc_pairs = pairs, PPH4_susie = pph4, susie_best_snp = snp,
      PPH3_abf = as.numeric(pp["PP.H3.abf"]), PPH4_abf = as.numeric(pp["PP.H4.abf"]),
      timed_out = timed_out)
    fwrite(rbindlist(res, fill = TRUE), OUT)
  }
  rm(M); gc(verbose = FALSE)
}
out <- rbindlist(res, fill = TRUE)
fwrite(out, OUT)
message("\n============ FINNISH LD SENSITIVITY ============")
print(out[, .(protein, config, n_samples, n_ld_matched, n_cs_pqtl, n_cs_gwas, PPH4_susie, published_PPH4_susie)])
