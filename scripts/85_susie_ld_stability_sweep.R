#!/usr/bin/env Rscript
## Script 85: LD-stability sweep for the SuSiE-dependent colocalization loci
## ─────────────────────────────────────────────────────────────────────────────
## QUESTION
## EFNA1 (published SuSiE PPH4 = 0.963) and ATRAID (0.996) are the only two Tier 1
## assignments that rest on coloc.susie rather than coloc.abf -- ABF calls both
## "distinct causal variants" (PPH3 ~0.91 and ~0.997). Script 84 could not
## reproduce either value when the 1000 Genomes EUR LD matrix was rebuilt, while
## TNFRSF6B (0.885) reproduced to 0.01. The original LD matrix was never saved, so
## the discrepancy cannot be resolved by inspection.
##
## This script asks the decisive question instead: is the published value
## recoverable under ANY defensible LD construction, or was it specific to one
## particular matrix? For each locus it holds the summary statistics fixed and
## varies only how LD is estimated, then reports the resulting spread of PPH4.
##
## DESIGN
## Genotypes for each locus are streamed from 1000G ONCE and cached in memory;
## every LD variant is then derived from that cache. This is what makes a 21-run
## sweep tractable -- re-streaming per configuration would dominate the runtime.
##
## Configurations per locus:
##   baseline      published filter, all EUR samples, missingness < 0.05
##   maf01         additionally require MAF >= 0.01 in both datasets
##   maf05         additionally require MAF >= 0.05
##   strict_miss   missingness < 0.01 (drops variants with any appreciable no-call)
##   no_palin      all palindromic variants dropped (script 84's arm, for continuity)
##   eurboot1..3   LD from three independent 80% subsamples of EUR individuals
##
## TNFRSF6B is carried through every configuration as a positive control: if the
## sweep is working, its PPH4 should stay near 0.885 throughout.
##
## Output: results/validation/susie_ld_stability_sweep.csv

suppressPackageStartupMessages({
  library(data.table); library(coloc); library(susieR)
  library(Rsamtools); library(GenomicRanges); library(dplyr); library(readr)
})
set.seed(2026)
proj <- "."
out_dir <- file.path(proj, "results/validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(out_dir, "susie_ld_stability_sweep.csv")

WINDOW_BP <- 500000L
MIN_SNPS  <- 50L
BASE_1KG <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
`%||%` <- function(a, b) if (is.null(a)) b else a

GW <- list(file = file.path(proj, "data/cancer_gwas/Breast_GCST90018757.h.tsv.gz"),
           n_cases = 122977L, n_controls = 105974L)
GW$n_total <- GW$n_cases + GW$n_controls
GW$s <- GW$n_cases / GW$n_total

PUBLISHED <- c(EFNA1 = 0.963, ATRAID = 0.996, TNFRSF6B = 0.885)
TARGETS <- c("EFNA1", "ATRAID", "TNFRSF6B")

CONFIGS <- list(
  list(name = "baseline",    maf = 0,    miss = 0.05, palin = FALSE, boot = 0L),
  list(name = "maf01",       maf = 0.01, miss = 0.05, palin = FALSE, boot = 0L),
  list(name = "maf05",       maf = 0.05, miss = 0.05, palin = FALSE, boot = 0L),
  list(name = "strict_miss", maf = 0,    miss = 0.01, palin = FALSE, boot = 0L),
  list(name = "no_palin",    maf = 0,    miss = 0.05, palin = TRUE,  boot = 0L),
  list(name = "eurboot1",    maf = 0,    miss = 0.05, palin = FALSE, boot = 1L),
  list(name = "eurboot2",    maf = 0,    miss = 0.05, palin = FALSE, boot = 2L),
  list(name = "eurboot3",    maf = 0,    miss = 0.05, palin = FALSE, boot = 3L)
)

message("Loading 1000G EUR panel...")
panel <- fread("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt")
eur_ids <- panel[Superpopulation == "EUR", SampleID]
message(sprintf("  %d EUR samples", length(eur_ids)))

vcf_info <- function(chr) {
  url <- sprintf("%s/1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                 BASE_1KG, chr)
  hdr <- headerTabix(TabixFile(url))
  samples <- strsplit(hdr$header[length(hdr$header)], "\t", fixed = TRUE)[[1]][-(1:9)]
  list(url = url, eur_col = which(samples %in% eur_ids) + 9L)
}

## stream the region once, return dosage matrix (samples x variants)
fetch_dosages <- function(chr, lo, hi, eur_col, url) {
  message(sprintf("  streaming 1000G chr%s:%d-%d (once)...", chr, lo, hi))
  lines <- tryCatch(scanTabix(TabixFile(url),
                    param = GRanges(paste0("chr", chr), IRanges(lo, hi)))[[1]],
                    error = function(e) { message("   tabix: ", e$message); character(0) })
  if (!length(lines)) return(NULL)
  ids <- character(0); cols <- list()
  for (ln in lines) {
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(f) < 10) next
    id <- paste0(f[2], ":", toupper(f[4]), ":", toupper(f[5]))
    gt <- f[eur_col]
    d <- vapply(gt, function(g) {
      a <- substr(g, 1, 1); b <- substr(g, 3, 3)
      if (a == "." || b == ".") return(NA_real_)
      as.numeric(a) + as.numeric(b) }, numeric(1), USE.NAMES = FALSE)
    ids <- c(ids, id); cols[[length(cols) + 1]] <- d
  }
  if (!length(ids)) return(NULL)
  M <- do.call(cbind, cols); colnames(M) <- ids
  M <- M[, !duplicated(colnames(M)), drop = FALSE]
  message(sprintf("   cached %d samples x %d variants", nrow(M), ncol(M)))
  M
}

ld_from_cache <- function(M, want, miss_thresh, boot_seed) {
  keep <- intersect(colnames(M), want)
  if (length(keep) < 2) return(NULL)
  S <- M[, keep, drop = FALSE]
  if (boot_seed > 0L) {
    set.seed(1000L + boot_seed)
    idx <- sample(seq_len(nrow(S)), size = floor(0.8 * nrow(S)))
    S <- S[idx, , drop = FALSE]
  }
  okv <- apply(S, 2, function(v) {
    mean(is.na(v)) <= miss_thresh && length(unique(v[!is.na(v)])) >= 2 })
  S <- S[, okv, drop = FALSE]
  if (ncol(S) < 2) return(NULL)
  ld <- suppressWarnings(cor(S, use = "pairwise.complete.obs"))
  ld[is.na(ld)] <- 0; diag(ld) <- 1
  list(ld = ld, ids = colnames(S))
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

harmonise <- function(pqtl, gwas, drop_palin, maf_floor) {
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
  h <- if (drop_palin) filter(h, !palin) else filter(h, !palin | (alt_freq > 0.1 & alt_freq < 0.9))
  h <- h |> filter(match_d | match_s | match_fl | match_fs) |>
    mutate(beta_g_h = if_else(match_d | match_fl, beta_g, -beta_g),
           eaf_g_h  = if_else(match_d | match_fl, eaf_g, 1 - eaf_g),
           maf_p = pmin(alt_freq, 1 - alt_freq),
           maf_g = pmin(eaf_g_h, 1 - eaf_g_h))
  if (maf_floor > 0) h <- filter(h, maf_p >= maf_floor, maf_g >= maf_floor)
  h
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
  message(sprintf("\n################ %s (published SuSiE PPH4 = %.3f) ################",
                  prot, PUBLISHED[[prot]]))
  pq <- load_pqtl(prot)
  lead <- pq[which.min(p)]
  chr <- as.character(lead$chr[1]); lp <- as.integer(lead$pos[1])
  lo <- lp - WINDOW_BP; hi <- lp + WINDOW_BP
  gw <- gwas_region(chr, lo, hi)
  pqw <- pq[pos >= lo & pos <= hi]
  vi <- vcf_info(chr)
  M <- fetch_dosages(chr, lo, hi, vi$eur_col, vi$url)
  if (is.null(M)) { message("  no genotypes; skipping"); next }

  for (cf in CONFIGS) {
    h <- harmonise(pqw, gw, drop_palin = cf$palin, maf_floor = cf$maf)
    if (nrow(h) < MIN_SNPS) {
      message(sprintf("  [%-11s] only %d harmonised variants; skipping", cf$name, nrow(h))); next
    }
    key <- paste0(h$pos, ":", toupper(h$ref), ":", toupper(h$alt))
    li <- ld_from_cache(M, key, cf$miss, cf$boot)
    if (is.null(li)) { message(sprintf("  [%-11s] no LD", cf$name)); next }
    idx <- match(li$ids, key); ok <- !is.na(idx)
    hs <- h[idx[ok], ]; ld <- li$ld[ok, ok, drop = FALSE]; ids <- li$ids[ok]
    if (nrow(hs) < MIN_SNPS) { message(sprintf("  [%-11s] too few LD-matched", cf$name)); next }

    D1 <- list(beta = hs$beta, varbeta = hs$se^2, snp = ids, type = "quant",
               N = 619L, MAF = hs$maf_p, LD = ld)
    D2 <- list(beta = hs$beta_g_h, varbeta = hs$se_g^2, snp = ids, type = "cc",
               N = GW$n_total, s = GW$s, MAF = hs$maf_g, LD = ld)

    ## ABF is LD-free: report it alongside as an internal anchor
    pp <- tryCatch(coloc.abf(D1[setdiff(names(D1), "LD")],
                             D2[setdiff(names(D2), "LD")])$summary,
                   error = function(e) c(PP.H3.abf = NA, PP.H4.abf = NA))

    s1 <- tryCatch(runsusie(D1, repeat_until_convergence = TRUE, maxit = 10000L),
                   error = function(e) NULL)
    s2 <- tryCatch(runsusie(D2, repeat_until_convergence = TRUE, maxit = 10000L),
                   error = function(e) NULL)
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
    message(sprintf("  [%-11s] n=%4d  CS pQTL=%d GWAS=%d  SuSiE PPH4=%s  (ABF PPH4=%.4f)",
                    cf$name, nrow(hs), n1, n2,
                    ifelse(is.na(pph4), "NA", sprintf("%.4f", pph4)),
                    as.numeric(pp["PP.H4.abf"])))

    res[[paste(prot, cf$name)]] <- data.table(
      protein = prot, published_PPH4_susie = PUBLISHED[[prot]], config = cf$name,
      maf_floor = cf$maf, miss_thresh = cf$miss, drop_palindromic = cf$palin,
      eur_subsample = cf$boot, n_ld_matched = nrow(hs),
      n_cs_pqtl = n1, n_cs_gwas = n2, n_coloc_pairs = pairs,
      PPH4_susie = pph4, susie_best_snp = snp,
      PPH3_abf = as.numeric(pp["PP.H3.abf"]), PPH4_abf = as.numeric(pp["PP.H4.abf"]))
    fwrite(rbindlist(res, fill = TRUE), OUT)
  }
  rm(M); gc(verbose = FALSE)
}

out <- rbindlist(res, fill = TRUE)
fwrite(out, OUT)
message("\n================ LD-STABILITY SWEEP ================")
print(out[, .(protein, config, n_ld_matched, n_cs_pqtl, n_cs_gwas,
              PPH4_susie, PPH4_abf, published_PPH4_susie)])
message("\n---- per-locus SuSiE PPH4 spread ----")
print(out[, .(n_configs = .N,
              n_recovering_published = sum(!is.na(PPH4_susie) &
                    abs(PPH4_susie - published_PPH4_susie) < 0.10),
              min_PPH4 = min(PPH4_susie, na.rm = TRUE),
              max_PPH4 = max(PPH4_susie, na.rm = TRUE),
              n_NA = sum(is.na(PPH4_susie))), by = protein])
