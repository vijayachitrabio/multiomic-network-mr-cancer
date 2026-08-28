#!/usr/bin/env Rscript
## Script 83: ABF-only palindromic-SNP sensitivity for all 17 prioritized proteins
## ─────────────────────────────────────────────────────────────────────────────
## Companion to script 82. coloc.abf needs no LD reference, so this covers all
## 17 loci in minutes, whereas the SuSiE arm (script 82) must stream 1000G LD
## per locus and run runsusie() to convergence.
##
## Everything is identical to script 39 except the palindromic filter: this run
## drops ALL palindromic (A/T, C/G) variants instead of keeping those in the
## ambiguous 0.1-0.9 frequency band. Priors are coloc defaults (p1=p2=1e-4,
## p12=1e-5), MAF comes from the input datasets, and the case fraction s and
## total N are the published values per outcome.
##
## Output: results/validation/protein_coloc_abf_palindromic_sensitivity.csv

suppressPackageStartupMessages({
  library(data.table); library(coloc); library(dplyr); library(readr)
})
set.seed(2026)
proj <- "."
out_dir <- file.path(proj, "results/validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
WINDOW_BP <- 500000L

mk_gwas <- function(file, n_cases, n_controls) {
  g <- list(file = file.path(proj, "data/cancer_gwas", file),
            n_cases = n_cases, n_controls = n_controls)
  g$n_total <- g$n_cases + g$n_controls
  g$s <- g$n_cases / g$n_total
  g
}
GWAS <- list(
  Breast      = mk_gwas("Breast_GCST90018757.h.tsv.gz", 122977L, 105974L),
  Endometrial = mk_gwas("Endometrial_GCST006464.h.tsv.gz", 12906L, 108979L)
)
PROTEINS <- c("EFNA1","TNFRSF6B","ATRAID","FGF5","ABO","SNX15","PM20D1","UMOD",
              "APOE","TSPAN8","IL34","ITIH3","SWAP70","KLB","FGFR4","CGREF1","INHBB")
outcome_of <- function(p) if (p == "ABO") "Endometrial" else "Breast"

load_pqtl <- function(protein) {
  x <- fread(file.path(proj, "data/pqtl/priority_regions",
                       paste0(protein, "_pqtl_regions.tsv.gz")))
  x[!is.na(beta) & !is.na(se) & se > 0 & !is.na(alt_freq) & alt_freq > 0 & alt_freq < 1]
}

## Cache each outcome GWAS once -- reading a 300-400 MB gz per protein is the
## dominant cost otherwise.
gwas_cache <- list()
get_gwas <- function(key) {
  if (!is.null(gwas_cache[[key]])) return(gwas_cache[[key]])
  message(sprintf("  reading %s GWAS (once)...", key))
  g <- read_tsv(GWAS[[key]]$file, show_col_types = FALSE,
        col_types = cols_only(chromosome = col_integer(),
          base_pair_location = col_integer(), effect_allele = col_character(),
          other_allele = col_character(), beta = col_double(),
          standard_error = col_double(), effect_allele_frequency = col_double(),
          p_value = col_double())) |>
    rename(pos = base_pair_location, ea = effect_allele, oa = other_allele,
           beta_g = beta, se_g = standard_error,
           eaf_g = effect_allele_frequency, p_g = p_value) |>
    filter(!is.na(beta_g), !is.na(se_g), se_g > 0, !is.na(eaf_g), eaf_g > 0, eaf_g < 1)
  g <- as.data.table(g)
  gwas_cache[[key]] <<- g
  g
}

harmonise <- function(pqtl, gwas, drop_palindromic) {
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
  qc <- list(n_overlap = nrow(h),
             n_palin = sum(h$palin, na.rm = TRUE),
             n_palin_ambig = sum(h$palin & h$alt_freq > 0.1 & h$alt_freq < 0.9, na.rm = TRUE))
  h <- if (drop_palindromic) filter(h, !palin)
       else filter(h, !palin | (alt_freq > 0.1 & alt_freq < 0.9))
  h <- h |> filter(match_d | match_s | match_fl | match_fs) |>
    mutate(beta_g_h = if_else(match_d | match_fl, beta_g, -beta_g),
           eaf_g_h  = if_else(match_d | match_fl, eaf_g, 1 - eaf_g))
  qc$n_after <- nrow(h)
  attr(h, "qc") <- qc
  h
}

run_abf <- function(h, gcfg) {
  D1 <- list(beta = h$beta, varbeta = h$se^2, snp = paste0(h$pos, ":", h$ref, ":", h$alt),
             type = "quant", N = 619L, MAF = pmin(h$alt_freq, 1 - h$alt_freq))
  D2 <- list(beta = h$beta_g_h, varbeta = h$se_g^2, snp = D1$snp,
             type = "cc", N = gcfg$n_total, s = gcfg$s,
             MAF = pmin(h$eaf_g_h, 1 - h$eaf_g_h))
  tryCatch(coloc.abf(D1, D2)$summary, error = function(e) {
    message("  ABF error: ", e$message)
    c(PP.H0.abf = NA, PP.H1.abf = NA, PP.H2.abf = NA, PP.H3.abf = NA, PP.H4.abf = NA)
  })
}

res <- list()
for (prot in PROTEINS) {
  message(sprintf("\n===== %s (%s) =====", prot, outcome_of(prot)))
  key <- outcome_of(prot); gcfg <- GWAS[[key]]
  pq <- tryCatch(load_pqtl(prot), error = function(e) { message("  ", e$message); NULL })
  if (is.null(pq) || !nrow(pq)) next
  lead <- pq[which.min(p)]
  chr <- as.integer(lead$chr[1]); lead_pos <- as.integer(lead$pos[1])
  lo <- lead_pos - WINDOW_BP; hi <- lead_pos + WINDOW_BP

  gw_all <- get_gwas(key)
  gw <- gw_all[chromosome == chr & pos >= lo & pos <= hi]
  pq_w <- pq[pos >= lo & pos <= hi]
  if (!nrow(gw) || !nrow(pq_w)) { message("  no regional overlap"); next }

  ## published behaviour (ambiguous palindromes retained) vs conservative (all dropped)
  h_pub <- harmonise(pq_w, gw, drop_palindromic = FALSE)
  h_new <- harmonise(pq_w, gw, drop_palindromic = TRUE)
  qc <- attr(h_new, "qc")
  message(sprintf("  overlap %d | palindromic %d (%d in ambiguous band) | kept %d vs %d",
                  qc$n_overlap, qc$n_palin, qc$n_palin_ambig, nrow(h_pub), nrow(h_new)))
  if (nrow(h_new) < 10) { message("  too few variants after filtering"); next }

  pp_pub <- run_abf(h_pub, gcfg)
  pp_new <- run_abf(h_new, gcfg)
  message(sprintf("  as-published : PPH3=%.4f PPH4=%.4f", pp_pub["PP.H3.abf"], pp_pub["PP.H4.abf"]))
  message(sprintf("  no-palindrome: PPH3=%.4f PPH4=%.4f", pp_new["PP.H3.abf"], pp_new["PP.H4.abf"]))

  res[[prot]] <- data.table(
    protein = prot, outcome = key, lead_chr = chr, lead_pos = lead_pos,
    n_overlap = qc$n_overlap, n_palindromic = qc$n_palin,
    n_palindromic_ambiguous_band = qc$n_palin_ambig,
    n_kept_as_published = nrow(h_pub), n_kept_no_palindromic = nrow(h_new),
    PPH3_abf_as_published = as.numeric(pp_pub["PP.H3.abf"]),
    PPH4_abf_as_published = as.numeric(pp_pub["PP.H4.abf"]),
    PPH3_abf_no_palindromic = as.numeric(pp_new["PP.H3.abf"]),
    PPH4_abf_no_palindromic = as.numeric(pp_new["PP.H4.abf"]),
    delta_PPH4 = as.numeric(pp_new["PP.H4.abf"]) - as.numeric(pp_pub["PP.H4.abf"]))
  fwrite(rbindlist(res, fill = TRUE),
         file.path(out_dir, "protein_coloc_abf_palindromic_sensitivity.csv"))
}

out <- rbindlist(res, fill = TRUE)
fwrite(out, file.path(out_dir, "protein_coloc_abf_palindromic_sensitivity.csv"))
message(sprintf("\nDone: %d proteins", nrow(out)))
print(out[, .(protein, n_palindromic, n_palindromic_ambiguous_band,
              PPH4_abf_as_published, PPH4_abf_no_palindromic, delta_PPH4)])
