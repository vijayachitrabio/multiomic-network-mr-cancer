#!/usr/bin/env Rscript

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

proj <- normalizePath(".")
out_dir <- file.path(proj, "results", "validation")
cache_dir <- file.path(proj, "data", "pqtl", "priority_regions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

target_table <- data.frame(
  protein = c("ABO", "CGREF1", "FGF5", "FGFR4", "INHBB", "KLB", "SWAP70", "TSPAN8", "UMOD"),
  cancer = c("Endometrial", rep("Breast", 8)),
  chr = c("9", "2", "4", "5", "2", "4", "11", "12", "16"),
  lead_pos = c(133261703L, 27101856L, 80263187L, 177093242L, 120548864L,
               39448609L, 9753544L, 71143549L, 20381010L),
  stringsAsFactors = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
protein <- if (length(args) >= 1) args[1] else "FGF5"
abf_only <- length(args) >= 2 && args[2] == "abf_only"
if (!protein %in% target_table$protein) {
  stop("Unknown protein: ", protein, ". Valid targets: ", paste(target_table$protein, collapse = ", "))
}

target <- target_table[target_table$protein == protein, ]
chr <- target$chr
lead_pos <- target$lead_pos
cancer_label <- target$cancer
window_bp <- 500000L
min_snps <- 50L
maf_floor <- 0.01

if (cancer_label == "Breast") {
  gwas_cfg <- list(
    file = file.path(proj, "data", "cancer_gwas", "Breast_GCST90018757.h.tsv.gz"),
    n_cases = 122977L,
    n_controls = 105974L
  )
} else {
  gwas_cfg <- list(
    file = file.path(proj, "data", "cancer_gwas", "Endometrial_GCST006464.h.tsv.gz"),
    n_cases = 12906L,
    n_controls = 108979L
  )
}
gwas_cfg$n_total <- gwas_cfg$n_cases + gwas_cfg$n_controls
gwas_cfg$s <- gwas_cfg$n_cases / gwas_cfg$n_total

base_pqtl <- "https://storage.googleapis.com/finngen-public-data-r10/omics/proteomics/release_2023_03_02/data/Olink/pQTL"
base_1kg <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
panel_url <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

cache_file <- file.path(cache_dir, paste0(protein, "_pqtl_regions.tsv.gz"))
out_file <- file.path(out_dir, paste0("protein_coloc_pending_pilot_", protein, ".csv"))

region_start <- lead_pos - window_bp
region_end <- lead_pos + window_bp
region_label <- sprintf("%s:%d-%d", chr, region_start, region_end)

cache_pqtl_region <- function() {
  if (file.exists(cache_file) && file.info(cache_file)$size > 0) {
    message("Using cached pQTL region: ", cache_file)
    return(invisible(cache_file))
  }

  url <- sprintf("%s/Olink_Batch1_%s.txt.gz", base_pqtl, protein)
  message("Fetching pQTL region from FinnGen Olink: ", url)
  region <- GRanges(chr, IRanges(region_start, region_end))
  lines <- scanTabix(TabixFile(url), param = region)[[1]]
  if (length(lines) == 0) stop("No pQTL lines returned for ", protein, " ", region_label)

  x <- fread(text = paste(lines, collapse = "\n"), header = FALSE, showProgress = FALSE)
  if (ncol(x) == 11) {
    setnames(x, c("chr", "pos", "variant_id", "ref", "alt", "alt_freq",
                  "beta", "se", "t_stat", "p", "log10_p"))
    x[, n := 619L]
  } else {
    setnames(x, c("chr", "pos", "variant_id", "ref", "alt", "alt_freq",
                  "beta", "se", "t_stat", "p", "log10_p", "n"))
  }
  x[, `:=`(protein = protein, region = region_label)]
  fwrite(x, cache_file)
  message("Cached ", nrow(x), " pQTL variants: ", cache_file)
  invisible(cache_file)
}

load_pqtl <- function() {
  cache_pqtl_region()
  x <- fread(cache_file)
  x[!is.na(beta) & !is.na(se) & se > 0 & !is.na(alt_freq) & alt_freq > 0 & alt_freq < 1]
}

load_gwas <- function(pos_min, pos_max) {
  read_tsv(gwas_cfg$file, show_col_types = FALSE,
           col_types = cols_only(
             chromosome = col_integer(),
             base_pair_location = col_integer(),
             effect_allele = col_character(),
             other_allele = col_character(),
             beta = col_double(),
             standard_error = col_double(),
             effect_allele_frequency = col_double(),
             p_value = col_double()
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

harmonise <- function(pqtl, gwas) {
  flip <- c(A = "T", T = "A", C = "G", G = "C")
  inner_join(
    pqtl |> mutate(pos = as.integer(pos)),
    gwas |> mutate(pos = as.integer(pos)),
    by = "pos"
  ) |>
    mutate(
      ea_p = toupper(alt), oa_p = toupper(ref),
      ea_g2 = toupper(ea), oa_g2 = toupper(oa),
      ea_pf = flip[ea_p], oa_pf = flip[oa_p],
      match_d = ea_p == ea_g2 & oa_p == oa_g2,
      match_s = ea_p == oa_g2 & oa_p == ea_g2,
      match_fl = !is.na(ea_pf) & ea_pf == ea_g2 & oa_pf == oa_g2,
      match_fs = !is.na(ea_pf) & ea_pf == oa_g2 & oa_pf == ea_g2,
      palin = ea_p == flip[oa_p]
    ) |>
    filter(!palin | (alt_freq > 0.1 & alt_freq < 0.9)) |>
    filter(match_d | match_s | match_fl | match_fs) |>
    mutate(
      beta_g_h = if_else(match_d | match_fl, beta_g, -beta_g),
      eaf_g_h = if_else(match_d | match_fl, eaf_g, 1 - eaf_g)
    )
}

get_eur_vcf <- function() {
  panel <- fread(panel_url)
  eur_ids <- panel[Superpopulation == "EUR", SampleID]
  vcf_url <- sprintf("%s/1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                     base_1kg, chr)
  hdr <- headerTabix(TabixFile(vcf_url))
  chrom_line <- hdr$header[length(hdr$header)]
  all_samples <- strsplit(chrom_line, "\t", fixed = TRUE)[[1]][-(1:9)]
  list(url = vcf_url, eur_col = which(all_samples %in% eur_ids) + 9L, n_eur = length(eur_ids))
}

build_ld_matrix <- function(pos_min, pos_max, eur_col, vcf_url, target_snp_ids) {
  region <- GRanges(paste0("chr", chr), IRanges(pos_min, pos_max))
  lines <- scanTabix(TabixFile(vcf_url), param = region)[[1]]
  if (length(lines) == 0) return(NULL)

  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  n_lines <- length(split_lines)
  pos_vec <- as.integer(vapply(split_lines, `[[`, character(1), 2))
  ref_vec <- vapply(split_lines, `[[`, character(1), 4)
  alt_vec <- vapply(split_lines, `[[`, character(1), 5)
  snp_key <- paste0(pos_vec, ":", ref_vec, ":", alt_vec)
  keep_snp <- nchar(ref_vec) == 1 & nchar(alt_vec) == 1 & !grepl(",", alt_vec, fixed = TRUE)

  geno <- matrix(NA_real_, nrow = n_lines, ncol = length(eur_col))
  for (i in seq_len(n_lines)) {
    flds <- split_lines[[i]][eur_col]
    geno[i, ] <- as.integer(substr(flds, 1L, 1L)) + as.integer(substr(flds, 3L, 3L))
  }

  af <- rowMeans(geno, na.rm = TRUE) / 2
  keep <- keep_snp & af > maf_floor & af < (1 - maf_floor)
  snp_key_f <- snp_key[keep]
  geno_f <- geno[keep, , drop = FALSE]
  rownames(geno_f) <- snp_key_f

  shared <- intersect(snp_key_f, target_snp_ids)
  message("LD shared variants: ", length(shared))
  if (length(shared) < min_snps) return(NULL)

  geno_s <- geno_f[shared, , drop = FALSE]
  ld <- cor(t(geno_s))
  ld[is.na(ld)] <- 0
  diag(ld) <- 1
  eig <- eigen(ld, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-4)
  ld_reg <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  d_inv <- 1 / sqrt(diag(ld_reg))
  ld_reg <- diag(d_inv) %*% ld_reg %*% diag(d_inv)
  diag(ld_reg) <- 1
  rownames(ld_reg) <- colnames(ld_reg) <- shared
  list(ld = ld_reg, snp_ids = shared)
}

clean_lbf <- function(s) {
  if (!is.null(s$lbf_variable)) {
    bad <- is.na(s$lbf_variable) | is.nan(s$lbf_variable)
    if (any(bad)) s$lbf_variable[bad] <- 0
  }
  s
}

extract_pph4 <- function(smry) {
  if (is.null(smry)) return(list(pph4 = 0, snp = NA_character_, n = 0L))
  if (is.numeric(smry) && "PP.H4.abf" %in% names(smry)) {
    return(list(pph4 = as.numeric(smry["PP.H4.abf"]), snp = NA_character_, n = 1L))
  }
  if ((is.data.frame(smry) || is.data.table(smry)) && nrow(smry) > 0 && "PP.H4.abf" %in% names(smry)) {
    best <- which.max(as.numeric(smry[["PP.H4.abf"]]))
    return(list(
      pph4 = as.numeric(smry[["PP.H4.abf"]][best]),
      snp = if ("hit1" %in% names(smry)) as.character(smry[["hit1"]][best]) else NA_character_,
      n = nrow(smry)
    ))
  }
  list(pph4 = 0, snp = NA_character_, n = 0L)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

run_coloc <- function(harm, ld_info) {
  snp_key_ref <- paste0(harm$pos, ":", toupper(harm$ref), ":", toupper(harm$alt))
  shared <- ld_info$snp_ids
  idx <- match(shared, snp_key_ref)
  ok <- !is.na(idx)
  n_ok <- sum(ok)

  harm_sub <- harm[idx[ok], ]
  ld_sub <- ld_info$ld[ok, ok, drop = FALSE]
  snp_sub <- shared[ok]

  make_d <- function(use_ld = TRUE) {
    d1 <- list(beta = harm_sub$beta, varbeta = harm_sub$se^2, snp = snp_sub,
               type = "quant", N = 619L,
               MAF = pmin(harm_sub$alt_freq, 1 - harm_sub$alt_freq))
    d2 <- list(beta = harm_sub$beta_g_h, varbeta = harm_sub$se_g^2, snp = snp_sub,
               type = "cc", N = gwas_cfg$n_total, s = gwas_cfg$s,
               MAF = pmin(harm_sub$eaf_g_h, 1 - harm_sub$eaf_g_h))
    if (use_ld) {
      d1$LD <- ld_sub
      d2$LD <- ld_sub
    }
    list(D1 = d1, D2 = d2)
  }

  ds <- make_d(use_ld = FALSE)
  pp_abf <- tryCatch(
    coloc.abf(ds$D1, ds$D2)$summary,
    error = function(e) {
      message("ABF error: ", e$message)
      c(PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
        PP.H3.abf = NA_real_, PP.H4.abf = NA_real_)
    }
  )

  pph4_s <- 0
  best_snp <- NA_character_
  n_pairs <- 0L
  n_cs1 <- 0L
  n_cs2 <- 0L

  if (n_ok >= min_snps && !abf_only) {
    ds_ld <- make_d(use_ld = TRUE)
    s1 <- tryCatch(runsusie(ds_ld$D1, repeat_until_convergence = TRUE, maxit = 10000L),
                   error = function(e) { message("pQTL SuSiE error: ", e$message); NULL })
    s2 <- tryCatch(runsusie(ds_ld$D2, repeat_until_convergence = TRUE, maxit = 10000L),
                   error = function(e) { message("GWAS SuSiE error: ", e$message); NULL })
    if (!is.null(s1)) {
      s1 <- clean_lbf(s1)
      n_cs1 <- length(s1$sets$cs %||% list())
    }
    if (!is.null(s2)) {
      s2 <- clean_lbf(s2)
      n_cs2 <- length(s2$sets$cs %||% list())
    }
    if (!is.null(s1) && !is.null(s2) && n_cs1 > 0 && n_cs2 > 0) {
      csr <- tryCatch(coloc.susie(s1, s2), error = function(e) {
        message("coloc.susie error: ", e$message)
        NULL
      })
      if (!is.null(csr)) {
        ex <- extract_pph4(csr$summary)
        pph4_s <- ex$pph4
        best_snp <- ex$snp
        n_pairs <- ex$n
      }
    }
  }

  pp4 <- as.numeric(pp_abf["PP.H4.abf"])
  pp3 <- as.numeric(pp_abf["PP.H3.abf"])
  interpretation <- dplyr::case_when(
    is.finite(pph4_s) & pph4_s >= 0.8 ~ "STRONG coloc (SuSiE)",
    is.finite(pph4_s) & pph4_s >= 0.5 ~ "MODERATE coloc (SuSiE)",
    !is.na(pp4) & pp4 >= 0.8 ~ "STRONG coloc (ABF only)",
    !is.na(pp4) & pp4 >= 0.5 ~ "MODERATE coloc (ABF only)",
    !is.na(pp3) & pp3 >= 0.5 ~ "DISTINCT causal variants",
    TRUE ~ "INSUFFICIENT evidence"
  )

  data.frame(
    protein = protein,
    cancer = cancer_label,
    n_harm = nrow(harm),
    n_ld_harm = as.integer(n_ok),
    n_cs_pqtl = as.integer(n_cs1),
    n_cs_gwas = as.integer(n_cs2),
    n_coloc_pairs = as.integer(n_pairs),
    PPH4_susie = round(as.numeric(pph4_s), 4),
    susie_best_snp = best_snp,
    PPH0_abf = round(as.numeric(pp_abf["PP.H0.abf"]), 4),
    PPH1_abf = round(as.numeric(pp_abf["PP.H1.abf"]), 4),
    PPH2_abf = round(as.numeric(pp_abf["PP.H2.abf"]), 4),
    PPH3_abf = round(as.numeric(pp_abf["PP.H3.abf"]), 4),
    PPH4_abf = round(as.numeric(pp_abf["PP.H4.abf"]), 4),
    interpretation = interpretation,
    stringsAsFactors = FALSE
  )
}

message("Starting ", protein, " -> ", cancer_label, " coloc pilot.")
pqtl <- load_pqtl()
message("pQTL variants: ", nrow(pqtl))
gwas <- load_gwas(min(pqtl$pos), max(pqtl$pos))
message("GWAS variants: ", nrow(gwas))
harm <- harmonise(pqtl, gwas)
message("Harmonised variants: ", nrow(harm))
if (nrow(harm) < min_snps) stop("Too few harmonised variants for coloc pilot.")

vcf <- get_eur_vcf()
target_ids <- union(
  paste0(harm$pos, ":", toupper(harm$ref), ":", toupper(harm$alt)),
  paste0(harm$pos, ":", toupper(harm$alt), ":", toupper(harm$ref))
)
ld_info <- build_ld_matrix(min(pqtl$pos), max(pqtl$pos), vcf$eur_col, vcf$url, target_ids)
if (is.null(ld_info)) stop("LD matrix unavailable or too few shared variants.")

res <- run_coloc(harm, ld_info)
fwrite(res, out_file)
print(res)
message("Saved: ", out_file)
