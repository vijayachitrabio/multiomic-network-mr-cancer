#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)

WINDOW_BP <- 500000
MIN_SNPS <- 50
MAX_LOCI_PER_PAIR <- as.integer(Sys.getenv("MAX_LOCI_PER_PAIR", "10"))

cancer_meta <- data.frame(
  outcome = c("Breast_GCST90018757", "Endometrial_GCST006464", "Ovarian_GCST90016665"),
  file = c(
    "data/cancer_gwas/Breast_GCST90018757.h.tsv.gz",
    "data/cancer_gwas/Endometrial_GCST006464.h.tsv.gz",
    "data/cancer_gwas/Ovarian_GCST90016665.h.tsv.gz"
  ),
  n_cases = c(122977, 12906, 25509),
  n_controls = c(105974, 108979, 40138),
  stringsAsFactors = FALSE
)
cancer_meta$n_total <- cancer_meta$n_cases + cancer_meta$n_controls
cancer_meta$s <- cancer_meta$n_cases / cancer_meta$n_total

queue <- fread("results/validation/metabolite_cancer_evidence_tiers.csv")
queue <- queue[fdr < 0.05]

instruments <- fread("data/mqtl/mqtl_gwas_instruments.csv")

read_cancer <- function(path) {
  x <- fread(path, select = c(
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid"
  ))
  setnames(x, c(
    "chr", "pos", "ea_cancer", "oa_cancer", "beta_cancer", "se_cancer",
    "eaf_cancer", "p_cancer", "SNP"
  ))
  x <- x[!is.na(SNP) & SNP != "" & !is.na(beta_cancer) & !is.na(se_cancer) & se_cancer > 0]
  x[, chr := as.integer(chr)]
  x
}

read_metabolite <- function(metabolite) {
  path <- file.path("data/mqtl/mqtl_full_gwas", paste0(metabolite, "_full_regenie.tsv.gz"))
  x <- fread(path, select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P"))
  setnames(x, c("chr", "pos", "SNP", "oa_met", "ea_met", "eaf_met", "N_met", "beta_met", "se_met", "log10p_met"))
  x <- x[!is.na(SNP) & SNP != "" & !is.na(beta_met) & !is.na(se_met) & se_met > 0]
  x[, chr := as.integer(chr)]
  x[, p_met := 10^(-log10p_met)]
  x
}

align_alleles <- function(x) {
  x[, strand := fifelse(ea_met == ea_cancer & oa_met == oa_cancer, "same",
    fifelse(ea_met == oa_cancer & oa_met == ea_cancer, "flip", NA_character_)
  )]
  x <- x[!is.na(strand)]
  x[strand == "flip", `:=`(
    beta_cancer = -beta_cancer,
    eaf_cancer = 1 - eaf_cancer
  )]
  x
}

results <- list()

for (outcome in unique(queue$outcome)) {
  meta <- cancer_meta[cancer_meta$outcome == outcome, ]
  if (nrow(meta) != 1) stop("Missing cancer metadata for ", outcome)
  cat("Loading cancer GWAS:", outcome, "\n")
  cancer <- read_cancer(meta$file)

  for (metabolite in unique(queue$exposure[queue$outcome == outcome])) {
    cat("  Metabolite:", metabolite, "\n")
    met <- read_metabolite(metabolite)
    leads <- instruments[instruments$metabolite == metabolite]
    if (nrow(leads) == 0) next
    setorder(leads, pval)
    leads <- leads[!duplicated(paste(SNP, chr, pos))]
    if (!is.na(MAX_LOCI_PER_PAIR) && MAX_LOCI_PER_PAIR > 0 && nrow(leads) > MAX_LOCI_PER_PAIR) {
      leads <- leads[seq_len(MAX_LOCI_PER_PAIR)]
    }

    merged <- merge(met, cancer, by = "SNP", suffixes = c("_met", "_cancer"))
    if ("chr_met" %in% names(merged)) {
      merged <- merged[chr_met == chr_cancer]
      merged[, `:=`(chr = chr_met, pos = pos_met)]
    }
    merged <- align_alleles(merged)

    for (i in seq_len(nrow(leads))) {
      lead <- leads[i]
      locus <- merged[chr == lead$chr & pos >= lead$pos - WINDOW_BP & pos <= lead$pos + WINDOW_BP]
      locus <- locus[!duplicated(SNP)]

      base <- data.table(
        exposure = metabolite,
        outcome = outcome,
        lead_snp = lead$SNP,
        chr = lead$chr,
        lead_pos = lead$pos,
        n_snps = nrow(locus),
        status = "not_run",
        PP.H0 = NA_real_,
        PP.H1 = NA_real_,
        PP.H2 = NA_real_,
        PP.H3 = NA_real_,
        PP.H4 = NA_real_
      )

      if (nrow(locus) < MIN_SNPS) {
        base$status <- "too_few_overlapping_snps"
        results[[length(results) + 1]] <- base
        next
      }

      locus <- locus[!is.na(eaf_met) & eaf_met > 0 & eaf_met < 1]

      d1 <- list(
        beta = locus$beta_met,
        varbeta = locus$se_met^2,
        snp = locus$SNP,
        position = locus$pos,
        type = "quant",
        N = median(locus$N_met, na.rm = TRUE),
        MAF = pmin(locus$eaf_met, 1 - locus$eaf_met)
      )
      d2 <- list(
        beta = locus$beta_cancer,
        varbeta = locus$se_cancer^2,
        snp = locus$SNP,
        position = locus$pos,
        type = "cc",
        N = meta$n_total,
        s = meta$s,
        MAF = pmin(locus$eaf_cancer, 1 - locus$eaf_cancer)
      )

      fit <- tryCatch({
        invisible(capture.output(ans <- coloc.abf(d1, d2)))
        ans
      }, error = function(e) e)
      if (inherits(fit, "error")) {
        base$status <- paste0("error: ", conditionMessage(fit))
        results[[length(results) + 1]] <- base
      } else {
        pp <- fit$summary
        base$status <- "ok"
        base$PP.H0 <- unname(pp["PP.H0.abf"])
        base$PP.H1 <- unname(pp["PP.H1.abf"])
        base$PP.H2 <- unname(pp["PP.H2.abf"])
        base$PP.H3 <- unname(pp["PP.H3.abf"])
        base$PP.H4 <- unname(pp["PP.H4.abf"])
        results[[length(results) + 1]] <- base
      }
    }
  }
}

out <- rbindlist(results, fill = TRUE)
setorder(out, exposure, outcome, -PP.H4)
fwrite(out, "results/validation/metabolite_cancer_coloc_loci.csv")

best <- out[status == "ok", .SD[which.max(PP.H4)], by = .(exposure, outcome)]
setorder(best, -PP.H4)
fwrite(best, "results/validation/metabolite_cancer_coloc_best_loci.csv")

cat("Wrote metabolite-cancer coloc results:\n")
cat("  results/validation/metabolite_cancer_coloc_loci.csv\n")
cat("  results/validation/metabolite_cancer_coloc_best_loci.csv\n")
print(best[, .(exposure, outcome, lead_snp, n_snps, PP.H3, PP.H4)])
