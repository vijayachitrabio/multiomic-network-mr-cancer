#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(rtracklayer)
  library(GenomicRanges)
})

dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)

MIN_SNPS <- as.integer(Sys.getenv("COLOC_MIN_SNPS", "50"))

cancer_meta <- data.table(
  cancer = "Breast_GCST90018757",
  file = "data/cancer_gwas/Breast_GCST90018757.h.tsv.gz",
  n_cases = 122977,
  n_controls = 105974
)
cancer_meta[, `:=`(
  n_total = n_cases + n_controls,
  s = n_cases / (n_cases + n_controls)
)]

top_paths <- fread("results/validation/top_mediation_paths_for_followup.csv")

make_variant_key <- function(chr, pos, allele1, allele2) {
  chr <- sub("^chr", "", as.character(chr))
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  a_min <- ifelse(allele1 <= allele2, allele1, allele2)
  a_max <- ifelse(allele1 <= allele2, allele2, allele1)
  paste(chr, pos, a_min, a_max, sep = ":")
}

chain_file <- Sys.getenv("CHAIN_FILE", unset = "/tmp/hg38ToHg19.over.chain")
if (!file.exists(chain_file)) {
  stop("Missing liftover chain file: ", chain_file)
}
chain <- import.chain(chain_file)

read_pqtl <- function(protein) {
  path <- file.path("data/pqtl/priority_regions", paste0(protein, "_pqtl_regions.tsv.gz"))
  if (!file.exists(path)) stop("Missing pQTL regional file for ", protein, ": ", path)
  x <- fread(path)
  x <- x[!is.na(beta) & !is.na(se) & se > 0 & !is.na(alt_freq) & alt_freq > 0 & alt_freq < 1]
  x[, `:=`(
    chr = as.integer(chr),
    ea_pqtl = alt,
    oa_pqtl = ref,
    beta_pqtl = beta,
    se_pqtl = se,
    eaf_pqtl = alt_freq,
    n_pqtl = fifelse(is.na(n), 619, as.numeric(n))
  )]
  x[, .(protein, chr, pos, variant_id, ea_pqtl, oa_pqtl, beta_pqtl, se_pqtl, eaf_pqtl, n_pqtl)]
}

lift_pqtl_to_gr37 <- function(pqtl) {
  gr38 <- GRanges(
    seqnames = paste0("chr", pqtl$chr),
    ranges = IRanges(start = pqtl$pos, width = 1)
  )
  names(gr38) <- seq_len(nrow(pqtl))
  gr37 <- unlist(liftOver(gr38, chain))
  if (length(gr37) == 0) return(data.table())

  lift <- data.table(
    row_id = as.integer(names(gr37)),
    chr37 = sub("^chr", "", as.character(seqnames(gr37))),
    pos37 = start(gr37)
  )
  out <- pqtl[lift$row_id]
  out[, `:=`(chr37 = as.integer(lift$chr37), pos37 = lift$pos37)]
  out[, SNP_gr37 := make_variant_key(chr37, pos37, ea_pqtl, oa_pqtl)]
  out
}

read_cancer <- function(path) {
  x <- fread(path, select = c(
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid"
  ))
  setnames(x, c(
    "chr", "pos", "ea_out", "oa_out", "beta_out", "se_out",
    "eaf_out", "p_out", "rsid"
  ))
  x <- x[!is.na(beta_out) & !is.na(se_out) & se_out > 0 & !is.na(eaf_out) & eaf_out > 0 & eaf_out < 1]
  x[, chr := as.integer(chr)]
  x
}

read_metabolite <- function(metabolite) {
  path <- file.path("data/mqtl/mqtl_full_gwas", paste0(metabolite, "_full_regenie.tsv.gz"))
  x <- fread(path, select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P"))
  setnames(x, c("chr", "pos", "rsid", "oa_out", "ea_out", "eaf_out", "n_out", "beta_out", "se_out", "log10p_out"))
  x <- x[!is.na(beta_out) & !is.na(se_out) & se_out > 0 & !is.na(eaf_out) & eaf_out > 0 & eaf_out < 1]
  x[, chr := as.integer(chr)]
  x[, SNP_gr37 := make_variant_key(chr, pos, ea_out, oa_out)]
  x
}

align_to_pqtl <- function(x) {
  x[, strand := fifelse(ea_pqtl == ea_out & oa_pqtl == oa_out, "same",
    fifelse(ea_pqtl == oa_out & oa_pqtl == ea_out, "flip", NA_character_)
  )]
  x <- x[!is.na(strand)]
  x[strand == "flip", `:=`(
    beta_out = -beta_out,
    eaf_out = 1 - eaf_out
  )]
  x
}

run_coloc <- function(locus, outcome_type, outcome_n = NA_real_, outcome_s = NA_real_) {
  locus <- unique(locus, by = c("chr", "pos", "ea_pqtl", "oa_pqtl"))
  if (nrow(locus) < MIN_SNPS) {
    return(list(status = "too_few_overlapping_snps", summary = NULL, n = nrow(locus)))
  }

  d1 <- list(
    beta = locus$beta_pqtl,
    varbeta = locus$se_pqtl^2,
    snp = paste(locus$chr, locus$pos, locus$ea_pqtl, locus$oa_pqtl, sep = ":"),
    position = locus$pos,
    type = "quant",
    N = median(locus$n_pqtl, na.rm = TRUE),
    MAF = pmin(locus$eaf_pqtl, 1 - locus$eaf_pqtl)
  )

  d2 <- list(
    beta = locus$beta_out,
    varbeta = locus$se_out^2,
    snp = paste(locus$chr, locus$pos, locus$ea_pqtl, locus$oa_pqtl, sep = ":"),
    position = locus$pos,
    type = outcome_type,
    MAF = pmin(locus$eaf_out, 1 - locus$eaf_out)
  )
  if (outcome_type == "cc") {
    d2$N <- outcome_n
    d2$s <- outcome_s
  } else {
    d2$N <- median(locus$n_out, na.rm = TRUE)
  }

  fit <- tryCatch({
    invisible(capture.output(ans <- coloc.abf(d1, d2)))
    ans
  }, error = function(e) e)

  if (inherits(fit, "error")) {
    list(status = paste0("error: ", conditionMessage(fit)), summary = NULL, n = nrow(locus))
  } else {
    list(status = "ok", summary = fit$summary, n = nrow(locus))
  }
}

format_result <- function(result, leg, protein, mediator, cancer) {
  out <- data.table(
    leg = leg,
    protein = protein,
    metabolite = mediator,
    cancer = cancer,
    n_snps = result$n,
    status = result$status,
    PP.H0 = NA_real_,
    PP.H1 = NA_real_,
    PP.H2 = NA_real_,
    PP.H3 = NA_real_,
    PP.H4 = NA_real_
  )
  if (!is.null(result$summary)) {
    out[, `:=`(
      PP.H0 = unname(result$summary["PP.H0.abf"]),
      PP.H1 = unname(result$summary["PP.H1.abf"]),
      PP.H2 = unname(result$summary["PP.H2.abf"]),
      PP.H3 = unname(result$summary["PP.H3.abf"]),
      PP.H4 = unname(result$summary["PP.H4.abf"])
    )]
  }
  out
}

cancer <- read_cancer(cancer_meta$file)
results <- list()

for (i in seq_len(nrow(top_paths))) {
  protein <- top_paths$protein[i]
  metabolite <- top_paths$metabolite[i]
  cancer_name <- top_paths$cancer[i]
  cat(sprintf("Coloc: %s -> %s -> %s\n", protein, metabolite, cancer_name))

  pqtl <- read_pqtl(protein)

  pc <- merge(pqtl, cancer, by = c("chr", "pos"))
  pc <- align_to_pqtl(pc)
  pc_fit <- run_coloc(pc, "cc", cancer_meta$n_total, cancer_meta$s)
  results[[length(results) + 1]] <- format_result(pc_fit, "protein_cancer", protein, metabolite, cancer_name)

  met <- read_metabolite(metabolite)
  pqtl_gr37 <- lift_pqtl_to_gr37(pqtl)
  pm <- merge(pqtl_gr37, met, by = "SNP_gr37")
  if ("chr.x" %in% names(pm)) pm[, chr := chr.x]
  if ("pos.x" %in% names(pm)) pm[, pos := pos.x]
  pm <- align_to_pqtl(pm)
  pm_fit <- run_coloc(pm, "quant")
  results[[length(results) + 1]] <- format_result(pm_fit, "protein_metabolite", protein, metabolite, cancer_name)
}

out <- rbindlist(results, fill = TRUE)
setorder(out, protein, metabolite, leg)
fwrite(out, "results/validation/priority_mediation_leg_coloc.csv")

wide <- dcast(
  out,
  protein + metabolite + cancer ~ leg,
  value.var = c("status", "n_snps", "PP.H4", "PP.H3")
)

wide[, evidence_note := fifelse(
  !is.na(PP.H4_protein_cancer) & !is.na(PP.H4_protein_metabolite) &
    PP.H4_protein_cancer >= 0.8 & PP.H4_protein_metabolite >= 0.8,
  "protein_legs_coloc_strong",
  fifelse(
    (!is.na(PP.H4_protein_cancer) & PP.H4_protein_cancer >= 0.8) |
      (!is.na(PP.H4_protein_metabolite) & PP.H4_protein_metabolite >= 0.8),
    "one_protein_leg_coloc_strong",
    fifelse(
      status_protein_metabolite == "too_few_overlapping_snps",
      "protein_metabolite_coloc_blocked_by_low_overlap",
      "protein_legs_not_coloc_strong"
    )
  )
)]
fwrite(wide, "results/validation/priority_mediation_leg_coloc_wide.csv")

cat("Wrote priority mediation leg coloc:\n")
cat("  results/validation/priority_mediation_leg_coloc.csv\n")
cat("  results/validation/priority_mediation_leg_coloc_wide.csv\n\n")
print(wide)
