#!/usr/bin/env Rscript

# Script 29: Prepare MAGMA SNP-location and p-value inputs
#
# Builds MAGMA-ready input files from the harmonised cancer GWAS used in the MR
# pipeline. This complements the gprofiler2 over-representation analysis with
# gene-level GWAS aggregation from the full breast and endometrial summary stats.

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- normalizePath(".")
out_dir <- file.path(project_dir, "results", "pathway", "magma")
input_dir <- file.path(out_dir, "inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

cancers <- data.table(
  trait = c("breast", "endometrial"),
  file = c(
    file.path(project_dir, "data", "cancer_gwas", "Breast_GCST90018757.h.tsv.gz"),
    file.path(project_dir, "data", "cancer_gwas", "Endometrial_GCST006464.h.tsv.gz")
  ),
  n_cases = c(122977L, 12906L),
  n_controls = c(105974L, 108979L)
)
cancers[, n_total := n_cases + n_controls]

read_trait <- function(path) {
  hdr <- names(fread(path, nrows = 0))
  rsid_col <- if ("rsid" %in% hdr) "rsid" else if ("hm_rsid" %in% hdr) "hm_rsid" else NA_character_
  chr_col <- if ("chromosome" %in% hdr) "chromosome" else if ("hm_chrom" %in% hdr) "hm_chrom" else NA_character_
  bp_col <- if ("base_pair_location" %in% hdr) "base_pair_location" else if ("hm_pos" %in% hdr) "hm_pos" else NA_character_
  p_col <- if ("p_value" %in% hdr) "p_value" else NA_character_

  needed <- c(rsid_col, chr_col, bp_col, p_col)
  if (anyNA(needed)) {
    stop("Could not detect the required MAGMA columns in: ", basename(path))
  }

  dt <- fread(path, select = needed, showProgress = FALSE)
  setnames(dt, old = needed, new = c("SNP", "CHR", "BP", "P"))

  dt[, `:=`(
    CHR = suppressWarnings(as.integer(CHR)),
    BP = suppressWarnings(as.integer(BP)),
    P = suppressWarnings(as.numeric(P))
  )]

  dt <- dt[
    !is.na(SNP) & SNP != "" & SNP != "NA" &
      !is.na(CHR) & !is.na(BP) &
      !is.na(P) & P > 0 & P <= 1
  ]
  dt <- dt[CHR %in% 1:22 & !is.na(BP) & BP > 0 & !is.na(P) & P > 0 & P <= 1]
  dt[P < 1e-300, P := 1e-300]

  # MAGMA expects unique SNP rows. Keep the strongest signal when duplicated.
  setorder(dt, SNP, P, BP)
  dt <- dt[, .SD[1], by = SNP]
  dt[]
}

summary_rows <- list()

for (i in seq_len(nrow(cancers))) {
  trait <- cancers$trait[[i]]
  path <- cancers$file[[i]]
  n_total <- cancers$n_total[[i]]

  message("Preparing MAGMA inputs for ", trait, " ...")
  dt <- read_trait(path)

  snploc <- dt[, .(SNP, CHR, BP)]
  pval <- dt[, .(SNP, P, N = n_total)]

  snploc_path <- file.path(input_dir, paste0(trait, ".snploc"))
  pval_path <- file.path(input_dir, paste0(trait, ".pval"))

  fwrite(snploc, snploc_path, sep = "\t")
  fwrite(pval, pval_path, sep = "\t")

  summary_rows[[trait]] <- data.table(
    trait = trait,
    input_file = basename(path),
    n_cases = cancers$n_cases[[i]],
    n_controls = cancers$n_controls[[i]],
    n_total = n_total,
    n_snps = nrow(dt),
    snploc_file = basename(snploc_path),
    pval_file = basename(pval_path)
  )

  message("  Wrote ", nrow(dt), " SNPs for ", trait)
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_dt, file.path(out_dir, "magma_input_summary.csv"))

message("Done. MAGMA inputs are in: ", input_dir)
