#!/usr/bin/env Rscript

# Script 04d: Protein -> Metabolite harmonisation (FinnGen R10 Olink 701 proteins)
#
# FinnGen pQTL instruments are GRCh38.
# Metabolite GWAS files (Rahu/Alasoo REGENIE output) are GRCh37.
#
# Strategy:
#   1. Liftover the 1062 FinnGen pQTL positions from GRCh38 → GRCh37
#      using rtracklayer + UCSC hg38ToHg19.over.chain.gz
#   2. Build GRCh37-based variant keys (sorted allele pairs) for both
#      the lifted pQTL positions and the metabolite GWAS rows
#   3. Match by GRCh37 variant key; format and harmonise as usual
#   4. Fix: after format_data(), apply toupper() to restore uppercase SNP column
#
# Chain file expected at /tmp/hg38ToHg19.over.chain.gz (or set CHAIN_FILE env var).

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(rtracklayer)
  library(GenomicRanges)
})

project_dir  <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
out_dir      <- file.path(project_dir, "data", "harmonised")
mqtl_gwas_dir <- file.path(project_dir, "data", "mqtl", "mqtl_full_gwas")

make_variant_key <- function(chr, pos, allele1, allele2) {
  chr    <- sub("^chr", "", as.character(chr))
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  a_min  <- ifelse(allele1 <= allele2, allele1, allele2)
  a_max  <- ifelse(allele1 <= allele2, allele2, allele1)
  paste(chr, pos, a_min, a_max, sep = ":")
}

# ---------------------------------------------------------------
# 1. Load FinnGen pQTL instruments
# ---------------------------------------------------------------
pqtl_file <- file.path(project_dir, "data", "pqtl", "pqtl_instruments.csv")
cat("Loading pQTL instruments...\n")
pqtls   <- fread(pqtl_file)
exp_dat_orig <- pqtls
cat(sprintf("  %d instruments across %d proteins (GRCh38)\n",
            nrow(exp_dat_orig), length(unique(exp_dat_orig$exposure))))
cat("  SNP sample:", paste(head(exp_dat_orig$SNP, 3), collapse = ", "), "\n\n")

# ---------------------------------------------------------------
# 2. Liftover GRCh38 → GRCh37 using rtracklayer
# ---------------------------------------------------------------
# Chain file: try decompressed first, then .gz; download if missing.
chain_file_gz <- Sys.getenv("CHAIN_FILE", unset = "/tmp/hg38ToHg19.over.chain.gz")
chain_file    <- sub("\\.gz$", "", chain_file_gz)
if (!file.exists(chain_file)) {
  if (file.exists(chain_file_gz)) {
    cat("Decompressing chain file...\n")
    R.utils::gunzip(chain_file_gz, destname = chain_file, remove = FALSE)
  } else {
    stop("Chain file not found. Download with:\n",
         "  curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
         " -o /tmp/hg38ToHg19.over.chain.gz && gunzip /tmp/hg38ToHg19.over.chain.gz")
  }
}

cat("Loading chain file:", chain_file, "\n")
chain <- import.chain(chain_file)

# Build GRanges for FinnGen positions using pre-parsed chr.exposure/pos.exposure columns
gr38 <- GRanges(
  seqnames = paste0("chr", exp_dat_orig$chr.exposure),
  ranges   = IRanges(start = exp_dat_orig$pos.exposure, width = 1)
)
names(gr38) <- exp_dat_orig$SNP

cat(sprintf("Lifting over %d pQTL positions GRCh38 → GRCh37...\n", length(gr38)))
gr37_list <- liftOver(gr38, chain)
gr37_flat <- unlist(gr37_list)   # only positions that lifted successfully

# Build lift table: original GRCh38 key → GRCh37 chr + pos
n_lifted  <- length(gr37_flat)
n_failed  <- length(gr38) - n_lifted
cat(sprintf("  Lifted: %d / %d (failed: %d)\n\n", n_lifted, n_lifted + n_failed, n_failed))

lift_tbl <- data.table(
  SNP_gr38 = names(gr37_flat),
  gr37_chr = as.character(seqnames(gr37_flat)),
  gr37_pos = start(gr37_flat)
)
lift_tbl[, gr37_chr := sub("^chr", "", gr37_chr)]

# Merge lifted positions back to exp_dat; use GRCh37 alleles from exp_dat (unchanged by liftover)
exp_dat_orig[, SNP_gr38 := SNP]   # keep original key
lift_tbl <- merge(lift_tbl, exp_dat_orig[, .(SNP_gr38=SNP, effect_allele.exposure, other_allele.exposure)],
                  by = "SNP_gr38", all.x = TRUE)
# Build GRCh37 variant key (sorted alleles, same convention as metabolite GWAS key)
lift_tbl[, SNP_gr37 := make_variant_key(gr37_chr, gr37_pos,
                                         effect_allele.exposure, other_allele.exposure)]
lift_tbl <- unique(lift_tbl[!is.na(SNP_gr37)])

cat(sprintf("GRCh37 variant keys generated: %d\n", nrow(lift_tbl)))
cat("Sample GRCh37 keys:", paste(head(lift_tbl$SNP_gr37, 5), collapse = ", "), "\n\n")

# Augment exp_dat with GRCh37 variant keys for matching
# Keep original exp_dat columns; add SNP_gr37 for outcome matching
exp_dat <- merge(
  exp_dat_orig,
  lift_tbl[, .(SNP_gr38, SNP_gr37)],
  by = "SNP_gr38", all.x = FALSE   # drop SNPs that failed liftover
)
setDF(exp_dat)
cat(sprintf("exp_dat after liftover filter: %d rows\n\n", nrow(exp_dat)))

# ---------------------------------------------------------------
# 3. Protein -> Metabolite harmonisation
# ---------------------------------------------------------------
mqtl_files <- list.files(mqtl_gwas_dir, pattern = "_full_regenie\\.tsv\\.gz$", full.names = TRUE)
cat(sprintf("Found %d metabolite GWAS files.\n\n", length(mqtl_files)))

harm_log <- list()
n_done   <- 0L

for (mf in mqtl_files) {
  met_name <- gsub("_full_regenie\\.tsv\\.gz$", "", basename(mf))
  n_done   <- n_done + 1L
  cat(sprintf("[%d/%d] %s\n", n_done, length(mqtl_files), met_name))

  out_dat_raw <- tryCatch(
    fread(mf, select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P")),
    error = function(e) { cat("  READ ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(out_dat_raw) || nrow(out_dat_raw) == 0) {
    harm_log[[met_name]] <- data.table(
      exposure="All_Proteins", outcome=met_name,
      type="Protein->Metabolite", n_before=nrow(exp_dat), n_after=0, n_dropped=nrow(exp_dat))
    next
  }

  # Build GRCh37 variant key from metabolite GWAS (CHROM, GENPOS are GRCh37)
  # ALLELE0=ref, ALLELE1=alt in REGENIE; sort alphabetically for key
  out_dat_raw[, SNP_gr37 := make_variant_key(CHROM, GENPOS, ALLELE0, ALLELE1)]

  # Filter to rows matching any lifted pQTL GRCh37 key
  out_dat_raw <- out_dat_raw[SNP_gr37 %in% exp_dat$SNP_gr37]
  cat(sprintf("  Matched (GRCh37 key): %d variants\n", nrow(out_dat_raw)))

  if (nrow(out_dat_raw) == 0) {
    saveRDS(data.frame(), file.path(out_dir, paste0("harmonised_protein_", met_name, ".rds")))
    harm_log[[met_name]] <- data.table(
      exposure="All_Proteins", outcome=met_name,
      type="Protein->Metabolite", n_before=nrow(exp_dat), n_after=0, n_dropped=nrow(exp_dat))
    next
  }

  setDF(out_dat_raw)
  out_dat <- format_data(
    out_dat_raw,
    type             = "outcome",
    snp_col          = "SNP_gr37",
    beta_col         = "BETA",
    se_col           = "SE",
    eaf_col          = "A1FREQ",
    effect_allele_col = "ALLELE1",
    other_allele_col  = "ALLELE0",
    pval_col         = "LOG10P",
    log_pval         = TRUE,
    samplesize_col   = "N",
    chr_col          = "CHROM",
    pos_col          = "GENPOS"
  )
  # format_data() lowercases the SNP col; restore uppercase to match exp_dat$SNP_gr37
  out_dat$SNP       <- toupper(out_dat$SNP)
  out_dat$outcome   <- met_name
  out_dat$id.outcome <- met_name

  # Temporarily set exp_dat$SNP to GRCh37 keys for this harmonisation round
  exp_dat_tmp       <- exp_dat
  exp_dat_tmp$SNP   <- toupper(exp_dat_tmp$SNP_gr37)
  # Fix id.exposure (TwoSampleMR uses this for grouping)
  exp_dat_tmp$id.exposure <- exp_dat_tmp$id.exposure  # keep as-is (protein ID)

  cat(sprintf("  SNP overlap (GRCh37): %d\n", length(intersect(exp_dat_tmp$SNP, out_dat$SNP))))

  harm_dat <- tryCatch(
    harmonise_data(exposure_dat = exp_dat_tmp, outcome_dat = out_dat, action = 2),
    error = function(e) { cat("  harmonise_data ERROR:", conditionMessage(e), "\n"); data.frame() }
  )
  n_aft  <- nrow(harm_dat)
  n_keep <- if (n_aft > 0) sum(harm_dat$mr_keep) else 0L
  cat(sprintf("  Harmonised: %d rows (%d mr_keep)\n", n_aft, n_keep))
  # Note: SNP column contains GRCh37 keys; MR analysis uses beta/se not SNP IDs, so this is fine.

  rds_path <- file.path(out_dir, paste0("harmonised_protein_", met_name, ".rds"))
  saveRDS(harm_dat, rds_path)

  harm_log[[met_name]] <- data.table(
    exposure="All_Proteins", outcome=met_name,
    type="Protein->Metabolite", n_before=nrow(exp_dat), n_after=n_aft, n_dropped=nrow(exp_dat) - n_aft)

  rm(out_dat_raw, out_dat, harm_dat); gc(verbose = FALSE)
}

# ---------------------------------------------------------------
# Update harmonisation log
# ---------------------------------------------------------------
log_file <- file.path(out_dir, "harmonisation_log.csv")
new_log  <- rbindlist(harm_log, fill = TRUE)
if (file.exists(log_file)) {
  old_log <- fread(log_file)
  old_log <- old_log[type != "Protein->Metabolite"]
  full_log <- rbind(old_log, new_log, fill = TRUE)
} else {
  full_log <- new_log
}
fwrite(full_log, log_file)
cat("\nHarmonisation log updated:", log_file, "\n")
cat("Done. Next: run scripts/09_protein_metabolite_mr.R\n")
sessionInfo()
