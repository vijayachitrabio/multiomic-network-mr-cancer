#!/usr/bin/env Rscript

# Script 03: Download cancer GWAS summary statistics
#
# Female hormone-sensitive cancers only (prostate excluded by design):
#   1. Endometrial cancer  — O'Mara et al 2018, GWAS Catalog GCST006464
#                            12,906 cases / 108,979 controls (ECAC)
#   2. Ovarian cancer      — Phelan et al, GWAS Catalog GCST90016665
#                            25,509 cases / 40,138 controls (OCAC)
#   3. Breast cancer       — Zhang et al 2020 Nature Genetics, GCST90018757
#                            122,977 cases / 105,974 controls (BCAC)
#
# All files are downloaded as GWAS Catalog harmonised summary statistics
# (.h.tsv.gz) from the EBI FTP. The harmonised format is detected
# automatically by 04_harmonise_all.R.
#
# Re-running this script is safe — existing files are not re-downloaded.

set.seed(42)

if (!require("here",       quietly = TRUE)) install.packages("here",       repos = "https://cloud.r-project.org")
if (!require("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")

library(here)
library(data.table)

project_dir <- here::here()
cancer_dir  <- file.path(project_dir, "data", "cancer_gwas")
if (!dir.exists(cancer_dir)) dir.create(cancer_dir, recursive = TRUE)

# ------------------------------------------------------------------
# Cancer GWAS download table
# ------------------------------------------------------------------
cancers <- data.frame(
  label      = c("Endometrial", "Ovarian", "Breast"),
  accession  = c("GCST006464", "GCST90016665", "GCST90018757"),
  ftp_range  = c("GCST006001-GCST007000",
                 "GCST90016001-GCST90017000",
                 "GCST90018001-GCST90019000"),
  n_cases    = c(12906,  25509,  122977),
  n_controls = c(108979, 40138,  105974),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------
# Download helper
# ------------------------------------------------------------------
download_gwas <- function(label, accession, ftp_range, cancer_dir) {
  dest <- file.path(cancer_dir, paste0(label, "_", accession, ".h.tsv.gz"))

  if (file.exists(dest)) {
    cat(sprintf("[%s] Already present: %s (%.1f MB) — skipping download.\n",
                label, basename(dest), file.size(dest) / 1e6))
    return(dest)
  }

  url <- paste0(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
    ftp_range, "/", accession,
    "/harmonised/", accession, ".h.tsv.gz"
  )

  cat(sprintf("[%s] Downloading %s ...\n", label, accession))
  cat(sprintf("  URL: %s\n", url))
  options(timeout = 3600)

  ok <- tryCatch({
    download.file(url, destfile = dest, method = "curl", mode = "wb",
                  extra = "-L --retry 3 --retry-delay 5 --silent --show-error")
    TRUE
  }, error = function(e) {
    message("  curl failed: ", conditionMessage(e), " — retrying with wget...")
    tryCatch({
      download.file(url, destfile = dest, method = "wget", mode = "wb",
                    extra = "--tries=3 --timeout=300 --quiet")
      TRUE
    }, error = function(e2) {
      message("  wget also failed: ", conditionMessage(e2))
      FALSE
    })
  })

  if (ok && file.exists(dest)) {
    cat(sprintf("  Done: %.1f MB\n", file.size(dest) / 1e6))
  } else {
    cat(sprintf("  FAILED — download manually from:\n  %s\n  Save as: %s\n", url, dest))
  }
  dest
}

# ------------------------------------------------------------------
# Run downloads
# ------------------------------------------------------------------
cat("=== Downloading female hormone-sensitive cancer GWAS ===\n\n")
for (i in seq_len(nrow(cancers))) {
  download_gwas(cancers$label[i], cancers$accession[i],
                cancers$ftp_range[i], cancer_dir)
}

# ------------------------------------------------------------------
# Column check — confirm format detected by 04_harmonise_all.R
# ------------------------------------------------------------------
detect_format <- function(fpath, label) {
  if (!file.exists(fpath)) {
    cat(sprintf("[%s] MISSING\n", label))
    return(invisible(NULL))
  }
  hdr <- tryCatch(fread(fpath, nrows = 0), error = function(e) NULL)
  if (is.null(hdr)) { cat(sprintf("[%s] unreadable\n", label)); return(invisible(NULL)) }
  fmt <- if ("hm_rsid" %in% colnames(hdr)) "hm_rsid (harmonised)"
         else if ("rsid" %in% colnames(hdr)) "rsid"
         else "SNP"
  cat(sprintf("[%s] %d cols | SNP key: %s\n", label, ncol(hdr), fmt))
}

cat("\n--- Format check ---\n")
for (i in seq_len(nrow(cancers))) {
  f <- file.path(cancer_dir,
                 paste0(cancers$label[i], "_", cancers$accession[i], ".h.tsv.gz"))
  detect_format(f, cancers$label[i])
}

# ------------------------------------------------------------------
# Final summary
# ------------------------------------------------------------------
cat("\n--- Cancer GWAS files present ---\n")
all_files <- list.files(cancer_dir, pattern = "\\.tsv\\.gz$|\\.csv$", full.names = TRUE)
for (f in all_files) {
  cat(sprintf("  %-45s  %.1f MB\n", basename(f), file.size(f) / 1e6))
}

cat("\nScript 03 complete.\n")
cat("Next: run 04_harmonise_all.R\n")
sessionInfo()
