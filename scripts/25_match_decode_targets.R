#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
listing_path <- file.path(project_dir, "data", "decode", "decode_folder_listing.csv")
priority_path <- file.path(project_dir, "results", "validation", "replication_priority_targets.csv")
out_dir <- file.path(project_dir, "results", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(listing_path)) {
  stop("Missing deCODE folder listing: ", listing_path,
       "\nRun: Rscript scripts/24_decode_token_helper.R list <TOKEN>")
}

listing <- fread(listing_path)
priority <- fread(priority_path)

listing[, base_name := basename(Key)]

parse_decode_name <- function(x) {
  x <- sub("\\.txt\\.gz$", "", x, ignore.case = TRUE)
  parts <- strsplit(x, "_", fixed = TRUE)
  data.table(
    seq_id = vapply(parts, function(p) if (length(p) >= 1) p[[1]] else NA_character_, character(1)),
    gene_name = vapply(parts, function(p) if (length(p) >= 2) p[[2]] else NA_character_, character(1)),
    protein_label = vapply(parts, function(p) if (length(p) >= 3) paste(p[3:length(p)], collapse = "_") else NA_character_, character(1))
  )
}

parsed <- parse_decode_name(listing$base_name)
listing <- cbind(listing, parsed)

targets <- unique(priority$exposure)

support_files <- listing[grepl("^assocvariants\\.(annotated|excluded)\\.txt\\.gz$", base_name, ignore.case = TRUE)]
target_files <- listing[gene_name %in% targets]

setorder(target_files, gene_name, base_name)

fwrite(target_files, file.path(out_dir, "decode_target_file_candidates.csv"))
fwrite(support_files, file.path(out_dir, "decode_support_file_candidates.csv"))

cat("Wrote:\n")
cat("  results/validation/decode_target_file_candidates.csv\n")
cat("  results/validation/decode_support_file_candidates.csv\n\n")

if (nrow(target_files)) {
  print(target_files[, .(gene_name, seq_id, protein_label, Key, Size)])
} else {
  cat("No target-gene files found in current deCODE listing.\n")
}

if (nrow(support_files)) {
  cat("\nSupport files detected:\n")
  print(support_files[, .(base_name, Key, Size)])
}

