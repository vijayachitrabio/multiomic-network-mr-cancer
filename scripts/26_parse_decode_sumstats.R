#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/26_parse_decode_sumstats.R <APTAMER_FILE> <ANNOTATED_FILE> <EXCLUDED_FILE> <OUTPUT_CSV>\n",
    sep = ""
  )
  quit(status = 1)
}

if (length(args) < 4) usage()

aptamer_file <- args[[1]]
annotated_file <- args[[2]]
excluded_file <- args[[3]]
output_csv <- args[[4]]

if (!file.exists(aptamer_file)) stop("Missing aptamer file: ", aptamer_file)
if (!file.exists(annotated_file)) stop("Missing annotated file: ", annotated_file)
if (!file.exists(excluded_file)) stop("Missing excluded file: ", excluded_file)

message("Reading aptamer summary statistics...")
x <- fread(aptamer_file)

required_cols <- c("Chrom", "Pos", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval", "SE", "N", "ImpMAF")
missing_cols <- setdiff(required_cols, names(x))
if (length(missing_cols)) {
  stop("Aptamer file missing required columns: ", paste(missing_cols, collapse = ", "))
}

message("Reading deCODE extra annotation...")
ann <- fread(annotated_file)
ann_required <- c("Name", "effectAllele", "otherAllele", "effectAlleleFreq")
ann_missing <- setdiff(ann_required, names(ann))
if (length(ann_missing)) {
  stop("Annotated file missing required columns: ", paste(ann_missing, collapse = ", "))
}

message("Reading excluded variants...")
ex <- fread(excluded_file)
if (!"Name" %in% names(ex)) stop("Excluded file missing Name column.")

x <- x[!Name %in% ex$Name]
ann <- unique(ann[, .(Name, effectAllele_ann = effectAllele, otherAllele_ann = otherAllele, effectAlleleFreq)])
x <- merge(x, ann, by = "Name", all.x = TRUE)

x[, effectAllele_final := fifelse(!is.na(effectAllele_ann), effectAllele_ann, effectAllele)]
x[, otherAllele_final := fifelse(!is.na(otherAllele_ann), otherAllele_ann, otherAllele)]

parsed_name <- strsplit(basename(aptamer_file), "_", fixed = TRUE)[[1]]
seq_id <- if (length(parsed_name) >= 1) parsed_name[[1]] else NA_character_
gene_name <- if (length(parsed_name) >= 2) parsed_name[[2]] else NA_character_
protein_label <- if (length(parsed_name) >= 3) paste(parsed_name[3:length(parsed_name)], collapse = "_") else NA_character_
protein_label <- sub("\\.txt\\.gz$", "", protein_label, ignore.case = TRUE)

out <- x[, .(
  seq_id = seq_id,
  gene_name = gene_name,
  protein_label = protein_label,
  chr = as.integer(Chrom),
  pos = as.integer(Pos),
  variant_name = Name,
  rsid = rsids,
  effect_allele = effectAllele_final,
  other_allele = otherAllele_final,
  effect_allele_freq = as.numeric(effectAlleleFreq),
  beta = as.numeric(Beta),
  se = as.numeric(SE),
  pval = as.numeric(Pval),
  n = as.numeric(N),
  imp_maf = as.numeric(ImpMAF),
  min_log10_pval = suppressWarnings(as.numeric(min_log10_pval))
)]

out <- out[!is.na(beta) & !is.na(se) & se > 0]

fwrite(out, output_csv)
message("Wrote parsed deCODE summary statistics to: ", output_csv)

