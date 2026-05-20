#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Rsamtools)
  library(GenomicRanges)
})

dir.create("data/pqtl/priority_regions", recursive = TRUE, showWarnings = FALSE)
dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)

WINDOW_BP <- as.integer(Sys.getenv("PQTL_WINDOW_BP", "500000"))
base_url <- "https://storage.googleapis.com/finngen-public-data-r10/omics/proteomics/release_2023_03_02/data/Olink/pQTL"

top_paths <- fread("results/validation/top_mediation_paths_for_followup.csv")
priority_proteins <- unique(top_paths$protein)

pqtl_inst <- fread("data/pqtl/pqtl_instruments.csv")
pqtl_inst <- pqtl_inst[protein %in% priority_proteins]

read_remote_region <- function(prot, chr, start, end) {
  url <- sprintf("%s/Olink_Batch1_%s.txt.gz", base_url, prot)
  region <- GRanges(as.character(chr), IRanges(start, end))
  lines <- scanTabix(TabixFile(url), param = region)[[1]]
  if (length(lines) == 0) return(data.table())

  x <- fread(
    text = paste(lines, collapse = "\n"),
    header = FALSE
  )
  cols <- c(
      "chr", "pos", "variant_id", "ref", "alt", "alt_freq", "beta",
      "se", "t_stat", "p", "log10_p", "n"
  )
  if (ncol(x) == 11) {
    cols <- cols[-12]
  }
  setnames(x, cols)
  if (!"n" %in% names(x)) x[, n := 619L]
  x[, protein := prot]
  x[, region := sprintf("%s:%s-%s", chr, start, end)]
  x
}

manifest <- list()

for (prot in priority_proteins) {
  leads <- pqtl_inst[protein == prot]
  if (nrow(leads) == 0) next

  protein_regions <- list()
  for (i in seq_len(nrow(leads))) {
    chr <- leads$chr.exposure[i]
    start <- max(1L, leads$pos.exposure[i] - WINDOW_BP)
    end <- leads$pos.exposure[i] + WINDOW_BP
    cat(sprintf("Extracting %s chr%s:%s-%s\n", prot, chr, start, end))
    region <- read_remote_region(prot, chr, start, end)
    if (nrow(region) > 0) protein_regions[[length(protein_regions) + 1]] <- region
  }

  if (length(protein_regions) == 0) next

  out <- unique(rbindlist(protein_regions, fill = TRUE), by = c("chr", "pos", "ref", "alt"))
  setorder(out, chr, pos)
  out_file <- file.path("data/pqtl/priority_regions", paste0(prot, "_pqtl_regions.tsv.gz"))
  fwrite(out, out_file, sep = "\t")

  manifest[[length(manifest) + 1]] <- data.table(
    protein = prot,
    n_regions_requested = nrow(leads),
    n_variants = nrow(out),
    file = out_file
  )
}

manifest <- rbindlist(manifest, fill = TRUE)
fwrite(manifest, "results/validation/priority_pqtl_region_manifest.csv")

cat("Wrote priority pQTL regional files:\n")
print(manifest)
