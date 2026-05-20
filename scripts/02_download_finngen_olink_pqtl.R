#!/usr/bin/env Rscript

# Script 01b: Download FinnGen R10 Olink pQTL cis instruments
#
# Uses the public FinnGen proteomics autoreporting "group_report" files, which
# contain grouped genome-wide significant pQTL lead variants. This is a compact
# fallback pQTL source while the larger UKB-PPP panel is unavailable.

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(curl)
})

project_dir <- "."
out_dir     <- file.path(project_dir, "data", "pqtl")
report_dir  <- file.path(out_dir, "finngen_olink_group_report")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
write_primary <- "--primary" %in% args

bucket <- "finngen-public-data-r10"
base_prefix <- "omics/proteomics/release_2023_03_02"
gcs_api <- sprintf("https://storage.googleapis.com/storage/v1/b/%s/o", bucket)
public_base <- sprintf("https://storage.googleapis.com/%s", bucket)

probe_map_url <- file.path(public_base, base_prefix, "data/Olink/probe_map.tsv")
probe_map_file <- file.path(out_dir, "finngen_olink_probe_map.tsv")

message("Downloading FinnGen Olink probe map...")
curl_download(probe_map_url, probe_map_file, quiet = TRUE, mode = "wb")

list_gcs_objects <- function(prefix, max_results = 1000) {
  out <- list()
  page_token <- NULL
  repeat {
    url <- paste0(
      gcs_api,
      "?prefix=", utils::URLencode(prefix, reserved = TRUE),
      "&maxResults=", max_results
    )
    if (!is.null(page_token)) {
      url <- paste0(url, "&pageToken=", utils::URLencode(page_token, reserved = TRUE))
    }
    tf <- tempfile(fileext = ".json")
    curl_download(url, tf, quiet = TRUE, mode = "wb")
    x <- fromJSON(tf, simplifyVector = FALSE)
    unlink(tf)
    if (!is.null(x$items)) out <- c(out, x$items)
    page_token <- x$nextPageToken
    if (is.null(page_token) || length(page_token) == 0) break
  }
  out
}

message("Listing FinnGen Olink grouped pQTL reports...")
group_prefix <- file.path(base_prefix, "data/Olink/autoreporting/group_report/")
objects <- list_gcs_objects(group_prefix)
manifest <- rbindlist(lapply(objects, function(x) {
  data.table(
    name = x$name,
    size = as.numeric(x$size),
    updated = x$updated
  )
}), fill = TRUE)

manifest <- manifest[grepl("\\.top\\.out$", name)]
manifest[, protein := sub("\\.top\\.out$", "", basename(name))]
manifest[, url := paste0(public_base, "/", name)]
manifest[, local_file := file.path(report_dir, basename(name))]
setorder(manifest, protein)

manifest_file <- file.path(out_dir, "finngen_olink_group_report_manifest.csv")
fwrite(manifest, manifest_file)
message(sprintf("Found %d grouped report files.", nrow(manifest)))

needs_download <- manifest[!file.exists(local_file) | file.info(local_file)$size == 0]
if (nrow(needs_download) > 0) {
  message(sprintf("Downloading %d missing grouped report files...", nrow(needs_download)))
  dl <- multi_download(
    urls = needs_download$url,
    destfiles = needs_download$local_file,
    resume = TRUE,
    progress = TRUE
  )
  dl <- as.data.table(dl)
  failed <- dl[success == FALSE]
  if (nrow(failed) > 0) {
    stop("Some FinnGen grouped report downloads failed; first failed URL: ", failed$url[1])
  }
} else {
  message("All grouped report files already present.")
}

message("Reading grouped pQTL reports...")
read_report <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  required <- c(
    "phenotype", "locus_id", "chrom", "pos", "ref", "alt", "pval",
    "lead_BETA", "lead_SE", "lead_ALT_FREQ"
  )
  if (!all(required %in% names(dt))) return(NULL)
  dt
}

reports <- rbindlist(lapply(manifest$local_file, read_report), fill = TRUE)
if (nrow(reports) == 0) stop("No valid grouped pQTL report rows were read.")

message(sprintf("Raw genome-wide significant grouped rows: %d", nrow(reports)))

probe_map <- fread(probe_map_file)
probe_map[, chr_clean := sub("^chr", "", chr)]
probe_map[, gene_start := as.integer(start)]
probe_map[, gene_end := as.integer(end)]
probe_map <- unique(probe_map[, .(geneName, chr_clean, gene_start, gene_end)])

reports[, chrom_clean := as.character(chrom)]
reports[, pos := as.integer(pos)]
reports[, pval := as.numeric(pval)]
reports[, lead_BETA := as.numeric(lead_BETA)]
reports[, lead_SE := as.numeric(lead_SE)]
reports[, lead_ALT_FREQ := as.numeric(lead_ALT_FREQ)]

pqtl <- merge(
  reports,
  probe_map,
  by.x = "phenotype",
  by.y = "geneName",
  all.x = FALSE,
  allow.cartesian = TRUE
)

cis_window_bp <- 1e6
pqtl <- pqtl[
  chrom_clean == chr_clean &
    pos >= (gene_start - cis_window_bp) &
    pos <= (gene_end + cis_window_bp)
]
pqtl <- pqtl[
  pval < 5e-8 &
    !is.na(lead_BETA) &
    !is.na(lead_SE) &
    lead_SE > 0
]
pqtl <- pqtl[!(chrom_clean == "6" & pos >= 28000000 & pos <= 34000000)]
pqtl[, F_stat := (lead_BETA^2) / (lead_SE^2)]
pqtl <- pqtl[F_stat > 10]

make_variant_key <- function(chr, pos, allele1, allele2) {
  chr <- sub("^chr", "", as.character(chr))
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  a_min <- ifelse(allele1 <= allele2, allele1, allele2)
  a_max <- ifelse(allele1 <= allele2, allele2, allele1)
  paste(chr, pos, a_min, a_max, sep = ":")
}

pqtl[, variant_key := make_variant_key(chrom_clean, pos, ref, alt)]
setorder(pqtl, phenotype, pval)
pqtl <- pqtl[!duplicated(paste(phenotype, variant_key))]

out <- pqtl[, .(
  protein = phenotype,
  SNP = variant_key,
  variant_id.exposure = locus_id,
  chr.exposure = chrom_clean,
  pos.exposure = pos,
  effect_allele.exposure = toupper(alt),
  other_allele.exposure = toupper(ref),
  eaf.exposure = lead_ALT_FREQ,
  beta.exposure = lead_BETA,
  se.exposure = lead_SE,
  pval.exposure = pval,
  exposure = phenotype,
  id.exposure = paste0("FINNGEN_OLINK_", phenotype),
  units.exposure = "IRNT_Olink_NPX",
  phenotype.col = "phenotype",
  data_source = "FinnGen_R10_Olink_autoreporting_group_report",
  mr_keep = TRUE,
  pQTL_source = "FinnGen_R10_Olink",
  genome_build.exposure = "GRCh38",
  cis_window_mb = 1,
  N = 619,
  samplesize.exposure = 619,
  F_stat = round(F_stat, 2),
  clump_source = "FinnGen_autoreporting_grouped_GWS_leads"
)]

setorder(out, protein, pval.exposure)

out_file <- file.path(out_dir, "pqtl_instruments_finngen_olink.csv")
fwrite(out, out_file)
saveRDS(out, file.path(out_dir, "pqtl_instruments_finngen_olink.rds"))

summary_file <- file.path(out_dir, "pqtl_instruments_finngen_olink_summary.csv")
summary <- out[, .(
  n_instruments = .N,
  min_p = min(pval.exposure, na.rm = TRUE),
  mean_F = round(mean(F_stat, na.rm = TRUE), 2),
  min_F = round(min(F_stat, na.rm = TRUE), 2),
  max_F = round(max(F_stat, na.rm = TRUE), 2)
), by = protein]
setorder(summary, protein)
fwrite(summary, summary_file)

message(sprintf("Cis instruments retained: %d across %d proteins.", nrow(out), uniqueN(out$protein)))
message(sprintf("Wrote: %s", out_file))

if (write_primary) {
  primary_file <- file.path(out_dir, "pqtl_instruments.csv")
  if (file.exists(primary_file)) {
    backup_file <- file.path(
      out_dir,
      sprintf("pqtl_instruments_backup_before_finngen_%s.csv", format(Sys.Date(), "%Y-%m-%d"))
    )
    if (!file.exists(backup_file)) file.copy(primary_file, backup_file)
    message(sprintf("Backed up previous primary pQTL file to: %s", backup_file))
  }
  fwrite(out, primary_file)
  message(sprintf("Updated primary pQTL instrument file: %s", primary_file))
}

message("Done.")
sessionInfo()
