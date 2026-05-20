#!/usr/bin/env Rscript

# Script 37: Prepare ARIC SomaScan pQTL instruments
#
# Source:
#   Sarnowski/Ma et al. ARIC plasma proteome pQTL supplement
#   PMC11469718.1/media-2.xlsx (Tables S5/S6)
#
# Outputs are kept separate from the active FinnGen pQTL file. Use --primary only
# if you intentionally want to replace data/pqtl/pqtl_instruments.csv.

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

project_dir <- "."
out_dir <- file.path(project_dir, "data", "pqtl")
aric_dir <- file.path(out_dir, "aric")
dir.create(aric_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
write_primary <- "--primary" %in% args
use_mr_selected_primary <- "--mr-selected-primary" %in% args

supplement_url <- "https://pmc-oa-opendata.s3.amazonaws.com/PMC11469718.1/media-2.xlsx"
supplement_file <- file.path(aric_dir, "sarnowski_2024_supplement2.xlsx")

if (!file.exists(supplement_file) || file.info(supplement_file)$size < 1e6) {
  message("Downloading ARIC pQTL supplementary workbook...")
  utils::download.file(supplement_url, supplement_file, mode = "wb", quiet = FALSE)
}

make_variant_key <- function(chr, pos, allele1, allele2) {
  chr <- sub("^chr", "", as.character(chr))
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  a_min <- ifelse(allele1 <= allele2, allele1, allele2)
  a_max <- ifelse(allele1 <= allele2, allele2, allele1)
  paste(chr, pos, a_min, a_max, sep = ":")
}

read_aric_table <- function(sheet, ancestry, n) {
  message("Reading ", sheet, " (", ancestry, ")...")
  x <- as.data.table(suppressWarnings(read_excel(supplement_file, sheet = sheet, skip = 1)))
  needed <- c(
    "SeqID", "Target", "UniProt", "EntrezGeneSymbol", "cis_trans",
    "Sentinel_variant", "rs_dbSNP150", "Conditional.variant",
    "b", "se", "p", "freq", "VariantType"
  )
  missing <- setdiff(needed, names(x))
  if (length(missing)) {
    stop(sheet, " missing required columns: ", paste(missing, collapse = ", "))
  }

  x[, ancestry := ancestry]
  x[, samplesize := n]
  if (!"SelectedForMR" %in% names(x)) x[, SelectedForMR := NA_character_]

  v <- tstrsplit(x$Conditional.variant, ":", fixed = TRUE)
  if (length(v) != 4) stop("Unexpected variant format in ", sheet)
  x[, chr := sub("^chr", "", v[[1]])]
  x[, pos := as.integer(v[[2]])]
  x[, ref := toupper(v[[3]])]
  x[, alt := toupper(v[[4]])]

  x[, beta := as.numeric(b)]
  x[, se_num := as.numeric(se)]
  x[, pval := as.numeric(p)]
  x[, eaf := as.numeric(freq)]
  x[, F_stat := (beta^2) / (se_num^2)]
  x[, is_mhc := chr == "6" & pos >= 28000000 & pos <= 34000000]
  x[, is_cis := tolower(cis_trans) == "cis"]
  x[, is_mr_selected := toupper(as.character(SelectedForMR)) == "YES"]
  x[is.na(is_mr_selected), is_mr_selected := FALSE]

  x[, SNP := make_variant_key(chr, pos, ref, alt)]
  x[, exposure := EntrezGeneSymbol]
  x <- x[
    !is.na(exposure) & exposure != "" &
      !is.na(beta) & !is.na(se_num) & se_num > 0 &
      !is.na(pos) & !is.na(ref) & !is.na(alt) &
      !is.na(pval)
  ]

  x[, .(
    protein = exposure,
    SNP,
    variant_id.exposure = Conditional.variant,
    chr.exposure = chr,
    pos.exposure = pos,
    effect_allele.exposure = alt,
    other_allele.exposure = ref,
    eaf.exposure = eaf,
    beta.exposure = beta,
    se.exposure = se_num,
    pval.exposure = pval,
    exposure,
    id.exposure = paste("ARIC", ancestry, "SOMASCAN", exposure, sep = "_"),
    units.exposure = "SomaScan_normalized_plasma_protein",
    phenotype.col = "EntrezGeneSymbol",
    data_source = paste0("ARIC_SomaScan_", ancestry, "_Sarnowski2024_supplement_Table_", sub("^Table_S", "S", sheet)),
    mr_keep = is_cis & !is_mhc & F_stat > 10,
    aric_mr_selected = is_mr_selected,
    pQTL_source = paste0("ARIC_SomaScan_", ancestry),
    genome_build.exposure = "GRCh38",
    cis_window_mb = fifelse(is_cis, 1, NA_real_),
    N = samplesize,
    samplesize.exposure = samplesize,
    F_stat = round(F_stat, 2),
    clump_source = "ARIC_GCTA_COJO_conditional_pQTL",
    ancestry,
    seq_id = SeqID,
    target_label = Target,
    uniprot = UniProt,
    rsid = rs_dbSNP150,
    sentinel_variant = Sentinel_variant,
    cis_trans = cis_trans,
    variant_type = VariantType,
    selected_for_mr = SelectedForMR
  )]
}

aa <- read_aric_table("Table_S5", "AA", 1871L)
ea <- read_aric_table("Table_S6", "EA", 7584L)

write_set <- function(x, path) {
  setorder(x, protein, pval.exposure, SNP)
  fwrite(x, path)
  saveRDS(x, sub("\\.csv$", ".rds", path))
  message("Wrote: ", path)
}

ea_all <- copy(ea)
aa_all <- copy(aa)
ea_cis <- ea_all[mr_keep == TRUE]
aa_cis <- aa_all[mr_keep == TRUE]
ea_mr_selected <- ea_cis[aric_mr_selected == TRUE]

write_set(ea_all, file.path(out_dir, "pqtl_instruments_aric_ea_all_conditional.csv"))
write_set(aa_all, file.path(out_dir, "pqtl_instruments_aric_aa_all_conditional.csv"))
write_set(ea_cis, file.path(out_dir, "pqtl_instruments_aric_ea_cis.csv"))
write_set(aa_cis, file.path(out_dir, "pqtl_instruments_aric_aa_cis.csv"))
write_set(ea_mr_selected, file.path(out_dir, "pqtl_instruments_aric_ea_cis_mr_selected.csv"))

summary <- rbindlist(list(
  ea_all[, .(source = "ARIC_EA_all_conditional", n_instruments = .N, n_proteins = uniqueN(protein), n_cis = sum(mr_keep))],
  ea_cis[, .(source = "ARIC_EA_cis", n_instruments = .N, n_proteins = uniqueN(protein), n_cis = .N)],
  ea_mr_selected[, .(source = "ARIC_EA_cis_mr_selected", n_instruments = .N, n_proteins = uniqueN(protein), n_cis = .N)],
  aa_all[, .(source = "ARIC_AA_all_conditional", n_instruments = .N, n_proteins = uniqueN(protein), n_cis = sum(mr_keep))],
  aa_cis[, .(source = "ARIC_AA_cis", n_instruments = .N, n_proteins = uniqueN(protein), n_cis = .N)]
), fill = TRUE)
fwrite(summary, file.path(out_dir, "pqtl_instruments_aric_summary.csv"))

if (write_primary) {
  primary <- if (use_mr_selected_primary) ea_mr_selected else ea_cis
  primary_file <- file.path(out_dir, "pqtl_instruments.csv")
  if (file.exists(primary_file)) {
    backup_file <- file.path(
      out_dir,
      sprintf("pqtl_instruments_backup_before_aric_%s.csv", format(Sys.Date(), "%Y-%m-%d"))
    )
    if (!file.exists(backup_file)) file.copy(primary_file, backup_file)
    message("Backed up previous primary pQTL file to: ", backup_file)
  }
  fwrite(primary, primary_file)
  saveRDS(primary, sub("\\.csv$", ".rds", primary_file))
  message("Updated primary pQTL instrument file: ", primary_file)
}

message("ARIC instrument preparation complete.")
print(summary)
sessionInfo()
