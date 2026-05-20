#!/usr/bin/env Rscript

# Script 38: ARIC pQTL replication MR for protein -> cancer hits
#
# This script is intentionally non-overwriting:
#   - reads an ARIC pQTL instrument CSV
#   - harmonises only the requested discovery hit proteins/cancers
#   - writes results under results/replication/aric_<source>/

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

project_dir <- "."

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit) || hit == length(args)) return(default)
  args[[hit + 1L]]
}

pqtl_file <- get_arg(
  "--pqtl-file",
  file.path(project_dir, "data", "pqtl", "pqtl_instruments_aric_ea_cis.csv")
)
target_file <- get_arg(
  "--targets",
  file.path(project_dir, "results", "tables", "STable2_17_FDR_hits_complete.csv")
)
source_label <- get_arg("--source-label", sub("^pqtl_instruments_", "", sub("\\.csv$", "", basename(pqtl_file))))

out_dir <- file.path(project_dir, "results", "replication", paste0("aric_", source_label))
harm_dir <- file.path(out_dir, "harmonised")
dir.create(harm_dir, recursive = TRUE, showWarnings = FALSE)

make_variant_key <- function(chr, pos, allele1, allele2) {
  chr <- sub("^chr", "", as.character(chr))
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  a_min <- ifelse(allele1 <= allele2, allele1, allele2)
  a_max <- ifelse(allele1 <= allele2, allele2, allele1)
  paste(chr, pos, a_min, a_max, sep = ":")
}

clean_name <- function(path) {
  x <- basename(path)
  x <- sub("\\.h\\.tsv\\.gz$", "", x)
  x <- sub("\\.tsv\\.gz$", "", x)
  x
}

normalise_cancer_name <- function(x) {
  x <- sub("\\.h\\.tsv\\.gz$", "", basename(x))
  x <- sub("\\.tsv\\.gz$", "", x)
  x
}

read_targets <- function(path) {
  if (!file.exists(path)) stop("Missing target file: ", path)
  x <- fread(path)
  if (all(c("protein", "cancer") %in% names(x))) {
    out <- unique(x[, .(protein, outcome = cancer)])
  } else if (all(c("exposure", "outcome") %in% names(x))) {
    out <- unique(x[, .(protein = exposure, outcome)])
  } else {
    stop("Target file must contain protein/cancer or exposure/outcome columns.")
  }
  out[!is.na(protein) & protein != "" & !is.na(outcome) & outcome != ""]
}

read_discovery <- function(path) {
  x <- fread(path)
  col_or_na <- function(col) {
    if (col %in% names(x)) x[[col]] else rep(NA_real_, nrow(x))
  }
  if (all(c("protein", "cancer") %in% names(x))) {
    setnames(x, c("protein", "cancer"), c("discovery_protein", "discovery_outcome"), skip_absent = TRUE)
    x[, .(
      protein = discovery_protein,
      outcome = discovery_outcome,
      discovery_method = method,
      discovery_beta = beta,
      discovery_se = se,
      discovery_pval = col_or_na("pvalue"),
      discovery_fdr = col_or_na("FDR"),
      discovery_or = col_or_na("OR")
    )]
  } else {
    data.table()
  }
}

if (!file.exists(pqtl_file)) {
  stop(
    "Missing ARIC pQTL file: ", pqtl_file, "\n",
    "Run: Rscript scripts/37_prepare_aric_pqtl_instruments.R"
  )
}

message("Loading ARIC pQTL instruments: ", pqtl_file)
pqtls <- fread(pqtl_file)
required_exp <- c(
  "SNP", "effect_allele.exposure", "other_allele.exposure",
  "beta.exposure", "se.exposure", "pval.exposure", "exposure", "id.exposure"
)
missing_exp <- setdiff(required_exp, names(pqtls))
if (length(missing_exp)) stop("pQTL file missing columns: ", paste(missing_exp, collapse = ", "))

pqtls <- pqtls[mr_keep == TRUE]
pqtls[, SNP := toupper(SNP)]
message(sprintf("  %d ARIC cis instruments across %d proteins", nrow(pqtls), uniqueN(pqtls$exposure)))

targets <- read_targets(target_file)
targets[, outcome_clean := normalise_cancer_name(outcome)]
message(sprintf("Loaded %d replication target protein/cancer pairs.", nrow(targets)))

discovery <- read_discovery(target_file)
if (nrow(discovery)) discovery[, outcome_clean := normalise_cancer_name(outcome)]

cancer_dir <- file.path(project_dir, "data", "cancer_gwas")
cancer_files <- list.files(cancer_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
cancer_map <- data.table(path = cancer_files)
cancer_map[, outcome_clean := clean_name(path)]

targets <- merge(targets, cancer_map, by = "outcome_clean", all.x = TRUE, sort = FALSE)
if (anyNA(targets$path)) {
  missing <- targets[is.na(path), unique(outcome)]
  stop("No local cancer GWAS file matched: ", paste(missing, collapse = ", "))
}

mr_results <- list()
coverage <- list()
harm_logs <- list()

for (outcome_name in unique(targets$outcome_clean)) {
  target_subset <- targets[outcome_clean == outcome_name]
  proteins <- unique(target_subset$protein)
  exp_dat <- pqtls[exposure %in% proteins]
  exp_dt <- copy(exp_dat)

  if (nrow(exp_dat) == 0) {
    coverage[[outcome_name]] <- data.table(
      protein = proteins,
      outcome = outcome_name,
      status = "no_aric_cis_instrument",
      n_aric_instruments = 0L,
      n_harmonised = 0L,
      n_mr_keep = 0L
    )
    next
  }

  gwas_path <- target_subset$path[1]
  message("\n=== ", outcome_name, " ===")
  message("Reading GWAS: ", basename(gwas_path))
  gwas <- fread(gwas_path)

  has_hm <- all(c("hm_chrom", "hm_pos", "hm_effect_allele", "hm_other_allele") %in% names(gwas))
  has_raw <- all(c("chromosome", "base_pair_location", "effect_allele", "other_allele") %in% names(gwas))
  if (has_hm) {
    gwas[, SNP_variant_key := make_variant_key(hm_chrom, hm_pos, hm_effect_allele, hm_other_allele)]
    ea_col <- "hm_effect_allele"
    oa_col <- "hm_other_allele"
    beta_col <- if ("hm_beta" %in% names(gwas)) "hm_beta" else "beta"
    se_col <- "standard_error"
    eaf_col <- "hm_effect_allele_frequency"
    pval_col <- "p_value"
  } else if (has_raw) {
    gwas[, SNP_variant_key := make_variant_key(chromosome, base_pair_location, effect_allele, other_allele)]
    ea_col <- "effect_allele"
    oa_col <- "other_allele"
    beta_col <- "beta"
    se_col <- "standard_error"
    eaf_col <- "effect_allele_frequency"
    pval_col <- "p_value"
  } else {
    stop("GWAS file lacks positional allele columns: ", gwas_path)
  }
  gwas[, SNP_variant_key := toupper(SNP_variant_key)]
  gwas <- gwas[SNP_variant_key %in% exp_dat$SNP]
  for (num_col in c(beta_col, se_col, eaf_col, pval_col)) {
    if (num_col %in% names(gwas)) gwas[, (num_col) := as.numeric(get(num_col))]
  }
  message(sprintf("  Matched outcome rows: %d", nrow(gwas)))

  if (nrow(gwas) == 0) {
    inst_counts <- exp_dt[, .(n_aric_instruments = .N), by = .(protein = exposure)]
    coverage[[outcome_name]] <- data.table(
      protein = proteins,
      outcome = outcome_name,
      status = fifelse(proteins %in% exp_dat$exposure, "no_gwas_overlap", "no_aric_cis_instrument"),
      n_aric_instruments = inst_counts[match(proteins, protein), n_aric_instruments],
      n_harmonised = 0L,
      n_mr_keep = 0L
    )
    coverage[[outcome_name]][is.na(n_aric_instruments), n_aric_instruments := 0L]
    next
  }

  setDF(gwas)
  out_dat <- format_data(
    gwas,
    type = "outcome",
    snp_col = "SNP_variant_key",
    beta_col = beta_col,
    se_col = se_col,
    eaf_col = eaf_col,
    effect_allele_col = ea_col,
    other_allele_col = oa_col,
    pval_col = pval_col
  )
  out_dat$SNP <- toupper(out_dat$SNP)
  out_dat$outcome <- outcome_name
  out_dat$id.outcome <- outcome_name

  setDF(exp_dat)
  harm <- tryCatch(
    harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2),
    error = function(e) {
      message("  harmonise_data error: ", conditionMessage(e))
      data.frame()
    }
  )
  harm_path <- file.path(harm_dir, paste0("harmonised_aric_", outcome_name, ".rds"))
  saveRDS(harm, harm_path)

  if (nrow(harm) == 0) {
    coverage[[outcome_name]] <- data.table(
      protein = proteins,
      outcome = outcome_name,
      status = "harmonise_empty",
      n_aric_instruments = as.integer(table(factor(exp_dat$exposure, levels = proteins))),
      n_harmonised = 0L,
      n_mr_keep = 0L
    )
    next
  }

  harm_dt <- as.data.table(harm)
  harm_dt <- harm_dt[mr_keep == TRUE & !is.na(mr_keep)]
  harm_dt <- harm_dt[!is.na(beta.exposure) & !is.na(se.exposure) & se.exposure > 0]
  harm_dt[, F_stat := beta.exposure^2 / se.exposure^2]
  harm_dt <- harm_dt[F_stat > 10]

  cov <- merge(
    data.table(protein = proteins),
    exp_dt[, .(n_aric_instruments = .N), by = .(protein = exposure)],
    by = "protein", all.x = TRUE
  )
  cov[is.na(n_aric_instruments), n_aric_instruments := 0L]
  cov <- merge(
    cov,
    as.data.table(harm)[, .(n_harmonised = .N), by = .(protein = exposure)],
    by = "protein", all.x = TRUE
  )
  cov <- merge(
    cov,
    harm_dt[, .(
      n_mr_keep = .N,
      aric_snps = paste(unique(SNP), collapse = ";"),
      variant_types = paste(unique(na.omit(variant_type)), collapse = ";")
    ), by = .(protein = exposure)],
    by = "protein", all.x = TRUE
  )
  cov[is.na(n_harmonised), n_harmonised := 0L]
  cov[is.na(n_mr_keep), n_mr_keep := 0L]
  cov[, outcome := outcome_name]
  cov[, status := fifelse(
    n_aric_instruments == 0, "no_aric_cis_instrument",
    fifelse(n_mr_keep == 0, "no_mr_keep_after_harmonise", "ready")
  )]
  coverage[[outcome_name]] <- cov

  if (nrow(harm_dt) == 0) next

  mr <- tryCatch(
    mr(
      as.data.frame(harm_dt),
      method_list = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median")
    ),
    error = function(e) {
      message("  MR error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(mr) || nrow(mr) == 0) next

  mr <- as.data.table(generate_odds_ratios(mr))
  qc <- harm_dt[, .(
    F_stat_mean = mean(F_stat, na.rm = TRUE),
    n_instruments = .N,
    aric_snps = paste(unique(SNP), collapse = ";"),
    variant_types = paste(unique(na.omit(variant_type)), collapse = ";")
  ), by = .(id.exposure, exposure)]
  mr <- merge(mr, qc, by = c("id.exposure", "exposure"), all.x = TRUE, sort = FALSE)
  mr[, source_label := source_label]
  mr_results[[outcome_name]] <- mr

  harm_logs[[outcome_name]] <- data.table(
    outcome = outcome_name,
    n_exposure_rows = nrow(exp_dat),
    n_outcome_matches = nrow(gwas),
    n_harmonised = nrow(harm),
    n_mr_keep = nrow(harm_dt)
  )
}

mr_all <- rbindlist(mr_results, fill = TRUE)
coverage_all <- rbindlist(coverage, fill = TRUE)
harm_log <- rbindlist(harm_logs, fill = TRUE)

if (nrow(mr_all)) {
  mr_all[, fdr := p.adjust(pval, method = "BH"), by = outcome]
  if (nrow(discovery)) {
    mr_all <- merge(
      mr_all,
      discovery[, .(
        exposure = protein,
        outcome_clean,
        discovery_method,
        discovery_beta,
        discovery_se,
        discovery_pval,
        discovery_fdr,
        discovery_or
      )],
      by.x = c("exposure", "outcome"),
      by.y = c("exposure", "outcome_clean"),
      all.x = TRUE,
      sort = FALSE
    )
    mr_all[, direction_concordant := sign(b) == sign(discovery_beta)]
  }
}

fwrite(mr_all, file.path(out_dir, "aric_protein_cancer_mr_results.csv"))
fwrite(coverage_all, file.path(out_dir, "aric_replication_coverage.csv"))
fwrite(harm_log, file.path(out_dir, "aric_harmonisation_log.csv"))

message("\nARIC replication complete.")
message("Wrote: ", file.path(out_dir, "aric_protein_cancer_mr_results.csv"))
message("Wrote: ", file.path(out_dir, "aric_replication_coverage.csv"))
if (nrow(coverage_all)) print(coverage_all[order(status, protein)])
if (nrow(mr_all)) print(mr_all[order(pval), .(exposure, outcome, method, nsnp, b, se, pval, fdr, or, direction_concordant)])
sessionInfo()
