#!/usr/bin/env Rscript

# OpenGWAS pQTL replication MR for prioritized proteins with available datasets.

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
inst_file <- file.path(project_dir, "results", "opengwas", "opengwas_7protein_instruments_p5e8.csv")
out_dir <- file.path(project_dir, "results", "opengwas")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

instrument_map <- data.table(
  id.exposure = c(
    "ebi-a-GCST90019469",
    "prot-a-3052",
    "prot-a-131",
    "prot-a-132",
    "prot-a-1524",
    "ebi-a-GCST90000469"
  ),
  protein = c("ITIH3", "TNFRSF6B", "APOE_E3", "APOE_E2", "IL34", "FGF5")
)

outcomes <- data.table(
  outcome_id = c("Breast_GCST90018757", "Endometrial_GCST006464", "Ovarian_GCST90016665"),
  outcome = c("Breast", "Endometrial", "Ovarian"),
  file = file.path(
    project_dir,
    "data",
    "cancer_gwas",
    c(
      "Breast_GCST90018757.h.tsv.gz",
      "Endometrial_GCST006464.h.tsv.gz",
      "Ovarian_GCST90016665.h.tsv.gz"
    )
  )
)

if (!file.exists(inst_file)) {
  stop("Missing OpenGWAS instrument file: ", inst_file)
}

inst <- fread(inst_file)
inst <- merge(inst, instrument_map, by = "id.exposure", all.x = TRUE, sort = FALSE)
inst <- inst[!is.na(protein)]
inst[, F_stat := beta.exposure^2 / se.exposure^2]
inst <- inst[F_stat > 10]

format_exposure <- function(x) {
  format_data(
    as.data.frame(x),
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    pval_col = "pval.exposure",
    samplesize_col = "samplesize.exposure",
    phenotype_col = "protein",
    id_col = "id.exposure"
  )
}

read_outcome_for_snps <- function(path, snps) {
  gwas <- fread(path)
  rs_col <- if ("rsid" %in% names(gwas)) "rsid" else if ("hm_rsid" %in% names(gwas)) "hm_rsid" else NA_character_
  if (is.na(rs_col)) stop("No rsid/hm_rsid column in ", path)
  gwas <- gwas[get(rs_col) %in% snps]
  if (nrow(gwas) == 0) return(gwas)
  gwas[, SNP_match := get(rs_col)]
  gwas
}

format_outcome <- function(x, outcome_name, outcome_id) {
  format_data(
    as.data.frame(x),
    type = "outcome",
    snp_col = "SNP_match",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value",
    phenotype_col = "outcome_name",
    id_col = "outcome_id"
  )
}

harm_all <- list()
mr_all <- list()
coverage <- list()

for (i in seq_len(nrow(outcomes))) {
  oc <- outcomes[i]
  message("Processing ", oc$outcome)
  out_raw <- read_outcome_for_snps(oc$file, unique(inst$SNP))
  coverage[[oc$outcome]] <- data.table(
    outcome = oc$outcome,
    requested_snps = uniqueN(inst$SNP),
    matched_snps = uniqueN(out_raw$SNP_match)
  )
  if (nrow(out_raw) == 0) next

  out_raw[, outcome_name := oc$outcome]
  out_raw[, outcome_id := oc$outcome_id]
  out_dat <- format_outcome(out_raw, oc$outcome, oc$outcome_id)

  for (pid in unique(inst$id.exposure)) {
    exp_raw <- inst[id.exposure == pid]
    exp_dat <- format_exposure(exp_raw)
    harm <- as.data.table(harmonise_data(exp_dat, out_dat, action = 2))
    if ("mr_keep" %in% names(harm)) {
      harm <- harm[mr_keep == TRUE]
    }
    if (nrow(harm) == 0) next
    harm_all[[paste(pid, oc$outcome, sep = "__")]] <- as.data.table(harm)

    methods <- if (nrow(harm) == 1) {
      "mr_wald_ratio"
    } else if (nrow(harm) == 2) {
      c("mr_ivw", "mr_weighted_median")
    } else {
      c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
    }
    res <- tryCatch(mr(harm, method_list = methods), error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    res <- as.data.table(generate_odds_ratios(res))
    res[, protein := unique(exp_raw$protein)]
    res[, F_stat_mean := mean(exp_raw$F_stat, na.rm = TRUE)]
    res[, exposure_sample_size := unique(exp_raw$samplesize.exposure)[1]]
    mr_all[[paste(pid, oc$outcome, sep = "__")]] <- res
  }
}

coverage_dt <- rbindlist(coverage, fill = TRUE)
harm_dt <- rbindlist(harm_all, fill = TRUE)
mr_dt <- rbindlist(mr_all, fill = TRUE)

if (nrow(mr_dt) > 0) {
  mr_dt[, fdr_within_outcome := p.adjust(pval, method = "BH"), by = outcome]
  setcolorder(
    mr_dt,
    intersect(
      c("protein", "id.exposure", "exposure", "id.outcome", "outcome", "method", "nsnp",
        "b", "se", "pval", "or", "or_lci95", "or_uci95", "fdr_within_outcome",
        "F_stat_mean", "exposure_sample_size"),
      names(mr_dt)
    )
  )
}

fwrite(coverage_dt, file.path(out_dir, "opengwas_5protein_replication_snp_coverage.csv"))
fwrite(harm_dt, file.path(out_dir, "opengwas_5protein_replication_harmonised.csv"))
fwrite(mr_dt, file.path(out_dir, "opengwas_5protein_replication_mr_results.csv"))

message("Saved:")
message("  ", file.path(out_dir, "opengwas_5protein_replication_snp_coverage.csv"))
message("  ", file.path(out_dir, "opengwas_5protein_replication_harmonised.csv"))
message("  ", file.path(out_dir, "opengwas_5protein_replication_mr_results.csv"))
message("MR rows: ", nrow(mr_dt))
