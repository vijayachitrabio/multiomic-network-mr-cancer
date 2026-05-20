#!/usr/bin/env Rscript

# Script 32: Trial selective protein + metabolite multivariable MR
#
# Pilot MVMR for the most interpretable breast-cancer pairs that are feasible
# with the current local data and public FinnGen pQTL files.

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(Rsamtools)
})

project_dir <- normalizePath(".")
out_dir <- file.path(project_dir, "results", "mvmr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pairs <- data.table(
  protein = c("PM20D1", "IL34"),
  metabolite = c("Gly", "Total_BCAA"),
  outcome = "Breast_GCST90018757"
)
TOP_MQTL_INSTRUMENTS <- 5L

make_variant_key <- function(chr, pos, allele1, allele2) {
  chr <- sub("^chr", "", as.character(chr))
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  a_min <- ifelse(allele1 <= allele2, allele1, allele2)
  a_max <- ifelse(allele1 <= allele2, allele2, allele1)
  paste(chr, pos, a_min, a_max, sep = ":")
}

query_local_gwas_by_pos <- function(path, keys, type = c("breast", "metabolite")) {
  type <- match.arg(type)
  if (length(keys) == 0) {
    return(data.table(
      SNP = character(), chr = integer(), pos = integer(),
      effect_allele = character(), other_allele = character(),
      eaf = numeric(), beta = numeric(), se = numeric(), pval = numeric()
    ))
  }

  key_string <- paste(unique(keys), collapse = ",")
  cmd <- sprintf(
    "gunzip -c %s | awk 'BEGIN{FS=OFS=\"\\t\"; n=split(\"%s\", arr, \",\"); for(i=1;i<=n;i++) keep[arr[i]]=1} NR==1 || (($1 \":\" $2) in keep)'",
    shQuote(path), key_string
  )
  dt <- fread(cmd = cmd, showProgress = FALSE)
  if (nrow(dt) == 0) {
    return(data.table(
      SNP = character(), chr = integer(), pos = integer(),
      effect_allele = character(), other_allele = character(),
      eaf = numeric(), beta = numeric(), se = numeric(), pval = numeric()
    ))
  }

  if (type == "breast") {
    setnames(
      dt,
      c("chromosome", "base_pair_location", "effect_allele", "other_allele",
        "effect_allele_frequency", "beta", "standard_error", "p_value"),
      c("chr", "pos", "effect_allele", "other_allele",
        "eaf", "beta", "se", "pval"),
      skip_absent = TRUE
    )
  } else {
    setnames(
      dt,
      c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
        "N", "BETA", "SE", "LOG10P"),
      c("chr", "pos", "rsid", "other_allele", "effect_allele", "eaf",
        "n", "beta", "se", "log10p"),
      skip_absent = TRUE
    )
    dt[, pval := 10^(-pmin(as.numeric(log10p), 300))]
  }

  dt[, `:=`(
    chr = suppressWarnings(as.integer(chr)),
    pos = suppressWarnings(as.integer(pos)),
    beta = suppressWarnings(as.numeric(beta)),
    se = suppressWarnings(as.numeric(se)),
    eaf = suppressWarnings(as.numeric(eaf)),
    pval = suppressWarnings(as.numeric(pval))
  )]
  dt <- dt[!is.na(chr) & !is.na(pos)]
  dt[, SNP := make_variant_key(chr, pos, effect_allele, other_allele)]
  unique(dt, by = "SNP")
}

query_remote_pqtl_for_pos <- function(protein, pos_dt) {
  if (nrow(pos_dt) == 0) {
    return(data.table(
      SNP = character(), chr = integer(), pos = integer(),
      effect_allele = character(), other_allele = character(),
      eaf = numeric(), beta = numeric(), se = numeric(), pval = numeric(), n = numeric()
    ))
  }

  base_url <- "https://storage.googleapis.com/finngen-public-data-r10/omics/proteomics/release_2023_03_02/data/Olink/pQTL"
  url <- sprintf("%s/Olink_Batch1_%s.txt.gz", base_url, protein)
  tbx <- TabixFile(url)

  out <- vector("list", nrow(pos_dt))
  for (i in seq_len(nrow(pos_dt))) {
    if (i %% 5 == 0 || i == nrow(pos_dt)) {
      message("  remote pQTL lookups for ", protein, ": ", i, "/", nrow(pos_dt))
    }
    region <- GRanges(as.character(pos_dt$chr[i]), IRanges(pos_dt$pos[i], pos_dt$pos[i]))
    lines <- tryCatch(scanTabix(tbx, param = region)[[1]], error = function(e) character())
    if (length(lines) == 0) next
    x <- fread(text = paste(lines, collapse = "\n"), header = FALSE, showProgress = FALSE)
    if (ncol(x) == 11) {
      setnames(x, c("chr", "pos", "variant_id", "ref", "alt", "alt_freq",
                    "beta", "se", "t_stat", "pval", "log10p"))
      x[, n := 619L]
    } else {
      setnames(x, c("chr", "pos", "variant_id", "ref", "alt", "alt_freq",
                    "beta", "se", "t_stat", "pval", "log10p", "n"))
    }
    x[, `:=`(
      chr = suppressWarnings(as.integer(chr)),
      pos = suppressWarnings(as.integer(pos)),
      beta = suppressWarnings(as.numeric(beta)),
      se = suppressWarnings(as.numeric(se)),
      pval = suppressWarnings(as.numeric(pval)),
      eaf = suppressWarnings(as.numeric(alt_freq)),
      effect_allele = toupper(as.character(alt)),
      other_allele = toupper(as.character(ref))
    )]
    x[, SNP := make_variant_key(chr, pos, effect_allele, other_allele)]
    out[[i]] <- x[, .(SNP, chr, pos, effect_allele, other_allele, eaf, beta, se, pval, n)]
  }
  ans <- rbindlist(out, fill = TRUE)
  if (nrow(ans) == 0) {
    return(data.table(
      SNP = character(), chr = integer(), pos = integer(),
      effect_allele = character(), other_allele = character(),
      eaf = numeric(), beta = numeric(), se = numeric(), pval = numeric(), n = numeric()
    ))
  }
  ans
}

format_exposure <- function(dt, exposure_name, exposure_id) {
  dt[, .(
    SNP,
    beta.exposure = beta,
    se.exposure = se,
    pval.exposure = pval,
    effect_allele.exposure = effect_allele,
    other_allele.exposure = other_allele,
    eaf.exposure = eaf,
    exposure = exposure_name,
    id.exposure = exposure_id
  )]
}

format_outcome <- function(dt, outcome_name, outcome_id) {
  dt[, .(
    SNP,
    beta.outcome = beta,
    se.outcome = se,
    pval.outcome = pval,
    effect_allele.outcome = effect_allele,
    other_allele.outcome = other_allele,
    eaf.outcome = eaf,
    outcome = outcome_name,
    id.outcome = outcome_id
  )]
}

pqtl <- fread(file.path(project_dir, "data", "pqtl", "pqtl_instruments_finngen_olink.csv"))
mqtl_inst <- fread(file.path(project_dir, "data", "mqtl", "mqtl_gwas_instruments.csv"))
uv_protein <- fread(file.path(project_dir, "results", "phase2_protein_cancer", "protein_cancer_mr_results_full.csv"))
uv_pm <- fread(file.path(project_dir, "results", "phase3_protein_metabolite", "protein_metabolite_mr_results_full.csv"))
uv_mc <- fread(file.path(project_dir, "results", "phase4_metabolite_cancer", "metabolite_cancer_mr_results_full.csv"))

breast_path <- file.path(project_dir, "data", "cancer_gwas", "Breast_GCST90018757.h.tsv.gz")

all_results <- list()
all_qc <- list()

for (i in seq_len(nrow(pairs))) {
  prot_name <- pairs$protein[i]
  met_name <- pairs$metabolite[i]
  pair_id <- paste(prot_name, met_name, "Breast", sep = "__")
  message("Running MVMR trial for ", prot_name, " + ", met_name, " -> Breast")

  p_inst <- pqtl[exposure == prot_name & mr_keep == TRUE,
    .(SNP, chr = chr.exposure, pos = pos.exposure,
      effect_allele = effect_allele.exposure, other_allele = other_allele.exposure,
      eaf = eaf.exposure, beta = beta.exposure, se = se.exposure,
      pval = pval.exposure)]
  m_inst <- mqtl_inst[metabolite == met_name,
    .(rsid = SNP, chr, pos, effect_allele, other_allele,
      eaf, beta, se, pval, n)]
  setorder(m_inst, pval)
  m_inst <- m_inst[seq_len(min(nrow(m_inst), TOP_MQTL_INSTRUMENTS))]
  m_inst[, SNP := make_variant_key(chr, pos, effect_allele, other_allele)]

  union_pos <- unique(rbindlist(list(
    p_inst[, .(chr, pos)],
    m_inst[, .(chr, pos)]
  )))
  union_keys <- unique(c(p_inst$SNP, m_inst$SNP))
  message("  union SNP positions: ", nrow(union_pos))

  # Metabolite exposure: start with instrument rows, add protein-instrument positions
  m_extra <- query_local_gwas_by_pos(
    file.path(project_dir, "data", "mqtl", "mqtl_full_gwas", paste0(met_name, "_full_regenie.tsv.gz")),
    keys = unique(paste(p_inst$chr, p_inst$pos, sep = ":")),
    type = "metabolite"
  )
  m_assoc <- rbindlist(list(
    m_inst[, .(SNP, chr, pos, effect_allele, other_allele, eaf, beta, se, pval)],
    m_extra[, .(SNP, chr, pos, effect_allele, other_allele, eaf, beta, se, pval)]
  ), fill = TRUE)
  m_assoc <- unique(m_assoc, by = "SNP")
  m_assoc <- m_assoc[SNP %in% union_keys]
  message("  metabolite associations assembled: ", uniqueN(m_assoc$SNP))

  # Protein exposure: use local instruments plus remote lookups for all union SNP positions
  p_remote <- query_remote_pqtl_for_pos(prot_name, union_pos)
  p_assoc <- rbindlist(list(p_inst, p_remote), fill = TRUE)
  p_assoc <- unique(p_assoc, by = "SNP")
  p_assoc <- p_assoc[SNP %in% union_keys]
  message("  protein associations assembled: ", uniqueN(p_assoc$SNP))

  # Outcome: breast associations at union SNP positions
  outcome_raw <- query_local_gwas_by_pos(
    breast_path,
    keys = unique(paste(union_pos$chr, union_pos$pos, sep = ":")),
    type = "breast"
  )
  outcome_assoc <- outcome_raw[SNP %in% union_keys, .(SNP, chr, pos, effect_allele, other_allele, eaf, beta, se, pval)]
  outcome_assoc <- unique(outcome_assoc, by = "SNP")
  message("  breast outcome associations assembled: ", uniqueN(outcome_assoc$SNP))

  exposure_dat <- rbindlist(list(
    format_exposure(p_assoc, prot_name, paste0("pQTL_", prot_name)),
    format_exposure(m_assoc, met_name, paste0("mQTL_", met_name))
  ), fill = TRUE)
  outcome_dat <- format_outcome(outcome_assoc, "Breast_cancer", "Breast_GCST90018757")

  common_snps <- Reduce(intersect, list(
    unique(p_assoc$SNP), unique(m_assoc$SNP), unique(outcome_assoc$SNP)
  ))

  qc_row <- data.table(
    pair = pair_id,
    protein = prot_name,
    metabolite = met_name,
    n_protein_inst = nrow(p_inst),
    n_metabolite_inst = nrow(m_inst),
    n_union_snps = length(union_keys),
    n_protein_assoc = uniqueN(p_assoc$SNP),
    n_metabolite_assoc = uniqueN(m_assoc$SNP),
    n_outcome_assoc = uniqueN(outcome_assoc$SNP),
    n_common_snps = length(common_snps)
  )

  if (length(common_snps) < 3) {
    qc_row[, note := "fewer_than_3_common_snps"]
    all_qc[[pair_id]] <- qc_row
    next
  }

  exposure_dat <- exposure_dat[SNP %in% common_snps]
  outcome_dat <- outcome_dat[SNP %in% common_snps]

  mvdat <- mv_harmonise_data(exposure_dat, outcome_dat, harmonise_strictness = 2)
  mvres <- mv_multiple(mvdat, intercept = FALSE, instrument_specific = FALSE)

  res <- as.data.table(mvres$result)
  res[, `:=`(
    pair = pair_id,
    protein = prot_name,
    metabolite = met_name,
    n_common_snps = length(common_snps)
  )]
  all_results[[pair_id]] <- res

  qc_row[, `:=`(
    note = "ok",
    n_harmonised_snps = nrow(mvdat$exposure_beta)
  )]
  all_qc[[pair_id]] <- qc_row
}

qc_dt <- rbindlist(all_qc, fill = TRUE)
res_dt <- rbindlist(all_results, fill = TRUE)

if (nrow(res_dt) > 0) {
  uni_summary <- rbindlist(list(
    uv_protein[exposure %in% pairs$protein & outcome == "Breast_GCST90018757",
      .(type = "protein_to_breast_uv", exposure, outcome, method, nsnp, b, se, pval, fdr)],
    uv_pm[exposure %in% pairs$protein & outcome %in% pairs$metabolite,
      .(type = "protein_to_metabolite_uv", exposure, outcome, method, nsnp, b, se, pval, fdr)],
    uv_mc[exposure %in% pairs$metabolite & outcome == "Breast_GCST90018757",
      .(type = "metabolite_to_breast_uv", exposure, outcome, method, nsnp, b, se, pval, fdr)]
  ), fill = TRUE)
  fwrite(res_dt, file.path(out_dir, "selective_mvmr_trial_results.csv"))
  fwrite(uni_summary, file.path(out_dir, "selective_mvmr_univariable_context.csv"))
}

fwrite(qc_dt, file.path(out_dir, "selective_mvmr_trial_qc.csv"))

summary_lines <- c(
  "# Selective MVMR Trial",
  "",
  "Pairs tested:",
  "- PM20D1 + Gly -> Breast cancer",
  "- IL34 + Total_BCAA -> Breast cancer",
  "",
  "This pilot uses local metabolite GWAS + breast GWAS plus public FinnGen Olink",
  "pQTL summary statistics queried by genomic position."
)
writeLines(summary_lines, file.path(out_dir, "SELECTIVE_MVMR_TRIAL_NOTE.md"))

message("Selective MVMR trial complete.")
