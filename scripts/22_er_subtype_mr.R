#!/usr/bin/env Rscript
# Script 18: ER-subtype MR for the 17 FDR<0.05 Phase 2 protein→breast hits
#
# Tests each of the 16 Phase 2 breast-significant proteins against:
#   - ER+  breast cancer (GCST90018758, same paper/format as overall BCAC GCST90018757)
#   - ER-  breast cancer (GCST90018759)
#
# Matching strategy: make_variant_key(chr, pos, allele1, allele2) on the cancer GWAS —
# identical to 04c_harmonise_protein_cancer.R (which got 679/1062 matches for overall breast).
# Both ER files are GRCh37, same column format as Breast_GCST90018757.h.tsv.gz.
#
# Outputs:
#   data/harmonised/harmonised_protein_BreastER_pos_GCST90018758.rds
#   data/harmonised/harmonised_protein_BreastER_neg_GCST90018759.rds
#   results/er_subtype/er_subtype_mr_results.csv
#   results/er_subtype/er_subtype_comparison.csv  — wide: overall | ER+ | ER-

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

project_dir <- "."
harm_dir    <- file.path(project_dir, "data", "harmonised")
out_dir     <- file.path(project_dir, "results", "er_subtype")
dir.create(out_dir, showWarnings = FALSE)

make_variant_key <- function(chr, pos, allele1, allele2) {
  chr     <- sub("^chr", "", as.character(chr))
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  a_min   <- ifelse(allele1 <= allele2, allele1, allele2)
  a_max   <- ifelse(allele1 <= allele2, allele2, allele1)
  paste(chr, pos, a_min, a_max, sep = ":")
}

# ── 16 breast FDR<0.05 proteins (ABO→EC excluded)
breast_hits <- c("SNX15","EFNA1","FGF5","UMOD","SWAP70","ATRAID","TNFRSF6B",
                 "ITIH3","KLB","PM20D1","TSPAN8","FGFR4","IL34","APOE",
                 "CGREF1","INHBB")

# ── Load pQTL instruments — filter to breast hits
cat("Loading pQTL instruments...\n")
pqtls   <- fread(file.path(project_dir, "data", "pqtl", "pqtl_instruments.csv"))
exp_dat <- pqtls[exposure %in% breast_hits & mr_keep == TRUE]
setDF(exp_dat)
cat(sprintf("  %d instruments across %d proteins\n",
            nrow(exp_dat), length(unique(exp_dat$exposure))))
cat("  SNP sample:", paste(head(exp_dat$SNP, 3), collapse=", "), "\n\n")

# ── ER GWAS files
er_files <- list(
  BreastER_pos_GCST90018758 = file.path(project_dir, "data", "cancer_gwas",
                                         "BreastER_pos_GCST90018758.h.tsv.gz"),
  BreastER_neg_GCST90018759 = file.path(project_dir, "data", "cancer_gwas",
                                         "BreastER_neg_GCST90018759.h.tsv.gz")
)

harm_results <- list()  # store harmonised data by subtype

for (label in names(er_files)) {
  gwas_path <- er_files[[label]]
  cat(sprintf("=== Harmonising %s ===\n", label))

  cat("  Loading GWAS...\n")
  out_dat_raw <- fread(
    cmd = sprintf("gunzip -c '%s'", gwas_path),
    select = c("chromosome","base_pair_location","effect_allele","other_allele",
               "beta","standard_error","effect_allele_frequency","p_value")
  )
  cat(sprintf("  Loaded %d variants\n", nrow(out_dat_raw)))

  # Build variant key — same as 04c
  out_dat_raw[, SNP_variant_key := make_variant_key(
    chromosome, base_pair_location, effect_allele, other_allele)]

  # Filter to pQTL positions
  out_dat_raw <- out_dat_raw[SNP_variant_key %in% exp_dat$SNP]
  cat(sprintf("  After variant-key filter: %d rows\n", nrow(out_dat_raw)))

  if (nrow(out_dat_raw) == 0) {
    cat("  No matching variants — skipping\n\n")
    harm_results[[label]] <- data.frame()
    next
  }

  setDF(out_dat_raw)
  out_dat <- format_data(
    out_dat_raw,
    type              = "outcome",
    snp_col           = "SNP_variant_key",
    beta_col          = "beta",
    se_col            = "standard_error",
    eaf_col           = "effect_allele_frequency",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    pval_col          = "p_value"
  )
  out_dat$SNP        <- toupper(out_dat$SNP)
  out_dat$outcome    <- label
  out_dat$id.outcome <- label

  cat(sprintf("  SNP overlap (exp vs out): %d\n",
              length(intersect(exp_dat$SNP, out_dat$SNP))))

  harm_dat <- tryCatch(
    harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2),
    error = function(e) { cat("  harmonise_data ERROR:", conditionMessage(e), "\n"); data.frame() }
  )
  n_keep <- if (nrow(harm_dat) > 0) sum(harm_dat$mr_keep) else 0L
  cat(sprintf("  Harmonised: %d rows (%d mr_keep)\n\n", nrow(harm_dat), n_keep))

  # Save harmonised RDS
  rds_path <- file.path(harm_dir, paste0("harmonised_protein_", label, ".rds"))
  saveRDS(harm_dat, rds_path)
  harm_results[[label]] <- harm_dat
}

# ── MR for each subtype
run_mr <- function(harm_all, label) {
  all_res <- list()
  prots   <- intersect(breast_hits, unique(harm_all$exposure[harm_all$mr_keep]))

  for (prot in prots) {
    harm <- harm_all[harm_all$exposure == prot & harm_all$mr_keep, ]
    if (nrow(harm) == 0) next

    method <- if (nrow(harm) >= 2) "mr_ivw" else "mr_wald_ratio"
    res <- tryCatch(mr(harm, method_list = method), error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next

    res$exposure <- prot
    res$outcome  <- label
    res$or       <- exp(res$b)
    res$or_lci95 <- exp(res$b - 1.96 * res$se)
    res$or_uci95 <- exp(res$b + 1.96 * res$se)
    all_res[[prot]] <- res
  }
  rbindlist(all_res, fill = TRUE)
}

cat("=== Running MR ===\n")
mr_list <- list()
for (label in names(harm_results)) {
  harm_all <- harm_results[[label]]
  if (is.null(harm_all) || nrow(harm_all) == 0) next
  cat(sprintf("\n--- %s ---\n", label))
  res <- run_mr(harm_all, label)
  if (nrow(res) > 0) {
    mr_list[[label]] <- res
    print(res[, c("exposure","method","nsnp","b","se","pval","or")])
  }
}

all_mr <- rbindlist(mr_list, fill = TRUE)
fwrite(all_mr, file.path(out_dir, "er_subtype_mr_results.csv"))
cat(sprintf("\nSaved %d MR result rows\n", nrow(all_mr)))

# ── Load overall breast results
overall <- fread(file.path(project_dir, "results", "phase2_protein_cancer",
                            "protein_cancer_mr_results_full.csv"))
overall_b <- overall[outcome == "Breast_GCST90018757" & exposure %in% breast_hits,
                     .(exposure, b_overall=round(b,4), se_overall=round(se,5),
                       pval_overall=signif(pval,3), fdr_overall=signif(fdr,3),
                       or_overall=round(or,4), nsnp_overall=nsnp)]

# ── Wide comparison table
er_pos <- all_mr[outcome == "BreastER_pos_GCST90018758",
                 .(exposure, b_ERpos=round(b,4), se_ERpos=round(se,5),
                   pval_ERpos=signif(pval,3), or_ERpos=round(or,4), nsnp_ERpos=nsnp)]
er_neg <- all_mr[outcome == "BreastER_neg_GCST90018759",
                 .(exposure, b_ERneg=round(b,4), se_ERneg=round(se,5),
                   pval_ERneg=signif(pval,3), or_ERneg=round(or,4), nsnp_ERneg=nsnp)]

comp <- merge(overall_b, er_pos, by="exposure", all.x=TRUE)
comp <- merge(comp, er_neg,      by="exposure", all.x=TRUE)
comp[, consistent_dir := fifelse(
  !is.na(b_ERpos) & !is.na(b_ERneg),
  sign(b_overall)==sign(b_ERpos) & sign(b_overall)==sign(b_ERneg), NA)]
comp[, er_pattern := fifelse(
  !is.na(pval_ERpos) & !is.na(pval_ERneg),
  fifelse(pval_ERpos < 0.05 & pval_ERneg >= 0.05, "ER_pos_specific",
  fifelse(pval_ERneg < 0.05 & pval_ERpos >= 0.05, "ER_neg_specific",
  fifelse(pval_ERpos < 0.05 & pval_ERneg < 0.05,  "both_subtypes",  "neither"))),
  NA_character_)]
comp <- comp[order(pval_overall)]
fwrite(comp, file.path(out_dir, "er_subtype_comparison.csv"))

cat("\n=== ER SUBTYPE COMPARISON ===\n")
print(comp[, .(exposure, or_overall, pval_overall,
               or_ERpos, pval_ERpos,
               or_ERneg, pval_ERneg,
               consistent_dir, er_pattern)])

cat("\nOutputs:\n")
cat("  ", file.path(out_dir, "er_subtype_mr_results.csv"), "\n")
cat("  ", file.path(out_dir, "er_subtype_comparison.csv"), "\n")
cat("Done.\n")
sessionInfo()
