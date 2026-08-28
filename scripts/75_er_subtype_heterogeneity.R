#!/usr/bin/env Rscript
## Script 75: Formal ER+ versus ER- subtype heterogeneity tests
## Output:
##   results/tables/STable_ER_subtype_heterogeneity.csv
##   results/er_subtype/er_subtype_comparison_with_heterogeneity.csv

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- "."
er_dir <- file.path(project_dir, "results", "er_subtype")
out_dir <- file.path(project_dir, "results", "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

mr_file <- file.path(er_dir, "er_subtype_mr_results.csv")
comp_file <- file.path(er_dir, "er_subtype_comparison.csv")

if (!file.exists(mr_file)) {
  stop("Missing ER-subtype MR results: ", mr_file)
}
if (!file.exists(comp_file)) {
  stop("Missing ER-subtype comparison table: ", comp_file)
}

mr <- fread(mr_file)
comp <- fread(comp_file)

er_pos <- mr[
  outcome == "BreastER_pos_GCST90018758",
  .(
    exposure,
    b_ERpos_raw = b,
    se_ERpos_raw = se,
    pval_ERpos_raw = pval,
    or_ERpos_raw = or,
    nsnp_ERpos_raw = nsnp,
    method_ERpos = method
  )
]

er_neg <- mr[
  outcome == "BreastER_neg_GCST90018759",
  .(
    exposure,
    b_ERneg_raw = b,
    se_ERneg_raw = se,
    pval_ERneg_raw = pval,
    or_ERneg_raw = or,
    nsnp_ERneg_raw = nsnp,
    method_ERneg = method
  )
]

het <- merge(er_pos, er_neg, by = "exposure", all = FALSE)
het[, beta_diff_ERpos_minus_ERneg := b_ERpos_raw - b_ERneg_raw]
het[, se_diff := sqrt(se_ERpos_raw^2 + se_ERneg_raw^2)]
het[, z_heterogeneity := beta_diff_ERpos_minus_ERneg / se_diff]
het[, p_heterogeneity := 2 * pnorm(abs(z_heterogeneity), lower.tail = FALSE)]
het[, fdr_heterogeneity := p.adjust(p_heterogeneity, method = "BH")]
het[, stronger_subtype := fifelse(
  abs(b_ERpos_raw) > abs(b_ERneg_raw),
  "ER-positive stronger",
  "ER-negative stronger"
)]
## "Subtype-shared association" must require that the ER+ and ER- estimates agree
## in direction with each other AND with the overall estimate (comp$consistent_dir).
## Fixed 2026-08-28 after external peer review: nominal p<0.05 in both subtypes was
## previously sufficient, which mislabeled TNFRSF6B as "shared" even though both
## subtype estimates were protective while the overall MR estimate was risk-increasing
## (consistent_dir = FALSE for that row). Proteins that are nominally significant in
## both subtypes but directionally discordant with the overall estimate, and do not
## reach FDR-significant heterogeneity either, get their own explicit label instead of
## being folded into "shared" or silently defaulting to "No formal heterogeneity".
comp_dir <- comp[, .(exposure, consistent_dir)]
het <- merge(het, comp_dir, by = "exposure", all.x = TRUE)
het[, subtype_class := fcase(
  fdr_heterogeneity < 0.05 & stronger_subtype == "ER-negative stronger",
    "ER-negative-enriched",
  fdr_heterogeneity < 0.05 & stronger_subtype == "ER-positive stronger",
    "ER-positive-enriched",
  pval_ERpos_raw < 0.05 & pval_ERneg_raw < 0.05 & consistent_dir == TRUE,
    "Subtype-shared association",
  pval_ERpos_raw < 0.05 & pval_ERneg_raw < 0.05 & consistent_dir == FALSE,
    "Directionally discordant with overall estimate",
  default = "No formal heterogeneity"
)]

out <- merge(
  comp,
  het[, .(
    exposure,
    b_ERpos_raw,
    se_ERpos_raw,
    pval_ERpos_raw,
    or_ERpos_raw,
    nsnp_ERpos_raw,
    method_ERpos,
    b_ERneg_raw,
    se_ERneg_raw,
    pval_ERneg_raw,
    or_ERneg_raw,
    nsnp_ERneg_raw,
    method_ERneg,
    beta_diff_ERpos_minus_ERneg,
    se_diff,
    z_heterogeneity,
    p_heterogeneity,
    fdr_heterogeneity,
    stronger_subtype,
    subtype_class
  )],
  by = "exposure",
  all.x = TRUE
)

setorder(out, p_heterogeneity)

pretty <- out[, .(
  protein = exposure,
  b_overall,
  se_overall,
  pval_overall,
  b_ERpos = round(b_ERpos_raw, 5),
  se_ERpos = round(se_ERpos_raw, 5),
  pval_ERpos = signif(pval_ERpos_raw, 3),
  b_ERneg = round(b_ERneg_raw, 5),
  se_ERneg = round(se_ERneg_raw, 5),
  pval_ERneg = signif(pval_ERneg_raw, 3),
  beta_diff_ERpos_minus_ERneg = round(beta_diff_ERpos_minus_ERneg, 5),
  se_diff = round(se_diff, 5),
  z_heterogeneity = round(z_heterogeneity, 3),
  p_heterogeneity = signif(p_heterogeneity, 3),
  fdr_heterogeneity = signif(fdr_heterogeneity, 3),
  stronger_subtype,
  subtype_class,
  er_pattern_nominal = er_pattern,
  consistent_dir
)]

table_file <- file.path(out_dir, "STable_ER_subtype_heterogeneity.csv")
aug_file <- file.path(er_dir, "er_subtype_comparison_with_heterogeneity.csv")

fwrite(pretty, table_file)
fwrite(out, aug_file)

cat("Saved ER subtype heterogeneity tables:\n")
cat("  ", table_file, "\n", sep = "")
cat("  ", aug_file, "\n\n", sep = "")

cat("Top heterogeneity results:\n")
print(pretty[, .(
  protein,
  beta_diff_ERpos_minus_ERneg,
  p_heterogeneity,
  fdr_heterogeneity,
  stronger_subtype,
  subtype_class
)])
