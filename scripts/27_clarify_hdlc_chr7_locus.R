#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"

mqtl_path <- file.path(project_dir, "data", "mqtl", "mqtl_full_gwas", "HDL_C_full_regenie.tsv.gz")
cancer_path <- file.path(project_dir, "data", "cancer_gwas", "Breast_GCST90018757.h.tsv.gz")
best_loci_path <- file.path(project_dir, "results", "validation", "metabolite_cancer_coloc_best_loci.csv")
out_dir <- file.path(project_dir, "results", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

best <- fread(best_loci_path)
hdl_best <- best[exposure == "HDL_C" & outcome == "Breast_GCST90018757"]
if (nrow(hdl_best) != 1) stop("Could not identify unique HDL_C best locus row.")

lead_snp <- hdl_best$lead_snp[[1]]
lead_chr <- hdl_best$chr[[1]]
lead_pos <- hdl_best$lead_pos[[1]]
window_bp <- 500000L

message("Reading HDL_C summary statistics in lead window...")
mqtl <- fread(
  cmd = sprintf("gunzip -c '%s'", mqtl_path),
  select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P")
)
setnames(mqtl, c("chr", "pos_met", "SNP", "oa_met", "ea_met", "eaf_met", "beta_met", "se_met", "log10p_met"))
mqtl <- mqtl[chr == lead_chr & pos_met >= (lead_pos - window_bp) & pos_met <= (lead_pos + window_bp)]

message("Reading breast cancer summary statistics...")
cancer <- fread(
  cmd = sprintf("gunzip -c '%s'", cancer_path),
  select = c("chromosome", "base_pair_location", "rsid", "effect_allele", "other_allele", "beta", "standard_error", "p_value")
)
setnames(cancer, c("chr_cancer", "pos_cancer", "SNP", "ea_cancer", "oa_cancer", "beta_cancer", "se_cancer", "p_cancer"))
cancer[, `:=`(
  chr_cancer = as.integer(chr_cancer),
  pos_cancer = as.integer(pos_cancer),
  beta_cancer = as.numeric(beta_cancer),
  se_cancer = as.numeric(se_cancer),
  p_cancer = as.numeric(p_cancer)
)]
setkey(mqtl, SNP)
setkey(cancer, SNP)

overlap <- cancer[mqtl, nomatch = 0]
overlap <- overlap[chr_cancer == lead_chr]
overlap[, cancer_minus_log10p := -log10(p_cancer)]
overlap[, abs_pos_delta := abs(pos_cancer - pos_met)]

top_by_cancer <- overlap[order(p_cancer, -log10p_met)][1:25,
  .(SNP, pos_met, pos_cancer, abs_pos_delta, p_cancer, cancer_minus_log10p, log10p_met, beta_met, beta_cancer)]

top_by_hdl <- overlap[order(-log10p_met, p_cancer)][1:25,
  .(SNP, pos_met, pos_cancer, abs_pos_delta, p_cancer, cancer_minus_log10p, log10p_met, beta_met, beta_cancer)]

lead_row <- overlap[SNP == lead_snp,
  .(SNP, pos_met, pos_cancer, abs_pos_delta, p_cancer, cancer_minus_log10p, log10p_met, beta_met, beta_cancer)]

summary_dt <- data.table(
  lead_snp = lead_snp,
  lead_pos_met = lead_pos,
  n_window_snps_mqtl = nrow(mqtl),
  n_overlap_snps = nrow(overlap),
  lead_snp_breast_p = if (nrow(lead_row)) lead_row$p_cancer[[1]] else NA_real_,
  lead_snp_hdl_log10p = if (nrow(lead_row)) lead_row$log10p_met[[1]] else NA_real_,
  top_overlap_snp_by_breast_p = top_by_cancer$SNP[[1]],
  top_overlap_snp_by_breast_pos_met = top_by_cancer$pos_met[[1]],
  top_overlap_snp_by_breast_pos_cancer = top_by_cancer$pos_cancer[[1]],
  top_overlap_snp_by_breast_p = top_by_cancer$p_cancer[[1]],
  top_overlap_snp_by_breast_hdl_log10p = top_by_cancer$log10p_met[[1]],
  top_overlap_snp_by_hdl = top_by_hdl$SNP[[1]],
  top_overlap_snp_by_hdl_pos_met = top_by_hdl$pos_met[[1]],
  top_overlap_snp_by_hdl_pos_cancer = top_by_hdl$pos_cancer[[1]],
  top_overlap_snp_by_hdl_p = top_by_hdl$p_cancer[[1]],
  top_overlap_snp_by_hdl_log10p = top_by_hdl$log10p_met[[1]]
)

fwrite(summary_dt, file.path(out_dir, "hdl_c_chr7_locus_summary.csv"))
fwrite(top_by_cancer, file.path(out_dir, "hdl_c_chr7_top_overlap_by_breast.csv"))
fwrite(top_by_hdl, file.path(out_dir, "hdl_c_chr7_top_overlap_by_hdl.csv"))
if (nrow(lead_row)) fwrite(lead_row, file.path(out_dir, "hdl_c_chr7_lead_snp_detail.csv"))

note <- c(
  "# HDL_C chr7 Locus Clarification",
  "",
  sprintf("Lead coloc window anchor from metabolite-cancer analysis: `%s` at metabolite position chr7:%d.", lead_snp, lead_pos),
  "",
  "## Main clarification",
  "",
  "The strong HDL_C -> Breast colocalization signal (PP.H4 = 0.967) is not driven by the anchor SNP `rs62463430` itself.",
  "",
  sprintf("- `%s` has strong HDL_C association in the metabolite GWAS (log10P = %.3f) but only weak breast cancer association (p = %.3g).",
          lead_snp, summary_dt$lead_snp_hdl_log10p[[1]], summary_dt$lead_snp_breast_p[[1]]),
  sprintf("- The strongest overlapping breast-cancer SNPs in the exact colocalization SNP set are `%s`, `%s`, `%s`, `%s`, and `%s`, all in the 72.82-72.86 Mb metabolite-position block and all with breast cancer p-values around 10^-7.",
          top_by_cancer$SNP[[1]], top_by_cancer$SNP[[2]], top_by_cancer$SNP[[3]], top_by_cancer$SNP[[4]], top_by_cancer$SNP[[5]]),
  sprintf("- The strongest HDL_C SNPs in the same overlapping set are `%s`, `%s`, `%s`, `%s`, and `%s`; these are also in the same 72.85-72.90 Mb block and show breast cancer p-values around 10^-7 to 10^-6.",
          top_by_hdl$SNP[[1]], top_by_hdl$SNP[[2]], top_by_hdl$SNP[[3]], top_by_hdl$SNP[[4]], top_by_hdl$SNP[[5]]),
  "",
  "## Coordinate note",
  "",
  "The lead rsID `rs62463430` has different recorded positions between the HDL_C file and the breast cancer GWAS file. The colocalization script matches by rsID and then anchors the tested window on the metabolite position. Therefore, `rs62463430` should be treated as the window anchor rather than the mechanistic driver of the shared signal.",
  "",
  "## Interpretation",
  "",
  "The chr7 HDL_C colocalization appears to reflect a broader shared haplotype block centered nearer 72.85-72.90 Mb, not a single weak breast-cancer signal at the anchor SNP itself. This supports keeping the HDL_C -> Breast colocalization result, but the manuscript should describe `rs62463430` as the lead window anchor and the 72.85-72.90 Mb block as the likely shared signal region.",
  "",
  "## Candidate regional biology",
  "",
  "This block lies just upstream of the GRCh37 MLXIPL / ChREBP locus (NCBI Gene: chr7:73,007,532-73,038,852 on GRCh37), which is compatible with a lipid-metabolism interpretation of the HDL_C signal.",
  "",
  "## Output files",
  "",
  "- `results/validation/hdl_c_chr7_locus_summary.csv`",
  "- `results/validation/hdl_c_chr7_top_overlap_by_breast.csv`",
  "- `results/validation/hdl_c_chr7_top_overlap_by_hdl.csv`",
  "- `results/validation/hdl_c_chr7_lead_snp_detail.csv`"
)

writeLines(note, file.path(out_dir, "HDL_C_chr7_locus_clarification.md"))

message("Wrote HDL_C chr7 clarification outputs to results/validation/")
