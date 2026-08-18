suppressMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

project_dir <- "."
decode_dir <- "../SUMMARY DATA-proteins"
out_dir <- file.path(project_dir, "results", "decode")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Validation Targets
targets <- list(
  # list(prot="EFNA1", outcome="Breast_GCST90018757", snps=c("rs4971066"), file="*EFNA1*.txt.gz"), # MISSING in deCODE
  # list(prot="ATRAID", outcome="Breast_GCST90018757", snps=c("rs6727617"), file="*ATRAID*.txt.gz"), # MISSING in deCODE
  list(prot="TNFRSF6B", outcome="Breast_GCST90018757", snps=c("rs6011040"), file="*TNFRSF6B*.txt.gz"),
  # list(prot="SNX15", outcome="Breast_GCST90018757", snps=c("rs12905762"), file="*SNX15*.txt.gz"), # MISSING in deCODE
  list(prot="UMOD", outcome="Breast_GCST90018757", snps=c("rs12917707"), file="*UMOD*.txt.gz"),
  list(prot="IL34", outcome="Breast_GCST90018757", snps=c("rs4985556", "rs4985558"), file="*IL34*.txt.gz"),
  # list(prot="PM20D1", outcome="Breast_GCST90018757", snps=c("rs823116", "rs823118"), file="*PM20D1*.txt.gz"), # MISSING in deCODE
  list(prot="CGREF1", outcome="Breast_GCST90018757", snps=c("rs28399512"), file="*CGREF1*.txt.gz"),
  list(prot="ITIH3", outcome="Breast_GCST90018757", snps=c("rs2535629"), file="*ITIH3*.txt.gz"),
  list(prot="SWAP70", outcome="Breast_GCST90018757", snps=c("rs415895"), file="*SWAP70*.txt.gz"),
  list(prot="INHBB", outcome="Breast_GCST90018757", snps=c("rs1982566"), file="*INHBB*.txt.gz")
)

results_list <- list()

for (t in targets) {
  cat("\n========================================\n")
  cat(sprintf("Processing %s -> %s\n", t$prot, t$outcome))
  
  # Find the decode file in either the iCloud folder or the local decode_pqtl folder
  decode_file <- Sys.glob(file.path(decode_dir, t$file))
  if (length(decode_file) == 0) {
    decode_file <- Sys.glob(file.path(project_dir, "data", "decode_pqtl", paste0("*", t$prot, "*.txt.gz")))
  }
  
  if (length(decode_file) == 0) {
    cat(sprintf("ERROR: Could not find deCODE file for %s\n", t$prot))
    next
  }
  decode_file <- decode_file[1] # Take first match if multiple
  
  cat(sprintf("Found deCODE file: %s\n", basename(decode_file)))
  
  # Outcome file
  outcome_file <- file.path(project_dir, "data", "cancer_gwas", paste0(t$outcome, ".h.tsv.gz"))
  if (!file.exists(outcome_file)) {
    cat(sprintf("ERROR: Outcome file not found: %s\n", outcome_file))
    next
  }
  
  # 1. Extract specifically target SNPs from deCODE file
  cat("Extracting specific target SNPs from deCODE file...\n")
  # Use awk to filter only the exact target SNPs to save memory
  snp_list <- paste(t$snps, collapse = "|")
  cmd <- sprintf("gzcat '%s' | awk 'NR==1 || $4 ~ /^(%s)$/ || $3 ~ /^(%s)$/'", decode_file, snp_list, snp_list)
  
  pqtl_raw <- fread(cmd = cmd)
  
  if (nrow(pqtl_raw) == 0) {
    cat(sprintf("WARNING: None of the target SNPs (%s) found in deCODE file for %s\n", paste(t$snps, collapse=", "), t$prot))
    next
  }
  
  pqtl_raw[, SNP := ifelse(rsids != "" & !is.na(rsids), rsids, Name)]
  
  exp_dat <- format_data(
    as.data.frame(pqtl_raw),
    type = "exposure",
    snp_col = "SNP",
    beta_col = "Beta",
    se_col = "SE",
    effect_allele_col = "effectAllele",
    other_allele_col = "otherAllele",
    eaf_col = "ImpMAF",
    pval_col = "Pval",
    chr_col = "Chrom",
    pos_col = "Pos"
  )
  exp_dat$exposure <- paste0(t$prot, " (deCODE)")
  
  # 2. Match with outcome data
  cat(sprintf("Extracting from Outcome data (%s)...\n", t$outcome))
  out_raw <- fread(outcome_file)
  
  # If Endometrial, use hm_rsid, else rsid
  rsid_col <- ifelse(grepl("Endometrial", t$outcome), "hm_rsid", "rsid")
  
  shared_snps <- intersect(exp_dat$SNP, out_raw[[rsid_col]])
  exp_dat <- exp_dat[exp_dat$SNP %in% shared_snps, ]
  
  if(nrow(exp_dat) == 0) {
    cat("WARNING: Target SNPs not present in outcome data.\n")
    next
  }
  
  out_subset <- out_raw[out_raw[[rsid_col]] %in% exp_dat$SNP, ]
  out_dat <- format_data(
    as.data.frame(out_subset),
    type = "outcome",
    snp_col = rsid_col,
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
  )
  out_dat$outcome <- gsub("_GCST.*", "", t$outcome)
  
  # 3. Harmonize and MR
  cat("Harmonizing data...\n")
  harm_dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2)
  
  if (nrow(harm_dat) == 0) {
    cat("WARNING: No SNPs remained after harmonization.\n")
    next
  }
  
  cat("Running MR...\n")
  mr_res <- mr(harm_dat, method_list = c("mr_wald_ratio", "mr_ivw"))
  mr_res <- generate_odds_ratios(mr_res)
  
  print(mr_res)
  results_list[[t$prot]] <- mr_res
  
  # 4. Cleanup memory
  gc()
}

if (length(results_list) > 0) {
  final_results <- bind_rows(results_list)
  write.csv(final_results, file.path(out_dir, "decode_batch_mr_results.csv"), row.names = FALSE)
  cat("\nBatch MR results saved to", file.path(out_dir, "decode_batch_mr_results.csv"), "\n")
} else {
  cat("\nNo results generated.\n")
}
