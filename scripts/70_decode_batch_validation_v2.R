suppressMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

project_dir <- "."
decode_dir <- "../SUMMARY DATA-proteins"
out_dir <- file.path(project_dir, "results", "decode")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

targets <- list(
  list(prot="FGF5", outcome="Breast_GCST90018757", snps=c("4:80263187"), file="*FGF5*.txt.gz"),
  list(prot="UMOD", outcome="Breast_GCST90018757", snps=c("16:20381010", "16:20344600"), file="*UMOD*.txt.gz"),
  list(prot="APOE", outcome="Breast_GCST90018757", snps=c("19:44908684"), file="*APOE*.txt.gz"),
  list(prot="KLB", outcome="Breast_GCST90018757", snps=c("4:39448609", "4:39450812"), file="*KLB*.txt.gz"),
  list(prot="FGFR4", outcome="Breast_GCST90018757", snps=c("5:177093242"), file="*FGFR4*.txt.gz"),
  list(prot="CGREF1", outcome="Breast_GCST90018757", snps=c("2:27101856"), file="*CGREF1*.txt.gz"),
  list(prot="INHBB_ActivinB", outcome="Breast_GCST90018757", snps=c("2:120548864"), file="*16746_12_INHBB_Activin_B*.txt.gz")
)

results_list <- list()

for (t in targets) {
  cat("\n========================================\n")
  cat(sprintf("Processing %s -> %s\n", t$prot, t$outcome))
  
  decode_file <- Sys.glob(file.path(decode_dir, t$file))
  if (length(decode_file) == 0) {
    cat(sprintf("ERROR: Could not find deCODE file for %s\n", t$prot))
    next
  }
  decode_file <- decode_file[1]
  
  outcome_file <- file.path(project_dir, "data", "cancer_gwas", paste0(t$outcome, ".h.tsv.gz"))
  
  # Build awk condition for chr/pos
  awk_conds <- c()
  for (s in t$snps) {
    parts <- strsplit(s, ":")[[1]]
    chr <- parts[1]
    pos <- parts[2]
    # Check $1 == chr or $1 == chr#
    cond <- sprintf("(($1==\"%s\" || $1==\"chr%s\") && $2==\"%s\")", chr, chr, pos)
    awk_conds <- c(awk_conds, cond)
  }
  awk_filter <- paste(awk_conds, collapse=" || ")
  
  cmd <- sprintf("gzcat '%s' | awk 'NR==1 || %s'", decode_file, awk_filter)
  
  pqtl_raw <- fread(cmd = cmd)
  
  if (nrow(pqtl_raw) == 0) {
    cat(sprintf("WARNING: Target positions not found in deCODE file for %s\n", t$prot))
    next
  }
  
  # Create a unified chr:pos SNP ID
  pqtl_raw[, Chrom := gsub("chr", "", as.character(Chrom), ignore.case = TRUE)]
  pqtl_raw[, SNP_cp := paste(Chrom, Pos, sep=":")]
  
  exp_dat <- format_data(
    as.data.frame(pqtl_raw),
    type = "exposure",
    snp_col = "SNP_cp",
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
  
  cat(sprintf("Extracting from Outcome data (%s)...\n", t$outcome))
  out_raw <- fread(outcome_file)
  
  # Create unified chr:pos for outcome
  out_raw[, SNP_cp := paste(chromosome, base_pair_location, sep=":")]
  
  shared_snps <- intersect(exp_dat$SNP, out_raw$SNP_cp)
  exp_dat <- exp_dat[exp_dat$SNP %in% shared_snps, ]
  
  if(nrow(exp_dat) == 0) {
    cat("WARNING: Target positions not present in outcome data.\n")
    next
  }
  
  out_subset <- out_raw[out_raw$SNP_cp %in% exp_dat$SNP, ]
  out_dat <- format_data(
    as.data.frame(out_subset),
    type = "outcome",
    snp_col = "SNP_cp",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
  )
  out_dat$outcome <- gsub("_GCST.*", "", t$outcome)
  
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
  
  gc()
}

if (length(results_list) > 0) {
  final_results <- bind_rows(results_list)
  write.csv(final_results, file.path(out_dir, "decode_batch_mr_results_v2.csv"), row.names = FALSE)
  cat("\nBatch MR results saved to", file.path(out_dir, "decode_batch_mr_results_v2.csv"), "\n")
} else {
  cat("\nNo results generated.\n")
}
