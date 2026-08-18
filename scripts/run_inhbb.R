suppressMessages({ library(TwoSampleMR); library(data.table); library(dplyr) })
decode_file <- "../SUMMARY DATA-proteins/16746_12_INHBB_Activin_B.txt.gz"
outcome_file <- "data/cancer_gwas/Breast_GCST90018757.h.tsv.gz"
cmd <- sprintf("gzcat '%s' | awk 'NR==1 || (($1==\"2\" || $1==\"chr2\") && $2==\"120548864\")'", decode_file)
pqtl_raw <- fread(cmd = cmd)
pqtl_raw[, Chrom := gsub("chr", "", as.character(Chrom), ignore.case = TRUE)]
pqtl_raw[, SNP_cp := paste(Chrom, Pos, sep=":")]
exp_dat <- format_data(as.data.frame(pqtl_raw), type="exposure", snp_col="SNP_cp", beta_col="Beta", se_col="SE", effect_allele_col="effectAllele", other_allele_col="otherAllele", eaf_col="ImpMAF", pval_col="Pval", chr_col="Chrom", pos_col="Pos")
exp_dat$exposure <- "INHBB_Activin_B"
out_raw <- fread(outcome_file)
out_raw[, SNP_cp := paste(chromosome, base_pair_location, sep=":")]
shared_snps <- intersect(exp_dat$SNP, out_raw$SNP_cp)
exp_dat <- exp_dat[exp_dat$SNP %in% shared_snps, ]
out_subset <- out_raw[out_raw$SNP_cp %in% exp_dat$SNP, ]
out_dat <- format_data(as.data.frame(out_subset), type="outcome", snp_col="SNP_cp", beta_col="beta", se_col="standard_error", effect_allele_col="effect_allele", other_allele_col="other_allele", eaf_col="effect_allele_frequency", pval_col="p_value")
out_dat$outcome <- "Breast"
harm_dat <- harmonise_data(exposure_dat=exp_dat, outcome_dat=out_dat, action=2)
if(nrow(harm_dat) > 0 && sum(harm_dat$mr_keep)>0) {
  mr_res <- mr(harm_dat, method_list=c("mr_wald_ratio", "mr_ivw"))
  print(generate_odds_ratios(mr_res))
} else { cat("Failed harmonization.\n") }
