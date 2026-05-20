#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(utils)
})

dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

add_global_fdr <- function(x) {
  x$global_fdr <- p.adjust(x$pval, method = "fdr")
  x$bonferroni_sig <- x$pval < 0.05 / nrow(x)
  x
}

instrument_tier <- function(nsnp) {
  ifelse(nsnp >= 3, "multi_snp_stronger",
    ifelse(nsnp == 2, "two_snp_limited", "single_snp_wald"))
}

pc <- add_global_fdr(read_csv("results/phase2_protein_cancer/protein_cancer_mr_results_full.csv"))
pm <- add_global_fdr(read_csv("results/phase3_protein_metabolite/protein_metabolite_mr_results_full.csv"))
mc <- add_global_fdr(read_csv("results/phase4_metabolite_cancer/metabolite_cancer_mr_results_full.csv"))
med <- read_csv("results/mediation_mr_results.csv")

pc_sig <- pc[pc$fdr < 0.05, ]
pc_sig$instrument_tier <- instrument_tier(pc_sig$nsnp)
pc_sig$validation_priority <- ifelse(
  pc_sig$nsnp >= 3 & pc_sig$global_fdr < 0.05,
  "priority_coloc_plus_sensitivity",
  ifelse(pc_sig$global_fdr < 0.05, "priority_coloc_single_or_two_snp", "suggestive")
)
pc_sig <- pc_sig[order(pc_sig$fdr, pc_sig$pval), ]

pm_sig <- pm[pm$fdr < 0.05, ]
pm_sig$instrument_tier <- instrument_tier(pm_sig$nsnp)
pm_sig$validation_priority <- ifelse(
  pm_sig$nsnp >= 3 & pm_sig$global_fdr < 0.05,
  "stronger_mediation_candidate",
  ifelse(pm_sig$global_fdr < 0.05, "candidate_needs_coloc", "suggestive_only")
)
pm_sig <- pm_sig[order(pm_sig$global_fdr, pm_sig$fdr), ]

mc_sig <- mc[mc$fdr < 0.05, ]
mc_sig$instrument_tier <- instrument_tier(mc_sig$nsnp)
mc_sig$validation_priority <- ifelse(
  mc_sig$global_fdr < 0.05,
  "priority_coloc_plus_sensitivity",
  "fdr_only_needs_confirmation"
)
mc_sig <- mc_sig[order(mc_sig$fdr, mc_sig$pval), ]

nearest_by_fdr <- function(query_fdr, candidates) {
  if (nrow(candidates) == 0) return(rep(NA_real_, 6))
  idx <- which.min(abs(candidates$fdr - query_fdr))
  row <- candidates[idx, ]
  c(row$nsnp, row$method, row$pval, row$fdr, row$global_fdr, row$instrument_tier)
}

pc_sig$instrument_tier <- instrument_tier(pc_sig$nsnp)
pm_sig$instrument_tier <- instrument_tier(pm_sig$nsnp)
mc_sig$instrument_tier <- instrument_tier(mc_sig$nsnp)

med_rows <- lapply(seq_len(nrow(med)), function(i) {
  row <- med[i, ]
  pc_match <- pc_sig[
    pc_sig$exposure == row$protein &
      pc_sig$outcome == row$cancer &
      pc_sig$method == row$p2_method,
  ]
  pm_match <- pm_sig[
    pm_sig$exposure == row$protein &
      pm_sig$outcome == row$metabolite,
  ]
  mc_match <- mc_sig[
    mc_sig$exposure == row$metabolite &
      mc_sig$outcome == row$cancer,
  ]

  p2 <- nearest_by_fdr(row$p2_fdr, pc_match)
  p3 <- nearest_by_fdr(row$p3_fdr, pm_match)
  p4 <- nearest_by_fdr(row$p4_fdr, mc_match)

  direct <- isTRUE(row$direct_consistent)
  indirect <- !is.na(row$p_indirect) && row$p_indirect < 0.05
  all_global <- all(as.numeric(c(p2[5], p3[5], p4[5])) < 0.05, na.rm = FALSE)
  any_single <- any(c(p2[6], p3[6], p4[6]) == "single_snp_wald", na.rm = TRUE)

  priority <- if (direct && indirect && all_global && !any_single) {
    "tier1_fast_track"
  } else if (direct && indirect && all_global) {
    "tier2_coloc_single_snp_caution"
  } else if (direct && indirect) {
    "tier3_fdr_only_or_direction_caution"
  } else {
    "lower_priority"
  }

  data.frame(
    protein = row$protein,
    metabolite = row$metabolite,
    cancer = row$cancer,
    direct_consistent = row$direct_consistent,
    p_indirect = row$p_indirect,
    prop_med_pct = row$prop_med_pct,
    p2_method = row$p2_method,
    p2_fdr = row$p2_fdr,
    p2_nsnp = as.numeric(p2[1]),
    p2_global_fdr = as.numeric(p2[5]),
    p2_tier = p2[6],
    p3_fdr = row$p3_fdr,
    p3_nsnp = as.numeric(p3[1]),
    p3_method = p3[2],
    p3_global_fdr = as.numeric(p3[5]),
    p3_tier = p3[6],
    p4_fdr = row$p4_fdr,
    p4_nsnp = as.numeric(p4[1]),
    p4_method = p4[2],
    p4_global_fdr = as.numeric(p4[5]),
    p4_tier = p4[6],
    validation_priority = priority,
    stringsAsFactors = FALSE
  )
})

med_queue <- do.call(rbind, med_rows)
priority_rank <- c(
  tier1_fast_track = 1,
  tier2_coloc_single_snp_caution = 2,
  tier3_fdr_only_or_direction_caution = 3,
  lower_priority = 4
)
med_queue$priority_rank <- priority_rank[med_queue$validation_priority]
med_queue <- med_queue[order(med_queue$priority_rank, med_queue$p_indirect), ]

write.csv(pc_sig, "results/validation/protein_cancer_evidence_tiers.csv", row.names = FALSE)
write.csv(pm_sig, "results/validation/protein_metabolite_evidence_tiers.csv", row.names = FALSE)
write.csv(mc_sig, "results/validation/metabolite_cancer_evidence_tiers.csv", row.names = FALSE)
write.csv(med_queue, "results/validation/mediation_validation_queue.csv", row.names = FALSE)

summary_lines <- c(
  "# Validation Queue Summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "This summary is derived from existing MR result CSVs only. It does not rerun MR.",
  "",
  "## Counts",
  paste0("- Protein -> cancer FDR<0.05 rows: ", nrow(pc_sig)),
  paste0("- Protein -> metabolite FDR<0.05 rows: ", nrow(pm_sig)),
  paste0("- Metabolite -> cancer FDR<0.05 rows: ", nrow(mc_sig)),
  paste0("- Mediation triplets: ", nrow(med_queue)),
  "",
  "## Mediation Priority Counts",
  paste(capture.output(print(table(med_queue$validation_priority))), collapse = "\n"),
  "",
  "## Recommended Fast Start",
  "1. Start coloc and sensitivity checks with `tier1_fast_track` rows.",
  "2. Treat single-SNP protein -> cancer paths as high-risk until coloc is positive.",
  "3. Keep FDR-only or direction-inconsistent mediation rows as suggestive."
)

writeLines(summary_lines, "results/validation/VALIDATION_QUEUE_SUMMARY.md")

cat("Wrote validation queue files to results/validation/\n")
cat("Mediation priorities:\n")
print(table(med_queue$validation_priority))
