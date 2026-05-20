#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)

queue_path <- "results/validation/mediation_validation_queue.csv"
coloc_path <- "results/validation/metabolite_cancer_coloc_best_loci.csv"

if (!file.exists(queue_path)) {
  stop("Missing mediation validation queue: ", queue_path)
}
if (!file.exists(coloc_path)) {
  stop("Missing metabolite-cancer coloc results: ", coloc_path)
}

queue <- fread(queue_path)
mc_coloc <- fread(coloc_path)

setnames(
  mc_coloc,
  old = c("exposure", "PP.H4", "PP.H3", "lead_snp", "n_snps"),
  new = c("metabolite", "met_cancer_coloc_pph4", "met_cancer_coloc_pph3",
          "met_cancer_coloc_lead_snp", "met_cancer_coloc_n_snps")
)

keep_cols <- c(
  "metabolite", "outcome", "met_cancer_coloc_lead_snp",
  "met_cancer_coloc_n_snps", "met_cancer_coloc_pph3", "met_cancer_coloc_pph4"
)
mc_coloc <- mc_coloc[, ..keep_cols]

status <- merge(
  queue,
  mc_coloc,
  by.x = c("metabolite", "cancer"),
  by.y = c("metabolite", "outcome"),
  all.x = TRUE
)

status[, met_cancer_coloc_tier := fifelse(
  is.na(met_cancer_coloc_pph4), "not_available",
  fifelse(met_cancer_coloc_pph4 >= 0.8, "strong_pph4_ge_0.8",
    fifelse(met_cancer_coloc_pph4 >= 0.7, "moderate_pph4_ge_0.7", "weak_pph4_lt_0.7")
  )
)]

status[, protein_coloc_status := "pending_full_regional_pqtl_sumstats"]

status[, manuscript_status := fifelse(
  validation_priority == "tier1_fast_track" & met_cancer_coloc_pph4 >= 0.8,
  "top_priority_but_needs_protein_coloc",
  fifelse(
    validation_priority %in% c("tier2_coloc_single_snp_caution", "tier3_fdr_only_or_direction_caution") &
      met_cancer_coloc_pph4 >= 0.8,
    "interesting_supportive_metabolite_leg_but_protein_side_caution",
    fifelse(
      met_cancer_coloc_pph4 < 0.8,
      "lower_confidence_until_coloc_improves",
      "needs_more_validation"
    )
  )
)]

setcolorder(status, c(
  "protein", "metabolite", "cancer", "validation_priority", "manuscript_status",
  "direct_consistent", "p_indirect", "prop_med_pct",
  "met_cancer_coloc_tier", "met_cancer_coloc_pph4", "met_cancer_coloc_pph3",
  "met_cancer_coloc_lead_snp", "protein_coloc_status",
  setdiff(names(status), c(
    "protein", "metabolite", "cancer", "validation_priority", "manuscript_status",
    "direct_consistent", "p_indirect", "prop_med_pct",
    "met_cancer_coloc_tier", "met_cancer_coloc_pph4", "met_cancer_coloc_pph3",
    "met_cancer_coloc_lead_snp", "protein_coloc_status"
  ))
))

setorder(status, priority_rank, -met_cancer_coloc_pph4, p_indirect)

fwrite(status, "results/validation/mediation_with_coloc_status.csv")

summary <- status[, .N, by = .(manuscript_status, met_cancer_coloc_tier)]
setorder(summary, manuscript_status, met_cancer_coloc_tier)
fwrite(summary, "results/validation/mediation_with_coloc_status_summary.csv")

top <- status[
  manuscript_status %in% c(
    "top_priority_but_needs_protein_coloc",
    "interesting_supportive_metabolite_leg_but_protein_side_caution"
  ),
  .(
    protein, metabolite, cancer, validation_priority, manuscript_status,
    p_indirect, prop_med_pct, met_cancer_coloc_pph4,
    met_cancer_coloc_lead_snp, protein_coloc_status
  )
]
fwrite(top, "results/validation/top_mediation_paths_for_followup.csv")

cat("Wrote mediation coloc status files:\n")
cat("  results/validation/mediation_with_coloc_status.csv\n")
cat("  results/validation/mediation_with_coloc_status_summary.csv\n")
cat("  results/validation/top_mediation_paths_for_followup.csv\n\n")
print(top)
