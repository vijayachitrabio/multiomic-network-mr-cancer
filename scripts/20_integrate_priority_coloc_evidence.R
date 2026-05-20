#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)

paths <- fread("results/validation/top_mediation_paths_for_followup.csv")
protein_legs <- fread("results/validation/priority_mediation_leg_coloc_wide.csv")
met_cancer <- fread("results/validation/metabolite_cancer_coloc_best_loci.csv")

setnames(
  met_cancer,
  old = c("exposure", "outcome", "PP.H4", "PP.H3", "lead_snp"),
  new = c("metabolite", "cancer", "PP.H4_metabolite_cancer", "PP.H3_metabolite_cancer", "metabolite_cancer_lead_snp")
)

out <- merge(paths, protein_legs, by = c("protein", "metabolite", "cancer"), all.x = TRUE)
out <- merge(
  out,
  met_cancer[, .(metabolite, cancer, PP.H4_metabolite_cancer, PP.H3_metabolite_cancer, metabolite_cancer_lead_snp)],
  by = c("metabolite", "cancer"),
  all.x = TRUE
)

strong <- function(x) !is.na(x) & x >= 0.8
moderate <- function(x) !is.na(x) & x >= 0.5 & x < 0.8

out[, coloc_evidence_class := fifelse(
  strong(PP.H4_protein_cancer) & strong(PP.H4_protein_metabolite) & strong(PP.H4_metabolite_cancer),
  "complete_three_leg_coloc_strong",
  fifelse(
    strong(PP.H4_protein_metabolite) & strong(PP.H4_metabolite_cancer),
    "mechanistic_pm_and_mc_coloc_strong_pc_not_strong",
    fifelse(
      strong(PP.H4_protein_cancer) & strong(PP.H4_metabolite_cancer),
      "pc_and_mc_coloc_strong_pm_not_strong",
      fifelse(
        strong(PP.H4_metabolite_cancer) &
          (moderate(PP.H4_protein_cancer) | moderate(PP.H4_protein_metabolite)),
        "metabolite_cancer_strong_one_protein_leg_moderate",
        "incomplete_or_weak_coloc"
      )
    )
  )
)]

setcolorder(out, c(
  "protein", "metabolite", "cancer", "validation_priority", "coloc_evidence_class",
  "p_indirect", "prop_med_pct",
  "PP.H4_protein_cancer", "PP.H4_protein_metabolite", "PP.H4_metabolite_cancer",
  "PP.H3_protein_cancer", "PP.H3_protein_metabolite", "PP.H3_metabolite_cancer",
  "n_snps_protein_cancer", "n_snps_protein_metabolite",
  "metabolite_cancer_lead_snp",
  setdiff(names(out), c(
    "protein", "metabolite", "cancer", "validation_priority", "coloc_evidence_class",
    "p_indirect", "prop_med_pct",
    "PP.H4_protein_cancer", "PP.H4_protein_metabolite", "PP.H4_metabolite_cancer",
    "PP.H3_protein_cancer", "PP.H3_protein_metabolite", "PP.H3_metabolite_cancer",
    "n_snps_protein_cancer", "n_snps_protein_metabolite",
    "metabolite_cancer_lead_snp"
  ))
))

setorder(out, coloc_evidence_class, validation_priority, -PP.H4_metabolite_cancer)
fwrite(out, "results/validation/priority_mediation_integrated_coloc_evidence.csv")

summary <- out[, .N, by = coloc_evidence_class][order(coloc_evidence_class)]
fwrite(summary, "results/validation/priority_mediation_integrated_coloc_summary.csv")

cat("Wrote integrated priority coloc evidence:\n")
cat("  results/validation/priority_mediation_integrated_coloc_evidence.csv\n")
cat("  results/validation/priority_mediation_integrated_coloc_summary.csv\n\n")
print(out[, .(
  protein, metabolite, validation_priority, coloc_evidence_class,
  PP.H4_protein_cancer, PP.H4_protein_metabolite, PP.H4_metabolite_cancer,
  p_indirect, prop_med_pct
)])
