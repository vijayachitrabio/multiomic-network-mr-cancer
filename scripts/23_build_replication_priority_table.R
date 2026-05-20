#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"

phase2_path <- file.path(project_dir, "results", "phase2_protein_cancer", "protein_cancer_mr_results_significant.csv")
story_path  <- file.path(project_dir, "results", "validation", "manuscript_story_tiers.csv")
out_dir     <- file.path(project_dir, "results", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

phase2 <- fread(phase2_path)
story  <- fread(story_path)

phase2[, cancer_label := fifelse(
  grepl("^Breast", outcome), "Breast",
  fifelse(grepl("^Endometrial", outcome), "Endometrial", outcome)
)]

phase2[, discovery_rank := frank(fdr, ties.method = "first"), by = cancer_label]

priority_map <- data.table(
  exposure = c("SNX15", "EFNA1", "FGF5", "TNFRSF6B", "APOE", "IL34", "ABO"),
  replication_priority = c(
    "very_high",
    "very_high",
    "very_high",
    "high",
    "high",
    "high",
    "secondary_endometrial"
  ),
  reason_for_replication = c(
    "Top breast hit; single-SNP; no protein-side coloc yet",
    "Top breast hit; single-SNP; biologically strong but needs external support",
    "Top breast hit; single-SNP; weak literature overlap",
    "Best partial protein-cancer plus metabolite-cancer support",
    "Most coherent partial APOE-glycine-breast mechanistic story",
    "Largest mediation proportion but weak current protein-side coloc",
    "Only significant endometrial protein hit"
  )
)

out <- merge(priority_map, phase2, by = "exposure", all.x = TRUE)
out <- merge(
  out,
  story[, .(protein, metabolite, manuscript_tier, coloc_evidence_class, interpretation)],
  by.x = "exposure",
  by.y = "protein",
  all.x = TRUE
)

setcolorder(out, c(
  "replication_priority", "exposure", "cancer_label", "outcome", "method", "nsnp",
  "or", "or_lci95", "or_uci95", "pval", "fdr", "F_stat_mean", "metabolite",
  "manuscript_tier", "coloc_evidence_class", "reason_for_replication", "interpretation"
))

out[, priority_order := fifelse(
  replication_priority == "very_high", 1L,
  fifelse(replication_priority == "high", 2L, 3L)
)]
setorder(out, priority_order, fdr)
out[, priority_order := NULL]

fwrite(out, file.path(out_dir, "replication_priority_targets.csv"))

cat("Wrote:\n")
cat("  results/validation/replication_priority_targets.csv\n")
print(out[, .(
  replication_priority, exposure, cancer_label, method, nsnp,
  or = round(or, 4), fdr = signif(fdr, 3), metabolite, manuscript_tier
)])
