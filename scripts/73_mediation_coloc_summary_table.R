#!/usr/bin/env Rscript

# 73_mediation_coloc_summary_table.R
# Formats the integrated mediation colocalization evidence into a clean supplementary table.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

in_file <- "results/validation/priority_mediation_integrated_coloc_evidence.csv"
out_file <- "results/tables/STable_mediation_coloc_3leg.csv"

if (!file.exists(in_file)) stop("Missing: ", in_file)

dat <- fread(in_file)

res <- dat %>%
  select(
    Protein = protein,
    Metabolite = metabolite,
    Cancer = cancer,
    `Proportion Mediated (%)` = prop_med_pct,
    `Colocalization PPH4 (Protein -> Cancer)` = PP.H4_protein_cancer,
    `Colocalization PPH4 (Protein -> Metabolite)` = PP.H4_protein_metabolite,
    `Colocalization PPH4 (Metabolite -> Cancer)` = PP.H4_metabolite_cancer,
    `Evidence Class` = coloc_evidence_class,
    `Interpretation Note` = evidence_note
  ) %>%
  mutate(across(starts_with("Colocalization PPH4"), ~ round(as.numeric(.), 3))) %>%
  arrange(desc(`Proportion Mediated (%)`))

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
fwrite(res, out_file)

cat("Successfully generated concise 3-leg mediation colocalization table:\n")
cat(" ->", out_file, "\n")
print(res)
