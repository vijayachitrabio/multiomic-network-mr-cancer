#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

proj <- "."
out_dir <- file.path(proj, "results", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fmt_p <- function(x) {
  ifelse(is.na(x), NA_character_,
         ifelse(x < 1e-3, sprintf("%.2e", x), sprintf("%.3f", x)))
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

master <- fread(file.path(proj, "results", "tables", "STable_master_evidence.csv"))

tcga_expr_path <- file.path(proj, "results", "tcga_immune", "tcga_brca_target_tumour_vs_normal.csv")
tcga_imm_path  <- file.path(proj, "results", "tcga_immune", "tcga_brca_target_immune_correlations.csv")
cptac_imm_path <- file.path(proj, "results", "cptac", "cptac_brca_target_protein_immune_correlations.csv")
scrna_path     <- file.path(proj, "results", "scrna", "tisch_brca_emtab8107_target_top_celltypes.csv")

tcga_expr <- if (file.exists(tcga_expr_path)) fread(tcga_expr_path) else data.table()
tcga_imm  <- if (file.exists(tcga_imm_path)) fread(tcga_imm_path) else data.table()
cptac_imm <- if (file.exists(cptac_imm_path)) fread(cptac_imm_path) else data.table()
scrna     <- if (file.exists(scrna_path)) fread(scrna_path) else data.table()

if (nrow(tcga_expr)) {
  tcga_expr[, tcga_tumour_expression := paste0(
    ifelse(logFC_tumour_vs_normal > 0, "higher in tumour", "lower in tumour"),
    " (logFC=", fmt_num(logFC_tumour_vs_normal, 2),
    ", FDR=", fmt_p(fdr_tumour_vs_normal), ")"
  )]
  tcga_expr <- tcga_expr[, .(protein = gene, tcga_tumour_expression)]
}

top_by_fdr <- function(dt, gene_col = "gene", feature_col = "signature", fdr_max = 0.05) {
  if (!nrow(dt)) return(data.table())
  dt <- dt[!is.na(fdr) & fdr <= fdr_max]
  if (!nrow(dt)) return(data.table())
  setorderv(dt, c(gene_col, "fdr"))
  dt[, .SD[1], by = gene_col][, {
    z <- paste0(get(feature_col), " rho=", fmt_num(rho, 2), ", FDR=", fmt_p(fdr))
    data.table(protein = get(gene_col), top_immune_signal = z)
  }]
}

tcga_imm_top <- top_by_fdr(tcga_imm, "gene", "signature", fdr_max = 0.05)
if (nrow(tcga_imm_top)) setnames(tcga_imm_top, "top_immune_signal", "tcga_top_immune_signal")

cptac_imm_top <- top_by_fdr(cptac_imm, "gene", "immune_feature", fdr_max = 0.10)
if (nrow(cptac_imm_top)) setnames(cptac_imm_top, "top_immune_signal", "cptac_top_immune_signal")

if (nrow(scrna)) {
  scrna_wide <- dcast(
    scrna[annotation_level %in% c("malignancy", "major_lineage", "minor_lineage")],
    gene ~ annotation_level,
    value.var = "top_celltype"
  )
  setnames(scrna_wide, "gene", "protein")
} else {
  scrna_wide <- data.table()
}

hpa <- data.table(
  protein = c("TNFRSF6B", "EFNA1", "FGF5", "UMOD", "ATRAID", "ABO"),
  hpa_context = c(
    "Secreted; immune-cell enriched in plasmacytoid DCs",
    "Secreted; trophoblast/migratory-cell enhancement; not breast-specific",
    "Secreted growth factor; not immune/breast-specific",
    "Secreted/membrane; renal tubular context, not breast-specific",
    "Broad membrane/intracellular expression; not cell-type specific",
    "Secreted/membrane/intracellular; platelet/immune/renal epithelial context"
  )
)

summary <- copy(master)
summary[, mr_summary := paste0(
  cancer_mr, " OR=", fmt_num(mr_or, 3),
  " (FDR=", fmt_p(mr_fdr), ", nsnp=", mr_nsnp, ")"
)]
summary[, coloc_summary := paste0(coloc_verdict, " PPH4=", fmt_num(coloc_PPH4_best, 3), " (", coloc_method_best, ")")]
summary[, magma_summary := fifelse(
  is.na(magma_breast_p), NA_character_,
  paste0("MAGMA p=", fmt_p(magma_breast_p), ifelse(magma_breast_bonf, " Bonferroni", ""))
)]
summary[, mediation_summary := fifelse(
  is.na(med_metabolite),
  "none retained",
  paste0(med_metabolite, "; WM p=", fmt_p(med_p_indirect_wm),
         "; prop=", fmt_num(med_prop_pct, 1), "%; ",
         ifelse(!is.na(med_supported) & med_supported == TRUE,
                "supported by MR sensitivity",
                "not supported by WM sensitivity"))
)]
summary[, final_claim := fifelse(
  tier_short == "T1",
  "Colocalization-supported protein-cancer candidate",
  fifelse(
    tier_short == "T2a",
    "MR and MAGMA-supported candidate without coloc support",
    fifelse(
      tier_short == "T2b",
      "Moderate coloc / secondary candidate",
      "Suggestive MR discovery finding"
    )
  )
)]

keep <- summary[, .(
  protein, cancer = cancer_mr, tier_short, final_claim,
  mr_summary, coloc_summary, magma_summary, er_pattern,
  mediation_summary, steiger_ok,
  druggability = fifelse(n_drugs > 0, paste0(n_drugs, " drug entries"), "no registered drug program")
)]

for (dt in list(tcga_expr, tcga_imm_top, cptac_imm_top, scrna_wide, hpa)) {
  if (nrow(dt)) keep <- merge(keep, dt, by = "protein", all.x = TRUE)
}

keep[, immune_context := mapply(function(tcga, cptac) {
  vals <- c(tcga, cptac)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  paste(vals, collapse = "; ")
}, tcga_top_immune_signal, cptac_top_immune_signal, USE.NAMES = FALSE)]

setorder(keep, tier_short, protein)

out_csv <- file.path(out_dir, "integrated_validation_summary_2026-05-20.csv")
fwrite(keep, out_csv)

md <- c(
  "# Integrated Validation Summary",
  "",
  "Status date: 2026-05-20",
  "",
  "This table integrates MR, colocalization, MAGMA, ER-subtype, mediation, druggability, TCGA, CPTAC, TISCH scRNA, and HPA context for the 17 FDR-significant protein-cancer associations.",
  "",
  "Interpretation rule: TCGA/CPTAC/scRNA/HPA evidence is used as biological context only. It is not treated as causal validation.",
  "",
  "## Main Claim Groups",
  "",
  paste0("- Colocalization-supported candidates: ", paste(keep[final_claim == "Colocalization-supported protein-cancer candidate", protein], collapse = ", "), "."),
  paste0("- MR and MAGMA-supported without colocalization: ", paste(keep[final_claim == "MR and MAGMA-supported candidate without coloc support", protein], collapse = ", "), "."),
  paste0("- Moderate/secondary candidates: ", paste(keep[final_claim == "Moderate coloc / secondary candidate", protein], collapse = ", "), "."),
  paste0("- Suggestive discovery findings: ", paste(keep[final_claim == "Suggestive MR discovery finding", protein], collapse = ", "), "."),
  "",
  "## Manuscript-Ready Interpretation",
  "",
  "The integrated validation layer supports a cautious prioritisation framework. EFNA1, FGF5, UMOD, TNFRSF6B, ATRAID, and ABO are the strongest protein-cancer candidates because their MR associations are supported by protein-cancer colocalization. SNX15 and PM20D1 remain important because they combine MR significance with MAGMA gene-level support, although colocalization indicates distinct causal variants. TCGA, CPTAC, TISCH, and HPA evidence provide biological context, especially immune/stromal support for TNFRSF6B, APOE, IL34, and ITIH3 and malignant-cell context for EFNA1/ATRAID, but these data should not be framed as causal replication.",
  "",
  "## Compact Table",
  "",
  "| Protein | Final claim | MR | Coloc | ER pattern | Mediation | TCGA expression | TCGA/CPTAC immune context | scRNA context | HPA context |",
  "|---|---|---|---|---|---|---|---|---|---|"
)

table_lines <- keep[, paste0(
  "| ", protein,
  " | ", final_claim,
  " | ", mr_summary,
  " | ", coloc_summary,
  " | ", fifelse(is.na(er_pattern), "", er_pattern),
  " | ", mediation_summary,
  " | ", fifelse(is.na(tcga_tumour_expression), "", tcga_tumour_expression),
  " | ", immune_context,
  " | ", fifelse(is.na(major_lineage), "", paste0(major_lineage, "; ", malignancy)),
  " | ", fifelse(is.na(hpa_context), "", hpa_context),
  " |"
)]
md <- c(md, table_lines)

out_md <- file.path(out_dir, "INTEGRATED_VALIDATION_SUMMARY_2026-05-20.md")
writeLines(md, out_md)

message("Wrote: ", out_csv)
message("Wrote: ", out_md)
