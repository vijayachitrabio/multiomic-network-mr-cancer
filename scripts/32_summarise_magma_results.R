#!/usr/bin/env Rscript

# Script 31: Summarise MAGMA gene-level results for manuscript use

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- normalizePath(".")
out_dir <- file.path(project_dir, "results", "pathway", "magma")
results_dir <- file.path(out_dir, "results")
gene_loc_path <- file.path(project_dir, "..", "uterine_fibroids", "magma_inputs", "NCBI37.3.gene.loc")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bonf <- 2.5e-6

gene_map <- fread(
  gene_loc_path,
  col.names = c("GENE", "CHR_LOC", "START_LOC", "STOP_LOC", "STRAND", "SYMBOL")
)[, .(GENE, SYMBOL)]
gene_map <- unique(gene_map)

load_magma <- function(trait) {
  path <- file.path(results_dir, paste0(trait, "_genes.genes.out"))
  dt <- fread(path)
  dt <- merge(dt, gene_map, by = "GENE", all.x = TRUE)
  dt[, SYMBOL := fifelse(is.na(SYMBOL) | SYMBOL == "", as.character(GENE), SYMBOL)]
  dt[, trait := trait]
  dt[, LOG10P := -log10(pmax(P, 1e-300))]
  dt[, bonferroni_sig := P < bonf]
  setorder(dt, P)
  dt[]
}

breast <- load_magma("breast")
endometrial <- load_magma("endometrial")

fwrite(breast, file.path(out_dir, "breast_magma_all_genes.csv"))
fwrite(endometrial, file.path(out_dir, "endometrial_magma_all_genes.csv"))
fwrite(breast[1:100], file.path(out_dir, "breast_magma_top100.csv"))
fwrite(endometrial[1:100], file.path(out_dir, "endometrial_magma_top100.csv"))

hit_map <- data.table(
  protein = c("SNX15","EFNA1","FGF5","UMOD","SWAP70","ATRAID","TNFRSF6B",
              "ITIH3","KLB","PM20D1","TSPAN8","FGFR4","IL34","APOE",
              "CGREF1","INHBB","ABO"),
  magma_symbol = c("SNX15","EFNA1","FGF5","UMOD","SWAP70","ATRAID","TNFRSF6B",
                   "ITIH3","RIMKLB","PM20D1","TSPAN8","FGFR4","IL34","APOE",
                   "CGREF1","INHBB","ABO")
)

lookup_trait <- function(dt, trait_name) {
  x <- merge(
    hit_map,
    dt[, .(magma_symbol = SYMBOL, NSNPS, NPARAM, ZSTAT, P, LOG10P, bonferroni_sig)],
    by = "magma_symbol",
    all.x = TRUE
  )
  setnames(
    x,
    old = c("NSNPS", "NPARAM", "ZSTAT", "P", "LOG10P", "bonferroni_sig"),
    new = paste0(trait_name, c("_nsnps", "_nparam", "_zstat", "_p", "_log10p", "_bonferroni_sig"))
  )
  x
}

lookup <- merge(
  lookup_trait(breast, "breast"),
  lookup_trait(endometrial, "endometrial"),
  by = c("protein", "magma_symbol"),
  all = TRUE
)
lookup[, breast_rank := frank(breast_p, ties.method = "min", na.last = "keep")]
lookup[, endometrial_rank := frank(endometrial_p, ties.method = "min", na.last = "keep")]
setorder(lookup, breast_p, endometrial_p)

fwrite(lookup, file.path(out_dir, "mr_hit_gene_magma_lookup.csv"))

summary_lines <- c(
  "# MAGMA Summary",
  "",
  sprintf("- Breast cancer: %d genes tested; %d Bonferroni-significant genes (P < %.1e).",
          nrow(breast), sum(breast$bonferroni_sig, na.rm = TRUE), bonf),
  sprintf("- Endometrial cancer: %d genes tested; %d Bonferroni-significant genes (P < %.1e).",
          nrow(endometrial), sum(endometrial$bonferroni_sig, na.rm = TRUE), bonf),
  "",
  "## MR-hit lookup",
  "",
  "The table `mr_hit_gene_magma_lookup.csv` records breast and endometrial MAGMA support",
  "for the 17 FDR-significant protein MR hits. `KLB` is mapped to the MAGMA gene symbol",
  "`RIMKLB` in `NCBI37.3.gene.loc`.",
  "",
  "## Key files",
  "",
  "- `results/pathway/magma/breast_magma_all_genes.csv`",
  "- `results/pathway/magma/endometrial_magma_all_genes.csv`",
  "- `results/pathway/magma/mr_hit_gene_magma_lookup.csv`",
  "- `results/pathway/magma/breast_magma_top100.csv`",
  "- `results/pathway/magma/endometrial_magma_top100.csv`"
)
writeLines(summary_lines, file.path(out_dir, "MAGMA_SUMMARY.md"))

message("MAGMA summaries written to: ", out_dir)
