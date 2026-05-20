#!/usr/bin/env Rscript

# TISCH2 BRCA_EMTAB8107 scRNA-seq target-gene cell-type localization.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
data_dir <- file.path(project_dir, "data", "scrna_tisch")
expr_dir <- file.path(data_dir, "BRCA_EMTAB8107_Expression")
out_dir <- file.path(project_dir, "results", "scrna")
fig_dir <- file.path(project_dir, "results", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

targets <- c("EFNA1", "TNFRSF6B", "ATRAID", "ITIH3", "IL34", "FGF5", "APOE")

files <- list(
  malignancy = file.path(expr_dir, "BRCA_EMTAB8107_expression_Celltype_malignancy.txt"),
  major_lineage = file.path(expr_dir, "BRCA_EMTAB8107_expression_Celltype_majorlineage.txt"),
  minor_lineage = file.path(expr_dir, "BRCA_EMTAB8107_expression_Celltype_minorlineage.txt"),
  cluster = file.path(expr_dir, "BRCA_EMTAB8107_expression_Cluster.txt")
)

read_avg_expr <- function(path, level) {
  x <- fread(path)
  setnames(x, 1, "gene")
  x <- x[gene %in% targets]
  long <- melt(x, id.vars = "gene", variable.name = "celltype", value.name = "avg_expression")
  long[, annotation_level := level]
  long[]
}

avg_long <- rbindlist(Map(read_avg_expr, files, names(files)), fill = TRUE)
avg_long[, gene_z := as.numeric(scale(avg_expression)), by = .(annotation_level, gene)]
avg_long[, rank_within_gene := frank(-avg_expression, ties.method = "min"), by = .(annotation_level, gene)]
fwrite(avg_long, file.path(out_dir, "tisch_brca_emtab8107_target_avg_expression_long.csv"))

top_localization <- avg_long[rank_within_gene == 1,
  .(top_celltype = paste(celltype, collapse = "; "),
    top_avg_expression = max(avg_expression, na.rm = TRUE)),
  by = .(annotation_level, gene)
]
fwrite(top_localization, file.path(out_dir, "tisch_brca_emtab8107_target_top_celltypes.csv"))

meta_file <- file.path(data_dir, "BRCA_EMTAB8107_CellMetainfo_table.tsv")
meta <- fread(meta_file)
cell_counts <- rbindlist(list(
  meta[, .N, by = .(celltype = `Celltype (malignancy)`)][, annotation_level := "malignancy"],
  meta[, .N, by = .(celltype = `Celltype (major-lineage)`)][, annotation_level := "major_lineage"],
  meta[, .N, by = .(celltype = `Celltype (minor-lineage)`)][, annotation_level := "minor_lineage"],
  meta[, .N, by = .(celltype = as.character(Cluster))][, annotation_level := "cluster"]
), fill = TRUE)
setcolorder(cell_counts, c("annotation_level", "celltype", "N"))
fwrite(cell_counts, file.path(out_dir, "tisch_brca_emtab8107_celltype_counts.csv"))

diff_file <- file.path(data_dir, "BRCA_EMTAB8107_AllDiffGenes_table.tsv")
diff <- fread(diff_file)
diff_target <- diff[Gene %in% targets]
fwrite(diff_target, file.path(out_dir, "tisch_brca_emtab8107_target_diffgene_hits.csv"))

plot_level <- "major_lineage"
pdat <- avg_long[annotation_level == plot_level]
pdat[, celltype := factor(celltype, levels = unique(celltype[order(-avg_expression)]))]
pdat[, gene := factor(gene, levels = targets[targets %in% pdat$gene])]

p <- ggplot(pdat, aes(x = celltype, y = gene, fill = gene_z)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#3568a6", mid = "white", high = "#b4443c", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Gene z-score") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "fig11_tisch_scrna_majorlineage_heatmap.png"), p, width = 7, height = 3.8, dpi = 300)
ggsave(file.path(fig_dir, "fig11_tisch_scrna_majorlineage_heatmap.pdf"), p, width = 7, height = 3.8)

plot_level <- "minor_lineage"
pdat <- avg_long[annotation_level == plot_level]
pdat[, gene := factor(gene, levels = targets[targets %in% pdat$gene])]
p <- ggplot(pdat, aes(x = celltype, y = gene, fill = gene_z)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#3568a6", mid = "white", high = "#b4443c", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Gene z-score") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "fig12_tisch_scrna_minorlineage_heatmap.png"), p, width = 7, height = 3.8, dpi = 300)
ggsave(file.path(fig_dir, "fig12_tisch_scrna_minorlineage_heatmap.pdf"), p, width = 7, height = 3.8)

message("Saved TISCH2 scRNA-seq validation outputs to: ", out_dir)
