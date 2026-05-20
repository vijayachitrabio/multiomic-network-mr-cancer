#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

out_fig <- "results/figures"
out_tab <- "results/network"
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)

priority <- c("EFNA1", "ATRAID", "TNFRSF6B", "APOE", "ITIH3", "IL34", "FGF5")

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(data.frame())
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

master <- read_csv_safe("results/tables/STable7_master_evidence_compact.csv")
ogw <- read_csv_safe("results/opengwas/opengwas_5protein_replication_mr_results.csv")
ogw_inst <- read_csv_safe("results/opengwas/opengwas_7protein_instruments_p5e8.csv")
tcga_expr <- read_csv_safe("results/tcga_immune/tcga_brca_target_tumour_vs_normal.csv")
tcga_imm <- read_csv_safe("results/tcga_immune/tcga_brca_target_immune_correlations.csv")
cptac <- read_csv_safe("results/cptac/cptac_brca_target_protein_immune_correlations.csv")
scrna <- read_csv_safe("results/scrna/tisch_brca_emtab8107_target_top_celltypes.csv")
med <- read_csv_safe("results/mediation_mr_results.csv")

clean_text <- function(x) {
  x <- gsub("<e2><80><93>", "-", x, fixed = TRUE)
  x <- gsub("<e2><80><94>", "-", x, fixed = TRUE)
  x <- gsub("<e2><86><92>", "to ", x, fixed = TRUE)
  x <- gsub("<e2><9c><93>", "yes", x, fixed = TRUE)
  x
}

extract_or <- function(x) {
  as.numeric(sub("^([0-9.]+).*", "\\1", x))
}

short_feature <- function(x) {
  x <- gsub("Immune_Checkpoint", "Checkpoint", x)
  x <- gsub("M1_Macrophage", "M1 macro", x)
  x <- gsub("M2_Macrophage", "M2 macro", x)
  x <- gsub("Macrophages M1", "M1 macro", x)
  x <- gsub("Macrophages M2", "M2 macro", x)
  x <- gsub("Macrophage", "Macro", x)
  x <- gsub("Macrophages", "Macro", x)
  x <- gsub("MicroenvironmentScore", "Microenv", x)
  x <- gsub("StromaScore", "Stroma", x)
  x <- gsub("ImmuneScore", "Immune", x)
  x <- gsub("CD8\\+ T-cells", "CD8 T", x)
  x <- gsub("T cells CD8", "CD8 T", x)
  x <- gsub("T cells regulatory \\(Tregs\\)", "Treg", x)
  x <- gsub("Dendritic cells activated", "Dendritic", x)
  x <- gsub("B-cells", "B cell", x)
  x <- gsub("B_cell", "B cell", x)
  x
}

protein_tier <- setNames(master$Tier, master$Protein)
row_labels <- paste0(priority, " (", protein_tier[priority], ")")
row_labels[is.na(protein_tier[priority])] <- priority[is.na(protein_tier[priority])]
protein_label <- setNames(row_labels, priority)

panel <- data.frame()
add_cell <- function(protein, domain, label, class, detail = "") {
  panel <<- rbind(panel, data.frame(
    protein = protein,
    domain = domain,
    label = label,
    class = class,
    detail = detail,
    stringsAsFactors = FALSE
  ))
}

domains <- c(
  "Discovery MR",
  "Coloc",
  "OpenGWAS\npQTL MR",
  "TCGA\nexpression",
  "TCGA\nimmune",
  "CPTAC\nprotein",
  "scRNA\nsource",
  "Metabolite\nmediation"
)

for (p in priority) {
  m <- master[master$Protein == p, , drop = FALSE]

  if (nrow(m)) {
    or <- extract_or(m$`MR OR (95% CI)`[1])
    mr_label <- if (is.na(or)) "MR\navailable" else if (or >= 1) {
      sprintf("Risk +\nOR %.3f", or)
    } else {
      sprintf("Protective\nOR %.3f", or)
    }
    add_cell(p, domains[1], mr_label, ifelse(!is.na(or) && or >= 1, "risk", "protective"))

    verdict <- clean_text(m$`Coloc verdict`[1])
    pph4 <- clean_text(m$`Coloc PPH4`[1])
    pph4_num <- suppressWarnings(as.numeric(sub(" .*", "", pph4)))
    if (grepl("STRONG", verdict, ignore.case = TRUE)) {
      lab <- if (!is.na(pph4_num)) sprintf("Strong\nPPH4 %.2f", pph4_num) else "Strong"
      cls <- "strong"
    } else if (grepl("MODERATE", verdict, ignore.case = TRUE)) {
      lab <- if (!is.na(pph4_num)) sprintf("Moderate\nPPH4 %.2f", pph4_num) else "Moderate"
      cls <- "moderate"
    } else if (grepl("DISTINCT", verdict, ignore.case = TRUE)) {
      lab <- "Distinct\nsignal"
      cls <- "negative"
    } else if (grepl("not tested", verdict, ignore.case = TRUE)) {
      lab <- "Not\ntested"
      cls <- "missing"
    } else {
      lab <- "Insufficient"
      cls <- "partial"
    }
    add_cell(p, domains[2], lab, cls)
  } else {
    add_cell(p, domains[1], "No MR", "missing")
    add_cell(p, domains[2], "No coloc", "missing")
  }

  ogw_p <- ogw[ogw$protein == p | grepl(p, ogw$protein), , drop = FALSE]
  if (p == "APOE") ogw_p <- ogw[grepl("^APOE", ogw$protein), , drop = FALSE]
  breast_sig <- nrow(ogw_p) && any(ogw_p$outcome == "Breast" &
    ogw_p$method %in% c("Inverse variance weighted", "Wald ratio") &
    ogw_p$fdr_within_outcome < 0.05, na.rm = TRUE)
  breast_tested <- nrow(ogw_p) && any(ogw_p$outcome == "Breast", na.rm = TRUE)
  exposure_text <- paste(ogw_inst$exposure, collapse = " ")
  inst_found <- grepl(p, exposure_text, ignore.case = TRUE) ||
    (p == "ITIH3" && grepl("trypsin inhibitor heavy chain H3", exposure_text, ignore.case = TRUE)) ||
    (p == "TNFRSF6B" && grepl("receptor superfamily member 6B", exposure_text, ignore.case = TRUE)) ||
    (p == "APOE" && grepl("Apolipoprotein E", exposure_text, ignore.case = TRUE)) ||
    (p == "IL34" && grepl("Interleukin-34", exposure_text, ignore.case = TRUE)) ||
    (p == "FGF5" && grepl("Fibroblast growth factor 5", exposure_text, ignore.case = TRUE))
  if (breast_sig) {
    add_cell(p, domains[3], "Replicated\nBRCA", "strong")
  } else if (breast_tested) {
    add_cell(p, domains[3], "Tested\nNS", "partial")
  } else if (inst_found) {
    add_cell(p, domains[3], "pQTL\nfound", "partial")
  } else {
    add_cell(p, domains[3], "Not\nfound", "missing")
  }

  te <- tcga_expr[tcga_expr$gene == p, , drop = FALSE]
  if (nrow(te)) {
    fc <- te$logFC_tumour_vs_normal[1]
    add_cell(
      p, domains[4],
      sprintf("%s\nlogFC %.2f", ifelse(fc >= 0, "Tumour up", "Tumour down"), fc),
      ifelse(fc >= 0, "expr_up", "expr_down")
    )
  } else {
    add_cell(p, domains[4], "No data", "missing")
  }

  ti <- tcga_imm[tcga_imm$gene == p & tcga_imm$fdr < 0.05, , drop = FALSE]
  if (nrow(ti)) {
    ti <- ti[order(-abs(ti$rho)), , drop = FALSE][1, ]
    add_cell(
      p, domains[5],
      sprintf("%s\nrho %.2f", short_feature(ti$signature), ti$rho),
      ifelse(ti$rho >= 0, "immune_pos", "immune_neg")
    )
  } else {
    add_cell(p, domains[5], "No sig.", "missing")
  }

  cp <- cptac[cptac$gene == p & cptac$fdr < 0.05, , drop = FALSE]
  if (nrow(cp)) {
    cp <- cp[order(cp$fdr), , drop = FALSE][1, ]
    add_cell(
      p, domains[6],
      sprintf("%s\nrho %.2f", short_feature(cp$immune_feature), cp$rho),
      ifelse(cp$rho >= 0, "protein_pos", "protein_neg")
    )
  } else if (p %in% unique(cptac$gene)) {
    add_cell(p, domains[6], "Protein\nNS", "partial")
  } else {
    add_cell(p, domains[6], "Not\ncovered", "missing")
  }

  sc <- scrna[scrna$gene == p & scrna$annotation_level == "major_lineage", , drop = FALSE]
  if (nrow(sc)) {
    source <- sc$top_celltype[1]
    cls <- if (grepl("Malignant", source, ignore.case = TRUE)) "cell_malignant" else
      if (grepl("Macro|Mono", source, ignore.case = TRUE)) "cell_myeloid" else "cell_stromal"
    add_cell(p, domains[7], source, cls)
  } else {
    add_cell(p, domains[7], "Not\ndetected", "missing")
  }

  md <- med[med$protein == p & med$direct_consistent == TRUE & med$p_indirect < 0.05, , drop = FALSE]
  if (nrow(md)) {
    md <- md[order(md$p_indirect), , drop = FALSE][1, ]
    add_cell(
      p, domains[8],
      sprintf("%s\n%.1f%%", md$metabolite, md$prop_med_pct),
      "mediation"
    )
  } else {
    add_cell(p, domains[8], "None", "missing")
  }
}

panel$protein <- factor(panel$protein, levels = rev(priority), labels = rev(protein_label[priority]))
panel$domain <- factor(panel$domain, levels = domains)
panel$class <- factor(panel$class, levels = c(
  "risk", "protective", "strong", "moderate", "partial", "negative", "expr_up",
  "expr_down", "immune_pos", "immune_neg", "protein_pos", "protein_neg",
  "cell_malignant", "cell_myeloid", "cell_stromal", "mediation", "missing"
))

palette <- c(
  risk = "#C75146",
  protective = "#3B78B4",
  strong = "#2E8B57",
  moderate = "#67A9CF",
  partial = "#D8A03D",
  negative = "#8A8F98",
  expr_up = "#D46A6A",
  expr_down = "#5B8CC0",
  immune_pos = "#8E5EA2",
  immune_neg = "#4FA39A",
  protein_pos = "#9966AA",
  protein_neg = "#58AFA5",
  cell_malignant = "#B85C38",
  cell_myeloid = "#7A6FAF",
  cell_stromal = "#8C9A45",
  mediation = "#C18F2D",
  missing = "#E7E8EA"
)

p <- ggplot(panel, aes(x = domain, y = protein, fill = class)) +
  geom_tile(color = "white", linewidth = 1.2, width = 0.98, height = 0.92) +
  geom_text(aes(label = label), size = 3.05, lineheight = 0.9, color = "#1F2328") +
  scale_fill_manual(values = palette, drop = FALSE, guide = "none") +
  scale_x_discrete(position = "top") +
  coord_equal(ratio = 0.72, clip = "off") +
  labs(
    title = "Integrated evidence map for priority circulating proteins",
    subtitle = "Compact summary of genetic, tumour, proteomic, immune, single-cell and mediation evidence",
    x = NULL,
    y = NULL,
    caption = "NS: not significant. BRCA: breast cancer. Only the strongest immune association per dataset is shown."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 10.5, color = "#555A60", margin = margin(b = 10)),
    plot.caption = element_text(size = 8.5, color = "#62676D", hjust = 0, margin = margin(t = 10)),
    axis.text.x = element_text(face = "bold", size = 9.5, color = "#2A2F35", margin = margin(b = 6)),
    axis.text.y = element_text(face = "bold", size = 10, color = "#2A2F35", margin = margin(r = 8)),
    panel.grid = element_blank(),
    plot.margin = margin(14, 18, 12, 14)
  )

ggsave(file.path(out_fig, "fig13_integrated_evidence_map.png"), p, width = 12.5, height = 5.8, dpi = 320, bg = "white")
ggsave(file.path(out_fig, "fig13_integrated_evidence_map.pdf"), p, width = 12.5, height = 5.8, bg = "white")
write.csv(panel, file.path(out_tab, "integrated_evidence_map_cells.csv"), row.names = FALSE)

cat("Saved integrated evidence map:\n")
cat(" - results/figures/fig13_integrated_evidence_map.png\n")
cat(" - results/figures/fig13_integrated_evidence_map.pdf\n")
cat(" - results/network/integrated_evidence_map_cells.csv\n")
