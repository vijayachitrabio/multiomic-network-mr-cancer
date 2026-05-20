#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

out_fig <- "results/figures"
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

cells_path <- "results/sensitivity/priority_pleiotropy_robustness_map_cells.csv"
summary_path <- "results/sensitivity/priority_pleiotropy_summary.csv"

if (!file.exists(cells_path) || !file.exists(summary_path)) {
  stop("Run scripts/55_priority_pleiotropy_robustness_map.R first.")
}

cells <- read.csv(cells_path, stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.csv(summary_path, stringsAsFactors = FALSE, check.names = FALSE)

priority <- c("EFNA1", "TNFRSF6B", "FGF5", "ATRAID", "APOE", "ITIH3", "IL34")
domains <- c("Discovery", "Steiger", "Coloc", "OpenGWAS", "Mediation")
domain_labels <- c(
  Discovery = "Instrument",
  Steiger = "Direction",
  Coloc = "Coloc",
  OpenGWAS = "External pQTL",
  Mediation = "Mediation"
)

clean_domain <- function(x) {
  x <- gsub("\\n", " ", x)
  x <- sub("Discovery.*", "Discovery", x)
  x <- sub("MR-Egger.*", "MR-Egger", x)
  x <- sub("Steiger.*", "Steiger", x)
  x <- sub("Coloc.*|Colocalization.*", "Coloc", x)
  x <- sub("OpenGWAS.*", "OpenGWAS", x)
  x <- sub("Mediation.*", "Mediation", x)
  x
}

short_label <- function(domain, label) {
  x <- gsub("\n", " ", label)
  if (domain == "Discovery") {
    x <- gsub("1 SNP Wald", "1 SNP\nWald", x)
    x <- gsub("2 SNP IVW", "2 SNP\nIVW", x)
  } else if (domain == "Steiger") {
    x <- "forward"
  } else if (domain == "Coloc") {
    x <- gsub("Strong PPH4 ", "PPH4\n", x)
    x <- gsub("Moderate PPH4 ", "PPH4\n", x)
    x <- gsub("Distinct locus", "distinct\nlocus", x)
    x <- gsub("Insufficient", "insufficient", x)
  } else if (domain == "OpenGWAS") {
    x <- gsub("No pQTL found", "no\npQTL", x)
    x <- gsub("pQTL found BRCA NS", "pQTL\nNS", x)
    x <- gsub("Replicated OR ", "rep\nOR ", x)
    x <- gsub("Tested NS", "tested\nNS", x)
  } else if (domain == "Mediation") {
    x <- gsub("No robust path", "none", x)
    x <- gsub("Total_BCAA", "BCAA", x)
    x <- gsub("TG_by_PG", "TG/PG", x)
    x <- gsub(" ", "\n", x, fixed = TRUE)
  }
  x
}

cells$protein_clean <- sub(" \\(.*", "", cells$protein)
cells$domain_clean <- clean_domain(cells$domain)
cells <- cells[cells$protein_clean %in% priority & cells$domain_clean %in% domains, ]
cells$label_short <- mapply(short_label, cells$domain_clean, cells$label, USE.NAMES = FALSE)

summary <- summary[match(priority, summary$protein), ]
summary$row_label <- sprintf("%s   %s   OR %.3f", summary$protein, summary$tier, summary$discovery_or)
names(summary$row_label) <- summary$protein
summary$y <- rev(seq_along(priority))

cells$row_label <- summary$row_label[match(cells$protein_clean, summary$protein)]
cells$y <- summary$y[match(cells$protein_clean, summary$protein)]
cells$domain_clean <- factor(cells$domain_clean, levels = domains, labels = domain_labels[domains])
cells$class <- factor(cells$class, levels = c("strong", "partial", "limited", "warning", "missing"))

palette <- c(
  strong = "#2E7D5B",
  partial = "#D9A441",
  limited = "#86AFCF",
  warning = "#C26358",
  missing = "#E2E4E7"
)

row_stripes <- data.frame(
  row_label = summary$row_label,
  y = summary$y,
  stripe = rep(c(TRUE, FALSE), length.out = length(summary$row_label))
)

p <- ggplot(cells, aes(x = domain_clean, y = y)) +
  geom_tile(
    data = subset(row_stripes, stripe),
    aes(x = 3, y = y),
    inherit.aes = FALSE,
    width = Inf,
    height = 0.92,
    fill = "#F7F8F9",
    color = NA
  ) +
  geom_tile(
    aes(fill = class),
    width = 0.88,
    height = 0.72,
    color = "white",
    linewidth = 0.95
  ) +
  geom_text(
    aes(label = label_short),
    size = 3.15,
    lineheight = 0.88,
    color = "#1F252B"
  ) +
  scale_fill_manual(
    values = palette,
    breaks = c("strong", "partial", "limited", "warning", "missing"),
    labels = c("supports", "partial", "limited", "concern", "not available"),
    name = NULL
  ) +
  scale_x_discrete(position = "top", expand = expansion(add = 0.42)) +
  scale_y_continuous(
    breaks = summary$y,
    labels = summary$row_label,
    expand = expansion(add = 0.42)
  ) +
  labs(
    title = "Priority protein robustness profile",
    subtitle = "Traffic-light summary of instrument strength, directionality, colocalization, external pQTL replication, and mediation evidence",
    x = NULL,
    y = NULL,
    caption = paste(strwrap(
      "Rows show protein, evidence tier, and primary MR odds ratio. MR-Egger is not shown because these priority signals have fewer than three discovery instruments. Mediation cells are hypothesis-generating and do not imply definitive direct/indirect effect decomposition.",
      width = 150
    ), collapse = "\n")
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15.5, color = "#20252B", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 9.6, color = "#5A6067", margin = margin(b = 14)),
    plot.caption = element_text(size = 8.2, color = "#62676D", hjust = 0, margin = margin(t = 12), lineheight = 1.05),
    axis.text.x = element_text(face = "bold", size = 10.2, color = "#252A31", margin = margin(b = 8)),
    axis.text.y = element_text(face = "bold", size = 10.2, color = "#252A31", margin = margin(r = 12)),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.justification = "left",
    legend.text = element_text(size = 9.0),
    legend.key.width = unit(0.42, "cm"),
    legend.key.height = unit(0.32, "cm"),
    plot.margin = margin(16, 22, 16, 20)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave(
  file.path(out_fig, "fig15_priority_pleiotropy_chronograph.png"),
  p,
  width = 10.8,
  height = 5.4,
  dpi = 340,
  bg = "white"
)
ggsave(
  file.path(out_fig, "fig15_priority_pleiotropy_chronograph.pdf"),
  p,
  width = 10.8,
  height = 5.4,
  bg = "white"
)

cat("Saved redesigned traffic-light robustness figure:\n")
cat(" - results/figures/fig15_priority_pleiotropy_chronograph.png\n")
cat(" - results/figures/fig15_priority_pleiotropy_chronograph.pdf\n")
