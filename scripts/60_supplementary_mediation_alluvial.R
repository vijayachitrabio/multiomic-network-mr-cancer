#!/usr/bin/env Rscript

# Supplementary mediation flow figure.
# Uses only the six tested two-step MR paths and keeps labels cautious.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

project_dir <- "."
fig_dir <- file.path(project_dir, "results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, name, w = 10.5, h = 6.2) {
  message("Writing figure files...")
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 320)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  message("Saved: ", file.path(fig_dir, paste0(name, ".png")))
  message("Saved: ", file.path(fig_dir, paste0(name, ".pdf")))
}

step2_file <- file.path(project_dir, "results", "mediation", "mediation_step2_sensitivity.csv")
if (!file.exists(step2_file)) {
  stop("Missing input file: ", step2_file)
}

dt <- fread(step2_file)

# Weighted median is the conservative step-2 sensitivity used in the manuscript
# for the metabolite -> breast cancer leg.
paths <- dt[step2_method == "Weighted median",
  .(protein, metabolite, cancer, p_indirect, prop_med_pct, coloc_supported, evidence_grade)
]

paths[, status := fifelse(
  protein == "ATRAID", "Not robust",
  fifelse(coloc_supported, "Coloc-supported", "Suggestive")
)]

paths[, status_detail := fifelse(
  protein == "ATRAID", "Not robust in sensitivity analysis",
  fifelse(coloc_supported, "Shared-signal support", "MR evidence only")
)]

paths[, line_width := fifelse(status == "Not robust", 0.7, pmax(0.9, sqrt(prop_med_pct) * 0.55))]
paths[, line_alpha := fifelse(status == "Not robust", 0.45, 0.88)]

protein_order <- c("IL34", "EFNA1", "TNFRSF6B", "APOE", "ITIH3", "ATRAID")
paths[, protein := factor(protein, levels = protein_order)]
setorder(paths, protein)

# Fixed positions: four reading columns from left to right.
node_pos <- rbindlist(list(
  data.table(stage = "Protein", label = as.character(paths$protein), x = 1, y = seq(6, 1, by = -1)),
  data.table(stage = "Metabolite", label = c("Total_BCAA", "Gly", "TG_by_PG"),
             x = 2, y = c(5, 3, 1)),
  data.table(stage = "Outcome", label = "Breast cancer", x = 3, y = 3.5),
  data.table(stage = "Interpretation",
             label = c("Shared-signal support", "MR evidence only", "Not robust in sensitivity analysis"),
             x = 4, y = c(5.4, 3.4, 1.4))
), fill = TRUE)

node_pos[, display := fifelse(
  label == "Total_BCAA", "Total BCAA",
  fifelse(label == "TG_by_PG", "Triglycerides in phosphoglycerides", label)
)]

y_lookup <- setNames(node_pos$y, node_pos$label)
paths[, y_protein := unname(y_lookup[as.character(protein)])]
paths[, y_metabolite := unname(y_lookup[metabolite])]
paths[, y_outcome := unname(y_lookup["Breast cancer"])]
paths[, y_status := unname(y_lookup[status_detail])]

segments <- rbindlist(list(
  paths[, .(path_id = paste0(protein, "_protein_metabolite"),
            x = 1.18, xend = 1.82, y = y_protein, yend = y_metabolite,
            status, line_width, line_alpha)],
  paths[, .(path_id = paste0(protein, "_metabolite_outcome"),
            x = 2.18, xend = 2.82, y = y_metabolite, yend = y_outcome,
            status, line_width, line_alpha)],
  paths[, .(path_id = paste0(protein, "_outcome_status"),
            x = 3.18, xend = 3.82, y = y_outcome, yend = y_status,
            status, line_width, line_alpha)]
))

curve_points <- segments[, {
  t <- seq(0, 1, length.out = 60)
  ease <- 3 * t^2 - 2 * t^3
  .(x = x + (xend - x) * t,
    y = y + (yend - y) * ease)
}, by = .(path_id, status, line_width, line_alpha)]

status_cols <- c(
  "Coloc-supported" = "#2F7D73",
  "Suggestive" = "#5C6F91",
  "Not robust" = "#B95F54"
)

node_fill <- c(
  "Protein" = "#F7F2E8",
  "Metabolite" = "#EAF3F0",
  "Outcome" = "#F2EEF6",
  "Interpretation" = "#EEF1F5"
)

column_headers <- data.table(
  x = c(1, 2, 3, 4),
  y = 6.65,
  label = c("Protein", "Metabolite", "Outcome", "Interpretation")
)

p <- ggplot() +
  geom_path(
    data = curve_points,
    aes(x = x, y = y, group = path_id,
        colour = status, linewidth = line_width, alpha = line_alpha),
    lineend = "round"
  ) +
  geom_rect(
    data = node_pos,
    aes(xmin = x - 0.20, xmax = x + 0.20, ymin = y - 0.24, ymax = y + 0.24, fill = stage),
    colour = "#2C2C2C",
    linewidth = 0.35
  ) +
  geom_text(
    data = node_pos,
    aes(x = x, y = y, label = display),
    size = 3.05,
    lineheight = 0.9,
    colour = "#202020"
  ) +
  geom_text(
    data = column_headers,
    aes(x = x, y = y, label = label),
    size = 3.7,
    fontface = "bold",
    colour = "#202020"
  ) +
  geom_text(
    data = paths,
    aes(x = 2.52, y = y_metabolite + 0.28,
        label = paste0(round(prop_med_pct, 1), "%")),
    size = 2.65,
    colour = "#4A4A4A"
  ) +
  scale_colour_manual(values = status_cols, name = NULL) +
  scale_fill_manual(values = node_fill, guide = "none") +
  scale_linewidth_identity() +
  scale_alpha_identity() +
  coord_cartesian(xlim = c(0.65, 4.35), ylim = c(0.45, 6.95), clip = "off") +
  labs(
    title = "Two-step MR mediation paths for breast cancer",
    subtitle = "Flow width reflects estimated proportion mediated; interpretation reflects sensitivity and shared-signal evidence"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.02, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10.5, colour = "#4A4A4A", hjust = 0.02, margin = margin(b = 12)),
    plot.margin = margin(16, 22, 14, 22),
    legend.position = "bottom",
    legend.text = element_text(size = 9.5),
    legend.key.width = unit(1.2, "cm")
  )

save_fig(p, "sfig_mediation_alluvial", w = 10.5, h = 6.2)
