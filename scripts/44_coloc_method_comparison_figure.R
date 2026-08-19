#!/usr/bin/env Rscript
## Script 45 (v3): coloc.abf vs coloc.SuSiE comparison figure
## Output: results/figures/fig_coloc_method_comparison.png/.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

proj    <- "."
dat_f   <- file.path(proj, "results/tables/STable8_protein_coloc.csv")
out_dir <- file.path(proj, "results/figures")

d <- fread(dat_f)
d[is.na(PPH4_abf),   PPH4_abf   := 0]
d[is.na(PPH4_susie), PPH4_susie := 0]

priority_proteins <- c(
  "ATRAID", "EFNA1", "TNFRSF6B", "APOE",
  "PM20D1", "IL34", "ITIH3", "SNX15"
)

d <- d[protein %in% priority_proteins]
d[, protein := factor(protein, levels = priority_proteins)]
d <- d[order(protein)]

d[, group := fcase(
  PPH4_susie >= 0.8 & PPH4_abf < 0.5,  "SuSiE identifies shared signal",
  PPH4_susie >= 0.8 & PPH4_abf >= 0.5, "Concordant support",
  PPH4_abf   >= 0.5 & PPH4_susie < 0.5,"ABF-supported",
  default = "Neither"
)]

colour_map <- c(
  "SuSiE identifies shared signal" = "#B4473E",
  "Concordant support"             = "#28607A",
  "ABF-supported"                  = "#D98A2B",
  "Neither"                        = "#9AA4A6"
)

shape_map <- c(
  "SuSiE identifies shared signal" = 17,
  "Concordant support"             = 16,
  "ABF-supported"                  = 15,
  "Neither"                        = 16
)

label_pos <- data.table(
  protein = factor(priority_proteins, levels = priority_proteins),
  label_x = c(0.060, 0.105, 0.880, 0.535, 0.055, 0.055, 0.055, 0.055),
  label_y = c(0.940, 0.900, 0.925, 0.040, 0.085, 0.045, 0.015, -0.015),
  hjust = c(0, 0, 0, 0, 0, 0, 0, 0)
)
d <- merge(d, label_pos, by = "protein", all.x = TRUE, sort = FALSE)

ref_lines <- data.frame(
  x = c(-0.04, -0.04, 0.8),
  y = c(-0.04, 0.8, -0.04),
  xend = c(1.02, 1.02, 0.8),
  yend = c(1.02, 0.8, 1.02),
  reference = c(
    "Equal posterior support",
    "PPH4 = 0.80 threshold",
    "PPH4 = 0.80 threshold"
  )
)

subtitle_txt <- sprintf(
  "%d priority proteins; breast cancer GWAS (N = 228,951); FinnGen pQTL (N = 619)",
  nrow(d)
)

p <- ggplot(d, aes(x = PPH4_abf, y = PPH4_susie, colour = group, shape = group)) +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0.8, ymax = 1,
           fill = "#B4473E", alpha = 0.055) +
  annotate("rect", xmin = 0.8, xmax = 1, ymin = 0.8, ymax = 1,
           fill = "#28607A", alpha = 0.055) +
  geom_segment(
    data = ref_lines,
    aes(x = x, y = y, xend = xend, yend = yend, linetype = reference),
    colour = "#8D8D8D",
    linewidth = 0.36,
    inherit.aes = FALSE
  ) +
  geom_point(size = 4.3, stroke = 0.15, alpha = 0.96) +
  geom_segment(
    aes(xend = label_x, yend = label_y),
    colour = "#A5A5A5",
    linewidth = 0.28,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = label_x, y = label_y, label = protein, hjust = hjust),
    fontface = "italic",
    size = 3.8,
    colour = "#222222",
    show.legend = FALSE
  ) +
  scale_colour_manual(values = colour_map, name = NULL) +
  scale_shape_manual(values = shape_map, name = NULL) +
  scale_linetype_manual(
    values = c("PPH4 = 0.80 threshold" = "dotted",
               "Equal posterior support" = "longdash"),
    name = NULL
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = 0.02)) +
  labs(
    x = "PPH4 from coloc.abf",
    y = "PPH4 from coloc.SuSiE",
    title = "Colocalization posterior agreement",
    subtitle = subtitle_txt
  ) +
  coord_cartesian(xlim = c(-0.04, 1.02), ylim = c(-0.04, 1.02), clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14.5, colour = "#222222"),
    plot.subtitle = element_text(size = 10.2, colour = "#4F4F4F",
                                 margin = margin(t = 3, b = 8)),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 9.2),
    legend.margin = margin(t = -2),
    panel.grid.major = element_line(colour = "#ECECEC", linewidth = 0.3),
    axis.line = element_line(colour = "#4C4C4C", linewidth = 0.35),
    axis.ticks = element_line(colour = "#4C4C4C", linewidth = 0.3),
    axis.title = element_text(size = 11.2, colour = "#222222"),
    axis.text = element_text(size = 9.6, colour = "#4C4C4C"),
    plot.margin = margin(12, 18, 10, 24)
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1),
    shape = guide_legend(order = 1, nrow = 1),
    linetype = guide_legend(order = 2, nrow = 1)
  )

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(out_dir, "fig_coloc_method_comparison.png"),
       p, width = 8.4, height = 7.2, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig_coloc_method_comparison.pdf"),
       p, width = 8.4, height = 7.2)

message("Saved fig_coloc_method_comparison.png/.pdf")
