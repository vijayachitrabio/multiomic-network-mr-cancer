#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

out_dir <- "results/decode"
batch_file <- file.path(out_dir, "decode_batch_mr_results_v2.csv")
abo_file <- file.path(out_dir, "decode_abo_mr_results.csv")
master_file <- "results/tables/STable_master_evidence.csv"

batch <- fread(batch_file)
abo <- fread(abo_file)
master <- fread(master_file)

clean_protein <- function(x) {
  x <- gsub(" \\(deCODE\\)", "", x)
  x
}

res <- bind_rows(
  as.data.frame(batch),
  as.data.frame(abo)
) %>%
  mutate(
    protein = clean_protein(exposure),
    protein_key = ifelse(protein == "INHBB_Activin_B", "INHBB", protein),
    outcome_label = ifelse(grepl("Endometrial", outcome), "Endometrial", "Breast"),
    p_plot = ifelse(is.na(pval), NA_real_, pmin(-log10(pval), 14)),
    p_label = case_when(
      is.na(pval) ~ "not evaluable",
      pval < 1e-99 ~ "p < 1e-99",
      pval < 0.001 ~ paste0("p = ", scientific(pval, digits = 2)),
      TRUE ~ paste0("p = ", signif(pval, 2))
    )
  )

res <- res %>%
  left_join(
    master %>% select(protein_key = protein, discovery_b = mr_b, discovery_or = mr_or, tier_short),
    by = "protein_key"
  ) %>%
  mutate(
    concordant = case_when(
      is.na(b) | is.na(discovery_b) ~ NA,
      b * discovery_b > 0 ~ TRUE,
      TRUE ~ FALSE
    ),
    status = case_when(
      method == "Harmonization Failed" ~ "Not evaluable",
      is.na(pval) ~ "Not evaluable",
      concordant == FALSE ~ "Discordant",
      pval < 0.05 ~ "Concordant, p < 0.05",
      TRUE ~ "Concordant, not significant"
    ),
    status = factor(
      status,
      levels = c("Concordant, p < 0.05", "Concordant, not significant", "Discordant", "Not evaluable")
    ),
    direction = case_when(
      is.na(b) ~ "NA",
      b > 0 ~ "Risk increasing",
      b < 0 ~ "Protective",
      TRUE ~ "Null"
    ),
    direction_symbol = case_when(
      is.na(b) ~ "NA",
      b > 0 ~ "↑",
      b < 0 ~ "↓",
      TRUE ~ "–"
    ),
    method_label = case_when(
      protein == "INHBB_Activin_B" ~ "secondary assay",
      method == "Inverse variance weighted" ~ "2 SNP IVW",
      method == "Harmonization Failed" ~ "failed",
      TRUE ~ "1 SNP Wald"
    ),
    display_name = case_when(
      protein == "INHBB_Activin_B" ~ "INHBB (Activin B)",
      TRUE ~ protein
    )
  )

priority <- c(
  "ABO", "TNFRSF6B", "FGF5", "ITIH3", "APOE", "KLB",
  "FGFR4", "IL34", "INHBB", "INHBB (Activin B)", "SWAP70", "UMOD", "CGREF1"
)

res <- res %>%
  mutate(
    display_name = factor(display_name, levels = rev(priority)),
    x_text = 14.75
  )

status_cols <- c(
  "Concordant, p < 0.05" = "#1B7F5A",
  "Concordant, not significant" = "#7A8793",
  "Discordant" = "#B63E3E",
  "Not evaluable" = "#A7A7A7"
)

status_shapes <- c(
  "Concordant, p < 0.05" = 16,
  "Concordant, not significant" = 1,
  "Discordant" = 17,
  "Not evaluable" = 4
)

plot_df <- res %>% mutate(p_plot2 = ifelse(is.na(p_plot), 0, p_plot))

p <- ggplot(plot_df, aes(y = display_name)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "#6B7280", linewidth = 0.45) +
  geom_segment(
    aes(x = 0, xend = p_plot2, yend = display_name, color = status),
    linewidth = 1.4,
    alpha = 0.55,
    lineend = "round"
  ) +
  geom_point(aes(x = p_plot2, color = status, shape = status), size = 3.4, stroke = 1.1) +
  geom_text(
    aes(x = x_text, label = paste0(direction_symbol, "  ", p_label, "  |  ", method_label)),
    hjust = 0,
    size = 3.15,
    color = "#263645"
  ) +
  annotate(
    "text",
    x = -log10(0.05),
    y = length(priority) + 0.55,
    label = "p = 0.05",
    hjust = -0.05,
    vjust = 0,
    size = 3,
    color = "#4B5563"
  ) +
  scale_x_continuous(
    limits = c(0, 22),
    breaks = c(0, 1.3, 3, 6, 9, 12, 14),
    labels = c("0", "1.3", "3", "6", "9", "12", "≥14"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_manual(values = status_cols, drop = FALSE) +
  scale_shape_manual(values = status_shapes, drop = FALSE) +
  labs(
    title = "deCODE cross-platform pQTL sensitivity analysis of manuscript instruments",
    subtitle = "Bars show −log10(P); colors indicate concordance. OR magnitudes omitted because deCODE protein scales differ by assay.",
    x = expression(-log[10](italic(P))),
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle = element_text(size = 10.5, color = "#4B5563", hjust = 0, margin = margin(b = 10)),
    axis.text.y = element_text(face = "bold", color = "#263645", size = 10.5),
    axis.text.x = element_text(color = "#4B5563"),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 8)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#E5E7EB", linewidth = 0.4),
    legend.position = "bottom",
    legend.justification = "left",
    legend.text = element_text(size = 9.5),
    plot.margin = margin(12, 235, 12, 12)
  )

png_file <- file.path(out_dir, "decode_forest_plot.png")
pdf_file <- file.path(out_dir, "decode_forest_plot.pdf")
ggsave(png_file, plot = p, width = 10.5, height = 6.8, dpi = 600, bg = "white")
ggsave(pdf_file, plot = p, width = 10.5, height = 6.8, bg = "white")

cat("Saved improved deCODE validation plot:\n")
cat(" - ", png_file, "\n", sep = "")
cat(" - ", pdf_file, "\n", sep = "")
