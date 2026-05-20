#!/usr/bin/env Rscript
## Script 46: ARIC replication forest plot (discovery vs ARIC side-by-side)
## Output: results/figures/fig7_aric_replication.png/.pdf

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr)
})

proj    <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
out_dir <- file.path(proj, "results/figures")

## ── Load data ────────────────────────────────────────────────────────────────
# Discovery (FinnGen Olink)
disc <- fread(file.path(proj, "results/phase2_protein_cancer/protein_cancer_mr_results_significant.csv")) |>
  select(protein = exposure, or_disc = or, lci_disc = or_lci95, uci_disc = or_uci95,
         pval_disc = pval, cancer = outcome) |>
  mutate(cancer_short = case_when(
    grepl("Breast", cancer) ~ "Breast",
    grepl("Endometrial|GCST006464", cancer) ~ "Endometrial",
    TRUE ~ cancer
  ))

# ARIC best-method — use ea_cis as primary; fall back to aa_cis
aric <- fread(file.path(proj, "results/replication/aric_replication_best_method_summary.csv")) |>
  filter(aric_analysis %in% c("ea_cis", "aa_cis")) |>
  mutate(analysis_rank = if_else(aric_analysis == "ea_cis", 1L, 2L)) |>
  group_by(exposure, outcome) |>
  slice_min(analysis_rank) |>
  ungroup() |>
  mutate(lci_aric = exp(b - 1.96 * se), uci_aric = exp(b + 1.96 * se)) |>
  select(protein = exposure, cancer_short_aric = outcome,
         or_aric = or, lci_aric, uci_aric,
         pval_aric = pval, dir_conc = direction_concordant) |>
  mutate(cancer_short_aric = case_when(
    grepl("Breast", cancer_short_aric) ~ "Breast",
    grepl("Endometrial|GCST006464", cancer_short_aric) ~ "Endometrial",
    TRUE ~ cancer_short_aric
  ))

# Join
combined <- inner_join(disc, aric, by = c("protein", "cancer_short" = "cancer_short_aric")) |>
  filter(dir_conc == TRUE) |>
  mutate(protein = factor(protein, levels = rev(sort(unique(protein)))))

# Tidy to long format for two-panel forest
long <- bind_rows(
  combined |> transmute(protein, cancer_short, OR = or_disc,  LCI = lci_disc,  UCI = uci_disc,
                         pval = pval_disc, study = "Discovery\n(FinnGen Olink, N=619)", dir_conc),
  combined |> transmute(protein, cancer_short, OR = or_aric,  LCI = lci_aric,  UCI = uci_aric,
                         pval = pval_aric, study = "Replication\n(ARIC SomaScan, N≈7,000)", dir_conc)
) |>
  mutate(
    study   = factor(study, levels = c("Discovery\n(FinnGen Olink, N=619)",
                                       "Replication\n(ARIC SomaScan, N≈7,000)")),
    sig_lab = case_when(pval < 0.001 ~ "***", pval < 0.01 ~ "**", pval < 0.05 ~ "*", TRUE ~ ""),
    cancer_label = if_else(cancer_short == "Breast", "", paste0(" [", cancer_short, "]"))
  )

col_disc  <- "#2471A3"
col_aric  <- "#1E8449"
study_colours <- c("Discovery\n(FinnGen Olink, N=619)"       = col_disc,
                   "Replication\n(ARIC SomaScan, N≈7,000)" = col_aric)

p <- ggplot(long, aes(x = OR, y = protein, colour = study, shape = study)) +
  facet_wrap(~ study, ncol = 2, scales = "free_x") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.28, linewidth = 0.7) +
  geom_point(size = 3.6) +
  geom_text(aes(x = UCI + 0.003, label = sig_lab),
            hjust = 0, size = 4.5, show.legend = FALSE, colour = "grey20") +
  scale_colour_manual(values = study_colours, name = NULL) +
  scale_shape_manual(values = c("Discovery\n(FinnGen Olink, N=619)" = 16,
                                "Replication\n(ARIC SomaScan, N≈7,000)" = 17), name = NULL) +
  scale_x_continuous(labels = function(x) sprintf("%.3f", x)) +
  labs(
    x       = "Odds ratio (95% CI) per SD increase in protein level",
    y       = NULL,
    title   = "Independent replication of MR-prioritised breast cancer protein associations",
    subtitle = paste0(
      "4 proteins directionally replicated in ARIC EA SomaScan  ·  ",
      "All discovery instruments cis-pQTLs (F > 30)"
    ),
    caption = "Significance: *** p<0.001, ** p<0.01, * p<0.05.  ARIC = Atherosclerosis Risk in Communities study."
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 13, colour = "#1A5276"),
    plot.subtitle   = element_text(size = 10, colour = "grey35"),
    plot.caption    = element_text(size = 8.5, colour = "grey55"),
    axis.text.y     = element_text(face = "italic", size = 11),
    strip.background= element_rect(fill = "grey95", colour = "grey70"),
    strip.text      = element_text(face = "bold", size = 10.5),
    legend.position = "none",
    panel.grid.minor= element_blank(),
    panel.grid.major.y = element_line(colour = "grey94")
  )

ggsave(file.path(out_dir, "fig7_aric_replication.png"),
       p, width = 9, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig7_aric_replication.pdf"),
       p, width = 9, height = 5.5)

message("✓ fig7_aric_replication.png/.pdf saved")
