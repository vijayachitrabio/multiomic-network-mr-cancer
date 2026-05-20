#!/usr/bin/env Rscript
## Script 50: OpenGWAS cross-platform replication figure
## Discovery (FinnGen Olink) vs Replication (OpenGWAS SomaScan)
## Output: results/figures/fig8_opengwas_replication.png/.pdf

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr)
})

proj    <- "."
out_dir <- file.path(proj, "results/figures")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# ── Discovery ─────────────────────────────────────────────────────────────────
disc_all <- fread(file.path(proj, "results/phase2_protein_cancer/protein_cancer_mr_results_significant.csv"))
disc <- disc_all |>
  filter(grepl("Inverse|Wald", method)) |>
  mutate(cancer_short = case_when(
    grepl("Breast",      outcome) ~ "Breast",
    grepl("Endometrial|GCST006464", outcome) ~ "Endometrial",
    grepl("Ovarian",     outcome) ~ "Ovarian",
    TRUE ~ outcome
  )) |>
  select(protein=exposure, cancer_short, b_disc=b, se_disc=se,
         or_disc=or, lci_disc=or_lci95, uci_disc=or_uci95, pval_disc=pval)

# ── OpenGWAS replication (IVW / Wald ratio only) ─────────────────────────────
rep_raw <- fread(file.path(proj, "results/opengwas/opengwas_5protein_replication_mr_results.csv"))
rep <- rep_raw |>
  filter(grepl("Inverse|Wald", method)) |>
  select(protein, cancer_short=outcome, b_rep=b, se_rep=se,
         or_rep=or, lci_rep=or_lci95, uci_rep=or_uci95, pval_rep=pval, nsnp_rep=nsnp)

# ── Join: only rows where BOTH discovery and replication exist ─────────────────
# clean APOE isoform names to match discovery
rep <- rep |> mutate(protein = sub("_E[0-9]+$", "", protein))

combined <- inner_join(disc, rep, by = c("protein", "cancer_short")) |>
  mutate(dir_conc = sign(b_disc) == sign(b_rep),
         protein_label = protein)

message(sprintf("Matched protein-cancer pairs: %d", nrow(combined)))
print(combined |> select(protein, cancer_short, b_disc, pval_disc, b_rep, pval_rep, dir_conc))

if (!nrow(combined)) {
  message("No overlapping discovery+replication pairs — check protein name matching")
  quit(save="no")
}

# ── Long format for side-by-side forest ───────────────────────────────────────
long <- bind_rows(
  combined |> transmute(protein=protein_label, cancer_short, dir_conc,
                         OR=or_disc, LCI=lci_disc, UCI=uci_disc, pval=pval_disc,
                         study="Discovery\n(FinnGen Olink, N=619)"),
  combined |> transmute(protein=protein_label, cancer_short, dir_conc,
                         OR=or_rep,  LCI=lci_rep,  UCI=uci_rep,  pval=pval_rep,
                         study="Replication\n(OpenGWAS SomaScan)")
) |>
  mutate(
    study   = factor(study, levels=c("Discovery\n(FinnGen Olink, N=619)",
                                     "Replication\n(OpenGWAS SomaScan)")),
    sig_lab = case_when(pval < 0.001 ~ "***", pval < 0.01 ~ "**", pval < 0.05 ~ "*", TRUE ~ ""),
    y_label = if_else(cancer_short=="Breast", protein,
                      paste0(protein, "\n[", cancer_short, "]"))
  )

col_disc <- "#2471A3"
col_rep  <- "#7D3C98"
study_colours <- c("Discovery\n(FinnGen Olink, N=619)" = col_disc,
                   "Replication\n(OpenGWAS SomaScan)"  = col_rep)

p <- ggplot(long, aes(x=OR, y=y_label, colour=study, shape=study)) +
  facet_wrap(~study, ncol=2, scales="free_x") +
  geom_vline(xintercept=1, linetype="dashed", colour="grey45", linewidth=0.4) +
  geom_errorbarh(aes(xmin=LCI, xmax=UCI), height=0.25, linewidth=0.75) +
  geom_point(size=4) +
  geom_text(aes(x=UCI + 0.002, label=sig_lab),
            hjust=0, size=4.5, show.legend=FALSE, colour="grey20") +
  scale_colour_manual(values=study_colours, name=NULL) +
  scale_shape_manual(values=c("Discovery\n(FinnGen Olink, N=619)"=16,
                               "Replication\n(OpenGWAS SomaScan)"=17), name=NULL) +
  scale_x_continuous(labels=function(x) sprintf("%.3f", x)) +
  labs(
    x       = "Odds ratio (95% CI) per SD increase in protein level",
    y       = NULL,
    title   = "Cross-platform replication: FinnGen Olink → OpenGWAS SomaScan",
    subtitle= paste0(
      nrow(combined), " protein–cancer pairs  ·  ",
      sum(combined$dir_conc & combined$pval_rep < 0.05), " directionally replicated (p<0.05)"
    ),
    caption = "Significance: *** p<0.001, ** p<0.01, * p<0.05.\nOpenGWAS instruments from Sun et al. 2018 SomaScan (INTERVAL study, N≈3,301)."
  ) +
  theme_bw(base_size=13) +
  theme(
    plot.title      = element_text(face="bold", size=13, colour="#1A5276"),
    plot.subtitle   = element_text(size=10, colour="grey35"),
    plot.caption    = element_text(size=8.5, colour="grey55"),
    axis.text.y     = element_text(face="italic", size=11),
    strip.background= element_rect(fill="grey95", colour="grey70"),
    strip.text      = element_text(face="bold", size=10.5),
    legend.position = "none",
    panel.grid.minor= element_blank(),
    panel.grid.major.y = element_line(colour="grey94")
  )

ggsave(file.path(out_dir, "fig8_opengwas_replication.png"),
       p, width=9, height=max(4, 1.2*nrow(combined)+2), dpi=300, bg="white")
ggsave(file.path(out_dir, "fig8_opengwas_replication.pdf"),
       p, width=9, height=max(4, 1.2*nrow(combined)+2))
message("✓ fig8_opengwas_replication.png/.pdf saved")
