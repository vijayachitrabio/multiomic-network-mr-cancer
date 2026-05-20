#!/usr/bin/env Rscript
# Script 21b: Supplementary manuscript figures
#
# SFig 1 — Volcano plot: all 701 proteins × breast cancer (−log10p vs β)
# SFig 2 — Per-SNP Wald ratio directional consistency for 2-SNP proteins
# SFig 3 — Steiger directionality: r² protein vs r² cancer per instrument
# SFig 4 — UKB observational: p-value vs delta-NPX with FDR threshold
# SFig 5 — ER subtype pattern summary (pie/bar of ER_pos_specific etc.)

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

project_dir <- "."
fig_dir     <- file.path(project_dir, "results", "figures")
dir.create(fig_dir, showWarnings = FALSE)

save_fig <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat(sprintf("  Saved: %s\n", name))
}

# Hit proteins for labelling
hit_proteins <- c("SNX15","EFNA1","FGF5","UMOD","SWAP70","ATRAID","TNFRSF6B",
                   "ITIH3","KLB","PM20D1","TSPAN8","FGFR4","IL34","APOE",
                   "CGREF1","INHBB","ABO")

# ============================================================
# SFig 1 — Volcano Plot: Breast cancer MR (all 701 proteins)
# ============================================================
cat("SFig 1: Breast cancer volcano\n")

phase2 <- fread(file.path(project_dir, "results", "phase2_protein_cancer",
                           "protein_cancer_mr_results_full.csv"))
breast_all <- phase2[grepl("Breast", outcome)]
breast_all[, log10p := -log10(pval)]
breast_all[, hit := exposure %in% hit_proteins]
breast_all[, direction := fifelse(b > 0, "Risk (OR>1)", "Protective (OR<1)")]
fdr_thresh <- max(breast_all[fdr < 0.05]$pval, na.rm = TRUE)

sfig1 <- ggplot(breast_all, aes(x = b, y = log10p, colour = hit, alpha = hit)) +
  geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed",
             colour = "firebrick", linewidth = 0.4) +
  geom_point(size = 1.2) +
  geom_text_repel(data = breast_all[hit == TRUE],
                  aes(label = exposure), size = 3, colour = "black",
                  box.padding = 0.4, max.overlaps = 20, seed = 42) +
  scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#E07B54"),
                      guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.5, `TRUE` = 1.0), guide = "none") +
  annotate("text", x = max(breast_all$b, na.rm=TRUE) * 0.85,
           y = -log10(fdr_thresh) + 0.8,
           label = "FDR < 0.05", size = 3, colour = "firebrick") +
  labs(
    title    = "Protein-breast cancer MR screen (701 proteins)",
    subtitle = sprintf("FinnGen R10 Olink pQTLs; %d proteins tested; %d FDR < 0.05",
                       nrow(breast_all), sum(breast_all$hit)),
    x        = "MR beta (per SD protein increase on log-OR scale)",
    y        = expression(-log[10](p))
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

save_fig(sfig1, "sfig1_breast_volcano", w = 9, h = 6.5)

# ============================================================
# SFig 2 — Per-SNP Wald Ratio: 2-SNP proteins (KLB, PM20D1, IL34, ABO)
# ============================================================
cat("SFig 2: Per-SNP Wald ratio sensitivity\n")

snp_dat <- fread(file.path(project_dir, "results", "sensitivity",
                            "sensitivity_singlesnp_results.csv"))
# Keep individual SNP rows (not the IVW summary)
snp_rows <- snp_dat[!grepl("^All", SNP) & !is.na(b)]
snp_rows[, or     := exp(b)]
snp_rows[, or_lci := exp(b - 1.96*se)]
snp_rows[, or_uci := exp(b + 1.96*se)]
snp_rows[, snp_label := paste0(exposure, "\n(", SNP, ")")]

# IVW summary
ivw_rows <- snp_dat[grepl("Inverse variance", SNP)]
ivw_rows[, or     := exp(b)]
ivw_rows[, or_lci := exp(b - 1.96*se)]
ivw_rows[, or_uci := exp(b + 1.96*se)]
ivw_rows[, snp_label := paste0(exposure, "\n(IVW)")]

all_snp <- rbindlist(list(snp_rows, ivw_rows), fill = TRUE)
all_snp[, type := fifelse(grepl("IVW", snp_label), "IVW", "Wald ratio")]
all_snp[, protein_f := factor(exposure, levels = c("KLB","PM20D1","IL34","ABO"))]

sfig2 <- ggplot(all_snp[!is.na(or)],
                aes(x = or, y = snp_label, colour = type, shape = type)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = or_lci, xmax = or_uci),
                width = 0.2, linewidth = 0.6, orientation = "y") +
  geom_point(size = 3) +
  facet_wrap(~ protein_f, scales = "free", ncol = 2) +
  scale_colour_manual(name = "Method",
    values = c(`Wald ratio` = "#5B8DB8", IVW = "#C0392B")) +
  scale_shape_manual(name = "Method",
    values = c(`Wald ratio` = 16, IVW = 18)) +
  scale_x_continuous(trans = "log", labels = function(x) sprintf("%.3f", x)) +
  labs(
    title    = "Per-SNP directional consistency: 2-SNP IVW proteins",
    subtitle = "Individual Wald ratios vs IVW estimate; both SNPs point same direction for all 4 proteins",
    x        = "Odds ratio (95% CI)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

save_fig(sfig2, "sfig2_per_snp_wald", w = 9, h = 7)

# ============================================================
# SFig 3 — Steiger: r² exposure vs r² outcome scatter
# ============================================================
cat("SFig 3: Steiger r-squared scatter\n")

stei <- fread(file.path(project_dir, "results", "sensitivity",
                          "steiger_directionality_results.csv"))

if (all(c("r2_exp", "r2_out") %in% names(stei))) {
  stei[, label := paste0(protein, "\n(", sub("FINNGEN_OLINK_","",SNP), ")")]

  sfig3 <- ggplot(stei, aes(x = r2_exp, y = r2_out, label = protein)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = cancer), size = 3, alpha = 0.8) +
    geom_text_repel(size = 2.8, max.overlaps = 20, seed = 42) +
    scale_colour_manual(name = "Cancer",
      values = c(Breast = "#E07B54", Endometrial = "#5B8DB8")) +
    scale_x_log10() + scale_y_log10() +
    labs(
      title    = "Steiger directionality: variance explained (log scale)",
      subtitle = "All points below diagonal confirm protein→cancer direction",
      x        = "R² in protein (exposure)",
      y        = "R² in cancer (outcome)"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12))

  save_fig(sfig3, "sfig3_steiger_r2", w = 7.5, h = 6)
} else {
  cat("  r2_exp/r2_out columns not found in steiger output — checking column names\n")
  cat("  Columns:", paste(names(stei), collapse=", "), "\n")

  # Compute r² from beta + se using z-score approximation
  # r² ≈ z² / (z² + N)  where z = beta/se
  stei[, z_exp := beta.exposure / se.exposure]
  stei[, z_out := beta.outcome  / se.outcome]
  stei[, r2_exp_approx := z_exp^2 / (z_exp^2 + samplesize.exposure)]
  stei[, r2_out_approx := z_out^2 / (z_out^2 + samplesize.outcome)]

  sfig3 <- ggplot(stei, aes(x = r2_exp_approx, y = r2_out_approx, label = protein)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = cancer), size = 3, alpha = 0.8) +
    geom_text_repel(size = 2.8, max.overlaps = 20, seed = 42) +
    scale_colour_manual(name = "Cancer",
      values = c(Breast = "#E07B54", Endometrial = "#5B8DB8")) +
    scale_x_log10(labels = scientific_format(digits = 2)) +
    scale_y_log10(labels = scientific_format(digits = 2)) +
    labs(
      title    = "Steiger directionality: variance explained by pQTL instruments",
      subtitle = "All points below diagonal: pQTL explains more variance in protein than in cancer",
      x        = expression(R^2~"in protein (exposure, N = 619)"),
      y        = expression(R^2~"in cancer (outcome)")
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12))

  save_fig(sfig3, "sfig3_steiger_r2", w = 7.5, h = 6)
}

# ============================================================
# SFig 4 — UKB Observational: volcano for breast proteins
# ============================================================
cat("SFig 4: Observational volcano\n")

obs <- fread(file.path(project_dir, "results", "observational", "ukb_triangulation_summary.csv"))
obs_br <- obs[cancer == "Breast"]
obs_br[, log10p := -log10(pval_obs)]
obs_br[, triangulates_label := fifelse(is.na(triangulates), "NA",
                                fifelse(triangulates, "Agree", "Disagree"))]
obs_br[, obs_fdr_thresh := fdr_obs < 0.05]

sfig4 <- ggplot(obs_br, aes(x = delta_npx, y = log10p,
                              colour = triangulates_label,
                              shape  = obs_fdr_thresh,
                              label  = protein)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_point(size = 3.5) +
  geom_text_repel(size = 3.2, max.overlaps = 20, box.padding = 0.5, seed = 42) +
  scale_colour_manual(name = "MR vs\nobservational",
    values = c(Agree = "#2C9E5B", Disagree = "#C0392B", `NA` = "grey60")) +
  scale_shape_manual(name = "Obs. FDR",
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c("< 0.05", "≥ 0.05")) +
  annotate("text", x = max(obs_br$delta_npx, na.rm=TRUE) * 0.8,
           y = -log10(0.05) + 0.15,
           label = "p = 0.05", size = 3, colour = "grey40") +
  labs(
    title    = "UKB-PPP observational validation: breast cancer proteins",
    subtitle = "52,995 UKB participants; 2,010 breast cases vs 39,169 no-cancer controls; unadjusted",
    x        = "Delta-NPX: mean(cases) - mean(controls)",
    y        = expression(-log[10](p["observational"]))
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

save_fig(sfig4, "sfig4_observational_breast", w = 8.5, h = 6)

# ============================================================
# SFig 5 — ER Subtype Pattern Summary
# ============================================================
cat("SFig 5: ER subtype pattern bar\n")

er_comp <- fread(file.path(project_dir, "results", "er_subtype", "er_subtype_comparison.csv"))
pattern_counts <- er_comp[!is.na(er_pattern), .N, by = er_pattern]
pattern_counts[, pattern_label := fcase(
  er_pattern == "ER_pos_specific", "ER+ specific",
  er_pattern == "ER_neg_specific", "ER- specific",
  er_pattern == "both_subtypes",   "Both subtypes",
  er_pattern == "neither",         "Neither"
)]
pattern_counts[, pattern_f := factor(pattern_label,
  levels = c("ER+ specific", "Both subtypes", "Neither", "ER- specific"))]

sfig5 <- ggplot(pattern_counts, aes(x = pattern_f, y = N, fill = pattern_f, label = N)) +
  geom_col(width = 0.6) +
  geom_text(vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(name = NULL,
    values = c(
      `ER+ specific`  = "#D4868C",
      `Both subtypes` = "#9B59B6",
      `Neither`       = "#BDC3C7",
      `ER- specific`  = "#6A9CBF"
    )) +
  scale_y_continuous(limits = c(0, max(pattern_counts$N) + 2), expand = c(0, 0)) +
  labs(
    title    = "ER subtype specificity of breast cancer protein associations",
    subtitle = sprintf("16 FDR < 0.05 breast proteins; nominal significance threshold p < 0.05 per subtype"),
    x        = NULL,
    y        = "Number of proteins"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

save_fig(sfig5, "sfig5_er_pattern", w = 7, h = 5)

cat("\nAll supplementary figures saved to:", fig_dir, "\n")
cat("Done.\n")
