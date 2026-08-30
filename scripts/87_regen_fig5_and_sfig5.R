#!/usr/bin/env Rscript
## Script 87: regenerate Figure 5 (mediation paths) and Supplementary Figure 5 (ER pattern)
## ─────────────────────────────────────────────────────────────────────────────
## Two audit items:
##
## 1. FIGURE 5 had two problems. Its subtitle said "consistent direction paths
##    only" while giving the statistically unsupported ATRAID -> TG_by_PG path
##    the same visual weight as supported paths, and it carried no confidence
##    intervals. Worse, it was drawn from `results/mediation_mr_results.csv`,
##    whose step-2 estimates are IVW-based -- so it showed ATRAID at 9.3% mediated
##    with p = 0.00097, contradicting main-text Table 2, which uses the weighted
##    median step-2 estimator (1.0% mediated, p = 0.627).
##    This script draws it from the authoritative source Table 2 itself uses,
##    `results/mediation/mediation_step2_sensitivity.csv` (weighted median rows),
##    adds 95% confidence intervals, and marks the unsupported path explicitly.
##
## 2. SUPPLEMENTARY FIGURE 5 was titled "ER subtype specificity of breast cancer
##    protein associations". "Specificity" overstates it: the categories are
##    nominal per-subtype significance, not formal effect-modification. Only four
##    proteins show FDR-supported heterogeneity. Retitled accordingly.
##
## Both are exported at 600 dpi at the dimensions already used in the submission
## package (Figure 5: 10 x 6 in; Supplementary Figure 5: 7 x 5 in).
##
## Outputs: results/figures/fig4_mediation_paths.{png,pdf}
##          results/figures/sfig5_er_pattern.{png,pdf}

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
proj <- "."
out_dir <- file.path(proj, "results/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 600, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h,
         bg = "white")
  message(sprintf("  wrote %s.png / .pdf at 600 dpi", name))
}

## ── FIGURE 5 ────────────────────────────────────────────────────────────────
message("Figure 5: mediation paths (weighted median, matching Table 2)")
med <- fread(file.path(proj, "results/mediation/mediation_step2_sensitivity.csv"))
wm <- med[grepl("median", step2_method, ignore.case = TRUE)]
stopifnot(nrow(wm) == 6)

wm[, path_label := paste0(protein, " -> ", metabolite)]
wm[, supported := sig_indirect %in% c(TRUE, "TRUE", "True")]
wm[, dir_label := fifelse(!supported, "Not supported (p >= 0.05)",
                   fifelse(b_indirect < 0, "Protective", "Risk"))]
wm[, lab := sprintf("%.1f%%\np = %s", prop_med_pct,
                    ifelse(p_indirect < 0.001,
                           formatC(p_indirect, format = "e", digits = 1),
                           formatC(p_indirect, format = "f", digits = 3)))]
setorder(wm, b_indirect)

fig5 <- ggplot(wm, aes(x = reorder(path_label, -abs(b_indirect)),
                       y = b_indirect * 100, fill = dir_label)) +
  geom_col(width = 0.6,
           aes(colour = dir_label, linetype = !supported),
           linewidth = 0.6, show.legend = c(fill = TRUE, colour = FALSE, linetype = FALSE)) +
  geom_errorbar(aes(ymin = ci95_lo * 100, ymax = ci95_hi * 100),
                width = 0.18, linewidth = 0.5, colour = "grey25") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = lab,
                vjust = ifelse(b_indirect < 0, 1.5, -0.6)),
            size = 2.9, lineheight = 0.95, colour = "grey20") +
  scale_fill_manual(name = NULL,
    values = c(Protective = "#4A90D9", Risk = "#E05A5A",
               `Not supported (p >= 0.05)` = "#BDC3C7")) +
  scale_colour_manual(values = c(Protective = "#2C6BA8", Risk = "#A83A3A",
                                 `Not supported (p >= 0.05)` = "#7B8794"),
                      guide = "none") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.18))) +
  labs(
    title    = "Protein -> metabolite -> breast cancer mediation: indirect effects",
    x = "Mediation path",
    y = "Indirect effect on ln(OR) x 100"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
    plot.title    = element_text(face = "bold", size = 11),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )
save_fig(fig5, "fig4_mediation_paths", w = 10, h = 6)

## ── SUPPLEMENTARY FIGURE 5 ──────────────────────────────────────────────────
message("Supplementary Figure 5: nominal ER-subtype association patterns")
er <- fread(file.path(proj, "results/tables/STable_ER_subtype_heterogeneity.csv"))
pat_col <- grep("pattern", names(er), value = TRUE)[1]
lut <- c(ER_neg_specific = "ER-negative only",
         ER_pos_specific = "ER-positive only",
         both_subtypes   = "Both subtypes",
         neither         = "Neither")
er[, pat := fifelse(get(pat_col) %in% names(lut), lut[get(pat_col)], as.character(get(pat_col)))]
counts <- er[, .(N = .N), by = pat]
n_fdr <- er[fdr_heterogeneity < 0.05, .N]
setorder(counts, -N)

sfig5 <- ggplot(counts, aes(x = reorder(pat, -N), y = N, fill = pat)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = N), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c(`ER-positive only` = "#D4868C", `Both subtypes` = "#9B59B6",
                               Neither = "#BDC3C7", `ER-negative only` = "#6A9CBF"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, max(counts$N) + 2), expand = c(0, 0)) +
  labs(
    title    = "Nominal ER-subtype association patterns among MR-prioritized proteins",
    subtitle = sprintf(paste0("16 breast cancer proteins at FDR < 0.05. Categories reflect nominal ",
                              "per-subtype significance (p < 0.05) only.\n",
                              "Exploratory ER-positive versus ER-negative heterogeneity is FDR-supported ",
                              "for %d of the 16 (UMOD, FGF5, ATRAID, INHBB);\n",
                              "these labels are not evidence of subtype specificity."), n_fdr),
    x = NULL, y = "Number of proteins"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11.5),
    plot.subtitle = element_text(size = 8.2, colour = "grey30", lineheight = 1.2),
    plot.margin   = margin(t = 8, r = 14, b = 8, l = 8),
    panel.grid.major.x = element_blank()
  )
save_fig(sfig5, "sfig5_er_pattern", w = 8.2, h = 5.2)

message("\nDone.")
