#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

out_fig <- "results/figures"
out_tab <- "results/sensitivity"
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)

priority <- c("EFNA1", "ATRAID", "TNFRSF6B", "APOE", "ITIH3", "IL34", "FGF5")

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(data.frame())
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

master <- read_csv_safe("results/tables/STable7_master_evidence_compact.csv")
steiger <- read_csv_safe("results/sensitivity/steiger_comparison_table.csv")
egger <- read_csv_safe("results/sensitivity/sensitivity_egger_intercepts.csv")
ogw <- read_csv_safe("results/opengwas/opengwas_5protein_replication_mr_results.csv")
ogw_inst <- read_csv_safe("results/opengwas/opengwas_7protein_instruments_p5e8.csv")
med_sens <- read_csv_safe("results/mediation/mediation_step2_sensitivity.csv")

clean_text <- function(x) {
  x <- gsub("<e2><80><93>", "-", x, fixed = TRUE)
  x <- gsub("<e2><80><94>", "-", x, fixed = TRUE)
  x <- gsub("<e2><86><92>", "to ", x, fixed = TRUE)
  x <- gsub("<e2><9c><93>", "yes", x, fixed = TRUE)
  x
}

extract_or <- function(x) suppressWarnings(as.numeric(sub("^([0-9.]+).*", "\\1", x)))

instrument_count <- function(p) {
  x <- master[master$Protein == p, "MR SNPs"]
  if (!length(x)) return(NA_integer_)
  suppressWarnings(as.integer(x[1]))
}

has_ogw_instrument <- function(p) {
  exposure_text <- paste(ogw_inst$exposure, collapse = " ")
  grepl(p, exposure_text, ignore.case = TRUE) ||
    (p == "ITIH3" && grepl("trypsin inhibitor heavy chain H3", exposure_text, ignore.case = TRUE)) ||
    (p == "TNFRSF6B" && grepl("receptor superfamily member 6B", exposure_text, ignore.case = TRUE)) ||
    (p == "APOE" && grepl("Apolipoprotein E", exposure_text, ignore.case = TRUE)) ||
    (p == "IL34" && grepl("Interleukin-34", exposure_text, ignore.case = TRUE)) ||
    (p == "FGF5" && grepl("Fibroblast growth factor 5", exposure_text, ignore.case = TRUE))
}

cells <- data.frame()
add_cell <- function(protein, domain, label, class, note = "") {
  cells <<- rbind(cells, data.frame(
    protein = protein,
    domain = domain,
    label = label,
    class = class,
    note = note,
    stringsAsFactors = FALSE
  ))
}

domains <- c(
  "Discovery",
  "MR-Egger",
  "Steiger",
  "Coloc",
  "OpenGWAS",
  "Mediation"
)

summary_rows <- data.frame()

for (p in priority) {
  m <- master[master$Protein == p, , drop = FALSE]
  nsnp <- instrument_count(p)
  or <- if (nrow(m)) extract_or(m$`MR OR (95% CI)`[1]) else NA_real_
  tier <- if (nrow(m)) m$Tier[1] else NA_character_

  if (is.na(nsnp)) {
    add_cell(p, domains[1], "No\nMR", "missing")
  } else if (nsnp == 1) {
    add_cell(p, domains[1], "1 SNP\nWald", "limited")
  } else if (nsnp == 2) {
    add_cell(p, domains[1], "2 SNP\nIVW", "partial")
  } else {
    add_cell(p, domains[1], paste0(nsnp, " SNP\nMR"), "strong")
  }

  eg <- egger[egger$exposure == p, , drop = FALSE]
  if (!is.na(nsnp) && nsnp < 3) {
    add_cell(p, domains[2], "Not testable\n<3 SNP", "limited")
  } else if (nrow(eg) && !is.na(eg$pval[1])) {
    add_cell(
      p, domains[2],
      sprintf("Intercept\np=%.2g", eg$pval[1]),
      ifelse(eg$pval[1] < 0.05, "warning", "strong")
    )
  } else {
    add_cell(p, domains[2], "Not\navailable", "missing")
  }

  st <- steiger[steiger$protein == p & steiger$cancer == "Breast", , drop = FALSE]
  if (nrow(st)) {
    add_cell(
      p, domains[3],
      ifelse(isTRUE(st$direction_consistent[1]), "Forward\ndirection", "Reversed\nvariant"),
      ifelse(isTRUE(st$direction_consistent[1]), "strong", "warning"),
      sprintf("%s reversed SNPs", st$n_snp_rev[1])
    )
  } else {
    add_cell(p, domains[3], "Not\navailable", "missing")
  }

  if (nrow(m)) {
    verdict <- clean_text(m$`Coloc verdict`[1])
    pph4 <- clean_text(m$`Coloc PPH4`[1])
    pph4_num <- suppressWarnings(as.numeric(sub(" .*", "", pph4)))
    if (grepl("STRONG", verdict, ignore.case = TRUE)) {
      lab <- if (!is.na(pph4_num)) sprintf("Strong\nPPH4 %.2f", pph4_num) else "Strong"
      cls <- "strong"
    } else if (grepl("MODERATE", verdict, ignore.case = TRUE)) {
      lab <- if (!is.na(pph4_num)) sprintf("Moderate\nPPH4 %.2f", pph4_num) else "Moderate"
      cls <- "partial"
    } else if (grepl("DISTINCT", verdict, ignore.case = TRUE)) {
      lab <- "Distinct\nlocus"
      cls <- "warning"
    } else if (grepl("not tested", verdict, ignore.case = TRUE)) {
      lab <- "Not\ntested"
      cls <- "missing"
    } else {
      lab <- "Insufficient"
      cls <- "limited"
    }
    add_cell(p, domains[4], lab, cls)
  } else {
    add_cell(p, domains[4], "No data", "missing")
  }

  ogw_p <- ogw[ogw$protein == p | grepl(p, ogw$protein), , drop = FALSE]
  if (p == "APOE") ogw_p <- ogw[grepl("^APOE", ogw$protein), , drop = FALSE]
  ogw_sig <- nrow(ogw_p) && any(ogw_p$outcome == "Breast" &
    ogw_p$method %in% c("Inverse variance weighted", "Wald ratio") &
    ogw_p$fdr_within_outcome < 0.05, na.rm = TRUE)
  ogw_tested <- nrow(ogw_p) && any(ogw_p$outcome == "Breast", na.rm = TRUE)
  if (ogw_sig) {
    best <- ogw_p[ogw_p$outcome == "Breast" &
      ogw_p$method %in% c("Inverse variance weighted", "Wald ratio"), , drop = FALSE][1, ]
    add_cell(p, domains[5], sprintf("Replicated\nOR %.3f", best$or), "strong")
  } else if (ogw_tested) {
    add_cell(p, domains[5], "Tested\nNS", "partial")
  } else if (has_ogw_instrument(p)) {
    add_cell(p, domains[5], "pQTL found\nBRCA NS", "partial")
  } else {
    add_cell(p, domains[5], "No pQTL\nfound", "missing")
  }

  ms <- med_sens[med_sens$protein == p & med_sens$sig_indirect == TRUE, , drop = FALSE]
  if (nrow(ms)) {
    wm <- ms[grepl("Weighted median", ms$step2_method, ignore.case = TRUE), , drop = FALSE]
    supported <- if (nrow(wm)) wm[order(wm$p_indirect), , drop = FALSE][1, ] else ms[order(ms$p_indirect), , drop = FALSE][1, ]
    cls <- if (grepl("Not supported", supported$evidence_grade)) "warning" else
      if (grepl("Moderate", supported$evidence_grade)) "partial" else "limited"
    add_cell(
      p, domains[6],
      sprintf("%s\np=%.2g", supported$metabolite, supported$p_indirect),
      cls,
      supported$evidence_grade
    )
  } else {
    add_cell(p, domains[6], "No robust\npath", "missing")
  }

  summary_rows <- rbind(summary_rows, data.frame(
    protein = p,
    tier = tier,
    discovery_nsnp = nsnp,
    discovery_or = or,
    steiger_ok = if (nrow(st)) st$direction_consistent[1] else NA,
    steiger_reversed_snps = if (nrow(st)) st$n_snp_rev[1] else NA,
    egger_testable = !is.na(nsnp) && nsnp >= 3,
    coloc_verdict = if (nrow(m)) clean_text(m$`Coloc verdict`[1]) else NA,
    opengwas_brca_replicated = ogw_sig,
    mediation_sensitivity = if (nrow(ms)) paste(unique(ms$evidence_grade), collapse = "; ") else "No robust path",
    stringsAsFactors = FALSE
  ))
}

row_labels <- paste0(summary_rows$protein, " (", summary_rows$tier, ")")
row_labels[is.na(summary_rows$tier)] <- summary_rows$protein[is.na(summary_rows$tier)]
names(row_labels) <- summary_rows$protein

cells$protein <- factor(cells$protein, levels = rev(priority), labels = rev(row_labels[priority]))
cells$domain <- factor(cells$domain, levels = domains)
cells$class <- factor(cells$class, levels = c("strong", "partial", "limited", "warning", "missing"))

palette <- c(
  strong = "#2E8B57",
  partial = "#D6A13A",
  limited = "#7FA6C7",
  warning = "#B76B5C",
  missing = "#E6E7E9"
)

p <- ggplot(cells, aes(x = domain, y = protein, fill = class)) +
  geom_tile(color = "white", linewidth = 1.35, width = 0.98, height = 0.9) +
  geom_text(aes(label = label), size = 3.0, lineheight = 0.9, color = "#20252B") +
  scale_fill_manual(
    values = palette,
    breaks = c("strong", "partial", "limited", "warning", "missing"),
    labels = c("supports robustness", "suggestive/partial", "limited by design", "potential concern", "not available"),
    name = NULL
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title = "Pleiotropy and robustness summary for priority proteins",
    subtitle = "Most discovery instruments are single-SNP cis-pQTLs; formal MR-Egger pleiotropy testing is therefore limited",
    x = NULL,
    y = NULL,
    caption = "Coloc helps distinguish shared causal signal from linkage. Steiger supports protein-to-cancer directionality for the priority breast cancer hits."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.5, color = "#555A60", margin = margin(b = 10)),
    plot.caption = element_text(size = 8.5, color = "#62676D", hjust = 0, margin = margin(t = 10)),
    axis.text.x = element_text(face = "bold", size = 10.5, color = "#2A2F35", margin = margin(b = 6)),
    axis.text.y = element_text(face = "bold", size = 10.5, color = "#2A2F35", margin = margin(r = 8)),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.55, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(14, 18, 12, 14)
  )

ggsave(file.path(out_fig, "fig14_priority_pleiotropy_robustness_map.png"), p, width = 12.8, height = 5.7, dpi = 320, bg = "white")
ggsave(file.path(out_fig, "fig14_priority_pleiotropy_robustness_map.pdf"), p, width = 12.8, height = 5.7, bg = "white")
write.csv(cells, file.path(out_tab, "priority_pleiotropy_robustness_map_cells.csv"), row.names = FALSE)
write.csv(summary_rows, file.path(out_tab, "priority_pleiotropy_summary.csv"), row.names = FALSE)

cat("Saved pleiotropy robustness outputs:\n")
cat(" - results/figures/fig14_priority_pleiotropy_robustness_map.png\n")
cat(" - results/figures/fig14_priority_pleiotropy_robustness_map.pdf\n")
cat(" - results/sensitivity/priority_pleiotropy_summary.csv\n")
cat(" - results/sensitivity/priority_pleiotropy_robustness_map_cells.csv\n")
