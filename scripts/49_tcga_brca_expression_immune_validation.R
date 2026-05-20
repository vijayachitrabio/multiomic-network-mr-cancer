#!/usr/bin/env Rscript

# TCGA-BRCA expression and immune-signature validation using lightweight UCSC Xena matrices.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

project_dir <- "."
cache_dir <- file.path(project_dir, "data", "tcga_brca")
out_dir <- file.path(project_dir, "results", "tcga_immune")
fig_dir <- file.path(project_dir, "results", "figures")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(cache_dir, "HiSeqV2.gz")
clinical_file <- file.path(cache_dir, "BRCA_clinicalMatrix")
if (!file.exists(expr_file)) stop("Missing Xena expression file: ", expr_file)
if (!file.exists(clinical_file)) stop("Missing Xena clinical file: ", clinical_file)

targets <- c("EFNA1", "TNFRSF6B", "ATRAID", "ITIH3", "IL34", "FGF5", "APOE")

message("Reading TCGA-BRCA Xena expression matrix...")
expr_dt <- fread(expr_file)
setnames(expr_dt, 1, "gene")
genes <- expr_dt$gene
expr_mat <- as.matrix(expr_dt[, -"gene"])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- genes

if (anyDuplicated(rownames(expr_mat))) {
  m <- rowMeans(expr_mat, na.rm = TRUE)
  keep <- !duplicated(rownames(expr_mat)[order(rownames(expr_mat), -m)])
  ord <- order(rownames(expr_mat), -m)
  expr_mat <- expr_mat[ord, , drop = FALSE]
  expr_mat <- expr_mat[keep, , drop = FALSE]
}

clinical <- fread(clinical_file)
setnames(clinical, 1, "sample")
fwrite(clinical, file.path(out_dir, "tcga_brca_clinical_xena.csv"))

sample_meta <- data.table(
  sample = colnames(expr_mat),
  patient = substr(colnames(expr_mat), 1, 12),
  sample_type_code = substr(colnames(expr_mat), 14, 15)
)
sample_meta[, tcga_sample_type := fifelse(sample_type_code == "01", "Primary Tumor",
                                   fifelse(sample_type_code == "11", "Solid Tissue Normal", "Other"))]
sample_meta <- merge(sample_meta, clinical, by = "sample", all.x = TRUE)
sample_type_by_expr <- sample_meta[match(colnames(expr_mat), sample), tcga_sample_type]

target_genes_present <- intersect(targets, rownames(expr_mat))
target_expr <- as.data.table(t(expr_mat[target_genes_present, , drop = FALSE]), keep.rownames = "sample")
target_expr <- merge(sample_meta, target_expr, by = "sample", all.y = TRUE)
fwrite(target_expr, file.path(out_dir, "tcga_brca_target_gene_xena_expression.csv"))

tn_results <- rbindlist(lapply(target_genes_present, function(g) {
  dat <- data.table(expr = as.numeric(expr_mat[g, ]), sample_type = sample_type_by_expr)
  dat <- dat[sample_type %in% c("Primary Tumor", "Solid Tissue Normal") & !is.na(expr)]
  if (uniqueN(dat$sample_type) < 2) return(NULL)
  fit <- lm(expr ~ relevel(factor(sample_type), ref = "Solid Tissue Normal"), data = dat)
  data.table(
    gene = g,
    n_tumour = sum(dat$sample_type == "Primary Tumor"),
    n_normal = sum(dat$sample_type == "Solid Tissue Normal"),
    mean_tumour = mean(dat$expr[dat$sample_type == "Primary Tumor"]),
    mean_normal = mean(dat$expr[dat$sample_type == "Solid Tissue Normal"]),
    logFC_tumour_vs_normal = unname(coef(fit)[2]),
    p_tumour_vs_normal = summary(fit)$coefficients[2, 4]
  )
}), fill = TRUE)
tn_results[, fdr_tumour_vs_normal := p.adjust(p_tumour_vs_normal, method = "BH")]
fwrite(tn_results, file.path(out_dir, "tcga_brca_target_tumour_vs_normal.csv"))

er_col <- if ("ER_Status_nature2012" %in% names(target_expr)) "ER_Status_nature2012" else
  if ("breast_carcinoma_estrogen_receptor_status" %in% names(target_expr)) "breast_carcinoma_estrogen_receptor_status" else NA_character_
if (!is.na(er_col)) {
  er_results <- rbindlist(lapply(target_genes_present, function(g) {
    dat <- copy(target_expr)[tcga_sample_type == "Primary Tumor", .(expr = get(g), er_raw = get(er_col))]
    dat[, er_status := fifelse(grepl("pos", er_raw, ignore.case = TRUE), "ER-positive",
                        fifelse(grepl("neg", er_raw, ignore.case = TRUE), "ER-negative", NA_character_))]
    dat <- dat[!is.na(er_status) & !is.na(expr)]
    if (uniqueN(dat$er_status) < 2) return(NULL)
    fit <- lm(expr ~ relevel(factor(er_status), ref = "ER-negative"), data = dat)
    data.table(
      gene = g,
      er_column = er_col,
      n_er_positive = sum(dat$er_status == "ER-positive"),
      n_er_negative = sum(dat$er_status == "ER-negative"),
      mean_er_positive = mean(dat$expr[dat$er_status == "ER-positive"]),
      mean_er_negative = mean(dat$expr[dat$er_status == "ER-negative"]),
      beta_er_positive_vs_negative = unname(coef(fit)[2]),
      p_er = summary(fit)$coefficients[2, 4]
    )
  }), fill = TRUE)
  er_results[, fdr_er := p.adjust(p_er, method = "BH")]
  fwrite(er_results, file.path(out_dir, "tcga_brca_target_er_status.csv"))
}

pam_col <- if ("PAM50Call_RNAseq" %in% names(target_expr)) "PAM50Call_RNAseq" else
  if ("PAM50_mRNA_nature2012" %in% names(target_expr)) "PAM50_mRNA_nature2012" else NA_character_
if (!is.na(pam_col)) {
  pam_results <- rbindlist(lapply(target_genes_present, function(g) {
    dat <- copy(target_expr)[tcga_sample_type == "Primary Tumor", .(expr = get(g), subtype = get(pam_col))]
    dat <- dat[!is.na(expr) & !is.na(subtype) & subtype != ""]
    if (uniqueN(dat$subtype) < 2) return(NULL)
    fit <- lm(expr ~ subtype, data = dat)
    a <- anova(fit)
    means <- dat[, .(mean_expr = mean(expr), n = .N), by = subtype]
    data.table(gene = g, pam50_column = pam_col, p_pam50 = a$`Pr(>F)`[1])[
      , paste(means$subtype, round(means$mean_expr, 3), sep = "=", collapse = "; "),
      by = .(gene, pam50_column, p_pam50)
    ]
  }), fill = TRUE)
  setnames(pam_results, "V1", "subtype_means")
  pam_results[, fdr_pam50 := p.adjust(p_pam50, method = "BH")]
  fwrite(pam_results, file.path(out_dir, "tcga_brca_target_pam50.csv"))
}

immune_signatures <- list(
  CD8_T = c("CD8A", "CD8B", "GZMB", "PRF1", "NKG7"),
  Cytotoxic = c("GZMA", "GZMB", "PRF1", "GNLY", "NKG7"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "IKZF2"),
  Macrophage = c("CD68", "CD163", "MRC1", "CSF1R", "MSR1"),
  M1_Macrophage = c("IL1B", "TNF", "CXCL9", "CXCL10", "CD86"),
  M2_Macrophage = c("CD163", "MRC1", "MSR1", "IL10", "TGFB1"),
  Dendritic = c("ITGAX", "CD1C", "CLEC9A", "LAMP3", "CCR7"),
  B_cell = c("MS4A1", "CD79A", "CD79B", "BANK1"),
  Immune_Checkpoint = c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT")
)

zexpr <- t(scale(t(expr_mat)))
score_dt <- rbindlist(lapply(names(immune_signatures), function(sig) {
  gs <- intersect(immune_signatures[[sig]], rownames(zexpr))
  if (length(gs) == 0) return(NULL)
  data.table(sample = colnames(zexpr), signature = sig, score = colMeans(zexpr[gs, , drop = FALSE], na.rm = TRUE))
}), fill = TRUE)
immune_wide <- dcast(score_dt, sample ~ signature, value.var = "score")
immune_wide <- merge(sample_meta, immune_wide, by = "sample", all.y = TRUE)
fwrite(immune_wide, file.path(out_dir, "tcga_brca_immune_signature_scores.csv"))

tumour_samples <- sample_meta[tcga_sample_type == "Primary Tumor", sample]
immune_cols <- setdiff(names(immune_wide), names(sample_meta))
immune_cor <- rbindlist(lapply(target_genes_present, function(g) {
  rbindlist(lapply(immune_cols, function(sig) {
    dt <- merge(
      data.table(sample = tumour_samples, expression = as.numeric(expr_mat[g, tumour_samples])),
      immune_wide[, .(sample, score = get(sig))],
      by = "sample"
    )
    ok <- complete.cases(dt[, .(expression, score)])
    if (sum(ok) < 20) return(NULL)
    ct <- suppressWarnings(cor.test(dt$expression[ok], dt$score[ok], method = "spearman"))
    data.table(gene = g, signature = sig, n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
  }))
}), fill = TRUE)
immune_cor[, fdr := p.adjust(p, method = "BH")]
fwrite(immune_cor, file.path(out_dir, "tcga_brca_target_immune_correlations.csv"))

top_cor <- immune_cor[order(fdr)][1:min(.N, 45)]
p <- ggplot(top_cor, aes(x = signature, y = gene, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#3568a6", mid = "white", high = "#b4443c", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Spearman rho") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "fig9_tcga_immune_correlations.png"), p, width = 7, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "fig9_tcga_immune_correlations.pdf"), p, width = 7, height = 4)

message("Saved TCGA/immune validation outputs to: ", out_dir)
