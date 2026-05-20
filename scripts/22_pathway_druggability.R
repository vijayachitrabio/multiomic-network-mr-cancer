#!/usr/bin/env Rscript
# Script 22: Pathway enrichment + druggability for Phase 2 hit proteins
#
# Part A — Pathway enrichment via gprofiler2
#   Input : 17 Phase 2 FDR<0.05 proteins (16 breast + 1 EC)
#   Background: all 701 FinnGen Olink proteins tested
#   Sources: GO (BP, MF, CC), KEGG, Reactome, WikiPathways, CORUM
#   Output: results/pathway/gprofiler_results.csv
#           results/pathway/gprofiler_manhattan.png
#
# Part B — OpenTargets druggability query (REST API)
#   Input : 17 proteins → gene symbols → ENSEMBL IDs via gprofiler
#   Queries: known drugs, tractability, clinical trials
#   Output: results/pathway/opentargets_druggability.csv
#
# Part C — GWAS Catalog overlap (via EBI REST API)
#   Input : rs1260326, rs1047891 — check cancer associations
#   Output: printed to console + logged in results/pathway/gwas_catalog_check.txt

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(gprofiler2)
  library(ggplot2)
  library(httr)
  library(jsonlite)
})

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && a != "") a else b

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
out_dir     <- file.path(project_dir, "results", "pathway")
dir.create(out_dir, showWarnings = FALSE)

# ── 17 Phase 2 hit proteins
hit_proteins <- c("SNX15","EFNA1","FGF5","UMOD","SWAP70","ATRAID","TNFRSF6B",
                   "ITIH3","KLB","PM20D1","TSPAN8","FGFR4","IL34","APOE",
                   "CGREF1","INHBB","ABO")
breast_hits  <- hit_proteins[hit_proteins != "ABO"]

# ── All 701 proteins tested (background)
pqtl <- fread(file.path(project_dir, "data", "pqtl", "pqtl_instruments.csv"))
all_proteins <- unique(pqtl$exposure)
cat(sprintf("Background: %d unique proteins\n", length(all_proteins)))

# ============================================================
# PART A — Pathway Enrichment (gprofiler2)
# ============================================================
cat("\n=== Part A: gprofiler2 pathway enrichment ===\n")

gp_res <- tryCatch({
  gost(
    query            = hit_proteins,
    organism         = "hsapiens",
    ordered_query    = FALSE,
    multi_query      = FALSE,
    significant      = TRUE,
    exclude_iea      = FALSE,          # include electronic GO annotations
    measure_underrepresentation = FALSE,
    evcodes          = TRUE,
    user_threshold   = 0.05,
    correction_method = "fdr",
    custom_bg        = all_proteins,   # use all tested proteins as background
    sources          = c("GO:BP","GO:MF","GO:CC","KEGG","REAC","WP","CORUM")
  )
}, error = function(e) {
  cat("gprofiler2 error:", conditionMessage(e), "\n")
  cat("Retrying without custom background...\n")
  gost(
    query             = hit_proteins,
    organism          = "hsapiens",
    significant       = TRUE,
    evcodes           = TRUE,
    user_threshold    = 0.05,
    correction_method = "fdr",
    sources           = c("GO:BP","GO:MF","GO:CC","KEGG","REAC","WP","CORUM")
  )
})

if (!is.null(gp_res) && !is.null(gp_res$result) && nrow(gp_res$result) > 0) {
  gp_dt <- as.data.table(gp_res$result)
  setorder(gp_dt, p_value)
  cat(sprintf("  %d significant terms (FDR < 0.05)\n", nrow(gp_dt)))
  cat("\n  Top 20 terms:\n")
  print(gp_dt[1:min(20, .N),
              .(source, term_name, term_size, query_size=intersection_size,
                p_value=signif(p_value,3), FDR=signif(p_value,3))])

  # Save full results
  gp_out <- gp_dt[, .(source, term_id, term_name, term_size,
                        intersection_size, p_value, significant,
                        intersection = intersection_size,
                        genes = intersection_input)]
  fwrite(gp_out, file.path(out_dir, "gprofiler_results.csv"))

  # ── Dot plot — top 20 terms coloured by source
  top20 <- gp_dt[1:min(20,.N)]
  top20[, term_short := ifelse(nchar(term_name) > 45,
                               paste0(substr(term_name,1,42),"..."), term_name)]
  source_cols <- c(
    "GO:BP" = "#4A90D9", "GO:MF" = "#7B68EE", "GO:CC" = "#20B2AA",
    "KEGG"  = "#E07B54", "REAC"  = "#2C9E5B", "WP"    = "#9B59B6",
    "CORUM" = "#C0392B"
  )

  p_enrich <- ggplot(top20, aes(x = -log10(p_value),
                                 y = reorder(term_short, -p_value),
                                 size = intersection_size,
                                 colour = source)) +
    geom_point(alpha = 0.85) +
    scale_colour_manual(name = "Database", values = source_cols,
                        drop = TRUE) +
    scale_size_continuous(name = "Proteins\nin term", range = c(3,10)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    labs(
      title    = "Pathway enrichment: 17 Phase 2 protein-cancer hits",
      subtitle = sprintf("gprofiler2 FDR < 0.05; background = %d FinnGen Olink proteins",
                         length(all_proteins)),
      x        = expression(-log[10](FDR)),
      y        = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "right")

  ggsave(file.path(out_dir, "gprofiler_dotplot.pdf"), p_enrich, width = 10, height = 7)
  ggsave(file.path(out_dir, "gprofiler_dotplot.png"), p_enrich, width = 10, height = 7, dpi = 300)
  cat("  Saved: gprofiler_dotplot.pdf|png\n")
} else {
  cat("  No significant terms found.\n")
}

# ============================================================
# PART B — OpenTargets druggability (REST API)
# ============================================================
cat("\n=== Part B: OpenTargets druggability query ===\n")

# Convert gene symbols to ENSEMBL IDs via gprofiler name conversion
name_conv <- tryCatch(
  gconvert(hit_proteins, organism = "hsapiens", target = "ENSG",
           mthreshold = 1, filter_na = TRUE),
  error = function(e) NULL
)

drug_results <- list()

if (!is.null(name_conv) && nrow(name_conv) > 0) {
  ensg_ids <- name_conv$target
  gene_map  <- setNames(name_conv$target, name_conv$input)
  cat(sprintf("  Converted %d/%d gene symbols to ENSEMBL IDs\n",
              length(ensg_ids), length(hit_proteins)))

  for (i in seq_along(ensg_ids)) {
    ensg <- ensg_ids[i]
    gene <- names(gene_map)[gene_map == ensg][1]

    url <- paste0("https://api.platform.opentargets.org/api/v4/graphql")
    query <- sprintf('{
      target(ensemblId: "%s") {
        id
        approvedSymbol
        approvedName
        tractability {
          label
          modality
          value
        }
        knownDrugs {
          count
          rows {
            drug { name }
            phase
            status
            disease { name }
          }
        }
      }
    }', ensg)

    resp <- tryCatch({
      req <- list(query = query)
      body_json <- jsonlite::toJSON(req, auto_unbox = TRUE)
      resp_raw <- httr::POST(url,
        httr::add_headers("Content-Type" = "application/json"),
        body = body_json,
        encode = "raw",
        httr::timeout(15))
      if (httr::status_code(resp_raw) == 200) {
        jsonlite::fromJSON(httr::content(resp_raw, "text", encoding="UTF-8"),
                           simplifyVector = FALSE)
      } else NULL
    }, error = function(e) NULL)

    if (!is.null(resp) && !is.null(resp$data$target)) {
      tgt   <- resp$data$target
      drugs <- tgt$knownDrugs
      n_drugs <- if (!is.null(drugs$count)) drugs$count else 0L

      # Tractability — pull SM and AB labels
      tract_sm <- tract_ab <- "unknown"
      if (!is.null(tgt$tractability) && length(tgt$tractability) > 0) {
        for (tr in tgt$tractability) {
          if (!is.null(tr$modality)) {
            if (grepl("SM", tr$modality, ignore.case=TRUE) && isTRUE(tr$value))
              tract_sm <- tr$label
            if (grepl("AB|Antibody", tr$modality, ignore.case=TRUE) && isTRUE(tr$value))
              tract_ab <- tr$label
          }
        }
      }

      # Top drug (highest clinical phase)
      top_drug <- top_phase <- top_indication <- NA_character_
      if (!is.null(drugs$rows) && length(drugs$rows) > 0) {
        drug_dt <- rbindlist(lapply(drugs$rows, function(r) {
          data.table(
            drug       = if (!is.null(r$drug$name)) r$drug$name else NA_character_,
            phase      = if (!is.null(r$phase))      as.integer(r$phase) else NA_integer_,
            status     = if (!is.null(r$status))     r$status else NA_character_,
            indication = if (!is.null(r$disease$name)) r$disease$name else NA_character_
          )
        }), fill = TRUE)
        drug_dt <- drug_dt[!is.na(phase)][order(-phase)]
        if (nrow(drug_dt) > 0) {
          top_drug       <- drug_dt$drug[1]
          top_phase      <- as.character(drug_dt$phase[1])
          top_indication <- drug_dt$indication[1]
        }
      }

      drug_results[[gene]] <- data.table(
        protein         = gene,
        ensembl_id      = ensg,
        approved_name   = tgt$approvedName %||% NA_character_,
        n_known_drugs   = n_drugs,
        tractability_SM = tract_sm,
        tractability_AB = tract_ab,
        top_drug        = top_drug,
        top_drug_phase  = top_phase,
        top_indication  = top_indication
      )
      cat(sprintf("  %s (%s): %d known drugs; top=%s (Phase %s)\n",
                  gene, ensg, n_drugs,
                  ifelse(is.na(top_drug),"none",top_drug),
                  ifelse(is.na(top_phase),"—",top_phase)))
    } else {
      cat(sprintf("  %s: no OpenTargets data\n", gene))
    }
    Sys.sleep(0.3)  # polite API rate
  }
} else {
  cat("  Gene symbol conversion failed — skipping OpenTargets query\n")
}

if (length(drug_results) > 0) {
  drug_dt <- rbindlist(drug_results, fill = TRUE)
  setorder(drug_dt, -n_known_drugs)
  cat("\n=== Druggability Summary ===\n")
  print(drug_dt[, .(protein, n_known_drugs, top_drug, top_drug_phase,
                     tractability_SM, tractability_AB)])
  fwrite(drug_dt, file.path(out_dir, "opentargets_druggability.csv"))
  cat("  Saved: opentargets_druggability.csv\n")
}

# ============================================================
# PART C — GWAS Catalog check for key coloc SNPs
# ============================================================
cat("\n=== Part C: GWAS Catalog association check ===\n")

check_snps <- c("rs1260326", "rs1047891")

gwas_log <- character(0)
for (snp in check_snps) {
  url <- paste0("https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/",
                snp, "/associations?projection=associationBySnp&size=200")
  resp <- tryCatch({
    raw <- httr::GET(url, httr::timeout(20))
    if (httr::status_code(raw) == 200)
      jsonlite::fromJSON(httr::content(raw, "text", encoding="UTF-8"),
                         simplifyVector = FALSE)
    else NULL
  }, error = function(e) NULL)

  if (!is.null(resp) && !is.null(resp[["_embedded"]])) {
    assoc <- resp[["_embedded"]][["associations"]]
    cat(sprintf("\n%s — %d GWAS Catalog associations:\n", snp, length(assoc)))

    rows <- lapply(assoc, function(a) {
      # Pull trait labels
      traits <- if (!is.null(a$efoTraits)) {
        paste(sapply(a$efoTraits, function(t) t$trait %||% "?"), collapse="; ")
      } else "?"
      pval_mant <- a$pvalueMantissa %||% NA
      pval_exp  <- a$pvalueExponent %||% NA
      pval_str  <- if (!is.na(pval_mant) && !is.na(pval_exp))
        sprintf("%se%s", pval_mant, pval_exp) else "?"
      data.table(snp=snp, trait=traits, pvalue=pval_str,
                 beta=a$betaNum %||% NA_real_,
                 or=a$orPerCopyNum %||% NA_real_)
    })
    rows_dt <- rbindlist(rows, fill=TRUE)

    # Flag cancer hits
    cancer_rows <- rows_dt[grepl("cancer|carcinoma|tumor|tumour|neoplasm",
                                  trait, ignore.case=TRUE)]
    if (nrow(cancer_rows) > 0) {
      cat("  *** CANCER ASSOCIATIONS FOUND ***\n")
      print(cancer_rows)
    } else {
      cat("  No cancer associations in GWAS Catalog for this SNP.\n")
    }
    cat("  All traits:", paste(unique(rows_dt$trait), collapse=" | "), "\n")
    gwas_log <- c(gwas_log, sprintf("\n%s:\n%s", snp,
                                     paste(unique(rows_dt$trait), collapse="\n")))
  } else {
    cat(sprintf("\n%s: No GWAS Catalog associations found (or API error)\n", snp))
    gwas_log <- c(gwas_log, sprintf("\n%s: no associations found", snp))
  }
}

writeLines(gwas_log, file.path(out_dir, "gwas_catalog_check.txt"))

cat("\nOutputs in:", out_dir, "\n")
cat("Done.\n")
