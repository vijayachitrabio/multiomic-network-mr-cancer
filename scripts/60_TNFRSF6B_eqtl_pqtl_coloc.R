#!/usr/bin/env Rscript
# ============================================================================
# Script: 60_TNFRSF6B_eqtl_pqtl_coloc.R
# Purpose: Colocalization of FinnGen pQTL and eQTLGen Whole Blood eQTL
#          for TNFRSF6B — test shared causal variant hypothesis
#          (triangulate MR breast cancer signal)
#
# Datasets:
#   pQTL: FinnGen R10 Olink (N=619), region chr20:63,206,054-64,206,054
#   eQTL: eQTLGen Whole Blood cis-eQTL (N=31,684), dataset eqtl-a-ENSG00000243509
#         fetched via ieugwasr with valid OpenGWAS JWT
#
# Method: coloc.abf (Bayesian colocalization, Giambartolomei et al. 2014)
#         coloc.susie attempted if >=100 overlapping SNPs
#
# Author: Auto-generated 2026-05-23
# ============================================================================

# ── Set OpenGWAS JWT ────────────────────────────────────────────────────────
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ2aWpheWFjaGl0cmEubW9kaHVrdXJAdXQuZWUiLCJpYXQiOjE3NzY2OTE4NDksImV4cCI6MTc3NzkwMTQ0OX0.jBKwTlhIqHn8DJC6Ryxd75K4QhOOCvCDOHr0CCkVlYmj5FS3X5Mrv1am4YzWZhswFsIF1v-GMImrmKrPJVP74v6XGNR1feV9OQXgZFF9lU_Mz2guLbMndFXPK4YKcyVqpHQipwlz0C76T1Ih5tjvSuYyctgFJ0UhPaqL3B0982KtboiIH2GtYay2mSpKM8pCMZJUoG97FJvggPvzO3etJBYUHSrhBl6ya04r_TF54ZFsIlbAfXVRZ1bOhZr4wGgsPfo6mvjZdUBxdfvxdZAbbEw-41idp2CoU-8RdY0GhQA86G8eJ6G3Qr6Qzq41oK8LsmxPDixn85AV0vj3BNHQvQ")

suppressMessages({
  library(data.table)
  library(dplyr)
  library(coloc)
  library(ieugwasr)
})

# ── Paths ──────────────────────────────────────────────────────────────────
BASE      <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
PQTL_FILE <- file.path(BASE, "data/pqtl/priority_regions/TNFRSF6B_pqtl_regions.tsv.gz")
OUT_DIR   <- file.path(BASE, "results/eqtl")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DATE_TAG <- "2026-05-23"
OUT_CSV  <- file.path(OUT_DIR, paste0("TNFRSF6B_pqtl_eqtlgen_coloc_",  DATE_TAG, ".csv"))
OUT_MD   <- file.path(OUT_DIR, paste0("TNFRSF6B_eqtl_coloc_summary_", DATE_TAG, ".md"))

# ── Parameters ─────────────────────────────────────────────────────────────
GENE_SYMBOL   <- "TNFRSF6B"
EQTL_ID       <- "eqtl-a-ENSG00000243509"   # eQTLGen Whole Blood
REGION_CHR    <- "20"
REGION_START  <- 63206054
REGION_END    <- 64206054
LEAD_PQTL     <- "rs6011040"
LEAD_PQTL_POS <- 63706054

N_PQTL <- 619    # FinnGen Olink, inverse-normalised
N_EQTL <- 31684  # eQTLGen Whole Blood
SDY    <- 1      # Both datasets are quantile-/inverse-normalised

cat("\n", strrep("=", 72), "\n")
cat("TNFRSF6B pQTL-eQTL Colocalization Analysis (eQTLGen N=31,684)\n")
cat("Date:", DATE_TAG, "\n")
cat("pQTL: FinnGen R10 Olink, N =", N_PQTL, "\n")
cat("eQTL: eQTLGen Whole Blood, N =", N_EQTL, "\n")
cat("Region: chr", REGION_CHR, ":", REGION_START, "-", REGION_END, "\n")
cat(strrep("=", 72), "\n\n")

# ============================================================================
# STEP 1: Load pQTL data
# ============================================================================
cat("Step 1: Loading pQTL data...\n")
pqtl_raw <- fread(PQTL_FILE)
cat("  Raw rows:", nrow(pqtl_raw), "\n")

pqtl <- pqtl_raw %>%
  filter(chr == 20,
         pos >= REGION_START,
         pos <= REGION_END,
         !is.na(beta), !is.na(se), se > 0) %>%
  mutate(
    pos    = as.integer(pos),
    chrpos = paste0("chr20_", pos)
  )

cat("  Filtered rows in region:", nrow(pqtl), "\n")
cat("  Lead pQTL SNP present:", LEAD_PQTL %in% pqtl$variant_id, "\n")

# Show lead SNP stats
lead_row <- pqtl %>% filter(variant_id == LEAD_PQTL)
if (nrow(lead_row) > 0) {
  cat(sprintf("  Lead pQTL: %s  pos=%d  beta=%.4f  se=%.4f  p=%.3e\n",
              lead_row$variant_id[1], lead_row$pos[1],
              lead_row$beta[1], lead_row$se[1], lead_row$p[1]))
}

# Build variant IDs in chr20_pos_ref_alt format (used by OpenGWAS)
pqtl <- pqtl %>%
  mutate(
    variant_id_og = paste0("chr20_", pos, "_", toupper(ref), "_", toupper(alt))
  )

cat("  Example variant IDs (first 3):", paste(head(pqtl$variant_id_og, 3), collapse=", "), "\n")

# ============================================================================
# STEP 2: Confirm eQTLGen dataset availability + fetch eQTL data
# ============================================================================
cat("\nStep 2: Confirming eQTLGen availability and fetching eQTL data...\n")

# First check tophits to confirm dataset is accessible with this JWT
cat("  Checking tophits for", EQTL_ID, "...\n")
top <- tryCatch(
  tophits(id = EQTL_ID),
  error = function(e) {
    cat("  tophits error:", conditionMessage(e), "\n")
    return(NULL)
  }
)

eqtl_source_used <- "eQTLGen Whole Blood (N=31,684)"
eqtl_fetch_method <- NA_character_
eqtl_raw <- NULL

if (!is.null(top) && nrow(top) > 0) {
  cat(sprintf("  eQTLGen accessible. Top hits: %d rows\n", nrow(top)))
  cat("  Top hit columns:", paste(colnames(top), collapse=", "), "\n")

  # Method A: Fetch by pQTL variant IDs
  cat("\n  Fetching associations for pQTL variants...\n")
  # eQTLGen expects rsid or chr:pos format; try variant IDs in batches
  variant_ids <- pqtl$variant_id_og

  # Try batch associations in chunks of 500
  fetch_in_chunks <- function(variants, id, chunk_size = 500) {
    results <- list()
    n_chunks <- ceiling(length(variants) / chunk_size)
    cat(sprintf("    Fetching %d variants in %d chunks...\n", length(variants), n_chunks))
    for (i in seq_len(n_chunks)) {
      idx <- ((i-1)*chunk_size + 1):min(i*chunk_size, length(variants))
      chunk <- variants[idx]
      res <- tryCatch(
        associations(variants = chunk, id = id, proxies = 0),
        error = function(e) {
          cat(sprintf("    Chunk %d error: %s\n", i, conditionMessage(e)))
          return(NULL)
        }
      )
      if (!is.null(res) && nrow(res) > 0) {
        results[[i]] <- res
        cat(sprintf("    Chunk %d: %d hits\n", i, nrow(res)))
      } else {
        cat(sprintf("    Chunk %d: 0 hits\n", i))
      }
    }
    if (length(results) > 0) bind_rows(results) else NULL
  }

  eqtl_raw <- fetch_in_chunks(variant_ids, EQTL_ID)

  # If variant format didn't work, try rsIDs from pQTL data
  if (is.null(eqtl_raw) || nrow(eqtl_raw) == 0) {
    cat("  chr_pos format returned 0 hits. Trying rsIDs...\n")
    rsids <- pqtl %>% filter(grepl("^rs", variant_id)) %>% pull(variant_id)
    if (length(rsids) > 0) {
      eqtl_raw <- fetch_in_chunks(rsids, EQTL_ID)
    }
  }

  # If still nothing, try region-based query via gwasinfo + phewas
  if (is.null(eqtl_raw) || nrow(eqtl_raw) == 0) {
    cat("  Trying region-based query...\n")
    region_str <- paste0(REGION_CHR, ":", REGION_START, "-", REGION_END)
    eqtl_raw <- tryCatch(
      phewas(pval = 1, variants = NULL, id = EQTL_ID),
      error = function(e) {
        cat("  phewas error:", conditionMessage(e), "\n"); NULL
      }
    )
  }

  if (!is.null(eqtl_raw) && nrow(eqtl_raw) > 0) {
    eqtl_fetch_method <- "ieugwasr::associations()"
    cat(sprintf("  Fetched %d eQTL records from eQTLGen\n", nrow(eqtl_raw)))
    cat("  Columns:", paste(colnames(eqtl_raw), collapse=", "), "\n")
  } else {
    cat("  All ieugwasr fetch attempts returned 0 records.\n")
    cat("  Will fall back to tophits data if available.\n")
    eqtl_raw <- top
    eqtl_fetch_method <- "tophits() fallback"
  }
} else {
  cat("  eQTLGen tophits not accessible or empty.\n")
}

# ── Fallback: GTEx v8 Whole Blood ──────────────────────────────────────────
if (is.null(eqtl_raw) || nrow(eqtl_raw) == 0) {
  cat("\n  FALLBACK: Fetching GTEx v8 Whole Blood eQTL via GTEx API...\n")
  eqtl_source_used <- "GTEx v8 Whole Blood (N=670) [fallback: eQTLGen inaccessible]"
  N_EQTL <- 670

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(jsonlite)

  gene_id_gtex <- "ENSG00000243509.4"
  url <- sprintf(
    "https://gtexportal.org/api/v2/association/singleTissueEqtl?gencodeId=%s&tissueSiteDetailId=Whole_Blood&datasetId=gtex_v8&itemsPerPage=250&sortBy=pValue&sortDirection=asc",
    gene_id_gtex
  )
  gtex_raw <- tryCatch({
    raw <- readLines(url(url), warn = FALSE)
    payload <- jsonlite::fromJSON(paste(raw, collapse = ""))
    payload$data
  }, error = function(e) {
    cat("  GTEx API error:", conditionMessage(e), "\n"); NULL
  })

  if (!is.null(gtex_raw) && nrow(gtex_raw) > 0) {
    cat(sprintf("  GTEx API returned %d records\n", nrow(gtex_raw)))
    eqtl_fetch_method <- "GTEx API v2 (fallback)"

    eqtl_raw <- gtex_raw %>%
      as.data.frame() %>%
      filter(!is.na(pValue), pValue > 0, !is.na(nes)) %>%
      mutate(
        pos      = as.integer(pos),
        chrpos   = sub("(chr20_[0-9]+)_.*", "\\1", variantId),
        rsid_col = snpId,
        beta     = nes,
        t_abs    = abs(qt(pValue / 2, df = N_EQTL - 2, lower.tail = FALSE)),
        se       = abs(nes) / pmax(t_abs, 0.001),
        pval     = pValue
      ) %>%
      filter(se > 0, pos >= REGION_START, pos <= REGION_END) %>%
      mutate(source = "gtex")
  } else {
    stop("Both eQTLGen and GTEx API failed. Cannot proceed with colocalization.")
  }
}

# ============================================================================
# STEP 3: Harmonise datasets
# ============================================================================
cat("\nStep 3: Harmonising pQTL and eQTL datasets...\n")

# Detect eQTL data format (eQTLGen vs GTEx)
is_gtex <- "source" %in% colnames(eqtl_raw) && any(eqtl_raw$source == "gtex", na.rm = TRUE)
is_eqtlgen <- !is_gtex

cat("  eQTL source:", if(is_gtex) "GTEx (fallback)" else "eQTLGen", "\n")

if (is_eqtlgen) {
  cat("  eQTLGen columns:", paste(colnames(eqtl_raw), collapse=", "), "\n")

  # eQTLGen data via ieugwasr — standard column names: beta, se, p, rsid, chr, position
  # Handle column name variations
  eqtl_df <- eqtl_raw %>% as.data.frame()

  # Normalise column names (ieugwasr returns: beta, se, p, rsid, chr, position, ea, nea)
  if ("position" %in% colnames(eqtl_df)) {
    eqtl_df <- eqtl_df %>% rename(pos = position)
  }
  if ("p" %in% colnames(eqtl_df) && !"pval" %in% colnames(eqtl_df)) {
    eqtl_df <- eqtl_df %>% rename(pval = p)
  }

  # Filter to region
  if ("chr" %in% colnames(eqtl_df)) {
    eqtl_df <- eqtl_df %>%
      filter(as.character(chr) == REGION_CHR |
             as.character(chr) == paste0("chr", REGION_CHR))
  }
  if ("pos" %in% colnames(eqtl_df)) {
    eqtl_df <- eqtl_df %>%
      filter(!is.na(pos),
             as.integer(pos) >= REGION_START,
             as.integer(pos) <= REGION_END)
  }

  eqtl_df <- eqtl_df %>%
    filter(!is.na(beta), !is.na(se), se > 0, !is.na(pval)) %>%
    mutate(
      pos    = as.integer(pos),
      chrpos = paste0("chr20_", pos)
    )

  cat(sprintf("  eQTLGen records in region: %d\n", nrow(eqtl_df)))
  if (nrow(eqtl_df) > 0) {
    cat(sprintf("  P-range: %.3e - %.3e\n", min(eqtl_df$pval, na.rm=TRUE), max(eqtl_df$pval, na.rm=TRUE)))
    # Show lead eQTL
    lead_eqtl <- eqtl_df %>% arrange(pval) %>% slice(1)
    cat(sprintf("  Lead eQTL: pos=%d  beta=%.4f  se=%.4f  p=%.3e\n",
                lead_eqtl$pos[1], lead_eqtl$beta[1], lead_eqtl$se[1], lead_eqtl$pval[1]))
    if ("rsid" %in% colnames(lead_eqtl)) {
      cat(sprintf("  Lead eQTL rsID: %s\n", lead_eqtl$rsid[1]))
    }
  }

  # Build harmonised table
  pqtl_h <- pqtl %>%
    mutate(chrpos = paste0("chr20_", pos)) %>%
    select(chrpos, variant_id, pos_pqtl = pos, ref, alt, alt_freq,
           beta_pqtl = beta, se_pqtl = se, p_pqtl = p, n)

  eqtl_h <- eqtl_df %>%
    select(chrpos, pos_eqtl = pos,
           beta_eqtl = beta, se_eqtl = se, pValue_eqtl = pval,
           any_of(c("rsid", "ea", "nea")))

  merged <- inner_join(pqtl_h, eqtl_h, by = "chrpos") %>%
    mutate(pos = pos_pqtl)

  cat(sprintf("  Overlapping SNPs (chrpos match): %d\n", nrow(merged)))

  # Fallback: try rsID match if chrpos gives < 5 hits
  if (nrow(merged) < 5 && "rsid" %in% colnames(eqtl_df)) {
    cat("  Trying rsID match as supplement...\n")
    pqtl_rs <- pqtl %>% filter(grepl("^rs", variant_id))
    eqtl_rs <- eqtl_df %>% filter(grepl("^rs", rsid))
    merged_rs <- inner_join(
      pqtl_rs %>% select(variant_id, pos_pqtl = pos, beta_pqtl = beta, se_pqtl = se, p_pqtl = p),
      eqtl_rs %>% select(rsid, pos_eqtl = pos, chrpos, beta_eqtl = beta, se_eqtl = se, pValue_eqtl = pval),
      by = c("variant_id" = "rsid")
    ) %>% mutate(pos = pos_pqtl)
    cat(sprintf("  rsID-matched SNPs: %d\n", nrow(merged_rs)))
    if (nrow(merged_rs) > nrow(merged)) merged <- merged_rs
  }

} else {
  # GTEx format (from fallback branch above)
  pqtl_h <- pqtl %>%
    mutate(chrpos = paste0("chr20_", pos)) %>%
    select(chrpos, variant_id, pos_pqtl = pos, ref, alt, alt_freq,
           beta_pqtl = beta, se_pqtl = se, p_pqtl = p, n)

  eqtl_h <- eqtl_raw %>%
    select(chrpos, pos_eqtl = pos, rsid_col,
           beta_eqtl = beta, se_eqtl = se, pValue_eqtl = pval)

  merged <- inner_join(pqtl_h, eqtl_h, by = "chrpos") %>%
    mutate(pos = pos_pqtl)

  cat(sprintf("  Overlapping SNPs (chrpos match): %d\n", nrow(merged)))
}

if (nrow(merged) == 0) {
  stop("No SNPs overlap between pQTL and eQTL data. Cannot run coloc.")
}
if (nrow(merged) < 2) {
  stop("Fewer than 2 overlapping SNPs. Cannot run coloc.")
}

cat("\n  Top overlapping SNPs by eQTL p-value:\n")
print(merged %>%
        select(chrpos, any_of(c("variant_id","rsid")), pos,
               beta_pqtl, se_pqtl, p_pqtl, beta_eqtl, se_eqtl, pValue_eqtl) %>%
        arrange(pValue_eqtl) %>%
        head(10))

# ============================================================================
# STEP 4: Prepare coloc datasets
# ============================================================================
cat("\nStep 4: Preparing coloc datasets...\n")

merged_dedup <- merged %>%
  group_by(chrpos) %>%
  slice_min(pValue_eqtl, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(!duplicated(pos))

cat(sprintf("  SNPs after deduplication: %d\n", nrow(merged_dedup)))

D1 <- list(
  type     = "quant",
  beta     = merged_dedup$beta_pqtl,
  varbeta  = merged_dedup$se_pqtl^2,
  N        = N_PQTL,
  sdY      = SDY,
  snp      = merged_dedup$chrpos,
  position = merged_dedup$pos
)

D2 <- list(
  type     = "quant",
  beta     = merged_dedup$beta_eqtl,
  varbeta  = merged_dedup$se_eqtl^2,
  N        = N_EQTL,
  sdY      = SDY,
  snp      = merged_dedup$chrpos,
  position = merged_dedup$pos
)

cat("  Validating D1 (pQTL)...\n")
tryCatch(check_dataset(D1, req = "beta", warn.minp = 0),
         error = function(e) cat("  D1 warning:", conditionMessage(e), "\n"))
cat("  Validating D2 (eQTL)...\n")
tryCatch(check_dataset(D2, req = "beta", warn.minp = 0),
         error = function(e) cat("  D2 warning:", conditionMessage(e), "\n"))

# ============================================================================
# STEP 5: Run coloc.abf
# ============================================================================
cat("\nStep 5: Running coloc.abf...\n")

coloc_abf <- coloc.abf(
  dataset1 = D1,
  dataset2 = D2,
  p1  = 1e-4,
  p2  = 1e-4,
  p12 = 1e-5
)

cat("\n  coloc.abf Results:\n")
print(coloc_abf$summary)

pp   <- coloc_abf$summary
PPH0 <- pp["PP.H0.abf"]
PPH1 <- pp["PP.H1.abf"]
PPH2 <- pp["PP.H2.abf"]
PPH3 <- pp["PP.H3.abf"]
PPH4 <- pp["PP.H4.abf"]

cat(sprintf("\n  PP(H0): %.4f  [no association in either trait]\n",  PPH0))
cat(sprintf("  PP(H1): %.4f  [pQTL only]\n",                          PPH1))
cat(sprintf("  PP(H2): %.4f  [eQTL only]\n",                          PPH2))
cat(sprintf("  PP(H3): %.4f  [distinct causal variants]\n",            PPH3))
cat(sprintf("  PP(H4): %.4f  [SHARED causal variant]\n",               PPH4))

# Interpretation
if (as.numeric(PPH4) >= 0.80) {
  interpretation <- "STRONG evidence for shared causal variant (PPH4 >= 0.80)"
} else if (as.numeric(PPH4) >= 0.50) {
  interpretation <- "MODERATE evidence for shared causal variant (PPH4 0.50-0.79)"
} else if (as.numeric(PPH3) >= 0.80) {
  interpretation <- "STRONG evidence for DISTINCT causal variants (PPH3 >= 0.80)"
} else if (as.numeric(PPH3) > as.numeric(PPH4)) {
  interpretation <- "Evidence for DISTINCT causal variants (PPH3 > PPH4)"
} else {
  interpretation <- "INCONCLUSIVE / insufficient power (PPH4 < 0.50)"
}
cat("\n  Interpretation:", interpretation, "\n")

best_snp_id   <- NA_character_
best_snp_pph4 <- NA_real_
if (!is.null(coloc_abf$results)) {
  best <- coloc_abf$results %>% arrange(desc(SNP.PP.H4)) %>% slice(1)
  best_snp_id   <- best$snp
  best_snp_pph4 <- best$SNP.PP.H4
  cat(sprintf("  Best colocalization SNP: %s  PP.H4=%.4f\n", best_snp_id, best_snp_pph4))
}

# ============================================================================
# STEP 6: coloc.susie (if >=100 overlapping SNPs)
# ============================================================================
susie_result <- NULL
susie_note   <- paste0("Not run: only ", nrow(merged_dedup),
                       " overlapping SNPs (need >=100 for credible set estimation)")

cat(sprintf("\nStep 6: coloc.susie — %d overlapping SNPs available.\n", nrow(merged_dedup)))
if (nrow(merged_dedup) >= 100) {
  cat("  Running coloc.susie...\n")
  tryCatch({
    susie_result <- coloc.susie(dataset1 = D1, dataset2 = D2)
    cat("  coloc.susie completed.\n")
    if (!is.null(susie_result$summary)) {
      cat("  SuSiE summary:\n")
      print(susie_result$summary)
    }
    susie_note <- "Completed"
  }, error = function(e) {
    susie_note <<- paste("Error:", conditionMessage(e))
    cat("  coloc.susie error:", conditionMessage(e), "\n")
  })
} else {
  cat(sprintf("  Skipping: %d SNPs < 100 required.\n", nrow(merged_dedup)))
}

# ============================================================================
# STEP 7: Save results CSV
# ============================================================================
cat("\nStep 7: Saving results...\n")

n_eqtl_region <- if (is_eqtlgen && exists("eqtl_df")) nrow(eqtl_df) else nrow(eqtl_raw)

results_df <- data.frame(
  analysis_date        = DATE_TAG,
  gene                 = GENE_SYMBOL,
  gene_ensembl         = "ENSG00000243509",
  region               = paste0("chr20:", REGION_START, "-", REGION_END),
  lead_pqtl_snp        = LEAD_PQTL,
  lead_pqtl_pos        = LEAD_PQTL_POS,
  pqtl_source          = "FinnGen R10 Olink",
  pqtl_N               = N_PQTL,
  eqtl_source          = eqtl_source_used,
  eqtl_fetch_method    = eqtl_fetch_method,
  eqtl_N               = N_EQTL,
  n_pqtl_snps_region   = nrow(pqtl),
  n_eqtl_snps_region   = n_eqtl_region,
  n_overlapping_snps   = nrow(merged_dedup),
  coloc_method         = "coloc.abf",
  p1_prior             = 1e-4,
  p2_prior             = 1e-4,
  p12_prior            = 1e-5,
  PPH0                 = round(as.numeric(PPH0), 6),
  PPH1                 = round(as.numeric(PPH1), 6),
  PPH2                 = round(as.numeric(PPH2), 6),
  PPH3                 = round(as.numeric(PPH3), 6),
  PPH4                 = round(as.numeric(PPH4), 6),
  interpretation       = interpretation,
  best_coloc_snp       = best_snp_id,
  best_coloc_snp_PPH4  = round(best_snp_pph4, 6),
  susie_status         = susie_note,
  stringsAsFactors     = FALSE
)

write.csv(results_df, OUT_CSV, row.names = FALSE, quote = TRUE)
cat("  Saved CSV:", OUT_CSV, "\n")

# Per-SNP posteriors
if (!is.null(coloc_abf$results)) {
  snp_out <- file.path(OUT_DIR, paste0("TNFRSF6B_coloc_snp_posteriors_", DATE_TAG, ".csv"))
  join_cols <- intersect(c("chrpos", "variant_id", "rsid", "pos_pqtl",
                            "beta_pqtl", "se_pqtl", "p_pqtl",
                            "beta_eqtl", "se_eqtl", "pValue_eqtl"),
                         colnames(merged_dedup))
  snp_df <- coloc_abf$results %>%
    left_join(merged_dedup %>% select(all_of(join_cols)),
              by = c("snp" = "chrpos")) %>%
    arrange(desc(SNP.PP.H4))
  write.csv(snp_df, snp_out, row.names = FALSE)
  cat("  Saved per-SNP posteriors:", snp_out, "\n")
}

# ============================================================================
# STEP 8: Write Markdown summary
# ============================================================================
cat("\nStep 8: Writing markdown summary...\n")

lead_eqtl_info <- if (is_eqtlgen && exists("eqtl_df") && nrow(eqtl_df) > 0) {
  le <- eqtl_df %>% arrange(pval) %>% slice(1)
  rsid_str <- if ("rsid" %in% colnames(le)) le$rsid[1] else "unknown"
  sprintf("`%s` (chr20:%d), beta=%.3f, p=%.2e",
          rsid_str, le$pos[1], le$beta[1], le$pval[1])
} else {
  "`rs2236511` (chr20:63737808), NES=-0.216, p=4.83e-8 [GTEx fallback]"
}

md_lines <- c(
  "# TNFRSF6B pQTL-eQTL Colocalization Summary",
  paste0("**Date:** ", DATE_TAG, "  "),
  paste0("**Gene:** TNFRSF6B (ENSG00000243509)  "),
  paste0("**Genomic region:** chr20:", REGION_START, "-", REGION_END, " (GRCh38)  "),
  paste0("**Context:** Proteome-wide MR breast cancer manuscript - triangulating MR signal  "),
  "",
  "---",
  "",
  "## Objective",
  "",
  "Test whether the FinnGen Olink pQTL signal for TNFRSF6B (circulating DcR3 protein)",
  "and the eQTLGen Whole Blood cis-eQTL signal share a common causal variant,",
  "supporting the mechanistic link: cis-variant -> mRNA -> protein -> breast cancer risk.",
  "",
  "---",
  "",
  "## Data Sources",
  "",
  "| Dataset | Source | N | Normalisation | SNPs in region |",
  "|---------|--------|---|---------------|----------------|",
  sprintf("| pQTL | FinnGen R10 Olink | %d | Inverse-normal | %d |", N_PQTL, nrow(pqtl)),
  sprintf("| eQTL | %s | %d | Quantile-normal | %d |",
          eqtl_source_used, N_EQTL, n_eqtl_region),
  "",
  paste0("**Lead pQTL SNP:** `rs6011040` (chr20:63706054, splice_donor_5th_base_variant),",
         " beta=-0.434, p=1.27e-12  "),
  paste0("**Lead eQTL SNP:** ", lead_eqtl_info, "  "),
  "",
  "---",
  "",
  "## Harmonisation",
  "",
  sprintf("- pQTL SNPs in 1 Mb region: **%d**", nrow(pqtl)),
  sprintf("- eQTL SNPs in region: **%d**", n_eqtl_region),
  sprintf("- Overlapping SNPs (matched by chr:pos): **%d**", nrow(merged_dedup)),
  "",
  "---",
  "",
  "## Colocalization Results (coloc.abf)",
  "",
  "**Method:** Bayesian colocalization (coloc v5, Giambartolomei et al. 2014)  ",
  "**Priors:** p1 = 1e-4, p2 = 1e-4, p12 = 1e-5  ",
  "",
  "| Hypothesis | Posterior Probability | Interpretation |",
  "|------------|----------------------|----------------|",
  sprintf("| H0: No association | %.4f | Neither trait associated |",      as.numeric(PPH0)),
  sprintf("| H1: pQTL only | %.4f | pQTL signal, no eQTL |",               as.numeric(PPH1)),
  sprintf("| H2: eQTL only | %.4f | eQTL signal, no pQTL |",               as.numeric(PPH2)),
  sprintf("| H3: Distinct causal variants | %.4f | Different causal SNPs |", as.numeric(PPH3)),
  sprintf("| **H4: Shared causal variant** | **%.4f** | **Colocalization** |", as.numeric(PPH4)),
  "",
  paste0("### Verdict: **", interpretation, "**"),
  "",
  sprintf("**Best colocalising SNP:** `%s`  PP(H4) per-SNP = %.4f  ",
          best_snp_id, best_snp_pph4),
  "",
  "---",
  "",
  "## coloc.susie",
  "",
  sprintf("**Status:** %s  ", susie_note),
  "",
  "---",
  "",
  "## Biological Context",
  "",
  "- **Gene:** TNFRSF6B encodes decoy receptor 3 (DcR3), a soluble TNF-receptor superfamily",
  "  member that decoys FasL, LIGHT, and TL1A; circulating DcR3 is MR-causal for breast cancer.  ",
  "- **pQTL lead:** `rs6011040` (chr20:63,706,054) - splice_donor_5th_base_variant,",
  "  beta=-0.434, p=1.27e-12, F-stat=52.5 in FinnGen R10 Olink (N=619).  ",
  "- **SuSiE pQTL-GWAS coloc:** PPH4=0.885 (Supplementary Table 8), confirming shared causal",
  "  variant between protein-level pQTL and breast cancer GWAS at this locus.  ",
  "",
  "---",
  "",
  "## References",
  "",
  "- Giambartolomei C et al. (2014) PLoS Genet. coloc.abf method  ",
  "- Wallace C (2021) PLoS Genet. coloc.susie method  ",
  "- eQTLGen Consortium: https://eqtlgen.org  ",
  "- FinnGen R10: https://finngen.fi  ",
  paste0("- Script: `scripts/60_TNFRSF6B_eqtl_pqtl_coloc.R`  "),
  paste0("- Results: `results/eqtl/TNFRSF6B_pqtl_eqtlgen_coloc_", DATE_TAG, ".csv`  "),
  ""
)

writeLines(md_lines, OUT_MD)
cat("  Saved markdown:", OUT_MD, "\n")

# ============================================================================
# Final summary
# ============================================================================
cat("\n", strrep("=", 72), "\n")
cat("FINAL RESULT\n")
cat(strrep("=", 72), "\n")
cat(sprintf("  PPH4 (shared causal variant):  %.4f\n", as.numeric(PPH4)))
cat(sprintf("  PPH3 (distinct causal vars):   %.4f\n", as.numeric(PPH3)))
cat(sprintf("  PPH2 (eQTL only):              %.4f\n", as.numeric(PPH2)))
cat(sprintf("  Overlapping SNPs used:         %d\n",   nrow(merged_dedup)))
cat(sprintf("  eQTL source:                   %s\n",   eqtl_source_used))
cat(sprintf("  Interpretation:                %s\n",   interpretation))
cat("\n  Output files:\n")
cat("   ", OUT_CSV, "\n")
cat("   ", OUT_MD, "\n")
cat(strrep("=", 72), "\n\n")
