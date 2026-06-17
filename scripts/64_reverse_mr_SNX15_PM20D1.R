#!/usr/bin/env Rscript
# Script 64: Reverse MR — Breast cancer liability -> SNX15 / PM20D1
# Extracts GWS breast cancer variants in ±1Mb cis window, clumps r²<0.001,
# harmonises against FinnGen pQTL as outcome, runs Wald ratio / IVW
# Date: 2026-05-25

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(ieugwasr)
})

PROJ <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"

# ── Cis windows (±1Mb around lead pQTL SNP) ─────────────────────────────────
windows <- list(
  PM20D1 = list(chr = "1",  start = 204843784, end = 206847765),
  SNX15  = list(chr = "11", start = 64025679,  end = 66025679)
)

# ── Load BCAC GWAS summary stats ─────────────────────────────────────────────
cat("Loading BCAC GWAS...\n")
gwas_file <- file.path(PROJ, "data/cancer_gwas/Breast_GCST90018757.h.tsv.gz")
bcac <- fread(gwas_file, sep = "\t", showProgress = FALSE)
cat("BCAC rows:", nrow(bcac), "\n")
cat("Columns:", paste(names(bcac), collapse=", "), "\n\n")

# Standardise column names
setnames(bcac,
  old = c("chromosome","base_pair_location","effect_allele","other_allele",
          "beta","standard_error","effect_allele_frequency","p_value","variant_id"),
  new = c("chr","pos","ea","oa","beta","se","eaf","pval","variant_id"),
  skip_absent = TRUE)

bcac[, chr := as.character(chr)]
bcac[, pos := as.integer(pos)]

# ── Load FinnGen pQTL outcome data ────────────────────────────────────────────
cat("Loading FinnGen pQTL instruments for outcome...\n")
pqtl_dat <- fread(file.path(PROJ, "data/pqtl/pqtl_instruments.csv"), showProgress = FALSE)
cat("pQTL instruments:", nrow(pqtl_dat), "\n\n")

# ── Function: extract, clump, harmonise, run MR ──────────────────────────────
run_reverse_mr <- function(gene, window) {
  cat(sprintf("\n════════════════════════════════════════\n"))
  cat(sprintf("Reverse MR: Breast cancer -> %s\n", gene))
  cat(sprintf("Cis window: chr%s:%d-%d\n", window$chr, window$start, window$end))
  cat(sprintf("════════════════════════════════════════\n"))

  # Step 1: Extract GWS SNPs in cis window from BCAC
  cis_snps <- bcac[chr == window$chr &
                     pos >= window$start &
                     pos <= window$end &
                     pval < 5e-8]

  cat(sprintf("GWS SNPs (p<5e-8) in cis window: %d\n", nrow(cis_snps)))

  if (nrow(cis_snps) == 0) {
    cat("No GWS breast cancer SNPs in cis window — reverse MR not feasible\n")
    return(data.frame(
      gene = gene, n_gws_snps = 0, n_after_clump = 0,
      method = "not_feasible", b = NA, se = NA, pval = NA,
      OR = NA, OR_lo = NA, OR_hi = NA,
      note = "No GWS breast cancer variants in cis window",
      stringsAsFactors = FALSE
    ))
  }

  # Format as TwoSampleMR exposure data frame
  exp_dat <- data.frame(
    SNP             = paste0(cis_snps$chr, ":", cis_snps$pos, "_",
                             cis_snps$ea, "_", cis_snps$oa),
    beta.exposure   = cis_snps$beta,
    se.exposure     = cis_snps$se,
    effect_allele.exposure = cis_snps$ea,
    other_allele.exposure  = cis_snps$oa,
    eaf.exposure    = cis_snps$eaf,
    pval.exposure   = cis_snps$pval,
    exposure        = "Breast_cancer_liability",
    id.exposure     = "Breast_GCST90018757",
    chr.exposure    = cis_snps$chr,
    pos.exposure    = cis_snps$pos,
    stringsAsFactors = FALSE
  )

  # Step 2: LD clumping r²<0.001, 1000G EUR reference via ieugwasr
  cat(sprintf("Running LD clumping (r²<0.001, 10Mb window)...\n"))

  # Try ieugwasr clumping
  clumped <- tryCatch({
    exp_dat_clumped <- clump_data(exp_dat,
                                   clump_r2 = 0.001,
                                   clump_kb = 10000,
                                   pop = "EUR")
    cat(sprintf("After clumping: %d independent SNPs\n", nrow(exp_dat_clumped)))
    exp_dat_clumped
  }, error = function(e) {
    cat("  Clumping via ieugwasr failed:", conditionMessage(e), "\n")
    cat("  Falling back to lead SNP only (most significant)\n")
    exp_dat[which.min(exp_dat$pval.exposure), ]
  })

  cat(sprintf("Instruments for reverse MR: %d\n", nrow(clumped)))
  if (nrow(clumped) > 0) {
    cat("  SNPs:", paste(clumped$SNP, collapse=", "), "\n")
  }

  if (nrow(clumped) == 0) {
    cat("No instruments after clumping\n")
    return(data.frame(
      gene = gene, n_gws_snps = nrow(cis_snps), n_after_clump = 0,
      method = "not_feasible", b = NA, se = NA, pval = NA,
      OR = NA, OR_lo = NA, OR_hi = NA,
      note = "No instruments after LD clumping",
      stringsAsFactors = FALSE
    ))
  }

  # Step 3: Set up outcome (FinnGen pQTL for this protein as outcome)
  # Extract pQTL SNPs for this gene
  gene_pqtl <- pqtl_dat[grepl(gene, pqtl_dat$protein, ignore.case = TRUE), ]
  cat(sprintf("pQTL rows for %s: %d\n", gene, nrow(gene_pqtl)))

  # Build outcome from pQTL data frame
  outcome_dat <- data.frame(
    SNP              = paste0(gene_pqtl$chr, ":", gene_pqtl$pos, "_",
                              gene_pqtl$effect_allele, "_", gene_pqtl$other_allele),
    beta.outcome     = gene_pqtl$beta,
    se.outcome       = gene_pqtl$se,
    effect_allele.outcome = gene_pqtl$effect_allele,
    other_allele.outcome  = gene_pqtl$other_allele,
    eaf.outcome      = gene_pqtl$eaf,
    pval.outcome     = gene_pqtl$pval,
    outcome          = gene,
    id.outcome       = paste0("FinnGen_", gene),
    samplesize.outcome = 619,
    stringsAsFactors = FALSE
  )

  # Step 4: Harmonise
  # We need matching SNPs between clumped BCAC instruments and pQTL region
  # Since we're in cis, look for any BCAC GWS SNP also present in pQTL data
  # OR: use chr:pos matching

  clumped$chr_pos <- paste0(clumped$chr.exposure, ":", clumped$pos.exposure)

  # Try extracting outcome data from OpenGWAS for the clumped SNPs
  cat("Extracting pQTL outcome data from OpenGWAS for clumped SNPs...\n")

  # Get rsIDs via variant_id column if available
  clumped_snp_ids <- clumped$SNP

  # Try to get outcome via FinnGen pQTL — extract from regional data
  # Use the harmonised breast cancer data to find matching pQTL betas
  harm_file <- file.path(PROJ, "data/harmonised/harmonised_protein_Breast_GCST90018757.rds")
  harm_dat <- readRDS(harm_file)
  gene_harm <- harm_dat[grepl(gene, harm_dat$exposure, ignore.case=TRUE), ]

  cat(sprintf("pQTL instruments in harmonised data: %d\n", nrow(gene_harm)))

  if (nrow(gene_harm) == 0) {
    cat("No harmonised pQTL data for", gene, "\n")
    return(data.frame(
      gene = gene, n_gws_snps = nrow(cis_snps), n_after_clump = nrow(clumped),
      method = "not_feasible", b = NA, se = NA, pval = NA,
      OR = NA, OR_lo = NA, OR_hi = NA,
      note = "No pQTL outcome data available",
      stringsAsFactors = FALSE
    ))
  }

  # Build full outcome dataset from FinnGen pQTL instruments file
  # (all SNPs in cis region, not just the instruments)
  # Check if we have regional pQTL sumstats
  pqtl_regional <- file.path(PROJ, "data/pqtl/finngen_olink_group_report")

  # Use the cis_snps positions to find if any BCAC GWS SNPs have pQTL data
  # Match by chr:pos between clumped BCAC SNPs and pQTL instruments
  pqtl_all <- pqtl_dat[grepl(gene, pqtl_dat$protein, ignore.case=TRUE), ]

  if (nrow(pqtl_all) == 0) {
    cat("No pQTL data for", gene, "in instruments file\n")
    return(data.frame(
      gene = gene, n_gws_snps = nrow(cis_snps), n_after_clump = nrow(clumped),
      method = "not_feasible", b = NA, se = NA, pval = NA,
      OR = NA, OR_lo = NA, OR_hi = NA,
      note = "No regional pQTL sumstats available",
      stringsAsFactors = FALSE
    ))
  }

  # Match clumped BCAC SNPs to pQTL by chr:pos
  pqtl_all$chr_pos <- paste0(pqtl_all$chr, ":", pqtl_all$pos)
  matched <- merge(
    clumped[, c("SNP","chr_pos","beta.exposure","se.exposure",
                "effect_allele.exposure","other_allele.exposure",
                "eaf.exposure","pval.exposure")],
    pqtl_all[, c("chr_pos","beta","se","effect_allele","other_allele","eaf")],
    by = "chr_pos"
  )

  cat(sprintf("Clumped BCAC SNPs matching pQTL by chr:pos: %d\n", nrow(matched)))

  if (nrow(matched) > 0) {
    # Manual Wald ratio / IVW from matched data
    # exposure = breast cancer, outcome = protein level
    # Need to harmonise alleles
    matched$flip <- matched$effect_allele.exposure != matched$effect_allele
    matched$beta_outcome_harm <- ifelse(matched$flip, -matched$beta, matched$beta)

    if (nrow(matched) == 1) {
      b   <- matched$beta_outcome_harm / matched$beta.exposure
      se  <- abs(matched$se / matched$beta.exposure)
      p   <- 2 * pnorm(-abs(b/se))
      meth <- "Wald ratio"
    } else {
      # IVW
      w    <- 1 / (matched$se^2)
      b    <- sum(w * (matched$beta_outcome_harm / matched$beta.exposure)) / sum(w)
      se   <- sqrt(1 / sum(w))
      p    <- 2 * pnorm(-abs(b/se))
      meth <- "IVW"
    }

    cat(sprintf("Result: %s  b=%.4f SE=%.4f p=%.4e OR=%.4f\n",
                meth, b, se, p, exp(b)))

    return(data.frame(
      gene = gene, n_gws_snps = nrow(cis_snps), n_after_clump = nrow(clumped),
      n_matched_pqtl = nrow(matched),
      method = meth, b = round(b,4), se = round(se,4), pval = round(p,6),
      OR = round(exp(b),4), OR_lo = round(exp(b-1.96*se),4),
      OR_hi = round(exp(b+1.96*se),4),
      note = "Reverse MR: breast cancer liability -> protein level",
      stringsAsFactors = FALSE
    ))
  } else {
    # No matched SNPs by position — instruments are different
    # This means BCAC GWS SNPs in cis window are not in the pQTL dataset
    # This is expected: pQTL instruments are only the lead/clumped pQTL hits
    # For proper reverse MR, we'd need the full regional pQTL sumstats
    cat("No positional overlap between clumped BCAC SNPs and pQTL instruments\n")
    cat("Need full regional pQTL summary statistics for proper reverse MR\n")

    return(data.frame(
      gene = gene, n_gws_snps = nrow(cis_snps), n_after_clump = nrow(clumped),
      n_matched_pqtl = 0,
      method = "not_feasible",
      b = NA, se = NA, pval = NA, OR = NA, OR_lo = NA, OR_hi = NA,
      note = paste0("No positional overlap between ", nrow(clumped),
                    " clumped BCAC SNPs and pQTL instrument SNPs; ",
                    "full regional pQTL sumstats needed"),
      stringsAsFactors = FALSE
    ))
  }
}

# ── Run for both genes ────────────────────────────────────────────────────────
results <- lapply(names(windows), function(g) {
  run_reverse_mr(g, windows[[g]])
})

results_df <- do.call(rbind, results)

# ── Save ─────────────────────────────────────────────────────────────────────
out_path <- file.path(PROJ, "results/tables/STable_ReverseMR_SNX15_PM20D1_2026-05-25.csv")
write.csv(results_df, out_path, row.names = FALSE)
cat("\n\nResults saved to:", out_path, "\n\n")

# ── Print summary ─────────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════\n")
cat("REVERSE MR SUMMARY: Breast cancer liability -> protein level\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
for (i in seq_len(nrow(results_df))) {
  r <- results_df[i,]
  cat(sprintf("%-8s: GWS_in_window=%d  clumped=%d  matched_pQTL=%s\n",
              r$gene, r$n_gws_snps, r$n_after_clump,
              ifelse(is.na(r$n_matched_pqtl), "N/A", r$n_matched_pqtl)))
  if (!is.na(r$pval)) {
    cat(sprintf("          Method=%s  OR=%.4f [%.4f-%.4f]  p=%.4e\n",
                r$method, r$OR, r$OR_lo, r$OR_hi, r$pval))
  } else {
    cat(sprintf("          Status: %s\n", r$note))
  }
  cat("\n")
}
