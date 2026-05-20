## Script 34: UKB-PPP Replication MR
## Purpose: Replicate FinnGen-MR hits using UKB-PPP pQTL summary stats (N≈34,000)
## Run incrementally as tar files are downloaded from Synapse (UKB-PPP OLINK NGS)
##
## UKB-PPP file format (regenie output):
##   CHROM GENPOS(GRCh38-bp) ID(chr:GRCh37pos:A0:A1:imp:v1) ALLELE0 ALLELE1
##   A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
##   ALLELE1 = effect allele; BETA = effect of ALLELE1
##   GRCh37 bp position is embedded in the ID field (field 2 when split by ":")
##
## Breast GWAS (GCST90018757): GRCh37, tab-sep
##
## Status [2026-05-07]:
##   APOE: DOWNLOADED — but lead instrument rs429358 (chr19:45411941) absent
##          from breast GWAS (likely excluded in GWAS QC due to complex APOE
##          haplotype region). Best proxy has LOG10P=5.19 in UKB-PPP — too weak.
##          APOE UKB-PPP replication blocked by GWAS coverage, not pQTL quality.
##   Other proteins: pending download from Synapse

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

set.seed(2026)
proj    <- "."
out_dir <- file.path(proj, "results/replication")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ─────────────────────────────────────────────────────────────────────────────
## Protein manifest: name → tar file, chromosomes, GRCh38 instrument position
## (cis window will be ±2 Mb around the gene, using GRCh38 as centre for
##  chromosome extraction; actual GRCh37 positions come from UKB-PPP ID field)
##
## Tar file naming pattern on Synapse:
##   GENE_UNIPROT_OID_v1_Panel.tar
##   e.g. APOE_P02649_OID30727_v1_Inflammation_II.tar
## ─────────────────────────────────────────────────────────────────────────────
protein_manifest <- tribble(
  ~protein,    ~tar_file,                                                    ~chr, ~pos_hg38_approx, ~finngen_cancer,
  "APOE",      "APOE_P02649_OID30727_v1_Inflammation_II.tar",               19,   44908684,          "Breast",
  ## ── Below: pending download ────────────────────────────────────────────────
  ## Synapse project: UKB-PPP (PPP3) → OLINK NGS → search for gene name in
  ## file list. Panel names vary (Inflammation_II, Cardiovascular_II, etc.)
  "SNX15",     "SNX15_Q9Y5X2_OID30xxx_v1_Panel.tar",                        11,   65025679,          "Breast",
  "EFNA1",     "EFNA1_P21709_OID30xxx_v1_Panel.tar",                         1,  155158049,          "Breast",
  "IL34",      "IL34_Q6ZMJ4_OID30xxx_v1_Panel.tar",                         16,   70660097,          "Breast",
  "PM20D1",    "PM20D1_Q9UKU1_OID30xxx_v1_Panel.tar",                        1,  205847765,          "Breast",
  "TNFRSF6B",  "TNFRSF6B_O95407_OID30xxx_v1_Panel.tar",                    20,   63706054,          "Breast",
  "FGF5",      "FGF5_P12034_OID30xxx_v1_Panel.tar",                          4,   80263187,          "Breast",
  "KLB",       "KLB_Q86Z14_OID30xxx_v1_Panel.tar",                           4,   39448609,          "Breast",
  "ITIH3",     "ITIH3_Q06033_OID30xxx_v1_Panel.tar",                         3,   52775509,          "Breast",
  "SWAP70",    "SWAP70_Q9UH65_OID30xxx_v1_Panel.tar",                       11,    9753544,          "Breast",
  "INHBB",     "INHBB_P09529_OID30xxx_v1_Panel.tar",                         2,  120548864,          "Breast",
  "UMOD",      "UMOD_P07911_OID30xxx_v1_Panel.tar",                         16,   20344600,          "Breast",
  "CGREF1",    "CGREF1_Q9Y382_OID30xxx_v1_Panel.tar",                        2,   27101856,          "Breast",
  "ATRAID",    "ATRAID_Q6UXG2_OID30xxx_v1_Panel.tar",                        2,   27163931,          "Breast",
  "TSPAN8",    "TSPAN8_P19075_OID30xxx_v1_Panel.tar",                       12,   71143549,          "Breast",
  "FGFR4",     "FGFR4_P22455_OID30xxx_v1_Panel.tar",                         5,  177093242,          "Breast",
  "ABO",       "ABO_P16442_OID30xxx_v1_Panel.tar",                           9,  133261703,          "Endometrial"
)

## FinnGen reference hits
finngen_hits <- read_csv(
  file.path(proj, "results/tables/STable2_17_FDR_hits_complete.csv"),
  show_col_types = FALSE
) |>
  select(protein, cancer,
         finngen_beta = beta, finngen_se = se,
         finngen_or   = OR, finngen_or_lci = OR_95CI_lower, finngen_or_uci = OR_95CI_upper,
         finngen_fdr  = FDR)

## ── GWAS loader ──────────────────────────────────────────────────────────────
gwas_files <- list(
  Breast      = file.path(proj, "data/cancer_gwas/Breast_GCST90018757.h.tsv.gz"),
  Endometrial = file.path(proj, "data/cancer_gwas/Endometrial_GCST006464.h.tsv.gz")
)

load_gwas_region <- function(cancer, chrom, start_hg19, end_hg19) {
  gf <- gwas_files[[cancer]]
  if (!file.exists(gf)) stop("GWAS file not found: ", gf)
  read_tsv(gf, show_col_types = FALSE) |>
    filter(chromosome == chrom,
           base_pair_location >= start_hg19,
           base_pair_location <= end_hg19) |>
    select(pos_hg19 = base_pair_location,
           ea_gwas  = effect_allele,
           oa_gwas  = other_allele,
           beta_gwas    = beta,
           se_gwas      = standard_error,
           pval_gwas    = p_value,
           rsid)
}

## ── UKB-PPP extractor ────────────────────────────────────────────────────────
extract_ukbppp_cis <- function(tar_path, chrom, pos_hg38, cis_mb = 1.5) {
  chr_str <- as.character(chrom)
  tf      <- untar(tar_path, list = TRUE)
  gz_file <- tf[grep(paste0("_chr", chr_str, "_"), tf)]
  if (length(gz_file) == 0) stop("No chr", chr_str, " gz in tar")
  gz_file <- gz_file[1]
  tmp_dir <- tempdir()
  untar(tar_path, files = gz_file, exdir = tmp_dir)
  tmp_gz  <- file.path(tmp_dir, gz_file)
  on.exit(unlink(tmp_gz), add = TRUE)

  dat <- read_table(tmp_gz, show_col_types = FALSE) |>
    ## Parse GRCh37 position from ID field
    mutate(
      pos_hg19 = as.integer(str_split_fixed(ID, ":", 6)[,2]),
      ea = ALLELE1,   # effect allele
      oa = ALLELE0    # other allele
    ) |>
    ## Cis window: ±cis_mb Mb around GRCh38 centre (generous for liftover diff)
    filter(
      GENPOS >= (pos_hg38 - cis_mb * 1e6),
      GENPOS <= (pos_hg38 + cis_mb * 1e6),
      INFO > 0.8,
      !is.na(LOG10P)
    ) |>
    select(pos_hg19, ea, oa,
           eaf  = A1FREQ, info = INFO, N = N,
           beta_pqtl = BETA, se_pqtl = SE, log10p_pqtl = LOG10P)
  dat
}

## ── Allele harmoniser ────────────────────────────────────────────────────────
flip_map <- c(A="T",T="A",C="G",G="C")
harmonise <- function(pqtl_row, gwas_dat) {
  ea_u  <- toupper(pqtl_row$ea)
  oa_u  <- toupper(pqtl_row$oa)
  ea_fl <- flip_map[ea_u]
  oa_fl <- flip_map[oa_u]
  palin <- (ea_u == flip_map[oa_u])

  hit <- gwas_dat |>
    filter(pos_hg19 == pqtl_row$pos_hg19) |>
    mutate(
      ea_g = toupper(ea_gwas), oa_g = toupper(oa_gwas),
      d  = (ea_u  == ea_g  & oa_u  == oa_g),
      s  = (ea_u  == oa_g  & oa_u  == ea_g),
      fl = (ea_fl == ea_g  & oa_fl == oa_g),
      fs = (ea_fl == oa_g  & oa_fl == ea_g)
    ) |>
    filter(d | s | fl | fs) |>
    slice(1) |>
    mutate(beta_gwas_h = ifelse(d | fl, beta_gwas, -beta_gwas))

  if (nrow(hit) == 0) return(NULL)
  if (palin && (pqtl_row$eaf < 0.05 | pqtl_row$eaf > 0.95)) return(NULL)
  hit
}

## ── Wald ratio + delta-method SE ─────────────────────────────────────────────
wald_ratio <- function(b_exp, se_exp, b_out, se_out) {
  b    <- b_out / b_exp
  se   <- sqrt(se_out^2/b_exp^2 + b_out^2 * se_exp^2 / b_exp^4)
  pval <- 2 * pnorm(abs(b/se), lower.tail = FALSE)
  tibble(beta_mr = b, se_mr = se, pval_mr = pval,
         or_mr = exp(b), or_lci = exp(b - 1.96*se), or_uci = exp(b + 1.96*se))
}

## ═════════════════════════════════════════════════════════════════════════════
## Main loop
## ═════════════════════════════════════════════════════════════════════════════
results <- list()

for (i in seq_len(nrow(protein_manifest))) {
  pm   <- protein_manifest[i,]
  prot <- pm$protein
  tar_full <- file.path(proj, pm$tar_file)

  if (!file.exists(tar_full)) {
    message(sprintf("⏳  %s: tar not yet downloaded — skipping", prot))
    next
  }
  message(sprintf("\n═══ %s ═══", prot))

  ## Extract cis-pQTLs
  pqtl <- tryCatch(
    extract_ukbppp_cis(tar_full, pm$chr, pm$pos_hg38_approx),
    error = function(e) { message("  ERR: ", e$message); NULL })
  if (is.null(pqtl) || nrow(pqtl) == 0) { message("  No cis variants"); next }
  message(sprintf("  %d cis variants (INFO>0.8)", nrow(pqtl)))

  ## Load breast/endometrial GWAS region
  ## Use generous ±3 Mb to cover GRCh37/38 coordinate differences (esp chr4, chr20)
  g37_approx <- pm$pos_hg38_approx  # approximate; real hg19 from ID field
  gwas_reg <- tryCatch(
    load_gwas_region(pm$finngen_cancer, pm$chr, g37_approx - 3e6, g37_approx + 3e6),
    error = function(e) { message("  ERR GWAS: ", e$message); NULL })
  if (is.null(gwas_reg) || nrow(gwas_reg) == 0) { message("  No GWAS variants"); next }
  message(sprintf("  %d GWAS variants in window", nrow(gwas_reg)))

  ## Full inner-join: find ALL overlapping positions between pQTL and GWAS,
  ## then pick the best by pQTL signal.  This handles cases where the lead
  ## pQTL (e.g. rs429358 for APOE) is absent from the GWAS.
  flip_map <- c(A="T",T="A",C="G",G="C")
  merged_all <- inner_join(pqtl, gwas_reg, by = "pos_hg19") |>
    mutate(
      ea_u = toupper(ea), oa_u = toupper(oa),
      ea_g = toupper(ea_gwas), oa_g = toupper(oa_gwas),
      ea_f = flip_map[ea_u], oa_f = flip_map[oa_u],
      match_d  = (ea_u == ea_g  & oa_u == oa_g),
      match_s  = (ea_u == oa_g  & oa_u == ea_g),
      match_f  = (!is.na(ea_f) & ea_f == ea_g  & oa_f == oa_g),
      match_fs = (!is.na(ea_f) & ea_f == oa_g  & oa_f == ea_g),
      palin    = (ea_u == flip_map[oa_u])
    ) |>
    filter(!palin | (eaf > 0.1 & eaf < 0.9)) |>
    filter(match_d | match_s | match_f | match_fs) |>
    mutate(beta_gwas_h = ifelse(match_d | match_f, beta_gwas, -beta_gwas),
           se_gwas_h   = se_gwas)

  if (nrow(merged_all) == 0) {
    message(sprintf("  ✗ No position-overlapping variants found in GWAS"))
    message(sprintf("    Lead pQTL LOG10P=%.1f at chr%d:%d",
                    max(pqtl$log10p_pqtl, na.rm=TRUE), pm$chr,
                    pqtl$pos_hg19[which.max(pqtl$log10p_pqtl)]))
    results[[prot]] <- tibble(
      protein = prot, cancer = pm$finngen_cancer,
      status  = "GWAS_MISSING",
      note    = sprintf("Lead pQTL LOG10P=%.1f not in GWAS; %d cis-pQTLs tested",
                        max(pqtl$log10p_pqtl, na.rm=TRUE), nrow(pqtl))
    )
    next
  }

  ## Select best overlapping variant by pQTL signal
  best_merged <- slice_max(merged_all, log10p_pqtl, n = 1)
  lead_log10p <- max(pqtl$log10p_pqtl, na.rm = TRUE)
  is_lead     <- abs(best_merged$log10p_pqtl - lead_log10p) < 0.5
  used_snp_rank <- which(pqtl$log10p_pqtl[order(-pqtl$log10p_pqtl)] == best_merged$log10p_pqtl[1])[1]
  harm_hit <- list(pqtl = best_merged, gwas = best_merged)

  pq <- harm_hit$pqtl
  gw <- harm_hit$gwas
  mr <- wald_ratio(pq$beta_pqtl, pq$se_pqtl, gw$beta_gwas_h, gw$se_gwas_h)
  fstat <- (pq$beta_pqtl / pq$se_pqtl)^2
  ## Status: LEAD if the best overlapping variant IS the top pQTL; PROXY otherwise
  status_label <- ifelse(is_lead, "LEAD",
                         sprintf("PROXY_rank%d_LOG10P%.1f_of_%.1f",
                                 used_snp_rank, pq$log10p_pqtl, lead_log10p))

  ## FinnGen cross-ref
  fg <- finngen_hits |> filter(protein == prot) |> slice(1)
  dir_ok <- if (nrow(fg) > 0) sign(mr$beta_mr) == sign(fg$finngen_beta) else NA

  row <- tibble(
    protein            = prot,
    cancer             = pm$finngen_cancer,
    status             = status_label,
    snp_id             = sprintf("chr%d:%d:%s:%s", pm$chr, pq$pos_hg19[1], pq$ea[1], pq$oa[1]),
    rsid               = pq$rsid[1],
    pos_hg19           = pq$pos_hg19[1],
    ea_ukb             = pq$ea[1],
    eaf_ukb            = round(pq$eaf[1], 4),
    beta_pqtl_ukb      = round(pq$beta_pqtl[1], 5),
    se_pqtl_ukb        = round(pq$se_pqtl[1], 6),
    log10p_pqtl_ukb    = round(pq$log10p_pqtl[1], 2),
    lead_log10p_pqtl   = round(lead_log10p, 2),
    fstat_ukb          = round(fstat, 1),
    N_ukb              = pq$N[1],
    beta_gwas_harm     = round(gw$beta_gwas_h[1], 6),
    se_gwas            = round(gw$se_gwas_h[1], 6),
    pval_gwas          = pq$pval_gwas[1],
    beta_mr_ukb        = round(mr$beta_mr, 5),
    se_mr_ukb          = round(mr$se_mr, 5),
    pval_mr_ukb        = signif(mr$pval_mr, 3),
    or_mr_ukb          = round(mr$or_mr, 4),
    or_lci_ukb         = round(mr$or_lci, 4),
    or_uci_ukb         = round(mr$or_uci, 4),
    finngen_or         = if (nrow(fg)>0) round(fg$finngen_or, 4) else NA,
    finngen_or_lci     = if (nrow(fg)>0) round(fg$finngen_or_lci, 4) else NA,
    finngen_or_uci     = if (nrow(fg)>0) round(fg$finngen_or_uci, 4) else NA,
    finngen_fdr        = if (nrow(fg)>0) signif(fg$finngen_fdr, 3) else NA,
    direction_consistent = dir_ok
  )
  results[[prot]] <- row

  message(sprintf("  ✓ [%s] chr%d:%d  pQTL LOG10P=%.1f  Fstat=%.0f",
                  row$status, pm$chr, pq$pos_hg19, pq$log10p_pqtl, fstat))
  message(sprintf("    UKB-PPP MR: OR=%.4f (%.4f–%.4f), p=%.3g",
                  mr$or_mr, mr$or_lci, mr$or_uci, mr$pval_mr))
  if (nrow(fg) > 0)
    message(sprintf("    FinnGen MR: OR=%.4f (%.4f–%.4f)",
                    fg$finngen_or, fg$finngen_or_lci, fg$finngen_or_uci))
  message(sprintf("    Direction consistent: %s", dir_ok))
}

## ── Save ─────────────────────────────────────────────────────────────────────
if (length(results) > 0) {
  final <- bind_rows(results)
  out_path <- file.path(out_dir, "ukbppp_replication_mr.csv")
  write_csv(final, out_path)
  message(sprintf("\n✓ Saved %d rows → %s", nrow(final), out_path))

  cat("\n═══ REPLICATION SUMMARY ═══\n")
  success_rows <- final[!is.na(final$or_mr_ukb), ]
  if (nrow(success_rows) > 0) {
    print(success_rows[, intersect(c("protein","status","finngen_or","or_mr_ukb",
                                     "or_lci_ukb","or_uci_ukb","pval_mr_ukb",
                                     "fstat_ukb","direction_consistent"),
                                   names(success_rows))])
  }
  cat("\n── Proteins with GWAS coverage issue ──\n")
  if ("note" %in% names(final)) {
    miss_rows <- final[!is.na(final$note), ]
    if (nrow(miss_rows) > 0) print(miss_rows[, c("protein","note")])
  } else {
    cat("(none)\n")
  }
} else {
  message("No proteins processed — check that tar files are downloaded")
}

## ─────────────────────────────────────────────────────────────────────────────
## NOTES ON APOE:
##   rs429358 (APOE ε4, GRCh37 chr19:45411941, LOG10P=2069 in UKB-PPP) is the
##   overwhelmingly dominant cis-pQTL for APOE protein levels (beta=-1.01).
##   This SNP is absent from the breast cancer GWAS (GCST90018757), likely
##   excluded during GWAS QC because the APOE ε2/ε3/ε4 haplotype block is
##   treated as a special region in many GWAS pipelines.
##   Alternative: use ER+ or ER-negative GWAS files, or IEU OpenGWAS
##   breast cancer GWAS datasets that may have rs429358 included.
##   The FinnGen APOE finding (OR=1.018, FDR=0.022) is supported by:
##   (i) MAGMA p=0.004, (ii) observational support in UKB-PPP, (iii) strong
##   coloc with glycine mediation. UKB-PPP pQTL replication blocked by GWAS.
## ─────────────────────────────────────────────────────────────────────────────
