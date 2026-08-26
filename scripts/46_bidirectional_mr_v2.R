#!/usr/bin/env Rscript
## Script 48 v3: Bidirectional MR (breast cancer liability → protein level)
## Strategy: for each protein, fetch FinnGen pQTL cis data and breast GWAS GWS
##   data, then match on position + allele concordance.
##
## GENOME BUILD: both sides are GRCh38. The BCAC file used here is the GWAS
##   Catalog *harmonised* release (Breast_GCST90018757.h.tsv.gz), in which
##   `base_pair_location` is already lifted to GRCh38 (hm_coordinate_conversion
##   column present). Verified empirically: GRCh38 cis-pQTL instrument positions
##   match `base_pair_location` exactly (e.g. 1:1023525 AGRN rs2710890).
##   An earlier version of this header incorrectly described the GWAS side as
##   hg19 and justified a ±200 bp tolerance as liftover slack. That was wrong:
##   no cross-build tolerance is required, so matching is now exact (tolerance 0).
##
## Proteins covered: those on chr 1,2,3,11,16,19,20 (1kG VCF available)
## Output: results/bidirectional/bidirectional_v2_*.csv
##
## PROVENANCE OF ARCHIVED RESULTS (read before citing results/bidirectional/):
##   The files currently in results/bidirectional/ were generated with the
##   earlier POS_TOL_BP = 200L setting, i.e. BEFORE the tolerance was corrected
##   to 0 above. Because both datasets are GRCh38, tightening the tolerance can
##   only ever REMOVE a match (a variant matched at 1-200 bp offset with
##   concordant alleles), never add one. The archived result is 7 proteins with
##   estimable reverse effects and 5 with no genome-wide significant breast
##   variant in the cis window, of which only ITIH3 reached nominal
##   significance - i.e. the reported "sparse, no consistent evidence"
##   conclusion cannot be strengthened by the tolerance change, only left
##   unchanged or made more sparse.
##   Rerunning is nonetheless advisable so script and output agree. It requires
##   TwoSampleMR, Rsamtools and GenomicRanges plus network access to the FinnGen
##   public tabix endpoint; none of these were available in the environment
##   where the tolerance fix was made.

suppressPackageStartupMessages({
  library(data.table); library(TwoSampleMR)
  library(Rsamtools);  library(GenomicRanges)
})

proj    <- "."
out_dir <- file.path(proj, "results", "bidirectional")

P_THRESH    <- 5e-8
WINDOW_KB   <- 500L
F_THRESH    <- 10
POS_TOL_BP  <- 0L        # both sides are GRCh38 - exact position matching

priority_proteins <- c("SNX15","EFNA1","UMOD","IL34","PM20D1","CGREF1","ATRAID",
                        "ITIH3","TNFRSF6B","SWAP70","INHBB","APOE")

probe_map <- fread(file.path(proj, "data/pqtl/finngen_olink_probe_map.tsv"))
probe_map[, chr_int := suppressWarnings(as.integer(chr2))]
probe_map <- probe_map[geneName %in% priority_proteins]

# ── Helpers ──────────────────────────────────────────────────────────────────
norm_alleles <- function(a1, a2) {
  # return sorted allele pair as c(lo, hi)
  a1 <- toupper(a1); a2 <- toupper(a2)
  ifelse(a1 <= a2, paste(a1,a2,sep="/"), paste(a2,a1,sep="/"))
}

# RC flip
rc <- function(x) chartr("ACGT","TGCA", x)

# Fetch FinnGen pQTL for a region (same URL as coloc scripts)
read_pqtl <- function(prot, chr, start, end) {
  url  <- sprintf(
    "https://storage.googleapis.com/finngen-public-data-r10/omics/proteomics/release_2023_03_02/data/Olink/pQTL/Olink_Batch1_%s.txt.gz",
    prot)
  prm  <- GRanges(as.character(chr), IRanges(start, end))
  lns  <- tryCatch(scanTabix(TabixFile(url), param=prm)[[1]],
                   error=function(e) character(0))
  if (!length(lns)) return(data.table())
  x <- fread(text=paste(lns,collapse="\n"), header=FALSE)
  cols <- c("chr","pos38","variant_id","ref","alt","alt_freq",
            "beta_out","se_out","t_stat","p_out","log10p","n_out")
  if (ncol(x)==11) cols <- cols[-12]
  setnames(x, cols[seq_len(ncol(x))])
  if (!"n_out" %in% names(x)) x[, n_out:=619L]
  x[, chr:=as.integer(chr)]
  x[, allele_pair := norm_alleles(alt, ref)]
  x
}

# ── Load breast GWAS GWS instruments (once) ──────────────────────────────────
message("Loading breast GWAS GWS variants...")
gwas_path <- file.path(proj,"data/cancer_gwas/Breast_GCST90018757.h.tsv.gz")
gw <- fread(cmd=sprintf("gunzip -c '%s'", gwas_path),
            select=c("chromosome","base_pair_location","effect_allele",
                     "other_allele","beta","standard_error",
                     "effect_allele_frequency","p_value","rsid"))
setnames(gw, c("chr","pos38","ea","oa","beta","se","eaf","pval","rsid"))
gw[, ':='(chr   = suppressWarnings(as.integer(chr)),
           pos38 = suppressWarnings(as.integer(pos38)),
           beta  = suppressWarnings(as.numeric(beta)),
           se    = suppressWarnings(as.numeric(se)),
           eaf   = suppressWarnings(as.numeric(eaf)),
           pval  = suppressWarnings(as.numeric(pval)))]
gw <- gw[!is.na(pval) & pval < P_THRESH & !is.na(beta) & !is.na(se) & se > 0]
gw[, ':='(F_stat = beta^2/se^2,
           allele_pair = norm_alleles(ea, oa),
           allele_pair_rc = norm_alleles(rc(ea), rc(oa)))]
gw <- gw[F_stat > F_THRESH]

# Clump per chromosome
clump_simple <- function(dt) {
  setorder(dt, pval)
  keep <- logical(nrow(dt)); kept_pos <- numeric()
  for (i in seq_len(nrow(dt))) {
    if (!length(kept_pos) || all(abs(dt$pos38[i]-kept_pos) > WINDOW_KB*1000L)) {
      keep[i] <- TRUE; kept_pos <- c(kept_pos, dt$pos38[i])
    }
  }
  dt[keep]
}
gw_c <- gw[, clump_simple(.SD), by=chr]
message(sprintf("Clumped GWS instruments: %d across %d chromosomes",
                nrow(gw_c), uniqueN(gw_c$chr)))

# ── Main loop ─────────────────────────────────────────────────────────────────
all_res <- list(); qc_list <- list()

for (i in seq_len(nrow(probe_map))) {
  prot  <- probe_map$geneName[i]
  chr_p <- probe_map$chr_int[i]
  gene_s<- as.integer(probe_map$start[i])
  gene_e<- as.integer(probe_map$end[i])
  p_s   <- max(1L, gene_s - 1000000L)
  p_e   <- gene_e + 1000000L

  message(sprintf("\n── %s (chr%d) ──", prot, chr_p))

  # Breast GWS instruments in the cis window (both GRCh38; ±1 Mb around gene)
  inst <- gw_c[chr==chr_p & pos38>=p_s & pos38<=p_e]
  message(sprintf("  Breast GWS in cis window: %d", nrow(inst)))
  if (!nrow(inst)) {
    qc_list[[prot]] <- data.table(protein=prot,n_inst=0,n_matched=0,
                                   n_harm=0,mr_success=FALSE,note="no_gws_in_cis_window")
    next
  }

  # Fetch pQTL
  pqtl <- tryCatch(read_pqtl(prot, chr_p, p_s, p_e), error=function(e) data.table())
  if (!nrow(pqtl)) {
    qc_list[[prot]] <- data.table(protein=prot,n_inst=nrow(inst),n_matched=0,
                                   n_harm=0,mr_success=FALSE,note="pqtl_fetch_failed")
    next
  }
  message(sprintf("  pQTL variants in region: %d", nrow(pqtl)))

  # Position + allele matching (both GRCh38; POS_TOL_BP = 0 -> exact)
  matched <- list()
  for (j in seq_len(nrow(inst))) {
    pos_j  <- inst$pos38[j]
    ap_j   <- inst$allele_pair[j]
    ap_j_rc<- inst$allele_pair_rc[j]
    hits   <- pqtl[abs(pos38 - pos_j) <= POS_TOL_BP &
                     (allele_pair == ap_j | allele_pair == ap_j_rc)]
    if (!nrow(hits)) next
    hits <- hits[1]   # take closest
    # flip sign if RC match
    if (hits$allele_pair == ap_j_rc) hits[, beta_out := -beta_out]
    row <- cbind(inst[j], hits)
    matched[[length(matched)+1]] <- row
  }

  if (!length(matched)) {
    # Retained as a diagnostic only: with both sides on GRCh38 an exact-position
    # miss means the variant is genuinely absent from the pQTL file, not shifted.
    message("  No exact-position matches; retrying with ±1000bp as diagnostic...")
    for (j in seq_len(nrow(inst))) {
      pos_j  <- inst$pos38[j]
      ap_j   <- inst$allele_pair[j]
      ap_j_rc<- inst$allele_pair_rc[j]
      hits   <- pqtl[abs(pos38 - pos_j) <= 1000L &
                       (allele_pair == ap_j | allele_pair == ap_j_rc)]
      if (!nrow(hits)) next
      hits <- hits[1]
      if (hits$allele_pair == ap_j_rc) hits[, beta_out := -beta_out]
      matched[[length(matched)+1]] <- cbind(inst[j], hits)
    }
  }

  if (!length(matched)) {
    qc_list[[prot]] <- data.table(protein=prot,n_inst=nrow(inst),n_matched=0,
                                   n_harm=0,mr_success=FALSE,note="no_pos_allele_match")
    next
  }
  mg <- rbindlist(matched, fill=TRUE)
  message(sprintf("  Matched variants: %d", nrow(mg)))

  # Harmonise
  harm <- tryCatch(harmonise_data(
    exposure_dat = data.frame(
      SNP                   = mg$rsid,
      beta.exposure         = mg$beta,
      se.exposure           = mg$se,
      effect_allele.exposure= mg$ea,
      other_allele.exposure = mg$oa,
      eaf.exposure          = mg$eaf,
      pval.exposure         = mg$pval,
      samplesize.exposure   = 228951L,
      exposure              = "Breast_cancer_liability",
      id.exposure           = "Breast_GCST90018757"
    ),
    outcome_dat = data.frame(
      SNP                   = mg$rsid,
      beta.outcome          = mg$beta_out,
      se.outcome            = mg$se_out,
      effect_allele.outcome = mg$alt,
      other_allele.outcome  = mg$ref,
      eaf.outcome           = mg$alt_freq,
      pval.outcome          = mg$p_out,
      samplesize.outcome    = mg$n_out,
      outcome               = prot,
      id.outcome            = paste0("FINNGEN_OLINK_",prot)
    ), action=2), error=function(e) NULL)

  if (is.null(harm)) {
    qc_list[[prot]] <- data.table(protein=prot,n_inst=nrow(inst),n_matched=nrow(mg),
                                   n_harm=0,mr_success=FALSE,note="harmonise_error")
    next
  }
  harm <- as.data.table(harm)[mr_keep==TRUE]
  if (!nrow(harm)) {
    qc_list[[prot]] <- data.table(protein=prot,n_inst=nrow(inst),n_matched=nrow(mg),
                                   n_harm=0,mr_success=FALSE,note="all_harmonised_excluded")
    next
  }

  methods <- if (nrow(harm)>=3) c("mr_ivw","mr_weighted_median","mr_egger_regression") else if (nrow(harm)==2) c("mr_ivw","mr_weighted_median") else "mr_wald_ratio"
  res <- tryCatch(mr(as.data.frame(harm), method_list=methods), error=function(e) NULL)

  if (is.null(res)||!nrow(res)) {
    qc_list[[prot]] <- data.table(protein=prot,n_inst=nrow(inst),n_matched=nrow(mg),
                                   n_harm=nrow(harm),mr_success=FALSE,note="mr_failed")
    next
  }

  res <- as.data.table(res)
  res[, ':='(protein=prot, direction="breast_to_protein",
             or=exp(b), or_lci=exp(b-1.96*se), or_uci=exp(b+1.96*se))]
  all_res[[prot]] <- res
  qc_list[[prot]] <- data.table(protein=prot,n_inst=nrow(inst),n_matched=nrow(mg),
                                 n_harm=nrow(harm),mr_success=TRUE,note="ok")
  message(sprintf("  ✓ %d instruments harmonised; methods: %s",
                  nrow(harm), paste(res$method, collapse=", ")))
  print(res[, .(protein, method, nsnp, b, se, pval, or)])
}

# ── Save ─────────────────────────────────────────────────────────────────────
qc_dt <- rbindlist(qc_list, fill=TRUE)
fwrite(qc_dt, file.path(out_dir, "bidirectional_v2_qc.csv"))
cat("\n=== QC ===\n"); print(qc_dt)

if (length(all_res)) {
  out_dt <- rbindlist(all_res, fill=TRUE)
  fwrite(out_dt, file.path(out_dir, "bidirectional_v2_results.csv"))
  cat("\n=== Results ===\n")
  print(out_dt[, .(protein, method, nsnp, b, se, pval, or, or_lci, or_uci)])
  message("✓ Bidirectional MR complete")
} else {
  message("No results — check QC")
}
