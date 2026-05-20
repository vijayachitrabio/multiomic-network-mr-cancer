#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(rtracklayer)
  library(GenomicRanges)
  library(Rsamtools)
})

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
out_dir <- file.path(project_dir, "results", "bidirectional")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chain_file <- Sys.getenv("CHAIN_FILE", unset = "/private/tmp/hg19ToHg38.over.chain")
if (!file.exists(chain_file)) stop("Missing hg19ToHg38 chain file: ", chain_file)
chain <- import.chain(chain_file)

WINDOW_BP <- as.integer(Sys.getenv("PQTL_WINDOW_BP", "500000"))
P_THRESHOLD <- 5e-8
WINDOW_KB <- 500
F_THRESHOLD <- 10

priority_proteins <- c("SNX15", "EFNA1", "FGF5", "UMOD", "IL34", "PM20D1", "CGREF1", "ATRAID", "ITIH3", "FGFR4")

probe_map <- fread(file.path(project_dir, "data", "pqtl", "finngen_olink_probe_map.tsv"))
probe_map <- probe_map[geneName %in% priority_proteins]

make_variant_key <- function(chr, pos, allele1, allele2) {
  chr <- sub("^chr", "", as.character(chr))
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  a_min <- ifelse(allele1 <= allele2, allele1, allele2)
  a_max <- ifelse(allele1 <= allele2, allele2, allele1)
  paste(chr, pos, a_min, a_max, sep = ":")
}

extract_breast_instruments <- function() {
  path <- file.path(project_dir, "data", "cancer_gwas", "Breast_GCST90018757.h.tsv.gz")
  x <- fread(
    cmd = sprintf("gunzip -c '%s'", path),
    select = c("chromosome", "base_pair_location", "effect_allele", "other_allele",
               "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid")
  )
  setnames(x, c("chr", "pos37", "ea", "oa", "beta", "se", "eaf", "pval", "rsid"))
  x[, `:=`(
    chr = suppressWarnings(as.integer(chr)),
    pos37 = suppressWarnings(as.integer(pos37)),
    beta = suppressWarnings(as.numeric(beta)),
    se = suppressWarnings(as.numeric(se)),
    eaf = suppressWarnings(as.numeric(eaf)),
    pval = suppressWarnings(as.numeric(pval))
  )]
  x <- x[
    !is.na(rsid) & rsid != "" &
      !is.na(chr) & !is.na(pos37) &
      !is.na(beta) & !is.na(se) & se > 0 &
      !is.na(pval)
  ]
  x <- x[pval < P_THRESHOLD]
  setorder(x, chr, pval, pos37)

  kept <- list()
  for (cur_chr in unique(x$chr)) {
    sub <- x[chr == cur_chr]
    keep_idx <- logical(nrow(sub))
    kept_pos <- numeric()
    for (i in seq_len(nrow(sub))) {
      if (length(kept_pos) == 0 || all(abs(sub$pos37[i] - kept_pos) > WINDOW_KB * 1000)) {
        keep_idx[i] <- TRUE
        kept_pos <- c(kept_pos, sub$pos37[i])
      }
    }
    kept[[length(kept) + 1]] <- sub[keep_idx]
  }
  if (length(kept) == 0) return(data.table())
  out <- rbindlist(kept)
  out[, F_stat := beta^2 / se^2]
  out <- out[F_stat > F_THRESHOLD]
  out
}

lift_to_hg38 <- function(dt) {
  if (nrow(dt) == 0) return(data.table())
  gr37 <- GRanges(
    seqnames = paste0("chr", dt$chr),
    ranges = IRanges(start = dt$pos37, width = 1)
  )
  names(gr37) <- seq_len(nrow(dt))
  gr38 <- unlist(liftOver(gr37, chain))
  if (length(gr38) == 0) return(data.table())

  map <- data.table(
    row_id = as.integer(names(gr38)),
    chr38 = sub("^chr", "", as.character(seqnames(gr38))),
    pos38 = start(gr38)
  )
  out <- dt[map$row_id]
  out[, `:=`(chr38 = as.integer(map$chr38), pos38 = map$pos38)]
  out[, variant_key38 := make_variant_key(chr38, pos38, ea, oa)]
  out
}

read_remote_region <- function(prot, chr, start, end) {
  base_url <- "https://storage.googleapis.com/finngen-public-data-r10/omics/proteomics/release_2023_03_02/data/Olink/pQTL"
  url <- sprintf("%s/Olink_Batch1_%s.txt.gz", base_url, prot)
  region <- GRanges(as.character(chr), IRanges(start, end))
  lines <- scanTabix(TabixFile(url), param = region)[[1]]
  if (length(lines) == 0) return(data.table())
  x <- fread(text = paste(lines, collapse = "\n"), header = FALSE)
  cols <- c("chr", "pos38", "variant_id", "ref", "alt", "alt_freq", "beta_outcome", "se_outcome", "t_stat", "p_outcome", "log10_p", "n_outcome")
  if (ncol(x) == 11) cols <- cols[-12]
  setnames(x, cols)
  if (!"n_outcome" %in% names(x)) x[, n_outcome := 619L]
  x[, `:=`(
    chr = as.integer(chr),
    variant_key38 = make_variant_key(chr, pos38, alt, ref)
  )]
  x
}

breast_inst <- extract_breast_instruments()
breast_inst38 <- lift_to_hg38(breast_inst)

all_results <- list()
qc <- list()

for (i in seq_len(nrow(probe_map))) {
  prot <- probe_map$geneName[i]
  chr <- as.integer(probe_map$chr2[i])
  start <- max(1L, as.integer(probe_map$start[i]) - 1000000L)
  end <- as.integer(probe_map$end[i]) + 1000000L

  message(sprintf("Protein %s chr%s:%s-%s", prot, chr, start, end))
  exp_dat <- breast_inst38[chr38 == chr & pos38 >= start & pos38 <= end]

  if (nrow(exp_dat) == 0) {
    qc[[prot]] <- data.table(protein = prot, n_breast_gws_instruments = 0, n_harmonised = 0, mr_success = FALSE, note = "no_breast_gws_in_cis_window")
    next
  }

  out_dat_raw <- tryCatch(read_remote_region(prot, chr, start, end), error = function(e) data.table())
  if (nrow(out_dat_raw) == 0) {
    qc[[prot]] <- data.table(protein = prot, n_breast_gws_instruments = nrow(exp_dat), n_harmonised = 0, mr_success = FALSE, note = "no_remote_pqtl_region_data")
    next
  }

  merged <- merge(exp_dat, out_dat_raw, by = "variant_key38")
  if (nrow(merged) == 0) {
    qc[[prot]] <- data.table(protein = prot, n_breast_gws_instruments = nrow(exp_dat), n_harmonised = 0, mr_success = FALSE, note = "no_variant_key_overlap_after_liftover")
    next
  }

  merged[, `:=`(
    effect_allele.exposure = ea,
    other_allele.exposure = oa,
    beta.exposure = beta,
    se.exposure = se,
    eaf.exposure = eaf,
    pval.exposure = pval,
    samplesize.exposure = 122977 + 105974,
    outcome = prot,
    id.exposure = "Breast_liability",
    exposure = "Breast_cancer_liability",
    SNP = variant_key38,
    effect_allele.outcome = alt,
    other_allele.outcome = ref,
    beta.outcome = beta_outcome,
    se.outcome = se_outcome,
    eaf.outcome = alt_freq,
    pval.outcome = p_outcome,
    samplesize.outcome = n_outcome,
    id.outcome = paste0("FINNGEN_OLINK_", prot)
  )]

  harm <- harmonise_data(
    exposure_dat = as.data.frame(merged[, .(
      SNP, beta.exposure, se.exposure, effect_allele.exposure, other_allele.exposure,
      eaf.exposure, pval.exposure, samplesize.exposure, exposure, id.exposure
    )]),
    outcome_dat = as.data.frame(merged[, .(
      SNP, beta.outcome, se.outcome, effect_allele.outcome, other_allele.outcome,
      eaf.outcome, pval.outcome, samplesize.outcome, outcome, id.outcome
    )]),
    action = 2
  )

  harm <- as.data.table(harm)
  harm <- harm[mr_keep == TRUE]

  if (nrow(harm) == 0) {
    qc[[prot]] <- data.table(protein = prot, n_breast_gws_instruments = nrow(exp_dat), n_harmonised = 0, mr_success = FALSE, note = "harmonised_zero")
    next
  }

  res <- tryCatch(
    mr(harm, method_list = if (nrow(harm) >= 2) c("mr_ivw", "mr_weighted_median", "mr_egger_regression") else "mr_wald_ratio"),
    error = function(e) NULL
  )

  qc[[prot]] <- data.table(
    protein = prot,
    n_breast_gws_instruments = nrow(exp_dat),
    n_variant_overlap = nrow(merged),
    n_harmonised = nrow(harm),
    mr_success = !is.null(res) && nrow(res) > 0,
    note = if (is.null(res) || nrow(res) == 0) "mr_failed_or_empty" else "ok"
  )

  if (!is.null(res) && nrow(res) > 0) {
    res <- as.data.table(res)
    res[, `:=`(
      protein = prot,
      direction = "Breast_liability_to_protein",
      single_snp_only = nsnp == 1
    )]
    all_results[[prot]] <- res
  }
}

qc_dt <- rbindlist(qc, fill = TRUE)
fwrite(qc_dt, file.path(out_dir, "bidirectional_breast_to_protein_qc.csv"))

if (length(all_results) > 0) {
  out <- rbindlist(all_results, fill = TRUE)
  fwrite(out, file.path(out_dir, "bidirectional_breast_to_protein_results.csv"))
  cat("Wrote:\n")
  cat("  results/bidirectional/bidirectional_breast_to_protein_results.csv\n")
  cat("  results/bidirectional/bidirectional_breast_to_protein_qc.csv\n")
  print(out[, .(protein, method, nsnp, b, se, pval)])
} else {
  fwrite(data.table(), file.path(out_dir, "bidirectional_breast_to_protein_results.csv"))
  cat("No bidirectional MR results generated.\n")
}
