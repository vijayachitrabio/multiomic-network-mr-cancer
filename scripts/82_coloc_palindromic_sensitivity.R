#!/usr/bin/env Rscript
## Script 82: Palindromic-SNP sensitivity analysis for protein-cancer colocalization
## ─────────────────────────────────────────────────────────────────────────────
## WHY THIS EXISTS
## External audit (2026-08-28) flagged the palindromic-variant handling in the
## colocalization harmonisation used by scripts 36/38/39. Those scripts contain:
##
##     filter(!palin | (alt_freq > 0.1 & alt_freq < 0.9))
##
## which KEEPS palindromic (A/T, C/G) variants whose frequency lies in the
## 0.1-0.9 band and DROPS those outside it. That is inverted: strand for a
## palindromic variant can only be inferred when the allele frequency is far
## from 0.5, so the retained band is exactly the ambiguous one. The filter also
## never compares the pQTL frequency against the GWAS frequency, so for every
## palindromic variant both the direct and the flipped allele match tests
## succeed and the code falls through to `beta_g` unflipped -- an unresolved
## orientation, i.e. a coin flip on the sign.
##
## This script re-runs colocalization for all 17 MR-prioritized proteins with
## ALL palindromic variants removed (the conservative choice), and reports the
## result beside the published values so tier stability can be judged.
##
## Everything else -- window, priors, MAF source, case fractions, LD reference,
## SuSiE settings, convergence cascade -- is byte-identical to script 39.
## Package versions match the manuscript: coloc 5.2.3, susieR 0.14.2.
##
## Output: results/validation/protein_coloc_palindromic_sensitivity.csv

suppressPackageStartupMessages({
  library(data.table); library(coloc); library(susieR)
  library(Rsamtools); library(GenomicRanges); library(dplyr); library(readr)
})

set.seed(2026)
proj    <- "."
out_dir <- file.path(proj, "results/validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

WINDOW_BP <- 500000L
MIN_SNPS  <- 50L
BASE_1KG  <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
`%||%` <- function(a,b) if (is.null(a)) b else a

## ── Outcome configs ──────────────────────────────────────────────────────────
mk_gwas <- function(file, n_cases, n_controls) {
  g <- list(file=file.path(proj,"data/cancer_gwas",file),
            n_cases=n_cases, n_controls=n_controls)
  g$n_total <- g$n_cases + g$n_controls
  g$s       <- g$n_cases / g$n_total
  g
}
GWAS <- list(
  Breast      = mk_gwas("Breast_GCST90018757.h.tsv.gz", 122977L, 105974L),
  Endometrial = mk_gwas("Endometrial_GCST006464.h.tsv.gz", 12906L, 108979L)
)

## 17 MR-prioritized proteins. ABO is the endometrial comparator; rest breast.
PROTEINS <- c("EFNA1","TNFRSF6B","ATRAID","FGF5","ABO","SNX15","PM20D1","UMOD",
              "APOE","TSPAN8","IL34","ITIH3","SWAP70","KLB","FGFR4","CGREF1","INHBB")
outcome_of <- function(p) if (p == "ABO") "Endometrial" else "Breast"

## ── 1000G EUR sample index ───────────────────────────────────────────────────
message("Loading 1000G population panel...")
panel   <- fread("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt")
eur_ids <- panel[Superpopulation == "EUR", SampleID]
message(sprintf("  %d EUR samples", length(eur_ids)))

vcf_cache <- list()
get_vcf_info <- function(chr_num) {
  key <- as.character(chr_num)
  if (!is.null(vcf_cache[[key]])) return(vcf_cache[[key]])
  vcf_url <- sprintf("%s/1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz", BASE_1KG, chr_num)
  hdr <- headerTabix(TabixFile(vcf_url))
  all_samples <- strsplit(hdr$header[length(hdr$header)], "\t", fixed=TRUE)[[1]][-(1:9)]
  info <- list(url=vcf_url, eur_col=which(all_samples %in% eur_ids) + 9L)
  info$n_eur <- length(info$eur_col)
  vcf_cache[[key]] <<- info
  info
}

## ── Loaders (identical to script 39) ────────────────────────────────────────
load_pqtl <- function(protein) {
  x <- fread(file.path(proj,"data/pqtl/priority_regions", paste0(protein,"_pqtl_regions.tsv.gz")))
  x[!is.na(beta) & !is.na(se) & se > 0 & !is.na(alt_freq) & alt_freq > 0 & alt_freq < 1]
}
load_gwas <- function(gcfg, chr, pos_min, pos_max) {
  read_tsv(gcfg$file, show_col_types=FALSE,
    col_types=cols_only(chromosome=col_integer(), base_pair_location=col_integer(),
      effect_allele=col_character(), other_allele=col_character(), beta=col_double(),
      standard_error=col_double(), effect_allele_frequency=col_double(), p_value=col_double())) |>
    filter(chromosome==as.integer(chr), base_pair_location>=pos_min, base_pair_location<=pos_max) |>
    rename(pos=base_pair_location, ea=effect_allele, oa=other_allele,
           beta_g=beta, se_g=standard_error, eaf_g=effect_allele_frequency, p_g=p_value) |>
    filter(!is.na(beta_g), !is.na(se_g), se_g>0, !is.na(eaf_g), eaf_g>0, eaf_g<1)
}

## ── Harmonise: DROP ALL PALINDROMIC (the only change vs script 39) ──────────
harmonise <- function(pqtl, gwas, drop_palindromic) {
  flip <- c(A="T",T="A",C="G",G="C")
  h <- inner_join(pqtl |> mutate(pos=as.integer(pos)),
                  gwas |> mutate(pos=as.integer(pos)), by="pos") |>
    mutate(ea_p=toupper(alt), oa_p=toupper(ref), ea_g2=toupper(ea), oa_g2=toupper(oa),
           ea_pf=flip[toupper(alt)], oa_pf=flip[toupper(ref)],
           match_d =ea_p==ea_g2 & oa_p==oa_g2,
           match_s =ea_p==oa_g2 & oa_p==ea_g2,
           match_fl=!is.na(ea_pf) & ea_pf==ea_g2 & oa_pf==oa_g2,
           match_fs=!is.na(ea_pf) & ea_pf==oa_g2 & oa_pf==ea_g2,
           palin   =ea_p==flip[oa_p])
  n_before <- nrow(h); n_palin <- sum(h$palin, na.rm=TRUE)
  n_palin_ambig <- sum(h$palin & h$alt_freq > 0.1 & h$alt_freq < 0.9, na.rm=TRUE)
  h <- if (drop_palindromic) filter(h, !palin) else filter(h, !palin | (alt_freq>0.1 & alt_freq<0.9))
  h <- h |> filter(match_d|match_s|match_fl|match_fs) |>
    mutate(beta_g_h=if_else(match_d|match_fl, beta_g, -beta_g),
           eaf_g_h =if_else(match_d|match_fl, eaf_g, 1-eaf_g))
  attr(h,"qc") <- list(n_before=n_before, n_palin=n_palin,
                       n_palin_ambig=n_palin_ambig, n_after=nrow(h))
  h
}

## ── LD matrix (identical to script 39) ──────────────────────────────────────
build_ld_matrix <- function(chr_num, pos_min, pos_max, eur_col, vcf_url, target_snp_ids) {
  region <- GRanges(paste0("chr",chr_num), IRanges(pos_min, pos_max))
  message(sprintf("  1000G tabix: chr%s:%d-%d", chr_num, pos_min, pos_max))
  lines <- tryCatch(scanTabix(TabixFile(vcf_url), param=region)[[1]],
                    error=function(e){message("  tabix err: ",e$message); character(0)})
  if (!length(lines)) return(NULL)
  keep_ids <- character(0); dos <- list()
  for (ln in lines) {
    f <- strsplit(ln, "\t", fixed=TRUE)[[1]]
    if (length(f) < 10) next
    id <- paste0(f[2],":",toupper(f[4]),":",toupper(f[5]))
    if (!(id %in% target_snp_ids)) next
    gt <- f[eur_col]
    d  <- vapply(gt, function(g){
      a <- substr(g,1,1); b <- substr(g,3,3)
      if (a=="." || b==".") return(NA_real_)
      as.numeric(a)+as.numeric(b)
    }, numeric(1), USE.NAMES=FALSE)
    if (mean(is.na(d)) > 0.05) next
    if (length(unique(d[!is.na(d)])) < 2) next
    keep_ids <- c(keep_ids, id); dos[[length(dos)+1]] <- d
  }
  if (length(keep_ids) < 2) return(NULL)
  M <- do.call(cbind, dos); colnames(M) <- keep_ids
  dup <- duplicated(colnames(M)); if (any(dup)) M <- M[,!dup,drop=FALSE]
  ld <- suppressWarnings(cor(M, use="pairwise.complete.obs"))
  ld[is.na(ld)] <- 0; diag(ld) <- 1
  list(ld=ld, snp_ids=colnames(M))
}

clean_lbf <- function(s){ if(!is.null(s$lbf_variable)){b<-is.na(s$lbf_variable)|is.nan(s$lbf_variable); if(any(b)) s$lbf_variable[b]<-0}; s }
extract_pph4 <- function(smry){
  if (is.null(smry)) return(list(pph4=0, snp=NA_character_, n=0L))
  if (is.numeric(smry) && "PP.H4.abf" %in% names(smry))
    return(list(pph4=as.numeric(smry["PP.H4.abf"]), snp=NA_character_, n=1L))
  if ((is.data.frame(smry)||is.data.table(smry)) && nrow(smry)>0 && "PP.H4.abf" %in% names(smry)){
    b <- which.max(as.numeric(smry[["PP.H4.abf"]]))
    return(list(pph4=as.numeric(smry[["PP.H4.abf"]][b]),
                snp=if("hit1" %in% names(smry)) as.character(smry[["hit1"]][b]) else NA_character_,
                n=nrow(smry)))
  }
  list(pph4=0, snp=NA_character_, n=0L)
}

## ── Main ─────────────────────────────────────────────────────────────────────
results <- list()
for (prot in PROTEINS) {
  message(sprintf("\n══════ %s (%s) ══════", prot, outcome_of(prot)))
  gcfg <- GWAS[[outcome_of(prot)]]
  pq <- tryCatch(load_pqtl(prot), error=function(e){message("  pQTL load err: ",e$message); NULL})
  if (is.null(pq) || !nrow(pq)) next

  lead <- pq[which.min(p)]
  chr <- as.character(lead$chr[1]); lead_pos <- as.integer(lead$pos[1])
  pos_min <- lead_pos - WINDOW_BP; pos_max <- lead_pos + WINDOW_BP
  message(sprintf("  lead chr%s:%d  p=%.3g", chr, lead_pos, lead$p[1]))

  gw <- tryCatch(load_gwas(gcfg, chr, pos_min, pos_max), error=function(e){message("  GWAS err: ",e$message); NULL})
  if (is.null(gw) || !nrow(gw)) next
  pq_w <- pq[pos >= pos_min & pos <= pos_max]

  harm <- harmonise(pq_w, gw, drop_palindromic = TRUE)
  qc <- attr(harm,"qc")
  message(sprintf("  harmonised: %d overlap -> %d after dropping %d palindromic (%d of them in the ambiguous 0.1-0.9 band)",
                  qc$n_before, qc$n_after, qc$n_palin, qc$n_palin_ambig))
  if (nrow(harm) < 10) { message("  too few variants, skipping"); next }

  snp_key <- paste0(harm$pos,":",toupper(harm$ref),":",toupper(harm$alt))
  vi <- tryCatch(get_vcf_info(chr), error=function(e){message("  vcf err: ",e$message); NULL})
  ld_info <- if (!is.null(vi)) build_ld_matrix(chr, pos_min, pos_max, vi$eur_col, vi$url, snp_key) else NULL

  if (!is.null(ld_info)) {
    idx <- match(ld_info$snp_ids, snp_key); ok <- !is.na(idx)
    harm_sub <- harm[idx[ok],]; ld_sub <- ld_info$ld[ok,ok,drop=FALSE]; snp_sub <- ld_info$snp_ids[ok]
  } else { harm_sub <- harm; ld_sub <- NULL; snp_sub <- snp_key }
  n_ok <- nrow(harm_sub)
  message(sprintf("  LD-matched variants: %d", n_ok))

  D1 <- list(beta=harm_sub$beta, varbeta=harm_sub$se^2, snp=snp_sub, type="quant",
             N=619L, MAF=pmin(harm_sub$alt_freq, 1-harm_sub$alt_freq))
  D2 <- list(beta=harm_sub$beta_g_h, varbeta=harm_sub$se_g^2, snp=snp_sub, type="cc",
             N=gcfg$n_total, s=gcfg$s, MAF=pmin(harm_sub$eaf_g_h, 1-harm_sub$eaf_g_h))

  pp <- tryCatch(coloc.abf(D1,D2)$summary, error=function(e){
    message("  ABF err: ",e$message)
    c(PP.H0.abf=NA,PP.H1.abf=NA,PP.H2.abf=NA,PP.H3.abf=NA,PP.H4.abf=NA)})
  message(sprintf("  [ABF] PPH3=%.4f  PPH4=%.4f", pp["PP.H3.abf"], pp["PP.H4.abf"]))

  pph4_s <- NA_real_; best_snp <- NA_character_; n_cs1 <- 0L; n_cs2 <- 0L; n_pairs <- 0L
  if (!is.null(ld_sub) && n_ok >= MIN_SNPS) {
    D1L <- c(D1, list(LD=ld_sub)); D2L <- c(D2, list(LD=ld_sub))
    s1 <- tryCatch(runsusie(D1L, repeat_until_convergence=TRUE, maxit=10000L), error=function(e){message("  D1 err: ",e$message);NULL})
    s2 <- tryCatch(runsusie(D2L, repeat_until_convergence=TRUE, maxit=10000L), error=function(e){message("  D2 err: ",e$message);NULL})
    if (!is.null(s1)) { s1 <- clean_lbf(s1); n_cs1 <- length(s1$sets$cs %||% list()) }
    if (!is.null(s2)) { s2 <- clean_lbf(s2); n_cs2 <- length(s2$sets$cs %||% list()) }
    message(sprintf("  CS: pQTL=%d  GWAS=%d", n_cs1, n_cs2))
    if (!is.null(s1) && !is.null(s2) && n_cs1>0 && n_cs2>0) {
      csr <- tryCatch(coloc.susie(s1,s2), error=function(e) NULL)
      if (is.null(csr)) csr <- tryCatch(coloc.susie(D1L,D2L,
        runsusie.args=list(repeat_until_convergence=TRUE,maxit=10000L)), error=function(e) NULL)
      if (is.null(csr)) csr <- tryCatch({
        bf1 <- s1$lbf_variable[s1$sets$cs_index,,drop=FALSE]
        bf2 <- s2$lbf_variable[s2$sets$cs_index,,drop=FALSE]
        bf1[is.nan(bf1)] <- 0; bf2[is.nan(bf2)] <- 0
        coloc:::coloc.bf_bf(bf1,bf2)}, error=function(e) NULL)
      if (!is.null(csr)) { ex <- extract_pph4(csr$summary); pph4_s <- ex$pph4; best_snp <- ex$snp; n_pairs <- ex$n }
    }
    message(sprintf("  [SuSiE] PPH4=%s", ifelse(is.na(pph4_s),"NA",sprintf("%.4f",pph4_s))))
  } else message("  SuSiE skipped (insufficient LD-matched variants)")

  results[[prot]] <- data.table(
    protein=prot, outcome=outcome_of(prot), lead_chr=chr, lead_pos=lead_pos,
    n_overlap=qc$n_before, n_palindromic_dropped=qc$n_palin,
    n_palindromic_in_ambiguous_band=qc$n_palin_ambig, n_harmonised=qc$n_after,
    n_ld_matched=n_ok,
    PPH3_abf_nopalin=as.numeric(pp["PP.H3.abf"]), PPH4_abf_nopalin=as.numeric(pp["PP.H4.abf"]),
    PPH4_susie_nopalin=pph4_s, susie_best_snp_nopalin=best_snp,
    n_cs_pqtl=n_cs1, n_cs_gwas=n_cs2, n_coloc_pairs=n_pairs)
  fwrite(rbindlist(results, fill=TRUE), file.path(out_dir,"protein_coloc_palindromic_sensitivity.csv"))
}

out <- rbindlist(results, fill=TRUE)
fwrite(out, file.path(out_dir,"protein_coloc_palindromic_sensitivity.csv"))
message(sprintf("\n✓ %d proteins -> %s", nrow(out), file.path(out_dir,"protein_coloc_palindromic_sensitivity.csv")))
print(out[, .(protein, n_palindromic_dropped, PPH3_abf_nopalin, PPH4_abf_nopalin, PPH4_susie_nopalin)])
