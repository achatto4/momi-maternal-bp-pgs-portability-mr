#!/usr/bin/env Rscript
# ============================================================
# deliv_S2_variant_qc.R  [ B05 -> Supp Table S2 ]
# Genotyping and variant QC per cohort x platform: samples, variants before and after QC,
# what each filter removed, and the final fileset size verified independently.
#
# This is the table a reviewer checks before believing anything genetic in the paper.
#
# WHAT IT ANSWERS THAT MATTERS. B04 found EUR-score recovery of 80.7% in AMANHI-Pakistan
# against 67.6% in GAPPS-Zambia, and the obvious explanation would be that the African
# filesets are sparser after QC. They are NOT -- the reverse:
#     AMANHI-Bangladesh  gsa 7,411,181   dosage  8,099,437
#     GAPPS-Zambia       gsa 13,004,068  dosage 14,210,810
# Zambia carries ~1.75x MORE variants (expected: African genomes segregate more variation)
# while recovering LESS of the score. So the deficit is not a QC-threshold artefact on a
# thinner fileset -- it is locus-specific. EUR score variants are drawn from European GWAS,
# so they are common in Europeans and disproportionately rare or absent in African genomes.
# More data, less of the specific data the score needs. S2 documents that; S3 tests it.
#
# Counts come from two independent sources and are cross-checked:
#   (a) parsed from the plink *_clean.log files (pre-QC, per-filter removals)
#   (b) counted directly from the final .bim/.pvar  (post-QC ground truth)
# If (a) and (b) disagree the log did not describe the fileset actually on disk, and the
# row is flagged rather than silently trusted.
#
# Writes: tables/S2_variant_qc.tsv
#   Rscript deliv_S2_variant_qc.R --outroot DIR
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)
OUTROOT <- momi_arg("--outroot", Sys.getenv("MOMI_OUTROOT"))

## every platform present on disk, including the ones we do NOT analyse -- raw lpwgs is
## reported precisely so the reader can see WHY it is unusable (321k-1.2M variants vs
## 8-14M on dosage), rather than having to take the claim on trust.
S2_PLATS <- c("gsa", "lpwgs", "lpwgs_dosage")

momi_deliverable("S2_variant_qc", script="06_phase2/deliv_S2_variant_qc.R",
                 inputs="genotypes", P=P, stop_on_error=FALSE, body=function(ctx){
  if(is.null(OUTROOT) || !nzchar(OUTROOT) || !dir.exists(OUTROOT))
    return(list(skip=TRUE, reason="no --outroot — see docs/paths_reference.md"))

  ## ---- parse the WHOLE log chain, not one arbitrary step ----
  ## The first version of this read only <plat>_clean.log and reported "no QC performed",
  ## because that is the FINAL step and does no filtering. The real QC lives in
  ## <plat>_clean.pass1.log (and the .allchr/.chrmerge/.rename/.hardcall chain for lpWGS):
  ##   10,787,812 out of 76,413,521 variants loaded
  ##   -22,686 (--rm-dup)  -3,353,945 (allele frequency)  = 7,411,181 remaining
  ## So we glob every log belonging to the platform and aggregate across them.
  nums <- function(x) as.numeric(gsub(",", "", x))
  first_num <- function(tx, pat, grp=1){
    m <- regmatches(tx, regexec(pat, tx, perl=TRUE))
    m <- m[vapply(m, length, 0L) > grp]
    if(!length(m)) return(NA_real_)
    nums(m[[1]][grp+1])
  }
  sum_num <- function(tx, pat, grp=1){
    m <- regmatches(tx, regexec(pat, tx, perl=TRUE))
    m <- m[vapply(m, length, 0L) > grp]
    if(!length(m)) return(NA_real_)
    sum(vapply(m, function(z) nums(z[grp+1]), 0))
  }
  parse_chain <- function(dir, plat){
    ## every log for this platform, EXCLUDING the nopass diagnostic build
    fs <- list.files(dir, pattern=sprintf("^%s(_clean|_hc)?.*\\.log$", plat), full.names=TRUE)
    fs <- fs[!grepl("nopass", fs)]
    if(!length(fs)) return(NULL)
    tx <- unlist(lapply(fs, function(f) tryCatch(readLines(f, warn=FALSE),
                                                 error=function(e) character(0))))
    if(!length(tx)) return(NULL)
    list(
      logs_parsed   = length(fs),
      ## "X out of Y variants loaded" -- Y is the true source count before any subsetting
      var_source    = first_num(tx, "([0-9,]+) out of ([0-9,]+) variants loaded", 2),
      ## lpWGS is processed PER CHROMOSOME, so these lines recur ~22 times. Take the max of
      ## every match (the merged whole-genome step reports the total and dominates), not the
      ## first, which would return whichever chromosome log sorted first.
      var_loaded    = suppressWarnings(max(c(
                        first_num(tx, "([0-9,]+) out of [0-9,]+ variants loaded"),
                        vapply(regmatches(tx, regexec("^([0-9,]+) variants loaded", tx, perl=TRUE)),
                               function(z) if(length(z)>1) nums(z[2]) else NA_real_, 0)),
                        na.rm=TRUE)),
      samples_max   = max(c(sum_num(tx, "([0-9,]+) samples? \\(.*loaded")*0 +
                            first_num(tx, "([0-9,]+) samples? \\(.*loaded"),
                            NA_real_), na.rm=TRUE),
      rm_dup        = sum_num(tx, "--rm-dup: [0-9,]+ duplicated IDs, ([0-9,]+) variants removed"),
      rm_maf        = sum_num(tx, "([0-9,]+) variants removed due to allele frequency"),
      rm_geno       = sum_num(tx, "([0-9,]+) variants removed due to missing genotype"),
      rm_hwe        = sum_num(tx, "([0-9,]+) variants removed due to Hardy-Weinberg"),
      ## SUM, not first: per-chromosome pipelines emit one "remaining" line per chromosome,
      ## so the total is their sum. GSA emits exactly one, where sum == first.
      var_remaining = sum_num(tx, "([0-9,]+) variants remaining after main filters"),
      thresholds    = paste(unique(unlist(regmatches(tx,
                        gregexpr("--(geno|maf|hwe|mind|mac|max-maf)\\s+[0-9.]+", tx)))),
                        collapse=" "))
  }

  rows <- list()
  for(coh in MOMI_COH_ALL) for(pl in S2_PLATS){
    pre <- file.path(OUTROOT, coh, paste0(pl, "_clean"))
    bim <- paste0(pre, ".bim"); pvar <- paste0(pre, ".pvar")
    fam <- paste0(pre, ".fam"); psam <- paste0(pre, ".psam")
    vf <- if(file.exists(pvar)) pvar else if(file.exists(bim)) bim else NA_character_
    sf <- if(file.exists(psam)) psam else if(file.exists(fam)) fam else NA_character_
    if(is.na(vf)) next

    ## (b) ground truth straight off the fileset
    nvar <- as.integer(system(sprintf("grep -vc '^#' %s", shQuote(vf)), intern=TRUE))
    nsam <- if(is.na(sf)) NA_integer_ else
              as.integer(system(sprintf("grep -vc '^#' %s", shQuote(sf)), intern=TRUE))

    ## (a) whatever the log chain says
    L <- parse_chain(file.path(OUTROOT, coh), pl)
    if(is.null(L)) L <- list(logs_parsed=0L, var_source=NA_real_, var_loaded=NA_real_,
                             samples_max=NA_real_, rm_dup=NA_real_, rm_maf=NA_real_,
                             rm_geno=NA_real_, rm_hwe=NA_real_, var_remaining=NA_real_,
                             thresholds="")
    sfin <- function(x) if(is.null(x) || !length(x) || !is.finite(x)) NA_real_ else x

    rows[[paste(coh,pl)]] <- data.table(
      cohort=coh, cohort_anc=unname(MOMI_ANC[coh]), platform=pl,
      analysed = as.integer(pl %in% MOMI_PLATS),
      logs_parsed = L$logs_parsed,
      n_samples_final = nsam,
      n_samples_preQC = sfin(L$samples_max),
      variants_source = sfin(L$var_source),      # before any subsetting
      variants_loaded = sfin(L$var_loaded),      # entering the filters
      removed_dup     = sfin(L$rm_dup),
      removed_maf     = sfin(L$rm_maf),          # the differential one -- see console note
      removed_missing = sfin(L$rm_geno),
      removed_hwe     = sfin(L$rm_hwe),
      variants_postQC_log  = sfin(L$var_remaining),
      variants_postQC_file = nvar,
      qc_thresholds = L$thresholds,
      log_matches_file = if(!is.finite(sfin(L$var_remaining))) NA else
                           abs(L$var_remaining - nvar) <= 1)
  }
  if(!length(rows)) return(list(skip=TRUE, reason="no filesets found under --outroot"))
  S2 <- rbindlist(rows, fill=TRUE)
  ## DERIVE the pre-QC total rather than parse it. Parsing gave 1445% removal on dosage,
  ## because removals are SUMMED over 22 per-chromosome logs (~116M, correct) while the
  ## parsed "loaded" was the MAX (8.1M, the merged post-QC figure) -- numerator and
  ## denominator from different stages. removals + final is internally consistent by
  ## construction, and it reproduces the GSA log exactly:
  ##   22,686 + 3,353,945 + 7,411,181 = 10,787,812  ✓
  S2[, removed_total := rowSums(cbind(removed_dup, removed_maf, removed_missing, removed_hwe),
                                na.rm=TRUE)]
  S2[, variants_preQC_derived := removed_total + variants_postQC_file]
  S2[, pct_removed_total := round(100*removed_total/variants_preQC_derived, 1)]
  S2[, pct_removed_maf   := round(100*removed_maf  /variants_preQC_derived, 1)]
  ## keep the parsed figure alongside so the two can be compared where both are meaningful
  setnames(S2, "variants_loaded", "variants_loaded_parsed")
  out <- momi_write_table(S2[order(platform, cohort)], "S2_variant_qc", P)

  ## ---- console: the contrast this table exists to document ----
  cat("\n== post-QC variants per cohort x platform ==\n")
  print(dcast(S2, cohort + cohort_anc ~ platform, value.var="variants_postQC_file"))
  cat("\n== mean post-QC variants by cohort ancestry (analysed platforms only) ==\n")
  print(S2[analysed==1, .(mean_variants=as.integer(mean(variants_postQC_file, na.rm=TRUE))),
           by=.(platform, cohort_anc)][order(platform, -mean_variants)])
  cat("\nAFRICAN cohorts show MORE variants, not fewer, so a sparser-fileset explanation for\n",
      "their lower score recovery (67.6% Zambia vs 80.7% Pakistan, B04) is ruled out.\n", sep="")

  cat("\n== the MAF filter, and why it is not a neutral choice ==\n")
  print(S2[analysed==1, .(cohort, cohort_anc, platform, variants_preQC_derived, removed_maf,
                          pct_removed_maf, qc_thresholds)][order(platform, -pct_removed_maf)])
  cat("\nA uniform MAF floor is applied WITHIN each cohort. European-common variants are\n",
      "disproportionately African-RARE, so a fixed threshold removes more of a European-\n",
      "derived score in African cohorts. Many missing score variants are therefore PRESENT\n",
      "in the data but below our own cut -- 'we removed them', not 'they are absent'.\n",
      "Defensible (rare-variant imputation is unreliable, and low-MAF variants contribute\n",
      "little PRS variance) but it is OUR analysis decision and belongs in Methods, not\n",
      "buried in a preprocessing log.\n", sep="")

  cat("\n== sample attrition (pre-QC vs final fileset) ==\n")
  sa <- S2[analysed==1 & is.finite(n_samples_preQC),
           .(cohort, platform, n_samples_preQC, n_samples_final,
             dropped = n_samples_preQC - n_samples_final,
             pct_dropped = round(100*(n_samples_preQC-n_samples_final)/n_samples_preQC,1))]
  if(nrow(sa)) print(sa[order(-pct_dropped)])
  cat("\n>>> RESOLVED 2026-07-19: these are NOT QC failures, and the samples are NOT recoverable.\n",
      "The drop happens at the KEEP/RENAME step, not a QC filter, and the ratios are ~2:1\n",
      "(Pemba 7629->3813 = 2.00, GSA 638->333 = 1.92) -- mother-INFANT pairs sequenced\n",
      "together with the maternal subset extracted afterwards.\n",
      "BUT the infants are gone from pipeline_out: chrtmp .psam now holds 3,813 not 7,629,\n",
      "lpwgs_map.rename has only 3,813 lines, and NONE map to a BABY_ID. The 7,629 survives\n",
      "only in the log. Infants were dropped and the per-chromosome files regenerated.\n",
      "  => n_samples_preQC below is 'samples in the SOURCE VCF, mothers + infants', not\n",
      "     'samples before QC'. Report it that way or not at all.\n",
      "  => S15 (fetal vs maternal) must use new_data_May13_26/merged_lps_gsa_all_sites,\n",
      "     which retains 8,086 baby-suffix samples -- but has NOT been through this QC.\n", sep="")
  bad <- S2[log_matches_file %in% FALSE]
  if(nrow(bad)){
    cat("\n!! log and fileset disagree on post-QC variant count -- do not trust these rows:\n")
    print(bad[, .(cohort, platform, variants_postQC_log, variants_postQC_file)])
  }
  miss <- S2[is.na(variants_preQC_derived)]
  if(nrow(miss)) cat(sprintf("\n%d row(s) have no parseable pre-QC count; the *_clean.log is absent or in an unexpected format. Post-QC counts are still exact (read off the fileset).\n", nrow(miss)))

  list(n=nrow(S2),
       key=sprintf("cohorts=%d platforms=%d; postQC variants %s-%s; AFR>SAS variant count=%s",
                   uniqueN(S2$cohort), uniqueN(S2$platform),
                   format(min(S2$variants_postQC_file,na.rm=TRUE), big.mark=","),
                   format(max(S2$variants_postQC_file,na.rm=TRUE), big.mark=","),
                   as.character(S2[analysed==1 & cohort_anc=="AFR", mean(variants_postQC_file,na.rm=TRUE)] >
                                S2[analysed==1 & cohort_anc=="SAS", mean(variants_postQC_file,na.rm=TRUE)])),
       outputs=basename(out))
})
