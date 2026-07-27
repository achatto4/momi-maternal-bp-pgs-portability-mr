#!/usr/bin/env Rscript
# ============================================================
# deliv_S1_prs_manifest.R  [ B04 -> Supp Table S1 (+ S1b) ]
# Polygenic-score panel provenance. One row per PGS in MOMI_PANEL:
#   id, trait, ancestry, source, role, whether scored on our platforms, the PUBLISHED
#   variant count, the variants USED in our data, and the resulting SCORE RECOVERY RATE.
#
# WHY RECOVERY IS THE POINT OF THIS TABLE: naming the source score is not enough for a
# portability paper. A 7.36M-variant published score of which we recover 5.4M is a
# DIFFERENT INSTRUMENT from the one the authors validated.
#
# REVISION 2026-07-19 — RECOVERY IS NOW REPORTED PER COHORT. The previous version pooled
# every cohort's .sscore files and took one median, reporting "dosage recovery 78-90%
# (median 87%)". That hid the finding audit_B02 turned up: recovery is strongly
# DIFFERENTIAL BY ANCESTRY. For the EUR SBP score it runs 80.7% in AMANHI-Pakistan down to
# 67.6% in GAPPS-Zambia -- a 13-point spread that attenuates African R2 independently of
# any ancestry-mismatch effect, and is therefore a competing explanation for the paper's
# African result. A single median makes that invisible.
#   S1  keeps one row per score, but carries min/median/max recovery across cohorts.
#   S1b is the full score x cohort x platform detail.
#
# Catalog metadata (published variant count, method, PMID, discovery/training N) is joined
# from ref/pgs_catalog_metadata.tsv -- externally sourced, with URL and retrieval date per
# row. See ref/README.md. Nothing here is filled from memory.
#
# Writes: tables/S1_prs_manifest.tsv, tables/S1b_recovery_by_cohort.tsv
#   Rscript deliv_S1_prs_manifest.R   (SSC from env)
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)
SSC <- momi_arg("--sscore-dir", Sys.getenv("SSC"))
CATALOG <- file.path(PIPE, "ref", "pgs_catalog_metadata.tsv")

momi_deliverable("S1_prs_manifest", script="05_prs/deliv_S1_prs_manifest.R",
                 inputs="MOMI_PANEL;sscore;ref/pgs_catalog_metadata.tsv", P=P, body=function(ctx){
  if(is.null(SSC) || !nzchar(SSC)) stop("no SSC: export SSC")
  if(!file.exists(CATALOG))
    stop(sprintf("catalog metadata missing: %s (externally sourced, must be committed)", CATALOG))

  ## ---- variants used, per score x cohort x platform ----
  ## DENOM is the allele-observation count; audit_B02 confirmed it is CONSTANT within every
  ## file (cv = 0 in all 80), so the median is exact rather than a summary.
  vc_rows <- list()
  for(pid in MOMI_PANEL$id) for(pl in MOMI_PLATS) for(coh in MOMI_COH_ALL){
    f <- file.path(SSC, sprintf("%s__%s__%s.sscore", pid, coh, pl))
    if(!file.exists(f)) next
    d <- tryCatch(fread(f), error=function(e) NULL)
    if(is.null(d) || !("DENOM" %in% names(d))) next
    vc_rows[[paste(pid,coh,pl)]] <- data.table(
      id=pid, cohort=coh, platform=pl, n_samples=nrow(d),
      variants_used=as.integer(round(median(as.numeric(d$DENOM), na.rm=TRUE)/2)))
  }
  VC <- rbindlist(vc_rows, fill=TRUE)

  ## ---- catalog metadata (externally sourced) ----
  CAT <- fread(CATALOG, colClasses=list(character="PMID"))
  miss <- setdiff(MOMI_PANEL$id, CAT$id)
  if(length(miss))
    stop(sprintf("catalog rows missing for: %s -- add them to ref/pgs_catalog_metadata.tsv",
                 paste(miss, collapse=", ")))

  ## ---- S1b: the full per-cohort detail ----
  S1b <- merge(VC, CAT[, .(id, n_snps_in_score)], by="id", all.x=TRUE)
  S1b <- merge(S1b, as.data.table(MOMI_PANEL)[, .(id, trait, anc, role)], by="id", all.x=TRUE)
  S1b[, cohort_anc := unname(MOMI_ANC[cohort])]
  S1b[, recovery_pct := round(100*variants_used/n_snps_in_score, 1)]
  setcolorder(S1b, c("id","trait","anc","role","cohort","cohort_anc","platform",
                     "n_samples","n_snps_in_score","variants_used","recovery_pct"))
  setorder(S1b, id, platform, -recovery_pct)
  out_b <- momi_write_table(S1b, "S1b_recovery_by_cohort", P)

  ## ---- S1: one row per score, recovery summarised ACROSS cohorts (min/median/max) ----
  agg <- function(pl) S1b[platform==pl, .(
      cohorts_scored   = .N,
      variants_med     = as.integer(median(variants_used, na.rm=TRUE)),
      recovery_min     = min(recovery_pct, na.rm=TRUE),
      recovery_median  = as.numeric(median(recovery_pct, na.rm=TRUE)),
      recovery_max     = max(recovery_pct, na.rm=TRUE)), by=id]
  Ag <- agg("gsa");   setnames(Ag, setdiff(names(Ag),"id"), paste0(setdiff(names(Ag),"id"),"_gsa"))
  Ad <- agg("lpwgs_dosage"); setnames(Ad, setdiff(names(Ad),"id"), paste0(setdiff(names(Ad),"id"),"_dosage"))

  M <- as.data.table(MOMI_PANEL)
  M <- merge(M, CAT, by="id", all.x=TRUE, sort=FALSE)
  M <- merge(M, Ag, by="id", all.x=TRUE, sort=FALSE)
  M <- merge(M, Ad, by="id", all.x=TRUE, sort=FALSE)
  M[, on_our_platforms := as.integer(!is.na(cohorts_scored_gsa) | !is.na(cohorts_scored_dosage))]

  M[, notes := fifelse(on_our_platforms==0,
                       "not scored on gsa/lpwgs_dosage (raw-lpwgs only or absent)", "")]
  ## PRSmix/PRSmixPlus scores run no new discovery GWAS -- flag so the N is not misread
  M[is.na(discovery_N) & !is.na(training_N),
    notes := paste0(notes, fifelse(nzchar(notes), "; ", ""),
                    "mixture score: no new discovery GWAS; N shown is training sample")]
  ## flag differential recovery, the reason this table was revised
  M[!is.na(recovery_max_dosage) & (recovery_max_dosage - recovery_min_dosage) >= 5,
    notes := paste0(notes, fifelse(nzchar(notes), "; ", ""),
                    sprintf("recovery varies %.1f-%.1f%% across cohorts (see S1b)",
                            recovery_min_dosage, recovery_max_dosage))]

  setcolorder(M, c("id","trait","anc","role","source","publication","PMID","method",
                   "n_snps_in_score","discovery_N","training_N","discovery_ancestry",
                   "on_our_platforms","cohorts_scored_gsa","cohorts_scored_dosage",
                   "variants_med_gsa","variants_med_dosage",
                   "recovery_min_gsa","recovery_median_gsa","recovery_max_gsa",
                   "recovery_min_dosage","recovery_median_dosage","recovery_max_dosage",
                   "catalog_url","retrieved","notes"))
  out <- momi_write_table(M, "S1_prs_manifest", P)

  ## ---- console: the differential this revision exists to expose ----
  cat("\nDosage-platform recovery by cohort (%), scored scores only:\n")
  w <- dcast(S1b[platform=="lpwgs_dosage"], id + trait + anc ~ cohort, value.var="recovery_pct")
  print(w)
  cat("\nmean recovery by COHORT ancestry:\n")
  print(S1b[platform=="lpwgs_dosage", .(scores=.N, mean_recovery=round(mean(recovery_pct),1)),
            by=cohort_anc][order(mean_recovery)])
  cat("\nA lower recovery in the AFR cohorts attenuates their R2 for reasons unrelated to\n",
      "ancestry mismatch, and is a competing explanation for the paper's African result.\n", sep="")

  rec <- S1b[platform=="lpwgs_dosage" & !is.na(recovery_pct)]
  list(n=nrow(M),
       key=sprintf("panel=%d scored=%d; dosage recovery %.0f-%.0f%% across cohorts (worst=%s); S1b rows=%d",
                   nrow(M), sum(M$on_our_platforms),
                   min(rec$recovery_pct), max(rec$recovery_pct),
                   rec[which.min(recovery_pct)]$cohort, nrow(S1b)),
       outputs=c(basename(out), basename(out_b)))
})
