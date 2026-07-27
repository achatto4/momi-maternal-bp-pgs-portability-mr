#!/usr/bin/env Rscript
# ============================================================
# deliv_T1_cohorts.R  [ B08 -> Table 1, cohort characteristics ]
# One row per genotyped cohort (consortium display names), + a THSTI phenotype-only note.
# Columns: country/ancestry; BP protocol & median GA at first BP reading (the AMANHI-late
# vs GAPPS-early contrast); N genotyped (GSA / lpWGS-dosage / union); N Part I / Part II;
# maternal age, gravidity, BMI (mean+-SD, genotyped Part I); PTB/LBW/SGA %; preeclampsia N(%).
# Reads analytic_mothers, bp_readings_long, prs_z_platform (all frozen intermediates).
#
# Writes: tables/T1_cohorts.tsv
#   Rscript deliv_T1_cohorts.R
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

momi_deliverable("T1_cohorts", script="01_phenotypes/deliv_T1_cohorts.R",
                 inputs="analytic_mothers;bp_readings_long;prs_z_platform", P=P, body=function(ctx){
  A   <- momi_read_intermediate("analytic_mothers", P)
  LNG <- momi_read_intermediate("bp_readings_long", P)
  PLZ <- momi_read_intermediate("prs_z_platform", P)
  A[, hasBP := is.finite(S_mean) | is.finite(D_mean)]

  msd <- function(x) if(all(is.na(x))) "--" else sprintf("%.1f±%.1f", mean(x,na.rm=TRUE), sd(x,na.rm=TRUE))
  pct <- function(x){ x<-x[!is.na(x)]; if(!length(x)) "--" else sprintf("%d (%.1f%%)", sum(x), 100*mean(x)) }

  # median GA (weeks) at each mother's FIRST valid-GA BP reading, per cohort
  firstga <- LNG[is.finite(GA_days), .(gfirst=min(GA_days)), by=.(cohort,IID)][
                 , .(ga_first_wk=round(median(gfirst)/7,1),
                     ga_p10_wk =round(quantile(gfirst,.10)/7,1)), by=cohort]
  # platform genotyping counts
  plc <- PLZ[, .(N=uniqueN(IID)), by=.(cohort,platform)]
  gsaN <- plc[platform=="gsa",         .(cohort, gsa=N)]
  dosN <- plc[platform=="lpwgs_dosage",.(cohort, dosage=N)]
  ## Dual-platform count, carried as a COLUMN rather than asserted in prose. The provenance
  ## note used to claim "GAPPS = same mothers on both platforms"; audit_B02 measured it and
  ## that is false (GAPPS-B 7.0%, Zambia 30.9%, Pakistan 0.0%). A table that carries the
  ## number cannot drift from the data the way a sentence can.
  bothN <- PLZ[, .(np=uniqueN(platform)), by=.(cohort,IID)][np==2, .(both=.N), by=cohort]

  protocol <- function(coh) if(grepl("^AMANHI", coh)) "AMANHI (late)" else "GAPPS (early)"

  ## descriptive block shared by the per-cohort rows and the Overall column
  desc <- function(pI){
    pII <- pI[!is.na(PTB)]; lb <- pII[livebirth==1]
    data.table(
      age = msd(pI$AGE), gravidity = msd(pI$GRAV), BMI = msd(pI$BMI),
      SBP = msd(pI$S_mean), DBP = msd(pI$D_mean),
      GAdel_wk = msd(lb$GAd/7), BWT_g = msd(lb$BWT),
      ## live-birth outcomes on the causal-analysis subset; STILLBIRTH is a pregnancy-level
      ## outcome and must be counted on the full sample (pI), not on the live-birth subset pII
      ## which excludes stillbirths by construction (their PTB status is missing).
      PTB_pct = pct(pII$PTB), LBW_pct = pct(pII$LBW), SGA_pct = pct(pII$SGA),
      still_pct = pct(pI$still), PE_n_pct = pct(pI$PE),
      cHTN_pct = if("CHRON_HTN" %in% names(pI)) pct(pI$CHRON_HTN) else "--")
  }

  rows <- lapply(MOMI_COH_ALL, function(coh){
    a  <- A[cohort==coh]
    pI <- a[genotyped==1 & hasBP==TRUE]                 # Part I analytic
    cbind(data.table(
      cohort_display = MOMI_DISPLAY[coh], ancestry = MOMI_ANC[coh],
      protocol = protocol(coh),
      GA_first_med_wk = firstga[cohort==coh, ga_first_wk],
      GA_first_p10_wk = firstga[cohort==coh, ga_p10_wk],
      N_gsa    = { v<-gsaN[cohort==coh, gsa]; if(length(v)) v else 0L },
      N_dosage = { v<-dosN[cohort==coh, dosage]; if(length(v)) v else 0L },
      N_both_platforms = { v<-bothN[cohort==coh, both]; if(length(v)) v else 0L },
      pct_both = { v<-bothN[cohort==coh, both]; u<-sum(a$genotyped==1, na.rm=TRUE)
                   if(length(v) && u>0) round(100*v/u,1) else 0 },
      N_genotyped_union = sum(a$genotyped==1, na.rm=TRUE),
      N_partI = nrow(pI), N_partII = nrow(pI[!is.na(PTB)])),
      desc(pI))
  })
  T1 <- rbindlist(rows)
  setorder(T1, ancestry, cohort_display)

  ## Overall column: all five genotyped cohorts pooled
  pIall <- A[cohort %in% MOMI_COH_ALL & genotyped==1 & (is.finite(S_mean)|is.finite(D_mean))]
  T1 <- rbind(T1, cbind(data.table(
    cohort_display="Overall", ancestry="all", protocol="",
    GA_first_med_wk=NA_real_, GA_first_p10_wk=NA_real_, N_gsa=0L, N_dosage=0L,
    N_both_platforms=0L, pct_both=0, N_genotyped_union=nrow(pIall),
    N_partI=nrow(pIall), N_partII=nrow(pIall[!is.na(PTB)])), desc(pIall)), fill=TRUE)

  # THSTI (phenotype-only) note row
  th <- A[cohort=="THSTI-India"]
  if(nrow(th)) T1 <- rbind(T1, cbind(data.table(
    cohort_display=MOMI_DISPLAY[["THSTI-India"]], ancestry="SAS (pheno-only)",
    protocol="GARBH-INi", GA_first_med_wk=NA_real_, GA_first_p10_wk=NA_real_,
    N_gsa=0L, N_dosage=0L, N_both_platforms=0L, pct_both=0, N_genotyped_union=0L,
    N_partI=0L, N_partII=0L), desc(th)), fill=TRUE)

  out <- momi_write_table(T1, "T1_cohorts", P)
  list(n=sum(T1$N_partI),
       key=sprintf("cohorts=%d genotyped-union=%d partI=%d partII=%d",
                   nrow(T1[N_genotyped_union>0]), sum(T1$N_genotyped_union),
                   sum(T1$N_partI), sum(T1$N_partII)),
       outputs=basename(out))
})
