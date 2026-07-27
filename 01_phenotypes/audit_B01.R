#!/usr/bin/env Rscript
# ============================================================
# audit_B01.R — independent, READ-ONLY audit of the B01 analytic build.
#
# Purpose: re-derive B01's inputs straight from the raw EPI file and check every
# assumption from scratch, WITHOUT trusting any previously locked number and WITHOUT
# modifying build_analytic.R. Anything B01 asserts, this recomputes independently.
#
# It answers, in order:
#   A  file shape, expected columns, PREGNANCY_ID filter
#   B  is the SITE_CODE -> cohort map really right? (crosstab vs participant-ID prefix)
#   C  does VISITDT parse to a number? (if not, "first reading" silently means "lowest GA")
#   D  gestational-age quality, out-of-range by site
#   E  BP readings per mother by cohort (drives ICC and the mean/median definitions)
#   F  WHEN are height and weight measured? (BMI comparability across cohorts)
#   G  missingness of every analysis variable, by cohort
#   H  missing-becomes-zero audit for PE / CHRON_HTN / twin
#   I  headline counts, recomputed
#   J  ICC by cohort x trait, recomputed
#   K  BP-definition coverage by cohort
#   L  genotyped flag: does score membership agree across PGS?
#
# Run:
#   Rscript 01_phenotypes/audit_B01.R          # EPI/SSC from environment
# Writes: results/current/qc/B01_audit_*.tsv   (console output is the main product)
# ============================================================
suppressMessages({library(data.table)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R"))
source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)
QC <- file.path(P$root, "qc"); dir.create(QC, showWarnings=FALSE, recursive=TRUE)

EPI <- momi_arg("--epi", Sys.getenv("EPI"))
SSC <- momi_arg("--sscore-dir", Sys.getenv("SSC"))
if(is.null(EPI) || !nzchar(EPI)) stop("no EPI: export EPI")

hdr <- function(x) cat(sprintf("\n\n========== %s ==========\n", x))
sec <- function(x) cat(sprintf("\n-- %s --\n", x))
sav <- function(dt, nm){ fwrite(dt, file.path(QC, paste0("B01_audit_",nm,".tsv")), sep="\t"); invisible(NULL) }
pr  <- function(dt, n=40) print(head(dt, n), nrows=n, class=FALSE)

cat("audit_B01.R —", format(Sys.time(), "%F %T"), "\ngit:", momi_git_sha(PIPE), "\nEPI:", EPI, "\n")

## ============================================================
hdr("A. FILE SHAPE")
full <- fread(EPI, na.strings=MOMI_NA)
cat(sprintf("raw dimensions: %d rows x %d cols\n", nrow(full), ncol(full)))
cat("column names:\n"); print(names(full))

need <- c("PARTICIPANT_ID","PREGNANCY_ID","SITE_CODE","SBP","DBP","GA_HDLK_NEW","VISITDT",
          "PW_AGE","GRAVIDITY","PARITY","PW_EDUCATION","WEALTH_INDEX","MAT_HEIGHT","MAT_WEIGHT",
          "CHRON_HTN","GAGEBRTH_NEW","BIRTH_WEIGHT","SGA_10_NEW","SGA_3_NEW","PTB_NEW",
          "SPONT_LABOUR","SINGLE_TWIN","BIRTH_OUTCOME","PE_CAT")
sec("expected columns present?")
print(data.table(column=need, present=need %in% names(full)))

sec("PREGNANCY_ID distribution (B01 keeps ==1)")
print(full[, .N, by=PREGNANCY_ID][order(PREGNANCY_ID)])

raw <- full[PREGNANCY_ID==1]
cat(sprintf("\nafter PREGNANCY_ID==1 filter: %d visit-rows, %d unique participants\n",
            nrow(raw), uniqueN(raw$PARTICIPANT_ID)))

## ============================================================
hdr("B. SITE_CODE -> COHORT MAP, VERIFIED BY PARTICIPANT-ID PREFIX")
raw[, IID := as.character(PARTICIPANT_ID)]
raw[, id_prefix := sub("[0-9].*$", "", IID)]
raw[, cohort := MOMI_SITE[as.character(SITE_CODE)]]
xt <- raw[, .(visit_rows=.N, mothers=uniqueN(IID)), by=.(SITE_CODE, cohort, id_prefix)][order(SITE_CODE)]
cat("If the map is correct each SITE_CODE has exactly ONE id_prefix:\n")
pr(xt)
bad <- xt[, .N, by=SITE_CODE][N>1]
cat(sprintf("\nSITE_CODEs with more than one prefix: %d %s\n", nrow(bad),
            if(nrow(bad)) paste0("<-- INVESTIGATE: ", paste(bad$SITE_CODE, collapse=",")) else "(clean)"))
sav(xt, "site_prefix")

## ============================================================
hdr("C. DOES VISITDT PARSE TO A NUMBER?")
cat("This decides what 'first reading' means. B01 does setorder(raw, IID, VIS, GA)\n",
    "with VIS <- as.numeric(VISITDT). If VISITDT is a character date, VIS is all NA\n",
    "and the ordering silently collapses to GA only.\n\n", sep="")
vr <- full$VISITDT
cat("class(VISITDT):", paste(class(vr), collapse="/"), "\n")
cat("first 8 raw values:", paste(head(as.character(vr),8), collapse=" | "), "\n")
## CORRECTION 2026-07-19: the earlier version of this test did as.numeric(as.character(vr)),
## which of course fails on a date string. B01 does as.numeric(vr) DIRECTLY, and an
## IDate/Date is numeric underneath (days since 1970-01-01), so it parses fine. Test what
## B01 actually does, not a strawman.
vnum_direct <- suppressWarnings(as.numeric(vr))                    # what B01 does
vnum_viachr <- suppressWarnings(as.numeric(as.character(vr)))      # the broken test
cat(sprintf("as.numeric(VISITDT)          -> non-NA: %d / %d (%.1f%%)   [this is what B01 does]\n",
            sum(!is.na(vnum_direct)), length(vnum_direct), 100*mean(!is.na(vnum_direct))))
cat(sprintf("as.numeric(as.character(...)) -> non-NA: %d / %d (%.1f%%)   [strawman, ignore]\n",
            sum(!is.na(vnum_viachr)), length(vnum_viachr), 100*mean(!is.na(vnum_viachr))))
if(sum(!is.na(vnum_direct))) cat(sprintf("numeric range: %.0f to %.0f (days since 1970-01-01) = %s to %s\n",
        min(vnum_direct,na.rm=TRUE), max(vnum_direct,na.rm=TRUE),
        as.Date(min(vnum_direct,na.rm=TRUE), origin="1970-01-01"),
        as.Date(max(vnum_direct,na.rm=TRUE), origin="1970-01-01")))
cat(if(mean(is.na(vnum_direct)) > 0.5)
      ">>> VERDICT: VISITDT does NOT parse. 'first' is effectively 'lowest valid GA'.\n"
    else ">>> VERDICT: VISITDT parses. Visit-date ordering in B01 is REAL and correct.\n")

## ============================================================
hdr("D. GESTATIONAL AGE — FULL TAIL BEHAVIOUR")
cat("Goal: decide the bounds from evidence, not convention. B01 currently nulls GA outside\n",
    "[28,315] d. Decision taken 2026-07-19: DROP the lower bound. The upper bound is still\n",
    "open -- these tables are what decides it.\n", sep="")
raw[, GAraw := suppressWarnings(as.numeric(GA_HDLK_NEW))]
raw[, GAok  := ifelse(GAraw >= 28 & GAraw <= 315, GAraw, NA_real_)]   # CURRENT rule, for reference

sec("D1. overall distribution (days)")
cat(sprintf("rows=%d  non-missing GA=%d (%.1f%%)\n",
            nrow(raw), sum(!is.na(raw$GAraw)), 100*mean(!is.na(raw$GAraw))))
qs <- quantile(raw$GAraw, probs=c(0,.001,.005,.01,.05,.25,.5,.75,.95,.99,.995,.999,1), na.rm=TRUE)
print(round(qs,1))
cat("\nsame in weeks:\n"); print(round(qs/7,1))

sec("D2. LOW tail — every distinct value below 28 days, by cohort")
cat("If these are real early visits we keep them; if they are 0/negative/codes we do not.\n\n")
low <- raw[!is.na(GAraw) & GAraw < 28, .N, by=.(cohort, GAraw)][order(GAraw, cohort)]
if(nrow(low)) pr(low, 60) else cat("(none below 28 d)\n")
cat(sprintf("\ntotal readings below 28 d: %d\n", sum(low$N)))
sav(low, "ga_low_tail")

sec("D3. HIGH tail — distinct values above 294 days (42 wk), by cohort")
high <- raw[!is.na(GAraw) & GAraw > 294, .N, by=.(cohort, GAraw)][order(-GAraw, cohort)]
if(nrow(high)) pr(high, 80) else cat("(none above 294 d)\n")
cat(sprintf("\ntotal readings above 294 d: %d ; above 315 d: %d\n",
            sum(high$N), sum(raw$GAraw > 315, na.rm=TRUE)))
sav(high, "ga_high_tail")

sec("D4. HIGH tail binned — is it a smooth tail or a cliff?")
hb <- raw[!is.na(GAraw) & GAraw > 280]
if(nrow(hb)){
  hb[, bin := cut(GAraw, breaks=c(280,294,300,308,315,330,350,400,500,1000,Inf),
                  right=TRUE, dig.lab=6)]
  hbt <- dcast(hb[, .N, by=.(cohort,bin)], bin ~ cohort, value.var="N", fill=0)
  pr(hbt, 30); sav(hbt, "ga_high_binned")
} else cat("(nothing above 280 d)\n")

sec("D5. most frequent GA values overall — spots sentinel/missing codes")
topv <- raw[!is.na(GAraw), .N, by=GAraw][order(-N)][1:25]
topv[, weeks := round(GAraw/7,1)]
pr(topv, 25); sav(topv, "ga_top_values")

sec("D6. are out-of-range-GA rows junk rows, or good BP with a bad date?")
cat("This is the decisive question for the upper bound. If a GA=400 row carries a perfectly\n",
    "normal SBP, only the DATE is wrong -- the BP reading is real and should be kept (it\n",
    "already is; only window/trimester definitions drop it). If those rows ALSO carry absurd\n",
    "BP, the whole row is junk and should go entirely.\n\n", sep="")
raw[, SBPn := suppressWarnings(as.numeric(SBP))]
raw[, DBPn := suppressWarnings(as.numeric(DBP))]
raw[, ga_class := fifelse(is.na(GAraw), "1_GA missing",
                   fifelse(GAraw < 28,  "2_GA <28d",
                     fifelse(GAraw > 315, "4_GA >315d", "3_GA in [28,315]")))]
bpc <- raw[, .(n_rows=.N,
               pct_SBP_present   = round(100*mean(is.finite(SBPn)),1),
               median_SBP        = round(median(SBPn, na.rm=TRUE),1),
               pct_SBP_plausible = round(100*mean(is.finite(SBPn) & SBPn>=60 & SBPn<=250),1),
               median_DBP        = round(median(DBPn, na.rm=TRUE),1),
               pct_DBP_plausible = round(100*mean(is.finite(DBPn) & DBPn>=30 & DBPn<=150),1)),
           by=ga_class][order(ga_class)]
pr(bpc); sav(bpc, "ga_class_vs_bp")
cat("\nRead across: if pct_SBP_plausible is similar in the >315d row and the in-range row,\n",
    "the BP data are fine and only the dates are broken.\n", sep="")

sec("D7. per-cohort summary under the CURRENT [28,315] rule")
gq <- raw[, .(n_rows=.N,
              n_GA_missing   = sum(is.na(GAraw)),
              n_below_28     = sum(!is.na(GAraw) & GAraw < 28),
              n_above_315    = sum(!is.na(GAraw) & GAraw > 315),
              pct_nulled     = round(100*mean(!is.na(GAraw) & (GAraw < 28 | GAraw > 315)),1),
              GA_min = suppressWarnings(min(GAraw, na.rm=TRUE)),
              GA_max = suppressWarnings(max(GAraw, na.rm=TRUE))), by=cohort][order(cohort)]
pr(gq); sav(gq, "ga_quality")
cat("\nn_below_28 vs n_above_315 splits the previously quoted '13% out of range' into the\n",
    "two tails, and shows which cohort contributes the 7,173.\n", sep="")

sec("D8. can GA be RECONSTRUCTED from the delivery date? (repair, not discard)")
cat("The file carries VISITDT, DEL_DATE, DATE_LMP and GAGEBRTH_NEW. If the dates parse, then\n",
    "   GA_at_visit  =  GA_at_birth  -  (delivery date - visit date)\n",
    "gives an INDEPENDENT estimate of gestational age at each visit. Two payoffs: it\n",
    "validates GA_HDLK_NEW where that is in range, and it may RECOVER the ~13% of readings\n",
    "we currently throw away -- which would beat any choice of upper bound.\n\n", sep="")

parse_date <- function(x, label){
  xc <- as.character(x)
  fmts <- c("%Y-%m-%d","%d/%m/%Y","%m/%d/%Y","%Y/%m/%d","%d-%m-%Y","%d%b%Y","%d-%b-%Y","%b %d %Y")
  best <- list(d=NULL, fmt="UNPARSED", rate=0)
  for(f in fmts){
    d <- suppressWarnings(as.Date(xc, format=f))
    r <- mean(!is.na(d))
    if(r > best$rate) best <- list(d=d, fmt=f, rate=r)
  }
  n <- suppressWarnings(as.numeric(xc))
  cat(sprintf("  %-18s class=%-10s best-format=%-12s parse-rate=%.1f%%  numeric-rate=%.1f%%\n",
              label, paste(class(x), collapse="/"), best$fmt, 100*best$rate, 100*mean(!is.na(n))))
  cat(sprintf("  %-18s examples: %s\n", "", paste(utils::head(xc[!is.na(xc)],4), collapse=" | ")))
  if(best$rate > 0.5) best$d else NULL
}
vd <- parse_date(raw$VISITDT,  "VISITDT")
dd <- parse_date(raw$DEL_DATE, "DEL_DATE")
ld <- parse_date(raw$DATE_LMP, "DATE_LMP")

raw[, GAbirth := suppressWarnings(as.numeric(GAGEBRTH_NEW))]
if(!is.null(vd) && !is.null(dd)){
  raw[, days_to_del := as.numeric(dd - vd)]
  raw[, GA_recon    := GAbirth - days_to_del]
  cat(sprintf("\n  reconstructable rows: %d of %d (%.1f%%)\n",
              sum(is.finite(raw$GA_recon)), nrow(raw), 100*mean(is.finite(raw$GA_recon))))

  cmpd <- raw[is.finite(GA_recon) & is.finite(GAraw),
              .(n=.N,
                median_diff = round(median(GAraw - GA_recon),1),
                pct_within_7d  = round(100*mean(abs(GAraw - GA_recon) <= 7),1),
                pct_within_14d = round(100*mean(abs(GAraw - GA_recon) <= 14),1)),
              by=ga_class][order(ga_class)]
  cat("\n  agreement between GA_HDLK_NEW and the reconstruction, by GA class:\n")
  pr(cmpd); sav(cmpd, "ga_reconstruction_agreement")
  cat("\n  Read the '3_GA in [28,315]' row first: high pct_within_7d there means the\n",
      "  reconstruction is trustworthy. Then read the '4_GA >315d' row -- if the\n",
      "  reconstruction gives sane values where GA_HDLK_NEW does not, those readings are\n",
      "  recoverable and we should repair rather than null them.\n", sep="")

  rec <- raw[ga_class %in% c("2_GA <28d","4_GA >315d","1_GA missing") & is.finite(GA_recon),
             .(n=.N,
               recon_min = round(min(GA_recon),0), recon_med = round(median(GA_recon),0),
               recon_max = round(max(GA_recon),0),
               pct_recon_plausible = round(100*mean(GA_recon >= 0 & GA_recon <= 315),1)),
             by=.(cohort, ga_class)][order(ga_class, cohort)]
  cat("\n  what the reconstruction says for currently-discarded readings:\n")
  pr(rec, 40); sav(rec, "ga_reconstruction_rescue")
} else {
  cat("\n  >>> cannot reconstruct: VISITDT and/or DEL_DATE did not parse as dates.\n")
  cat("  If they are numeric date serials, tell me the origin (Stata=1960-01-01,\n")
  cat("  Excel=1899-12-30) and I will redo this section.\n")
}

sec("D9. INDEPENDENT GA check from DATE_LMP (D8 proved circular)")
cat("The 2026-07-19 run showed D8 agreeing with GA_HDLK_NEW to the DAY in every class,\n",
    "including the out-of-range one (median_diff=0, 100% within 7 d). That means\n",
    "GAGEBRTH_NEW and GA_HDLK_NEW are two views of the SAME arithmetic -- the reconstruction\n",
    "is not independent and validates nothing. DATE_LMP is a separate source:\n",
    "   GA_at_visit = VISITDT - DATE_LMP\n",
    "If THAT is sane where GA_HDLK_NEW is not, the broken readings are repairable.\n\n", sep="")
if(!is.null(vd) && !is.null(ld)){
  ## Attach the parsed dates as COLUMNS. Using the bare vectors inside a `by=` j-expression
  ## does NOT subset them per group -- dd[1] would return row 1 of the whole table for every
  ## group. That bug made D9b report an identical 250 d gestation in all six cohorts.
  raw[, VD := vd]
  if(!is.null(dd)) raw[, DD := dd]
  raw[, LD := ld]
  raw[, GA_lmp := as.numeric(vd - ld)]
  cat(sprintf("  rows with an LMP-based GA: %d of %d (%.1f%%)\n",
              sum(is.finite(raw$GA_lmp)), nrow(raw), 100*mean(is.finite(raw$GA_lmp))))

  cmpl <- raw[is.finite(GA_lmp) & is.finite(GAraw),
              .(n=.N,
                median_diff    = round(median(GAraw - GA_lmp),1),
                pct_within_7d  = round(100*mean(abs(GAraw - GA_lmp) <= 7),1),
                pct_within_14d = round(100*mean(abs(GAraw - GA_lmp) <= 14),1)),
              by=ga_class][order(ga_class)]
  cat("\n  GA_HDLK_NEW vs LMP-based GA, by class (a NON-zero median_diff here is GOOD --\n",
      "  it means the two sources are genuinely independent):\n", sep="")
  pr(cmpl); sav(cmpl, "ga_lmp_agreement")

  resc <- raw[ga_class %in% c("1_GA missing","2_GA <28d","4_GA >315d") & is.finite(GA_lmp),
              .(n=.N, lmp_min=round(min(GA_lmp)), lmp_med=round(median(GA_lmp)),
                lmp_max=round(max(GA_lmp)),
                pct_lmp_plausible = round(100*mean(GA_lmp >= 0 & GA_lmp <= 315),1)),
              by=.(cohort, ga_class)][order(ga_class, cohort)]
  cat("\n  what the LMP-based GA says for currently-discarded readings:\n")
  cat("  (a high pct_lmp_plausible = those readings are RECOVERABLE)\n")
  pr(resc, 40); sav(resc, "ga_lmp_rescue")

  sec("D9b. gestation length from dates vs recorded GAGEBRTH_NEW")
  cat("  DEL_DATE - DATE_LMP should approximate GAGEBRTH_NEW. Disagreement localises the\n",
      "  broken field: if gestation-from-dates is sane but GAGEBRTH_NEW is not, the outcome\n",
      "  variable is at fault; if the dates are wild, the dates are.\n\n", sep="")
  if(!is.null(dd)){
    mo <- raw[, .(gest_dates = as.numeric(DD[1] - LD[1]), gest_rec = GAbirth[1]), by=.(cohort,IID)]
    ms <- mo[is.finite(gest_dates) & is.finite(gest_rec),
             .(mothers=.N,
               med_gest_from_dates = round(median(gest_dates)),
               med_gest_recorded   = round(median(gest_rec)),
               pct_dates_plausible = round(100*mean(gest_dates>=140 & gest_dates<=320),1),
               pct_rec_plausible   = round(100*mean(gest_rec  >=140 & gest_rec  <=320),1),
               pct_agree_within_7d = round(100*mean(abs(gest_dates-gest_rec)<=7),1)),
             by=cohort][order(cohort)]
    pr(ms); sav(ms, "gestation_dates_vs_recorded")
  }
} else {
  cat("  >>> DATE_LMP or VISITDT unavailable as dates; cannot run the independent check.\n")
}

## ============================================================
hdr("E. BP READINGS PER MOTHER, BY COHORT")
raw[, SBPn := suppressWarnings(as.numeric(SBP))]
raw[, DBPn := suppressWarnings(as.numeric(DBP))]
## NB: median() on an integer vector returns integer for odd n and double for even n, so
## data.table sees inconsistent column types across groups and aborts. Force numeric.
rp <- raw[is.finite(SBPn), .N, by=.(cohort, IID)][
        , .(mothers=.N,
            med_readings = as.numeric(median(N)),
            q1           = as.numeric(quantile(N,.25)),
            q3           = as.numeric(quantile(N,.75)),
            max_readings = as.numeric(max(N)),
            pct_with_ge2 = round(100*mean(N>=2),1)), by=cohort][order(cohort)]
cat("SBP readings per mother (drives ICC, and why 'mean' beats a single reading):\n")
pr(rp); sav(rp, "readings_per_mother")

sec("GA at FIRST valid BP reading, by cohort (the timing contrast)")
setorder(raw, IID, GAok, na.last=TRUE)
ft <- raw[is.finite(SBPn) & is.finite(GAok), .(ga_first=GAok[1]), by=.(cohort,IID)][
        , .(mothers=.N, median_wk=round(median(ga_first)/7,1),
            q1_wk=round(quantile(ga_first,.25)/7,1), q3_wk=round(quantile(ga_first,.75)/7,1),
            pct_before_20wk=round(100*mean(ga_first < 140),1)), by=cohort][order(cohort)]
pr(ft); sav(ft, "ga_at_first_bp")

## ============================================================
hdr("F. WHEN ARE HEIGHT AND WEIGHT MEASURED? (BMI COMPARABILITY)")
cat("B01 takes the first non-missing HT and WT INDEPENDENTLY. Height is constant so that\n",
    "is harmless, but weight rises through pregnancy. If cohorts differ in when weight is\n",
    "first recorded, 'BMI' is not the same quantity across cohorts -- and BMI is a core\n",
    "confounder in every adjusted model, and the subject of B18.\n\n", sep="")
raw[, HTn := suppressWarnings(as.numeric(MAT_HEIGHT))]
raw[, WTn := suppressWarnings(as.numeric(MAT_WEIGHT))]
wt_ga <- raw[is.finite(WTn) & is.finite(GAok), .(ga_wt=GAok[1]), by=.(cohort,IID)]
ht_ga <- raw[is.finite(HTn) & is.finite(GAok), .(ga_ht=GAok[1]), by=.(cohort,IID)]
wsum <- wt_ga[, .(mothers=.N, median_wk_at_weight=round(median(ga_wt)/7,1),
                  q1=round(quantile(ga_wt,.25)/7,1), q3=round(quantile(ga_wt,.75)/7,1)),
              by=cohort][order(cohort)]
cat("GA (weeks) at the FIRST recorded maternal WEIGHT:\n"); pr(wsum); sav(wsum,"ga_at_weight")

both <- merge(wt_ga, ht_ga, by=c("cohort","IID"))
both[, same_visit := abs(ga_wt - ga_ht) < 1e-9]
bsum <- both[, .(n=.N, pct_ht_wt_same_visit=round(100*mean(same_visit),1),
                 median_gap_days=round(median(abs(ga_wt-ga_ht)),1)), by=cohort][order(cohort)]
cat("\nAre height and weight taken at the same visit?\n"); pr(bsum); sav(bsum,"ht_wt_same_visit")

bmi <- merge(raw[is.finite(WTn), .(WT=WTn[1]), by=.(cohort,IID)],
             raw[is.finite(HTn), .(HT=HTn[1]), by=.(cohort,IID)], by=c("cohort","IID"))
bmi[, BMI := { b <- WT/(HT/100)^2; ifelse(b>=12 & b<=60, b, NA_real_) }]
bsum2 <- bmi[, .(n=.N, BMI_ok=sum(!is.na(BMI)),
                 median_BMI=round(median(BMI,na.rm=TRUE),1),
                 mean_BMI=round(mean(BMI,na.rm=TRUE),1)), by=cohort][order(cohort)]
cat("\nComputed BMI by cohort:\n"); pr(bsum2); sav(bsum2,"bmi_by_cohort")

## ============================================================
hdr("F2. IS THE BP-TIMING CONTRAST REAL, OR AN ARTEFACT OF GA DATA QUALITY?")
cat("THE most important question in this audit. Part I's backbone is that AMANHI measures BP\n",
    "at ~24 wk and GAPPS at ~12-16 wk. But section F just showed AMANHI records maternal\n",
    "WEIGHT at 12.4-13.9 wk -- so those women ARE attending early. Two explanations:\n",
    "  (a) REAL: they attend early, but BP is not taken/recorded until later.\n",
    "  (b) ARTEFACT: early BP readings exist, but their GA is one of the 14-26%% that is\n",
    "      broken and gets nulled, pushing 'GA at first BP' spuriously late.\n",
    "(b) would mean the timing contrast -- and therefore much of F3 and the S8/S16 story --\n",
    "is partly manufactured by data quality. These two tests separate them.\n", sep="")

setorder(raw, IID, VISITDT, na.last=TRUE)
have_lmp <- "GA_lmp" %in% names(raw)

sec("F2a. at the visit where WEIGHT was first recorded, was BP also taken?")
cat("If AMANHI shows high weight-coverage but low BP-coverage at that same early visit,\n",
    "explanation (a) holds: the visit happened, BP simply was not recorded.\n\n", sep="")
wf <- raw[is.finite(WTn), .SD[1], by=.(cohort, IID), .SDcols=c("SBPn","GAraw","GAok")]
wfs <- wf[, .(mothers=.N,
              pct_SBP_present_at_that_visit = round(100*mean(is.finite(SBPn)),1),
              pct_GA_valid_at_that_visit    = round(100*mean(is.finite(GAok)),1),
              median_GA_wk = round(median(GAok, na.rm=TRUE)/7,1)), by=cohort][order(cohort)]
pr(wfs); sav(wfs, "bp_at_first_weight_visit")

sec("F2b. GA at first BP — current rule vs an INDEPENDENT (LMP-based) GA")
cat("Column 1 is what the paper currently reports. Column 2 recomputes it using LMP dates,\n",
    "which are unaffected by the GA_HDLK_NEW corruption. If the AMANHI/GAPPS gap SHRINKS in\n",
    "column 2, part of the contrast was data quality, not protocol.\n\n", sep="")
if(have_lmp){
  fb <- raw[is.finite(SBPn), .(
            ga_first_current = { v <- GAok[is.finite(GAok)];  if(!length(v)) NA_real_ else v[1] },
            ga_first_lmp     = { v <- GA_lmp[is.finite(GA_lmp) & GA_lmp>=0 & GA_lmp<=315]
                                 if(!length(v)) NA_real_ else v[1] }),
            by=.(cohort, IID)]
  fbs <- fb[, .(mothers=.N,
                median_wk_current = round(median(ga_first_current, na.rm=TRUE)/7,1),
                median_wk_lmp     = round(median(ga_first_lmp,     na.rm=TRUE)/7,1),
                shift_wk = round((median(ga_first_current,na.rm=TRUE) -
                                  median(ga_first_lmp,   na.rm=TRUE))/7,1),
                pct_before20_current = round(100*mean(ga_first_current < 140, na.rm=TRUE),1),
                pct_before20_lmp     = round(100*mean(ga_first_lmp     < 140, na.rm=TRUE),1)),
            by=cohort][order(cohort)]
  pr(fbs); sav(fbs, "bp_timing_current_vs_lmp")
  cat("\nshift_wk > 0 means the current rule reports BP as being measured LATER than an\n",
      "independent GA says it was. A large positive shift in AMANHI and ~0 in GAPPS is the\n",
      "signature of explanation (b).\n", sep="")
} else cat("  (GA_lmp unavailable — D9 did not run)\n")

sec("F2c. how many BP readings are lost purely to a broken GA?")
lost <- raw[is.finite(SBPn), .(
          n_bp_readings   = .N,
          n_GA_valid      = sum(is.finite(GAok)),
          n_GA_broken     = sum(!is.finite(GAok)),
          pct_BP_lost_to_GA = round(100*mean(!is.finite(GAok)),1)), by=cohort][order(cohort)]
pr(lost); sav(lost, "bp_lost_to_broken_ga")
cat("\nThese readings still count toward mean/median/first/last (GA-free definitions) but are\n",
    "invisible to trimester, <20/>=20 wk and the GA-residual. That is the exact set of\n",
    "definitions S7 and S8 are built on.\n", sep="")

## ============================================================
hdr("F3. ARE THE >315 d READINGS POSTNATAL VISITS RATHER THAN ERRORS?")
cat("New hypothesis from the 2026-07-19 run. The >315 d readings are NOT random noise:\n",
    "  * they bunch just above the cutoff -- (315,330] and (330,350] hold most of them\n",
    "  * the LMP-based GA is ALSO implausible there, so BOTH clocks say 'late'\n",
    "  * they occur in AMANHI (which does postnatal follow-up) and are absent in GAPPS\n",
    "  * their BP values are perfectly normal (median SBP 111 vs 108)\n",
    "A visit 1-6 weeks AFTER delivery has days-since-LMP of roughly 287-322 -- exactly this\n",
    "pattern. If so they are not corrupt at all; they are postpartum visits, and the real\n",
    "problem is the opposite of what we assumed: nulling the GA removes them from the\n",
    "trimester/window definitions, but their BP STILL ENTERS mean/median/first/last, which\n",
    "are GA-free. That would contaminate the PRIMARY exposure, differentially by cohort.\n",
    "Decisive test: is the visit date after the delivery date?\n", sep="")

if(exists("dd") && !is.null(dd) && "DD" %in% names(raw)){
  raw[, days_after_del := as.numeric(VD - DD)]
  raw[, visit_class := fifelse(is.na(days_after_del), "unknown",
                        fifelse(days_after_del > 0, "POSTNATAL", "antenatal"))]

  sec("F3a. is a >315 d GA the same thing as a post-delivery visit?")
  xtab <- dcast(raw[, .N, by=.(ga_class, visit_class)], ga_class ~ visit_class,
                value.var="N", fill=0)
  pr(xtab); sav(xtab, "ga_class_vs_postnatal")
  cat("\nIf the '4_GA >315d' row is overwhelmingly POSTNATAL, the hypothesis holds.\n")

  sec("F3b. postnatal visits and postnatal BP readings, by cohort")
  pn <- raw[, .(visits=.N,
                postnatal_visits = sum(visit_class=="POSTNATAL"),
                pct_visits_postnatal = round(100*mean(visit_class=="POSTNATAL"),1),
                bp_readings = sum(is.finite(SBPn)),
                postnatal_bp = sum(is.finite(SBPn) & visit_class=="POSTNATAL"),
                pct_BP_postnatal = round(100*sum(is.finite(SBPn) & visit_class=="POSTNATAL")/
                                          max(sum(is.finite(SBPn)),1),1)),
            by=cohort][order(cohort)]
  pr(pn); sav(pn, "postnatal_by_cohort")
  cat("\npct_BP_postnatal is the share of the PRIMARY EXPOSURE that is postpartum rather\n",
      "than antenatal. Non-zero in AMANHI and zero in GAPPS = differential contamination.\n", sep="")

  sec("F3c. how different is postnatal BP from antenatal BP?")
  bpd <- raw[is.finite(SBPn) & visit_class %in% c("antenatal","POSTNATAL"),
             .(n=.N, mean_SBP=round(mean(SBPn),1), mean_DBP=round(mean(DBPn, na.rm=TRUE),1)),
             by=.(cohort, visit_class)][order(cohort, visit_class)]
  pr(bpd, 20); sav(bpd, "bp_antenatal_vs_postnatal")

  sec("F3d. DOES IT MATTER? mean BP with vs without postnatal readings")
  cat("The number that decides whether B01 must change: per-mother mean SBP computed the\n",
      "current way (all readings) versus antenatal-only.\n\n", sep="")
  mm <- merge(
    raw[is.finite(SBPn), .(mean_all = mean(SBPn), n_all=.N), by=.(cohort,IID)],
    raw[is.finite(SBPn) & visit_class=="antenatal",
        .(mean_ante = mean(SBPn), n_ante=.N), by=.(cohort,IID)],
    by=c("cohort","IID"), all.x=TRUE)
  mms <- mm[, .(mothers=.N,
                mothers_affected = sum(!is.na(n_ante) & n_all > n_ante),
                pct_affected = round(100*mean(!is.na(n_ante) & n_all > n_ante),1),
                mean_shift_mmHg = round(mean(mean_all - mean_ante, na.rm=TRUE),2),
                max_shift_mmHg  = round(max(abs(mean_all - mean_ante), na.rm=TRUE),1),
                lost_all_readings = sum(is.na(n_ante))),
            by=cohort][order(cohort)]
  pr(mms); sav(mms, "mean_bp_with_without_postnatal")
  cat("\nmean_shift_mmHg is the average bias in the primary exposure. A shift that is\n",
      "materially non-zero in AMANHI and exactly zero in GAPPS would mean the cohort\n",
      "comparison in Table 2 / F2 / F3 is partly comparing contaminated to clean.\n",
      "lost_all_readings = mothers whose ONLY BP readings are postnatal.\n", sep="")
} else {
  cat("\n  DEL_DATE unavailable — cannot classify visits.\n")
}

## ============================================================
hdr("G. MISSINGNESS BY COHORT (per mother, first non-missing value)")
## NB: return a TYPED NA (v[NA_integer_]), not bare NA. Bare NA is logical, so a group where
## every value is missing yields a logical column while other groups yield integer/character,
## and data.table aborts on the type mismatch. build_analytic.R already does this correctly.
fnn <- function(v){ w <- which(!is.na(v)); if(!length(w)) v[NA_integer_] else v[w[1]] }
## NB: SMOK_FREQ was dropped as ~99% missing, but the file carries TWO other tobacco
## variables that were never checked (SNIFF_TOBA, PASSIVE_SMOK). BABY_SEX is needed for
## S12's negative control, and BABY_ID is the mother-infant link S15 needs. Audit them all.
vars <- c("PW_AGE","GRAVIDITY","PARITY","PW_EDUCATION","HUS_EDUC","PW_OCCUPATION",
          "WEALTH_INDEX","MAT_HEIGHT","MAT_WEIGHT","CHRON_HTN","DIABETES",
          "SMOK_FREQ","SNIFF_TOBA","PASSIVE_SMOK",
          "GAGEBRTH_NEW","BIRTH_WEIGHT","SGA_10_NEW","PTB_NEW",
          "SPONT_LABOUR","SINGLE_TWIN","BIRTH_OUTCOME","BABY_SEX","BABY_ID","PE_CAT")
vars <- vars[vars %in% names(raw)]
permo <- raw[, lapply(.SD, fnn), by=.(cohort,IID), .SDcols=vars]
miss <- permo[, c(list(mothers=.N),
                  lapply(.SD, function(v) round(100*mean(is.na(v)),1))), by=cohort, .SDcols=vars]
cat("percent MISSING per mother, by cohort:\n"); pr(miss); sav(miss,"missingness")

## ============================================================
hdr("H. MISSING-BECOMES-ZERO AUDIT")
cat("B01 builds PE, CHRON_HTN and twin with patterns that turn 'no data' into 0, not NA.\n",
    "So a 0 can mean 'confirmed absent' OR 'never measured'. This quantifies the mixture.\n\n", sep="")
pe <- raw[, .(any_pecat = any(!is.na(PE_CAT)),
              pe_pos    = any(as.character(PE_CAT) %in% c("EOPE","LOPE")),
              pe_tbd    = any(as.character(PE_CAT) %in% c("TBD"))), by=.(cohort,IID)]
pesum <- pe[, .(mothers=.N, PE_positive=sum(pe_pos),
                PE_TBD_only=sum(pe_tbd & !pe_pos),
                PE_CAT_entirely_missing=sum(!any_pecat),
                pct_no_PE_data=round(100*mean(!any_pecat),1)), by=cohort][order(cohort)]
cat("PE_CAT ascertainment (matters for the S12 positive control):\n"); pr(pesum); sav(pesum,"pe_ascertainment")

ch <- raw[, .(any_chtn=any(!is.na(CHRON_HTN)), chtn_pos=any(CHRON_HTN==1, na.rm=TRUE)),
          by=.(cohort,IID)]
chsum <- ch[, .(mothers=.N, chronHTN_positive=sum(chtn_pos),
                CHRON_HTN_entirely_missing=sum(!any_chtn),
                pct_no_chtn_data=round(100*mean(!any_chtn),1)), by=cohort][order(cohort)]
cat("\nCHRON_HTN ascertainment:\n"); pr(chsum); sav(chsum,"chtn_ascertainment")

tw <- raw[, .(any_twin=any(!is.na(SINGLE_TWIN)),
              twin_pos=any(!is.na(SINGLE_TWIN) & SINGLE_TWIN!=1)), by=.(cohort,IID)]
twsum <- tw[, .(mothers=.N, twin_positive=sum(twin_pos),
                SINGLE_TWIN_missing=sum(!any_twin),
                pct_no_twin_data=round(100*mean(!any_twin),1)), by=cohort][order(cohort)]
cat("\nSINGLE_TWIN ascertainment:\n"); pr(twsum); sav(twsum,"twin_ascertainment")

## ============================================================
hdr("I. HEADLINE COUNTS, RECOMPUTED FROM SCRATCH")
o <- raw[, .(PTB=fnn(PTB_NEW), GAd=fnn(GAGEBRTH_NEW), BWTraw=fnn(BIRTH_WEIGHT),
             BOUT=fnn(BIRTH_OUTCOME)), by=.(cohort,IID)]
o[, GAd := ifelse(GAd>=140 & GAd<=320, GAd, NA_real_)]
o[, BWT := ifelse(BWTraw>=500 & BWTraw<=6500, BWTraw, NA_real_)]
o[, LBW := as.integer(BWT < 2500)]
tot <- data.table(
  quantity = c("first-pregnancy mothers","BMI computable","PTB (PTB_NEW==1)",
               "PE (EOPE|LOPE)","chronic HTN","twins","live births","stillbirths","LBW"),
  recomputed = c(uniqueN(raw$IID), sum(!is.na(bmi$BMI)), sum(o$PTB==1, na.rm=TRUE),
                 sum(pe$pe_pos), sum(ch$chtn_pos), sum(tw$twin_pos),
                 sum(o$BOUT==1, na.rm=TRUE), sum(o$BOUT==2, na.rm=TRUE),
                 sum(o$LBW==1, na.rm=TRUE)),
  previously_recorded = c(21685, 21118, 2262, 413, 488, 209, NA, NA, NA))
tot[, matches := fifelse(is.na(previously_recorded), NA, recomputed==previously_recorded)]
pr(tot); sav(tot,"headline_counts")
cat("\nAny FALSE in `matches` means the build or the source file has changed.\n")

## ============================================================
hdr("J. ICC BY COHORT x TRAIT, RECOMPUTED")
lng <- rbindlist(list(
  raw[is.finite(SBPn), .(IID, cohort, trait="SBP", bp=SBPn)],
  raw[is.finite(DBPn), .(IID, cohort, trait="DBP", bp=DBPn)]))
ic <- lng[, .(icc=round(momi_icc(.SD, id="IID", y="bp"),3),
              n_readings=.N, n_mothers=uniqueN(IID)), by=.(cohort,trait)][order(cohort,trait)]
cat("ICC is the reliability of a single reading. Low ICC => a single reading is noisy =>\n",
    "the mean of several should transfer better. This is the mechanism F3 found.\n\n", sep="")
pr(ic); sav(ic,"icc")

## ============================================================
hdr("K. BP-DEFINITION COVERAGE BY COHORT")
cat("How many mothers actually HAVE each definition? A definition that exists for 8% of a\n",
    "cohort is not comparable to one that exists for 99%.\n\n", sep="")
setorder(raw, IID, GAok, na.last=TRUE)
wmn <- function(v, keep){ sel <- keep & !is.na(v); sel[is.na(sel)] <- FALSE
                          if(!any(sel)) NA_real_ else mean(v[sel]) }
defs <- raw[, .(first=fnn(SBPn), mean=if(all(is.na(SBPn))) NA_real_ else mean(SBPn,na.rm=TRUE),
                lt20=wmn(SBPn, GAok < MOMI_GA$wk20), ge20=wmn(SBPn, GAok >= MOMI_GA$wk20),
                tri1=wmn(SBPn, GAok < MOMI_GA$tri1_end),
                tri3=wmn(SBPn, GAok >= MOMI_GA$tri2_end)), by=.(cohort,IID)]
cov <- defs[, .(mothers=.N,
                pct_first=round(100*mean(!is.na(first)),1),
                pct_mean =round(100*mean(!is.na(mean)),1),
                pct_lt20 =round(100*mean(!is.na(lt20)),1),
                pct_ge20 =round(100*mean(!is.na(ge20)),1),
                pct_tri1 =round(100*mean(!is.na(tri1)),1),
                pct_tri3 =round(100*mean(!is.na(tri3)),1)), by=cohort][order(cohort)]
pr(cov); sav(cov,"definition_coverage")

## ============================================================
hdr("L. GENOTYPED — DEFINED BY PRESENCE OF ANY GENETIC DATA")
cat("B01 currently flags 'genotyped' using ONE score (", MOMI_EUR$SBP, ", EUR SBP). Decision\n",
    "2026-07-19: it should mean simply 'this mother has genetic data'. This section scans\n",
    "EVERY .sscore file rather than assuming, so we can see whether the definitions differ\n",
    "and by how much.\n", sep="")
if(is.null(SSC) || !nzchar(SSC)){
  cat("\nSSC not set — skipping.\n")
} else {
  allf <- list.files(SSC, pattern="\\.sscore$", full.names=FALSE)
  cat(sprintf("\n.sscore files found: %d\n", length(allf)))

  ## ---- filename -> score / cohort / platform / infant ----
  parts <- strsplit(sub("\\.sscore$","",allf), "__", fixed=TRUE)
  FI <- data.table(
    file     = allf,
    score    = vapply(parts, function(p) if(length(p)>=1) p[1] else NA_character_, ""),
    cohort   = vapply(parts, function(p) if(length(p)>=2) p[2] else NA_character_, ""),
    platform = vapply(parts, function(p) if(length(p)>=3) p[3] else NA_character_, ""),
    infant   = grepl("__infant", allf, fixed=TRUE))

  sec("L1. inventory — what actually exists on disk")
  print(FI[, .N, by=.(platform, infant)][order(platform, infant)])
  cat("\nscores present:\n"); print(sort(unique(FI$score)))
  cat("\ncohorts present:\n"); print(sort(unique(FI$cohort)))

  sec("L2. INFANT files — does S15 (fetal vs maternal) actually have data?")
  if(any(FI$infant)){
    cat(">>> INFANT SCORE FILES EXIST. S15 may NOT be blocked after all.\n\n")
    print(FI[infant==TRUE, .N, by=.(cohort, platform)][order(cohort, platform)])
  } else cat("no __infant files — S15 remains blocked on infant genotypes.\n")

  ## ---- read IIDs from every maternal file ----
  read_ids <- function(f){
    d <- tryCatch(fread(file.path(SSC,f), select=1), error=function(e) NULL)
    if(is.null(d)) return(character(0)); unique(as.character(d[[1]]))
  }
  MAT <- FI[infant==FALSE]
  cat(sprintf("\nreading %d maternal .sscore files ...\n", nrow(MAT)))
  idlist <- lapply(MAT$file, read_ids)
  MAT[, n_samples := vapply(idlist, length, integer(1))]

  sec("L3. sample count per cohort x platform (should be identical across scores)")
  cp <- MAT[, .(files=.N, n_min=min(n_samples), n_max=max(n_samples),
                consistent = min(n_samples)==max(n_samples)),
            by=.(cohort, platform)][order(cohort, platform)]
  pr(cp); sav(cp, "sscore_cohort_platform")
  cat("\nconsistent==FALSE anywhere means different scores cover different samples --\n",
      "which is exactly the fragility the single-score flag was exposed to.\n", sep="")

  sec("L4. competing definitions of 'genotyped'")
  usable   <- MAT$platform %in% MOMI_PLATS
  ids_eur  <- unique(unlist(idlist[MAT$score==MOMI_EUR$SBP & usable]))
  ids_use  <- unique(unlist(idlist[usable]))            # any score, usable platforms
  ids_any  <- unique(unlist(idlist))                    # any score, ANY platform incl. raw lpwgs
  D <- data.table(
    definition = c(sprintf("current B01: %s only, usable platforms", MOMI_EUR$SBP),
                   "ANY score, usable platforms (gsa + lpwgs_dosage)",
                   "ANY score, ANY platform (includes raw lpwgs)"),
    n_IIDs = c(length(ids_eur), length(ids_use), length(ids_any)))
  pr(D); sav(D, "genotyped_definitions")
  cat(sprintf("\ngained by moving single-score -> any-score (usable): %d\n",
              length(setdiff(ids_use, ids_eur))))
  cat(sprintf("gained by additionally allowing raw lpwgs        : %d\n",
              length(setdiff(ids_any, ids_use))))
  cat("\nNOTE a real distinction: a mother on raw lpwgs only HAS genetic data but has NO\n",
      "usable polygenic score, so she cannot contribute to transferability or MR. If that\n",
      "count is non-zero, F1's flow needs TWO steps -- 'genotyped' then 'has usable score' --\n",
      "rather than one.\n", sep="")

  sec("L5. do score IIDs match EPI mothers?")
  epi_ids <- unique(raw$IID)
  for(nm in c("ids_eur","ids_use","ids_any")){
    v <- get(nm)
    cat(sprintf("%-8s: %6d IIDs | matched to EPI first-preg: %6d | unmatched: %6d\n",
                nm, length(v), length(intersect(v, epi_ids)), length(setdiff(v, epi_ids))))
  }
  cat("\nUnmatched should be ~0. A large number means an ID-format mismatch between the\n",
      "genotype files and the EPI file, which would silently deflate every genotyped count.\n", sep="")
  um <- setdiff(ids_any, epi_ids)
  if(length(um)) { cat("\nexample unmatched IIDs:\n"); print(utils::head(um, 10)) }
}

hdr("AUDIT COMPLETE")
cat("Per-section TSVs written to:", QC, "\n")
