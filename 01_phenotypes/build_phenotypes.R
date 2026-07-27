#!/usr/bin/env Rscript
# ============================================================
# build_phenotypes.R
# Stage 01: from merged.fam + Epi file, build the analytic-mother set and all
# phenotype/covariate files the downstream expects. Keys entirely on
# PARTICIPANT_ID (stage-00 already renamed every sample to its Epi PID), so the
# old ZAPPS/-M/ORIG_ID normalization is no longer needed.
#
# Outputs (in --out-dir):
#   mothers.keep                              FID IID  (analytic mothers)
#   analytic_mothers.txt                      PIDs in geno AND pheno (pre-filter)
#   covariate_table_analytic_mothers.txt      PARTICIPANT_ID + covars + SITE + ANCESTRY
#   mothers_<TRAIT>_<window>_final.pheno      FID IID PHENO  (SBP/DBP x 5 windows)
#   mothers_PTB_NEW_final.pheno               FID IID PHENO  (1=control,2=case)
#   phenotype_build_report.txt                exclusion counts + summary
#
# Analytic filters (restricted to PREGNANCY_ID==1; consistency fix vs the old
# pipeline which restricted mother-selection but not BP-window extraction):
#   first pregnancy, SINGLE_TWIN==1, PREG_OUTCOME==2, BIRTH_OUTCOME==1,
#   PE_CAT=="NA" (no preeclampsia), and >=1 valid BP visit at GA>20wk (vis<=del).
#
# Usage:
#   Rscript build_phenotypes.R --merged-fam FILE --epi FILE --out-dir DIR
#           [--ptb-case-code 2]
# ============================================================
suppressMessages(library(data.table))

## ---- args ----
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default=NULL){ i <- which(args==flag); if(length(i)) args[i+1] else default }
FAM    <- getarg("--merged-fam")
EPI    <- getarg("--epi")
OUTDIR <- getarg("--out-dir")
PTB_CASE <- as.integer(getarg("--ptb-case-code", "2"))
stopifnot(!is.null(FAM), !is.null(EPI), !is.null(OUTDIR))
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

MISS <- c("-88","-77","-99","NA","na","",".")
is_miss <- function(x) is.na(x) | x %in% MISS

## prefix -> cohort/ancestry (matches config.cohorts; extend here if cohorts change)
PFX <- data.table(
  prefix   = c("AMANHIB","AMANHIP","AMANHIT","GAPPSB","ZAPPS"),
  cohort   = c("AMANHI-Bangladesh","AMANHI-Pakistan","AMANHI-Pemba","GAPPS-Bangladesh","GAPPS-Zambia"),
  ancestry = c("SAS","SAS","AFR","SAS","AFR"))
cohort_of <- function(pid){
  p <- sub("[0-9].*$","",pid); p <- sub("-$","",p)
  # longest-prefix match
  out <- rep(NA_character_, length(pid)); anc <- out
  for(i in seq_len(nrow(PFX))){
    hit <- startsWith(pid, paste0(PFX$prefix[i]))
    out[hit] <- PFX$cohort[i]; anc[hit] <- PFX$ancestry[i]
  }
  list(cohort=out, ancestry=anc)
}

## ---- robust, locale-independent date parser ----
## handles YYYY-MM-DD  and  DD-MON-YYYY / DDMONYYYY (e.g. 05JAN2018, 5-Jan-2018)
parse_date <- function(x){
  x <- toupper(gsub('[" \t]',"",as.character(x)))
  out <- as.Date(rep(NA_real_, length(x)), origin="1970-01-01")
  iso <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x)
  out[iso] <- as.Date(x[iso], "%Y-%m-%d")
  mon <- c(JAN="01",FEB="02",MAR="03",APR="04",MAY="05",JUN="06",
           JUL="07",AUG="08",SEP="09",OCT="10",NOV="11",DEC="12")
  mm <- regmatches(x, regexec("^([0-9]{1,2})-?([A-Z]{3})-?([0-9]{4})$", x))
  for(i in which(is.na(out))){
    g <- mm[[i]]
    if(length(g)==4 && g[3] %in% names(mon))
      out[i] <- as.Date(paste(g[4], mon[g[3]], sprintf("%02d",as.integer(g[2])), sep="-"), "%Y-%m-%d")
  }
  out
}

## ---- load ----
fam <- fread(FAM, header=FALSE)
setnames(fam, 1:2, c("FID","IID"))
geno_pids <- unique(fam$IID)

# na.strings=NULL: keep literal "NA" as a string — PE_CAT=="NA" MEANS "no preeclampsia"
# (the category to keep). Missingness is handled explicitly via is_miss().
epi <- fread(EPI, sep="\t", header=TRUE, colClasses="character", quote="", na.strings=NULL)
need <- c("PARTICIPANT_ID","PREGNANCY_ID","SINGLE_TWIN","PREG_OUTCOME","BIRTH_OUTCOME",
          "PE_CAT","DEL_DATE","VISITDT","SBP","DBP","GA_HDLK_NEW","PTB_NEW","SITE_CODE",
          "PW_AGE","PW_EDUCATION","WEALTH_INDEX","PARITY","CHRON_HTN","DIABETES",
          "MAT_HEIGHT","SMOK_FREQ","BABY_SEX")
miss_cols <- setdiff(need, names(epi))
if(length(miss_cols)) stop("Epi missing columns: ", paste(miss_cols, collapse=", "))

epi <- epi[PARTICIPANT_ID %in% geno_pids]            # restrict to genotyped mothers
epi1 <- epi[PREGNANCY_ID=="1"]                        # first pregnancy only
analytic <- sort(unique(epi1$PARTICIPANT_ID))
writeLines(analytic, file.path(OUTDIR,"analytic_mothers.txt"))

## ---- per-mother analytic filters ----
epi1[, ga_weeks := as.integer((suppressWarnings(as.numeric(GA_HDLK_NEW))+6)/7)]
epi1[, vis := parse_date(VISITDT)]
epi1[, del := parse_date(DEL_DATE)]
epi1[, sbp_ok := !is_miss(SBP) & suppressWarnings(as.numeric(SBP))>0]
epi1[, dbp_ok := !is_miss(DBP) & suppressWarnings(as.numeric(DBP))>0]
epi1[, validBPvisit := !is.na(vis) & !is.na(del) & vis<=del & !is.na(ga_weeks) & ga_weeks>20 & (sbp_ok|dbp_ok)]

agg <- epi1[, .(
  single    = any(SINGLE_TWIN=="1", na.rm=TRUE),
  livepreg  = any(PREG_OUTCOME=="2", na.rm=TRUE),
  livebirth = any(BIRTH_OUTCOME=="1", na.rm=TRUE),
  nope      = any(is.na(PE_CAT) | PE_CAT=="NA", na.rm=TRUE),   # "NA"/missing PE_CAT = no preeclampsia (keep)
  validBP   = any(validBPvisit, na.rm=TRUE)
), by=PARTICIPANT_ID]

excl <- list()
keep <- agg$PARTICIPANT_ID
# apply exclusions in order on the shrinking set
a <- agg[match(keep, PARTICIPANT_ID)]
keep <- a[single==TRUE,   PARTICIPANT_ID]; excl$non_singleton <- a[single!=TRUE,PARTICIPANT_ID]
a <- a[PARTICIPANT_ID %in% keep]
keep <- a[livepreg==TRUE, PARTICIPANT_ID]; excl$non_livepreg  <- a[livepreg!=TRUE,PARTICIPANT_ID]
a <- a[PARTICIPANT_ID %in% keep]
keep <- a[livebirth==TRUE,PARTICIPANT_ID]; excl$non_livebirth <- a[livebirth!=TRUE,PARTICIPANT_ID]
a <- a[PARTICIPANT_ID %in% keep]
keep <- a[nope==TRUE,     PARTICIPANT_ID]; excl$preeclampsia  <- a[nope!=TRUE,PARTICIPANT_ID]
a <- a[PARTICIPANT_ID %in% keep]
keep <- a[validBP==TRUE,  PARTICIPANT_ID]; excl$no_validBP    <- a[validBP!=TRUE,PARTICIPANT_ID]
mothers <- sort(unique(keep))

for(nm in names(excl)) if(length(excl[[nm]])) writeLines(excl[[nm]], file.path(OUTDIR, paste0("excluded_",nm,".txt")))

## ---- mothers.keep (FID IID) ----
keepdt <- fam[IID %in% mothers, .(FID, IID)]
fwrite(keepdt, file.path(OUTDIR,"mothers.keep"), sep="\t", col.names=FALSE)

## ---- covariate table (first preg==1 record per mother) ----
cov <- epi1[PARTICIPANT_ID %in% mothers]
cov <- cov[order(PARTICIPANT_ID)][, .SD[1], by=PARTICIPANT_ID]
ca <- cohort_of(cov$PARTICIPANT_ID)
covout <- data.table(
  PARTICIPANT_ID = cov$PARTICIPANT_ID,
  SITE   = ca$cohort, ANCESTRY = ca$ancestry, SITE_CODE = cov$SITE_CODE,
  PW_AGE = cov$PW_AGE, PW_EDUCATION = cov$PW_EDUCATION, WEALTH_INDEX = cov$WEALTH_INDEX,
  PARITY = cov$PARITY, CHRON_HTN = cov$CHRON_HTN, DIABETES = cov$DIABETES,
  MAT_HEIGHT = cov$MAT_HEIGHT, SMOK_FREQ = cov$SMOK_FREQ, BABY_SEX = cov$BABY_SEX)
for(j in names(covout)) covout[is_miss(get(j)), (j):=NA]
fwrite(covout, file.path(OUTDIR,"covariate_table_analytic_mothers.txt"), sep="\t", na="NA")

## ---- BP-window phenotypes ----
fammap <- setNames(fam$FID, fam$IID)
write_pheno <- function(iids, val, fn){
  dt <- data.table(FID=fammap[iids], IID=iids,
                    PHENO=ifelse(is.na(val),"NA",sprintf("%.1f",val)))
  setnames(dt, c("FID","IID","PHENO"))
  fwrite(dt, file.path(OUTDIR,fn), sep="\t")
}
bp <- epi1[PARTICIPANT_ID %in% mothers & !is.na(vis) & !is.na(del) & (vis<=del)]
bp[, days_before := as.numeric(del - vis)]
bp <- bp[days_before>=0]

windows <- c("earliest","5mo","7mo","8mo_and_later","latest")
for(TR in c("SBP","DBP")){
  bp[, v := suppressWarnings(as.numeric(get(TR)))]
  bp[is_miss(get(TR)), v := NA]
  d <- bp[!is.na(v)]
  for(w in windows){
    val <- sapply(mothers, function(m){
      vv <- d[PARTICIPANT_ID==m]
      if(nrow(vv)==0) return(NA_real_)
      if(w=="earliest") return(vv[which.min(vis), v])
      if(w=="latest")   return(vv[which.max(vis), v])
      if(w=="5mo"){ s<-vv[days_before>=90 & days_before<=210]; if(nrow(s)) return(s[which.min(abs(days_before-150)),v]) else return(NA_real_)}
      if(w=="7mo"){ s<-vv[days_before>=180 & days_before<=240]; if(nrow(s)) return(s[which.min(abs(days_before-210)),v]) else return(NA_real_)}
      if(w=="8mo_and_later"){ s<-vv[days_before>=0 & days_before<=30]; if(nrow(s)) return(s[which.max(vis),v]) else return(NA_real_)}
      NA_real_
    })
    write_pheno(mothers, as.numeric(val), sprintf("mothers_%s_%s_final.pheno",TR,w))
  }
}

## ---- PTB phenotype (plink: 1=control, 2=case) ----
ptb <- epi1[PARTICIPANT_ID %in% mothers][order(PARTICIPANT_ID)][, .SD[1], by=PARTICIPANT_ID]
ptb[, pv := PTB_NEW]
ptb[, code := ifelse(is_miss(pv), "NA", ifelse(pv=="1", as.character(PTB_CASE), as.character(3-PTB_CASE)))]
ptbout <- data.table(FID=fammap[ptb$PARTICIPANT_ID], IID=ptb$PARTICIPANT_ID, PHENO=ptb$code)
setnames(ptbout, c("FID","IID","PHENO"))
fwrite(ptbout, file.path(OUTDIR,"mothers_PTB_NEW_final.pheno"), sep="\t")

## ---- report ----
rep <- file.path(OUTDIR,"phenotype_build_report.txt")
con <- file(rep,"w")
wl <- function(...) writeLines(sprintf(...), con)
wl("=== Stage 01 phenotype build ===")
wl("genotyped mothers (merged.fam IIDs)   : %d", length(geno_pids))
wl("in Epi (first pregnancy)               : %d", length(analytic))
wl("excluded non-singleton                 : %d", length(excl$non_singleton))
wl("excluded non-livepreg                  : %d", length(excl$non_livepreg))
wl("excluded non-livebirth                 : %d", length(excl$non_livebirth))
wl("excluded preeclampsia (PE_CAT!=NA)     : %d", length(excl$preeclampsia))
wl("excluded no valid BP after 20wk        : %d", length(excl$no_validBP))
wl("ANALYTIC MOTHERS retained (mothers.keep): %d", length(mothers))
ptb_tab <- table(factor(ptbout$PHENO, levels=c("1","2","NA")))
wl("PTB controls(1)/cases(2)/NA            : %d / %d / %d", ptb_tab["1"], ptb_tab["2"], ptb_tab["NA"])
ca2 <- cohort_of(mothers)
wl("by cohort:")
for(co in sort(unique(ca2$cohort))) wl("  %-20s %d", co, sum(ca2$cohort==co, na.rm=TRUE))
nlatest <- function(TR){ p<-fread(file.path(OUTDIR,sprintf("mothers_%s_latest_final.pheno",TR))); sum(p$PHENO!="NA") }
wl("non-missing SBP(latest)/DBP(latest)    : %d / %d", nlatest("SBP"), nlatest("DBP"))
close(con)
cat(readLines(rep), sep="\n"); cat("\n")
