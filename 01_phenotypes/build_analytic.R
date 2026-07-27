#!/usr/bin/env Rscript
# ============================================================
# build_analytic.R  [ build step B01 ]
# Canonical per-mother analytic table for the whole paper. EVERY deliverable reads
# the intermediates written here; no downstream module re-reads raw EPI. This is the
# synchronization point that stops Table 2 and Figure 2 (etc.) from drifting apart.
#
# Schema is grounded in the real EPI file (41 cols, 173,113 rows of which 167,460 are
# first-pregnancy visits, 21,685 first-preg mothers over 6 sites). Notable data facts:
#   * NO BMI column -> compute from MAT_HEIGHT(cm) & MAT_WEIGHT(kg), earliest visit.
#   * Use precomputed PTB_NEW (agrees 100% with GAGEBRTH<259, adds QC NAs).
#   * PE = EOPE|LOPE (PE_CAT); drop TBD. Carry CHRON_HTN (chronic-HTN proxy).
#   * livebirth = BIRTH_OUTCOME==1 ; stillbirth = 2.
#   * confounders: AGE, GRAVIDITY, PW_EDUCATION, computed BMI (+ PARITY/WEALTH kept
#     for reference/sensitivity).
#
# ------------------------------------------------------------------
# MAJOR REVISION 2026-07-19, from audit_B01.R. Two changes, both consequential:
#
# 1. ANTENATAL RESTRICTION. ~40% of every AMANHI cohort's BP readings are POSTPARTUM
#    (visit date after delivery date) against 0% in GAPPS. 21,032 of the 21,087 readings
#    formerly discarded as "GA > 315 d errors" were not errors at all -- they were
#    postnatal visits (99.7%). Worse, a further 22,889 postnatal visits carry a GA
#    INSIDE [28,315] and so looked perfectly valid: a visit 2 wk after a 34-wk delivery
#    computes to 252 d and was being used as a third-trimester antenatal measurement.
#    Postpartum BP is 2.8-3.9 mmHg higher (pregnancy vasodilation resolving), so this
#    inflated mean SBP by 1.15-1.45 mmHg in AMANHI and 0.00 in GAPPS -- a DIFFERENTIAL
#    bias in the primary exposure, in exactly the cohorts with lower transferability.
#    => All BP definitions are now computed on ANTENATAL visits only.
#    => Parallel POSTNATAL columns (S_*_pn / D_*_pn) are retained, because the BP scores
#       were trained on non-pregnant adults and postpartum BP is closer to that state.
#       If the PRS predicts postpartum BP better than antenatal BP, that is direct
#       evidence that pregnancy physiology -- not ancestry -- degrades portability.
#
# 2. GA LOWER BOUND DROPPED. Only 50 readings fell below 28 d, so the floor bought
#    nothing; it is now [0, 315], excluding only impossible negatives (THSTI has values
#    to -596). The 315 ceiling is retained as a sanity bound but is nearly non-binding
#    once postnatal visits are removed.
#
# Delivery date is resolved through a fallback cascade before classifying visits, and a
# reading is discarded only if all fallbacks fail (decision: "try hard, then delete").
# ------------------------------------------------------------------
#
# Writes (results/current/intermediates/):
#   analytic_mothers.rds, bp_readings_long.rds, bp_icc.rds
#   Rscript build_analytic.R --epi EPI [--sscore-dir DIR] [--pipe PIPE]
# ============================================================
suppressMessages({library(data.table); library(splines)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R"))
source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))   # momi_icc()
P <- momi_paths(PIPE)

# EPI/SSC from CLI flag, else the environment (driver exports them; env is space-safe).
EPI <- momi_arg("--epi", Sys.getenv("EPI")); SSC <- momi_arg("--sscore-dir", Sys.getenv("SSC"))
if(is.null(EPI) || !nzchar(EPI)) stop("no EPI: pass --epi or export EPI")
if(is.null(SSC) || !nzchar(SSC)) SSC <- NULL
num <- function(x) suppressWarnings(as.numeric(x))
firstnn <- function(v){ w <- which(!is.na(v)); if(!length(w)) v[NA_integer_][1] else v[w[1]] }
asdate <- function(x){
  if(inherits(x, c("Date","IDate"))) return(as.Date(x))
  suppressWarnings(as.Date(as.character(x), format="%Y-%m-%d"))
}

momi_deliverable("B01_build_analytic", script="01_phenotypes/build_analytic.R",
                 inputs=paste0("epi=",basename(EPI)), P=P, body=function(ctx){

  raw <- fread(EPI, na.strings=MOMI_NA)[PREGNANCY_ID==1]
  getcol <- function(nm) if(nm %in% names(raw)) num(raw[[nm]]) else rep(NA_real_, nrow(raw))
  raw[, IID    := as.character(PARTICIPANT_ID)]
  raw[, SBP    := getcol("SBP")]
  raw[, DBP    := getcol("DBP")]
  garead       <- getcol("GA_HDLK_NEW")
  ## GA sanity bound: negatives are impossible; 315 d (45 wk) is a ceiling no real
  ## antenatal visit reaches. NB this NULLS the GA only -- the BP value survives and is
  ## still used by the GA-free definitions (mean/first/last/median).
  raw[, GA     := ifelse(garead >= 0 & garead <= 315, garead, NA_real_)]
  raw[, VIS    := getcol("VISITDT")]
  raw[, AGEv   := getcol("PW_AGE")]
  raw[, GRAVv  := getcol("GRAVIDITY")]
  raw[, PARv   := getcol("PARITY")]
  raw[, EDUv   := getcol("PW_EDUCATION")]
  raw[, WEALv  := getcol("WEALTH_INDEX")]
  raw[, HTv    := getcol("MAT_HEIGHT")]
  raw[, WTv    := getcol("MAT_WEIGHT")]
  raw[, CHTNv  := getcol("CHRON_HTN")]
  raw[, GAdv   := getcol("GAGEBRTH_NEW")]
  raw[, BWTv   := getcol("BIRTH_WEIGHT")]
  raw[, SGAv   := getcol("SGA_10_NEW")]
  raw[, SGA3v  := getcol("SGA_3_NEW")]
  raw[, PTBv   := getcol("PTB_NEW")]
  raw[, SPONTv := getcol("SPONT_LABOUR")]
  raw[, TWINv  := getcol("SINGLE_TWIN")]
  raw[, BOUTv  := getcol("BIRTH_OUTCOME")]
  raw[, PEv    := if("PE_CAT" %in% names(raw)) as.character(PE_CAT) else NA_character_]
  raw[, cohort := MOMI_SITE[as.character(SITE_CODE)]]
  setorder(raw, IID, VIS, GA)

  ## ---------- ANTENATAL / POSTNATAL CLASSIFICATION (see header, change 1) ----------
  raw[, VISd := asdate(VISITDT)]
  raw[, DELd := if("DEL_DATE"  %in% names(raw)) asdate(DEL_DATE)  else as.Date(NA)]
  raw[, LMPd := if("DATE_LMP"  %in% names(raw)) asdate(DATE_LMP)  else as.Date(NA)]
  raw[, BWTd := if("BWT_MEASURE_DATE" %in% names(raw)) asdate(BWT_MEASURE_DATE) else as.Date(NA)]

  ## fallback cascade for the delivery date, most reliable first
  raw[, DELeff := DELd]
  n_direct <- sum(!is.na(raw$DELeff))
  raw[is.na(DELeff) & !is.na(LMPd) & is.finite(GAdv), DELeff := LMPd + GAdv]
  n_lmp    <- sum(!is.na(raw$DELeff)) - n_direct
  raw[is.na(DELeff) & !is.na(BWTd), DELeff := BWTd]
  n_bwt    <- sum(!is.na(raw$DELeff)) - n_direct - n_lmp
  ## delivery date is a MOTHER-level fact: propagate the first known value within IID
  raw[, DELeff := { v <- DELeff[!is.na(DELeff)]; if(length(v)) v[1] else as.Date(NA) }, by=IID]
  n_prop   <- sum(!is.na(raw$DELeff)) - n_direct - n_lmp - n_bwt

  raw[, visit_class := fifelse(is.na(VISd) | is.na(DELeff), "unknown",
                        fifelse(VISd > DELeff, "postnatal", "antenatal"))]

  ## ---- Rule 4: LMP rescue (decision 2026-07-19) ----
  ## The cascade above provably failed for 442 mothers -- diagnostics showed 0 of them have
  ## GAGEBRTH_NEW or a birth-weight date, i.e. NO delivery record of any kind (lost to
  ## follow-up; 244 of the 442 are Zambian, 9.8% of that cohort). But 1,131 of their 1,141
  ## rows DO carry an LMP date, and a postnatal visit presupposes a delivery -- so with no
  ## delivery on record these visits are very likely antenatal. Classify them on LMP-based
  ## GA rather than discarding genuine antenatal BP from the paper's weakest cohort.
  ## Ceiling 294 d (42 wk) excludes the obvious postnatal cases; residual misclassification
  ## risk is real but small, and the logged distribution below makes it visible.
  raw[, rescued_lmp := FALSE]
  raw[, GA_lmp := as.numeric(VISd - LMPd)]
  n_unk_before <- sum(raw$visit_class == "unknown")
  raw[visit_class=="unknown" & is.finite(GA_lmp) & GA_lmp >= 0 & GA_lmp <= 294,
      `:=`(visit_class = "antenatal", rescued_lmp = TRUE)]
  n_resc <- sum(raw$rescued_lmp)
  cat(sprintf("LMP rescue: %d of %d unknown rows -> antenatal (%d mothers); %d rows still unknown\n",
              n_resc, n_unk_before, uniqueN(raw[rescued_lmp==TRUE]$IID),
              sum(raw$visit_class=="unknown")))
  if(n_resc){
    cat("  LMP-GA (days) of rescued rows -- values near 294 are the misclassification risk:\n")
    print(round(quantile(raw[rescued_lmp==TRUE]$GA_lmp, c(0,.25,.5,.75,.95,1), na.rm=TRUE), 0))
    print(raw[rescued_lmp==TRUE, .(rows=.N, mothers=uniqueN(IID)), by=cohort][order(-rows)])
  }

  vc <- raw[, .N, by=visit_class][order(visit_class)]
  cat("visit classification:\n"); print(vc)
  cat(sprintf("delivery date resolved: direct=%d +LMP+GAbirth=%d +BWTdate=%d +propagated=%d ; unresolved rows=%d\n",
              n_direct, n_lmp, n_bwt, n_prop, sum(is.na(raw$DELeff))))
  ## Prove the cascade actually tried. A fallback yielding 0 is ambiguous: it looks the same
  ## whether the source data is absent or the fallback silently failed. Show the inputs.
  if(sum(is.na(raw$DELeff))){
    u <- raw[is.na(DELeff)]
    cat(sprintf("  unresolved diagnostics: rows=%d mothers=%d | has LMP date=%d | has GAGEBRTH=%d | has BOTH (LMP fallback eligible)=%d | has BWT date=%d\n",
                nrow(u), uniqueN(u$IID), sum(!is.na(u$LMPd)), sum(is.finite(u$GAdv)),
                sum(!is.na(u$LMPd) & is.finite(u$GAdv)), sum(!is.na(u$BWTd))))
    cat("  if 'has BOTH' > 0 the LMP fallback is BROKEN, not merely unlucky.\n")
    print(u[, .(rows=.N, mothers=uniqueN(IID)), by=cohort][order(-rows)])
  }
  cat("BP readings by class (SBP present):\n")
  print(raw[is.finite(SBP), .N, by=.(cohort, visit_class)][order(cohort, visit_class)])

  ## Decision 2026-07-19: unknown-class readings are DELETED from every BP definition
  ## (only positively-classified visits are used). Covariates and outcomes below still
  ## read all rows, because those are mother-level facts, not visit-level measurements.
  ANTE <- raw[visit_class == "antenatal"]
  POST <- raw[visit_class == "postnatal"]

  ## ---------- per-visit long BP (ICC + BP-by-GA diagnostics) ----------
  ## carries visit_class so downstream figures (e.g. F3 panel A's BP-by-GA trajectory)
  ## cannot silently plot postpartum readings as antenatal.
  long <- rbindlist(list(
    raw[is.finite(SBP), .(IID, cohort, GA_days=GA, trait="SBP", bp=SBP, visit_class)],
    raw[is.finite(DBP), .(IID, cohort, GA_days=GA, trait="DBP", bp=DBP, visit_class)]))
  momi_save_intermediate(long, "bp_readings_long", P)

  ## ---------- reliability ICC per cohort x trait (ANTENATAL only) ----------
  icc <- long[visit_class=="antenatal",
              .(icc=momi_icc(.SD, id="IID", y="bp"),
                n_readings=.N, n_mothers=uniqueN(IID)), by=.(cohort, trait)]
  momi_save_intermediate(icc, "bp_icc", P)

  ## ---------- BP definitions per mother ----------
  fv  <- function(v){w<-which(!is.na(v)); if(!length(w)) NA_real_ else v[w[1]]}
  lv  <- function(v){w<-which(!is.na(v)); if(!length(w)) NA_real_ else v[w[length(w)]]}
  m2  <- function(v){h<-head(v[!is.na(v)],2); if(!length(h)) NA_real_ else mean(h)}
  mnv <- function(v) if(all(is.na(v))) NA_real_ else mean(v,na.rm=TRUE)
  mdv <- function(v) if(all(is.na(v))) NA_real_ else as.numeric(median(v,na.rm=TRUE))
  wm  <- function(v,keep){ sel<-keep & !is.na(v); sel[is.na(sel)]<-FALSE
                           if(!any(sel)) NA_real_ else mean(v[sel]) }
  g1<-MOMI_GA$tri1_end; g2<-MOMI_GA$tri2_end; w20<-MOMI_GA$wk20

  ## ANTENATAL definitions (the paper's exposures)
  bp_by <- function(D, bpcol) D[, .(
      first=fv(get(bpcol)), mean=mnv(get(bpcol)), mn2=m2(get(bpcol)),
      median=mdv(get(bpcol)), last=lv(get(bpcol)),
      tri1=wm(get(bpcol),GA<g1), tri2=wm(get(bpcol),GA>=g1&GA<g2), tri3=wm(get(bpcol),GA>=g2),
      lt20=wm(get(bpcol),GA<w20), ge20=wm(get(bpcol),GA>=w20)), by=IID]
  Sd<-bp_by(ANTE,"SBP"); setnames(Sd,setdiff(names(Sd),"IID"),paste0("S_",setdiff(names(Sd),"IID")))
  Dd<-bp_by(ANTE,"DBP"); setnames(Dd,setdiff(names(Dd),"IID"),paste0("D_",setdiff(names(Dd),"IID")))

  ## POSTNATAL companion columns (see header, change 1): postpartum BP is closer to the
  ## non-pregnant state the BP scores were trained in. Suffix _pn throughout.
  pn_by <- function(D, bpcol) D[, .(
      mean=mnv(get(bpcol)), first=fv(get(bpcol)), median=mdv(get(bpcol)), n=sum(!is.na(get(bpcol)))),
      by=IID]
  Sp<-pn_by(POST,"SBP"); setnames(Sp,setdiff(names(Sp),"IID"),paste0("S_",setdiff(names(Sp),"IID"),"_pn"))
  Dp<-pn_by(POST,"DBP"); setnames(Dp,setdiff(names(Dp),"IID"),paste0("D_",setdiff(names(Dp),"IID"),"_pn"))

  ## ---------- GA-residual (GA-standardized) definition, ANTENATAL only ----------
  ga_resid_mean <- function(tr){
    dd <- long[trait==tr & visit_class=="antenatal" & is.finite(bp) & is.finite(GA_days)]
    if(!nrow(dd)) return(data.table(IID=character(0), resid=numeric(0)))
    out <- dd[, {
      m <- tryCatch(lm(bp ~ ns(GA_days,3)), error=function(e) NULL)
      r <- if(is.null(m)) rep(NA_real_,.N) else resid(m)+mean(bp,na.rm=TRUE)
      .(IID=IID, r=r) }, by=cohort]
    out[, .(resid=mean(r,na.rm=TRUE)), by=IID]
  }
  Sr<-ga_resid_mean("SBP"); setnames(Sr,"resid","S_resid")
  Dr<-ga_resid_mean("DBP"); setnames(Dr,"resid","D_resid")

  ## ---------- outcomes + covariates (per mother; all rows, mother-level facts) --------
  cov <- raw[, .(
    cohort=cohort[1], SITE=SITE_CODE[1],
    AGE=firstnn(AGEv), GRAV=firstnn(GRAVv), PAR=firstnn(PARv),
    EDU=firstnn(EDUv), WEALTH=firstnn(WEALv),
    HT=firstnn(HTv), WT=firstnn(WTv),
    CHRON_HTN=as.integer(any(CHTNv==1, na.rm=TRUE)),
    PTB=firstnn(PTBv),
    GAd=firstnn(GAdv), BWTraw=firstnn(BWTv),
    SGA=as.integer(firstnn(SGAv)==1), SGA3=as.integer(firstnn(SGA3v)==1),
    SPONT=firstnn(SPONTv), TWIN=firstnn(TWINv), BOUT=firstnn(BOUTv),
    PE=as.integer(any(PEv %in% c("EOPE","LOPE"))),
    n_visits=.N, n_ante=sum(visit_class=="antenatal"),
    n_post=sum(visit_class=="postnatal"), n_unknown=sum(visit_class=="unknown"),
    n_rescued_lmp=sum(rescued_lmp)   # antenatal only via the LMP rule; keeps the rescue auditable
  ), by=IID]
  cov[, ancestry := MOMI_ANC[cohort]]
  cov[, BMI := { b <- WT/(HT/100)^2; ifelse(b>=12 & b<=60, b, NA_real_) }]        # computed
  cov[, GAd := ifelse(GAd>=140 & GAd<=320, GAd, NA_real_)]                        # clip birth GA
  cov[, BWT := ifelse(BWTraw>=500 & BWTraw<=6500, BWTraw, NA_real_)]              # clip birthweight
  cov[, LBW := as.integer(BWT < 2500)]
  cov[, livebirth := as.integer(BOUT==1)]
  cov[, still := as.integer(BOUT==2)]
  cov[, twin := as.integer(!is.na(TWIN) & TWIN!=1)]
  cov[, subtype := ifelse(is.na(PTB), NA_character_,
                    ifelse(PTB==0,"term",
                      ifelse(is.na(SPONT),"unk",
                        ifelse(SPONT==1,"spont","indic"))))]

  ## ---------- genotyped flag ----------
  ## Audit 2026-07-19 (section L) confirmed all three candidate definitions -- one score,
  ## any score on usable platforms, any score on any platform -- give exactly 14,032 with
  ## zero unmatched IDs, so the single-score probe is not fragile. Kept as-is.
  cov[, genotyped := NA_integer_]
  if(!is.null(SSC)){
    pgs0 <- MOMI_EUR$SBP
    gids <- unique(unlist(lapply(MOMI_COH_ALL, function(c)
      unlist(lapply(MOMI_PLATS, function(pl){
        f <- file.path(SSC, sprintf("%s__%s__%s.sscore", pgs0, c, pl))
        if(file.exists(f)) as.character(fread(f)[[1]]) else NULL})))))
    cov[, genotyped := as.integer(IID %in% gids)]
  }

  A <- Reduce(function(a,b) merge(a,b,by="IID",all.x=TRUE), list(cov,Sd,Dd,Sp,Dp,Sr,Dr))

  ## ---------- analytic-set exclusions (Anagh 2026-07-26) ----------
  ## Harmonised with the MOMI observational/metabolomics papers: EXCLUDE non-singleton
  ## pregnancies (twins/multiples). Stillbirths are RETAINED (birth-weight-based outcomes are
  ## livebirth-restricted downstream). No post-20-week BP requirement. Applying the singleton
  ## restriction here means every downstream deliverable inherits it via analytic_mothers.
  n_first_preg_all <- nrow(A)
  n_nonsingleton   <- sum(A$twin == 1, na.rm = TRUE)
  A <- A[is.na(twin) | twin == 0]                                   # singletons only
  cat(sprintf("analytic-set exclusion: dropped %d non-singleton pregnancies; %d singletons remain\n",
              n_nonsingleton, nrow(A)))
  momi_save_intermediate(
    data.table(n_first_preg_all = n_first_preg_all, n_nonsingleton = n_nonsingleton),
    "flow_pre", P)
  momi_save_intermediate(A, "analytic_mothers", P)

  gk <- if(all(is.na(A$genotyped))) "NA" else as.character(sum(A$genotyped,na.rm=TRUE))
  list(n=nrow(A),
       key=sprintf("mothers=%d genotyped=%s BMIok=%d PTB=%d PE=%d chronHTN=%d twins=%d | visits ante=%d post=%d unk=%d | S_mean(ante)=%d S_mean_pn=%d",
                   nrow(A), gk, sum(!is.na(A$BMI)), sum(A$PTB,na.rm=TRUE),
                   sum(A$PE,na.rm=TRUE), sum(A$CHRON_HTN,na.rm=TRUE), sum(A$twin,na.rm=TRUE),
                   nrow(ANTE), nrow(POST), sum(raw$visit_class=="unknown"),
                   sum(!is.na(A$S_mean)), sum(!is.na(A$S_mean_pn))),
       outputs=c("analytic_mothers.rds","bp_readings_long.rds","bp_icc.rds"))
})
