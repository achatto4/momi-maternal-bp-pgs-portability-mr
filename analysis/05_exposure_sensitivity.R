## -------------------------------------------------------------------------------
## 05_exposure_sensitivity.R
##
## How the Mendelian-randomization results respond to the definition of the antenatal
## blood-pressure exposure. This is a sensitivity analysis of the exposure and not a new
## instrument selection: the polygenic score chosen for each cohort and trait is the same
## one used in 04 and is never reconsidered here.
##
## Four definitions, all built from the same valid antenatal readings of the cleaned sample:
##
##   mean      the primary exposure: the mean of a woman's valid antenatal readings.
##   ge20      the mean of her valid antenatal readings taken at or after 140 gestational
##             days (20 weeks). A woman with no reading in that window leaves this sample.
##   resid     the gestational-age-standardised definition. Within each cohort the
##             reading-level pressure is regressed on a natural cubic spline of gestational
##             age in days with three degrees of freedom; the residual is taken and the
##             cohort's mean reading is added back, which returns the value to millimetres
##             of mercury; each woman's exposure is the mean of her re-centred readings.
##   ge20last  the last antenatal reading at or after 140 gestational days, adjusted for
##             gestational age. Among a woman's valid readings at or after 140 days, the
##             reading at her latest gestational age is taken, averaging any readings tied
##             at that gestational age; within cohort and trait the selected reading is
##             regressed on a natural cubic spline of gestational age with three degrees of
##             freedom, and the exposure is the residual plus the cohort mean of the
##             selected readings, again in millimetres of mercury.
##
## Because the units decide whether a 10 mmHg scaling means anything, the re-centring
## identity is checked rather than assumed for both standardised definitions, and the
## script stops if either fails. Everything else is held at the specification of 04: the
## same outcomes and denominators, the same covariates, the same single prespecified
## two-cell technology recoding, the same Wald ratio with its delta-method standard error,
## the same 10 mmHg scaling, and the same fixed-effect pooling. For every definition and
## outcome the complete-case sample is rebuilt and both stages are re-estimated in it.
## -------------------------------------------------------------------------------
if (!exists("MOMI_ROOT")) MOMI_ROOT <- getwd()
if (!exists("config")) source(file.path(MOMI_ROOT,
  if (file.exists(file.path(MOMI_ROOT, "config.R"))) "config.R" else "config.example.R"))
source(file.path(MOMI_ROOT, "analysis", "00_functions.R"))
suppressMessages({library(splines); library(metafor)})

say  <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
POOL <- TRUE

WK20   <- 140
BASE_DEFS <- c("mean","ge20","resid")
NEWDEF <- "ge20last"
DEFS   <- c(BASE_DEFS, NEWDEF)
DEFLAB <- c(mean  = "overall mean valid antenatal BP (primary)",
            ge20  = "mean valid antenatal BP at 20 gestational weeks or later",
            resid = "gestational-age-standardised BP (per-cohort spline residual, re-centred on the cohort mean, mmHg)",
            ge20last = "last valid antenatal reading at 20 gestational weeks or later (readings tied at that gestational age averaged), gestational-age-standardised within cohort and re-centred on the cohort mean of the selected readings, mmHg")

## The cleaned sample, the outcomes, the covariates and the scores.

S3 <- eligible_sample(EPI, SSC)
AN <- S3$AN[cohort %chin% COH5]
keep_ids <- unique(normid(fread(KEEPF, header=FALSE, colClasses="character")[[1]]))
AN <- AN[IID %chin% keep_ids]

OUT_COLS <- c("SITE_CODE","PARTICIPANT_ID","PREGNANCY_ID","GA_HDLK_NEW","VISITDT","PW_AGE",
              "SBP","DBP","DEL_DATE","DATE_LMP","GAGEBRTH_NEW","BWT_MEASURE_DATE",
              "PTB_NEW","BIRTH_WEIGHT","SGA_10_NEW","BIRTH_OUTCOME")
hdr <- names(fread(EPI, nrows=0)); mi <- setdiff(OUT_COLS, hdr)
if(length(mi)) stop("EPI lacks columns: ", jn(mi))
raw <- fread(EPI, na.strings=MOMI_NA, select=OUT_COLS)[PREGNANCY_ID == 1]
raw[, IID := as.character(PARTICIPANT_ID)]
for(cc in c("GA_HDLK_NEW","VISITDT","PW_AGE","SBP","DBP","GAGEBRTH_NEW","PTB_NEW",
            "BIRTH_WEIGHT","SGA_10_NEW","BIRTH_OUTCOME"))
  set(raw, j=paste0(cc,"_n"), value=num(raw[[cc]]))
PH <- raw[, .(AGE=firstnn(PW_AGE_n), PTB=firstnn(PTB_NEW_n), BWTraw=firstnn(BIRTH_WEIGHT_n),
              SGA=as.integer(firstnn(SGA_10_NEW_n) == 1), BOUT=firstnn(BIRTH_OUTCOME_n)), by=IID]
PH[, BWT := ifelse(BWTraw >= 500 & BWTraw <= 6500, BWTraw, NA_real_)]
PH[, livebirth := as.integer(BOUT == 1)]
AN <- merge(AN, PH[, .(IID, AGE, PTB, BWT, SGA, livebirth)], by="IID", all.x=TRUE)

PL <- if(grepl("\\.rds$", PCSF, ignore.case=TRUE)) as.data.table(readRDS(PCSF)) else
      fread(PCSF, colClasses=list(character="IID"))
PL <- as.data.table(PL); PL[, IID := normid(IID)]
if(length(setdiff(PCS5, names(PL)))) stop("PC table lacks: ", jn(setdiff(PCS5, names(PL))))
AN <- merge(AN, unique(PL[, c("IID", PCS5), with=FALSE], by="IID"), by="IID", all.x=TRUE)
AN[, has_pc5 := Reduce(`&`, lapply(PCS5, function(p) is.finite(get(p))))]

DROP <- if(!is.null(DROPF) && file.exists(DROPF)) fread(DROPF, colClasses="character") else
        data.table(ID=character(0), platform=character(0))
if(nrow(DROP)){ setnames(DROP, 1:2, c("ID","platform")); DROP[, key := paste(ID, platform)] } else DROP[, key := character(0)]
rd_sscore <- function(f, plat){ s <- fread(f); d <- data.table(IID=as.character(s[[1]]), raw=as.numeric(s$SCORE1_AVG))
  if(nrow(DROP)) d <- d[!(paste(IID, plat) %chin% DROP$key)]; d }
getz <- function(pid, coh){
  L <- lapply(PLATS, function(p){ f <- file.path(SSC, sprintf("%s__%s__%s.sscore", pid, coh, p))
    if(!file.exists(f)) return(NULL); s <- rd_sscore(f, p); data.table(IID=s$IID, z=as.numeric(zin(s$raw))) })
  L <- L[!vapply(L, is.null, TRUE)]; if(!length(L)) return(NULL)
  rbindlist(L)[is.finite(z)][, .(z=mean(z)), by=IID][, z := zin(z)][]
}
ZL <- list()
for(coh in COH5) for(sid in unique(c(INSTR[COH5]$SBP, INSTR[COH5]$DBP))){
  z <- getz(sid, coh); if(!is.null(z)) ZL[[paste(coh, sid)]] <- z
}

## The four exposure definitions.

mk_readings <- function(){
  r <- copy(raw)
  r[, SBPv := SBP_n][, DBPv := DBP_n]
  r[is.finite(SBPv) & is.finite(DBPv) & DBPv >= SBPv, `:=`(SBPv = NA_real_, DBPv = NA_real_)]
  r[, GA := ifelse(GA_HDLK_NEW_n >= 0 & GA_HDLK_NEW_n <= 315, GA_HDLK_NEW_n, NA_real_)]
  asd <- function(x){ if(inherits(x, c("Date","IDate"))) return(as.Date(x))
                      suppressWarnings(as.Date(as.character(x), format="%Y-%m-%d")) }
  r[, VISd := asd(VISITDT)][, DELd := asd(DEL_DATE)][, LMPd := asd(DATE_LMP)][, BWTd := asd(BWT_MEASURE_DATE)]
  r[, DELeff := DELd]
  r[is.na(DELeff) & !is.na(LMPd) & is.finite(GAGEBRTH_NEW_n), DELeff := LMPd + GAGEBRTH_NEW_n]
  r[is.na(DELeff) & !is.na(BWTd), DELeff := BWTd]
  r[, DELeff := { v <- DELeff[!is.na(DELeff)]; if(length(v)) v[1] else as.Date(NA) }, by=IID]
  r[, visit_class := fifelse(is.na(VISd) | is.na(DELeff), "unknown",
                      fifelse(VISd > DELeff, "postnatal", "antenatal"))]
  r[, GA_lmp := as.numeric(VISd - LMPd)]
  r[visit_class == "unknown" & is.finite(GA_lmp) & GA_lmp >= 0 & GA_lmp <= 294, visit_class := "antenatal"]
  r[visit_class == "antenatal"]
}
ANTE <- mk_readings()
ANTE <- merge(ANTE[, .(IID, GA, SBPv, DBPv)], AN[, .(IID, cohort)], by="IID")
LONG <- rbindlist(list(ANTE[is.finite(SBPv), .(IID, cohort, GA, trait="SBP", bp=SBPv)],
                       ANTE[is.finite(DBPv), .(IID, cohort, GA, trait="DBP", bp=DBPv)]))

wm <- function(v, keep){ sel <- keep & !is.na(v); sel[is.na(sel)] <- FALSE
                         if(!any(sel)) NA_real_ else mean(v[sel]) }
EXPO <- list(); SELINFO <- list()
for(tr in TRAITS){
  L <- LONG[trait == tr]
  a <- L[, .(mean = mean(bp), ge20 = wm(bp, GA >= WK20)), by=IID]

  KN <- list()
  rr <- L[is.finite(bp) & is.finite(GA), {
    m <- tryCatch(lm(bp ~ ns(GA, 3)), error=function(e) NULL)
    if(is.null(m)) .(IID=IID, r=rep(NA_real_, .N))
    else { b <- ns(GA, 3)
           KN[[paste(tr, .BY$cohort)]] <<- data.table(bp_trait=tr, cohort=.BY$cohort,
                 n_readings=.N, df=3L,
                 interior_knots=jn(sprintf("%.1f", attr(b, "knots"))),
                 boundary_knots=jn(sprintf("%.1f", attr(b, "Boundary.knots"))),
                 cohort_mean_reading=mean(bp), cohort_sd_reading=sd(bp))
           .(IID=IID, r=resid(m) + mean(bp)) } }, by=cohort]
  b <- rr[, .(resid = mean(r, na.rm=TRUE)), by=IID]

  SEL <- L[is.finite(bp) & is.finite(GA) & GA >= WK20]
  if(nrow(SEL)){
    SEL[, gmax := max(GA), by=IID]
    SEL <- SEL[GA == gmax, .(bp = mean(bp), n_readings_at_last = .N), by=.(IID, cohort, GA)]
  } else SEL <- data.table(IID=character(0), cohort=character(0), GA=numeric(0), bp=numeric(0),
                           n_readings_at_last=integer(0))
  SELINFO[[tr]] <- copy(SEL)[, bp_trait := tr]

  KN4 <- list()
  if(nrow(SEL)){
    r4 <- SEL[, { m <- tryCatch(lm(bp ~ ns(GA, 3)), error=function(e) NULL)
      if(is.null(m)) .(IID=IID, r=rep(NA_real_, .N))
      else { bs <- ns(GA, 3)
             KN4[[paste(tr, .BY$cohort)]] <<- data.table(bp_trait=tr, cohort=.BY$cohort,
                   exposure_definition=NEWDEF, n_selected_readings=.N, df=3L,
                   n_distinct_GA=uniqueN(GA),
                   interior_knots=jn(sprintf("%.1f", attr(bs, "knots"))),
                   boundary_knots=jn(sprintf("%.1f", attr(bs, "Boundary.knots"))),
                   cohort_mean_selected_reading=mean(bp), cohort_sd_selected_reading=sd(bp),
                   min_GA=min(GA), median_GA=as.numeric(median(GA)), max_GA=max(GA))
             .(IID=IID, r=resid(m) + mean(bp)) } }, by=cohort]
    d4 <- r4[, .(ge20last = r[1]), by=IID]
    if(nrow(d4) != uniqueN(SEL$IID))
      stop("the last-reading definition produced ", nrow(d4), " values for ", uniqueN(SEL$IID), " women")
  } else d4 <- data.table(IID=character(0), ge20last=numeric(0))
  assign(paste0("KNOTS4_", tr), if(length(KN4)) rbindlist(KN4) else
         data.table(bp_trait=character(0), cohort=character(0)))

  EXPO[[tr]] <- Reduce(function(x, y) merge(x, y, by="IID", all=TRUE), list(a, b, d4))
  assign(paste0("KNOTS_", tr), rbindlist(KN))
}
KNOTS <- rbindlist(list(KNOTS_SBP, KNOTS_DBP))
W_full(KNOTS, "ga_spline_knots.tsv")

recentred <- rbindlist(lapply(TRAITS, function(tr){
  L <- LONG[trait == tr & is.finite(bp) & is.finite(GA)]
  L[, { m <- tryCatch(lm(bp ~ ns(GA, 3)), error=function(e) NULL)
        .(bp_trait = tr, mean_raw_reading = mean(bp),
          mean_recentred_reading = if(is.null(m)) NA_real_ else mean(resid(m) + mean(bp)),
          n_readings = .N) }, by=cohort]
}))
recentred[, construction_diff_mmHg := abs(mean_raw_reading - mean_recentred_reading)]
UCHK <- rbindlist(lapply(TRAITS, function(tr){
  E <- merge(EXPO[[tr]][, .(IID, resid)], AN[, .(IID, cohort)], by="IID")
  E[is.finite(resid), .(bp_trait = tr, n_women = .N,
                        mean_woman_exposure = mean(resid), sd_woman_exposure = sd(resid)), by=cohort]
}))
UCHK <- merge(UCHK, recentred, by=c("cohort","bp_trait"))
UCHK[, `:=`(mean_shift_mmHg = abs(mean_raw_reading - mean_woman_exposure),
            looks_like_a_z_score = abs(mean_woman_exposure) < 5 & sd_woman_exposure < 1.5)]
W_full(UCHK, "ga_standardised_units_check.tsv")
if(any(!is.finite(UCHK$construction_diff_mmHg) | UCHK$construction_diff_mmHg > 1e-8))
  stop("the gestational-age-standardised exposure was not re-centred onto the reading scale: the ",
       "reading-level mean differs from the raw mean by ",
       sprintf("%.3g", max(UCHK$construction_diff_mmHg, na.rm=TRUE)), " mmHg")
if(any(UCHK$looks_like_a_z_score))
  stop("the gestational-age-standardised exposure looks like a standardised score, not millimetres of ",
       "mercury (a cohort mean near zero with a spread below 1.5). A 10-mmHg scaling would not be ",
       "interpretable; stopping rather than forcing a result.")
if(any(UCHK$mean_shift_mmHg > 8))
  stop("the gestational-age-standardised exposure sits ", sprintf("%.1f", max(UCHK$mean_shift_mmHg)),
       " mmHg from its cohort's mean reading, which is more than the weighting of women by their ",
       "number of readings can explain; stopping rather than forcing a result.")
say(sprintf("units check passed: re-centring exact to %.1e mmHg; woman-level exposure spread %.1f-%.1f mmHg",
            max(UCHK$construction_diff_mmHg), min(UCHK$sd_woman_exposure), max(UCHK$sd_woman_exposure)))

KNOTS4 <- rbindlist(list(KNOTS4_SBP, KNOTS4_DBP), fill=TRUE)
W_full(KNOTS4, "last_reading_spline_knots.tsv")
SELALL <- rbindlist(SELINFO)
SELSUM <- SELALL[, .(n_women_with_a_selected_reading = .N,
                     n_women_with_tied_readings_averaged = sum(n_readings_at_last > 1L),
                     max_readings_tied_at_last = max(n_readings_at_last),
                     mean_selected_reading = mean(bp), sd_selected_reading = sd(bp),
                     min_GA_days = min(GA), median_GA_days = as.numeric(median(GA)), max_GA_days = max(GA)),
                 by=.(cohort, bp_trait)]
SELSUM <- merge(SELSUM, AN[, .(n_in_cleaned_sample = .N), by=cohort], by="cohort", all.x=TRUE)
SELSUM[, n_women_without_a_reading_in_window := n_in_cleaned_sample - n_women_with_a_selected_reading]
W_full(SELSUM, "last_reading_selection_summary.tsv")

U4 <- rbindlist(lapply(TRAITS, function(tr){
  S <- SELALL[bp_trait == tr]
  if(!nrow(S)) return(NULL)
  S[, { m <- tryCatch(lm(bp ~ ns(GA, 3)), error=function(e) NULL)
        .(bp_trait = tr, n_selected = .N, mean_selected_reading = mean(bp),
          mean_recentred = if(is.null(m)) NA_real_ else mean(resid(m) + mean(bp)),
          sd_recentred   = if(is.null(m)) NA_real_ else sd(resid(m) + mean(bp))) }, by=cohort] }))
U4[, construction_diff_mmHg := abs(mean_selected_reading - mean_recentred)]
U4[, looks_like_a_z_score := abs(mean_recentred) < 5 & sd_recentred < 1.5]
U4[, exposure_definition := NEWDEF]
W_full(U4, "last_reading_units_check.tsv")
if(!nrow(U4) || any(!is.finite(U4$construction_diff_mmHg) | U4$construction_diff_mmHg > 1e-8))
  stop("the last-reading exposure was not re-centred onto the reading scale: the mean of the re-centred ",
       "selected readings differs from their raw mean by ",
       sprintf("%.3g", max(U4$construction_diff_mmHg, na.rm=TRUE)), " mmHg")
if(any(U4$looks_like_a_z_score))
  stop("the last-reading exposure looks like a standardised score rather than millimetres of mercury; a ",
       "10-mmHg scaling would not be interpretable, so it is not run")
say(sprintf("last-reading definition: %d women selected, %d had readings tied at the last gestational age; ",
            sum(SELSUM$n_women_with_a_selected_reading), sum(SELSUM$n_women_with_tied_readings_averaged)),
    sprintf("re-centring exact to %.1e mmHg", max(U4$construction_diff_mmHg)))

SM <- merge(EXPO$SBP[, .(IID, mean)], S3$ALLW[, .(IID, S_mean)], by="IID")
DM <- merge(EXPO$DBP[, .(IID, mean)], S3$ALLW[, .(IID, D_mean)], by="IID")
SM <- SM[IID %chin% AN$IID]; DM <- DM[IID %chin% AN$IID]
bad_mean <- sum(!is.na(SM$mean) & SM$mean != SM$S_mean) + sum(!is.na(DM$mean) & DM$mean != DM$D_mean)
if(bad_mean) stop("the rebuilt mean exposure differs from the primary one in ", bad_mean, " women")
say("the rebuilt mean exposure equals the primary S_mean / D_mean in every woman")

## The cells, rebuilt for every definition.

RES <- list(); CNT <- list()
for(dfn in DEFS) for(coh in COH5) for(tr in TRAITS){
  sid <- INSTR[coh][[tr]]
  zt  <- ZL[[paste(coh, sid)]]
  D <- merge(AN[cohort == coh], EXPO[[tr]][, c("IID", dfn), with=FALSE], by="IID", all.x=TRUE)
  setnames(D, dfn, "bp")
  D <- merge(D, if(is.null(zt)) data.table(IID=character(0), z=numeric(0)) else zt, by="IID", all.x=TRUE)
  for(oc in OUTSPEC$outcome){
    osp <- OUTSPEC[oc]; binary <- osp$type == "binary"
    lb <- D$livebirth == 1; lb[is.na(lb)] <- FALSE
    y <- switch(oc, PTB = as.integer(D$PTB == 1), LBW = as.integer(D$BWT < 2500),
                    SGA = as.integer(D$SGA),      BWT = as.numeric(D$BWT))
    den <- switch(oc, PTB = !is.na(D$PTB), LBW = lb & !is.na(D$BWT),
                      SGA = lb & !is.na(D$SGA), BWT = lb & !is.na(D$BWT))
    den[is.na(den)] <- FALSE
    ok <- den & is.finite(D$z) & is.finite(D$bp) & is.finite(D$AGE) &
          D$technology %chin% TECH_ALL & D$has_pc5 & is.finite(y)
    ok[is.na(ok)] <- FALSE
    idx <- which(ok)
    Dc <- as.data.frame(D[idx, c("technology","AGE","z","bp", PCS5), with=FALSE]); Dc$y <- y[idx]
    rownames(Dc) <- paste0("w", idx)

    tech_exc <- paste(coh, tr, oc, sep="|") %chin% TECH_EXCEPT$key
    tech_raw <- as.character(Dc$technology)
    tech_use <- if(tech_exc) ifelse(tech_raw %chin% TECH_COLLAPSE_FROM, TECH_COLLAPSE_TO, tech_raw) else tech_raw
    lv_all <- if(tech_exc) TECH_EXCEPT_LEVELS else TECH_ALL
    tabt <- table(factor(tech_use, levels=lv_all)); lv <- lv_all[tabt > 0]
    ref <- if(TECH_REF %in% lv) TECH_REF else names(tabt)[which.max(as.integer(tabt))]
    lv  <- c(ref, setdiff(lv, ref)); Dc$tech <- factor(tech_use, levels=lv)
    tech_in <- length(lv) > 1
    rhs  <- paste(c("z","AGE", if(tech_in) "tech", PCS5), collapse=" + ")
    n    <- nrow(Dc)
    npar <- 2L + 1L + (if(tech_in) length(lv) - 1L else 0L) + length(PCS5)
    nev  <- if(binary) sum(Dc$y == 1L) else NA_integer_

    CNT[[length(CNT)+1]] <- data.table(
      exposure_definition=dfn, exposure_label=unname(DEFLAB[dfn]), cohort=coh,
      cohort_display=unname(DISPLAY[coh]), ancestry_group=unname(GROUP[coh]), bp_trait=tr, outcome=oc,
      pgs_id=sid, N=n, n_events=nev,
      n_with_exposure=sum(is.finite(D$bp)), n_in_cohort=nrow(D),
      technology_levels_modelled=length(lv), technology_coding=if(tech_exc) "two-level exception" else "three-level")

    stops <- character(0); fitted_both <- FALSE
    fs_b <- fs_se <- fs_F <- rf_b <- rf_se <- NA_real_
    if(n < npar + 2L) stops <- c(stops, "too few complete cases to fit the prespecified model")
    if(binary && !length(stops) && (nev == 0L || nev == n)) stops <- c(stops, "no variation in the outcome")
    zero_ev_level <- FALSE
    if(binary && tech_in){
      tb <- as.data.table(Dc)[, .(n=.N, ev=sum(y == 1L)), by=tech]
      zero_ev_level <- any(tb$ev == 0L | tb$ev == tb$n)
    }
    if(!length(stops)){
      fitted_both <- TRUE
      m_fs <- lm(as.formula(paste("bp ~", rhs)), data=Dc)
      gg <- if(binary) fit_glm(as.formula(paste("y ~", rhs)), Dc) else
                       list(m=lm(as.formula(paste("y ~", rhs)), data=Dc), warn=character(0))
      m_rf <- gg$m
      mf1 <- model.frame(m_fs); mf2 <- model.frame(m_rf)
      X1 <- model.matrix(m_fs);  X2 <- model.matrix(m_rf)
      same <- identical(rownames(mf1), rownames(mf2)) && identical(colnames(X1), colnames(X2)) &&
              identical(dim(X1), dim(X2)) && max(abs(X1 - X2)) == 0 &&
              identical(attr(terms(m_fs), "term.labels"), attr(terms(m_rf), "term.labels"))
      s1 <- summary(m_fs)$coefficients
      if(any(is.na(coef(m_fs))) || !("z" %in% rownames(s1))) stops <- c(stops, "first stage is rank-deficient")
      else { fs_b <- s1["z",1]; fs_se <- s1["z",2]; fs_F <- (fs_b/fs_se)^2 }
      s2 <- tryCatch(summary(m_rf)$coefficients, error=function(e) NULL)
      if(!is.null(s2) && "z" %in% rownames(s2)){ rf_b <- s2["z",1]; rf_se <- s2["z",2] }
      if(binary){
        st <- logit_status(m_rf, gg$warn)
        if(!st$converged) stops <- c(stops, "logistic model did not converge")
        if(st$separation || zero_ev_level) stops <- c(stops, "logistic separation")
      } else if(any(is.na(coef(m_rf))) || !all(is.finite(coef(m_rf))))
        stops <- c(stops, "reduced-form OLS is rank-deficient")
      if(is.finite(fs_b) && !(fs_b > 0))
        stops <- c(stops, "first-stage PGS coefficient is not positive (orientation)")
      if(!same) stops <- c(stops, "the two stages were not fitted in the same women or the same specification")
    }
    status <- if(length(stops)) paste0("stopped: ", jn(unique(stops))) else "ok"

    th <- th_se <- e10 <- se10 <- or10 <- or_lo <- or_hi <- g10 <- g_lo <- g_hi <- wp <- NA_real_
    if(status == "ok" && is.finite(fs_b) && fs_b > 0 && is.finite(fs_se) && is.finite(rf_b) && is.finite(rf_se)){
      th <- rf_b / fs_b
      th_se <- sqrt(rf_se^2 / fs_b^2 + rf_b^2 * fs_se^2 / fs_b^4)
      e10 <- SCALE * th; se10 <- SCALE * th_se
      wp <- 2 * pnorm(-abs(th / th_se))
      if(binary){ or10 <- exp(e10); or_lo <- exp(e10 - ZCRIT*se10); or_hi <- exp(e10 + ZCRIT*se10) }
      else      { g10 <- e10; g_lo <- e10 - ZCRIT*se10; g_hi <- e10 + ZCRIT*se10 }
    }
    RES[[length(RES)+1]] <- data.table(
      exposure_definition=dfn, exposure_label=unname(DEFLAB[dfn]),
      cohort=coh, cohort_display=unname(DISPLAY[coh]), ancestry_group=unname(GROUP[coh]),
      bp_trait=tr, pgs_id=sid, outcome=oc, outcome_label=osp$label, outcome_type=osp$type,
      N=n, n_events=nev,
      fs_beta_mmHg_per_SD=fs_b, fs_se=fs_se, fs_F=fs_F, weak_instrument=isTRUE(fs_F < WEAK_F),
      rf_beta_per_SD=rf_b, rf_se=rf_se,
      wald_estimate_per_10mmHg=e10, wald_se_per_10mmHg=se10,
      or_per_10mmHg=or10, or_lo95=or_lo, or_hi95=or_hi,
      grams_per_10mmHg=g10, grams_lo95=g_lo, grams_hi95=g_hi,
      wald_p=wp, fitted_both_stages=fitted_both, status=status)
  }
}
RES <- rbindlist(RES); CNT <- rbindlist(CNT)

W_full(RES, "exposure_four_definitions_cohort.tsv")
W_full(CNT, "exposure_four_definitions_sample_counts.tsv")

## Pooling, the same fixed-effect model as in 04.

if(POOL && requireNamespace("metafor", quietly=TRUE)){
  suppressMessages(library(metafor))
  PO <- rbindlist(lapply(DEFS, function(dfn) rbindlist(lapply(TRAITS, function(tr)
    rbindlist(lapply(OUTSPEC$outcome, function(oc){
      D <- RES[exposure_definition == dfn & bp_trait == tr & outcome == oc & status == "ok" &
               is.finite(wald_estimate_per_10mmHg) & is.finite(wald_se_per_10mmHg) & wald_se_per_10mmHg > 0]
      bin <- oc != "BWT"
      if(nrow(D) < 2L)
        return(data.table(exposure_definition=dfn, bp_trait=tr, outcome=oc,
                          outcome_type=if(bin) "binary" else "continuous",
                          k_cohorts=nrow(D), N_total=sum(D$N),
                          events_total=if(bin) sum(as.numeric(D$n_events)) else NA_real_,
                          pooled_estimate=NA_real_, pooled_se=NA_real_, ci_lo=NA_real_, ci_hi=NA_real_,
                          or_per_10mmHg=NA_real_, or_lo95=NA_real_, or_hi95=NA_real_,
                          grams_per_10mmHg=NA_real_, grams_lo95=NA_real_, grams_hi95=NA_real_,
                          z=NA_real_, p=NA_real_, tau2=NA_real_, I2=NA_real_, Q=NA_real_,
                          Q_df=NA_integer_, Q_p=NA_real_,
                          status=sprintf("not pooled: only %d estimable cohort(s)", nrow(D))))
      m <- rma(yi=D$wald_estimate_per_10mmHg, sei=D$wald_se_per_10mmHg, method="FE", test="z")
      data.table(exposure_definition=dfn, bp_trait=tr, outcome=oc,
                 outcome_type=if(bin) "binary" else "continuous",
                 k_cohorts=m$k, N_total=sum(D$N),
                 events_total=if(bin) sum(as.numeric(D$n_events)) else NA_real_,
                 pooled_estimate=as.numeric(m$beta), pooled_se=m$se, ci_lo=m$ci.lb, ci_hi=m$ci.ub,
                 or_per_10mmHg=if(bin) exp(as.numeric(m$beta)) else NA_real_,
                 or_lo95=if(bin) exp(m$ci.lb) else NA_real_, or_hi95=if(bin) exp(m$ci.ub) else NA_real_,
                 grams_per_10mmHg=if(bin) NA_real_ else as.numeric(m$beta),
                 grams_lo95=if(bin) NA_real_ else m$ci.lb, grams_hi95=if(bin) NA_real_ else m$ci.ub,
                 z=m$zval, p=m$pval, tau2=m$tau2, I2=m$I2, Q=m$QE, Q_df=m$k - 1L, Q_p=m$QEp,
                 status="ok")
    }))))))
  W_full(PO, "exposure_four_definitions_pooled.tsv")
  say(sprintf("pooled: %d rows (%d definitions x %d traits x %d outcomes)",
              nrow(PO), length(DEFS), length(TRAITS), nrow(OUTSPEC)))
  cat("\n==== pooled, by exposure definition ====\n")
  print(PO[, .(exposure_definition, bp_trait, outcome, k=k_cohorts, N=N_total,
               est=signif(fifelse(outcome_type == "binary", or_per_10mmHg, grams_per_10mmHg), 3),
               lo=signif(fifelse(outcome_type == "binary", or_lo95, grams_lo95), 3),
               hi=signif(fifelse(outcome_type == "binary", or_hi95, grams_hi95), 3),
               p=signif(p, 2), I2=round(I2, 1))])
} else say("metafor unavailable or pooling suppressed; cohort results written, pooled table skipped")

W_full(data.table(
  item  = c("definitions","wk20_days","spline","units_check","R","metafor","run_at"),
  value = c(jn(sprintf("%s = %s", DEFS, DEFLAB[DEFS])), WK20,
            "per cohort, lm(bp ~ ns(GA_days, 3)); residual plus the cohort mean reading",
            "re-centring exact at reading level (1e-8 mmHg); the exposure is in mmHg, not a z score",
            R.version.string, as.character(packageVersion("metafor")),
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))),
  "exposure_sensitivity_run_info.tsv")
cat(sprintf("05: %d cells across %d exposure definitions\n", nrow(RES), length(DEFS)))
