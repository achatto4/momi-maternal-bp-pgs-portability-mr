suppressMessages(library(data.table))
eligible_sample <- function(EPI, SSC){
  MOMI_NA <- c("-88","-77","-99","NA","na",".","")
  SITE <- c("1"="GAPPS-Zambia","2"="GAPPS-Bangladesh","3"="AMANHI-Pakistan","4"="AMANHI-Bangladesh","5"="AMANHI-Pemba","6"="THSTI-India")
  COH5 <- c("AMANHI-Bangladesh","AMANHI-Pakistan","GAPPS-Bangladesh","AMANHI-Pemba","GAPPS-Zambia")
  EPI_COLS <- c("SITE_CODE","PARTICIPANT_ID","PREGNANCY_ID","SBP","DBP","GA_HDLK_NEW","VISITDT","PW_AGE",
                "SINGLE_TWIN","DEL_DATE","DATE_LMP","GAGEBRTH_NEW","BWT_MEASURE_DATE")
  num <- function(x) suppressWarnings(as.numeric(x))
  firstnn <- function(v){ w <- which(!is.na(v)); if(!length(w)) v[NA_integer_][1] else v[w[1]] }
  asdate <- function(x){ if(inherits(x, c("Date","IDate"))) return(as.Date(x)); suppressWarnings(as.Date(as.character(x), format="%Y-%m-%d")) }
  raw <- fread(EPI, na.strings=MOMI_NA, select=EPI_COLS)[PREGNANCY_ID == 1]
  getcol <- function(nm) if(nm %in% names(raw)) num(raw[[nm]]) else rep(NA_real_, nrow(raw))
  raw[, IID := as.character(PARTICIPANT_ID)]
  raw[, SBP := getcol("SBP")]; raw[, DBP := getcol("DBP")]
  raw[, pair_invalid0 := is.finite(SBP) & is.finite(DBP) & DBP >= SBP]
  raw[pair_invalid0 == TRUE, `:=`(SBP = NA_real_, DBP = NA_real_)]
  garead <- getcol("GA_HDLK_NEW"); raw[, GA := ifelse(garead >= 0 & garead <= 315, garead, NA_real_)]
  raw[, VIS := getcol("VISITDT")]; raw[, AGEv := getcol("PW_AGE")]; raw[, GAdv := getcol("GAGEBRTH_NEW")]; raw[, TWINv := getcol("SINGLE_TWIN")]
  raw[, cohort := SITE[as.character(SITE_CODE)]]
  setorder(raw, IID, VIS, GA)
  raw[, VISd := asdate(VISITDT)]; raw[, DELd := asdate(DEL_DATE)]; raw[, LMPd := asdate(DATE_LMP)]; raw[, BWTd := asdate(BWT_MEASURE_DATE)]
  raw[, DELeff := DELd]
  raw[is.na(DELeff) & !is.na(LMPd) & is.finite(GAdv), DELeff := LMPd + GAdv]
  raw[is.na(DELeff) & !is.na(BWTd), DELeff := BWTd]
  raw[, DELeff := { v <- DELeff[!is.na(DELeff)]; if(length(v)) v[1] else as.Date(NA) }, by=IID]
  raw[, visit_class := fifelse(is.na(VISd) | is.na(DELeff), "unknown", fifelse(VISd > DELeff, "postnatal", "antenatal"))]
  raw[, GA_lmp := as.numeric(VISd - LMPd)]
  raw[visit_class == "unknown" & is.finite(GA_lmp) & GA_lmp >= 0 & GA_lmp <= 294, visit_class := "antenatal"]
  mnv <- function(v) if(all(is.na(v))) NA_real_ else mean(v, na.rm=TRUE)
  ANTE <- raw[visit_class == "antenatal"]
  Sd <- ANTE[, .(S_mean=mnv(SBP)), by=IID]; Dd <- ANTE[, .(D_mean=mnv(DBP)), by=IID]
  cv <- raw[, .(cohort=cohort[1], TWIN=firstnn(TWINv)), by=IID]; cv[, twin := as.integer(!is.na(TWIN) & TWIN != 1)]
  ALLW <- Reduce(function(a, b) merge(a, b, by="IID", all.x=TRUE), list(cv, Sd, Dd))
  ALLW[, hasBP := is.finite(S_mean) | is.finite(D_mean)]
  probe <- function(pl) unique(unlist(lapply(COH5, function(c){ f <- file.path(SSC, sprintf("PGS004603__%s__%s.sscore", c, pl))
    if(file.exists(f)) as.character(fread(f)[[1]]) else character(0) })))
  gsa_ids <- probe("gsa"); lp_ids <- probe("lpwgs_dosage")
  ALLW[, technology := fifelse(IID %chin% gsa_ids & IID %chin% lp_ids, "both",
                       fifelse(IID %chin% gsa_ids, "GSA only", fifelse(IID %chin% lp_ids, "low-pass WGS only", "none")))]
  AN <- ALLW[twin == 0 & technology != "none" & hasBP == TRUE & cohort %in% COH5, .(IID, cohort, technology)]
  GEN <- ALLW[technology != "none" & cohort %in% COH5, .(IID, cohort, technology)]
  GEN[, f03_eligible := IID %chin% AN$IID]
  list(ALLW=ALLW, AN=AN, GEN=GEN) }
