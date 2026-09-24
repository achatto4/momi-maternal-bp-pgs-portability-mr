## -------------------------------------------------------------------------------
## 00_functions.R -- shared definitions: cohort labels, the prespecified polygenic-score
## instruments, outcome definitions, covariate coding, and the helper functions used by
## every later script.  Sourced by 01-06; it fits nothing on its own.
## -------------------------------------------------------------------------------
suppressMessages(library(data.table))

## Directories, all taken from the configuration file.
EPI     <- config$epi_file
SSC     <- config$sscore_dir
PCSF    <- config$pc_file
F4B     <- config$genotype_qc_dir
DERIVED <- config$derived_dir
RESULTS <- config$results_dir
FIGURES <- config$figures_dir
TABLES  <- config$tables_dir
for (d in c(DERIVED, RESULTS, FIGURES, TABLES))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

## Identifier lists written by 01 and read by 02, 03, 04 and 05.
KEEPF <- file.path(DERIVED, "cleaned_sample_ids.txt")
DROPF <- file.path(DERIVED, "dropped_genotype_records_ids.tsv")

## Participant identifiers are normalised the same way everywhere, so that identifiers
## read from genotype files, relatedness tables and the phenotype extract compare equal.

suppressMessages(library(data.table))
normid <- function(x) { x <- as.character(x); if(!length(x)) return(character(0))
  as.character(fread(text=paste(c("IID", x), collapse="\n"), sep="\t", header=TRUE)$IID) }
rd_psam <- function(f){ x <- fread(f, colClasses="character"); setnames(x, 1, sub("^#", "", names(x)[1]))
  if(!"IID" %in% names(x)) setnames(x, 1, "IID"); x[, nid := normid(IID)]; x }
wr_keep <- function(psam, nids, f){ cols <- intersect(c("FID", "IID"), names(psam)); y <- psam[nid %chin% nids, ..cols]
  setnames(y, 1, paste0("#", names(y)[1])); fwrite(y, f, sep="\t"); nrow(y) }
rd_tab <- function(f, idcols){ idcols <- intersect(idcols, names(fread(f, nrows=0))); x <- fread(f, colClasses=list(character=idcols))
  setnames(x, 1, sub("^#", "", names(x)[1]))
  for(c in sub("^#", "", idcols)) if(c %in% names(x)) set(x, j=c, value=normid(x[[c]])); x }
rd_kin0 <- function(f){ h <- names(fread(f, nrows=0)); idc <- intersect(h, c("#FID1", "FID1", "#IID1", "IID1", "FID2", "IID2"))
  x <- fread(f, colClasses=list(character=idc)); names(x) <- sub("^#", "", names(x))
  if(!nrow(x)) return(data.table(ID1=character(0), ID2=character(0), NSNP=numeric(0), IBS0_FRAC=numeric(0), KINSHIP=numeric(0)))
  x[, IBS0 := as.numeric(IBS0)]; if(any(x$IBS0 > 1, na.rm=TRUE)) x[, IBS0 := IBS0 / as.numeric(NSNP)]
  x[, .(ID1=normid(IID1), ID2=normid(IID2), NSNP=as.numeric(NSNP), IBS0_FRAC=IBS0, KINSHIP=as.numeric(KINSHIP))] }

## The analysis sample: first recorded pregnancy, not a known multiple pregnancy, at
## least one valid antenatal blood-pressure reading, and at least one genotype record.
## Returns the full per-woman table (ALLW), the analysis sample (AN) and the genotyped
## sample (GEN), together with the genotyping technology available for each woman.

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

## Prespecified study constants: cohort names and display labels, the polygenic score
## chosen in advance for each cohort and trait, the outcome definitions and their
## denominators, the genotyping-technology covariate and its single prespecified
## two-cell recoding, and the estimation settings.

MOMI_NA <- c("-88","-77","-99","NA","na",".","")
COH5    <- config$cohorts
DISPLAY <- c("AMANHI-Bangladesh"="AMANHI-Sylhet, Bangladesh", "AMANHI-Pakistan"="AMANHI-Karachi, Pakistan",
             "GAPPS-Bangladesh"="PreSSMat-Matlab, Bangladesh", "AMANHI-Pemba"="AMANHI-Pemba, Tanzania",
             "GAPPS-Zambia"="ZAPPS-Lusaka, Zambia")
SITE_SHORT <- c("AMANHI-Bangladesh"="Sylhet", "AMANHI-Pakistan"="Karachi", "GAPPS-Bangladesh"="Matlab",
                "AMANHI-Pemba"="Pemba", "GAPPS-Zambia"="Lusaka")
GROUP   <- c("AMANHI-Bangladesh"="South Asian","AMANHI-Pakistan"="South Asian","GAPPS-Bangladesh"="South Asian",
             "AMANHI-Pemba"="African","GAPPS-Zambia"="African")
PLATS   <- c("gsa","lpwgs_dosage")

INSTR <- data.table(cohort = c("AMANHI-Bangladesh","AMANHI-Pakistan","GAPPS-Bangladesh","AMANHI-Pemba","GAPPS-Zambia"),
                    SBP    = c("PGS004830","PGS004830","PGS004830","PGS004603","PGS004603"),
                    DBP    = c("PGS004758","PGS004758","PGS004758","PGS004604","PGS004604"))
setkey(INSTR, cohort)
if(length(setdiff(COH5, INSTR$cohort)))
  stop("no prespecified instrument for: ", paste(setdiff(COH5, INSTR$cohort), collapse=", "))
BPCOL  <- c(SBP="S_mean", DBP="D_mean")
TRAITS <- c("SBP","DBP")

OUTSPEC <- data.table(
  outcome = c("PTB","LBW","SGA","BWT"),
  type    = c("binary","binary","binary","continuous"),
  label   = c("Preterm birth", "Low birth weight (<2,500 g)",
              "Small for gestational age (<10th centile)", "Birth weight, g"),
  denom   = c("women with a non-missing preterm-birth indicator",
              "live births with a birth weight of 500-6,500 g",
              "live births with a small-for-gestational-age classification",
              "live births with a birth weight of 500-6,500 g"))
setkey(OUTSPEC, outcome)

TECH_ALL <- c("GSA only","low-pass WGS only","both")
TECH_REF <- "low-pass WGS only"

TECH_EXCEPT <- data.table(cohort   = c("GAPPS-Zambia", "GAPPS-Zambia"),
                          bp_trait = c("SBP",          "DBP"),
                          outcome  = c("SGA",          "SGA"))
TECH_EXCEPT[, key := paste(cohort, bp_trait, outcome, sep="|")]
TECH_COLLAPSE_FROM <- c("GSA only", "both")
TECH_COLLAPSE_TO   <- "any GSA"
TECH_EXCEPT_LEVELS <- c(TECH_REF, TECH_COLLAPSE_TO)
TECH_CODE_STD <- sprintf("three-level (%s)", paste(TECH_ALL, collapse=" / "))
TECH_CODE_EXC <- sprintf("two-level prespecified exception (%s; %s = %s)",
                         paste(TECH_EXCEPT_LEVELS, collapse=" / "), TECH_COLLAPSE_TO,
                         paste(TECH_COLLAPSE_FROM, collapse=" + "))

PCS5     <- paste0("PC", 1:5)
SCALE    <- 10
ZCRIT    <- 1.96
WEAK_F   <- 10
GLM_EPS  <- 1e-14
GLM_MAXIT<- 100L
SEP_P    <- 1e-8
SEP_BSD  <- 10
SEP_SESD <- 5
OUT_COLS <- c("SITE_CODE","PARTICIPANT_ID","PREGNANCY_ID","GA_HDLK_NEW","VISITDT","PW_AGE",
              "PTB_NEW","BIRTH_WEIGHT","SGA_10_NEW","BIRTH_OUTCOME")

num     <- function(x) suppressWarnings(as.numeric(x))
firstnn <- function(v){ w <- which(!is.na(v)); if(!length(w)) v[NA_integer_][1] else v[w[1]] }
md5     <- function(f) tryCatch(unname(tools::md5sum(f)), error=function(e) NA_character_)
W_      <- function(dt, name) fwrite(dt, file.path(RESULTS, name), sep="\t")

fmt_rt <- function(x){
  v <- rep(NA_character_, length(x))
  nf <- !is.na(x) & !is.finite(x); v[nf] <- as.character(x[nf])
  for(k in which(is.finite(x))){
    s <- NULL
    for(d in 15:17){ s <- sprintf(paste0("%.", d, "g"), x[k]); if(as.numeric(s) == x[k]) break }
    v[k] <- s
  }
  v
}
W_full  <- function(dt, name){ d <- copy(dt)
  for(cc in names(d)) if(is.double(d[[cc]])) set(d, j=cc, value=fmt_rt(d[[cc]]))
  fwrite(d, file.path(RESULTS, name), sep="\t") }
jn      <- function(v) paste(v, collapse="; ")

## Polygenic score construction.  Within cohort and platform the raw score is
## standardised; a woman scored on both platforms takes the mean of her two standardised
## values; the merged score is then re-standardised within cohort.  Genotype records
## excluded by 01 are removed before any standardisation.

zin <- function(x){ s <- sd(x, na.rm=TRUE); if(is.na(s) || s == 0) x*NA else (x - mean(x, na.rm=TRUE))/s }
DROP <- if(!is.null(DROPF) && file.exists(DROPF)) fread(DROPF, colClasses="character") else
        data.table(ID=character(0), platform=character(0))
if(nrow(DROP)){ setnames(DROP, 1:2, c("ID", "platform")); DROP[, key := paste(ID, platform)] } else DROP[, key := character(0)]
rd_sscore <- function(f, plat){ s <- fread(f); d <- data.table(IID=as.character(s[[1]]), raw=as.numeric(s$SCORE1_AVG))
  if(nrow(DROP)) d <- d[!(paste(IID, plat) %chin% DROP$key)]; d }
getz <- function(pid, coh){
  L <- lapply(PLATS, function(p){ f <- file.path(SSC, sprintf("%s__%s__%s.sscore", pid, coh, p))
    if(!file.exists(f)) return(NULL); s <- rd_sscore(f, p)
    data.table(IID=s$IID, z=as.numeric(zin(s$raw)), plat=p) })
  L <- L[!vapply(L, is.null, TRUE)]; if(!length(L)) return(NULL)
  X <- rbindlist(L)[is.finite(z)]
  X[, .(z=mean(z), n_plat=.N), by=IID][, z := zin(z)][]
}

## A second, independent derivation of the same merged score, written as explicit sums
## rather than with the grouping helpers above.  04 compares the two and reports the
## largest disagreement, so the score that enters the models is shown to be reproducible
## from the score files alone by a different code path.

z_alt <- function(pid, coh){
  L <- list()
  for(p in PLATS){
    f <- file.path(SSC, sprintf("%s__%s__%s.sscore", pid, coh, p))
    if(!file.exists(f)) next
    s  <- fread(f, showProgress=FALSE)
    id <- as.character(s[[1]]); rw <- as.numeric(s$SCORE1_AVG)
    if(nrow(DROP)){ k <- !(paste(id, p) %chin% DROP$key); id <- id[k]; rw <- rw[k] }
    ok <- is.finite(rw); id <- id[ok]; rw <- rw[ok]
    if(!length(rw)) next
    mu <- sum(rw)/length(rw); s2 <- sum((rw - mu)^2)/(length(rw) - 1)
    if(!(s2 > 0)) next
    L[[length(L) + 1]] <- data.table(IID=id, v=(rw - mu)/sqrt(s2))
  }
  if(!length(L)) return(NULL)
  A <- rbindlist(L)[, .(v=sum(v)/.N), by=IID]
  mu <- sum(A$v)/nrow(A); s2 <- sum((A$v - mu)^2)/(nrow(A) - 1)
  A[, .(IID, z_alt=(v - mu)/sqrt(s2))]
}
score_files <- function(pid, coh){
  f <- vapply(PLATS, function(p) file.path(SSC, sprintf("%s__%s__%s.sscore", pid, coh, p)), character(1))
  unname(f[file.exists(f)])
}

## Logistic-model diagnostics.  Separation is judged on the log-odds change per standard
## deviation of each covariate, so a covariate on a small scale cannot appear separated
## merely for being on a small scale.  A cell that separates or fails to converge is
## stopped and reported as such; no penalised estimator is substituted.

logit_status <- function(m, wmsg){
  cf <- coef(m)
  X  <- model.matrix(m)
  sdx <- apply(X, 2, stats::sd)
  se <- rep(NA_real_, length(cf)); names(se) <- names(cf)
  sm <- tryCatch(summary(m)$coefficients, error=function(e) NULL)
  if(!is.null(sm)) se[rownames(sm)] <- sm[, 2]
  sdx <- sdx[names(cf)]; sdx[!is.finite(sdx) | sdx <= 0] <- NA_real_
  fv <- fitted(m)
  na_coef   <- any(is.na(cf))
  converged <- isTRUE(m$converged) && !isTRUE(m$boundary) && !na_coef && all(is.finite(cf[!is.na(cf)]))
  sep_fit   <- any(fv < SEP_P | fv > 1 - SEP_P)
  sep_coef  <- any(abs(cf) * sdx > SEP_BSD, na.rm=TRUE)
  sep_se    <- any(se * sdx > SEP_SESD, na.rm=TRUE)
  sep_warn  <- any(grepl("fitted probabilities numerically 0 or 1", wmsg, fixed=TRUE)) ||
               any(grepl("did not converge", wmsg, fixed=TRUE))
  list(converged=converged, boundary=isTRUE(m$boundary), na_coef=na_coef, iter=m$iter,
       separation = sep_fit || sep_coef || sep_se || sep_warn,
       sep_fitted=sep_fit, sep_coef=sep_coef, sep_se=sep_se, sep_warning=sep_warn,
       min_fitted=min(fv), max_fitted=max(fv),
       max_scaled_coef=suppressWarnings(max(abs(cf)*sdx, na.rm=TRUE)),
       max_scaled_se=suppressWarnings(max(se*sdx, na.rm=TRUE)))
}
fit_glm <- function(fml, d){
  wmsg <- character(0)
  m <- withCallingHandlers(glm(fml, data=d, family=binomial(),
                               control=glm.control(epsilon=GLM_EPS, maxit=GLM_MAXIT)),
         warning=function(w){ wmsg <<- c(wmsg, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(m=m, warn=wmsg)
}

sep_selftest <- local({
  x  <- c(seq(-2, -0.1, length.out=20), seq(0.1, 2, length.out=20))
  ds <- data.frame(y=c(rep(0L, 20), rep(1L, 20)), x=x)
  dn <- data.frame(y=rep(c(0L, 1L), 20),          x=x)
  gs <- fit_glm(y ~ x, ds); gn <- fit_glm(y ~ x, dn)
  list(separated_flagged = isTRUE(logit_status(gs$m, gs$warn)$separation),
       clean_not_flagged = isTRUE(!logit_status(gn$m, gn$warn)$separation))
})

stopifnot(sep_selftest$separated_flagged, sep_selftest$clean_not_flagged)
