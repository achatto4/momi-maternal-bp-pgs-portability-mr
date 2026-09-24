## -------------------------------------------------------------------------------
## 03_pgs_portability.R
##
## How well published blood-pressure polygenic scores transfer to these cohorts.
##
## For each cohort, trait and score, two nested linear models are fitted in the cleaned
## sample -- blood pressure on the standardised score and maternal age, against blood
## pressure on maternal age alone -- and the incremental R-squared, the coefficient in
## millimetres of mercury per standard deviation of the score, its standard error and
## t-based 95% confidence interval, the F statistic and the p value are reported. The
## primary panel of five score families is fitted on the women who have all five scores
## of that trait, so the families are compared in the same women; score-specific sample
## sizes are reported alongside. Confidence intervals for the incremental R-squared come
## from mother-level bootstrap resampling within cohort, with one index matrix per cohort
## and trait shared by all five families of the primary panel.
##
## A separate diagnostic compares the two genotyping platforms for the women who carry
## both, per cohort and score. It has no effect on the models above.
## -------------------------------------------------------------------------------
if (!exists("MOMI_ROOT")) MOMI_ROOT <- getwd()
if (!exists("config")) source(file.path(MOMI_ROOT,
  if (file.exists(file.path(MOMI_ROOT, "config.R"))) "config.R" else "config.example.R"))
source(file.path(MOMI_ROOT, "analysis", "00_functions.R"))

t0    <- Sys.time()
B     <- config$bootstrap_replicates
SEED0 <- config$bootstrap_seed

SITE    <- c("1"="GAPPS-Zambia","2"="GAPPS-Bangladesh","3"="AMANHI-Pakistan",
             "4"="AMANHI-Bangladesh","5"="AMANHI-Pemba","6"="THSTI-India")
GENO_PROBE <- "PGS004603"
FAM     <- data.table(family=c("European","South Asian","East Asian","Diverse ancestry","MVP"),
                      SBP=c("PGS004603","PGS004830","PGS002376","PGS003968","PGS002238"),
                      DBP=c("PGS004604","PGS004758","PGS002362","PGS003964","PGS002239"))
SUPP    <- "PGS005008"
BPCOL   <- c(SBP="S_mean", DBP="D_mean")

EPI_COLS <- c("SITE_CODE","PARTICIPANT_ID","PREGNANCY_ID","SBP","DBP","GA_HDLK_NEW","VISITDT","PW_AGE",
              "SINGLE_TWIN","DEL_DATE","DATE_LMP","GAGEBRTH_NEW","BWT_MEASURE_DATE")
OUTCOME_COLS <- c("PTB_NEW","BIRTH_WEIGHT","SGA_10_NEW","SGA_3_NEW","BIRTH_OUTCOME","PREG_OUTCOME","PE_CAT","CHRON_HTN")

asdate  <- function(x){
  if(inherits(x, c("Date","IDate"))) return(as.Date(x))
  suppressWarnings(as.Date(as.character(x), format="%Y-%m-%d"))
}
W_  <- function(dt, name) fwrite(dt, file.path(RESULTS, name), sep="\t")

## The exposure: each woman's mean of her valid antenatal readings, built from the same
## rules as elsewhere, reading only the fields the exposure and the sample need.

hdr <- names(fread(EPI, nrows=0))
miss <- setdiff(EPI_COLS, hdr); if(length(miss)) stop("EPI lacks columns: ", paste(miss, collapse=", "))
raw <- fread(EPI, na.strings=MOMI_NA, select=EPI_COLS)[PREGNANCY_ID == 1]
getcol <- function(nm) if(nm %in% names(raw)) num(raw[[nm]]) else rep(NA_real_, nrow(raw))
raw[, IID := as.character(PARTICIPANT_ID)]
raw[, SBP := getcol("SBP")]; raw[, DBP := getcol("DBP")]

raw[, pair_invalid0 := is.finite(SBP) & is.finite(DBP) & DBP >= SBP]
PAIR_RULE_COUNTS <- raw[, .(invalid_pairs=sum(pair_invalid0)), by=.(SITE_CODE)]
raw[pair_invalid0 == TRUE, `:=`(SBP = NA_real_, DBP = NA_real_)]
garead <- getcol("GA_HDLK_NEW")
raw[, GA := ifelse(garead >= 0 & garead <= 315, garead, NA_real_)]
raw[, VIS := getcol("VISITDT")]
raw[, AGEv := getcol("PW_AGE")]; raw[, GAdv := getcol("GAGEBRTH_NEW")]; raw[, TWINv := getcol("SINGLE_TWIN")]
raw[, cohort := SITE[as.character(SITE_CODE)]]
setorder(raw, IID, VIS, GA)
raw[, VISd := asdate(VISITDT)]
raw[, DELd := if("DEL_DATE" %in% names(raw)) asdate(DEL_DATE) else as.Date(NA)]
raw[, LMPd := if("DATE_LMP" %in% names(raw)) asdate(DATE_LMP) else as.Date(NA)]
raw[, BWTd := if("BWT_MEASURE_DATE" %in% names(raw)) asdate(BWT_MEASURE_DATE) else as.Date(NA)]
raw[, DELeff := DELd]
raw[is.na(DELeff) & !is.na(LMPd) & is.finite(GAdv), DELeff := LMPd + GAdv]
raw[is.na(DELeff) & !is.na(BWTd), DELeff := BWTd]
raw[, DELeff := { v <- DELeff[!is.na(DELeff)]; if(length(v)) v[1] else as.Date(NA) }, by=IID]
raw[, visit_class := fifelse(is.na(VISd) | is.na(DELeff), "unknown",
                      fifelse(VISd > DELeff, "postnatal", "antenatal"))]
raw[, GA_lmp := as.numeric(VISd - LMPd)]
raw[visit_class == "unknown" & is.finite(GA_lmp) & GA_lmp >= 0 & GA_lmp <= 294, visit_class := "antenatal"]
mnv  <- function(v) if(all(is.na(v))) NA_real_ else mean(v, na.rm=TRUE)
ANTE <- raw[visit_class == "antenatal"]
Sd   <- ANTE[, .(S_mean=mnv(SBP)), by=IID]; Dd <- ANTE[, .(D_mean=mnv(DBP)), by=IID]
cov  <- raw[, .(cohort=cohort[1], AGE=firstnn(AGEv), TWIN=firstnn(TWINv)), by=IID]
cov[, twin := as.integer(!is.na(TWIN) & TWIN != 1)]
ALLW <- Reduce(function(a, b) merge(a, b, by="IID", all.x=TRUE), list(cov, Sd, Dd))
ALLW[, hasBP := is.finite(S_mean) | is.finite(D_mean)]

## Genotype availability and the merged scores.

DROP <- if(!is.null(DROPF) && file.exists(DROPF)) fread(DROPF, colClasses="character") else
        data.table(ID=character(0), platform=character(0))
if(nrow(DROP)){ setnames(DROP, 1:2, c("ID", "platform")); DROP[, key := paste(ID, platform)] } else DROP[, key := character(0)]
rd_sscore <- function(f, plat){ s <- fread(f); d <- data.table(IID=as.character(s[[1]]), raw=as.numeric(s$SCORE1_AVG))
  if(nrow(DROP)) d <- d[!(paste(IID, plat) %chin% DROP$key)]; d }
gids <- character(0)
for(coh in COH5) for(pl in PLATS){
  f <- file.path(SSC, sprintf("%s__%s__%s.sscore", GENO_PROBE, coh, pl))
  if(file.exists(f)) gids <- union(gids, rd_sscore(f, pl)$IID)
}
ALLW[, genotyped := as.integer(IID %in% gids)]
AN <- ALLW[twin == 0 & genotyped == 1 & hasBP == TRUE, .(IID, cohort, AGE, S_mean, D_mean)]
ELIG_N <- nrow(AN)
KEEP <- if(!is.null(KEEPF) && file.exists(KEEPF)) as.character(fread(KEEPF, header=FALSE, colClasses="character")[[1]]) else character(0)
if(length(KEEP)) AN <- AN[IID %chin% KEEP]
CLEAN_N <- nrow(AN)
getz <- function(pid, coh){
  L <- lapply(PLATS, function(p){ f <- file.path(SSC, sprintf("%s__%s__%s.sscore", pid, coh, p))
    if(!file.exists(f)) return(NULL); s <- rd_sscore(f, p)
    data.table(IID=s$IID, z=as.numeric(zin(s$raw)), plat=p) })
  L <- L[!vapply(L, is.null, TRUE)]; if(!length(L)) return(NULL)
  X <- rbindlist(L)[is.finite(z)]
  X[, .(z=mean(z), n_plat=.N), by=IID][, z := zin(z)][]
}
SCORES <- c(FAM$SBP, FAM$DBP, SUPP)
ZL <- list(); ZQC <- list()
for(coh in COH5) for(sid in SCORES){
  z <- getz(sid, coh); if(is.null(z)) stop("no .sscore files for ", sid, " / ", coh)
  ZL[[paste(coh, sid)]] <- z[, .(IID, z)]
  za <- merge(AN[cohort == coh, .(IID)], z, by="IID")
  ZQC[[length(ZQC) + 1]] <- data.table(cohort=coh, score_id=sid, n_scored=nrow(z), mean_scored=mean(z$z), sd_scored=sd(z$z),
                                        n_dual_platform_scored=sum(z$n_plat == 2), n_in_sample=nrow(za),
                                        mean_in_sample=mean(za$z), sd_in_sample=sd(za$z))
}
ZQC <- rbindlist(ZQC)
W_(ZQC, "pgs_portability_score_qc.tsv")

## The models and the bootstrap.

transfer <- function(z, y, age){
  m <- lm(y ~ z + age); m0 <- lm(y ~ age); sm <- summary(m)
  cc <- sm$coefficients; ci <- confint(m, "z", level=0.95)
  list(incR2=sm$r.squared - summary(m0)$r.squared, beta=cc["z", 1], se=cc["z", 2], beta_lo=ci[1], beta_hi=ci[2],
       F=(cc["z", 1]/cc["z", 2])^2, p=cc["z", 4])
}
inc_fast <- function(y, z, age){
  f0 <- .lm.fit(cbind(1, age), y); f1 <- .lm.fit(cbind(1, z, age), y)
  tss <- sum((y - mean(y))^2)
  if(f0$rank < 2 || f1$rank < 3 || !(tss > 0)) return(NA_real_)
  (sum(f0$residuals^2) - sum(f1$residuals^2)) / tss
}
boot_inc <- function(y, z, age, IDX) vapply(seq_len(ncol(IDX)), function(b){ i <- IDX[, b]; inc_fast(y[i], z[i], age[i]) }, numeric(1))
make_idx <- function(n, seed){ set.seed(seed); matrix(sample.int(n, n * B, replace=TRUE), nrow=n) }
summ <- function(bt) list(boot_B=B, boot_n_success=sum(is.finite(bt)),
                          boot_lo=unname(quantile(bt, 0.025, na.rm=TRUE, type=7)), boot_hi=unname(quantile(bt, 0.975, na.rm=TRUE, type=7)),
                          boot_mean=mean(bt, na.rm=TRUE), boot_sd=sd(bt, na.rm=TRUE))

RES <- list(); BOOT <- list(); SAMP <- list()
for(ic in seq_along(COH5)) for(it in 1:2){
  coh <- COH5[ic]; tr <- c("SBP","DBP")[it]; seed <- SEED0 + 10L*ic + it
  sids <- FAM[[tr]]; allsids <- if(tr == "SBP") c(sids, SUPP) else sids
  D <- AN[cohort == coh, .(IID, AGE, y=get(BPCOL[[tr]]))]
  for(s in allsids){ D <- merge(D, ZL[[paste(coh, s)]], by="IID", all.x=TRUE); setnames(D, "z", paste0("z_", s)) }
  setkey(D, IID)
  base <- is.finite(D$y) & !is.na(D$AGE)
  spec <- sapply(allsids, function(s) base & is.finite(D[[paste0("z_", s)]]))
  common <- base & Reduce(`&`, lapply(sids, function(s) is.finite(D[[paste0("z_", s)]])))
  Dc <- D[common]; n <- nrow(Dc)
  IDX <- make_idx(n, seed)
  SAMP[[length(SAMP) + 1]] <- as.data.table(c(list(cohort=coh, trait=tr, n_bp_age=sum(base), n_common=n, seed=seed),
                                              as.list(setNames(as.integer(colSums(spec)), paste0("n_", allsids)))))
  for(k in seq_along(sids)){
    s <- sids[k]; z <- Dc[[paste0("z_", s)]]
    est <- transfer(z, Dc$y, Dc$AGE); bt <- boot_inc(Dc$y, z, Dc$AGE, IDX)
    RES[[length(RES) + 1]] <- data.table(cohort=coh, cohort_group=GROUP[[coh]], trait=tr, family=FAM$family[k], score_id=s,
      role="Primary portability panel", exposure=BPCOL[[tr]], sample="common (all five primary scores)",
      N=n, N_score_specific=sum(spec[, s]), z_mean_in_sample=mean(z), z_sd_in_sample=sd(z), as.data.table(est),
      as.data.table(summ(bt)), boot_seed=seed, boot_indices="shared by the five primary families")
    BOOT[[length(BOOT) + 1]] <- data.table(cohort=coh, trait=tr, score_id=s, replicate=seq_len(B), incR2=bt)
  }
  if(tr == "SBP"){
    Ds <- D[spec[, SUPP]]; same <- identical(Ds$IID, Dc$IID)
    IDXs <- if(same) IDX else make_idx(nrow(Ds), seed + 5L)
    z <- Ds[[paste0("z_", SUPP)]]; est <- transfer(z, Ds$y, Ds$AGE); bt <- boot_inc(Ds$y, z, Ds$AGE, IDXs)
    RES[[length(RES) + 1]] <- data.table(cohort=coh, cohort_group=GROUP[[coh]], trait=tr, family="Multi-ancestry", score_id=SUPP,
      role="Supplementary SBP-only portability analysis", exposure=BPCOL[[tr]],
      sample=if(same) "own complete cases (identical to the SBP common sample)" else "own complete cases",
      N=nrow(Ds), N_score_specific=nrow(Ds), z_mean_in_sample=mean(z), z_sd_in_sample=sd(z), as.data.table(est),
      as.data.table(summ(bt)), boot_seed=if(same) seed else seed + 5L,
      boot_indices=if(same) "same matrix as the SBP primary families" else "own matrix")
    BOOT[[length(BOOT) + 1]] <- data.table(cohort=coh, trait=tr, score_id=SUPP, replicate=seq_len(B), incR2=bt)
  }
}
RES <- rbindlist(RES); BOOT <- rbindlist(BOOT); SAMP <- rbindlist(SAMP, fill=TRUE)
W_(RES, "pgs_portability_results.tsv"); W_(BOOT, "pgs_portability_bootstrap_replicates.tsv"); W_(SAMP, "pgs_portability_samples.tsv")
SUMM <- AN[, .(n=.N, mean_S_mean=mean(S_mean, na.rm=TRUE), n_S_mean=sum(is.finite(S_mean)),
               mean_D_mean=mean(D_mean, na.rm=TRUE), n_D_mean=sum(is.finite(D_mean)),
               mean_AGE=mean(AGE, na.rm=TRUE), n_AGE=sum(!is.na(AGE))), by=cohort]
W_(SUMM, "pgs_portability_sample_summary.tsv")

## A score whose coefficient is negative is reported rather than adjusted.

neg <- RES[beta < 0, .(cohort, trait, score_id, beta)]
W_(if(nrow(neg)) neg else data.table(cohort=character(0), trait=character(0), score_id=character(0), beta=numeric(0)),
   "pgs_portability_negative_betas.tsv")

## Cross-platform agreement diagnostic.

TRAIT  <- c(setNames(rep("SBP", 5), FAM$SBP), setNames(rep("DBP", 5), FAM$DBP), setNames("SBP", SUPP))
FAMILY <- c(setNames(FAM$family, FAM$SBP), setNames(FAM$family, FAM$DBP), setNames("Multi-ancestry", SUPP))
ROLE   <- c(setNames(rep("Primary portability panel", 10), c(FAM$SBP, FAM$DBP)),
            setNames("Supplementary SBP-only portability analysis", SUPP))
EXPECT_DUAL <- NULL
AGR <- list()
for(coh in COH5) for(sid in SCORES){
  zp <- lapply(PLATS, function(p){ f <- file.path(SSC, sprintf("%s__%s__%s.sscore", sid, coh, p))
    if(!file.exists(f)) return(data.table(IID=character(0), z=numeric(0)))
    s <- rd_sscore(f, p); data.table(IID=s$IID, z=as.numeric(zin(s$raw))) })
  m <- merge(merge(AN[cohort == coh, .(IID)], zp[[1]], by="IID"), zp[[2]], by="IID", suffixes=c("_gsa", "_wgs"))
  m <- m[is.finite(z_gsa) & is.finite(z_wgs)]; n <- nrow(m); d <- m$z_gsa - m$z_wgs
  r  <- if(n >= 3) cor(m$z_gsa, m$z_wgs) else NA_real_
  ci <- if(n >= 4 && is.finite(r) && abs(r) < 1) tanh(atanh(r) + c(-1, 1) * qnorm(0.975) / sqrt(n - 3)) else c(NA_real_, NA_real_)
  AGR[[length(AGR) + 1]] <- data.table(cohort=coh, cohort_group=GROUP[[coh]], score_id=sid, trait=TRAIT[[sid]], family=FAMILY[[sid]],
    role=ROLE[[sid]], n_dual=n, r=r, r_lo=ci[1], r_hi=ci[2], mean_diff_gsa_minus_wgs=if(n) mean(d) else NA_real_,
    sd_diff=if(n > 1) sd(d) else NA_real_)
}
AGR <- rbindlist(AGR)
W_(AGR, "pgs_portability_platform_agreement.tsv")
VA <- data.table(check=c("agreement rows (5 cohorts x 11 scores)",
                         "cohorts whose 11 scores all have the same dual-platform N"),
                 observed=c(nrow(AGR), sum(vapply(COH5, function(coh) uniqueN(AGR[cohort == coh]$n_dual) == 1L, logical(1)))),
                 expected=c(55L, 5L))
VA[, match := observed == expected]
W_(VA, "pgs_portability_platform_agreement_checks.tsv")

W_(data.table(
  item = c("eligible_women","cleaned_women","dropped_genotype_records","bootstrap_replicates",
           "seed_base","seed_rule","R","data.table","run_seconds"),
  value = c(ELIG_N, CLEAN_N, nrow(DROP), B, SEED0,
            "base + 10*cohort_index + trait_index (+5 for a separate supplementary-score matrix)",
            R.version.string, as.character(packageVersion("data.table")),
            sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "secs"))))),
   "pgs_portability_run_info.tsv")
cat(sprintf("03: portability results for %d cohort x trait x score cells\n", nrow(RES)))
