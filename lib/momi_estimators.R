#!/usr/bin/env Rscript
# ============================================================
# momi_estimators.R — the SINGLE source of estimation truth for the paper.
# Every deliverable computes transferability, reliability, MR, and observational
# effects by calling THESE functions, so Table 2 / Figure 2 / S4 / S5 can never
# disagree on how an R2 or a Wald ratio was formed. Sources momi_config.R helpers.
#   source(file.path(PIPE,"lib/momi_estimators.R"))
# ============================================================
suppressMessages(library(data.table))

## ---------- (1) TRANSFERABILITY ----------
# Incremental R2 of a standardized PRS on measured BP, adjusting for covariates
# (default age). Returns incR2, beta (mmHg per SD PRS), se, first-stage F, N.
#   z    : standardized PRS (numeric)
#   bp   : measured BP (numeric)
#   cov  : optional data.frame/data.table of covariates (e.g. AGE)
momi_transfer <- function(z, bp, cov=NULL, min_n=30){
  d <- data.table(z=as.numeric(z), bp=as.numeric(bp))
  if(!is.null(cov)) d <- cbind(d, as.data.table(cov))
  d <- d[is.finite(z) & is.finite(bp)]
  d <- d[complete.cases(d)]
  if(nrow(d) < min_n)
    return(data.table(n=nrow(d), incR2=NA_real_, beta=NA_real_, se=NA_real_, F=NA_real_, p=NA_real_))
  covn <- setdiff(names(d), c("z","bp"))
  rhs  <- if(length(covn)) paste(c("z",covn), collapse="+") else "z"
  rhs0 <- if(length(covn)) paste(covn, collapse="+") else "1"
  m  <- lm(as.formula(paste("bp ~", rhs)),  data=d)
  m0 <- lm(as.formula(paste("bp ~", rhs0)), data=d)
  cc <- summary(m)$coefficients
  b  <- cc["z",1]; se <- cc["z",2]
  data.table(n=nrow(d),
             incR2 = summary(m)$r.squared - summary(m0)$r.squared,
             beta  = b, se = se, F = (b/se)^2, p = cc["z",4])
}

## ---------- (2) RELIABILITY (ICC) for disattenuation ----------
# One-way random-effects ICC of repeated antenatal BP readings within mother,
# used to disattenuate transferability R2 (R2_true ~= R2_obs / ICC).
# Closed-form variance components (unbalanced) — NO per-mother design matrix, so it
# is fast and memory-light on thousands of mothers.
#   long : data.table with columns id, y (one row per reading)
momi_icc <- function(long, id="id", y="y"){
  d <- data.table(id=as.character(long[[id]]), y=as.numeric(long[[y]]))[is.finite(y)]
  if(nrow(d) < 5) return(NA_real_)
  grp <- d[, .(ni=.N, mi=mean(y)), by=id]
  k <- nrow(grp); N <- nrow(d)
  if(k < 5 || max(grp$ni) < 2 || N==k) return(NA_real_)   # need repeats
  grand <- mean(d$y)
  SSB <- grp[, sum(ni*(mi-grand)^2)]                       # between-mother
  d2  <- merge(d, grp[, .(id, mi)], by="id")
  SSW <- d2[, sum((y-mi)^2)]                               # within-mother
  MSB <- SSB/(k-1); MSW <- SSW/(N-k)
  n0  <- (N - sum(grp$ni^2)/N)/(k-1)                       # avg cluster size (unbalanced)
  icc <- (MSB - MSW)/(MSB + (n0-1)*MSW)
  max(min(icc, 1), 0)
}
## ---------- (2b) reliability of a MEAN of k readings (Spearman-Brown) ----------
# The ICC is the reliability of a SINGLE reading. Our exposures are means of ~4-5 readings,
# and averaging cancels noise, so the mean is markedly more reliable than any one reading.
#   rel(mean of k) = k*ICC / (1 + (k-1)*ICC)
# e.g. ICC 0.391 with k = 4.5 gives 0.74, not 0.39.
momi_rel_mean <- function(icc, k){
  if(is.na(icc) || icc <= 0 || is.na(k) || k < 1) return(icc)
  min(k*icc / (1 + (k-1)*icc), 1)
}

## ---------- (2c) DISATTENUATED R2 ----------
# Corrects R2 for measurement error in the OUTCOME of the transferability regression
# (bp ~ z). Standard correction for attenuation, Spearman 1904; for blood pressure
# specifically the landmark is MacMahon et al. Lancet 1990 (PMID 1969518), which showed a
# single BP reading dilutes observed associations by ~40%. See ref/methods_citations.tsv.
#
# NB measurement error in an OUTCOME does not bias the slope -- beta (mmHg per SD of PRS)
# and F are unbiased and must NOT be "corrected". Only R2 is attenuated.
#
# DEFECT FIXED 2026-07-19: this divided by the raw ICC, i.e. treated the exposure as a
# single reading. That over-corrected by roughly 2x -- GAPPS-Bangladesh SBP came out at a
# disattenuated 16.5%, higher than published European estimates, which was the tell. Pass k
# (mean readings per mother) to use the reliability of the mean instead. k = 1 reproduces
# the old behaviour, so any caller that does not supply k is unchanged.
momi_disattenuate <- function(r2, icc, k = 1){
  rel <- momi_rel_mean(icc, k)
  if(is.na(rel) || rel <= 0) return(r2)
  pmin(r2 / rel, 1)
}

## ---------- (3) MENDELIAN RANDOMIZATION (Wald ratio) ----------
# Reduced form (PRS->outcome) unchanged; first stage (PRS->BP) sets the denominator.
# Returns per-cohort building blocks; combine across cohorts with momi_meta_iv().
#   d      : per-mother data.table with columns z, <bpcol>, <outcol>, covariates
#   type   : "bin" (glm binomial, OR) or "lin" (lm)
#   cov    : character vector of covariate names (default "AGE")
momi_wald <- function(d, bpcol, outcol, type=c("bin","lin"), cov="AGE",
                      min_n=40, min_case=8){
  type <- match.arg(type)
  keep <- c("z", bpcol, outcol, cov)
  d <- d[, ..keep]
  setnames(d, c("z","bp","y", cov))
  d <- d[is.finite(z) & is.finite(y)]
  d <- d[complete.cases(d[, c("z","y",..cov)])]
  ## C1 (2026-07-26): both stages must be estimated on the SAME sample. The exposure `bp`
  ## (the GA-standardised residual) is missing for women without a valid GA, and that
  ## missingness is strongly cohort-differential. Restricting to finite(bp) here means the
  ## reduced form and the first stage share one sample, so the Wald ratio bo/be is a ratio of
  ## commensurable quantities. (Previously the reduced form used all finite-y rows and only the
  ## first stage was restricted to finite(bp).) This also makes the reported n the true sample.
  d <- d[is.finite(bp)]
  covrhs <- if(length(cov)) paste0("+", paste(cov, collapse="+")) else ""
  if(type=="bin" && (sum(d$y,na.rm=TRUE) < min_case || nrow(d) < min_n))
    return(NULL)
  if(type=="lin" && nrow(d) < min_n) return(NULL)
  # reduced form
  rf <- tryCatch(
    if(type=="bin") summary(glm(as.formula(paste0("y ~ z", covrhs)), binomial, d))$coef
    else            summary(lm(as.formula(paste0("y ~ z", covrhs)), d))$coef,
    error=function(e) NULL)
  if(is.null(rf) || !("z" %in% rownames(rf))) return(NULL)
  bo <- rf["z",1]; so <- rf["z",2]
  # first stage (PRS -> BP) on the SAME finite-BP sample used above
  fs <- tryCatch(summary(lm(as.formula(paste0("bp ~ z", covrhs)), d))$coef,
                 error=function(e) NULL)
  if(is.null(fs) || !("z" %in% rownames(fs))) return(NULL)
  be <- fs["z",1]; see <- fs["z",2]
  theta   <- bo/be
  se_th   <- sqrt(so^2/be^2 + bo^2*see^2/be^4)   # delta method
  data.table(n=nrow(d), bo=bo, so=so, be=be, see=see,
             Ffs=(be/see)^2, theta=theta, se_theta=se_th)
}

## ---------- (4) inverse-variance meta over cohorts ----------
momi_meta_iv <- function(b, se){
  ok <- is.finite(b) & is.finite(se) & se>0; b<-b[ok]; se<-se[ok]
  if(!length(b)) return(NULL)
  w <- 1/se^2; mu <- sum(w*b)/sum(w); s <- sqrt(1/sum(w))
  list(k=length(b), b=mu, se=s, z=mu/s, p=2*pnorm(-abs(mu/s)))
}

## ---------- (5) OBSERVATIONAL effect (adjusted) ----------
# OR (bin) or beta (lin) per `per` units of BP, adjusting for a confounder formula.
#   d    : per-mother data.table with <bpcol>, <outcol>, confounders
#   conf : character vector of confounder columns (e.g. MOMI_CONF_CORE)
momi_obs <- function(d, bpcol, outcol, type=c("bin","lin"), conf=character(0),
                     per=10, min_n=40, min_case=8){
  type <- match.arg(type)
  keep <- unique(c(bpcol, outcol, conf))
  x <- d[, ..keep]
  setnames(x, c("bp","y", conf))
  x <- x[is.finite(bp) & is.finite(y)]
  x <- x[complete.cases(x)]
  x[, bpu := bp/per]
  if(type=="bin" && (sum(x$y,na.rm=TRUE) < min_case || nrow(x) < min_n)) return(NULL)
  if(type=="lin" && nrow(x) < min_n) return(NULL)
  covrhs <- if(length(conf)) paste0("+", paste(conf, collapse="+")) else ""
  fit <- tryCatch(
    if(type=="bin") glm(as.formula(paste0("y ~ bpu", covrhs)), binomial, x)
    else            lm(as.formula(paste0("y ~ bpu", covrhs)), x),
    error=function(e) NULL)
  if(is.null(fit)) return(NULL)
  cc <- summary(fit)$coef
  if(!("bpu" %in% rownames(cc))) return(NULL)
  b <- cc["bpu",1]; se <- cc["bpu",2]
  data.table(n=nrow(x), est=b, se=se,
             lo=b-1.96*se, hi=b+1.96*se, p=cc["bpu",4],
             OR=exp(b), OR_lo=exp(b-1.96*se), OR_hi=exp(b+1.96*se))
}

invisible(TRUE)
