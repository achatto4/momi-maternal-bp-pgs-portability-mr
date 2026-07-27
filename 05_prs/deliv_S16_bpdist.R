#!/usr/bin/env Rscript
# ============================================================
# deliv_S16_bpdist.R  [ B15 -> Supp Table S16 + Supp Figure SF3 ]
# BP distribution comparison, AMANHI-Bangladesh vs GAPPS-Bangladesh, per definition x trait:
# mean, SD, skew, N per cohort, plus standardized mean difference, variance ratio,
# and Kolmogorov-Smirnov D/p. Also emits SF3 density overlays across all cohorts.
#
# FINDING -- SUPERSEDED TWICE, current statement below (2026-07-20 numbers).
#
# The 2026-07-18 run described "both extremes diverge while the central summaries converge",
# with DBP median reaching genuine concordance (KS p = 0.134). THAT WAS AN ARTEFACT OF
# POSTPARTUM CONTAMINATION: ~40% of AMANHI's BP readings were postnatal, and removing them
# (B01's antenatal restriction, 2026-07-19) moved the mean SBP SMD from -0.179 to -0.309 and
# eliminated the one non-significant row. No definition now reaches concordance.
#
# CURRENT FINDING, and both halves must be stated because either alone misleads:
#
#   1. AVERAGING GENUINELY ATTENUATES THE GAP. SBP SMD runs -0.503 (`last`), -0.336
#      (`early`), -0.320 (`residual`), -0.309 (`mean`), -0.235 (`median`), -0.211 (`mean-2`).
#      From the worst single reading to the mean of the first two is a 58% reduction. So the
#      old intuition that summarising over readings brings the cohorts closer was RIGHT.
#
#   2. IT NEVER CLOSES THE GAP. Every definition retains at least a fifth of a standard
#      deviation, and every KS test is significant. GAPPS-Bangladesh mothers simply have
#      higher measured BP than AMANHI-Bangladesh mothers at every summary.
#
# Saying only (1) reproduces the retracted "averaging reconciles the cohorts" claim; saying
# only (2) hides that the choice of definition matters a great deal to how large the gap looks.
#
# WHAT THIS RULES OUT, AND IT IS THE REASON S16 MATTERS FOR F3. The variance ratio is ~1.0 at
# every definition (SBP mean 1.049, residual 1.044) and is 0.828 for `last`, i.e. AMANHI has
# MORE BP variance there while having LOWER R2. So "GAPPS has more BP variance for the score
# to explain" -- the obvious competing explanation for F3's transferability gap -- is dead.
# Combined with S3's Fst of 0.00017 showing the two cohorts are genetically indistinguishable,
# F3's gap can be neither genetic nor variance-driven, which leaves measurement timing.
#
# DENOMINATOR (decision 2026-07-18): the PRIMARY comparison uses ALL mothers in the two
# cohorts, because this is a BP-measurement-protocol contrast that does not depend on
# genotyping. A Part I-restricted consistency check (genotyped + valid antenatal BP -- the
# exact F1/B07 rule) is reported alongside in the *_pI columns. The 2026-07-18 run showed
# the two agree to <=0.01 SMD on every well-covered definition, so genotyping selection is
# not driving the contrast.
#
# COVERAGE FLAG: a definition can clear the n>=30 bar while representing a tiny, selected
# slice of its cohort. `<20wk` is exactly this -- 225/2934 (7.7%) of AMANHI-Bangladesh
# against 3631/3633 (99.9%) of GAPPS-Bangladesh, i.e. self-selected early presenters vs
# essentially the whole cohort. Such rows are flagged LOW_COVERAGE and must not be read as
# a like-for-like comparison.
#
#   Writes: tables/S16_bpdist.tsv, figures/SF3_bpdist.{png,pdf}
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

MOMI_S16_MINCOV <- 0.25   # flag a row if either arm covers <25% of its cohort

momi_deliverable("S16_bpdist", script="05_prs/deliv_S16_bpdist.R",
                 inputs="analytic_mothers", P=P, body=function(ctx){
  A <- momi_read_intermediate("analytic_mothers", P)
  FOCUS <- c("AMANHI-Bangladesh","GAPPS-Bangladesh")

  ## ---- Part I flag: EXACTLY the F1/B07 rule (genotyped + any valid mean BP) ----
  ## NA-safe: genotyped is NA if the sscore dir was unavailable at build time.
  A[, partI := {
      g <- !is.na(genotyped) & genotyped == 1
      b <- is.finite(S_mean) | is.finite(D_mean)
      g & b }]

  skew <- function(x){ x<-x[is.finite(x)]; if(length(x)<3) return(NA_real_)
                       m<-mean(x); s<-sd(x); if(s==0) return(NA_real_); mean((x-m)^3)/s^3 }
  defs <- c(early="first", mean="mean", `mean-2`="mn2", median="median",
            last="last", `<20wk`="lt20", `>=20wk`="ge20", residual="resid")

  ## p-values from ks.test underflow to exactly 0 at these sample sizes; reporting "0" in a
  ## published table is wrong. Render the floor explicitly.
  fmt_p <- function(p){
    if(!is.finite(p)) return(NA_character_)
    if(p < 2.2e-16) "<2.2e-16" else format(signif(p,3), scientific=TRUE)
  }

  ## ---- one two-sample comparison; NULL if either arm is too small ----
  cmp <- function(x, y, min_n=30){
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    if(length(x) < min_n || length(y) < min_n) return(NULL)
    sp <- sqrt(((length(x)-1)*var(x) + (length(y)-1)*var(y)) / (length(x)+length(y)-2))
    ks <- suppressWarnings(stats::ks.test(x, y))
    list(nx=length(x), mx=mean(x), sx=sd(x), skx=skew(x),
         ny=length(y), my=mean(y), sy=sd(y), sky=skew(y),
         smd=if(is.finite(sp) && sp>0) (mean(x)-mean(y))/sp else NA_real_,
         vr =if(var(y)>0) var(x)/var(y) else NA_real_,
         ksD=unname(ks$statistic), ksp=ks$p.value)
  }

  rows <- list()
  for(tr in c("SBP","DBP")){
    pre <- substr(tr,1,1)
    ## cohort denominator for coverage = mothers with a valid MEAN reading of this trait
    den <- vapply(FOCUS, function(co)
             sum(is.finite(A[cohort==co, get(paste0(pre,"_mean"))])), numeric(1))

    for(nm in names(defs)){
      col <- paste0(pre,"_",defs[[nm]])
      if(!col %in% names(A)) next

      ## primary: all mothers in the two cohorts
      a <- cmp(A[cohort==FOCUS[1], get(col)], A[cohort==FOCUS[2], get(col)])
      if(is.null(a)) next

      ## consistency check: Part I analytic sample only
      p <- cmp(A[cohort==FOCUS[1] & partI, get(col)],
               A[cohort==FOCUS[2] & partI, get(col)])

      cov_x <- a$nx/den[[1]]; cov_y <- a$ny/den[[2]]
      flag  <- if(min(cov_x, cov_y) < MOMI_S16_MINCOV) "LOW_COVERAGE" else ""

      rows[[paste(tr,nm)]] <- data.table(
        trait=tr, definition=nm, flag=flag,
        ## ---- primary (all mothers) ----
        N_AMANHI_B=a$nx, cov_AMANHI_B=round(100*cov_x,1),
        mean_AMANHI_B=round(a$mx,2), sd_AMANHI_B=round(a$sx,2), skew_AMANHI_B=round(a$skx,2),
        N_GAPPS_B=a$ny, cov_GAPPS_B=round(100*cov_y,1),
        mean_GAPPS_B=round(a$my,2), sd_GAPPS_B=round(a$sy,2), skew_GAPPS_B=round(a$sky,2),
        SMD=round(a$smd,3), var_ratio=round(a$vr,3),
        KS_D=round(a$ksD,3), KS_p=fmt_p(a$ksp),
        ## ---- consistency check (Part I analytic sample) ----
        N_AMANHI_B_pI    = if(is.null(p)) NA_integer_ else p$nx,
        mean_AMANHI_B_pI = if(is.null(p)) NA_real_    else round(p$mx,2),
        N_GAPPS_B_pI     = if(is.null(p)) NA_integer_ else p$ny,
        mean_GAPPS_B_pI  = if(is.null(p)) NA_real_    else round(p$my,2),
        SMD_pI           = if(is.null(p)) NA_real_    else round(p$smd,3),
        KS_D_pI          = if(is.null(p)) NA_real_    else round(p$ksD,3))
    }
  }
  T <- rbindlist(rows)
  out1 <- momi_write_table(T, "S16_bpdist", P)

  ## ---- SF3: density overlays, all cohorts, earliest reading vs mean ----
  dd <- rbindlist(list(
    A[is.finite(S_first), .(cohort, trait="SBP", which="earliest reading", v=S_first)],
    A[is.finite(S_mean),  .(cohort, trait="SBP", which="mean of readings", v=S_mean)],
    A[is.finite(D_first), .(cohort, trait="DBP", which="earliest reading", v=D_first)],
    A[is.finite(D_mean),  .(cohort, trait="DBP", which="mean of readings", v=D_mean)]))
  dd <- dd[cohort %in% MOMI_COH_ALL]
  dd[, Cohort := unname(MOMI_DISPLAY[cohort])]
  ## export densities on a grid so report/build_figures.py can redraw SF3 in the shared theme
  denS <- dd[is.finite(v), {
    g <- density(v, n=256); .(sbp=round(g$x,2), density=signif(g$y,6))
  }, by=.(cohort, trait, which)]
  momi_write_table(denS, "SF3_bp_density", P)
  g <- ggplot(dd, aes(v, colour=Cohort)) + geom_density(linewidth=0.6) +
    facet_grid(trait ~ which, scales="free") +
    labs(x="Blood pressure (mmHg)", y="Density",
         title="SF3. BP distributions by cohort: earliest reading vs mean",
         subtitle="All mothers with a valid reading (not restricted to the genotyped Part I sample)") +
    theme_minimal(base_size=10) + theme(legend.position="bottom", legend.title=element_blank())
  outs <- momi_save_fig(g, "SF3_bpdist", width=9.5, height=6, P=P)

  ## ---- key: the SMD profile across definitions, both traits, plus flag count ----
  prof <- function(tr) paste(sprintf("%s=%.2f", T[trait==tr]$definition, T[trait==tr]$SMD),
                             collapse=",")
  list(n=nrow(T),
       key=sprintf("SMD SBP[%s] DBP[%s]; flagged=%d",
                   prof("SBP"), prof("DBP"), sum(nzchar(T$flag))),
       outputs=c(basename(out1), basename(outs)))
})
