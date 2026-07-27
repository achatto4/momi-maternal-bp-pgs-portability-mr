#!/usr/bin/env Rscript
# ============================================================
# deliv_F3_diagnostic.R  [ B16 -> Figure 3, why portability differs ]
# AMANHI-Bangladesh vs GAPPS-Bangladesh (same country & ancestry, different BP protocol),
# so the contrast isolates measurement TIMING rather than genetics:
#   A  BP-by-gestational-age trajectory (mechanism: when BP is measured)
#   B  R2 by BP definition for the SAME score (SAS SBP), with bootstrap 95% CI,
#      PLUS a bootstrap CI on the BETWEEN-COHORT DIFFERENCE dR2 = GAPPS - AMANHI
#   C  BP distribution overlay, earliest reading vs mean
#
# WHY dR2 (decision 2026-07-18): panel B originally showed two separate CIs and invited the
# reader to judge the gap by whether they overlap. That is not a valid test -- non-overlap
# implies a difference, but overlap does NOT imply no difference, and at the `mean`
# definition the two intervals nearly touched (4.86 vs 4.91). The cohorts are INDEPENDENT
# samples, so the difference of their bootstrap replicate vectors is a valid bootstrap
# distribution for dR2. We therefore report dR2 with a 95% CI per definition and make the
# claim directly. B raised 400 -> 1000 for a stable difference interval (~65 s).
#
# FINDING (2026-07-18 run): R2 rises monotonically with the number of readings averaged in
# BOTH cohorts (AMANHI 1.55 -> 2.56 -> 3.13 -> 3.45; PreSSMat 4.94 -> 5.94 -> 6.05 -> 6.26
# for early -> mean-2 -> median -> mean). That is attenuation by measurement error, and it
# explains B14: a single early reading loses more to noise than it gains from being closer
# to a pre-pregnancy baseline. B16 therefore CONFIRMS B14 rather than contradicting it.
#
# Reads bp_readings_long, prs_z, analytic_mothers (frozen intermediates).
#   Writes: figures/F3_diagnostic.{png,pdf} (or F3_panelA/B/C if patchwork unavailable)
#           tables/F3_r2_by_definition.tsv, tables/F3_delta_r2.tsv
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

MOMI_F3_B <- 1000L   # bootstrap replicates

momi_deliverable("F3_diagnostic", script="05_prs/deliv_F3_diagnostic.R",
                 inputs="bp_readings_long;prs_z;analytic_mothers", P=P, body=function(ctx){
  A   <- momi_read_intermediate("analytic_mothers", P)
  Z   <- momi_read_intermediate("prs_z", P)
  LNG <- momi_read_intermediate("bp_readings_long", P)
  FOCUS <- c("AMANHI-Bangladesh","GAPPS-Bangladesh")
  disp  <- function(c) unname(MOMI_DISPLAY[c])

  ## ---- Panel A: BP-by-GA trajectory ----
  traj <- LNG[cohort %in% FOCUS & trait=="SBP" & is.finite(GA_days),
              .(bp=mean(bp), n=.N), by=.(cohort, wk=floor(GA_days/7))][n>=20]
  traj[, Cohort := disp(cohort)]
  ## export panel-A data so report/build_figures.py (Python) can redraw it in the shared theme
  momi_write_table(traj[, .(cohort, ga_week=wk, mean_sbp=round(bp,2), n)], "F3_bp_by_ga", P)
  pA <- ggplot(traj, aes(wk, bp, colour=Cohort)) +
    geom_vline(xintercept=20, linetype=2, colour="grey55") +
    geom_line(linewidth=0.7) + geom_point(size=1.1) +
    annotate("text", x=20.4, y=-Inf, label="20 wk", hjust=0, vjust=-0.6, size=2.7, colour="grey40") +
    labs(x="Gestational age at reading (weeks)", y="Mean SBP (mmHg)",
         title="A. When BP is measured (SBP by gestational age)",
         subtitle="All antenatal readings; GA-weeks with <20 readings suppressed") +
    theme_minimal(base_size=10) + theme(legend.position="top", legend.title=element_blank())

  ## ---- Panel B: R2 by definition, SAME score, bootstrap CI + between-cohort dR2 ----
  sid  <- MOMI_SAS$SBP                       # one score for both cohorts -> isolates timing
  defmap <- c(early="S_first", mean="S_mean", `mean-2`="S_mn2",
              median="S_median", last="S_last", `>=20wk`="S_ge20")

  ## return the FULL replicate vector (in R2 %), so the difference distribution can be built
  boot_reps <- function(z, bp, age, B=MOMI_F3_B){
    d <- data.table(z=z, bp=bp, age=age); d <- d[is.finite(z) & is.finite(bp) & is.finite(age)]
    if(nrow(d) < 60) return(rep(NA_real_, B))
    replicate(B, { i <- sample.int(nrow(d), nrow(d), replace=TRUE)
      100 * momi_transfer(d$z[i], d$bp[i], cov=data.table(AGE=d$age[i]))$incR2 })
  }
  qci <- function(v) as.numeric(quantile(v, c(.025,.975), na.rm=TRUE))

  set.seed(1)
  rows <- list(); reps <- list()
  for(coh in FOCUS){
    z  <- Z[score_id==sid & cohort==coh, .(IID,z)]
    ac <- A[cohort==coh, c("IID","AGE", unname(defmap)), with=FALSE]
    d  <- merge(z, ac, by="IID")
    for(nm in names(defmap)){
      col <- defmap[[nm]]
      tt  <- momi_transfer(d$z, d[[col]], cov=d[, .(AGE)])
      v   <- boot_reps(d$z, d[[col]], d$AGE)
      reps[[paste(coh,nm,sep="|")]] <- v
      ci  <- qci(v)
      rows[[paste(coh,nm)]] <- data.table(cohort=coh, Cohort=disp(coh), definition=nm,
        R2pct=round(100*tt$incR2,3), lo=round(ci[1],3), hi=round(ci[2],3), N=tt$n)
    }
  }
  RB <- rbindlist(rows)
  RB[, definition := factor(definition, levels=names(defmap))]
  out_rb <- momi_write_table(RB, "F3_r2_by_definition", P)

  ## ---- between-cohort difference: GAPPS - AMANHI, per definition ----
  ## The two cohorts are independent samples, so differencing their bootstrap replicate
  ## vectors elementwise yields a valid bootstrap distribution for dR2.
  drows <- lapply(names(defmap), function(nm){
    vA <- reps[[paste(FOCUS[1],nm,sep="|")]]
    vG <- reps[[paste(FOCUS[2],nm,sep="|")]]
    if(is.null(vA) || is.null(vG) || all(is.na(vA)) || all(is.na(vG)))
      return(data.table(definition=nm, R2_AMANHI_B=NA_real_, R2_GAPPS_B=NA_real_,
                        dR2=NA_real_, dR2_lo=NA_real_, dR2_hi=NA_real_, excludes_zero=NA))
    dv <- vG - vA; ci <- qci(dv)
    data.table(definition=nm,
      R2_AMANHI_B = RB[cohort==FOCUS[1] & definition==nm, R2pct],
      R2_GAPPS_B  = RB[cohort==FOCUS[2] & definition==nm, R2pct],
      dR2 = round(mean(dv, na.rm=TRUE),3),
      dR2_lo = round(ci[1],3), dR2_hi = round(ci[2],3),
      excludes_zero = !(ci[1] <= 0 && ci[2] >= 0))
  })
  DR <- rbindlist(drows)
  out_dr <- momi_write_table(DR, "F3_delta_r2", P)

  sub_b <- sprintf("Genotyped mothers with a %s score; %d bootstrap replicates. Between-cohort dR2 in tables/F3_delta_r2.tsv",
                   sid, MOMI_F3_B)
  pB <- ggplot(RB, aes(definition, R2pct, fill=Cohort)) +
    geom_col(position=position_dodge(0.8), width=0.72) +
    geom_errorbar(aes(ymin=lo, ymax=hi), position=position_dodge(0.8), width=0.18, linewidth=0.4) +
    labs(x="BP definition", y="Incremental R² (%)",
         title=sprintf("B. Transferability by BP definition (same score: %s), bootstrap 95%% CI", sid),
         subtitle=sub_b) +
    theme_minimal(base_size=10) + theme(legend.position="top", legend.title=element_blank())

  ## ---- Panel C: distribution overlay, earliest vs mean ----
  dd <- rbind(
    A[cohort %in% FOCUS & is.finite(S_first), .(Cohort=disp(cohort), which="earliest reading", v=S_first)],
    A[cohort %in% FOCUS & is.finite(S_mean),  .(Cohort=disp(cohort), which="mean of readings", v=S_mean)])
  ## export panel-C densities (mean-of-readings) on a grid for report/build_figures.py
  denC <- rbindlist(lapply(FOCUS, function(coh){
    v <- A[cohort==coh & is.finite(S_mean), S_mean]
    if(length(v) < 30) return(NULL)
    g <- density(v, n=256)
    data.table(cohort=coh, sbp=round(g$x,2), density=signif(g$y,5))
  }))
  momi_write_table(denC, "F3_bp_density", P)
  pC <- ggplot(dd, aes(v, colour=Cohort, fill=Cohort)) +
    geom_density(alpha=0.16, linewidth=0.6) + facet_wrap(~which) +
    ## CAPTION CORRECTED 2026-07-20. It previously read "extremes diverge, central summaries
    ## converge". S16 refutes the second half: averaging ATTENUATES the between-cohort gap
    ## (standardised mean difference -0.50 on the last reading, -0.21 on the mean of the first
    ## two) but never closes it -- every definition retains a gap of at least a fifth of a
    ## standard deviation, and the mean definition sits at -0.31. "Converge" overstated it and
    ## was an artefact of the pre-2026-07-19 postpartum contamination.
    labs(x="SBP (mmHg)", y="Density",
         title="C. BP distributions: averaging narrows the gap but does not close it",
         subtitle="Between-cohort SMD: -0.50 (last reading), -0.34 (earliest), -0.31 (mean), -0.21 (mean of first two). All mothers with a valid reading; panel B is genotyped-only") +
    theme_minimal(base_size=10) + theme(legend.position="top", legend.title=element_blank())

  ## ---- combine (patchwork if available, else save panels separately) ----
  if(requireNamespace("patchwork", quietly=TRUE)){
    g <- patchwork::wrap_plots(list(pA,pB,pC), ncol=1)
    outs <- momi_save_fig(g, "F3_diagnostic", width=8.5, height=11, P=P)
  } else {
    outs <- c(momi_save_fig(pA,"F3_panelA",width=8,height=3.6,P=P),
              momi_save_fig(pB,"F3_panelB",width=8,height=3.6,P=P),
              momi_save_fig(pC,"F3_panelC",width=8,height=3.6,P=P))
    message("patchwork not installed — saved F3 panels A/B/C separately")
  }

  dm <- DR[definition=="mean"]
  list(n=nrow(RB),
       key=sprintf("mean-def R² AMANHI=%.2f%% GAPPS=%.2f%%; dR2=%.2f [%.2f,%.2f] excl0=%s; sig_defs=%s",
                   dm$R2_AMANHI_B, dm$R2_GAPPS_B, dm$dR2, dm$dR2_lo, dm$dR2_hi,
                   as.character(dm$excludes_zero),
                   paste(DR[excludes_zero %in% TRUE, definition], collapse=",")),
       outputs=c(basename(outs), basename(out_rb), basename(out_dr)))
})
