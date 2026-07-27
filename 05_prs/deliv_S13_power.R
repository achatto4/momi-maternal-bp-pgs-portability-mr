#!/usr/bin/env Rscript
# ============================================================
# deliv_S13_power.R  [ B25 -> Supp Table S13 ]  ** BROUGHT FORWARD 2026-07-19 **
# Minimum detectable effect and power, per cohort, for the observational and MR arms.
#
# WHY IT MOVED AHEAD OF B19-B24. S9b established that the causal effect is heterogeneous
# across cohorts (I2 = 70-92%, and 89-91% on continuous birth weight), and the paper's thesis
# is now that the DIRECTION of effect transports while the MAGNITUDE does not. That thesis
# rests on COMPARING effects between cohorts. It works on the observational arm -- those
# estimates are precise enough to reject homogeneity. Whether it works on the MR arm is
# unknown, and the arithmetic is not encouraging: the PRS explains only 1-6% of BP variance,
# cohorts run 1,448-3,950, and a Wald ratio inherits imprecision from BOTH stages.
#
# If per-cohort MR intervals span everything from protective to strongly harmful, then the
# heterogeneity finding stands on the observational arm ALONE. That is still a paper, but a
# different one from the triangulation currently planned -- and it is far cheaper to learn it
# now than after B19-B24 are built.
#
# THE HEADLINE NUMBER is power_MR_at_obs: if the true causal effect equals what we observe
# observationally in that cohort, what is our power to detect it by MR there? Anything under
# ~50% means a null MR result in that cohort is uninformative and must not be read as absence
# of effect.
#
# Reads analytic_mothers, transfer_grid, tables/S9_confounding.tsv.
# Writes: tables/S13_power.tsv
#   Rscript deliv_S13_power.R
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

ALPHA <- 0.05
CRIT  <- qchisq(1-ALPHA, df=1)          # 3.841
NCP80 <- (qnorm(1-ALPHA/2) + qnorm(0.80))^2   # 7.849, the NCP giving 80% power

momi_deliverable("S13_power", script="05_prs/deliv_S13_power.R",
                 inputs="analytic_mothers;transfer_grid", P=P, stop_on_error=FALSE,
                 body=function(ctx){
  A  <- momi_read_intermediate("analytic_mothers", P)
  G  <- momi_read_intermediate("transfer_grid", P)
  pd <- MOMI_BPDEF_DEFAULT
  f9 <- file.path(P$tables, "S9_confounding.tsv")
  S9 <- if(file.exists(f9)) fread(f9) else NULL

  OUT <- list(PTB=list(col="PTB", type="bin", lb=FALSE),
              LBW=list(col="LBW", type="bin", lb=TRUE),
              SGA=list(col="SGA", type="bin", lb=TRUE),
              BWT=list(col="BWT", type="lin", lb=TRUE))

  rows <- list()
  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    bpcol <- paste0(substr(tr,1,1), "_", pd)
    if(!bpcol %in% names(A)) next
    sid <- momi_instrument(coh, tr)
    ## first-stage R2 for the CHOSEN instrument at the primary definition
    gr <- G[cohort==coh & score_id==sid & definition==pd & trait==tr]
    R2 <- if(nrow(gr) && is.finite(gr$incR2[1])) gr$incR2[1] else NA_real_

    for(onm in names(OUT)){
      O <- OUT[[onm]]
      d <- A[cohort==coh & genotyped==1]
      if(O$lb) d <- d[livebirth==1]
      d <- d[is.finite(get(bpcol)) & !is.na(get(O$col))]
      if(nrow(d) < 50 || !is.finite(R2)) next

      N   <- nrow(d)
      x10 <- d[[bpcol]]/10                 # exposure in 10-mmHg units
      vx  <- var(x10, na.rm=TRUE)

      ## observed observational effect for this cell (per 10 mmHg, log scale)
      obs <- if(!is.null(S9)) S9[scope==coh & trait==tr & definition==pd & outcome==onm &
                                 sample=="common" & model=="M3_plus_EDU", est] else numeric(0)
      obs <- if(length(obs) && is.finite(obs[1])) obs[1] else NA_real_

      if(O$type=="bin"){
        p  <- mean(d[[O$col]]==1, na.rm=TRUE)
        k  <- sum(d[[O$col]]==1, na.rm=TRUE)
        if(!is.finite(p) || p<=0 || p>=1 || k < 10) next
        ## OBSERVATIONAL: logistic on x10. NCP = beta^2 * var(x10) * N * p(1-p)
        b_obs_min <- sqrt(NCP80 / (N * p*(1-p) * vx))
        ## MR (Wald): the instrument explains R2 of the exposure, so the effective
        ## exposure variance available to the instrument is R2 * var(x10).
        b_mr_min  <- sqrt(NCP80 / (N * p*(1-p) * vx * R2))
        ## power of the MR arm IF the truth equals the observed observational effect
        ncp_at    <- if(is.finite(obs)) obs^2 * N * p*(1-p) * vx * R2 else NA_real_
        eff_obs   <- if(is.finite(obs)) exp(obs) else NA_real_
        unit      <- "OR per 10 mmHg"
        mdeo      <- exp(b_obs_min); mdem <- exp(b_mr_min)
      } else {
        k  <- NA_integer_; p <- NA_real_
        sy <- sd(d[[O$col]], na.rm=TRUE)
        ## continuous: SE(beta) = sy / (sqrt(N) * sd(x10)); MR divides by sqrt(R2)
        b_obs_min <- sqrt(NCP80) * sy / (sqrt(N) * sqrt(vx))
        b_mr_min  <- sqrt(NCP80) * sy / (sqrt(N) * sqrt(vx * R2))
        ncp_at    <- if(is.finite(obs)) (obs^2 * N * vx * R2) / (sy^2) else NA_real_
        eff_obs   <- obs
        unit      <- "grams per 10 mmHg"
        mdeo      <- b_obs_min; mdem <- b_mr_min
      }

      rows[[paste(coh,tr,onm)]] <- data.table(
        cohort=coh, cohort_anc=unname(MOMI_ANC[coh]), trait=tr, outcome=onm,
        instrument=sid, unit=unit,
        N=N, cases=k, prevalence=round(100*p,1),
        firstStage_R2pct = round(100*R2,3),
        observed_effect  = round(eff_obs,3),
        MDE_observational = round(mdeo,3),
        MDE_MR            = round(mdem,3),
        power_MR_at_obs   = round(pchisq(CRIT, df=1, ncp=ncp_at, lower.tail=FALSE),3),
        MR_adequate       = fifelse(is.na(ncp_at), NA,
                              pchisq(CRIT, 1, ncp=ncp_at, lower.tail=FALSE) >= 0.80))
    }
  }
  if(!length(rows)) return(list(skip=TRUE, reason="no cells computable — check transfer_grid"))
  S13 <- rbindlist(rows, fill=TRUE)
  setorder(S13, trait, outcome, -power_MR_at_obs)
  out <- momi_write_table(S13, "S13_power", P)
  ## ADDED 2026-07-19: also publish as an INTERMEDIATE so B22 can join power onto every MR
  ## cell. A deliverable that only writes a .tsv cannot be depended on -- momi_fingerprint
  ## hashes intermediates, so a TSV-only dependency would leave B22 silently un-stale when
  ## the power model changes. Making it an intermediate makes the dependency real.
  momi_save_intermediate(S13, "power_grid", P)

  ## ---- console ----
  cat("\n== MR power per cohort, IF the truth equals the observed observational effect ==\n")
  pw <- dcast(S13[outcome %in% c("PTB","LBW","SGA","BWT")],
              cohort + cohort_anc ~ trait + outcome, value.var="power_MR_at_obs")
  print(pw)
  cat("\nPower below ~0.50 means a null MR result in that cohort is UNINFORMATIVE -- it cannot\n",
      "distinguish 'no effect' from 'we could not have seen it'. It must not be reported as\n",
      "evidence of absence, and it cannot contribute to a heterogeneity claim.\n", sep="")

  cat("\n== minimum detectable effect, observational vs MR (primary outcomes, SBP) ==\n")
  print(S13[trait=="SBP" & outcome %in% c("LBW","PTB"),
            .(cohort, outcome, N, cases, firstStage_R2pct,
              observed=observed_effect, MDE_obs=MDE_observational, MDE_MR)])
  cat("\nMDE_MR is inflated relative to MDE_observational by 1/sqrt(R2). At R2 = 2%, that is a\n",
      "~7-fold penalty; at 6%, ~4-fold. This is the cost of instrumenting a weakly-predicted\n",
      "exposure and it is unavoidable, not a modelling choice.\n", sep="")

  ok <- S13[MR_adequate %in% TRUE]
  cat(sprintf("\ncells with >=80%% MR power: %d of %d (%.0f%%)\n",
              nrow(ok), nrow(S13), 100*nrow(ok)/nrow(S13)))
  if(nrow(ok)) print(ok[, .(cohort, trait, outcome, power_MR_at_obs)])
  cat("\nIf this is empty or near-empty, per-cohort MR cannot support the heterogeneity claim\n",
      "and the thesis rests on the observational arm alone. B19-B24 need restructuring, not\n",
      "just rerunning.\n", sep="")

  list(n=nrow(S13),
       key=sprintf("cells=%d; MR power median=%.2f range %.2f-%.2f; adequate(>=80%%)=%d",
                   nrow(S13), median(S13$power_MR_at_obs, na.rm=TRUE),
                   min(S13$power_MR_at_obs, na.rm=TRUE), max(S13$power_MR_at_obs, na.rm=TRUE),
                   nrow(ok)),
       outputs=c(basename(out), "power_grid.rds"))
})
