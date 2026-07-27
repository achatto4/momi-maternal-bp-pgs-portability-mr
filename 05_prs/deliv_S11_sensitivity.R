#!/usr/bin/env Rscript
# ============================================================
# deliv_S11_sensitivity.R  [ B23 -> S11, is the pooled MR result real or one cohort? ]
#
# S10 produced two pooled results that look like findings: SBP -> birthweight
# (-79 g per 10 mmHg, p=0.005) and, more weakly, SBP -> LBW (OR 1.39, p=0.045).
# This module exists to attack them.
#
# ------------------------------------------------------------------
# WHY THIS IS THE DECIDING ANALYSIS, NOT A FOOTNOTE
#
# In the per-cohort panel, FOUR of the six reduced-form signals at p<0.05 come from
# GAPPS-Bangladesh alone. That cohort also has by far the strongest instrument (F 166-247
# vs 17-81 elsewhere) and the highest power in S13 (0.28-0.30 vs ~0.05-0.17). So the
# question is not "is the pooled estimate significant" -- it is "is the pooled estimate
# ANYTHING MORE THAN GAPPS-Bangladesh with four noisy cohorts averaged in". Those are
# different papers and the leave-one-out below is what distinguishes them.
#
# WHAT THE CONVENTIONAL DIAGNOSTICS MISSED, AND WHY WE DO NOT RELY ON THEM
#
#   I2 IS USELESS AT THIS SAMPLE SIZE. S10 reports I2 = 0 for SBP->BWT, which reads as
#   "perfectly homogeneous". But GAPPS-Zambia's point estimate is +163 g while the other
#   four run -56 to -449 g -- a cohort on the OPPOSITE SIDE OF THE NULL. I2 fails to see it
#   because the standard errors are so wide that every cohort is consistent with every
#   other. I2 = 0 here means "cannot detect heterogeneity with k=5", NOT "no heterogeneity".
#   S11 therefore reports the SIGN CONCORDANCE and the leave-one-out spread, which are
#   descriptive and cannot be talked out of.
#
#   F > 10 DID NOT PROTECT US. Every cell in S10 clears the conventional threshold, yet
#   AMANHI-Pemba returns -449 g per 10 mmHg for DBP (F=22) against -89 g for SBP (F=56) --
#   same cohort, same outcome, estimate quadrupling as the first stage weakens. The ratio
#   estimator inflates smoothly as F falls; 10 is not a cliff edge below which there is
#   trouble and above which there is none. S11 reports estimates restricted to progressively
#   stronger instruments so the reader can watch the estimate move.
# ------------------------------------------------------------------
# WHAT IS COMPUTED, for every trait x outcome cell in the pooled table:
#
#   1. LEAVE-ONE-OUT      -- k refits, each dropping one cohort. Reports the largest shift
#                            and which cohort caused it, and whether p<0.05 survives all k.
#   2. SIGN CONCORDANCE   -- how many cohorts share the pooled estimate's sign. Immune to
#                            the SE inflation that defeats I2.
#   3. ANCESTRY SUBGROUP  -- SAS-only and AFR-only pooled estimates. The portability gap
#                            (Part I) predicts the AFR arm is near-useless; this quantifies
#                            how much of the pooled result is really the SAS cohorts.
#   4. F-THRESHOLD LADDER -- pooled estimate restricted to cohorts with F >= 10, 25, 50, 100.
#                            A monotone drift as the threshold rises is the signature of
#                            weak-instrument inflation contributing to the pooled estimate.
#
# ALL FOUR ARE PRE-SPECIFIED HERE AS A SET AND REPORTED FOR EVERY CELL, INCLUDING THE NULL
# ONES. That matters: running leave-one-out only on the significant cells, or reporting only
# the exclusion that happens to help, is how sensitivity analyses become advocacy. The
# F-ladder in particular is rule-based rather than naming GAPPS-Zambia, even though Zambia
# is the cohort we noticed -- a threshold rule cannot be accused of being drawn around the
# inconvenient point after the fact.
#
#   Writes: tables/S11_sensitivity.tsv, tables/S11b_leave_one_out.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

## inverse-variance pool + Cochran's Q on an arbitrary subset; returns a one-row data.table
pool_of <- function(s, label){
  s <- s[is.finite(theta_per10) & is.finite(se_per10) & se_per10 > 0]
  if(!nrow(s)) return(NULL)
  if(nrow(s) == 1)
    return(data.table(subset=label, k=1L, cohorts=s$cohort[1], N=s$N[1],
                      theta=s$theta_per10[1], se=s$se_per10[1],
                      p=2*pnorm(-abs(s$theta_per10[1]/s$se_per10[1])),
                      Q=NA_real_, I2=NA_real_, n_neg=sum(s$theta_per10<0), n_pos=sum(s$theta_per10>0)))
  fe <- momi_meta_iv(s$theta_per10, s$se_per10)
  w  <- 1/s$se_per10^2
  Q  <- sum(w*(s$theta_per10 - fe$b)^2); dfQ <- nrow(s)-1
  data.table(subset=label, k=fe$k, cohorts=paste(sort(s$cohort), collapse=","),
             N=sum(s$N), theta=fe$b, se=fe$se, p=fe$p,
             Q=Q, I2=max(0, 100*(Q-dfQ)/Q),
             n_neg=sum(s$theta_per10<0), n_pos=sum(s$theta_per10>0))
}

momi_deliverable("S11_sensitivity", script="05_prs/deliv_S11_sensitivity.R",
                 inputs="mr_panel;mr_pooled;ptb_subtype", P=P, stop_on_error=FALSE,
                 body=function(ctx){

  MR <- tryCatch(momi_read_intermediate("mr_panel", P), error=function(e) NULL)
  if(is.null(MR) || !nrow(MR))
    return(list(skip=TRUE, reason="no mr_panel.rds -- run B22 first"))
  MR[, panel := "perinatal"]

  ## ---- EXTENDED 2026-07-20: the PTB-subtype result needs this MORE than the main panel ----
  ## B19 found indicated PTB at OR 2.19 (p=0.0025) against spontaneous 0.67 (ns), difference
  ## p=0.0026 -- the paper's strongest result. But GAPPS-Bangladesh contributes 265 of the
  ## 448 indicated cases (59%), the same single-cohort concentration that cost SBP->BWT its
  ## significance on deletion. A finding that rests on one cohort's obstetric practice is a
  ## different claim from one five cohorts agree on, and only leave-one-out separates them.
  ## The subtype panel is folded in here (subtype renamed to `outcome`) so it goes through
  ## the identical machinery rather than a parallel implementation that could drift.
  SUBP <- tryCatch(momi_read_intermediate("ptb_subtype", P), error=function(e) NULL)
  if(!is.null(SUBP) && nrow(SUBP)){
    sp <- copy(SUBP)
    setnames(sp, "subtype", "outcome")
    sp[, `:=`(panel = "PTB subtype", type = "bin")]
    keep <- intersect(names(MR), names(sp))
    MR <- rbind(MR[, ..keep], sp[, ..keep], fill=TRUE)
  } else cat("\nNB ptb_subtype.rds not found -- S11 covers the perinatal panel only.\n")

  cells <- unique(MR[, .(panel, trait, outcome, type)])
  main <- list(); loo <- list()

  for(i in seq_len(nrow(cells))){
    pn <- cells$panel[i]; tr <- cells$trait[i]; oc <- cells$outcome[i]; ty <- cells$type[i]
    s  <- MR[panel==pn & trait==tr & outcome==oc &
             is.finite(theta_per10) & is.finite(se_per10) & se_per10>0]
    if(nrow(s) < 2) next
    tag <- function(d, lab) if(is.null(d)) NULL else
             cbind(data.table(panel=pn, trait=tr, outcome=oc, type=ty), d)

    ## ---- 0. the full pooled estimate, repeated here so every row is self-contained ----
    full <- pool_of(s, "all cohorts")
    main[[paste(pn,tr,oc,"all")]] <- tag(full)

    ## ---- 1. leave-one-out ----
    lo <- rbindlist(lapply(s$cohort, function(cc)
      tag(pool_of(s[cohort != cc], paste0("drop ", cc)))), fill=TRUE)
    if(nrow(lo)){
      lo[, dropped := sub("^drop ", "", subset)]
      ## influence = how far the pooled estimate moves when this cohort is removed
      lo[, shift := theta - full$theta]
      lo[, pct_shift := 100*shift/abs(full$theta)]
      lo[, sig := p < 0.05]
      ## share of the cell's events carried by the dropped cohort -- the number that says
      ## whether an influential cohort is influential because it is large or because it
      ## disagrees. n_case is absent for continuous outcomes, hence the guard.
      if("n_case" %in% names(s) && any(is.finite(s$n_case)))
        lo[, pct_cases_dropped :=
             round(100*s$n_case[match(dropped, s$cohort)]/sum(s$n_case, na.rm=TRUE), 1)]
      loo[[paste(pn,tr,oc)]] <- lo
    }

    ## ---- 2. ancestry subgroups ----
    s[, anc := MOMI_ANC[cohort]]
    main[[paste(pn,tr,oc,"SAS")]] <- tag(pool_of(s[anc=="SAS"], "SAS only"))
    main[[paste(pn,tr,oc,"AFR")]] <- tag(pool_of(s[anc=="AFR"], "AFR only"))

    ## ---- 3. F-threshold ladder ----
    for(thr in c(10, 25, 50, 100))
      main[[paste(pn,tr,oc,"F",thr)]] <- tag(pool_of(s[F >= thr], sprintf("F >= %d", thr)))
  }

  S11  <- rbindlist(main, fill=TRUE)
  LOO  <- rbindlist(loo,  fill=TRUE)
  if(!nrow(S11)) return(list(skip=TRUE, reason="no cell had >=2 estimable cohorts"))

  ## ---- verdict per cell: does the headline survive every single-cohort deletion? ----
  V <- LOO[, .(loo_worst_p = max(p),               # worst case across the k refits
               loo_all_sig = all(p < 0.05),
               loo_theta_min = min(theta), loo_theta_max = max(theta),
               most_influential = dropped[which.max(abs(shift))],
               max_abs_pct_shift = round(max(abs(pct_shift)),1)),
           by=.(panel, trait, outcome)]
  S11 <- merge(S11, V, by=c("panel","trait","outcome"), all.x=TRUE)
  setorder(S11, panel, outcome, trait, subset)

  out  <- momi_write_table(S11, "S11_sensitivity", P)
  out2 <- momi_write_table(LOO, "S11b_leave_one_out", P)

  ## ---------------- console ----------------
  cat("\n=== pooled estimate under each restriction (theta per 10 mmHg) ===\n")
  print(S11[, .(panel, trait, outcome, subset, k, theta=round(theta,3), p=signif(p,3),
                sign_split=paste0(n_neg,"-/",n_pos,"+"))])

  cat("\n=== leave-one-out: does each cell survive dropping any ONE cohort? ===\n")
  print(V[, .(panel, trait, outcome, loo_all_sig, worst_p=signif(loo_worst_p,3),
              theta_range=sprintf("%.1f to %.1f", loo_theta_min, loo_theta_max),
              most_influential, max_abs_pct_shift)])

  cat("\n=== the PTB-subtype result, cohort by cohort (the 59%-of-cases question) ===\n")
  sub <- LOO[panel=="PTB subtype"]
  if(nrow(sub)) print(sub[, .(trait, outcome, dropped, k, theta=round(theta,3),
                              OR=round(exp(theta),3), p=signif(p,3),
                              pct_shift=round(pct_shift,1),
                              pct_cases = if("pct_cases_dropped" %in% names(sub))
                                            pct_cases_dropped else NA_real_)])
  cat("\nGAPPS-Bangladesh carries 265 of 448 indicated cases. If `indicated` stays significant\n",
      "with it dropped, the finding is not one cohort's obstetric practice. If it does not,\n",
      "say so plainly -- an OR of 2.2 resting on one site is still worth reporting, but as\n",
      "that site's result, and the health-system-dependence caveat becomes central.\n", sep="")

  cat("\nHOW TO READ THIS. loo_all_sig=FALSE means the pooled p<0.05 depends on the presence\n",
      "of at least one particular cohort, and the cell should NOT be reported as a pooled\n",
      "finding -- report it as that cohort's result. A large max_abs_pct_shift with\n",
      "loo_all_sig=TRUE is milder: the estimate is unstable but the direction holds.\n", sep="")

  cat("\n=== the cells S10 flagged as findings ===\n")
  key <- S11[outcome %in% c("BWT","LBW")]
  if(nrow(key)) print(key[, .(trait, outcome, subset, k, theta=round(theta,2), p=signif(p,3))])

  cat("\nSIGN CONCORDANCE is the diagnostic to trust here, not I2. A cell where one cohort\n",
      "sits on the far side of the null still returns I2=0 when the SEs are this wide.\n", sep="")

  ## The summary line previously read "LOO-robust cells = 0 of 12", counting cells that were
  ## NEVER significant as leave-one-out failures. That is meaningless -- a null cell cannot
  ## fail a robustness test it never passed -- and it twice made results look worse than they
  ## were. Restrict the denominator to cells whose FULL pooled estimate reached p<0.05, which
  ## is the only population for which "does it survive deletion" is a question.
  fullsig <- S11[subset=="all cohorts" & is.finite(p) & p < 0.05,
                 .(panel, trait, outcome)]
  Vs <- merge(V, fullsig, by=c("panel","trait","outcome"))
  cat(sprintf("\ncells significant in the full pool: %d; of those, robust to dropping any one cohort: %d\n",
              nrow(Vs), sum(Vs$loo_all_sig, na.rm=TRUE)))
  if(nrow(Vs)) print(Vs[, .(panel, trait, outcome, loo_all_sig,
                            worst_p=signif(loo_worst_p,3), most_influential,
                            max_abs_pct_shift)])

  list(n=nrow(S11),
       key=sprintf("cells=%d; of %d cells significant in full pool, %d survive LOO; subsets per cell=%d",
                   nrow(S11), nrow(Vs), sum(Vs$loo_all_sig, na.rm=TRUE),
                   uniqueN(S11$subset)),
       outputs=c(basename(out), basename(out2)))
})
