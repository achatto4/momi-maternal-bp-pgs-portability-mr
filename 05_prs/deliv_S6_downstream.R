#!/usr/bin/env Rscript
# ============================================================
# deliv_S6_downstream.R  [ B20 -> Supp Table S6, does the BP definition change the answer? ]
#
# ------------------------------------------------------------------
# THE QUESTION THIS ANSWERS, AND WHY IT IS NOT THE SAME AS S7
#
# S7 showed the FIRST STAGE is insensitive to how blood pressure is summarised: across ten
# cohort x trait cells, GA-standardising moves incremental R2 by at most 0.20 percentage
# points. That is reassuring about the instrument but says nothing about the causal estimates,
# because the Wald ratio divides by the first stage -- a definition could leave R2 untouched
# while shifting the exposure's MEANING (a late-pregnancy reading is not the same quantity as
# a first-trimester one) and so move theta.
#
# This module therefore re-runs the ENTIRE MR PANEL under every BP definition in the grid and
# asks whether the paper's conclusions survive. The exposure choice (`resid`) was made on
# conceptual grounds; S6 is what makes that choice defensible rather than arbitrary.
#
# WHAT WOULD FALSIFY THE PAPER'S HEADLINE. SBP -> birthweight is reported at roughly -80 to
# -90 g per 10 mmHg. If that estimate flips sign, or collapses toward zero, under a
# reasonable alternative definition, then the finding is an artefact of one summarisation
# choice and must be reported as such. The console prints the spread across definitions for
# exactly this reason.
#
# A PREDICTION WORTH RECORDING BEFORE LOOKING. The definitions differ enormously in what they
# measure: `last` is a single late reading (and S16 shows it is where the two Bangladeshi
# cohorts diverge MOST, SMD -0.50), `lt20` is early pregnancy in the few cohorts that have it,
# `mean` and `resid` average everything. If the causal estimate is stable across all of them,
# that is genuinely reassuring. If it moves systematically with gestational timing -- larger
# for late definitions, smaller for early ones -- that is not noise, it is evidence that the
# effect on birthweight accumulates through pregnancy, which would be a finding rather than a
# nuisance and should be reported as one.
#
# DEFINITIONS WITH TINY N ARE EXCLUDED, NOT SHOWN AS WIDE INTERVALS. `tri1` has n=6 in
# AMANHI-Bangladesh and `lt20` n=195; including them would fill the table with uninterpretable
# rows. A cohort must clear MIN_N for that definition to contribute -- the same discipline
# applied to the PTB subtypes, and for the same reason: a wide interval from an absent
# exposure reads as a weak effect when it is really no measurement.
# ------------------------------------------------------------------
#   Writes: tables/S6_definition_sensitivity.tsv, figures/S6_definition_spread.{png,pdf}
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

NPC_ADJ <- as.integer(momi_arg("--npc", "5"))
MIN_N   <- as.integer(momi_arg("--min-n", "300"))   # per cohort, per definition

momi_deliverable("S6_downstream", script="05_prs/deliv_S6_downstream.R",
                 inputs="analytic_mothers;prs_z;pcs", P=P, stop_on_error=FALSE,
                 body=function(ctx){

  A  <- momi_read_intermediate("analytic_mothers", P)
  Z  <- momi_read_intermediate("prs_z", P)
  PC <- tryCatch(momi_read_intermediate("pcs", P), error=function(e) NULL)
  if(is.null(PC) || !nrow(PC)) return(list(skip=TRUE, reason="no pcs.rds -- run B06b first"))
  pccols <- intersect(paste0("PC", seq_len(NPC_ADJ)), names(PC))
  covall <- c("AGE", pccols)

  A[, hasBP := is.finite(S_mean) | is.finite(D_mean)]
  A2 <- A[genotyped==1 & hasBP==TRUE & !is.na(PTB)]
  A2 <- merge(A2, PC[, c("IID","cohort",pccols), with=FALSE], by=c("IID","cohort"), all.x=TRUE)

  OUTS <- c(MOMI_OUTCOMES_MAIN, "BWT")
  rows <- list()

  for(dfn in MOMI_BPDEFS) for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    bpcol <- paste0(substr(tr,1,1), "_", dfn)
    if(!(bpcol %in% names(A2))) next
    sid <- momi_instrument(coh, tr)
    zc  <- Z[score_id==sid & cohort==coh, .(IID, z)]
    if(!nrow(zc)) next
    d0 <- merge(A2[cohort==coh], zc, by="IID")
    ## the definition must actually be measured in this cohort (see header)
    if(sum(is.finite(d0[[bpcol]])) < MIN_N) next

    for(oc in OUTS){
      if(!(oc %in% names(d0))) next
      d  <- if(oc %in% MOMI_LIVEBIRTH_OUTCOMES) d0[livebirth==1] else d0
      ty <- if(oc=="BWT") "lin" else "bin"
      w  <- momi_wald(d, bpcol=bpcol, outcol=oc, type=ty, cov=covall)
      if(is.null(w)) next
      rows[[paste(dfn,coh,tr,oc)]] <- data.table(
        definition=dfn, cohort=coh, ancestry=MOMI_ANC[[coh]], trait=tr, outcome=oc, type=ty,
        n=w$n, F=w$Ffs, theta_per10=w$theta*10, se_per10=w$se_theta*10,
        rf_p=2*pnorm(-abs(w$bo/w$so)))
    }
  }
  if(!length(rows)) return(list(skip=TRUE, reason="no definition x cohort cell cleared MIN_N"))
  D <- rbindlist(rows)

  ## ---- pool per definition x trait x outcome, SAS and all ----
  pool <- list()
  for(dfn in unique(D$definition)) for(tr in c("SBP","DBP")) for(oc in OUTS)
    for(st in c("SAS","all")){
      s <- D[definition==dfn & trait==tr & outcome==oc & is.finite(theta_per10) & se_per10>0]
      if(st!="all") s <- s[ancestry==st]
      if(nrow(s) < 2) next
      fe <- momi_meta_iv(s$theta_per10, s$se_per10)
      pool[[paste(dfn,tr,oc,st)]] <- data.table(
        definition=dfn, trait=tr, outcome=oc, type=s$type[1], stratum=st,
        k=fe$k, N=sum(s$n), theta_per10=fe$b, se=fe$se, p=fe$p,
        lo=fe$b-1.96*fe$se, hi=fe$b+1.96*fe$se,
        primary = dfn == MOMI_BPDEF_DEFAULT)
    }
  S6 <- rbindlist(pool, fill=TRUE)
  if(!nrow(S6)) return(list(skip=TRUE, reason="nothing poolable"))
  S6[, effect := fifelse(type=="bin", exp(theta_per10), theta_per10)]

  out <- momi_write_table(S6, "S6_definition_sensitivity", P)

  ## ---- figure: the spread of the causal estimate across definitions ----
  pl <- S6[stratum=="SAS"]
  g <- ggplot(pl, aes(x=theta_per10, y=definition, colour=primary)) +
    geom_vline(xintercept=0, linetype="dashed", linewidth=.3, colour="grey40") +
    geom_errorbarh(aes(xmin=lo, xmax=hi), height=.15, linewidth=.35) +
    geom_point(size=1.7) +
    facet_grid(trait ~ outcome, scales="free_x") +
    scale_colour_manual(values=c(`FALSE`="grey45", `TRUE`="#C44E52"), guide="none") +
    labs(x="effect per 10 mmHg (log-OR for binary, grams for BWT)", y="BP definition",
         title="Causal estimates under every BP definition (South Asian cohorts pooled)",
         subtitle=sprintf("Primary definition (%s) in red. Stability across definitions is what makes the exposure choice defensible.",
                          MOMI_BPDEF_DEFAULT)) +
    theme_minimal(base_size=8) +
    theme(panel.grid.minor=element_blank(),
          plot.subtitle=element_text(size=6.5, colour="grey35"))
  fig <- momi_save_fig(g, "S6_definition_spread", width=9, height=5.5, P=P)

  ## ---- console ----
  cat(sprintf("\ndefinitions retained (>= %d measured per cohort): %s\n",
              MIN_N, paste(sort(unique(D$definition)), collapse=", ")))
  cat("\n=== does the headline survive every definition? SBP -> BWT, SAS pooled ===\n")
  hl <- S6[trait=="SBP" & outcome=="BWT" & stratum=="SAS",
           .(definition, k, theta=round(theta_per10,1), lo=round(lo,1), hi=round(hi,1),
             p=signif(p,3), primary)]
  if(nrow(hl)){
    print(hl[order(theta)])
    cat(sprintf("\n  range across definitions: %.1f to %.1f g per 10 mmHg; all same sign: %s\n",
                min(hl$theta), max(hl$theta), all(sign(hl$theta)==sign(hl$theta[1]))))
  }

  cat("\n=== spread across definitions, every outcome (SAS pooled) ===\n")
  sp <- S6[stratum=="SAS", .(k_defs=.N,
                             min=round(min(theta_per10),3), max=round(max(theta_per10),3),
                             primary=round(theta_per10[primary][1],3),
                             same_sign=all(sign(theta_per10)==sign(theta_per10[1])),
                             n_sig=sum(p<0.05)), by=.(trait, outcome)]
  print(sp)
  cat("\nIf `same_sign` is TRUE and the primary sits inside the range, the conclusion does not\n",
      "depend on the summarisation choice. If the estimate moves MONOTONICALLY with gestational\n",
      "timing (lt20 -> tri1 -> tri2 -> tri3 -> last), that is a finding about when in pregnancy\n",
      "BP matters, not a robustness failure -- report it rather than averaging it away.\n", sep="")

  list(n=nrow(S6),
       key=sprintf("defs=%d; pooled rows=%d; SBP->BWT sign stable=%s",
                   uniqueN(S6$definition), nrow(S6),
                   as.character(all(sign(S6[trait=="SBP" & outcome=="BWT" & stratum=="SAS"]$theta_per10) < 0))),
       outputs=c(basename(out), basename(fig)))
})
