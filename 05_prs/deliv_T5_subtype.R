#!/usr/bin/env Rscript
# ============================================================
# deliv_T5_subtype.R  [ B19 -> T5, preterm birth by subtype: spontaneous vs indicated ]
#
# ------------------------------------------------------------------
# WHY THIS IS A MECHANISM TEST, NOT A SUBGROUP FISHING EXPEDITION
#
# S10 found combined PTB null (pooled OR 1.28 per 10 mmHg, p=0.21). It would be easy to read
# that as "BP does not cause preterm birth" and stop. That reading is wrong, because
# combined PTB pools two outcomes with DIFFERENT CAUSAL STRUCTURES:
#
#   INDICATED (provider-initiated) preterm birth is frequently IATROGENIC WITH RESPECT TO
#     BLOOD PRESSURE -- the clinician delivers early BECAUSE the mother is hypertensive or
#     preeclamptic. If maternal BP has any effect on preterm birth at all, this is the
#     pathway it must travel. The "effect" is partly a clinical decision rule, which is a
#     real causal pathway even though it is mediated by human judgement.
#
#   SPONTANEOUS preterm labour has no such pathway. BP would have to act through
#     placental/vascular mechanisms, and the prior for a large effect is much weaker.
#
# Combining them therefore DILUTES: a real effect on the indicated subtype gets averaged
# against a null on the spontaneous one. The null on combined PTB is exactly what a diluted
# indicated-only effect looks like.
#
# THIS MAKES T5 A SHARPER TEST THAN COMBINED PTB, NOT A WEAKER ONE. It predicts a SPECIFIC
# PATTERN in advance: indicated moves, spontaneous does not. That is falsifiable in a way
# "some effect somewhere" is not --
#   * indicated moves, spontaneous null  -> supports the iatrogenic pathway
#   * BOTH move by similar amounts       -> suspicious. A shared effect on both subtypes is
#                                           more consistent with confounding or with the
#                                           subtype coding being unreliable than with the
#                                           mechanism, and should be reported as such.
#   * neither moves                      -> combined PTB null was a true null, not dilution
#
# COMPETING-RISKS HANDLING. These subtypes compete: a woman delivering spontaneously preterm
# cannot also deliver indicated preterm. Analysing "indicated vs everyone else" would put
# spontaneous preterm births in the comparison group, which is wrong -- they are neither the
# outcome nor a clean control. Each model therefore uses ONLY that subtype as cases and TERM
# births as controls, dropping the other subtype entirely. Ns differ between the two models
# by design and are reported per row.
#
# B01 builds `subtype` as term / spont / indic / unk from SPONT_LABOUR, where unk means PTB
# with SPONT missing. `unk` is EXCLUDED from both models rather than assigned -- assigning it
# would manufacture the very contrast being tested. Its size is reported because if `unk` is
# large the whole analysis is fragile.
# ------------------------------------------------------------------
#   Writes: tables/T5_ptb_subtype.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

NPC_ADJ   <- as.integer(momi_arg("--npc", "5"))
MIN_CASES <- as.integer(momi_arg("--min-cases", "50"))

momi_deliverable("T5_subtype", script="05_prs/deliv_T5_subtype.R",
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

  ## ---- how much of the PTB set is uncodable? if this is large, T5 is fragile ----
  comp <- A2[, .(term=sum(subtype=="term", na.rm=TRUE),
                 spont=sum(subtype=="spont", na.rm=TRUE),
                 indic=sum(subtype=="indic", na.rm=TRUE),
                 unk=sum(subtype=="unk", na.rm=TRUE)), by=cohort][order(cohort)]
  comp[, pct_ptb_uncoded := round(100*unk/(spont+indic+unk), 1)]
  cat("\nPTB subtype composition (unk = preterm with SPONT_LABOUR missing):\n")
  print(comp)
  cat("\nIf pct_ptb_uncoded is high in a cohort, its subtype contrast is driven by whichever\n",
      "subtype happened to be recorded, and the row should not be interpreted.\n", sep="")
  cat(sprintf("\nCohorts reaching the %d-case floor -- indicated: %s | spontaneous: %s\n",
              MIN_CASES,
              paste(comp[indic >= MIN_CASES, cohort], collapse=", "),
              paste(comp[spont >= MIN_CASES, cohort], collapse=", ")))
  cat("THIS COMPOSITION IS ITSELF A RESULT. 265 indicated preterm births at Matlab against\n",
      "18 at Pemba is a statement about obstetric practice across these settings, and it is\n",
      "descriptive fact rather than inference. It also explains why only some sites can\n",
      "contribute an estimate -- the outcome barely occurs in the others.\n", sep="")

  BPDEF <- MOMI_BPDEF_DEFAULT
  SUB   <- c(spont="spontaneous", indic="indicated")

  rows <- list()
  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")) for(sb in names(SUB)){
    sid   <- momi_instrument(coh, tr)
    bpcol <- paste0(substr(tr,1,1), "_", BPDEF)
    zc <- Z[score_id==sid & cohort==coh, .(IID, z)]
    if(!nrow(zc)) next
    d <- merge(A2[cohort==coh], zc, by="IID")
    if(!nrow(d) || !(bpcol %in% names(d))) next

    ## competing risks: cases = this subtype, controls = TERM only, other subtype dropped
    d <- d[subtype %in% c("term", sb)]
    d[, y := as.integer(subtype == sb)]
    ## MINIMUM CASE THRESHOLD, raised 10 -> 50 on 2026-07-20.
    ##
    ## At 10 the table admitted cells with 18, 37 and 44 indicated preterm births. Those
    ## cannot estimate an odds ratio to any useful precision, and including them created a
    ## POOLED estimate (OR 2.19) that was really GAPPS-Bangladesh's number -- it holds 265 of
    ## 448 indicated cases -- averaged with four cohorts contributing noise. The leave-one-out
    ## then appeared to "reveal" fragility, but dropping GAPPS-B gives OR 1.04 [0.39, 2.76],
    ## an interval that COMFORTABLY CONTAINS 2.19. The other cohorts do not contradict the
    ## effect; they are uninformative about it. Presenting a five-cohort pooled row implied a
    ## breadth of evidence that does not exist.
    ##
    ## 50 events is the conventional rule-of-thumb floor for a stable logistic fit with a
    ## handful of covariates. Cohorts below it are reported in the COMPOSITION table (where
    ## their case counts are the informative content) but contribute no estimate.
    if(sum(d$y, na.rm=TRUE) < MIN_CASES) next

    w <- momi_wald(d, bpcol=bpcol, outcol="y", type="bin", cov=covall)
    if(is.null(w)) next
    rows[[paste(coh,tr,sb)]] <- data.table(
      cohort=coh, cohort_display=MOMI_DISPLAY[[coh]], ancestry=MOMI_ANC[[coh]],
      trait=tr, subtype=SUB[[sb]], instrument=sid,
      N=w$n, n_case=sum(d$y, na.rm=TRUE),
      rf_beta=w$bo, rf_se=w$so, rf_p=2*pnorm(-abs(w$bo/w$so)),
      F=w$Ffs, theta_per10=w$theta*10, se_per10=w$se_theta*10)
  }
  if(!length(rows)) return(list(skip=TRUE, reason="no subtype cell had >=10 cases"))
  T5 <- rbindlist(rows)
  T5[, `:=`(OR=exp(theta_per10),
            OR_lo=exp(theta_per10-1.96*se_per10),
            OR_hi=exp(theta_per10+1.96*se_per10))]

  ## ---- pooled per subtype x trait, and the SAS/AFR strata ----
  ## Pooling now runs only over cohorts that CLEARED the case threshold, so a pooled row can
  ## no longer be one cohort's estimate diluted by four uninformative ones. k is reported on
  ## every row and must be read alongside the estimate: a "pooled" OR with k=1 is not pooled,
  ## and the console labels it as a single-site result rather than printing it as a meta-analysis.
  pool <- list()
  for(tr in c("SBP","DBP")) for(sb in unname(SUB)) for(st in c("all","SAS","AFR")){
    s <- T5[trait==tr & subtype==sb & is.finite(theta_per10) & se_per10>0]
    if(st != "all") s <- s[ancestry==st]
    if(nrow(s) < 1) next
    if(nrow(s) == 1){
      pool[[paste(tr,sb,st)]] <- data.table(
        trait=tr, subtype=sb, stratum=st, k=1L, N=s$N[1], cases=s$n_case[1],
        theta_per10=s$theta_per10[1], se=s$se_per10[1],
        p=2*pnorm(-abs(s$theta_per10[1]/s$se_per10[1])),
        OR=exp(s$theta_per10[1]),
        OR_lo=exp(s$theta_per10[1]-1.96*s$se_per10[1]),
        OR_hi=exp(s$theta_per10[1]+1.96*s$se_per10[1]),
        I2=NA_real_, F_min=s$F[1], F_max=s$F[1],
        n_neg=as.integer(s$theta_per10[1]<0), n_pos=as.integer(s$theta_per10[1]>0),
        single_site=s$cohort[1])
      next
    }
    fe <- momi_meta_iv(s$theta_per10, s$se_per10)
    w  <- 1/s$se_per10^2
    Q  <- sum(w*(s$theta_per10-fe$b)^2); dfQ <- nrow(s)-1
    pool[[paste(tr,sb,st)]] <- data.table(
      trait=tr, subtype=sb, stratum=st, k=fe$k, N=sum(s$N), cases=sum(s$n_case),
      theta_per10=fe$b, se=fe$se, p=fe$p,
      OR=exp(fe$b), OR_lo=exp(fe$b-1.96*fe$se), OR_hi=exp(fe$b+1.96*fe$se),
      I2=max(0, 100*(Q-dfQ)/Q), F_min=min(s$F), F_max=max(s$F),
      n_neg=sum(s$theta_per10<0), n_pos=sum(s$theta_per10>0))
  }
  PL <- rbindlist(pool, fill=TRUE)

  out  <- momi_write_table(T5, "T5_ptb_subtype", P)
  out2 <- momi_write_table(PL, "T5b_ptb_subtype_pooled", P)
  momi_save_intermediate(T5, "ptb_subtype", P)

  ## ---------------- console: the pre-specified contrast ----------------
  cat("\n=== OR per 10 mmHg, by subtype and stratum (k=1 rows are SINGLE-SITE, not pooled) ===\n")
  if(nrow(PL)) print(PL[, .(trait, subtype, stratum, k, cases,
                            OR=round(OR,3), lo=round(OR_lo,3), hi=round(OR_hi,3),
                            p=signif(p,3), I2=round(I2), Fmin=round(F_min),
                            site = if("single_site" %in% names(PL)) single_site else NA_character_)])

  cat("\n=== THE TEST: indicated vs spontaneous, same trait and stratum ===\n")
  for(st in c("all","SAS")) for(tr in c("SBP","DBP")){
    a <- PL[trait==tr & stratum==st & subtype=="indicated"]
    b <- PL[trait==tr & stratum==st & subtype=="spontaneous"]
    if(!nrow(a) || !nrow(b)) next
    ## difference of two independent log-ORs
    dd <- a$theta_per10 - b$theta_per10
    sd_ <- sqrt(a$se^2 + b$se^2)
    cat(sprintf("  %s %-4s  indicated OR %.3f (p=%.3g) vs spontaneous OR %.3f (p=%.3g); diff p=%.3g\n",
                st, tr, a$OR, a$p, b$OR, b$p, 2*pnorm(-abs(dd/sd_))))
  }
  cat("\nSupports the iatrogenic pathway ONLY if indicated moves and spontaneous does not.\n",
      "If BOTH move similarly, that is evidence AGAINST the mechanism -- report it as such\n",
      "rather than as 'an effect on preterm birth'.\n", sep="")

  list(n=nrow(T5),
       key=sprintf("cells=%d; pooled rows=%d; cases indicated=%d spontaneous=%d",
                   nrow(T5), nrow(PL),
                   sum(T5[subtype=="indicated"]$n_case), sum(T5[subtype=="spontaneous"]$n_case)),
       outputs=c(basename(out), basename(out2), "ptb_subtype.rds"))
})
