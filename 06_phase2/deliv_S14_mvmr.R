#!/usr/bin/env Rscript
# ============================================================
# deliv_S14_mvmr.R  [ B31 -> Supp Table S14, is the BP effect independent of adiposity? ]
#
# ------------------------------------------------------------------
# THE OBJECTION. Part II reports that genetically predicted maternal SBP lowers birthweight
# by roughly 80-90 g per 10 mmHg. Blood pressure and body mass share substantial biology and
# substantial genetics, so a reviewer will reasonably ask whether the BP score is standing in
# for adiposity -- i.e. whether we have measured a BMI effect wearing a blood-pressure label.
#
# THE TEST. Multivariable MR. Put BOTH standardised scores in the model at once:
#     first stage A:   BP  ~ zBP + zBMI + covariates
#     first stage B:   BMI ~ zBP + zBMI + covariates
#     reduced form:    Y   ~ zBP + zBMI + covariates
# The coefficient on zBP in the reduced form is now the association of BP-raising genotype
# with the outcome HOLDING ADIPOSITY GENOTYPE FIXED. If the BP effect survives, it is not
# adiposity. If it collapses toward zero once zBMI is included, it largely was.
#
# WHAT WOULD MAKE THIS TEST UNINFORMATIVE, AND IS CHECKED BEFORE ANYTHING IS REPORTED:
#   1. If zBP and zBMI are strongly correlated in our data, the two cannot be separated and
#      the conditional estimates will simply be unstable -- wide, and sensitive to sample.
#      The correlation and the variance inflation factor are printed FIRST for that reason.
#   2. If the BMI score does not predict measured BMI in these cohorts, it is not a valid
#      instrument for adiposity here and conditioning on it accomplishes nothing. Given
#      Part I found EUR scores transfer poorly to these populations -- and PGS000027 is a
#      EUR BMI score -- this is a live possibility, not a formality. Its first-stage R2 on
#      measured BMI is reported per cohort and must be read before the MVMR estimates.
#
# HONEST FRAMING IF THE BMI INSTRUMENT IS WEAK. A weak zBMI cannot adjust away a confounder
# it does not capture, so "the BP effect survived MVMR" would be a hollow reassurance. In
# that case the correct statement is that we CANNOT rule out adiposity, not that we have.
# The console says so explicitly rather than leaving it to interpretation.
#
# COVARIATES: age + within-cohort PCs, matching B22 exactly (Burgess; see
# ref/methods_citations.tsv). Same set in every stage, or the ratio is incoherent.
# ------------------------------------------------------------------
#   Writes: tables/S14_mvmr.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

NPC_ADJ <- as.integer(momi_arg("--npc", "5"))
BMI_PGS <- momi_arg("--bmi-pgs", "PGS000027")

momi_deliverable("S14_mvmr", script="06_phase2/deliv_S14_mvmr.R",
                 inputs="analytic_mothers;prs_z;pcs", P=P, stop_on_error=FALSE,
                 body=function(ctx){

  A  <- momi_read_intermediate("analytic_mothers", P)
  Z  <- momi_read_intermediate("prs_z", P)
  PC <- tryCatch(momi_read_intermediate("pcs", P), error=function(e) NULL)
  if(is.null(PC) || !nrow(PC)) return(list(skip=TRUE, reason="no pcs.rds -- run B06b first"))
  if(!nrow(Z[score_id==BMI_PGS]))
    return(list(skip=TRUE,
                reason=sprintf("%s not in prs_z -- fetch it (bin/fetch_pgs.sh) and score it first", BMI_PGS)))

  pccols <- intersect(paste0("PC", seq_len(NPC_ADJ)), names(PC))
  covall <- c("AGE", pccols)
  A[, hasBP := is.finite(S_mean) | is.finite(D_mean)]
  A2 <- A[genotyped==1 & hasBP==TRUE & !is.na(PTB)]
  A2 <- merge(A2, PC[, c("IID","cohort",pccols), with=FALSE], by=c("IID","cohort"), all.x=TRUE)

  BPDEF <- MOMI_BPDEF_DEFAULT
  OUTS  <- c(MOMI_OUTCOMES_MAIN, "BWT")
  rows <- list(); diag <- list()

  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    sid   <- momi_instrument(coh, tr)
    bpcol <- paste0(substr(tr,1,1), "_", BPDEF)
    zb <- Z[score_id==sid     & cohort==coh, .(IID, zBP=z)]
    zm <- Z[score_id==BMI_PGS & cohort==coh, .(IID, zBMI=z)]
    if(!nrow(zb) || !nrow(zm)) next
    d <- merge(merge(A2[cohort==coh], zb, by="IID"), zm, by="IID")
    if(!nrow(d) || !(bpcol %in% names(d))) next

    ## ---- diagnostics that decide whether the MVMR means anything (see header) ----
    dd <- d[is.finite(zBP) & is.finite(zBMI)]
    r_scores <- if(nrow(dd) > 10) cor(dd$zBP, dd$zBMI, use="complete.obs") else NA_real_
    bmi_fs <- if("BMI" %in% names(d)) momi_transfer(d$zBMI, d$BMI, cov=d[, .(AGE)]) else NULL
    diag[[paste(coh,tr)]] <- data.table(
      cohort=coh, trait=tr, n=nrow(dd),
      cor_zBP_zBMI = round(r_scores, 4),
      vif = round(1/(1 - r_scores^2), 3),
      bmi_score_R2pct_on_BMI = if(is.null(bmi_fs)) NA_real_ else round(100*bmi_fs$incR2, 3),
      bmi_score_F            = if(is.null(bmi_fs)) NA_real_ else round(bmi_fs$F, 1))

    for(oc in OUTS){
      if(!(oc %in% names(d))) next
      x  <- if(oc %in% MOMI_LIVEBIRTH_OUTCOMES) d[livebirth==1] else d
      ty <- if(oc=="BWT") "lin" else "bin"
      keep <- c("zBP","zBMI", bpcol, "BMI", oc, covall)
      x <- x[, intersect(keep, names(x)), with=FALSE]
      setnames(x, c(bpcol, oc), c("bp","y"))
      x <- x[is.finite(zBP) & is.finite(zBMI) & is.finite(y)]
      x <- x[complete.cases(x[, c("zBP","zBMI","y",..covall)])]
      if(nrow(x) < 100) next
      if(ty=="bin" && sum(x$y, na.rm=TRUE) < 25) next

      crhs <- paste(covall, collapse="+")
      ## univariable reduced form (what B22 reports), then the multivariable one
      f_uni <- as.formula(paste("y ~ zBP +", crhs))
      f_mv  <- as.formula(paste("y ~ zBP + zBMI +", crhs))
      m_uni <- tryCatch(if(ty=="bin") glm(f_uni, binomial, x) else lm(f_uni, x), error=function(e) NULL)
      m_mv  <- tryCatch(if(ty=="bin") glm(f_mv,  binomial, x) else lm(f_mv,  x), error=function(e) NULL)
      if(is.null(m_uni) || is.null(m_mv)) next
      cu <- summary(m_uni)$coef; cm <- summary(m_mv)$coef
      if(!("zBP" %in% rownames(cu)) || !("zBP" %in% rownames(cm))) next

      ## first stage for the BP exposure, conditional on zBMI -- the MVMR denominator
      fs <- tryCatch(summary(lm(as.formula(paste("bp ~ zBP + zBMI +", crhs)),
                                x[is.finite(bp)]))$coef, error=function(e) NULL)
      be <- if(is.null(fs) || !("zBP" %in% rownames(fs))) NA_real_ else fs["zBP",1]
      se_be <- if(is.null(fs) || !("zBP" %in% rownames(fs))) NA_real_ else fs["zBP",2]

      rows[[paste(coh,tr,oc)]] <- data.table(
        cohort=coh, cohort_display=MOMI_DISPLAY[[coh]], ancestry=MOMI_ANC[[coh]],
        trait=tr, outcome=oc, type=ty, N=nrow(x),
        rf_zBP_uni   = cu["zBP",1], rf_zBP_uni_se = cu["zBP",2], rf_zBP_uni_p = cu["zBP",4],
        rf_zBP_mv    = cm["zBP",1], rf_zBP_mv_se  = cm["zBP",2], rf_zBP_mv_p  = cm["zBP",4],
        rf_zBMI_mv   = cm["zBMI",1], rf_zBMI_mv_se= cm["zBMI",2], rf_zBMI_mv_p = cm["zBMI",4],
        pct_change_zBP = round(100*(cm["zBP",1]-cu["zBP",1])/abs(cu["zBP",1]), 1),
        fs_beta_cond = be, fs_F_cond = if(is.finite(be)) (be/se_be)^2 else NA_real_,
        theta_mv_per10 = if(is.finite(be) && be!=0) 10*cm["zBP",1]/be else NA_real_)
    }
  }
  if(!length(rows)) return(list(skip=TRUE, reason="no estimable MVMR cell"))
  S14 <- rbindlist(rows); DG <- rbindlist(diag, fill=TRUE)
  out  <- momi_write_table(S14, "S14_mvmr", P)
  out2 <- momi_write_table(DG,  "S14b_mvmr_diagnostics", P)

  ## ---------------- console: diagnostics FIRST ----------------
  cat("\n=== is this test informative? read before the estimates ===\n")
  print(DG)
  wk <- DG[is.finite(bmi_score_F) & bmi_score_F < 10]
  cat(sprintf("\nBMI-score first stage on measured BMI: median R2 = %.3f%%, cohorts with F<10: %d of %d\n",
              median(DG$bmi_score_R2pct_on_BMI, na.rm=TRUE), nrow(wk), nrow(DG)))
  if(nrow(wk) >= nrow(DG)/2){
    cat("\n*** THE BMI INSTRUMENT IS WEAK IN MOST COHORTS. ***\n",
        "Conditioning on a score that does not capture adiposity cannot adjust adiposity away.\n",
        "If the BP effect survives below, that is NOT evidence it is independent of BMI -- the\n",
        "correct statement is that this analysis cannot rule adiposity out. Report it that way.\n",
        "This is the expected result if PGS000027 (a EUR BMI score) fails to transfer here,\n",
        "which is exactly what Part I documents for the EUR BP scores.\n", sep="")
  }
  cat(sprintf("\nzBP-zBMI correlation: median %.3f, max |r| %.3f (VIF max %.2f)\n",
              median(DG$cor_zBP_zBMI, na.rm=TRUE), max(abs(DG$cor_zBP_zBMI), na.rm=TRUE),
              max(DG$vif, na.rm=TRUE)))

  cat("\n=== BP effect before vs after conditioning on adiposity genotype ===\n")
  key <- S14[outcome %in% c("BWT","LBW"),
             .(cohort, trait, outcome, N,
               uni=round(rf_zBP_uni,4), mv=round(rf_zBP_mv,4),
               pct_change=pct_change_zBP, mv_p=signif(rf_zBP_mv_p,3),
               zBMI_p=signif(rf_zBMI_mv_p,3))]
  if(nrow(key)) print(key[order(outcome, trait, cohort)])
  cat("\nA small pct_change means the BP association does not run through adiposity genetics.\n",
      "A large one means it substantially does. Judge this per cohort -- and only where the\n",
      "BMI instrument above is actually strong enough for the question to be answerable.\n", sep="")

  list(n=nrow(S14),
       key=sprintf("cells=%d; BMI-score median R2 on BMI=%.3f%%; weak-BMI cohorts=%d/%d; median |pct change| in zBP=%.1f%%",
                   nrow(S14), median(DG$bmi_score_R2pct_on_BMI, na.rm=TRUE),
                   nrow(wk), nrow(DG), median(abs(S14$pct_change_zBP), na.rm=TRUE)),
       outputs=c(basename(out), basename(out2)))
})
