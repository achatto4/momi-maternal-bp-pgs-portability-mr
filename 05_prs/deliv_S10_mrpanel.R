#!/usr/bin/env Rscript
# ============================================================
# deliv_S10_mrpanel.R  [ B22 -> S10, the Mendelian randomization panel ]
#
# The Part-II core analysis: does GENETICALLY PREDICTED maternal blood pressure affect
# perinatal outcomes? One row per cohort x trait x outcome, plus pooled rows.
#
# ------------------------------------------------------------------
# WHAT IS ESTIMATED, AND WHICH NUMBER IS THE TEST
#
# For a single instrument the ratio (Wald) estimator is
#       theta = beta_reduced_form / beta_first_stage
# where the reduced form is (outcome ~ PRS) and the first stage is (BP ~ PRS).
#
# The REDUCED FORM p-value is the honest test of the causal null. Under the null of no
# causal effect of BP on the outcome, the PRS is unassociated with the outcome, so
# beta_reduced_form = 0 -- and that test does not involve the first stage at all. The Wald
# ratio then only RESCALES that same evidence into mmHg units. Dividing by a small, noisily
# estimated first-stage coefficient inflates both the point estimate and its SE, which is why
# a weak instrument produces spectacular-looking effect sizes that mean nothing. S13 found
# 0 of 40 cells reach 80% power, so that failure mode is live here, not hypothetical.
#
# CONSEQUENCE FOR THE PAPER: p_reduced_form is the column to read for significance, and
# theta is the column to read for magnitude, and the two must be reported together. The
# script therefore never reports a Wald estimate without its reduced-form p and its F.
# ------------------------------------------------------------------
# COVARIATES: AGE + WITHIN-COHORT PCs. NOTHING ELSE.
#
# This is deliberately NOT the observational confounder set (MOMI_CONF_CORE = AGE, GRAV,
# EDU, BMI). Burgess et al. (Wellcome Open Res 2019;4:186, ref/methods_citations.tsv)
# recommend restricting MR covariates to age, sex, genomic PCs and technical covariates,
# because adjusting further risks conditioning on a MEDIATOR (BMI plausibly sits on the
# path from BP-raising genotype to outcome) or inducing COLLIDER bias. Gravidity and
# education are post-randomization or non-genetic and have no business in a first stage.
#
# The same covariate set is used in the reduced form AND the first stage. Burgess is
# explicit that the two 2SLS stages take the same covariates; mismatched stages make the
# ratio incoherent, since numerator and denominator would condition on different things.
#
# WHY WITHIN-COHORT PCs SPECIFICALLY (from B06b): the confounding an MR analysis must
# control is population structure INSIDE the analysis sample. Projected 1000G PCs describe
# where a cohort sits globally, which is a DESCRIPTIVE question (SF1), not the adjustment.
# ------------------------------------------------------------------
# POOLING: TWO METHODS, REPORTED SIDE BY SIDE, BECAUSE THEY DISAGREE ABOUT WHAT IS ASSUMED
#
#   pooled_meta   -- inverse-variance meta-analysis of the five per-cohort Wald ratios,
#                    with Cochran's Q, I2 and a random-effects (DerSimonian-Laird) companion.
#                    Assumes nothing about cross-cohort commensurability. PREFERRED.
#   pooled_stack  -- one regression on the stacked data with cohort fixed effects.
#                    More efficient IF the effect is homogeneous, which S9b's I2 of 70-92%
#                    on the observational side gives us reason to doubt.
#
# A TRAP IN THE STACKED MODEL. Within-cohort PCs come from FIVE SEPARATE PCAs, so "PC1" in
# GAPPS-Zambia and "PC1" in AMANHI-Pakistan are different axes that merely share a name.
# Entering PC1 as a single pooled term would impose a common coefficient on incommensurable
# coordinates. The stacked model therefore interacts every PC with cohort (PC1:cohort etc),
# which is algebraically the same as allowing each cohort its own PC adjustment. The meta
# route avoids the issue entirely, which is one reason it is preferred.
#
# A SECOND TRAP, NOT FIXABLE HERE. The instrument is not the same score in every cohort:
# T3 locked SAS (PGS004830/PGS004758) for GAPPS-Bangladesh and EUR elsewhere. Pooling
# therefore averages effects estimated with DIFFERENT instruments. That is defensible for a
# ratio estimate -- theta is on the mmHg scale regardless of which score got you there, and
# the ratio divides out instrument strength -- but it is an assumption and is flagged in the
# output as `instrument_mixed`.
# ------------------------------------------------------------------
# SAMPLE. Part II (genotyped & has antenatal BP & non-missing PTB), per F1/B07. Outcomes
# derived from birthweight (LBW, SGA, BWT) are additionally restricted to LIVE BIRTHS
# (MOMI_LIVEBIRTH_OUTCOMES) -- a stillbirth has a birthweight but conditioning an analysis
# of fetal growth on survival is a different question and a collider.
#
# PE and twins are NOT excluded. They are descriptors in F1, and preeclampsia is a MEDIATOR
# on the BP-to-PTB path (Fig 1B DAG), so excluding it would remove part of the effect we are
# trying to estimate. S11/B23 tests sensitivity to that choice.
#
#   Writes: tables/S10_mr_panel.tsv, tables/S10b_mr_pooled.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

NPC_ADJ <- as.integer(momi_arg("--npc", "5"))   # PCs used for adjustment; 5 is conventional

momi_deliverable("S10_mrpanel", script="05_prs/deliv_S10_mrpanel.R",
                 inputs="analytic_mothers;prs_z;pcs;power_grid", P=P,
                 stop_on_error=FALSE, body=function(ctx){

  A <- momi_read_intermediate("analytic_mothers", P)
  Z <- momi_read_intermediate("prs_z", P)
  PC <- tryCatch(momi_read_intermediate("pcs", P), error=function(e) NULL)
  if(is.null(PC) || !nrow(PC))
    return(list(skip=TRUE, reason="no pcs.rds -- run B06b first; MR must not be run unadjusted for structure"))

  pccols <- intersect(paste0("PC", seq_len(NPC_ADJ)), names(PC))
  if(!length(pccols))
    return(list(skip=TRUE, reason="pcs.rds has no within-cohort PC columns"))

  ## ---- analysis sample ----
  ## hasBP is NOT a column of analytic_mothers -- F1/B07 and T1/B08 each derive it locally as
  ## `is.finite(S_mean) | is.finite(D_mean)`. Repeated verbatim here so Part II is the SAME
  ## 13,212 women F1 reports; any divergence would make the flow diagram wrong about the MR
  ## sample. (Three copies of one definition is a latent inconsistency -- it belongs in B01
  ## as a column. Noted, not fixed here, because moving it re-locks B01.)
  A[, hasBP := is.finite(S_mean) | is.finite(D_mean)]
  ## C2 (2026-07-26): do NOT gate the whole analysis on non-missing preterm birth. That dropped
  ## women whose PTB was missing from the LBW/SGA/birth-weight analyses even though PTB is not
  ## their outcome, biasing those samples and desynchronising them from S13's power sample. The
  ## sample is now defined per outcome -- momi_wald() filters to non-missing values of the
  ## outcome actually being analysed, so each cell uses its own complete-case sample.
  A2 <- A[genotyped==1 & hasBP==TRUE]

  ## NB the exposure is `resid`, which additionally requires a valid GA at the visit, so some
  ## Part II women have S_mean but no S_resid and drop out of the first stage. That is why
  ## the per-cell N below can sit under 13,212, and why N is reported on every row.
  A2 <- merge(A2, PC[, c("IID","cohort",pccols), with=FALSE], by=c("IID","cohort"), all.x=TRUE)

  ## Mothers with no PCs cannot enter an adjusted model. Report the loss rather than let it
  ## happen silently inside complete.cases() -- if a cohort loses many here, its MR row is
  ## not comparable to the others and that must be visible.
  covall <- c("AGE", pccols)
  A2[, has_pc := Reduce(`&`, lapply(pccols, function(p) is.finite(get(p))))]
  loss <- A2[, .(N_partII=.N, N_with_PC=sum(has_pc),
                 pct_lost=round(100*(1-sum(has_pc)/.N),1)), by=cohort][order(cohort)]
  cat("\nPC coverage within Part II (mothers dropped here are dropped from every MR cell):\n")
  print(loss)

  BPDEF <- MOMI_BPDEF_DEFAULT
  OUTS  <- c(MOMI_OUTCOMES_MAIN, "BWT")

  ## ---------------- per-cohort estimates ----------------
  rows <- list()
  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    sid   <- momi_instrument(coh, tr)
    bpcol <- paste0(substr(tr,1,1), "_", BPDEF)
    zc <- Z[score_id==sid & cohort==coh, .(IID, z)]
    if(!nrow(zc)) next
    d0 <- merge(A2[cohort==coh], zc, by="IID")
    if(!nrow(d0) || !(bpcol %in% names(d0))) next

    for(oc in OUTS){
      if(!(oc %in% names(d0))) next
      d <- if(oc %in% MOMI_LIVEBIRTH_OUTCOMES) d0[livebirth==1] else d0
      ty <- if(oc=="BWT") "lin" else "bin"
      w  <- momi_wald(d, bpcol=bpcol, outcol=oc, type=ty, cov=covall)
      if(is.null(w)) next
      ncase <- if(ty=="bin") sum(d[[oc]], na.rm=TRUE) else NA_integer_
      rows[[paste(coh,tr,oc)]] <- data.table(
        cohort=coh, cohort_display=MOMI_DISPLAY[[coh]], ancestry=MOMI_ANC[[coh]],
        trait=tr, exposure_def=BPDEF, instrument=sid, outcome=oc, type=ty,
        N=w$n, n_case=ncase,
        ## reduced form: the TEST. beta is per SD of PRS, i.e. NOT on the mmHg scale.
        rf_beta=w$bo, rf_se=w$so, rf_p=2*pnorm(-abs(w$bo/w$so)),
        ## first stage: instrument strength in this exact sample
        fs_beta=w$be, fs_se=w$see, F=w$Ffs,
        ## Wald ratio: the same evidence rescaled to per-1-mmHg, then reported per 10 mmHg
        theta_per10=w$theta*10, se_per10=w$se_theta*10)
    }
  }
  if(!length(rows)) return(list(skip=TRUE, reason="no estimable cohort x trait x outcome cells"))
  MR <- rbindlist(rows)

  ## effect scale: OR per 10 mmHg for binary, grams per 10 mmHg for BWT
  MR[, `:=`(est_lo = theta_per10 - 1.96*se_per10,
            est_hi = theta_per10 + 1.96*se_per10)]
  MR[type=="bin", `:=`(OR=exp(theta_per10), OR_lo=exp(est_lo), OR_hi=exp(est_hi))]
  MR[, weak := F < 10]   # conventional threshold; see S13 for the power consequence

  ## ---------------- pooling ----------------
  pool <- list()
  for(tr in c("SBP","DBP")) for(oc in OUTS){
    s <- MR[trait==tr & outcome==oc & is.finite(theta_per10) & is.finite(se_per10) & se_per10>0]
    if(nrow(s) < 2) next

    ## (a) fixed-effect inverse-variance meta over per-cohort Wald ratios
    fe <- momi_meta_iv(s$theta_per10, s$se_per10)

    ## (a2) RATIO-OF-POOLED-COEFFICIENTS. Added 2026-07-20 after methodological review.
    ##
    ## (a) meta-analyses the five per-cohort Wald RATIOS. Burgess, Small & Thompson
    ## (Stat Methods Med Res 2017;26:2333) advise against exactly this: "a study-level
    ## meta-analysis can accentuate weak instrument bias, and so it is preferable to combine
    ## either individual-level or summarized data from each study."
    ##
    ## The mechanism is specific and worth stating: within a cohort, when the first-stage
    ## coefficient comes out large by chance, the Wald ratio is BOTH pulled toward the
    ## confounded observational association AND given a smaller standard error. Inverse-
    ## variance weighting therefore systematically upweights the cohorts that are most biased.
    ## Burgess & Thompson (IJE 2011;40:755) put it directly: "when the observed F-statistic is
    ## larger than expected in a particular study, the causal estimate is more biased towards
    ## the observational association and its standard error is smaller."
    ##
    ## The recommended alternative pools the two REGRESSIONS separately and forms the ratio
    ## once, at the end -- so no per-cohort ratio is ever weighted by its own noise.
    ##
    ## BOTH ARE REPORTED. This is a check, not a replacement: our first stages are strong
    ## (F 17-250, median ~49) and the distortion scales roughly as 1/F, so the two should
    ## agree closely. If they DISAGREE materially, the ratio-of-pooled version is the one to
    ## trust and the divergence itself needs explaining.
    rfp <- momi_meta_iv(s$rf_beta, s$rf_se)      # pooled instrument -> outcome
    fsp <- momi_meta_iv(s$fs_beta, s$fs_se)      # pooled instrument -> exposure
    rop_theta <- rop_se <- rop_p <- NA_real_
    if(!is.null(rfp) && !is.null(fsp) && is.finite(fsp$b) && fsp$b != 0){
      rop_theta <- 10 * rfp$b / fsp$b
      ## delta method, same form as momi_wald
      rop_se <- 10 * sqrt(rfp$se^2/fsp$b^2 + rfp$b^2*fsp$se^2/fsp$b^4)
      rop_p  <- 2*pnorm(-abs(rop_theta/rop_se))
    }
    ## Cochran's Q -> I2 -> DerSimonian-Laird tau2 -> random-effects estimate
    w   <- 1/s$se_per10^2
    Q   <- sum(w*(s$theta_per10 - fe$b)^2)
    dfQ <- nrow(s)-1
    I2  <- max(0, 100*(Q-dfQ)/Q)
    tau2<- max(0, (Q-dfQ)/(sum(w) - sum(w^2)/sum(w)))
    wr  <- 1/(s$se_per10^2 + tau2)
    re_b<- sum(wr*s$theta_per10)/sum(wr); re_se <- sqrt(1/sum(wr))

    ## (b) stacked model with cohort fixed effects and COHORT-INTERACTED PCs (see header).
    ## Built by hand rather than via momi_wald because the covariate structure is a formula,
    ## not a column list.
    sd0 <- list()
    for(coh in unique(s$cohort)){
      sid <- momi_instrument(coh, tr); bpcol <- paste0(substr(tr,1,1),"_",BPDEF)
      zc  <- Z[score_id==sid & cohort==coh, .(IID, z)]
      dd  <- merge(A2[cohort==coh], zc, by="IID")
      if(oc %in% MOMI_LIVEBIRTH_OUTCOMES) dd <- dd[livebirth==1]
      if(!nrow(dd)) next
      sd0[[coh]] <- dd[, c("IID","cohort","AGE",pccols,"z",bpcol,oc), with=FALSE]
      setnames(sd0[[coh]], bpcol, "bp"); setnames(sd0[[coh]], oc, "y")
    }
    st <- rbindlist(sd0, fill=TRUE)
    st[, cohort := factor(cohort)]
    pcterm <- paste(sprintf("%s:cohort", pccols), collapse="+")
    rhs    <- paste0("z + cohort + AGE + ", pcterm)
    ty     <- if(oc=="BWT") "lin" else "bin"
    stk <- tryCatch({
      rf <- if(ty=="bin") summary(glm(as.formula(paste("y ~", rhs)), binomial, st))$coef
            else          summary(lm(as.formula(paste("y ~", rhs)), st))$coef
      fs <- summary(lm(as.formula(paste("bp ~", rhs)), st[is.finite(bp)]))$coef
      bo<-rf["z",1]; so<-rf["z",2]; be<-fs["z",1]; see<-fs["z",2]
      list(theta=10*bo/be, se=10*sqrt(so^2/be^2 + bo^2*see^2/be^4),
           rf_p=2*pnorm(-abs(bo/so)), F=(be/see)^2, n=nrow(st))
    }, error=function(e) NULL)

    pool[[paste(tr,oc)]] <- data.table(
      trait=tr, outcome=oc, type=ty, k_cohorts=fe$k, N_total=sum(s$N),
      instrument_mixed = uniqueN(s$instrument) > 1,
      ## ---- PRIMARY: ratio of pooled coefficients (decision 2026-07-20) ----
      ## Named theta_per10/se/p without a prefix so downstream consumers read the PRIMARY
      ## estimate by default. Burgess, Small & Thompson (Stat Methods Med Res 2017;26:2333):
      ## study-level meta-analysis of causal estimates "can accentuate weak instrument bias,
      ## and so it is preferable to combine either individual-level or summarized data".
      ## Two reasons it is the right default here:
      ##  (i) inverse-variance weights for a Wald ratio go as SE(theta)^-2 ~ betaX^2, i.e. the
      ##      weight depends on the SAME random draw that sits in the ratio's denominator, so
      ##      weights are correlated with what they weight and IVW loses its optimality;
      ## (ii) a ratio of normals has non-finite mean and Cauchy-like tails, so its SE does not
      ##      summarise its uncertainty -- pooling coefficients defers the single ill-behaved
      ##      division to the end, where one delta-method interval handles it.
      ## Empirically in this study the two agree to <8% in the SAS stratum (F 39-248) and
      ## differ by up to 13-fold with a sign flip in AFR (F 17-26) -- exactly the predicted
      ## failure mode, and itself evidence for the paper's argument about where the
      ## instrument works.
      theta_per10=rop_theta, se=rop_se, p=rop_p,
      ## ---- SECONDARY, retained for transparency ----
      meta_of_ratios_theta_per10=fe$b, meta_of_ratios_se=fe$se, meta_of_ratios_p=fe$p,
      pooling_diff_pct = if(is.finite(rop_theta) && fe$b != 0)
                        round(100*(rop_theta-fe$b)/abs(fe$b), 1) else NA_real_,
      Q=Q, df=dfQ, I2=I2, tau2=tau2,
      re_theta_per10=re_b, re_se=re_se, re_p=2*pnorm(-abs(re_b/re_se)),
      stack_theta_per10 = if(is.null(stk)) NA_real_ else stk$theta,
      stack_se          = if(is.null(stk)) NA_real_ else stk$se,
      stack_rf_p        = if(is.null(stk)) NA_real_ else stk$rf_p,
      stack_F           = if(is.null(stk)) NA_real_ else stk$F,
      min_F=min(s$F), max_F=max(s$F), any_weak=any(s$weak))
  }
  PL <- rbindlist(pool, fill=TRUE)

  ## ---------------- join S13 power ----------------
  ## S13 is per cohort x trait x outcome; attaching it here means no reader can look at a
  ## null cell without simultaneously seeing whether it COULD have been anything else.
  PW <- tryCatch(momi_read_intermediate("power_grid", P), error=function(e) NULL)
  if(!is.null(PW)){
    kk   <- intersect(c("cohort","trait","outcome"), names(PW))
    pcol <- intersect(c("power_MR_at_obs","MR_adequate","MDE_MR","MDE_observational",
                        "observed_effect"), names(PW))
    if(length(kk)==3 && length(pcol))
      MR <- merge(MR, unique(PW[, c(kk,pcol), with=FALSE]), by=kk, all.x=TRUE)
  } else cat("\nNB power_grid not found -- S10 written without the power columns.\n")

  setorder(MR, outcome, trait, ancestry, cohort)
  out  <- momi_write_table(MR, "S10_mr_panel", P)
  out2 <- momi_write_table(PL, "S10b_mr_pooled", P)
  ## published as intermediates so B23 (S11 sensitivity) can DEPEND on them and be marked
  ## stale when the MR model changes -- momi_fingerprint hashes intermediates, not tables.
  momi_save_intermediate(MR, "mr_panel", P)
  momi_save_intermediate(PL, "mr_pooled", P)

  ## ---------------- console ----------------
  cat(sprintf("\nMR panel: %d cells (%d cohorts x 2 traits x %d outcomes), exposure = %s\n",
              nrow(MR), uniqueN(MR$cohort), uniqueN(MR$outcome), BPDEF))
  cat(sprintf("covariates: %s (within-cohort PCs from B06b)\n", paste(covall, collapse=", ")))
  cat(sprintf("\nfirst-stage F: median %.1f, range %.1f-%.1f; cells with F<10: %d of %d\n",
              median(MR$F,na.rm=TRUE), min(MR$F,na.rm=TRUE), max(MR$F,na.rm=TRUE),
              sum(MR$weak,na.rm=TRUE), nrow(MR)))
  cat("\nreduced-form tests (the honest test of the causal null) surviving p<0.05:\n")
  sig <- MR[rf_p < 0.05, .(cohort, trait, outcome, N, rf_p=signif(rf_p,3), F=round(F,1))]
  if(nrow(sig)) print(sig) else cat("  none.\n")
  cat("\nExpect none-to-few: with 0 of 40 cells at 80% power (S13), this panel is designed to\n",
      "bound the effect, not to detect it. Report it as a bounding exercise.\n", sep="")

  cat("\npooled, per 10 mmHg (PRIMARY = ratio of pooled coefficients):\n")
  if(nrow(PL)) print(PL[, .(trait, outcome, k=k_cohorts,
                            primary=round(theta_per10,3), p=signif(p,3),
                            I2=round(I2), RE=round(re_theta_per10,3),
                            stack=round(stack_theta_per10,3))])
  cat("\nIf I2 is high the FE row is not interpretable as a common effect -- prefer RE, and\n",
      "note that with k=5 the between-cohort variance tau2 is itself badly estimated.\n",
      "von Hippel (BMC Med Res Methodol 2015;15:35) shows I2 is biased by ~12 percentage\n",
      "points at k=7 and worse below, so do NOT read I2=0 as evidence of homogeneity.\n", sep="")

  cat("\n=== POOLING CHECK: meta-analysing ratios vs pooling coefficients then dividing ===\n")
  if(nrow(PL)) print(PL[, .(trait, outcome,
                            meta_of_ratios = round(meta_of_ratios_theta_per10,3),
                            ratio_of_pooled = round(theta_per10,3),
                            diff_pct = pooling_diff_pct,
                            p_primary = signif(p,3),
                            p_meta_of_ratios = signif(meta_of_ratios_p,3))])
  cat("\nBurgess et al. prefer the second: meta-analysing per-cohort RATIOS upweights cohorts\n",
      "whose first stage was large by chance, and those are exactly the cohorts whose estimate\n",
      "is most pulled toward the confounded observational association. With our first stages\n",
      "(F 17-250) the two should agree to a few percent. A large diff_pct means the headline\n",
      "should be taken from ratio_of_pooled, and the divergence explained.\n", sep="")

  list(n=nrow(MR),
       key=sprintf("cells=%d; F median=%.1f weak(F<10)=%d; reduced-form p<0.05: %d; pooled rows=%d",
                   nrow(MR), median(MR$F,na.rm=TRUE), sum(MR$weak,na.rm=TRUE),
                   sum(MR$rf_p<0.05, na.rm=TRUE), nrow(PL)),
       outputs=c(basename(out), basename(out2), "mr_panel.rds", "mr_pooled.rds"))
})
