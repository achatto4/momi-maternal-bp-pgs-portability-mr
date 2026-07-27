#!/usr/bin/env Rscript
# ============================================================
# deliv_T4_triangulation.R  [ B21 -> Table 4 + Figure 4, the paper's headline ]
#
# ------------------------------------------------------------------
# WHAT TRIANGULATION MEANS HERE, AND WHY IT IS THE RIGHT FRAME FOR THIS DATASET
#
# Three lines of evidence on the same question, each with DIFFERENT and largely NON-OVERLAPPING
# sources of bias. Agreement between them is informative precisely because a single artefact
# is unlikely to move all three the same way:
#
#   1. OBSERVATIONAL (S9/B18)   -- measured BP vs outcome, adjusted for age, gravidity, BMI,
#        education. Well powered. Vulnerable to residual confounding and reverse causation.
#   2. MR, THIS STUDY (S10/B22) -- genetic instrument, immune to reverse causation and to
#        confounding by anything not correlated with genotype. Badly underpowered here:
#        S13 found 0 of 40 cells reach 80% power.
#   3. MR, EXTERNAL (ref/external_mr_estimates.tsv) -- Morales-Berstein 2026, two-sample IVW
#        in European samples. Well powered. Vulnerable to whatever does not transfer from a
#        European population to South Asian and African pregnancy cohorts -- which is exactly
#        what Part I of this paper measures.
#
# THE ARGUMENT THE TABLE MAKES. Our MR is too small to reject the null for most outcomes. But
# "underpowered" and "null" are different claims, and the distinction is testable: if our
# POINT ESTIMATES land on the external ones despite our p-values, that is an underpowered
# REPLICATION, not an absence of effect. Reporting "no evidence of an effect on SGA" when our
# OR is 1.17 and the external is 1.16 would be a straightforward misreading of our own data.
#
# THIS IS WHY THE TABLE LEADS WITH ESTIMATES AND INTERVALS, NOT SIGNIFICANCE. Every row shows
# the effect and its CI from all three lines, and the p-value is present but not the organising
# principle. A reader should be able to see agreement or disagreement in the point estimates
# without consulting a single p.
#
# STRATA. SAS, AFR and overall are reported for our MR, per Anagh 2026-07-20: the strata may
# reflect genuine causal heterogeneity, not merely differing precision, so neither is treated
# as a nuisance. The AFR rows carry their F range because Part I established the instrument
# does not transfer there -- WITHOUT THAT ANNOTATION A WIDE AFR NULL READS AS EVIDENCE OF
# ABSENCE, WHICH IT IS NOT. "No effect in African cohorts" and "no instrument in African
# cohorts" are different findings and only the second is supported.
#
# SCALE. Two panels rather than one. PTB/LBW/SGA are odds ratios per 10 mmHg and go directly
# against the external estimates, which are published in the same units -- no conversion, no
# assumptions. Birthweight is grams per 10 mmHg on its own axis; a log scale for grams is
# meaningless, and standardising everything to per-SD would require assuming the external
# study's BP distribution, which we cannot verify.
#
# NO EXTERNAL COMPARATOR EXISTS FOR CONTINUOUS BIRTHWEIGHT in the reference file (the source
# reports HBW 0.76 and LGA 0.87, implying a growth effect, but not grams). That is a FILLABLE
# gap, marked as such rather than silently omitted -- our most robust finding is the one with
# no external anchor, and a reader deserves to know that.
# ------------------------------------------------------------------
#   Writes: tables/T4_triangulation.tsv, figures/F4_triangulation.{png,pdf}
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

EXT <- file.path(PIPE, "ref", "external_mr_estimates.tsv")

momi_deliverable("T4_triangulation", script="05_prs/deliv_T4_triangulation.R",
                 inputs="mr_panel;mr_pooled;power_grid", P=P, stop_on_error=FALSE,
                 body=function(ctx){

  MR <- tryCatch(momi_read_intermediate("mr_panel",  P), error=function(e) NULL)
  PL <- tryCatch(momi_read_intermediate("mr_pooled", P), error=function(e) NULL)
  if(is.null(MR) || is.null(PL)) return(list(skip=TRUE, reason="no mr_panel/mr_pooled -- run B22"))

  ## ---------- line 1: observational, from the S9 table ----------
  ## S9 is a table not an intermediate, so it is read defensively and its absence degrades
  ## the figure rather than failing the build.
  s9f <- file.path(P$tables, "S9_confounding.tsv")
  OBS <- if(file.exists(s9f)) fread(s9f) else NULL
  if(is.null(OBS)) cat("\nNB S9_confounding.tsv absent -- observational arm omitted.\n")

  ## ---------- line 3: external ----------
  EX <- if(file.exists(EXT)) fread(EXT) else NULL
  if(is.null(EX)) cat("\nNB ref/external_mr_estimates.tsv absent -- external arm omitted.\n")

  OUTS <- c("PTB","LBW","SGA","BWT")
  rows <- list()

  for(tr in c("SBP","DBP")) for(oc in OUTS){
    ty <- if(oc=="BWT") "lin" else "bin"

    ## --- our MR, three strata ---
    ## POOLING METHOD, changed 2026-07-20 after methodological review.
    ##
    ## These strata previously inverse-variance meta-analysed the per-cohort WALD RATIOS.
    ## Burgess, Small & Thompson (Stat Methods Med Res 2017;26:2333) advise against that:
    ## "a study-level meta-analysis can accentuate weak instrument bias, and so it is
    ## preferable to combine either individual-level or summarized data from each study."
    ## The reason it matters here is not only bias: with first stages spanning F 17-250,
    ## meta-of-ratios weights each cohort by the precision of its RATIO, which mixes the two
    ## regressions together, whereas pooling each regression separately and dividing once
    ## does not. In S10 the two differed by 7-29% -- far more than the ~1/F bias argument
    ## alone predicts -- so this is a real choice of estimator, not a rounding detail.
    ##
    ## We now lead with RATIO-OF-POOLED and carry meta-of-ratios alongside for comparison,
    ## because the former is what the methodological literature recommends and because it
    ## does not let a cohort whose first stage was large by chance dominate the pooled result.
    for(st in c("SAS","AFR","all")){
      s <- MR[trait==tr & outcome==oc & is.finite(theta_per10) & se_per10 > 0]
      s[, anc := MOMI_ANC[cohort]]
      if(st != "all") s <- s[anc==st]
      if(nrow(s) < 1) next
      fe  <- momi_meta_iv(s$theta_per10, s$se_per10)          # legacy: meta of ratios
      rfp <- momi_meta_iv(s$rf_beta, s$rf_se)                 # pooled instrument -> outcome
      fsp <- momi_meta_iv(s$fs_beta, s$fs_se)                 # pooled instrument -> exposure
      if(is.null(fe)) next
      rop <- rop_se <- NA_real_
      if(!is.null(rfp) && !is.null(fsp) && is.finite(fsp$b) && fsp$b != 0){
        rop    <- 10 * rfp$b / fsp$b
        rop_se <- 10 * sqrt(rfp$se^2/fsp$b^2 + rfp$b^2*fsp$se^2/fsp$b^4)
      }
      ## primary = ratio-of-pooled where estimable, else fall back to meta-of-ratios
      est <- if(is.finite(rop)) rop else fe$b
      se  <- if(is.finite(rop)) rop_se else fe$se
      rows[[paste("mr",tr,oc,st)]] <- data.table(
        trait=tr, outcome=oc, type=ty,
        line = paste0("MR (this study, ", st, ")"),
        arm  = "MR_ours", stratum = st,
        k = fe$k, N = sum(s$N), n_case = sum(s$n_case, na.rm=TRUE),
        est = est, lo = est-1.96*se, hi = est+1.96*se,
        p = 2*pnorm(-abs(est/se)),
        pooling = if(is.finite(rop)) "ratio_of_pooled" else "meta_of_ratios",
        est_meta_of_ratios = fe$b, p_meta_of_ratios = fe$p,
        pooling_diff_pct = if(is.finite(rop) && fe$b != 0)
                             round(100*(rop-fe$b)/abs(fe$b),1) else NA_real_,
        F_min = min(s$F, na.rm=TRUE), F_max = max(s$F, na.rm=TRUE))
    }

    ## --- observational, pooled, fully adjusted, own-sample ---
    if(!is.null(OBS)){
      o <- OBS[trait==tr & outcome==oc & scope=="pooled" &
               model=="M3_plus_EDU" & sample=="own" &
               definition==MOMI_BPDEF_DEFAULT]
      if(nrow(o)){
        e <- if(ty=="bin") log(o$OR[1]) else o$est[1]
        l <- if(ty=="bin") log(o$lo[1])  else o$lo[1]
        h <- if(ty=="bin") log(o$hi[1])  else o$hi[1]
        rows[[paste("obs",tr,oc)]] <- data.table(
          trait=tr, outcome=oc, type=ty, line="Observational (adjusted)",
          arm="observational", stratum="all", k=5L, N=o$n[1], n_case=NA_integer_,
          est=e, lo=l, hi=h, p=o$p[1], F_min=NA_real_, F_max=NA_real_)
      }
    }

    ## --- external MR ---
    ## Retrieved in full from Additional File 3, Supplementary Table 6 on 2026-07-20, so BOTH
    ## traits and all four MOMI outcomes now have a comparator. Two scales are involved and
    ## they must not be mixed up:
    ##   binary outcomes    -> odds ratio per 10 mmHg   (stored as OR, logged here)
    ##   continuous BWT/GAd -> SD of the outcome per 10 mmHg, NOT grams. The source paper
    ##                         never uses grams. Our estimate is converted to SD below using
    ##                         the SD of birthweight in our own live-birth sample.
    if(!is.null(EX)){
      x <- EX[exposure==tr & outcome==oc & is.finite(estimate)]
      if(nrow(x)){
        is_sd <- grepl("^SD", x$units[1])
        ## for binary we store on the log-OR scale (exponentiated later with everything else);
        ## for SD-scale continuous the estimate is already on the reporting scale, but `est`
        ## must stay in OUR units (grams) so it plots with our rows -- multiply by our SD.
        if(ty=="bin"){
          e <- log(x$estimate[1]); l <- log(x$ci_lo[1]); h <- log(x$ci_hi[1])
        } else if(is_sd){
          e <- x$estimate[1]; l <- x$ci_lo[1]; h <- x$ci_hi[1]   # rescaled to grams below
        } else { e <- x$estimate[1]; l <- x$ci_lo[1]; h <- x$ci_hi[1] }
        rows[[paste("ext",tr,oc)]] <- data.table(
          trait=tr, outcome=oc, type=ty, line="MR (external, European)",
          arm="MR_external", stratum="EUR", k=NA_integer_, N=NA_integer_, n_case=NA_integer_,
          est=e, lo=l, hi=h, p=NA_real_, F_min=NA_real_, F_max=NA_real_,
          ext_units=x$units[1], ext_on_sd_scale=is_sd)
      }
    }
  }
  if(!length(rows)) return(list(skip=TRUE, reason="no triangulation rows assembled"))
  T4 <- rbindlist(rows, fill=TRUE)

  ## ---- power, joined onto our MR rows: a null must never be read without it ----
  PW <- tryCatch(momi_read_intermediate("power_grid", P), error=function(e) NULL)
  if(!is.null(PW) && all(c("trait","outcome","power_MR_at_obs") %in% names(PW))){
    pw <- PW[, .(power_max = max(power_MR_at_obs, na.rm=TRUE)), by=.(trait, outcome)]
    T4 <- merge(T4, pw, by=c("trait","outcome"), all.x=TRUE)
    T4[arm != "MR_ours", power_max := NA_real_]
  }

  ## reporting scale: OR for binary, grams for continuous
  T4[, `:=`(effect = fifelse(type=="bin", exp(est), est),
            eff_lo = fifelse(type=="bin", exp(lo),  lo),
            eff_hi = fifelse(type=="bin", exp(hi),  hi))]
  T4[, units := fifelse(type=="bin", "OR per 10 mmHg", "grams per 10 mmHg")]

  ## ---- birthweight ALSO in SD units, so it is comparable to the external study ----------
  ## Established 2026-07-20 by reading the source paper: Morales-Berstein et al. report
  ## continuous birth weight "in standard deviations as a secondary outcome". The word
  ## "grams" does not appear anywhere in that paper. So a grams-scale external comparator
  ## does NOT exist and cannot be obtained -- the fix is to convert OUR estimate to SD units
  ## rather than to keep waiting for theirs in grams.
  ##
  ## The SD used is that of birthweight in OUR analytic sample (live births, Part II), which
  ## is the correct denominator for putting our effect on a per-SD scale. NB this makes the
  ## two comparable in SCALE but not identical in DEFINITION: their SD comes from their own
  ## cohorts, and birthweight distributions differ between European and South Asian
  ## populations, so a like-for-like reading still carries that caveat. Recorded here rather
  ## than left for a reader to assume equivalence.
  A_bw <- momi_read_intermediate("analytic_mothers", P)
  sd_bwt <- sd(A_bw[livebirth==1 & is.finite(BWT), BWT], na.rm=TRUE)

  ## The external continuous rows arrive on the SD scale. Put them into GRAMS using OUR SD so
  ## every row of the continuous panel shares one axis. This is the only defensible direction:
  ## converting ours to their SD would need the SD of birthweight in THEIR cohorts, which the
  ## paper does not report. Flagged on the row (`ext_on_sd_scale`) so the assumption is visible.
  if("ext_on_sd_scale" %in% names(T4)){
    T4[arm=="MR_external" & ext_on_sd_scale==TRUE,
       `:=`(effect = est*sd_bwt, eff_lo = lo*sd_bwt, eff_hi = hi*sd_bwt)]
  }
  T4[type=="lin", `:=`(effect_sd = effect/sd_bwt,
                       eff_lo_sd = eff_lo/sd_bwt,
                       eff_hi_sd = eff_hi/sd_bwt)]
  T4[type=="lin", units := sprintf("grams per 10 mmHg (SD of BWT = %.0f g)", sd_bwt)]

  ## flag the gap rather than leaving it silent (see header)
  T4[, external_available := any(arm=="MR_external"), by=.(trait, outcome)]

  setorder(T4, outcome, trait, arm, stratum)
  out <- momi_write_table(T4, "T4_triangulation", P)

  ## ---------------- figure: two panels ----------------
  mk <- function(D, id, title, xlab, expo){
    if(!nrow(D)) return(character(0))
    D <- copy(D)
    D[, lab := fifelse(is.finite(F_min) & arm=="MR_ours",
                       sprintf("%s  [F %.0f-%.0f]", line, F_min, F_max), line)]
    D[, row := factor(paste(outcome, lab, sep=" | "),
                      levels=rev(unique(paste(outcome, lab, sep=" | "))))]
    g <- ggplot(D, aes(x=effect, y=row, colour=arm)) +
      geom_vline(xintercept=if(expo) 1 else 0, linetype="dashed",
                 linewidth=.3, colour="grey40") +
      geom_errorbarh(aes(xmin=eff_lo, xmax=eff_hi), height=.18, linewidth=.4) +
      geom_point(size=1.9) +
      facet_grid(trait ~ ., scales="free_y", space="free_y") +
      scale_colour_manual(values=c(observational="#4C72B0", MR_ours="#C44E52",
                                   MR_external="#55A868"), name=NULL) +
      labs(x=xlab, y=NULL, title=title,
           subtitle="Three lines of evidence with largely non-overlapping biases. Read the estimates and intervals, not the p-values.") +
      theme_minimal(base_size=8) +
      theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
            legend.position="bottom",
            plot.subtitle=element_text(size=6.5, colour="grey35"))
    if(expo) g <- g + scale_x_continuous(trans="log10")
    momi_save_fig(g, id, width=8.5, height=max(4.5, 0.24*nrow(D)+2), P=P)
  }

  f1 <- mk(T4[type=="bin"], "F4_triangulation",
           "Maternal blood pressure and perinatal outcomes: triangulation",
           "OR per 10 mmHg (log scale)", expo=TRUE)
  f2 <- mk(T4[type=="lin"], "F4b_triangulation_bwt",
           "Maternal blood pressure and birthweight: triangulation",
           "grams per 10 mmHg", expo=FALSE)

  ## ---------------- console ----------------
  cat("\n=== T4: ours vs external, binary outcomes (the replication argument) ===\n")
  cmp <- dcast(T4[type=="bin" & (arm=="MR_external" | (arm=="MR_ours" & stratum %in% c("all","SAS")))],
               trait + outcome ~ line, value.var="effect")
  print(cmp)
  cat("\nIf our point estimates sit near the external ones while our CIs are wide, that is an\n",
      "UNDERPOWERED REPLICATION. Describing such a cell as 'no evidence of an effect' would\n",
      "misstate our own data.\n", sep="")

  cat(sprintf("\n=== birthweight in SD units (external study reports per-SD, not grams; SD here = %.0f g) ===\n", sd_bwt))
  bw <- T4[type=="lin" & arm=="MR_ours",
           .(trait, stratum, k, grams=round(effect,1),
             per_SD=round(effect_sd,4), lo_SD=round(eff_lo_sd,4), hi_SD=round(eff_hi_sd,4),
             p=signif(p,3))]
  if(nrow(bw)) print(bw)
  cat("Morales-Berstein report continuous birthweight in SD (Additional File 3, Suppl Table 6).\n",
      "Once those values are entered in ref/external_mr_estimates.tsv the per_SD column above is\n",
      "the like-for-like comparison. Caveat: their SD comes from European cohorts and ours from\n",
      "South Asian and African ones, so the scales are comparable but not identical.\n", sep="")

  cat("\n=== outcomes with no external comparator ===\n")
  gap <- unique(T4[external_available==FALSE, .(trait, outcome)])
  if(nrow(gap)) print(gap) else cat("  none\n")
  cat("NB continuous birthweight is our most robust finding and has no external anchor in\n",
      "ref/external_mr_estimates.tsv. Retrieving a grams-per-mmHg estimate is a FILLABLE gap.\n", sep="")

  cat("\n=== POOLING: ratio-of-pooled (primary) vs meta-of-ratios (legacy) ===\n")
  pc <- T4[arm=="MR_ours", .(trait, outcome, stratum, k,
                             primary=round(effect,3),
                             meta_of_ratios=round(fifelse(type=="bin",
                                 exp(est_meta_of_ratios), est_meta_of_ratios),3),
                             diff_pct=pooling_diff_pct,
                             p=signif(p,3), p_legacy=signif(p_meta_of_ratios,3))]
  if(nrow(pc)) print(pc[order(outcome, trait, stratum)])
  cat("\nThe SAS rows are the ones the paper leads with, so their diff_pct is what matters.\n",
      "A large difference means the earlier headline was partly an artefact of weighting\n",
      "cohorts by the precision of their Wald ratio rather than of the underlying regressions.\n", sep="")

  cat("\n=== AFR strata: check F before reading any null ===\n")
  print(T4[arm=="MR_ours" & stratum=="AFR",
           .(trait, outcome, k, N, effect=round(effect,3),
             lo=round(eff_lo,3), hi=round(eff_hi,3), F_min=round(F_min))])

  list(n=nrow(T4),
       key=sprintf("rows=%d; outcomes=%d; arms=%d; external comparators=%d",
                   nrow(T4), uniqueN(T4$outcome), uniqueN(T4$arm),
                   nrow(T4[arm=="MR_external"])),
       outputs=c(basename(out), basename(c(f1,f2))))
})
