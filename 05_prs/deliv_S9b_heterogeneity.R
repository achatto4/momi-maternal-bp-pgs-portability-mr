#!/usr/bin/env Rscript
# ============================================================
# deliv_S9b_heterogeneity.R  [ B18b -> Supp Table S9b + early per-site forest ]
# Is the observational BP -> outcome effect poolable across cohorts?
#
# WHY THIS EXISTS. S9 produced fully-adjusted per-cohort estimates spanning OR 0.94 to 1.70
# for SBP -> LBW, with AMANHI-Bangladesh and AMANHI-Pakistan null across ALL six
# outcome x trait combinations while Pemba, GAPPS-Bangladesh and Zambia show strong effects.
# T4, F4 and the whole MR panel are designed around a single pooled estimate per outcome.
# If that pooling is not defensible, the design changes -- and it is far cheaper to find out
# now than after four more deliverables are built on the assumption.
#
# WHAT IT DOES NOT DO. It does not attempt to explain the two null cohorts. Three separate
# explanations have been floated and refuted today (dual-platform precision, GSA/dosage mix,
# MAF filtering), so this module measures the heterogeneity and reports the two MECHANICAL
# candidates -- outcome prevalence and exposure variance -- without interpreting them.
# A null association is unremarkable if a cohort has almost no cases or almost no BP spread.
#
# Reads tables/S9_confounding.tsv (B18) and analytic_mothers.
# Writes: tables/S9b_heterogeneity.tsv, figures/SF4a_forest_observational.{png,pdf}
#   Rscript deliv_S9b_heterogeneity.R
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

momi_deliverable("S9b_heterogeneity", script="05_prs/deliv_S9b_heterogeneity.R",
                 inputs="analytic_mothers", P=P, body=function(ctx){
  f9 <- file.path(P$tables, "S9_confounding.tsv")
  if(!file.exists(f9)) stop("S9_confounding.tsv not found — run B18 first")
  S9 <- fread(f9)
  A  <- momi_read_intermediate("analytic_mothers", P)
  pd <- MOMI_BPDEF_DEFAULT

  ## fully-adjusted, fixed-sample, primary definition -- the estimates T4 would pool
  ## BWT is included DELIBERATELY and is the decisive row. LBW and SGA are binary outcomes
  ## defined against thresholds and an international reference; SGA prevalence runs from
  ## 10.7% (Pemba) to 44.0% (AMANHI-Bangladesh), which means the reference is badly
  ## mismatched in the South-Asian cohorts and "SGA" there is not measuring what it measures
  ## elsewhere. Birth weight in grams has no threshold and no reference population, so:
  ##   heterogeneous in BWT too  -> the CAUSAL EFFECT genuinely varies across populations
  ##   homogeneous in BWT only   -> the heterogeneity is an artefact of the binary definitions
  ## The paper's thesis depends on which.
  OUTS <- c("PTB","LBW","SGA","BWT")
  PC <- S9[scope != "pooled" & sample == "common" & model == "M3_plus_EDU" &
           definition == pd & outcome %in% OUTS & is.finite(est) & se > 0]
  PL <- S9[scope == "pooled" & sample == "common" & model == "M3_plus_EDU" &
           definition == pd & outcome %in% OUTS]

  ## ---- Cochran's Q, I^2, tau^2, and the inverse-variance pooled estimate ----
  het <- PC[, {
    b <- est; s <- se; w <- 1/s^2
    mu <- sum(w*b)/sum(w); semu <- sqrt(1/sum(w))
    Q  <- sum(w*(b-mu)^2); df <- .N - 1
    I2 <- if(Q > df) 100*(Q-df)/Q else 0
    C  <- sum(w) - sum(w^2)/sum(w)
    tau2 <- if(Q > df && C > 0) (Q-df)/C else 0
    ## Q and I^2 are scale-free, so they work on the raw estimate for both outcome types.
    ## Only the REPORTED effect differs: exp() for binary (OR), identity for BWT (grams).
    isbin <- .BY$outcome != "BWT"
    tf <- function(z) if(isbin) exp(z) else z
    .(k=.N, scale = if(isbin) "OR per 10 mmHg" else "grams per 10 mmHg",
      pooled_eff = round(tf(mu),3),
      pooled_lo  = round(tf(mu-1.96*semu),3),
      pooled_hi  = round(tf(mu+1.96*semu),3),
      Q = round(Q,2), df = df,
      Q_p = signif(pchisq(Q, df, lower.tail=FALSE), 3),
      I2_pct = round(I2,1), tau2 = round(tau2,4),
      eff_min = round(tf(min(b)),3), eff_max = round(tf(max(b)),3),
      ## random-effects pooled, for comparison when tau2 > 0
      re_eff = { wr <- 1/(s^2 + tau2); round(tf(sum(wr*b)/sum(wr)),3) })
  }, by=.(trait, outcome)]

  ## direct pooled model (one regression with a site term) vs meta-analytic pooling.
  ## They answer different questions: the direct model assumes ONE common effect, the
  ## meta-analysis allows the effect to differ by cohort. Divergence is itself diagnostic.
  PL[, direct_eff := fifelse(outcome=="BWT", est, OR)]
  het <- merge(het, PL[, .(trait, outcome, direct_pooled_eff = direct_eff)],
               by=c("trait","outcome"), all.x=TRUE)
  het[, poolable := fifelse(I2_pct < 50 & Q_p > 0.05, "yes",
                     fifelse(I2_pct < 75, "borderline", "NO"))]
  setorder(het, trait, outcome)
  out <- momi_write_table(het, "S9b_heterogeneity", P)

  ## ---- the two mechanical explanations for a null cohort ----
  A[, livebirth_ok := livebirth == 1]
  diag <- rbindlist(lapply(MOMI_COH_ALL, function(coh){
    a <- A[cohort == coh]
    lb <- a[livebirth_ok %in% TRUE]
    data.table(cohort = coh, ancestry = unname(MOMI_ANC[coh]),
      n_mothers   = nrow(a),
      sd_SBP      = round(sd(a[[paste0("S_",pd)]], na.rm=TRUE), 2),
      sd_DBP      = round(sd(a[[paste0("D_",pd)]], na.rm=TRUE), 2),
      PTB_n       = sum(a$PTB == 1, na.rm=TRUE),
      PTB_pct     = round(100*mean(a$PTB == 1, na.rm=TRUE), 1),
      LBW_n       = sum(lb$LBW == 1, na.rm=TRUE),
      LBW_pct     = round(100*mean(lb$LBW == 1, na.rm=TRUE), 1),
      SGA_n       = sum(lb$SGA == 1, na.rm=TRUE),
      SGA_pct     = round(100*mean(lb$SGA == 1, na.rm=TRUE), 1),
      pct_BWT_miss= round(100*mean(is.na(lb$BWT)), 1))
  }))
  momi_write_table(diag, "S9b_cohort_diagnostics", P)

  ## ---- forest: per-cohort with the pooled estimate, faceted ----
  FP <- PC[outcome != "BWT", .(cohort = scope, trait, outcome, OR, lo, hi, n)]
  FP[, Cohort := unname(MOMI_DISPLAY[cohort])]
  ## Show BOTH pooled estimates (decision 2026-07-19: report aggregate AND component-wise).
  ## Random-effects is the defensible summary here -- fixed-effect assumes one common effect,
  ## which I2 = 70-92% rejects. Plotting them together makes the gap between the two visible
  ## rather than leaving the reader to trust whichever we chose.
  ## RE interval uses the RE weights: se_RE = sqrt(1/sum(1/(se^2 + tau2))).
  reint <- PC[outcome != "BWT"][het[outcome != "BWT", .(trait, outcome, tau2)],
                                on=.(trait, outcome)][
              , .(se_re = sqrt(1/sum(1/(se^2 + tau2[1])))), by=.(trait, outcome)]
  RE <- merge(het[outcome != "BWT", .(trait, outcome, re_eff)], reint,
              by=c("trait","outcome"))
  RE[, `:=`(Cohort = "POOLED (random-effects)",
            OR = re_eff, lo = round(exp(log(re_eff) - 1.96*se_re),3),
            hi = round(exp(log(re_eff) + 1.96*se_re),3), n = NA_integer_)]

  PO <- het[outcome != "BWT", .(trait, outcome, Cohort = "POOLED (fixed-effect)",
                OR = pooled_eff, lo = pooled_lo, hi = pooled_hi, n = NA_integer_)]
  FP <- rbind(FP[, .(Cohort, trait, outcome, OR, lo, hi, n)], PO,
              RE[, .(Cohort, trait, outcome, OR, lo, hi, n)], fill=TRUE)
  pooled_labs <- c("POOLED (fixed-effect)", "POOLED (random-effects)")
  FP[, Cohort := factor(Cohort, levels = c(sort(setdiff(unique(Cohort), pooled_labs)),
                                           pooled_labs))]
  FP[, is_pooled := Cohort %in% pooled_labs]

  g <- ggplot(FP, aes(x = OR, y = Cohort, colour = is_pooled)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.5) +
    geom_point(aes(shape = is_pooled), size = 2.1) +
    facet_grid(trait ~ outcome) +
    scale_x_log10() +
    scale_colour_manual(values = c(`FALSE` = "grey25", `TRUE` = "#b2182b"), guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
    labs(x = "Adjusted OR per 10 mmHg (log scale)", y = NULL,
         title = "Observational BP-outcome association by cohort",
         subtitle = sprintf("Exposure = %s; fully adjusted (age, gravidity, BMI, education); fixed analytic sample", pd)) +
    theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank())
  outs <- momi_save_fig(g, "SF4a_forest_observational", width = 9, height = 6, P = P)

  ## ---- console ----
  cat("\n== heterogeneity across the five cohorts ==\n")
  print(het[, .(trait, outcome, scale, k, eff_min, eff_max, pooled_eff,
                direct_pooled_eff, I2_pct, Q_p, poolable)])
  cat("\nI2 above 75% means a single pooled number misrepresents the cohorts. Compare\n",
      "pooled_OR (meta-analytic) against direct_pooled_OR (one regression with a site term):\n",
      "they diverge exactly when the common-effect assumption fails.\n", sep="")

  cat("\n== mechanical checks: is a null cohort simply short of cases or BP spread? ==\n")
  print(diag[, .(cohort, ancestry, sd_SBP, sd_DBP, PTB_n, PTB_pct, LBW_n, LBW_pct,
                 SGA_n, SGA_pct, pct_BWT_miss)])
  cat("\nAMANHI-Bangladesh and AMANHI-Pakistan were null on all six estimates. If their case\n",
      "counts and BP standard deviations are comparable to the others, that is NOT explained\n",
      "and must be reported as heterogeneity rather than attributed to power.\n", sep="")

  worst <- het[which.max(I2_pct)]
  list(n=nrow(het),
       key=sprintf("I2 %.0f-%.0f%%; worst %s %s I2=%.0f%% Q_p=%s; poolable: %s",
                   min(het$I2_pct), max(het$I2_pct), worst$trait, worst$outcome,
                   worst$I2_pct, format(worst$Q_p),
                   paste(sprintf("%s=%s", het$outcome, het$poolable)[!duplicated(het$outcome)],
                         collapse=",")),
       outputs=c(basename(out), "S9b_cohort_diagnostics.tsv", basename(outs)))
})
