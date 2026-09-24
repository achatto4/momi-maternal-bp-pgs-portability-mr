## -------------------------------------------------------------------------------
## 04_mr_cohort_models_and_pooling.R
##
## One-sample Mendelian randomization, cohort by cohort, and the pooling of those
## cohort estimates.
##
##   Instrument   one polygenic score per cohort and trait, chosen in advance and never
##                re-selected, re-oriented or flipped here.
##   Exposure     the mean of a woman's valid antenatal systolic or diastolic readings.
##   Outcomes     preterm birth, low birth weight, small for gestational age, and birth
##                weight in grams.
##   Sample       one outcome-specific complete-case sample per cell. A woman enters a
##                cell if she is in the cleaned sample and the cohort and has the score,
##                the corresponding mean antenatal pressure, the outcome and its
##                denominator condition, maternal age, genotyping technology and
##                principal components 1-5. Both stages are fitted in exactly those women.
##   Covariates   maternal age, genotyping technology, principal components 1-5, and
##                nothing else.
##   First stage  blood pressure on the score and the covariates, by ordinary least
##                squares, re-estimated in every outcome-specific sample.
##   Reduced form logistic regression for the binary outcomes and ordinary least squares
##                in grams for birth weight, on the same right-hand side. A cell that
##                separates or fails to converge is stopped and reported as stopped; no
##                penalised estimator is substituted.
##   Estimate     the Wald ratio of the reduced-form to the first-stage coefficient, with
##                a delta-method standard error, scaled to a 10 mmHg difference.
##   Pooling      random-effects meta-analysis by restricted maximum likelihood with a z
##                test, using the same model for every pooled estimate: the South Asian
##                cohorts, the African cohorts, and all five together. Binary outcomes are
##                pooled on the log-odds scale and exponentiated afterwards; birth weight
##                is pooled in grams. The subgroup estimates are descriptive: no
##                meta-regression and no formal ancestry-interaction test is fitted.
## -------------------------------------------------------------------------------
if (!exists("MOMI_ROOT")) MOMI_ROOT <- getwd()
if (!exists("config")) source(file.path(MOMI_ROOT,
  if (file.exists(file.path(MOMI_ROOT, "config.R"))) "config.R" else "config.example.R"))
source(file.path(MOMI_ROOT, "analysis", "00_functions.R"))
suppressMessages(library(metafor))

D0  <- readRDS(file.path(DERIVED, "analysis_sample_with_scores.rds"))
AN   <- D0$AN
ZL   <- D0$ZL
ZALT <- D0$ZALT

RES <- list(); FS <- list(); RF <- list(); TL <- list(); ZCHK <- list()
SAMPLE_IDX <- list()
for(coh in COH5){
  base <- AN[cohort == coh]
  for(tr in TRAITS){
    sid <- INSTR[coh][[tr]]
    D   <- merge(base, ZL[[paste(coh, sid)]][, .(IID, z)], by="IID", all.x=TRUE)
    setorder(D, IID)
    D[, bp := get(BPCOL[[tr]])]

    zz <- merge(D[, .(IID, z)], ZALT[[paste(coh, sid)]], by="IID")
    zo <- merge(D[, .(IID, z)], ZL[[paste(coh, INSTR[coh][[setdiff(TRAITS, tr)]])]][, .(IID, z_other=z)], by="IID")
    ZCHK[[length(ZCHK) + 1]] <- data.table(cohort=coh, bp_trait=tr, pgs_id=sid, n_compared=nrow(zz),
      max_abs_diff_vs_independent_rebuild = if(nrow(zz)) max(abs(zz$z - zz$z_alt)) else NA_real_,
      max_abs_diff_vs_other_trait_score   = if(nrow(zo)) max(abs(zo$z - zo$z_other)) else NA_real_)

    for(oc in OUTSPEC$outcome){
      osp    <- OUTSPEC[oc]
      binary <- osp$type == "binary"
      lb <- D$livebirth == 1; lb[is.na(lb)] <- FALSE
      if(oc == "PTB"){ y <- as.integer(D$PTB == 1);    den <- !is.na(D$PTB) }
      if(oc == "LBW"){ y <- as.integer(D$BWT < 2500);  den <- lb & !is.na(D$BWT) }
      if(oc == "SGA"){ y <- as.integer(D$SGA);         den <- lb & !is.na(D$SGA) }
      if(oc == "BWT"){ y <- as.numeric(D$BWT);         den <- lb & !is.na(D$BWT) }
      den[is.na(den)] <- FALSE
      ok <- den & is.finite(D$z) & is.finite(D$bp) & is.finite(D$AGE) &
            D$technology %chin% TECH_ALL & D$has_pc5 & is.finite(y)
      ok[is.na(ok)] <- FALSE
      idx <- which(ok)
      SAMPLE_IDX[[paste(coh, tr, oc)]] <- idx

      Dc <- as.data.frame(D[idx, c("technology","AGE","z","bp", PCS5), with=FALSE])
      Dc$y <- y[idx]
      rownames(Dc) <- paste0("w", idx)

      tech_exc <- paste(coh, tr, oc, sep="|") %chin% TECH_EXCEPT$key
      tech_raw <- as.character(Dc$technology)
      tech_use <- if(tech_exc) ifelse(tech_raw %chin% TECH_COLLAPSE_FROM, TECH_COLLAPSE_TO, tech_raw)
                  else         tech_raw
      lv_all    <- if(tech_exc) TECH_EXCEPT_LEVELS else TECH_ALL
      tech_code <- if(tech_exc) TECH_CODE_EXC      else TECH_CODE_STD
      tabt_raw <- table(factor(tech_raw, levels=TECH_ALL))
      tabt <- table(factor(tech_use, levels=lv_all))
      lv   <- lv_all[tabt > 0]
      ref  <- if(TECH_REF %in% lv) TECH_REF else names(tabt)[which.max(as.integer(tabt))]
      lv   <- c(ref, setdiff(lv, ref))
      Dc$tech <- factor(tech_use, levels=lv)
      tech_in <- length(lv) > 1
      rhs  <- paste(c("z", "AGE", if(tech_in) "tech", PCS5), collapse=" + ")
      n    <- nrow(Dc)
      npar <- 2L + 1L + (if(tech_in) length(lv) - 1L else 0L) + length(PCS5)
      nev  <- if(binary) sum(Dc$y == 1L) else NA_integer_

      mk_tl <- function(L, sel, inmod, role, into)
        data.table(cohort=coh, cohort_display=unname(DISPLAY[coh]), bp_trait=tr, pgs_id=sid, outcome=oc,
                   technology_level=L, level_role=role, in_model_levels=inmod,
                   is_reference=inmod && identical(L, ref),
                   n=sum(sel), n_events=if(binary) sum(Dc$y[sel] == 1L) else NA_integer_,
                   technology_coding=tech_code,
                   collapsed_from=if(tech_exc && identical(L, TECH_COLLAPSE_TO)) jn(TECH_COLLAPSE_FROM) else "",
                   collapsed_into=into)
      lvcnt <- rbindlist(c(
        lapply(lv_all, function(L) mk_tl(L, tech_use == L, L %chin% lv,
                                         if(L %chin% lv) "modelled" else "absent", "")),
        if(tech_exc) lapply(TECH_COLLAPSE_FROM, function(L)
                              mk_tl(L, tech_raw == L, FALSE, "collapsed component", TECH_COLLAPSE_TO))
        else list()))
      TL[[length(TL) + 1]] <- lvcnt
      zero_ev_level <- binary && tech_in &&
                       any(lvcnt[in_model_levels == TRUE & n > 0, n_events == 0L | n_events == n])

      stops <- character(0); notes <- character(0); wmsg <- character(0)
      fitted_both <- FALSE
      fs_b <- fs_se <- fs_p <- fs_F <- NA_real_; rf_b <- rf_se <- rf_p <- NA_real_
      fs_n <- rf_n <- NA_integer_
      same_rows <- same_cols <- same_vals <- same_terms <- NA
      conv <- NA; sepinfo <- NULL
      if(n < npar + 2L) stops <- c(stops, "too few complete cases to fit the prespecified model")
      if(binary && length(stops) == 0L && (nev == 0L || nev == n)) stops <- c(stops, "no variation in the outcome")

      if(length(stops) == 0L){
        fitted_both <- TRUE
        m_fs <- lm(as.formula(paste("bp ~", rhs)), data=Dc)
        gg   <- if(binary) fit_glm(as.formula(paste("y ~", rhs)), Dc) else
                           list(m=lm(as.formula(paste("y ~", rhs)), data=Dc), warn=character(0))
        m_rf <- gg$m; wmsg <- gg$warn

        mf1 <- model.frame(m_fs); mf2 <- model.frame(m_rf)
        X1  <- model.matrix(m_fs); X2 <- model.matrix(m_rf)
        same_rows  <- identical(rownames(mf1), rownames(mf2))
        same_cols  <- identical(colnames(X1), colnames(X2))
        same_vals  <- same_cols && identical(dim(X1), dim(X2)) && max(abs(X1 - X2)) == 0
        same_terms <- identical(attr(terms(m_fs), "term.labels"), attr(terms(m_rf), "term.labels"))
        fs_n <- nrow(mf1); rf_n <- nrow(mf2)
        s1 <- summary(m_fs)$coefficients
        if(any(is.na(coef(m_fs))) || !("z" %in% rownames(s1))) stops <- c(stops, "first stage is rank-deficient")
        else { fs_b <- s1["z",1]; fs_se <- s1["z",2]; fs_p <- s1["z",4]; fs_F <- (fs_b/fs_se)^2 }
        s2 <- tryCatch(summary(m_rf)$coefficients, error=function(e) NULL)
        if(!is.null(s2) && "z" %in% rownames(s2)){ rf_b <- s2["z",1]; rf_se <- s2["z",2]; rf_p <- s2["z",4] }
        if(binary){
          sepinfo <- logit_status(m_rf, wmsg)
          conv <- sepinfo$converged && !sepinfo$separation && !zero_ev_level
          if(!sepinfo$converged) stops <- c(stops, "logistic model did not converge")
          if(sepinfo$separation || zero_ev_level) stops <- c(stops, "logistic separation")
        } else {
          conv <- !any(is.na(coef(m_rf))) && all(is.finite(coef(m_rf)))
          if(!conv) stops <- c(stops, "reduced-form OLS is rank-deficient")
        }

        if(is.finite(fs_b) && !(fs_b > 0)) stops <- c(stops, "first-stage PGS coefficient is not positive (orientation)")
        if(!(isTRUE(same_rows) && isTRUE(same_vals) && isTRUE(same_terms)))
          stops <- c(stops, "the two stages were not fitted in the same women or the same specification")
      }
      stop_orientation <- any(grepl("orientation", stops))
      stop_separation  <- any(grepl("separation",  stops))
      stop_convergence <- any(grepl("did not converge", stops))
      status <- if(length(stops)) paste0("stopped: ", jn(unique(stops))) else "ok"
      if(zero_ev_level) notes <- c(notes, "a technology level in the model has no events or no non-events")
      if(tech_exc) notes <- c(notes, sprintf(
        "prespecified two-cell exception: genotyping technology recoded to %s, identically in both stages, in the same women",
        tech_code))
      if(!tech_in) notes <- c(notes, sprintf("only one technology level present (%s); the term was dropped", lv[1]))
      if(length(wmsg)) notes <- c(notes, unique(wmsg))

      th <- th_se <- e10 <- se10 <- or10 <- or_lo <- or_hi <- g10 <- g_lo <- g_hi <- wp <- NA_real_
      if(status == "ok" && is.finite(fs_b) && fs_b > 0 && is.finite(fs_se) && is.finite(rf_b) && is.finite(rf_se)){
        th    <- rf_b / fs_b
        th_se <- sqrt(rf_se^2 / fs_b^2 + rf_b^2 * fs_se^2 / fs_b^4)
        e10   <- SCALE * th
        se10  <- SCALE * th_se
        wp    <- 2 * pnorm(-abs(th / th_se))
        if(binary){ or10 <- exp(e10); or_lo <- exp(e10 - ZCRIT*se10); or_hi <- exp(e10 + ZCRIT*se10) }
        else      { g10  <- e10;      g_lo  <- e10 - ZCRIT*se10;      g_hi  <- e10 + ZCRIT*se10 }
      }
      disp <- if(!is.finite(e10)) "--" else if(binary)
                sprintf("%.2f (%.2f, %.2f)", or10, or_lo, or_hi) else
                sprintf("%.1f (%.1f, %.1f)", g10, g_lo, g_hi)
      pdisp <- if(!is.finite(wp)) "--" else if(wp < 1e-3) formatC(wp, format="e", digits=1) else sprintf("%.3f", wp)

      RES[[length(RES) + 1]] <- data.table(
        cohort=coh, cohort_display=unname(DISPLAY[coh]), site=unname(SITE_SHORT[coh]),
        ancestry_group=unname(GROUP[coh]), bp_trait=tr, exposure=unname(BPCOL[[tr]]), pgs_id=sid,
        outcome=oc, outcome_label=osp$label, outcome_type=osp$type, denominator=osp$denom,
        N=n, n_events=nev,
        fs_beta_mmHg_per_SD=fs_b, fs_se=fs_se, fs_p=fs_p, fs_F=fs_F, weak_instrument=isTRUE(fs_F < WEAK_F),
        rf_beta_per_SD=rf_b, rf_se=rf_se, rf_p=rf_p,
        wald_theta_per_mmHg=th, wald_se_theta_per_mmHg=th_se,
        wald_estimate_per_10mmHg=e10, wald_se_per_10mmHg=se10,
        or_per_10mmHg=or10, or_lo95=or_lo, or_hi95=or_hi,
        grams_per_10mmHg=g10, grams_lo95=g_lo, grams_hi95=g_hi,
        wald_p=wp, estimate_display=disp, wald_p_display=pdisp,
        model_convergence=if(is.na(conv)) "not fitted" else if(isTRUE(conv)) "converged" else "not converged",
        status=status, fitted_both_stages=fitted_both, stop_orientation=stop_orientation,
        stop_separation=stop_separation, stop_convergence=stop_convergence,
        n_technology_levels=length(lv), technology_levels=jn(lv), technology_reference=ref,
        technology_term_in_model=tech_in, notes=jn(notes),
        technology_exception_applied=tech_exc, technology_coding=tech_code)

      FS[[length(FS) + 1]] <- data.table(
        cohort=coh, cohort_display=unname(DISPLAY[coh]), ancestry_group=unname(GROUP[coh]),
        bp_trait=tr, exposure=unname(BPCOL[[tr]]), pgs_id=sid, outcome=oc, outcome_type=osp$type,
        denominator=osp$denom, N=n, n_model_frame=fs_n, model_rhs=rhs,
        fs_beta_mmHg_per_SD=fs_b, fs_se=fs_se, fs_p=fs_p, fs_F=fs_F, weak_instrument=isTRUE(fs_F < WEAK_F),
        bp_mean=if(n) mean(Dc$bp) else NA_real_, bp_sd=if(n > 1) sd(Dc$bp) else NA_real_,
        z_mean=if(n) mean(Dc$z) else NA_real_, z_sd=if(n > 1) sd(Dc$z) else NA_real_,
        age_mean=if(n) mean(Dc$AGE) else NA_real_,
        n_technology_levels=length(lv), technology_levels=jn(lv), technology_reference=ref,
        technology_term_in_model=tech_in,
        n_GSA_only=as.integer(tabt_raw[["GSA only"]]),
        n_lowpass_WGS_only=as.integer(tabt_raw[["low-pass WGS only"]]),
        n_both=as.integer(tabt_raw[["both"]]),
        technology_counts=jn(sprintf("%s=%d", lv, as.integer(tabt[lv]))), status=status,
        technology_exception_applied=tech_exc, technology_coding=tech_code,
        raw_technology_counts=jn(sprintf("%s=%d", TECH_ALL, as.integer(tabt_raw))))

      RF[[length(RF) + 1]] <- data.table(
        cohort=coh, cohort_display=unname(DISPLAY[coh]), ancestry_group=unname(GROUP[coh]),
        bp_trait=tr, pgs_id=sid, outcome=oc, outcome_label=osp$label, outcome_type=osp$type,
        model=if(binary) "logistic (binomial, logit link)" else "OLS (grams)",
        N=n, n_model_frame=rf_n, n_events=nev, model_rhs=rhs,
        rf_beta_per_SD=rf_b, rf_se=rf_se, rf_p=rf_p,
        converged=if(is.null(sepinfo)) conv else sepinfo$converged,
        separation=if(is.null(sepinfo)) NA else sepinfo$separation,
        sep_fitted_at_bound=if(is.null(sepinfo)) NA else sepinfo$sep_fitted,
        sep_large_scaled_coefficient=if(is.null(sepinfo)) NA else sepinfo$sep_coef,
        sep_large_scaled_se=if(is.null(sepinfo)) NA else sepinfo$sep_se,
        sep_fit_warning=if(is.null(sepinfo)) NA else sepinfo$sep_warning,
        zero_event_technology_level=zero_ev_level,
        glm_iterations=if(is.null(sepinfo)) NA_integer_ else as.integer(sepinfo$iter),
        min_fitted=if(is.null(sepinfo)) NA_real_ else sepinfo$min_fitted,
        max_fitted=if(is.null(sepinfo)) NA_real_ else sepinfo$max_fitted,
        max_scaled_coefficient=if(is.null(sepinfo)) NA_real_ else sepinfo$max_scaled_coef,
        max_scaled_coefficient_se=if(is.null(sepinfo)) NA_real_ else sepinfo$max_scaled_se,
        events_by_technology_level=if(binary) jn(sprintf("%s: %d/%d",
              lvcnt[in_model_levels == TRUE & n > 0]$technology_level,
              lvcnt[in_model_levels == TRUE & n > 0]$n_events,
              lvcnt[in_model_levels == TRUE & n > 0]$n)) else "",
        events_by_raw_technology_level=if(binary) jn(sprintf("%s: %d/%d", TECH_ALL[tabt_raw > 0],
              vapply(TECH_ALL[tabt_raw > 0], function(L) sum(Dc$y[tech_raw == L] == 1L), 0L),
              as.integer(tabt_raw[tabt_raw > 0]))) else "",
        same_sample_as_first_stage=same_rows, same_design_as_first_stage=same_vals,
        same_terms_as_first_stage=same_terms, status=status, warnings=jn(unique(wmsg)),
        technology_exception_applied=tech_exc, technology_coding=tech_code)
    }
  }
}
RES <- rbindlist(RES); FS <- rbindlist(FS); RF <- rbindlist(RF)
TL  <- rbindlist(TL); ZCHK <- rbindlist(ZCHK)

W_full(RES, "cohort_mr_results.tsv")
W_full(FS,  "first_stage_by_outcome.tsv")
W_full(RF,  "reduced_form_results.tsv")
W_full(TL,  "technology_levels.tsv")
cat(sprintf("04: %d cells fitted; %d reached a Wald estimate\n",
            nrow(RES), sum(RES$status == "ok")))

## ---- pooling -----------------------------------------------------------------

ex <- function(x) vapply(x, function(v){
  if(is.na(v)) return(NA_character_); if(!is.finite(v)) return(as.character(v))
  for(d in 1:17){ s <- sprintf(paste0("%.", d, "g"), v); if(identical(as.numeric(s), v)) return(s) }
  sprintf("%.17g", v) }, character(1), USE.NAMES=FALSE)
exact <- function(dt){ dt <- copy(dt); for(cn in names(dt)) if(is.double(dt[[cn]])) set(dt, j=cn, value=ex(dt[[cn]])); dt }
W_ <- function(dt, f) fwrite(exact(dt), file.path(RESULTS, f), sep="\t")

R <- fread(file.path(RESULTS, "cohort_mr_results.tsv"))
stopifnot(nrow(R) == 40L)
if (any(R$status != "ok"))
  stop("the cohort results contain ", sum(R$status != "ok"), " stopped cell(s)")

COH_ORDER <- c("AMANHI-Bangladesh","AMANHI-Pakistan","GAPPS-Bangladesh","AMANHI-Pemba","GAPPS-Zambia")
SA  <- c("AMANHI-Bangladesh","AMANHI-Pakistan","GAPPS-Bangladesh")
AFR <- c("AMANHI-Pemba","GAPPS-Zambia")
ANALYSES <- data.table(bp_trait = rep(c("SBP","DBP"), each = 4),
                       outcome  = rep(c("PTB","LBW","SGA","BWT"), 2))
OUTLAB <- c(PTB="Preterm birth", LBW="Low birth weight (<2,500 g)",
            SGA="Small for gestational age (<10th centile)", BWT="Birth weight, g")
SCALE  <- c(PTB="log odds ratio", LBW="log odds ratio", SGA="log odds ratio", BWT="grams")

IN <- R[, .(bp_trait, outcome, outcome_label, outcome_type, cohort, cohort_display, site, ancestry_group,
            pgs_id, N, n_events,
            fs_beta_mmHg_per_SD, fs_se, fs_p, fs_F,
            rf_beta_per_SD, rf_se, rf_p,
            yi = wald_estimate_per_10mmHg, sei = wald_se_per_10mmHg,
            or_per_10mmHg, or_lo95, or_hi95, grams_per_10mmHg, grams_lo95, grams_hi95,
            wald_p, analysis_scale = SCALE[outcome], status, model_convergence)]
IN[, cohort := factor(cohort, levels = COH_ORDER)]
setorderv(IN, c("bp_trait","outcome","cohort")); IN[, cohort := as.character(cohort)]
if(any(!is.finite(IN$yi)) || any(!is.finite(IN$sei)) || any(IN$sei <= 0))
  stop("a Wald estimate or SE is missing, infinite or non-positive")
W_(IN, "pooling_input.tsv")

GROUPS <- list(`South Asian` = SA, `African` = AFR, `Overall` = COH_ORDER)

POOL <- rbindlist(lapply(seq_len(nrow(ANALYSES)), function(i){
  tr <- ANALYSES$bp_trait[i]; oc <- ANALYSES$outcome[i]; bin <- oc != "BWT"
  rbindlist(lapply(names(GROUPS), function(g){
    D <- IN[bp_trait == tr & outcome == oc & cohort %in% GROUPS[[g]]]
    if(nrow(D) != length(GROUPS[[g]]))
      stop("group ", g, " for ", tr, " -> ", oc, " has ", nrow(D), " cohorts, expected ", length(GROUPS[[g]]))

    m <- rma(yi = D$yi, sei = D$sei, method = "REML", test = "z")
    data.table(
      bp_trait = tr, outcome = oc, outcome_label = OUTLAB[[oc]],
      outcome_type = if(bin) "binary" else "continuous", analysis_scale = SCALE[[oc]],
      population_group = g, model = "random effects (REML), z test",
      k_cohorts = m$k, N_total = sum(D$N),
      events_total = if(bin) sum(as.numeric(D$n_events)) else NA_real_,
      pooled_estimate = as.numeric(m$beta), pooled_se = m$se, ci_lo = m$ci.lb, ci_hi = m$ci.ub,
      or_per_10mmHg = if(bin) exp(as.numeric(m$beta)) else NA_real_,
      or_lo95 = if(bin) exp(m$ci.lb) else NA_real_, or_hi95 = if(bin) exp(m$ci.ub) else NA_real_,
      grams_per_10mmHg = if(bin) NA_real_ else as.numeric(m$beta),
      grams_lo95 = if(bin) NA_real_ else m$ci.lb, grams_hi95 = if(bin) NA_real_ else m$ci.ub,
      z = m$zval, p = m$pval, tau2 = m$tau2, I2 = m$I2, H2 = m$H2,
      Q = m$QE, Q_df = m$k - 1L, Q_p = m$QEp,
      cohorts = paste(D$cohort_display, collapse = "; "))
  }))
}))
W_(POOL, "ancestry_pooled_mr_results.tsv")

cat(sprintf("04: %d pooled rows (%d analyses x 3 population groups)\n",
            nrow(POOL), nrow(ANALYSES)))
