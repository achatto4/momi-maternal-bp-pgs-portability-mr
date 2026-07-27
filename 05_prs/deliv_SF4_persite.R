#!/usr/bin/env Rscript
# ============================================================
# deliv_SF4_persite.R  [ B27 -> SF4, per-cohort MR forests ]
#
# ------------------------------------------------------------------
# THIS FIGURE CARRIES THE HONESTY THAT S11 USED TO CARRY.
#
# The leave-one-out sensitivity analysis (S11/B23) was removed from the paper on 2026-07-20
# as too much apparatus for what it conveyed. But the CONCLUSIONS it produced still govern
# how results are described, and this figure is why that is defensible: a forest showing
# every cohort's estimate and interval makes the same point directly, without asking the
# reader to understand what leave-one-out means.
#
# Concretely, a reader looking at these panels can see for themselves that:
#   * for INDICATED PTB, GAPPS-Bangladesh's interval is the only one excluding 1 -- which is
#     why the text calls it a Matlab result rather than a five-cohort effect;
#   * for BIRTHWEIGHT, four of five cohort estimates sit below zero and the fifth
#     (GAPPS-Zambia) has an interval so wide it is uninformative -- which is why the pooled
#     estimate is reported as robust;
#   * for SBP -> LBW, AMANHI-Pakistan is doing the work.
# None of that needs a sensitivity framework. It needs the intervals drawn.
#
# DESIGN DECISIONS THAT MATTER FOR HOW THIS IS READ
#
#   CASE COUNTS ARE PRINTED ON EVERY ROW. Without them a wide interval looks like a weak
#   result rather than an absent outcome, and the two are constantly confused in this
#   dataset -- AMANHI-Pemba has 18 indicated preterm births. The n on the row is what stops
#   a reader concluding "no effect at Pemba" when the honest statement is "not estimable".
#
#   THE FIRST-STAGE F IS PRINTED TOO. Part I established the instrument does not transfer to
#   the African cohorts. A reader seeing GAPPS-Zambia's F of 17 next to its wide interval can
#   tell that the imprecision is an instrument problem, not evidence of a smaller effect.
#
#   SAS AND AFR SUBTOTALS ARE SHOWN SEPARATELY, plus the overall. Anagh's call 2026-07-20:
#   the strata may reflect GENUINE CAUSAL HETEROGENEITY, not merely differing precision, so
#   both are displayed rather than one being designated a nuisance. The overall diamond
#   carries a heterogeneity annotation.
#
#   NO POOLED DIAMOND IS DRAWN FOR A STRATUM WITH k=1. A single cohort is not a
#   meta-analysis and must not be drawn as one.
# ------------------------------------------------------------------
#   Writes: figures/SF4_forest_mr.{png,pdf}, figures/SF4b_forest_subtype.{png,pdf}
#           tables/SF4_forest_data.tsv
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

## pooled row for a set of per-cohort estimates; NULL if fewer than 2 (see header)
pooled_row <- function(s, label){
  s <- s[is.finite(theta_per10) & is.finite(se_per10) & se_per10 > 0]
  if(nrow(s) < 2) return(NULL)
  fe <- momi_meta_iv(s$theta_per10, s$se_per10)
  w  <- 1/s$se_per10^2
  Q  <- sum(w*(s$theta_per10-fe$b)^2); dfQ <- nrow(s)-1
  data.table(label=label, is_pooled=TRUE, k=fe$k,
             theta=fe$b, lo=fe$b-1.96*fe$se, hi=fe$b+1.96*fe$se,
             N=sum(s$N), n_case=sum(s$n_case, na.rm=TRUE), F=NA_real_,
             I2=max(0, 100*(Q-dfQ)/Q), p=fe$p)
}

momi_deliverable("SF4_persite", script="05_prs/deliv_SF4_persite.R",
                 inputs="mr_panel;ptb_subtype", P=P, stop_on_error=FALSE,
                 body=function(ctx){

  MR <- tryCatch(momi_read_intermediate("mr_panel", P), error=function(e) NULL)
  if(is.null(MR) || !nrow(MR))
    return(list(skip=TRUE, reason="no mr_panel.rds -- run B22 first"))

  ## ---- assemble one long frame: cohorts, then SAS / AFR / overall summaries ----
  build_rows <- function(D, panel_name){
    out <- list()
    for(tr in c("SBP","DBP")) for(oc in unique(D$outcome)){
      s <- D[trait==tr & outcome==oc & is.finite(theta_per10) & se_per10 > 0]
      if(!nrow(s)) next
      s[, anc := MOMI_ANC[cohort]]
      coh <- s[, .(label = MOMI_DISPLAY[cohort], is_pooled = FALSE, k = 1L,
                   theta = theta_per10, lo = theta_per10-1.96*se_per10,
                   hi = theta_per10+1.96*se_per10,
                   N, n_case, F, I2 = NA_real_, p = rf_p, anc)]
      setorder(coh, anc, label)
      rows <- rbindlist(list(
        coh[, -"anc"],
        pooled_row(s[anc=="SAS"], "  South Asian (pooled)"),
        pooled_row(s[anc=="AFR"], "  African (pooled)"),
        pooled_row(s,             "  All cohorts (pooled)")), fill=TRUE)
      out[[paste(tr,oc)]] <- cbind(data.table(panel=panel_name, trait=tr, outcome=oc), rows)
    }
    rbindlist(out, fill=TRUE)
  }

  FD <- build_rows(MR, "perinatal")

  SUBP <- tryCatch(momi_read_intermediate("ptb_subtype", P), error=function(e) NULL)
  if(!is.null(SUBP) && nrow(SUBP)){
    sp <- copy(SUBP); setnames(sp, "subtype", "outcome")
    FD <- rbind(FD, build_rows(sp, "PTB subtype"), fill=TRUE)
  }
  if(!nrow(FD)) return(list(skip=TRUE, reason="nothing estimable to plot"))

  ## annotation strings drawn beside each row -- the point of the figure (see header)
  FD[, ann := fifelse(is_pooled,
                      fifelse(is.finite(I2), sprintf("k=%d, I2=%.0f%%", k, I2), sprintf("k=%d", k)),
                      fifelse(is.finite(n_case),
                              sprintf("n=%d, %d events, F=%.0f", N, n_case, F),
                              sprintf("n=%d, F=%.0f", N, F)))]
  out_t <- momi_write_table(FD, "SF4_forest_data", P)

  ## ---- plotting helper. Binary outcomes are drawn on the OR (exp) scale, continuous on
  ## the natural scale, because a log axis for grams is meaningless. ----
  draw <- function(D, id, title, xlab, expo){
    if(!nrow(D)) return(character(0))
    D <- copy(D)
    if(expo) D[, `:=`(theta=exp(theta), lo=exp(lo), hi=exp(hi))]
    D[, row := factor(paste(outcome, label, sep=" | "),
                      levels=rev(unique(paste(outcome, label, sep=" | "))))]
    null_at <- if(expo) 1 else 0
    g <- ggplot(D, aes(x=theta, y=row)) +
      geom_vline(xintercept=null_at, linetype="dashed", linewidth=.3, colour="grey40") +
      geom_errorbarh(aes(xmin=lo, xmax=hi), height=.18, linewidth=.35) +
      geom_point(aes(shape=is_pooled, size=is_pooled)) +
      geom_text(aes(label=ann), x=Inf, hjust=1.02, size=2.1, colour="grey30") +
      scale_shape_manual(values=c(`FALSE`=16, `TRUE`=18), guide="none") +
      scale_size_manual(values=c(`FALSE`=1.6, `TRUE`=2.8), guide="none") +
      facet_grid(trait ~ ., scales="free_y", space="free_y") +
      labs(x=xlab, y=NULL, title=title,
           subtitle="Diamonds are pooled; circles are individual cohorts. Event counts and first-stage F shown at right.") +
      theme_minimal(base_size=8) +
      theme(panel.grid.minor=element_blank(),
            panel.grid.major.y=element_blank(),
            plot.subtitle=element_text(size=6.5, colour="grey35"))
    if(expo) g <- g + scale_x_continuous(trans="log10")
    momi_save_fig(g, id, width=8.5, height=max(4, 0.22*nrow(D)+2), P=P)
  }

  f1 <- draw(FD[panel=="perinatal" & outcome != "BWT"], "SF4_forest_mr",
             "Mendelian randomization: maternal blood pressure and perinatal outcomes",
             "OR per 10 mmHg (log scale)", expo=TRUE)
  f2 <- draw(FD[panel=="perinatal" & outcome == "BWT"], "SF4c_forest_bwt",
             "Mendelian randomization: maternal blood pressure and birthweight",
             "grams per 10 mmHg", expo=FALSE)
  f3 <- draw(FD[panel=="PTB subtype"], "SF4b_forest_subtype",
             "Preterm birth by subtype: spontaneous vs provider-initiated",
             "OR per 10 mmHg (log scale)", expo=TRUE)

  ## ---- console ----
  cat("\n=== rows drawn ===\n")
  print(FD[, .N, by=.(panel, outcome)][order(panel, outcome)])
  cat("\n=== cohort intervals excluding the null, by outcome ===\n")
  ex <- FD[is_pooled==FALSE & ((lo>0 & hi>0) | (lo<0 & hi<0)),
           .(panel, trait, outcome, label, theta=round(theta,3), n_case)]
  if(nrow(ex)) print(ex) else cat("  none\n")
  cat("\nThis is the list a reader reconstructs visually. If one cohort dominates it for a\n",
      "given outcome, the text must describe that outcome as that cohort's result.\n", sep="")

  list(n=nrow(FD),
       key=sprintf("rows=%d; panels=%d; figures=%d",
                   nrow(FD), uniqueN(FD$panel),
                   length(c(f1,f2,f3))/2),
       outputs=c(basename(out_t), basename(c(f1,f2,f3))))
})
