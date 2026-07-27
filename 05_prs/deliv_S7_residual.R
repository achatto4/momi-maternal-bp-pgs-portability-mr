#!/usr/bin/env Rscript
# ============================================================
# deliv_S7_residual.R  [ B13 -> Supp Table S7, GA-residual transferability ]
# For the chosen instrument, R2 of the GA-standardised residual definition vs mean BP,
# per cohort x trait. View of transfer_grid (definitions mean & resid).
#
# ------------------------------------------------------------------
# WHAT THIS TABLE DOES AND DOES NOT ESTABLISH -- rewritten 2026-07-20.
#
# It has been framed two ways, and BOTH WERE WRONG:
#   * before 2026-07-19: "GA-standardisation adds nothing" -- used to dismiss `resid`;
#   * after the re-tiering:  "it justifies the headline exposure" -- used to promote `resid`.
#
# The numbers support neither. Across all ten cohort x trait cells the largest change in R2
# is 0.20 PERCENTAGE POINTS (PreSSMat SBP, 6.258 -> 6.459); most move by under 0.05. On the
# scale of the finding this table exists alongside -- a 2.6-fold portability gap between
# cohorts -- that is nothing.
#
# THE CORRECT READING: transferability is INSENSITIVE to whether blood pressure is
# GA-standardised. That makes S7 a ROBUSTNESS check, not an argument for the exposure. The
# choice of `resid` as primary rests on a conceptual claim -- that removing the gestational
# trajectory gives a cleaner measure of a mother's underlying BP -- and this table's job is
# to show the headline result does not depend on that claim being right.
#
# ONE CURIOSITY, DELIBERATELY NOT INTERPRETED: the direction splits perfectly by ancestry.
# GA-standardisation raises R2 in all six South Asian cells and lowers it in all four African
# ones. With a maximum effect of 0.2 percentage points this is far too small to build on, and
# it is recorded here only so that a future reader who notices the pattern knows it was seen
# and judged uninterpretable rather than missed.
# ------------------------------------------------------------------
#   Writes: tables/S7_residual.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

momi_deliverable("S7_residual", script="05_prs/deliv_S7_residual.R",
                 inputs="transfer_grid", P=P, body=function(ctx){
  G <- momi_read_intermediate("transfer_grid", P)
  rows <- list()
  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    sid <- momi_instrument(coh,tr)
    rm <- G[cohort==coh & trait==tr & score_id==sid & definition=="mean"]
    rr <- G[cohort==coh & trait==tr & score_id==sid & definition=="resid"]
    if(!nrow(rm)) next
    rows[[paste(coh,tr)]] <- data.table(
      cohort_display=MOMI_DISPLAY[coh], ancestry=MOMI_ANC[coh], trait=tr, chosen_PGS=sid,
      R2pct_mean=rm$R2pct, R2pct_resid=if(nrow(rr)) rr$R2pct else NA_real_,
      N=rm$n)
  }
  T <- rbindlist(rows); setorder(T, ancestry, cohort_display, trait)
  out <- momi_write_table(T, "S7_residual", P)
  list(n=nrow(T), key="GA-residual R² vs mean-BP R² (chosen instrument)", outputs=basename(out))
})
