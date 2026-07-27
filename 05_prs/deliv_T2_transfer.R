#!/usr/bin/env Rscript
# ============================================================
# deliv_T2_transfer.R  [ B10 -> Table 2, transferability of the CHOSEN instrument ]
# The Part-I core table. Exposure = MOMI_BPDEF_DEFAULT. For each cohort x trait, the chosen
# PRS (EUR primary; SAS for GAPPS-Bangladesh, via momi_instrument) and its OBSERVED
# incremental R2(%), first-stage F, and N. No cross-validation column.
# The measurement-error-corrected R2 is deliberately NOT here — see the note at R2pct below.
# A filtered view of transfer_grid — identical numbers to S4/Fig2 by construction.
#   Writes: tables/T2_transfer.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

momi_deliverable("T2_transfer", script="05_prs/deliv_T2_transfer.R",
                 inputs="transfer_grid", P=P, body=function(ctx){
  G <- momi_read_intermediate("transfer_grid", P)
  ## Read the headline exposure from config, never a literal. Tiering changed 2026-07-19
  ## (primary = resid, then ge20; mean demoted to auxiliary) and a hardcoded "mean" here
  ## would have silently kept Table 2 on the old exposure while config said otherwise.
  M <- G[definition==MOMI_BPDEF_DEFAULT]
  rows <- list()
  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    sid <- momi_instrument(coh, tr)
    r <- M[cohort==coh & trait==tr & score_id==sid]
    if(!nrow(r)) next
    rows[[paste(coh,tr)]] <- data.table(
      cohort_display=MOMI_DISPLAY[coh], ancestry=MOMI_ANC[coh], trait=tr,
      chosen_instrument=MOMI_INSTR_ANC[[coh]], chosen_PGS=sid,
      ## DECISION 2026-07-19 (user): the disattenuated R2 is a SUPPLEMENTARY quantity and
      ## does not appear in Table 2. Table 2 reports only what was measured.
      ##
      ## Rationale: R2pct_disatt is not an observation, it is an inference that depends on a
      ## reliability model (ICC + Spearman-Brown on k readings) whose k differs across
      ## cohorts by a factor of ~1.8 (2.51 in AMANHI-Karachi vs 4.44 in PreSSMat). Since the
      ## paper's central claim is a CROSS-COHORT CONTRAST in portability, putting a
      ## cohort-specific model correction into the headline table would make the contrast
      ## partly an artefact of how often each site happened to measure BP. The observed R2 is
      ## the like-for-like comparison. The correction, with all its inputs, lives in S4.
      R2pct=r$R2pct,
      F=round(r$F,1), N=r$n)
  }
  T2 <- rbindlist(rows)
  setorder(T2, ancestry, cohort_display, trait)
  out <- momi_write_table(T2, "T2_transfer", P)
  list(n=nrow(T2),
       key=sprintf("SAS R2 range %.1f-%.1f%%; AFR %.1f-%.1f%%",
                   min(T2[ancestry=="SAS",R2pct]), max(T2[ancestry=="SAS",R2pct]),
                   min(T2[ancestry=="AFR",R2pct]), max(T2[ancestry=="AFR",R2pct])),
       outputs=basename(out))
})
