#!/usr/bin/env Rscript
# ============================================================
# deliv_S5_platform.R  [ B12 -> Supp Table S5, transferability by platform ]
# GSA vs lpWGS-dosage R2 for the chosen instrument (mean BP), per cohort x trait, with N.
# Shows dosage stability and small-N GSA caveats. Computed from prs_z_platform (per-platform
# z) + analytic_mothers via the shared estimator (same R2 definition as the merged grid).
#   Writes: tables/S5_platform.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

momi_deliverable("S5_platform", script="05_prs/deliv_S5_platform.R",
                 inputs="prs_z_platform;analytic_mothers", P=P, body=function(ctx){
  A   <- momi_read_intermediate("analytic_mothers", P)
  PLZ <- momi_read_intermediate("prs_z_platform", P)
  rows <- list()
  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    sid <- momi_instrument(coh,tr); bpcol <- paste0(substr(tr,1,1),"_mean")
    ac <- A[cohort==coh, c("IID","AGE",bpcol), with=FALSE]
    rr <- data.table(cohort_display=MOMI_DISPLAY[coh], ancestry=MOMI_ANC[coh],
                     trait=tr, chosen_PGS=sid)
    for(pl in MOMI_PLATS){
      z <- PLZ[score_id==sid & cohort==coh & platform==pl, .(IID,z)]
      tag <- if(pl=="gsa") "gsa" else "dosage"
      if(nrow(z)){
        d <- merge(z, ac, by="IID")
        tt <- momi_transfer(d$z, d[[bpcol]], cov=d[, .(AGE)])
        rr[[paste0("R2pct_",tag)]] <- round(100*tt$incR2,3); rr[[paste0("N_",tag)]] <- tt$n
      } else { rr[[paste0("R2pct_",tag)]] <- NA_real_; rr[[paste0("N_",tag)]] <- 0L }
    }
    rows[[paste(coh,tr)]] <- rr
  }
  T <- rbindlist(rows, fill=TRUE); setorder(T, ancestry, cohort_display, trait)
  out <- momi_write_table(T, "S5_platform", P)
  list(n=nrow(T), key="GSA vs lpWGS-dosage R² (chosen instrument, mean BP)", outputs=basename(out))
})
