#!/usr/bin/env Rscript
# ============================================================
# deliv_S8_wk20.R  [ B14 -> Supp Table S8, <20 vs >=20 wk transferability ]
# For the chosen instrument, R2 of the chronic proxy (BP <20 wk) vs gestational (>=20 wk),
# per cohort x trait, WITH N so sparse cells are visible (AMANHI has ~no <20 wk BP, so
# lt20 is only meaningful for GAPPS-Bangladesh & Zambia). View of transfer_grid.
#   Writes: tables/S8_wk20.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

momi_deliverable("S8_wk20", script="05_prs/deliv_S8_wk20.R",
                 inputs="transfer_grid", P=P, body=function(ctx){
  G <- momi_read_intermediate("transfer_grid", P)
  rows <- list()
  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    sid <- momi_instrument(coh,tr)
    lt <- G[cohort==coh & trait==tr & score_id==sid & definition=="lt20"]
    ge <- G[cohort==coh & trait==tr & score_id==sid & definition=="ge20"]
    rows[[paste(coh,tr)]] <- data.table(
      cohort_display=MOMI_DISPLAY[coh], ancestry=MOMI_ANC[coh], trait=tr, chosen_PGS=sid,
      R2pct_lt20=if(nrow(lt)) lt$R2pct else NA, N_lt20=if(nrow(lt)) lt$n else NA_integer_,
      R2pct_ge20=if(nrow(ge)) ge$R2pct else NA, N_ge20=if(nrow(ge)) ge$n else NA_integer_,
      note=if(nrow(lt) && lt$n < 200) "sparse <20wk (interpret with caution)" else "")
  }
  T <- rbindlist(rows); setorder(T, ancestry, cohort_display, trait)
  out <- momi_write_table(T, "S8_wk20", P)
  list(n=nrow(T),
       key="chronic(<20wk) vs gestational(>=20wk) R²; usable mainly for GAPPS-B & Zambia",
       outputs=basename(out))
})
