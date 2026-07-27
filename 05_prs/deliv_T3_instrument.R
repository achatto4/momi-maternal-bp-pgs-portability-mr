#!/usr/bin/env Rscript
# ============================================================
# deliv_T3_instrument.R  [ B17 -> Table 3, instrument decision & first-stage strength ]
# The locked Part-I decision, written as a result. Rows = population group; columns =
# selected SBP PGS, selected DBP PGS, rationale, first-stage F range, R2 range.
# Derived from transfer_grid (mean def) + momi_instrument.
#   Writes: tables/T3_instrument.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

momi_deliverable("T3_instrument", script="05_prs/deliv_T3_instrument.R",
                 inputs="transfer_grid", P=P, body=function(ctx){
  ## headline exposure from config, not a literal (see note in deliv_T2_transfer.R)
  G <- momi_read_intermediate("transfer_grid", P)[definition==MOMI_BPDEF_DEFAULT]
  chosen <- function(cohs){
    rows <- rbindlist(lapply(cohs, function(coh) rbindlist(lapply(c("SBP","DBP"), function(tr){
      sid <- momi_instrument(coh,tr); r <- G[cohort==coh & trait==tr & score_id==sid]
      if(nrow(r)) data.table(coh=coh, tr=tr, sid=sid, F=r$F, R2=r$R2pct) }))))
    rows
  }
  grp <- function(name, cohs, rationale){
    r <- chosen(cohs)
    sbp <- unique(r[tr=="SBP", sid]); dbp <- unique(r[tr=="DBP", sid])
    data.table(group=name,
               cohorts=paste(vapply(cohs, function(c) MOMI_DISPLAY[[c]], ""), collapse="; "),
               selected_SBP_PGS=paste(sbp,collapse="/"), selected_DBP_PGS=paste(dbp,collapse="/"),
               rationale=rationale,
               F_range=sprintf("%.0f–%.0f", min(r$F), max(r$F)),
               R2_range_pct=sprintf("%.1f–%.1f", min(r$R2), max(r$R2)))
  }
  T3 <- rbindlist(list(
    grp("South-Asian (EUR primary)", c("AMANHI-Bangladesh","AMANHI-Pakistan"),
        "EUR PGS best/near-best; used for consistency across populations"),
    grp("South-Asian (GAPPS-Bangladesh)", "GAPPS-Bangladesh",
        "South-Asian PGS best-fitting here (early-BP cohort); SAS chosen"),
    grp("African (EUR primary)", c("AMANHI-Pemba","GAPPS-Zambia"),
        "EUR ≥ other panels but low R² — distinct genealogical branch; EUR for consistency")))
  out <- momi_write_table(T3, "T3_instrument", P)
  list(n=nrow(T3), key="EUR primary; SAS for GAPPS-Bangladesh (locked Part-I decision)",
       outputs=basename(out))
})
