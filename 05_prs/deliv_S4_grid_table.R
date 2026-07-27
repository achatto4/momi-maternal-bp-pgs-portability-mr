#!/usr/bin/env Rscript
# ============================================================
# deliv_S4_grid_table.R  [ B09 -> Supp Table S4, full transferability grid ]
# Reader-facing wide table: R2 (%) for every score x cohort x BP definition, with F
# and N alongside. A pivot of the frozen transfer_grid (no re-estimation).
#   Writes: tables/S4_transferability.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

momi_deliverable("S4_transferability", script="05_prs/deliv_S4_grid_table.R",
                 inputs="transfer_grid", P=P, body=function(ctx){
  G <- momi_read_intermediate("transfer_grid", P)
  defs <- intersect(c("early","mean","mn2","median","last","tri1","tri2","tri3","lt20","ge20","resid"),
                    unique(G$definition))
  W <- dcast(G, score_id+score_anc+trait+cohort+cohort_anc ~ definition, value.var="R2pct")
  setcolorder(W, c("score_id","score_anc","trait","cohort","cohort_anc", defs))
  setorder(W, trait, cohort_anc, cohort, -mean)
  out <- momi_write_table(W, "S4_transferability", P)
  list(n=nrow(W), key=sprintf("rows=%d defs=%d", nrow(W), length(defs)), outputs=basename(out))
})
