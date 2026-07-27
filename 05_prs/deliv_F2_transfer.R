#!/usr/bin/env Rscript
# ============================================================
# deliv_F2_transfer.R  [ B11 -> Figure 2, transferability heatmap ]
# Heatmap of incremental R2(%) for each score ancestry x cohort, mean BP, platforms
# merged, faceted by trait (SBP/DBP). Shows EUR/SAS strong in SAS cohorts, weak in AFR,
# MVP/EAS poor everywhere. A pivot of transfer_grid (mean def).
#   Writes: figures/F2_transfer.{png,pdf}
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

momi_deliverable("F2_transfer", script="05_prs/deliv_F2_transfer.R",
                 inputs="transfer_grid", P=P, body=function(ctx){
  G <- momi_read_intermediate("transfer_grid", P)
  ## headline exposure from config, not a literal (see note in deliv_T2_transfer.R)
  M <- G[definition==MOMI_BPDEF_DEFAULT]
  M[, coh_disp := MOMI_DISPLAY[cohort]]
  # order cohorts SAS then AFR; score ancestries EUR/SAS/EAS/MVP
  coh_order <- MOMI_DISPLAY[c(MOMI_COH_SAS, MOMI_COH_AFR)]
  M[, coh_disp := factor(coh_disp, levels=coh_order)]
  M[, score_anc := factor(score_anc, levels=c("EUR","SAS","EAS","MVP"))]
  M[, trait := factor(trait, levels=c("SBP","DBP"))]
  M <- M[!is.na(score_anc)]

  g <- ggplot(M, aes(coh_disp, score_anc, fill=R2pct)) +
    geom_tile(colour="white", linewidth=0.6) +
    geom_text(aes(label=sprintf("%.1f", R2pct)), size=3) +
    facet_wrap(~trait) +
    scale_fill_gradient(low="#f7fbff", high="#08519c", name="R² (%)") +
    labs(x=NULL, y="PRS ancestry",
         title="BP-PRS transferability (mean BP, incremental R² %)",
         subtitle="EUR/SAS strong in South-Asian cohorts; weak in African; MVP/EAS poor everywhere") +
    theme_minimal(base_size=11) +
    theme(axis.text.x=element_text(angle=35, hjust=1),
          panel.grid=element_blank(), plot.subtitle=element_text(size=8.5))
  outs <- momi_save_fig(g, "F2_transfer", width=9, height=4.4, P=P)
  list(n=nrow(M), key=sprintf("cells=%d max R2=%.1f%%", nrow(M), max(M$R2pct,na.rm=TRUE)),
       outputs=basename(outs))
})
