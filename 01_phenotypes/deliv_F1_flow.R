#!/usr/bin/env Rscript
# ============================================================
# deliv_F1_flow.R  [ B07 -> Figure 1, Panel A (analytic sample flow) ]
# CONSORT-style cascade computed from analytic_mothers (PE kept as mediator):
#   first-preg mothers -> genotyped -> Part I (genotyped + valid antenatal BP)
#   -> Part II (+ PTB outcome), with per-cohort Part I counts and PE/twin descriptors.
# Emits the counts table, a LaTeX macro file (so the figure/manuscript numbers can
# never drift from the data), and the ggplot flow figure. (Panel B DAG is a static
# conceptual diagram assembled in LaTeX.)
#
# Writes: tables/F1_flow_counts.tsv, figures/flow_counts.tex, figures/F1_flow.{png,pdf}
#   Rscript deliv_F1_flow.R
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

momi_deliverable("F1_flow", script="01_phenotypes/deliv_F1_flow.R",
                 inputs="analytic_mothers", P=P, body=function(ctx){
  A <- momi_read_intermediate("analytic_mothers", P)
  A[, hasBP := is.finite(S_mean) | is.finite(D_mean)]

  ## pre-exclusion counts (enrolled first pregnancies, and non-singletons dropped in B01)
  FP <- tryCatch(momi_read_intermediate("flow_pre", P), error=function(e) NULL)
  n_enrolled <- if(!is.null(FP)) as.integer(FP$n_first_preg_all[1]) else nrow(A)
  n_nonsing  <- if(!is.null(FP)) as.integer(FP$n_nonsingleton[1])   else 0L

  N1     <- nrow(A)   # singleton first pregnancies (analytic_mothers is singleton-only)
  geno   <- A[genotyped==1]
  Ngeno  <- nrow(geno)
  partI  <- geno[hasBP==TRUE]
  NpartI <- nrow(partI)
  partII <- partI[!is.na(PTB)]
  NpartII<- nrow(partII)
  exGeno <- N1 - Ngeno; exBP <- Ngeno - NpartI; exPTB <- NpartI - NpartII

  ## Split the BP exclusion into its TWO real causes. Since B01's antenatal restriction
  ## (2026-07-19) "no valid antenatal BP" conflates two quite different women: those with no
  ## BP recorded at all, and those whose ONLY readings were postpartum. The second group is
  ## an artefact of AMANHI's follow-up protocol rather than missing data, and lumping them
  ## together hides why the Part I count moved from 13,896 to 13,717.
  exBPall  <- geno[hasBP==FALSE]
  hasPN    <- function(d) is.finite(d$S_mean_pn) | is.finite(d$D_mean_pn)
  exBP_pn  <- if(all(c("S_mean_pn","D_mean_pn") %in% names(A))) sum(hasPN(exBPall)) else NA_integer_
  exBP_non <- if(is.na(exBP_pn)) NA_integer_ else exBP - exBP_pn

  # per-cohort Part I, with consortium display names
  bc <- partI[, .(n=.N), by=cohort][order(-n)]
  bc[, display := MOMI_DISPLAY[cohort]]
  pe <- sum(partI$PE, na.rm=TRUE); tw <- sum(partI$twin, na.rm=TRUE)

  # ---- counts table ----
  counts <- rbindlist(list(
    data.table(step="0_enrolled_first_preg",     n=n_enrolled),
    data.table(step="  excluded_non_singleton",  n=n_nonsing),
    data.table(step="1_first_preg_mothers",      n=N1),
    data.table(step="  excluded_not_genotyped",  n=exGeno),
    data.table(step="2_genotyped",               n=Ngeno),
    data.table(step="  excluded_no_valid_BP",     n=exBP),
    data.table(step="    of_which_postnatal_BP_only", n=exBP_pn),
    data.table(step="    of_which_no_BP_at_all",      n=exBP_non),
    data.table(step="3_partI_analytic",          n=NpartI),
    data.table(step="  excluded_missing_PTB",     n=exPTB),
    data.table(step="4_partII_MR",               n=NpartII),
    data.table(step="descriptor_preeclampsia",   n=pe),
    data.table(step="descriptor_twins",          n=tw)))
  counts <- rbind(counts, bc[, .(step=paste0("partI_",cohort), n=n)])
  momi_write_table(counts, "F1_flow_counts", P)

  # ---- LaTeX macros (single source for figure/manuscript numbers) ----
  tex <- c(
    sprintf("\\newcommand{\\NfirstPreg}{%s}",  formatC(N1,big.mark=",",format="d")),
    sprintf("\\newcommand{\\NexGeno}{%s}",     formatC(exGeno,big.mark=",",format="d")),
    sprintf("\\newcommand{\\Ngenotyped}{%s}",  formatC(Ngeno,big.mark=",",format="d")),
    sprintf("\\newcommand{\\NexBP}{%s}",       formatC(exBP,big.mark=",",format="d")),
    sprintf("\\newcommand{\\NpartI}{%s}",      formatC(NpartI,big.mark=",",format="d")),
    sprintf("\\newcommand{\\NexPTB}{%s}",      formatC(exPTB,big.mark=",",format="d")),
    sprintf("\\newcommand{\\NpartII}{%s}",     formatC(NpartII,big.mark=",",format="d")),
    sprintf("\\newcommand{\\Npreeclampsia}{%s}", formatC(pe,big.mark=",",format="d")),
    sprintf("\\newcommand{\\Ntwins}{%s}",      formatC(tw,big.mark=",",format="d")))
  writeLines(tex, file.path(P$figures, "flow_counts.tex"))

  # ---- flow figure (ggplot) ----
  fmt <- function(x) formatC(x, big.mark=",", format="d")
  main <- data.table(
    x=3, y=c(9,6.5,4,1.5),
    lab=c(sprintf("First-pregnancy mothers\n(n = %s)", fmt(N1)),
          sprintf("Genotyped\n(GSA and/or low-pass WGS)\n(n = %s)", fmt(Ngeno)),
          sprintf("Part I analytic\n(genotyped + valid antenatal BP)\n(n = %s)", fmt(NpartI)),
          sprintf("Part II MR\n(+ preterm-birth outcome)\n(n = %s)", fmt(NpartII))))
  excl <- data.table(
    x=7.2, y=c(7.75,5.25,2.75),
    lab=c(sprintf("Excluded: not genotyped\n(n = %s)", fmt(exGeno)),
          if(is.na(exBP_pn)) sprintf("Excluded: no valid\nantenatal BP (n = %s)", fmt(exBP))
          else sprintf("Excluded: no valid antenatal BP (n = %s)\n(%s postpartum readings only, %s no BP)",
                       fmt(exBP), fmt(exBP_pn), fmt(exBP_non)),
          sprintf("Excluded: missing\nPTB outcome (n = %s)", fmt(exPTB))))
  bw<-3.4; bh<-1.15; ew<-3.0; eh<-0.95
  seg <- data.table(x=3, xend=3, y=main$y[-4]-bh/2, yend=main$y[-1]+bh/2)
  bra <- data.table(x=3, xend=excl$x-ew/2, y=excl$y, yend=excl$y)
  cohlab <- paste0("Part I by cohort:  ",
                   paste(sprintf("%s = %s", bc$display, fmt(bc$n)), collapse="   |   "),
                   sprintf("\nDescriptors kept: preeclampsia = %s, twins/multiples = %s", fmt(pe), fmt(tw)))

  g <- ggplot() +
    geom_segment(data=seg, aes(x,y,xend=xend,yend=yend),
                 arrow=arrow(length=unit(0.18,"cm"),type="closed"), linewidth=0.5) +
    geom_segment(data=bra, aes(x,y,xend=xend,yend=yend),
                 arrow=arrow(length=unit(0.15,"cm"),type="closed"), linewidth=0.4, colour="grey40") +
    geom_rect(data=main, aes(xmin=x-bw/2,xmax=x+bw/2,ymin=y-bh/2,ymax=y+bh/2),
              fill="#e8eef7", colour="#2c3e6b", linewidth=0.6) +
    geom_text(data=main, aes(x,y,label=lab), size=3.1, lineheight=0.95) +
    geom_rect(data=excl, aes(xmin=x-ew/2,xmax=x+ew/2,ymin=y-eh/2,ymax=y+eh/2),
              fill="grey96", colour="grey55", linewidth=0.4) +
    geom_text(data=excl, aes(x,y,label=lab), size=2.7, lineheight=0.95, colour="grey25") +
    annotate("text", x=0.2, y=0.15, hjust=0, vjust=1, size=2.5, colour="grey30", label=cohlab) +
    coord_cartesian(xlim=c(0,9), ylim=c(-0.4,10)) + theme_void()
  outs <- momi_save_fig(g, "F1_flow", width=8, height=6.2, P=P)

  list(n=NpartII,
       key=sprintf("firstPreg=%d -> geno=%d (-%d) -> partI=%d (-%d) -> partII=%d (-%d); PE=%d twins=%d",
                   N1, Ngeno, exGeno, NpartI, exBP, NpartII, exPTB, pe, tw),
       outputs=c(basename(outs), "F1_flow_counts.tsv","flow_counts.tex"))
})
