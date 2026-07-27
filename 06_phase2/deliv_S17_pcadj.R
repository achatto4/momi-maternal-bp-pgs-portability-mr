#!/usr/bin/env Rscript
# ============================================================
# deliv_S17_pcadj.R  [ B30 -> Supp Table S17, is Part I an artefact of population structure? ]
#
# ------------------------------------------------------------------
# THE OBJECTION THIS ANSWERS
#
# Part I reports incremental R2 of a polygenic score on measured blood pressure, adjusting
# only for age. A reviewer will ask whether those R2 values are inflated by POPULATION
# STRUCTURE: if a cohort contains subgroups that differ in both allele frequencies and mean
# blood pressure -- through diet, altitude, ethnicity, migration history, anything -- then a
# score built from allele frequencies can predict blood pressure WITHOUT any of the causal
# variants being relevant in that population. The score would be tagging structure, not
# biology, and the "transferability" we report would be partly spurious.
#
# This is a serious objection and it bites hardest exactly where our headline sits: the
# cohort with the highest R2 (GAPPS-Bangladesh, 6.46%) is also the largest and most
# geographically concentrated, and Matlab has well-documented internal population structure.
#
# THE TEST. Re-fit every transferability model with the within-cohort principal components
# from B06b added to the covariates, and compare. Within-cohort PCs are precisely the right
# control: they describe structure INSIDE the analysis sample, which is the confounding at
# issue. (Projected 1000G PCs would not do -- they describe where a cohort sits globally,
# which is a different question and is SF1's job.)
#
# HOW TO READ THE RESULT:
#   * R2 essentially unchanged  -> the score is not tagging structure; Part I stands as
#     reported, and this table is the evidence for that sentence in the Methods.
#   * R2 falls substantially    -> some of the reported transferability WAS structure. The
#     PC-adjusted values become the ones to report, and the portability gap must be
#     re-examined under them.
#   * R2 falls in some cohorts but not others -> the most awkward outcome, because it would
#     mean the CROSS-COHORT CONTRAST that Part I rests on is partly a contrast in how much
#     structure each cohort contains. Watch GAPPS-Bangladesh vs AMANHI-Bangladesh specifically:
#     if PC adjustment narrows that particular gap, F3's timing explanation loses ground and
#     structure gains it.
#
# THIS ALSO VALIDATES THE MR COVARIATE SET. B22 adjusts the MR for age plus these same
# within-cohort PCs, following Burgess (ref/methods_citations.tsv). If the PCs turn out to
# matter here, that choice is doing real work; if they change nothing, the MR estimates are
# insensitive to it and that is worth stating rather than leaving implicit.
# ------------------------------------------------------------------
#   Writes: tables/S17_pc_adjusted.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

NPC_ADJ <- as.integer(momi_arg("--npc", "5"))

momi_deliverable("S17_pcadj", script="06_phase2/deliv_S17_pcadj.R",
                 inputs="analytic_mothers;prs_z;pcs", P=P, stop_on_error=FALSE,
                 body=function(ctx){

  A  <- momi_read_intermediate("analytic_mothers", P)
  Z  <- momi_read_intermediate("prs_z", P)
  PC <- tryCatch(momi_read_intermediate("pcs", P), error=function(e) NULL)
  if(is.null(PC) || !nrow(PC)) return(list(skip=TRUE, reason="no pcs.rds -- run B06b first"))
  pccols <- intersect(paste0("PC", seq_len(NPC_ADJ)), names(PC))
  if(!length(pccols)) return(list(skip=TRUE, reason="pcs.rds has no within-cohort PCs"))

  A2 <- merge(A, PC[, c("IID","cohort",pccols), with=FALSE], by=c("IID","cohort"), all.x=TRUE)
  BPDEF <- MOMI_BPDEF_DEFAULT

  rows <- list()
  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    sid   <- momi_instrument(coh, tr)
    bpcol <- paste0(substr(tr,1,1), "_", BPDEF)
    if(!(bpcol %in% names(A2))) next
    zc <- Z[score_id==sid & cohort==coh, .(IID, z)]
    if(!nrow(zc)) next
    d <- merge(A2[cohort==coh], zc, by="IID")
    ## SAME SAMPLE for both models, or the comparison is confounded by who has PCs.
    ## Without this, a cohort losing mothers to missing PCs would show an R2 change that is
    ## really a change of sample -- the single easiest way to get this table wrong.
    need <- c("z", bpcol, "AGE", pccols)
    d <- d[complete.cases(d[, ..need])]
    if(nrow(d) < 100) next

    base <- momi_transfer(d$z, d[[bpcol]], cov=d[, .(AGE)])
    padj <- momi_transfer(d$z, d[[bpcol]], cov=d[, c("AGE", pccols), with=FALSE])
    if(is.na(base$incR2) || is.na(padj$incR2)) next

    rows[[paste(coh,tr)]] <- data.table(
      cohort_display=MOMI_DISPLAY[[coh]], cohort=coh, ancestry=MOMI_ANC[[coh]],
      trait=tr, chosen_PGS=sid, N=base$n, n_pcs=length(pccols),
      R2pct_age_only = round(100*base$incR2, 3),
      R2pct_age_plus_PCs = round(100*padj$incR2, 3),
      delta_pp = round(100*(padj$incR2 - base$incR2), 3),
      pct_change = round(100*(padj$incR2 - base$incR2)/base$incR2, 1),
      F_age_only = round(base$F, 1), F_age_plus_PCs = round(padj$F, 1),
      beta_age_only = round(base$beta, 3), beta_age_plus_PCs = round(padj$beta, 3))
  }
  if(!length(rows)) return(list(skip=TRUE, reason="no cohort x trait cell estimable"))
  S17 <- rbindlist(rows)
  setorder(S17, ancestry, cohort_display, trait)
  out <- momi_write_table(S17, "S17_pc_adjusted", P)

  ## ---- console ----
  cat(sprintf("\n=== transferability with and without %d within-cohort PCs (same sample in both) ===\n",
              length(pccols)))
  print(S17[, .(cohort_display, trait, N, R2pct_age_only, R2pct_age_plus_PCs,
                delta_pp, pct_change)])

  cat(sprintf("\nlargest absolute change: %.3f pp (%s %s); median %%change: %.1f%%\n",
              max(abs(S17$delta_pp)),
              S17[which.max(abs(delta_pp)), cohort_display],
              S17[which.max(abs(delta_pp)), trait],
              median(S17$pct_change)))

  ## ---- the specific comparison Part I depends on ----
  cat("\n=== the contrast Part I rests on: the two Bangladeshi cohorts ===\n")
  bd <- S17[cohort %in% c("AMANHI-Bangladesh","GAPPS-Bangladesh")]
  if(nrow(bd)){
    print(bd[, .(cohort_display, trait, R2pct_age_only, R2pct_age_plus_PCs, pct_change)])
    for(tr in c("SBP","DBP")){
      a <- bd[cohort=="AMANHI-Bangladesh" & trait==tr]
      g <- bd[cohort=="GAPPS-Bangladesh"  & trait==tr]
      if(nrow(a) && nrow(g))
        cat(sprintf("  %s gap: %.2f-fold unadjusted -> %.2f-fold PC-adjusted\n", tr,
                    g$R2pct_age_only/a$R2pct_age_only,
                    g$R2pct_age_plus_PCs/a$R2pct_age_plus_PCs))
    }
    cat("\n  If the fold-gap barely moves, population structure is NOT the explanation for the\n",
        "  portability difference, and F3's measurement-timing account survives this challenge.\n",
        "  If it narrows materially, structure is doing part of the work and Part I must say so.\n", sep="")
  }

  list(n=nrow(S17),
       key=sprintf("cells=%d npc=%d; max |dR2|=%.3f pp; median change=%.1f%%",
                   nrow(S17), length(pccols), max(abs(S17$delta_pp)), median(S17$pct_change)),
       outputs=basename(out))
})
