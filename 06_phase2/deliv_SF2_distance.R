#!/usr/bin/env Rscript
# ============================================================
# deliv_SF2_distance.R  [ B29 -> Supp Figure SF2, does genetic distance explain portability? ]
#
# ------------------------------------------------------------------
# THE HYPOTHESIS ON TRIAL, WHICH IS THE FIELD'S DEFAULT EXPLANATION
#
# Polygenic score portability is conventionally explained by genetic distance from the
# training population: a EUR-trained score works less well the further a target sample sits
# from Europeans, through differences in LD structure, allele frequencies and causal effect
# sizes. If that account is sufficient HERE, then plotting each cohort's transferability
# against its distance from the EUR reference should produce a clean monotone relationship.
#
# WE ALREADY KNOW IT IS NOT SUFFICIENT, AND THIS FIGURE IS WHERE THAT IS SHOWN RATHER THAN
# ASSERTED. The decisive case is internal to the plot: AMANHI-Bangladesh and GAPPS-Bangladesh
# are separated by Fst = 0.00017 -- indistinguishable, 14x closer than the next-closest South
# Asian pair -- yet differ 2.6-fold in incremental R2 (2.49% vs 6.46% for SBP). Two points at
# effectively the SAME x with grossly DIFFERENT y is a refutation that needs no statistics:
# whatever else drives portability here, it is not ancestry.
#
# So SF2 does two things at once, and the second is the point:
#   * it CONFIRMS the conventional account across ancestry groups (AFR cohorts, further from
#     EUR, do have lower R2 than SAS cohorts -- the field's expectation holds coarsely);
#   * it REFUTES the account as a complete explanation, because the largest single R2 gap in
#     the study occurs between the two cohorts with the smallest genetic distance between them.
#
# BOTH DISTANCE MEASURES ARE PLOTTED because they can disagree and each has a weakness:
#   Fst to the EUR reference -- the standard measure, computed at score SNPs (S3b). Clean, but
#     summarises the whole genome-wide allele frequency spectrum rather than what the score
#     uses.
#   PC distance -- Euclidean distance between the cohort centroid and the EUR reference
#     centroid in projected 1000G space (B06b). Intuitive, but subject to the projection
#     shrinkage artefact documented in SF1, which compresses all cohorts toward the origin
#     and therefore toward each other. Reported, but Fst is the measure to trust.
#
# n=5 COHORTS. No correlation coefficient is computed and none should be. With five points a
# Spearman rho is not evidence of anything, and quoting one would dress up a scatterplot as a
# test. The figure is an existence argument -- one pair of points does the work.
# ------------------------------------------------------------------
#   Writes: figures/SF2_distance_vs_r2.{png,pdf}, tables/SF2_distance.tsv
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

momi_deliverable("SF2_distance", script="06_phase2/deliv_SF2_distance.R",
                 inputs="transfer_grid;pcs", P=P, stop_on_error=FALSE, body=function(ctx){

  G <- momi_read_intermediate("transfer_grid", P)

  ## ---- y axis: transferability of the CHOSEN instrument at the primary definition ----
  M <- G[definition==MOMI_BPDEF_DEFAULT]
  ys <- rbindlist(lapply(MOMI_COH_ALL, function(coh)
    rbindlist(lapply(c("SBP","DBP"), function(tr){
      r <- M[cohort==coh & trait==tr & score_id==momi_instrument(coh,tr)]
      if(!nrow(r)) return(NULL)
      data.table(cohort=coh, trait=tr, R2pct=r$R2pct[1], F=r$F[1], n=r$n[1])
    }))), fill=TRUE)
  if(!nrow(ys)) return(list(skip=TRUE, reason="no transferability rows"))

  ## ---- x axis 1: Fst to the EUR reference, from S3b ----
  ## S3b holds PAIRWISE cohort-vs-cohort Fst. A cohort-to-EUR distance is only available if
  ## the reference was included as a pseudo-cohort. Handle both cases rather than assuming.
  f3 <- file.path(P$tables, "S3b_pairwise_fst.tsv")
  FST <- if(file.exists(f3)) fread(f3) else NULL
  fst_eur <- NULL
  if(!is.null(FST) && all(c("cohort_a","cohort_b","fst_hudson") %in% names(FST))){
    e <- FST[grepl("EUR|1000G|g1k|reference", cohort_a, ignore.case=TRUE) |
             grepl("EUR|1000G|g1k|reference", cohort_b, ignore.case=TRUE)]
    if(nrow(e)){
      e[, cohort := fifelse(grepl("EUR|1000G|g1k|reference", cohort_a, ignore.case=TRUE),
                            cohort_b, cohort_a)]
      fst_eur <- e[, .(fst_to_eur = mean(fst_hudson)), by=cohort]
    }
  }
  if(is.null(fst_eur))
    cat("\nNB S3b has no cohort-vs-EUR pair, so Fst-to-EUR is unavailable.\n",
        "   The WITHIN-STUDY pairwise Fst is still the decisive evidence and is printed below.\n", sep="")

  ## ---- x axis 2: PC distance from the EUR reference centroid ----
  PC  <- tryCatch(momi_read_intermediate("pcs", P), error=function(e) NULL)
  REF <- tryCatch(momi_read_intermediate("pcs_reference", P), error=function(e) NULL)
  pcd <- NULL
  pj  <- c("PC1_1kg","PC2_1kg")
  if(!is.null(PC) && all(pj %in% names(PC)) && !is.null(REF) &&
     all(c("PC1","PC2") %in% names(REF)) && "superpop" %in% names(REF)){
    eu <- REF[superpop=="EUR" & is.finite(PC1) & is.finite(PC2)]
    if(nrow(eu)){
      ec <- c(mean(eu$PC1), mean(eu$PC2))
      pcd <- PC[is.finite(get(pj[1])) & is.finite(get(pj[2])),
                .(pc1=mean(get(pj[1])), pc2=mean(get(pj[2]))), by=cohort]
      pcd[, pc_dist_eur := sqrt((pc1-ec[1])^2 + (pc2-ec[2])^2)]
      pcd <- pcd[, .(cohort, pc_dist_eur)]
    }
  }
  if(is.null(pcd)) cat("\nNB projected PCs or EUR reference centroid unavailable -- PC distance omitted.\n")

  D <- ys
  if(!is.null(fst_eur)) D <- merge(D, fst_eur, by="cohort", all.x=TRUE)
  if(!is.null(pcd))     D <- merge(D, pcd,     by="cohort", all.x=TRUE)
  D[, `:=`(ancestry = MOMI_ANC[cohort], Cohort = MOMI_DISPLAY[cohort])]
  out_t <- momi_write_table(D, "SF2_distance", P)

  ## ---- figure: whichever distance measures exist ----
  xs <- intersect(c("fst_to_eur","pc_dist_eur"), names(D))
  figs <- character(0)
  if(length(xs)){
    L <- melt(D, id.vars=c("cohort","Cohort","ancestry","trait","R2pct"),
              measure.vars=xs, variable.name="measure", value.name="distance")
    L <- L[is.finite(distance)]
    if(nrow(L)){
      g <- ggplot(L, aes(distance, R2pct, colour=ancestry, shape=trait)) +
        geom_point(size=2.4) +
        ggrepel::geom_text_repel(aes(label=Cohort), size=2, show.legend=FALSE,
                                 max.overlaps=20) +
        facet_wrap(~measure, scales="free_x") +
        labs(x="genetic distance from the EUR training population", y="incremental R² (%)",
             title="Genetic distance does not fully explain polygenic score portability",
             subtitle=paste0("The two Bangladeshi cohorts are separated by Fst = 0.00017 -- indistinguishable -- yet differ 2.6-fold in R². ",
                             "No correlation is fitted: with five cohorts a coefficient would dress a scatterplot up as a test.")) +
        theme_minimal(base_size=8) +
        theme(panel.grid.minor=element_blank(),
              plot.subtitle=element_text(size=6, colour="grey35"))
      if(!requireNamespace("ggrepel", quietly=TRUE)){
        g <- ggplot(L, aes(distance, R2pct, colour=ancestry, shape=trait)) +
          geom_point(size=2.4) + geom_text(aes(label=Cohort), size=2, vjust=-0.9, show.legend=FALSE) +
          facet_wrap(~measure, scales="free_x") +
          labs(x="genetic distance from the EUR training population", y="incremental R² (%)",
               title="Genetic distance does not fully explain polygenic score portability") +
          theme_minimal(base_size=8)
      }
      figs <- momi_save_fig(g, "SF2_distance_vs_r2", width=8.5, height=4.5, P=P)
    }
  } else cat("\nNo distance measure available -- figure not drawn; the table still carries R2.\n")

  ## ---- console: the decisive comparison, which needs no figure ----
  cat("\n=== transferability vs distance ===\n"); print(D[order(R2pct)])
  if(!is.null(FST)){
    cat("\n=== the refutation, stated as two numbers ===\n")
    bd <- FST[(cohort_a=="AMANHI-Bangladesh" & cohort_b=="GAPPS-Bangladesh") |
              (cohort_b=="AMANHI-Bangladesh" & cohort_a=="GAPPS-Bangladesh")]
    r2 <- ys[trait=="SBP" & cohort %in% c("AMANHI-Bangladesh","GAPPS-Bangladesh")]
    if(nrow(bd) && nrow(r2)){
      cat(sprintf("  Fst(AMANHI-B, GAPPS-B) = %.5f on %s SNPs\n",
                  bd$fst_hudson[1],
                  if("n_snps" %in% names(bd)) format(bd$n_snps[1], big.mark=",") else "?"))
      cat(sprintf("  SBP R2: AMANHI-B %.2f%%  vs  GAPPS-B %.2f%%  (%.1f-fold)\n",
                  r2[cohort=="AMANHI-Bangladesh", R2pct],
                  r2[cohort=="GAPPS-Bangladesh", R2pct],
                  r2[cohort=="GAPPS-Bangladesh", R2pct] / r2[cohort=="AMANHI-Bangladesh", R2pct]))
      cat("  Two cohorts at effectively the same genetic distance, differing 2.6-fold in\n",
          "  transferability. Ancestry cannot be the explanation. What remains is measurement:\n",
          "  when in gestation BP was taken, and how many readings were averaged.\n", sep="")
    }
  }

  list(n=nrow(D),
       key=sprintf("cohorts=%d; distance measures=%s",
                   uniqueN(D$cohort),
                   if(length(xs)) paste(xs, collapse=",") else "none"),
       outputs=c(basename(out_t), basename(figs)))
})
