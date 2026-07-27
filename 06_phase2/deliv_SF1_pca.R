#!/usr/bin/env Rscript
# ============================================================
# deliv_SF1_pca.R  [ B28 -> Supp Figure SF1, ancestry of the five cohorts ]
#
# ------------------------------------------------------------------
# WHAT THIS FIGURE IS FOR. momi_config.R asserts an ancestry label for every cohort --
# three SAS, two AFR -- and those labels do real work: they define the strata in T4/F4, they
# justify reporting the African MR arm separately, and they underpin the claim that the
# instrument fails to transfer to African cohorts. A reviewer is entitled to ask whether the
# labels are supported by the genotypes. This figure answers that, and it is the ONLY place
# in the paper where the labels are checked rather than assumed.
#
# THE SHRINKAGE ARTEFACT, WHICH MUST BE STATED RATHER THAN QUIETLY IGNORED.
# Projected PCs are regression predictions, so they have LESS VARIANCE than the reference
# coordinates they are projected onto. Every projected sample is therefore pulled toward the
# origin relative to the 1000 Genomes cloud. The consequence for reading this figure is
# specific and serious: A COHORT SITTING "BETWEEN" TWO REFERENCE CLUSTERS MAY BE AN ARTEFACT
# OF SHRINKAGE RATHER THAN EVIDENCE OF ADMIXTURE. B06b flagged this and SF1 is where it has
# to be handled.
#
# We handle it two ways. First, the reference samples are plotted with the SAME projection
# applied to a held-out portion is not available, so instead the caption states the artefact
# explicitly and the reference cloud is drawn at low alpha as context, not as a measuring
# stick. Second, and more usefully, the QUANTITATIVE claim rests on Fst (S3b), which has no
# such artefact -- the figure is descriptive and Fst is the evidence.
#
# THE SPECIFIC QUESTION THIS WAS BUILT TO TEST. AMANHI-Pemba is labelled AFR. Zanzibar has
# documented Omani and South Asian admixture from centuries of Indian Ocean trade, so the
# label may be doing work the data does not support. If Pemba sits between the AFR and SAS
# reference clusters BY MORE THAN SHRINKAGE EXPLAINS, the label is wrong and the AFR stratum
# in T4 is not a coherent group. Note S3b already gives the cleaner answer: Pemba vs Zambia
# Fst = 0.00434, which is 25x the divergence between the two Bangladeshi cohorts (0.00017) --
# i.e. the two "AFR" cohorts are NOT interchangeable, whatever the PCA shows.
# ------------------------------------------------------------------
#   Writes: figures/SF1_pca.{png,pdf}, tables/SF1_pca_centroids.tsv
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

momi_deliverable("SF1_pca", script="06_phase2/deliv_SF1_pca.R",
                 inputs="pcs;pcs_reference", P=P, stop_on_error=FALSE, body=function(ctx){

  PC  <- tryCatch(momi_read_intermediate("pcs", P), error=function(e) NULL)
  REF <- tryCatch(momi_read_intermediate("pcs_reference", P), error=function(e) NULL)
  if(is.null(PC) || !nrow(PC)) return(list(skip=TRUE, reason="no pcs.rds -- run B06b"))

  pj <- c("PC1_1kg","PC2_1kg")
  if(!all(pj %in% names(PC)))
    return(list(skip=TRUE, reason="pcs.rds has no projected PCs; SF1 needs the 1000G projection"))

  D <- PC[is.finite(get(pj[1])) & is.finite(get(pj[2]))]
  D[, Cohort := MOMI_DISPLAY[cohort]]
  D[, label  := sprintf("%s [%s]", Cohort, MOMI_ANC[cohort])]

  ## ---- centroids: the quantitative content of the figure ----
  cen <- D[, .(n=.N,
               PC1=mean(get(pj[1])), PC2=mean(get(pj[2])),
               PC1_sd=sd(get(pj[1])), PC2_sd=sd(get(pj[2]))),
           by=.(cohort, Cohort, ancestry=MOMI_ANC[cohort])][order(PC1)]
  out_t <- momi_write_table(cen, "SF1_pca_centroids", P)

  ## ---- per-individual points for a standard PCA scatter (report/build_figures.py) ----
  ## Sample-wide in-sample PCs preferred (per Anagh: no need to project to 1000G); fall
  ## back to the projected PCs if plain PC1/PC2 are absent. Subsample for a light figure.
  ppc <- if(all(c("PC1","PC2") %in% names(PC))) c("PC1","PC2") else pj
  PTS <- PC[is.finite(get(ppc[1])) & is.finite(get(ppc[2])),
            .(cohort, PC1=get(ppc[1]), PC2=get(ppc[2]))]
  set.seed(1)
  PTS <- PTS[, .SD[sample(.N, min(.N, 1500))], by=cohort]   # <=1500 points/cohort
  momi_write_table(PTS, "SF1_pca_points", P)

  ## ---- reference backdrop, if the panel file gave superpopulations ----
  refok <- !is.null(REF) && all(c("PC1","PC2") %in% names(REF))
  rp <- NULL
  if(refok){
    rp <- copy(REF)
    if(!("superpop" %in% names(rp))) rp[, superpop := "1000G"]
    rp <- rp[is.finite(PC1) & is.finite(PC2)]
  }

  g <- ggplot()
  if(!is.null(rp) && nrow(rp))
    g <- g + geom_point(data=rp, aes(PC1, PC2, colour=superpop),
                        size=.5, alpha=.18, show.legend=TRUE)
  g <- g +
    geom_point(data=D, aes(get(pj[1]), get(pj[2]), shape=label),
               size=.7, alpha=.45, colour="grey20") +
    geom_point(data=cen, aes(PC1, PC2), size=3.4, shape=4, stroke=1.1, colour="black") +
    labs(x="PC1 (projected onto 1000 Genomes)", y="PC2 (projected)",
         shape=NULL, colour="1000G superpop",
         title="Cohort ancestry: MOMI samples projected onto 1000 Genomes axes",
         subtitle=paste0("Crosses are cohort centroids. NOTE: projected PCs shrink toward the origin relative to the reference ",
                         "(a projection has less variance than what it predicts), so a cohort lying between reference clusters ",
                         "is NOT by itself evidence of admixture. Quantitative divergence is in S3b (Fst).")) +
    theme_minimal(base_size=8) +
    theme(panel.grid.minor=element_blank(),
          legend.position="right",
          plot.subtitle=element_text(size=6, colour="grey35"))
  fig <- momi_save_fig(g, "SF1_pca", width=9, height=6, P=P)

  ## ---- console ----
  cat("\n=== projected PC centroids per cohort ===\n")
  print(cen)
  cat("\nThis is what tests whether the config's ancestry labels survive contact with the\n",
      "genotypes. AMANHI-Pemba is the case to watch: Zanzibar has documented Omani and South\n",
      "Asian admixture, so if it sits appreciably toward the SAS cohorts the AFR stratum in\n",
      "T4 is not a coherent group. Read this ALONGSIDE S3b -- Pemba vs Zambia Fst = 0.00434,\n",
      "25x the two Bangladeshi cohorts' 0.00017, so the two AFR cohorts already differ more\n",
      "from each other than our two SAS cohorts do.\n", sep="")

  if(!is.null(rp) && "superpop" %in% names(rp)){
    cat("\n=== 1000G reference centroids, for comparison ===\n")
    print(rp[, .(n=.N, PC1=round(mean(PC1),4), PC2=round(mean(PC2),4)), by=superpop][order(PC1)])
    cat("\nCompare each cohort's centroid to these. Remember the shrinkage: cohort centroids are\n",
        "systematically closer to the origin than reference centroids on the same axes, so\n",
        "compare DIRECTION and RANK ORDER rather than absolute distance.\n", sep="")
  }

  list(n=nrow(D),
       key=sprintf("mothers plotted=%d cohorts=%d; reference samples=%s",
                   nrow(D), uniqueN(D$cohort),
                   if(is.null(rp)) "none" else as.character(nrow(rp))),
       outputs=c(basename(out_t), basename(fig)))
})
