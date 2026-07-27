#!/usr/bin/env Rscript
# ============================================================
# plot_pca.R  (Stage 02 QC)
# Plot genetic PCs colored by cohort and by ancestry, plus a scree plot.
# Inputs: plink2 pca_results.eigenvec/.eigenval + covariate table (SITE, ANCESTRY).
# Outputs (in --out-dir): pca_by_cohort.png, pca_by_ancestry.png, pca_scree.png,
#                         pca_with_meta.txt
# Usage:
#   Rscript plot_pca.R --eigenvec FILE --eigenval FILE --covar FILE --out-dir DIR
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
getarg <- function(f,d=NULL){i<-which(args==f); if(length(i)) args[i+1] else d}
EVEC<-getarg("--eigenvec"); EVAL<-getarg("--eigenval"); COV<-getarg("--covar"); OUT<-getarg("--out-dir")
stopifnot(!is.null(EVEC),!is.null(EVAL),!is.null(COV),!is.null(OUT))
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

ev <- fread(EVEC)
setnames(ev, 1:2, c("FID","IID"))
pcs <- grep("^PC", names(ev), value=TRUE)
if(!length(pcs)){ # plink2 may name them V3.. ; rename trailing cols
  setnames(ev, 3:ncol(ev), paste0("PC", seq_len(ncol(ev)-2)))
  pcs <- grep("^PC", names(ev), value=TRUE)
}
cov <- fread(COV)
m <- merge(ev, cov[,.(PARTICIPANT_ID, SITE, ANCESTRY)], by.x="IID", by.y="PARTICIPANT_ID", all.x=TRUE)
m[is.na(SITE), SITE:="UNKNOWN"]; m[is.na(ANCESTRY), ANCESTRY:="UNKNOWN"]
fwrite(m, file.path(OUT,"pca_with_meta.txt"), sep="\t")

gg <- function(colvar, fn, title){
  p <- ggplot(m, aes(PC1, PC2, color=.data[[colvar]])) +
    geom_point(size=1, alpha=0.7) + theme_bw() +
    labs(title=title, color=colvar)
  ggsave(file.path(OUT,fn), p, width=7, height=5, dpi=150)
}
gg("SITE","pca_by_cohort.png","Genetic PCs by cohort")
gg("ANCESTRY","pca_by_ancestry.png","Genetic PCs by ancestry")

eval <- fread(EVAL, header=FALSE)$V1
scree <- data.table(PC=seq_along(eval), eigenvalue=eval, pct=100*eval/sum(eval))
p <- ggplot(scree, aes(PC, pct)) + geom_col(fill="steelblue") + geom_line() + geom_point() +
  theme_bw() + labs(title="Scree plot", y="% variance", x="PC") + scale_x_continuous(breaks=scree$PC)
ggsave(file.path(OUT,"pca_scree.png"), p, width=7, height=4, dpi=150)

cat("PCA plots written to", OUT, "\n")
cat("variance % by PC:", sprintf("%.1f", scree$pct), "\n")
