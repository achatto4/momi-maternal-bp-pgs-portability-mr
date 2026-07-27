#!/usr/bin/env Rscript
# ============================================================
# make_covar.R  (Stage 03 helper)
# Build a PLINK2 covariate file: FID IID PW_AGE PC1..PCn + SITE dummies.
# (Fixes the old-pipeline gap where PCs were computed but never used in GWAS.)
# Usage:
#   Rscript make_covar.R --covar COVTABLE --eigenvec EVEC --out FILE [--n-pcs 10]
# ============================================================
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
ga <- function(f,d=NULL){i<-which(args==f); if(length(i)) args[i+1] else d}
COV<-ga("--covar"); EVEC<-ga("--eigenvec"); OUT<-ga("--out"); NPC<-as.integer(ga("--n-pcs","10"))
stopifnot(!is.null(COV),!is.null(EVEC),!is.null(OUT))

cov <- fread(COV)
ev  <- fread(EVEC); setnames(ev,1:2,c("FID","IID"))
pcn <- grep("^PC", names(ev), value=TRUE)
if(!length(pcn)){ setnames(ev,3:ncol(ev),paste0("PC",seq_len(ncol(ev)-2))); pcn<-grep("^PC",names(ev),value=TRUE) }
pcn <- pcn[seq_len(min(NPC,length(pcn)))]

m <- merge(ev[, c("FID","IID",pcn), with=FALSE],
           cov[, .(IID=PARTICIPANT_ID, PW_AGE=suppressWarnings(as.numeric(PW_AGE)), SITE)],
           by="IID")
# one-hot SITE dummies, dropping the first (reference) level
# NOTE: no SITE dummies — site is collinear with the top PCs (PC1 = SAS/AFR ancestry),
# which causes plink2 VIF_TOO_HIGH. The genetic PCs already control site/ancestry structure.
m[, SITE := NULL]
setcolorder(m, c("FID","IID"))
fwrite(m, OUT, sep="\t", na="NA", quote=FALSE)   # plink needs UNquoted FID/IID
cat(sprintf("covar written: %s  (%d samples, cols: %s)\n", OUT, nrow(m),
            paste(setdiff(names(m),c("FID","IID")), collapse=", ")))
