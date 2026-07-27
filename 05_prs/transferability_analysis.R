#!/usr/bin/env Rscript
# ============================================================
# transferability_analysis.R
# How well does each ancestry's BP PGS predict measured BP in each cohort/ancestry?
# For each (panel PGS x cohort x platform): standardize the score within that group,
# regress measured BP (the PGS's trait) on it (+ PW_AGE), and record incremental R2,
# beta, F, and N. Then summarise into a PGS-ancestry x cohort-ancestry matrix.
#
# Inputs:
#   --panel transfer_panel.tsv         (PGS_id, trait, ancestry, source)
#   --sscore-dir DIR                   (<PGS_id>__<COHORT>__<PLATFORM>.sscore)
#   --pheno-dir DIR                    (mothers_{SBP,DBP}_<win>_final.pheno, covariate_table_*)
#   --out-dir DIR  [--bp-window latest]
# Outputs: transferability_long.txt (one row per PGS x cohort x platform),
#          transferability_matrix_<trait>.txt (mean incR2: score-ancestry x cohort-ancestry)
# ============================================================
suppressMessages(library(data.table))
ga <- function(f,d=NULL){a<-commandArgs(TRUE);i<-which(a==f); if(length(i)) a[i+1] else d}
PANEL<-ga("--panel"); SSC<-ga("--sscore-dir"); PHD<-ga("--pheno-dir"); OUT<-ga("--out-dir"); WIN<-ga("--bp-window","latest")
stopifnot(!is.null(PANEL),!is.null(SSC),!is.null(PHD),!is.null(OUT)); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)

ANC <- c("AMANHI-Bangladesh"="SAS","AMANHI-Pakistan"="SAS","AMANHI-Pemba"="AFR",
         "GAPPS-Bangladesh"="SAS","GAPPS-Zambia"="AFR")
panel <- fread(PANEL)                      # PGS_id trait ancestry source
setkey(panel, PGS_id)

## phenotypes + covariate
rd <- function(fn,nm){ d<-fread(file.path(PHD,fn)); d[,.(IID, v=suppressWarnings(as.numeric(ifelse(PHENO=="NA",NA,PHENO))))][, setNames(.(IID,v),c("IID",nm))] }
ph <- merge(rd(sprintf("mothers_SBP_%s_final.pheno",WIN),"SBP"),
            rd(sprintf("mothers_DBP_%s_final.pheno",WIN),"DBP"), by="IID", all=TRUE)
cov <- fread(file.path(PHD,"covariate_table_analytic_mothers.txt"))[,.(IID=PARTICIPANT_ID, PW_AGE=suppressWarnings(as.numeric(PW_AGE)))]
ph <- merge(ph, cov, by="IID", all.x=TRUE)

zin <- function(x){ s<-sd(x,na.rm=TRUE); if(is.na(s)||s==0) rep(NA_real_,length(x)) else (x-mean(x,na.rm=TRUE))/s }

files <- list.files(SSC, pattern="__.*__.*\\.sscore$", full.names=TRUE)
res <- list()
for(f in files){
  b   <- sub("\\.sscore$","",basename(f)); p <- strsplit(b,"__",fixed=TRUE)[[1]]
  pid <- p[1]; coh <- p[2]; plat <- p[3]
  if(!(pid %in% panel$PGS_id)) next
  trait <- panel[pid, trait]; anc <- panel[pid, ancestry]; src <- panel[pid, source]
  s <- fread(f); setnames(s,1,"IID"); s <- s[,.(IID, score=SCORE1_AVG)]
  d <- merge(s, ph, by="IID")
  d <- d[!is.na(get(trait)) & !is.na(PW_AGE)]
  if(nrow(d) < 30) next
  d[, sz := zin(score)]; d <- d[!is.na(sz)]
  m  <- lm(as.formula(sprintf("%s ~ sz + PW_AGE", trait)), data=d); cc <- summary(m)$coefficients
  m0 <- lm(as.formula(sprintf("%s ~ PW_AGE", trait)), data=d)
  res[[b]] <- data.table(PGS_id=pid, score_ancestry=anc, source=src, trait=trait,
                         cohort=coh, cohort_ancestry=unname(ANC[coh]), platform=plat,
                         n=nrow(d), incR2=summary(m)$r.squared-summary(m0)$r.squared,
                         beta=cc["sz",1], se=cc["sz",2], F=(cc["sz",1]/cc["sz",2])^2, p=cc["sz",4])
}
long <- rbindlist(res)
fwrite(long, file.path(OUT,"transferability_long.txt"), sep="\t")

## ---- matrices: mean incremental R2, score-ancestry (rows) x cohort-ancestry (cols) ----
## Platform matters as much as ancestry: hard-call lpWGS loses most of the PRS signal,
## so a matrix pooled across platforms is misleading. Report PER PLATFORM, plus a pooled
## "_ALL" for reference (interpret _ALL with care — it mixes good and lossy genotypes).
mkmat <- function(d, tag){
  if(nrow(d)==0) return(invisible())
  mat <- dcast(d, score_ancestry ~ cohort_ancestry, value.var="incR2",
               fun.aggregate=function(x) mean(x,na.rm=TRUE))
  ncel <- dcast(d, score_ancestry ~ cohort_ancestry, value.var="incR2", fun.aggregate=length)
  fwrite(mat,  file.path(OUT, sprintf("transferability_matrix_%s.txt",tag)),  sep="\t")
  fwrite(ncel, file.path(OUT, sprintf("transferability_ncells_%s.txt",tag)), sep="\t")
  cat(sprintf("\n=== %s: mean incremental R2 (score ancestry x cohort ancestry) ===\n", tag))
  print(mat)
}
for(tr in sort(unique(long$trait))){
  for(pf in sort(unique(long$platform))) mkmat(long[trait==tr & platform==pf], sprintf("%s_%s",tr,pf))
  mkmat(long[trait==tr], sprintf("%s_ALL",tr))   # pooled across platforms — reference only
}
cat("\nper-cohort detail -> ", file.path(OUT,"transferability_long.txt"), "\n")
cat("Reading: compare WITHIN a platform. On GSA, expect EUR/SAS strong in SAS cohorts,\n",
    "weaker in AFR. lpWGS hard-calls lose signal (motivates the dosage re-run); the\n",
    "lpwgs_dosage matrices will appear here once dosage filesets are scored.\n")
