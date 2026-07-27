#!/usr/bin/env Rscript
# ============================================================
# pooled_prs_analysis.R  (combined PRS-based MR across ALL cohorts x platforms)
# Pools per-sample PGS computed separately on each cohort x platform bed
# (score_all_cohorts.sh), standardizes WITHIN cohort x platform, then runs:
#   (1) instrument strength : BP  ~ PRS_z + PW_AGE + factor(GROUP)
#   (2) PRS -> PTB          : PTB ~ PRS_z + PW_AGE + factor(GROUP)  (binomial)
#   (3) PRS-based MR (Wald) : BP->PTB = beta_ptb/beta_bp
# GROUP = cohort__platform fixed effect controls site/platform/ancestry structure
# (no cross-cohort PCs, since genotypes aren't pooled; within-group standardization
#  + GROUP fixed effects absorb between-group differences).
#
# Usage:
#   Rscript pooled_prs_analysis.R --prs-dir DIR --pheno-dir DIR --out-dir DIR
#           [--bp-window latest] [--ptb-case-code 2]
# ============================================================
suppressMessages(library(data.table))
have_proc <- requireNamespace("pROC", quietly=TRUE)
args <- commandArgs(trailingOnly=TRUE)
ga <- function(f,d=NULL){i<-which(args==f); if(length(i)) args[i+1] else d}
PRSDIR<-ga("--prs-dir"); PHDIR<-ga("--pheno-dir"); OUT<-ga("--out-dir")
WIN<-ga("--bp-window","latest"); PTBCASE<-as.integer(ga("--ptb-case-code","2"))
stopifnot(!is.null(PRSDIR),!is.null(PHDIR),!is.null(OUT)); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)

ANC <- c("AMANHI-Bangladesh"="SAS","AMANHI-Pakistan"="SAS","AMANHI-Pemba"="AFR",
         "GAPPS-Bangladesh"="SAS","GAPPS-Zambia"="AFR")
TAGS <- c("SBP_EUR","DBP_EUR","SBP_SAS","DBP_SAS")

## ---- read all sscores: <COHORT>__<PLATFORM>__<TAG>.sscore ----
files <- list.files(PRSDIR, pattern="__.*\\.sscore$", full.names=TRUE)
stopifnot(length(files) > 0)
meta <- function(f){ b<-sub("\\.sscore$","",basename(f)); p<-strsplit(b,"__",fixed=TRUE)[[1]]
                     list(cohort=p[1], platform=p[2], tag=p[3]) }
groups <- unique(sapply(files, function(f){m<-meta(f); paste(m$cohort,m$platform,sep="__")}))

read_group <- function(grp){
  cp <- strsplit(grp,"__",fixed=TRUE)[[1]]; coh<-cp[1]; plat<-cp[2]
  dt <- NULL
  for(tg in TAGS){
    f <- file.path(PRSDIR, paste0(grp,"__",tg,".sscore")); if(!file.exists(f)) next
    d <- fread(f); setnames(d,1,"IID"); d <- d[,.(IID, score=SCORE1_AVG)]; setnames(d,"score",tg)
    dt <- if(is.null(dt)) d else merge(dt,d,by="IID",all=TRUE)
  }
  if(is.null(dt)) return(NULL)
  dt[, `:=`(COHORT=coh, PLATFORM=plat, GROUP=grp, ANCESTRY=ANC[[coh]])]; dt
}
prs <- rbindlist(lapply(groups, read_group), fill=TRUE)
cat(sprintf("pooled samples scored: %d across %d groups\n", nrow(prs), length(groups)))

## ---- phenotypes + covariates (built on all-cohort mothers) ----
read_ph <- function(fn,name){ d<-fread(file.path(PHDIR,fn))
  d[,.(IID, v=suppressWarnings(as.numeric(ifelse(PHENO=="NA",NA,PHENO))))][, setNames(.(IID,v),c("IID",name))] }
sbp <- read_ph(sprintf("mothers_SBP_%s_final.pheno",WIN),"SBP")
dbp <- read_ph(sprintf("mothers_DBP_%s_final.pheno",WIN),"DBP")
ptbf<- fread(file.path(PHDIR,"mothers_PTB_NEW_final.pheno"))
ptb <- ptbf[,.(IID, PTB=ifelse(PHENO=="NA",NA, as.integer(as.integer(PHENO)==PTBCASE)))]
cov <- fread(file.path(PHDIR,"covariate_table_analytic_mothers.txt"))[,.(IID=PARTICIPANT_ID, PW_AGE=suppressWarnings(as.numeric(PW_AGE)))]

d <- Reduce(function(a,b) merge(a,b,by="IID",all.x=TRUE), list(prs,sbp,dbp,ptb,cov))
d <- d[IID %in% cov$IID]                 # keep analytic mothers only (those in covariate table)
fwrite(d, file.path(OUT,"pooled_prs_dataset.txt"), sep="\t", na="NA")
cat(sprintf("analytic pooled mothers: %d\n", nrow(d)))

## ---- standardize each score WITHIN group ----
zin <- function(x){ s<-sd(x,na.rm=TRUE); if(is.na(s)||s==0) return(rep(NA_real_,length(x))); (x-mean(x,na.rm=TRUE))/s }
for(tg in intersect(TAGS,names(d))) d[, (paste0(tg,"_z")) := zin(get(tg)), by=GROUP]

analyze <- function(dd, trait, anc, scope){
  tg <- paste0(trait,"_",anc,"_z"); if(!(tg %in% names(dd))) return(NULL)
  dd <- dd[!is.na(get(tg)) & !is.na(get(trait)) & !is.na(PTB) & !is.na(PW_AGE)]
  if(nrow(dd) < 50) return(NULL)
  ncase<-sum(dd$PTB==1); nctrl<-sum(dd$PTB==0); if(ncase<10||nctrl<10) return(NULL)
  ng <- length(unique(dd$GROUP))
  base <- if(scope=="combined" && ng>1) "PW_AGE + factor(GROUP)" else "PW_AGE"
  m_bp <- lm(as.formula(sprintf("%s ~ %s + %s", trait, tg, base)), data=dd); cb<-summary(m_bp)$coefficients
  b_bp<-cb[tg,1]; se_bp<-cb[tg,2]
  m0  <- lm(as.formula(sprintf("%s ~ %s", trait, base)), data=dd)
  incR2 <- summary(m_bp)$r.squared - summary(m0)$r.squared; Fstat <- (b_bp/se_bp)^2
  m_ptb <- glm(as.formula(sprintf("PTB ~ %s + %s", tg, base)), data=dd, family=binomial()); cp<-summary(m_ptb)$coefficients
  b_ptb<-cp[tg,1]; se_ptb<-cp[tg,2]; p_ptb<-cp[tg,4]
  auc <- if(have_proc) tryCatch(as.numeric(pROC::auc(dd$PTB, predict(m_ptb,type="response"),quiet=TRUE)),error=function(e)NA_real_) else NA_real_
  ratio<-b_ptb/b_bp; se_ratio<-sqrt(se_ptb^2/b_bp^2 + b_ptb^2*se_bp^2/b_bp^4); p_mr<-2*pnorm(-abs(ratio/se_ratio))
  data.table(scope=scope,trait=trait,ancestry=anc,n=nrow(dd),n_case=ncase,n_ctrl=nctrl,
    instr_incR2=incR2,instr_F=Fstat,prs_ptb_OR=exp(b_ptb),prs_ptb_p=p_ptb,prs_ptb_AUC=auc,
    mr_OR=exp(ratio),mr_se=se_ratio,mr_p=p_mr)
}

res<-list()
for(trait in c("SBP","DBP")) for(anc in c("EUR","SAS")){
  res[[paste("comb",trait,anc)]] <- analyze(d,trait,anc,"combined")
  for(g in sort(unique(d$GROUP))){ r<-analyze(d[GROUP==g],trait,anc,g); if(!is.null(r)) res[[paste(g,trait,anc)]]<-r }
}
out<-rbindlist(Filter(Negate(is.null),res))
fwrite(out, file.path(OUT,"pooled_prs_mr_results.txt"), sep="\t")
cat("\n=== POOLED (all cohorts x platforms) PRS-based MR ===\n")
print(out[scope=="combined", .(trait,ancestry,n,n_case,incR2=round(instr_incR2,4),F=round(instr_F,1),
      PTB_OR=round(prs_ptb_OR,3),PTB_p=signif(prs_ptb_p,2),MR_OR=round(mr_OR,3),MR_p=signif(mr_p,2))])
cat("\nfull table ->", file.path(OUT,"pooled_prs_mr_results.txt"), "\n")
