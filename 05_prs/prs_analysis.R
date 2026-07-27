#!/usr/bin/env Rscript
# ============================================================
# prs_analysis.R  (Stage 05: PRS association + PRS-based MR)
# For each BP trait x ancestry PGS:
#   (1) instrument strength : BP ~ PRS_std + covars   -> beta_bp, incremental R2, F
#   (2) PRS -> PTB          : PTB ~ PRS_std + covars   -> logOR, OR, p, AUC
#   (3) PRS-based MR (Wald) : BP->PTB = beta_ptb/beta_bp, delta-method SE, p, OR
# Run COMBINED (all cohorts, PRS standardized within site, +SITE+PCs) and PER-COHORT.
#
# Inputs (from earlier stages):
#   --prs-dir   : has {SBP,DBP}_{EUR,SAS}.sscore  (col SCORE1_AVG)
#   --pheno-dir : has mothers_{SBP,DBP}_<window>_final.pheno + mothers_PTB_NEW_final.pheno
#   --covar     : covariate_table_analytic_mothers.txt (PARTICIPANT_ID, SITE, PW_AGE, ...)
#   --eigenvec  : pca_results.eigenvec
#   --out-dir   : outputs
#   [--bp-window latest] [--n-pcs 10] [--ptb-case-code 2]
# Outputs: prs_analytic_dataset.txt, prs_mr_results.txt (one row per trait x ancestry x scope)
# ============================================================
suppressMessages(library(data.table))
have_proc <- requireNamespace("pROC", quietly=TRUE)

args <- commandArgs(trailingOnly=TRUE)
ga <- function(f,d=NULL){i<-which(args==f); if(length(i)) args[i+1] else d}
PRSDIR<-ga("--prs-dir"); PHDIR<-ga("--pheno-dir"); COV<-ga("--covar")
EVEC<-ga("--eigenvec"); OUT<-ga("--out-dir"); WIN<-ga("--bp-window","latest")
NPC<-as.integer(ga("--n-pcs","10")); PTBCASE<-as.integer(ga("--ptb-case-code","2"))
stopifnot(!is.null(PRSDIR),!is.null(PHDIR),!is.null(COV),!is.null(EVEC),!is.null(OUT))
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

## ---- load PRS scores ----
read_score <- function(tag){
  f <- file.path(PRSDIR, paste0(tag,".sscore"))
  if(!file.exists(f)) return(NULL)
  d <- fread(f); setnames(d, 1, "IID")
  d[, .(IID, score=SCORE1_AVG)][, setNames(.(IID, score), c("IID", tag))]
}
tags <- c("SBP_EUR","DBP_EUR","SBP_SAS","DBP_SAS")
prs <- Reduce(function(a,b) merge(a,b,by="IID",all=TRUE),
              Filter(Negate(is.null), lapply(tags, read_score)))

## ---- phenotypes ----
read_ph <- function(fn,name){
  f<-file.path(PHDIR,fn); d<-fread(f)
  d[,.(IID, val=suppressWarnings(as.numeric(ifelse(PHENO=="NA",NA,PHENO))))][, setNames(.(IID,val),c("IID",name))]
}
sbp <- read_ph(sprintf("mothers_SBP_%s_final.pheno",WIN),"SBP")
dbp <- read_ph(sprintf("mothers_DBP_%s_final.pheno",WIN),"DBP")
ptb <- fread(file.path(PHDIR,"mothers_PTB_NEW_final.pheno"))
ptb <- ptb[,.(IID, PTB=ifelse(PHENO=="NA",NA, as.integer(as.integer(PHENO)==PTBCASE)))]  # 1=case,0=control

## ---- covariates + PCs ----
cov <- fread(COV)
cov[, PW_AGE := suppressWarnings(as.numeric(PW_AGE))]
ev <- fread(EVEC); setnames(ev,1:2,c("FID","IID"))
pcn <- grep("^PC", names(ev), value=TRUE); if(!length(pcn)){setnames(ev,3:ncol(ev),paste0("PC",seq_len(ncol(ev)-2))); pcn<-grep("^PC",names(ev),value=TRUE)}
pcn <- pcn[seq_len(min(NPC,length(pcn)))]
ev <- ev[, c("IID",pcn), with=FALSE]

## ---- assemble ----
d <- Reduce(function(a,b) merge(a,b,by="IID",all.x=TRUE),
            list(prs, sbp, dbp, ptb, cov[,.(IID=PARTICIPANT_ID,SITE,ANCESTRY,PW_AGE)], ev))
fwrite(d, file.path(OUT,"prs_analytic_dataset.txt"), sep="\t", na="NA")
cat(sprintf("assembled %d mothers; cohorts: %s\n", nrow(d), paste(sort(unique(d$SITE)),collapse=", ")))

## ---- standardize each PRS within SITE ----
zin <- function(x){ if(sd(x,na.rm=TRUE)==0||is.na(sd(x,na.rm=TRUE))) return(rep(NA_real_,length(x))); (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE) }
for(tg in intersect(tags,names(d))) d[, (paste0(tg,"_z")) := zin(get(tg)), by=SITE]

pcterm <- paste(pcn, collapse=" + ")

analyze <- function(dd, trait, anc, scope){
  tg <- paste0(trait,"_",anc,"_z")
  if(!(tg %in% names(dd))) return(NULL)
  dd <- dd[!is.na(get(tg)) & !is.na(get(trait)) & !is.na(PTB) & !is.na(PW_AGE)]
  if(nrow(dd) < 30) return(NULL)
  ncase<-sum(dd$PTB==1); nctrl<-sum(dd$PTB==0)
  if(ncase<10 || nctrl<10) return(NULL)
  base <- if(scope=="combined" && length(unique(dd$SITE))>1) paste0("PW_AGE + ",pcterm," + factor(SITE)") else paste0("PW_AGE + ",pcterm)
  # (1) instrument: BP ~ PRS
  m_bp <- lm(as.formula(sprintf("%s ~ %s + %s", trait, tg, base)), data=dd)
  cb <- summary(m_bp)$coefficients
  b_bp<-cb[tg,1]; se_bp<-cb[tg,2]
  m0  <- lm(as.formula(sprintf("%s ~ %s", trait, base)), data=dd)
  r2_full<-summary(m_bp)$r.squared; r2_null<-summary(m0)$r.squared; incR2<-r2_full-r2_null
  Fstat <- (b_bp/se_bp)^2
  # (2) PRS -> PTB
  m_ptb <- glm(as.formula(sprintf("PTB ~ %s + %s", tg, base)), data=dd, family=binomial())
  cp <- summary(m_ptb)$coefficients
  b_ptb<-cp[tg,1]; se_ptb<-cp[tg,2]; p_ptb<-cp[tg,4]
  auc <- NA_real_
  if(have_proc) auc <- tryCatch(as.numeric(pROC::auc(dd$PTB, predict(m_ptb,type="response"), quiet=TRUE)), error=function(e) NA_real_)
  # (3) MR Wald ratio (BP -> PTB)
  ratio <- b_ptb/b_bp
  se_ratio <- sqrt(se_ptb^2/b_bp^2 + b_ptb^2*se_bp^2/b_bp^4)   # delta method
  z <- ratio/se_ratio; p_mr <- 2*pnorm(-abs(z))
  data.table(scope=scope, trait=trait, ancestry=anc, n=nrow(dd), n_case=ncase, n_ctrl=nctrl,
    instr_beta_bp=b_bp, instr_se=se_bp, instr_incR2=incR2, instr_F=Fstat,
    prs_ptb_logOR=b_ptb, prs_ptb_OR=exp(b_ptb), prs_ptb_p=p_ptb, prs_ptb_AUC=auc,
    mr_BPtoPTB_logOR=ratio, mr_se=se_ratio, mr_OR=exp(ratio), mr_p=p_mr)
}

res <- list()
for(trait in c("SBP","DBP")) for(anc in c("EUR","SAS")){
  res[[paste("comb",trait,anc)]] <- analyze(d, trait, anc, "combined")
  for(co in sort(unique(d$SITE))){
    r <- analyze(d[SITE==co], trait, anc, co); if(!is.null(r)) res[[paste(co,trait,anc)]] <- r
  }
}
out <- rbindlist(Filter(Negate(is.null), res))
fwrite(out, file.path(OUT,"prs_mr_results.txt"), sep="\t")

cat("\n=== COMBINED results (instrument strength, PRS->PTB, MR) ===\n")
print(out[scope=="combined", .(trait,ancestry,n,incR2=round(instr_incR2,4),F=round(instr_F,1),
        PTB_OR=round(prs_ptb_OR,3),PTB_p=signif(prs_ptb_p,2),AUC=round(prs_ptb_AUC,3),
        MR_OR=round(mr_OR,3),MR_p=signif(mr_p,2))])
cat("\nFull results -> ", file.path(OUT,"prs_mr_results.txt"), "\n")
