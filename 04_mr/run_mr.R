#!/usr/bin/env Rscript
# ============================================================
# run_mr.R  (Stage 04b)  -- SNP-based two-sample MR: BP (exposure) -> PTB (outcome)
# Exposure = external BP GWAS instruments (from prep_instruments.sh).
# Outcome  = internal PTB GWAS (plink2 .glm.logistic.hybrid), matched on chr:pos.
# Methods  = IVW, MR-Egger, weighted median, weighted mode + sensitivity + plots.
#
# Usage:
#   Rscript run_mr.R --instruments instruments_SBP.txt --ptb-gwas <glm.logistic.hybrid> \
#                    --trait SBP --out-dir DIR
# ============================================================
suppressMessages(library(data.table))
if(!requireNamespace("TwoSampleMR", quietly=TRUE))
  stop("TwoSampleMR not installed. Install (remotes::install_github('MRCIEU/TwoSampleMR')) or load a module that provides it.")
suppressMessages(library(TwoSampleMR))

args <- commandArgs(trailingOnly=TRUE)
ga <- function(f,d=NULL){i<-which(args==f); if(length(i)) args[i+1] else d}
INSTR<-ga("--instruments"); PTB<-ga("--ptb-gwas"); TRAIT<-ga("--trait","BP"); OUT<-ga("--out-dir")
stopifnot(!is.null(INSTR),!is.null(PTB),!is.null(OUT))
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

## ---- exposure ----
exp_df <- fread(INSTR)
exp_df[, Phenotype := TRAIT]
exposure <- format_data(as.data.frame(exp_df), type="exposure",
  snp_col="SNP", beta_col="BETA", se_col="SE",
  effect_allele_col="EA", other_allele_col="OA", eaf_col="EAF", pval_col="P",
  phenotype_col="Phenotype")
cat(sprintf("exposure instruments: %d\n", nrow(exposure)))

## ---- outcome: internal PTB GWAS at the instrument SNPs ----
ptb <- fread(PTB)
setnames(ptb, names(ptb), sub("^#","",names(ptb)))   # drop leading # on CHROM
# keep additive test rows if a TEST column exists
if("TEST" %in% names(ptb)) ptb <- ptb[TEST=="ADD"]
ptb <- ptb[ID %in% exp_df$SNP]
# beta = log(OR); SE column is LOG(OR)_SE
orcol <- intersect(c("OR"), names(ptb)); betacol <- intersect(c("BETA"), names(ptb))
secol  <- grep("LOG.OR._SE|^SE$", names(ptb), value=TRUE)[1]
ptb[, beta_out := if(length(betacol)) as.numeric(get(betacol)) else log(as.numeric(get(orcol)))]
ptb[, se_out := as.numeric(get(secol))]
ptb[, other := ifelse(A1==ALT, REF, ALT)]
out_df <- ptb[, .(SNP=ID, EA=A1, OA=other, BETA=beta_out, SE=se_out, P=as.numeric(P))]
out_df[, Phenotype := "PTB"]
cat(sprintf("outcome SNPs matched: %d\n", nrow(out_df)))
if(nrow(out_df) < 3) stop("Too few matched outcome SNPs for MR (need >=3).")

outcome <- format_data(as.data.frame(out_df), type="outcome",
  snp_col="SNP", beta_col="BETA", se_col="SE",
  effect_allele_col="EA", other_allele_col="OA", pval_col="P", phenotype_col="Phenotype")

## ---- harmonise + MR ----
dat <- harmonise_data(exposure, outcome, action=2)
fwrite(dat, file.path(OUT, sprintf("harmonised_%s_PTB.txt",TRAIT)), sep="\t")
nsnp <- sum(dat$mr_keep)
cat(sprintf("harmonised SNPs kept: %d\n", nsnp))

res  <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
het  <- tryCatch(mr_heterogeneity(dat), error=function(e) NULL)
ple  <- tryCatch(mr_pleiotropy_test(dat), error=function(e) NULL)
sin  <- tryCatch(mr_singlesnp(dat), error=function(e) NULL)
loo  <- tryCatch(mr_leaveoneout(dat), error=function(e) NULL)

res_or <- generate_odds_ratios(res)
fwrite(res_or, file.path(OUT, sprintf("mr_results_%s_PTB.txt",TRAIT)), sep="\t")
if(!is.null(het)) fwrite(het, file.path(OUT, sprintf("mr_heterogeneity_%s.txt",TRAIT)), sep="\t")
if(!is.null(ple)) fwrite(ple, file.path(OUT, sprintf("mr_pleiotropy_%s.txt",TRAIT)), sep="\t")
if(!is.null(sin)) fwrite(sin, file.path(OUT, sprintf("mr_singlesnp_%s.txt",TRAIT)), sep="\t")
if(!is.null(loo)) fwrite(loo, file.path(OUT, sprintf("mr_leaveoneout_%s.txt",TRAIT)), sep="\t")

## ---- plots ----
sv <- function(p, fn, w=6,h=6) tryCatch(ggplot2::ggsave(file.path(OUT,fn), p[[1]], width=w, height=h), error=function(e) NULL)
tryCatch(sv(mr_scatter_plot(res, dat),        sprintf("mr_scatter_%s_PTB.png",TRAIT)), error=function(e) NULL)
if(!is.null(sin)) tryCatch(sv(mr_forest_plot(sin),       sprintf("mr_forest_%s_PTB.png",TRAIT)), error=function(e) NULL)
if(!is.null(loo)) tryCatch(sv(mr_leaveoneout_plot(loo),  sprintf("mr_loo_%s_PTB.png",TRAIT)), error=function(e) NULL)
if(!is.null(sin)) tryCatch(sv(mr_funnel_plot(sin),       sprintf("mr_funnel_%s_PTB.png",TRAIT)), error=function(e) NULL)

cat("\n=== MR results:", TRAIT, "-> PTB ===\n")
print(res_or[, c("method","nsnp","b","se","pval","or","or_lci95","or_uci95")])
if(!is.null(ple)) cat(sprintf("Egger intercept p (pleiotropy): %.3g\n", ple$pval))
cat("outputs ->", OUT, "\n")
