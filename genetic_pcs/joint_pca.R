## joint_pca.R
##
## Within-cohort joint-platform PCA: relationship matrix accumulated over variant blocks, axes from the unrelated
## reference, relatives projected with a relatedness-aware shrinkage correction, and the platform diagnostics.
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(here, "helpers", "ids.R")); source(file.path(here, "helpers", "tracy_widom.R")); source(file.path(here, "pca_lib.R"))
args <- commandArgs(TRUE); arg <- function(f, d=NULL){ i <- which(args == f); if(length(i)) args[i + 1] else d }
W <- arg("--work"); COH <- arg("--cohort"); ENC <- arg("--encoding"); AGG <- arg("--agg"); ACC <- arg("--previous-pcs")
NPC <- as.integer(arg("--npc", "10")); SEED <- as.integer(arg("--seed", "20260925"))
stopifnot(ENC %in% c("dosage", "hardcall")); dir.create(AGG, FALSE, TRUE)
TW_P_MAX <- 0.001; CC_MIN <- 0.70; NEFF_MIN <- 30; CVR_MIN <- 0.75
pcs <- paste0("PC", seq_len(NPC)); t0 <- Sys.time()
A_ <- function(dt, f) fwrite(cbind(data.table(cohort=COH), dt), file.path(AGG, f), sep="\t")
say <- function(...) cat(sprintf("[%s] %s: ", format(Sys.time(), "%H:%M:%S"), COH), ..., "\n", sep="")

Wm <- fread(file.path(W, "women.tsv"), colClasses=list(character="IID"))
Wm[, IID := normid(IID)]
Wm[, key := paste(IID, fifelse(record == "low-pass", "lp", "gsa"), sep="|")]
EX <- Wm[extra_gsa == TRUE, .(IID, key=paste(IID, "gsax", sep="|"), kappa=kappa_extra_main)]
keys <- c(Wm$key, EX$key)
if(anyDuplicated(Wm$IID)) stop("a woman appears twice")
plat_of <- setNames(c(ifelse(Wm$record == "low-pass", "low-pass", "GSA"), rep("GSA", nrow(EX))), keys)
ref_keys <- Wm[in_ref == TRUE, key]
isref <- keys %chin% ref_keys

blocks <- list.files(file.path(W, "blocks"), pattern="^lp_[0-9]+\\.raw$")
blocks <- blocks[order(as.integer(sub("^lp_([0-9]+)\\.raw$", "\\1", blocks)))]
if(!length(blocks)) stop("no exported blocks")
n <- length(keys); G <- matrix(0, n, n, dimnames=list(keys, keys)); Ghc <- if(ENC == "dosage") G else NULL
m_used <- 0L; m_used_hc <- 0L; VQ <- list(); PF <- list(); row_miss <- setNames(numeric(n), keys); n_var_seen <- 0L; nmiss_mean <- 0
isl <- plat_of[keys] == "low-pass"
for(b in blocks){
  k <- sub("^lp_", "", sub("\\.raw$", "", b))
  Xl <- read_block(file.path(W, "blocks", b)); rownames(Xl) <- paste(normid(rownames(Xl)), "lp", sep="|")
  fg <- file.path(W, "blocks", paste0("gsa_", k, ".raw"))
  Xg <- NULL
  if(file.exists(fg)){ Xg <- read_block(fg); gid <- normid(rownames(Xg))
    rownames(Xg) <- ifelse(paste(gid, "gsa", sep="|") %chin% keys, paste(gid, "gsa", sep="|"), paste(gid, "gsax", sep="|")) }
  vl <- sub("_[^_]+$", "", colnames(Xl)); al <- sub("^.*_", "", colnames(Xl))
  if(!is.null(Xg)){ vg <- sub("_[^_]+$", "", colnames(Xg)); ag <- sub("^.*_", "", colnames(Xg))
    if(!setequal(vl, vg)) stop("block ", k, ": the technologies carry different variants")
    o <- match(vl, vg); if(!identical(al, ag[o])) stop("block ", k, ": counted alleles differ between technologies (alignment failure)")
    Xg <- Xg[, o, drop=FALSE]; colnames(Xg) <- colnames(Xl) }
  X <- rbind(Xl, Xg)
  if(!all(keys %chin% rownames(X))) stop("block ", k, ": ", sum(!keys %chin% rownames(X)), " records missing")
  X <- X[keys, , drop=FALSE]
  if(ENC == "hardcall") X <- hardcall(X)
  row_miss <- row_miss + rowSums(is.na(X)); n_var_seen <- n_var_seen + ncol(X)
  s <- std_block(X, isref, plat_of[keys])

  PF[[k]] <- data.table(p_lp=colMeans(X[isl, , drop=FALSE], na.rm=TRUE) / 2, n_lp=colSums(!is.na(X[isl, , drop=FALSE])),
                        p_gsa=colMeans(X[!isl, , drop=FALSE], na.rm=TRUE) / 2, n_gsa=colSums(!is.na(X[!isl, , drop=FALSE])), used=s$keep)
  VQ[[k]] <- data.table(block=k, exported=ncol(X), kept=sum(s$keep), maf_fail=sum(pmin(s$p_ref, 1 - s$p_ref) < 0.05, na.rm=TRUE),
                        miss_fail=sum(s$miss_ref > 0.05 | s$miss_max_platform > 0.05))
  if(!is.null(s$Z)){ G <- G + tcrossprod(s$Z); m_used <- m_used + ncol(s$Z); nmiss_mean <- nmiss_mean + s$n_missing_set_to_mean }
  if(ENC == "dosage"){ sh <- std_block(hardcall(X), isref, plat_of[keys]); if(!is.null(sh$Z)){ Ghc <- Ghc + tcrossprod(sh$Z); m_used_hc <- m_used_hc + ncol(sh$Z) } }
  rm(X, Xl, Xg, s); gc(verbose=FALSE)
}
if(m_used == 0L) stop("no variant survived QC")
G <- G / m_used; if(ENC == "dosage") Ghc <- Ghc / m_used_hc
say(sprintf("G accumulated over %d variants (%d blocks), %d records, %.0f s", m_used, length(blocks), n, as.numeric(difftime(Sys.time(), t0, units="secs"))))

kap <- setNames(c(Wm$kappa_main, EX$kappa), keys)
F <- fit_project(G, ref_keys, kap, K=NPC)
S <- F$scores[Wm$key, , drop=FALSE]
PCW <- data.table(IID=Wm$IID, cohort=COH, technology=Wm$technology, record=Wm$record, in_reference=Wm$in_ref,
                  kappa=fifelse(Wm$in_ref, 0, Wm$kappa_main), as.data.table(S))
fwrite(PCW, file.path(W, "pcs_joint.tsv"), sep="\t")
lam <- F$lambda; tot <- F$trace
self_r <- sapply(seq_len(NPC), function(j) cor(F$in_sample[, j], (G[ref_keys, ref_keys] %*% F$U[, j]) / sqrt(lam[j])))

tw <- tw_seq(lam, F$trace, F$frob2, F$n_ref)
nf <- t(sapply(seq_len(NPC), function(j) neff(F$U[, j])))
set.seed(SEED + 1L); fold <- sample(rep_len(1:5, length(ref_keys))); set.seed(SEED + 2L); half <- sample(rep_len(1:2, length(ref_keys)))
CV <- rbindlist(lapply(1:5, function(f){ tr <- ref_keys[fold != f]; te <- ref_keys[fold == f]
  e <- top_eigen(G[tr, tr], NPC); Str <- sweep(e$vectors, 2, sqrt(e$values), "*"); Pte <- G[te, tr] %*% sweep(e$vectors, 2, 1 / sqrt(e$values), "*")
  data.table(fold=f, pc=pcs, n_train=length(tr), n_test=length(te), eigenvalue=e$values, d_emp=apply(Pte, 2, sd) / apply(Str, 2, sd),
             d_theory=1 - mean(diag(G[tr, tr])) / e$values) }))
cvs <- CV[, .(d_cv=mean(d_emp), d_theory_cv=mean(d_theory), cv_ratio=mean(d_emp / d_theory)), by=pc]
hs <- lapply(1:2, function(h){ hh <- ref_keys[half == h]; e <- top_eigen(G[hh, hh], NPC)
  G[ref_keys, hh] %*% sweep(e$vectors, 2, 1 / sqrt(e$values), "*") })
cc <- sapply(seq_len(NPC), function(k) min_cancor(hs[[1]][, 1:k, drop=FALSE], hs[[2]][, 1:k, drop=FALSE]))
AX <- data.table(pc=pcs, eigenvalue=lam[seq_len(NPC)], pct_variance=100 * lam[seq_len(NPC)] / tot, tw_stat=tw$tw_stat[seq_len(NPC)],
                 tw_p=tw$tw_p[seq_len(NPC)], split_half_cc=cc, n_eff=nf[, "n_eff"], top5_share=nf[, "top5_share"],
                 d_theory_full=1 - F$gbar / lam[seq_len(NPC)], self_projection_r=self_r)
AX <- merge(AX, cvs, by="pc", sort=FALSE); AX <- AX[match(pcs, pc)]
AX[, qualifies := tw_p < TW_P_MAX & split_half_cc >= CC_MIN & n_eff >= NEFF_MIN & cv_ratio >= CVR_MIN]
kstar <- { q <- AX$qualifies; if(!q[1]) 0L else { w <- which(!q); if(length(w)) w[1] - 1L else length(q) } }
AX[, `:=`(informative=seq_len(.N) <= kstar, n_informative=kstar, n_reference=F$n_ref, mean_grm_diag=F$gbar, eig_method=F$eig_method)]
A_(AX, "agg_eigen.tsv"); A_(CV, "agg_cv_folds.tsv")
A_(data.table(pc=paste0("PC", seq_along(lam)), eigenvalue=lam, pct_variance=100 * lam / tot, tw), "agg_tracy_widom_1_20.tsv")

sepr <- function(dt, scope){ rbindlist(lapply(pcs, function(p){
  a <- separation(dt[[p]], dt$record); b <- separation(dt[[p]], dt$technology)
  rbind(cbind(data.table(scope=scope, pc=p, grouping="record platform", eta2=a$eta2), a$by_level),
        cbind(data.table(scope=scope, pc=p, grouping="technology", eta2=b$eta2), b$by_level)) })) }
SEP <- rbind(sepr(PCW, "all women"), sepr(PCW[in_reference == TRUE], "reference women (in-sample)"))
A_(SEP, "agg_separation.tsv")

sp <- if(nrow(EX)) diag(G[paste(EX$IID, "lp", sep="|"), EX$key, drop=FALSE]) else numeric(0)
A_(data.table(n_dual=length(sp), median_same_person_grm=if(length(sp)) median(sp) else NA_real_,
              p05=if(length(sp)) unname(quantile(sp, 0.05)) else NA_real_, min=if(length(sp)) min(sp) else NA_real_,
              median_self_grm_lowpass=median(diag(G)[Wm[record == "low-pass", key]]),
              median_self_grm_gsa=if(Wm[record == "GSA", .N]) median(diag(G)[Wm[record == "GSA", key]]) else NA_real_), "agg_same_person_grm.tsv")
DU <- NULL
if(nrow(EX)){
  xm <- F$scores[EX$key, , drop=FALSE]; lm_ <- F$scores[paste(EX$IID, "lp", sep="|"), , drop=FALSE]
  DU <- rbindlist(lapply(seq_len(NPC), function(j) cbind(data.table(fit="main (extra GSA record projected; kappa includes own low-pass record)", pc=pcs[j]),
                                                     agree(lm_[, j], xm[, j], sd(F$in_sample[, j])))))
  refnd <- Wm[in_ref_nodual == TRUE, key]
  if(length(refnd) > NPC + 10){
    kl <- setNames(c(Wm$kappa_ldo, Wm[extra_gsa == TRUE, kappa_ldo]), keys)
    Fl <- fit_project(G[c(refnd, paste(EX$IID, "lp", sep="|"), EX$key), c(refnd, paste(EX$IID, "lp", sep="|"), EX$key)], refnd, kl, K=NPC)
    DU <- rbind(DU, rbindlist(lapply(seq_len(NPC), function(j) cbind(data.table(fit="leave-dual-out (both records projected and corrected)", pc=pcs[j]),
                agree(Fl$scores[paste(EX$IID, "lp", sep="|"), j], Fl$scores[EX$key, j], sd(Fl$in_sample[, j]))))))
  }
  fwrite(data.table(IID=EX$IID, xm), file.path(W, "extra_gsa_projection.tsv"), sep="\t")
}
A_(if(is.null(DU)) data.table(fit=character(0)) else DU, "agg_dual_agreement.tsv")

Gs <- G[Wm$key, Wm$key]; diag(Gs) <- NA
hi <- which(Gs >= 0.80, arr.ind=TRUE); hi <- hi[hi[, 1] < hi[, 2], , drop=FALSE]
A_(data.table(selected_records=nrow(Wm), pairs_grm_ge_0.80=nrow(hi), max_offdiag_grm=max(Gs, na.rm=TRUE),
              same_person_pairs_ge_0.80=sum(sp >= 0.80), same_person_pairs=length(sp), duplicated_women=anyDuplicated(Wm$IID) > 0), "agg_duplicates.tsv")
rm(Gs); gc(verbose=FALSE)
VQd <- rbindlist(VQ)
PFd <- rbindlist(PF)[used == TRUE & is.finite(p_lp) & is.finite(p_gsa)]; dpf <- PFd$p_lp - PFd$p_gsa
pbf <- (PFd$p_lp * PFd$n_lp + PFd$p_gsa * PFd$n_gsa) / (PFd$n_lp + PFd$n_gsa)
A_(data.table(encoding=ENC, variants=nrow(PFd), lp_records=sum(isl), gsa_records=sum(!isl), r_frequency=cor(PFd$p_lp, PFd$p_gsa),
              pct_beyond_0.10=100 * mean(abs(dpf) > 0.10),
              pct_expected_by_sampling=100 * mean(2 * pnorm(-0.10 / sqrt(pbf * (1 - pbf) * (1 / (2 * PFd$n_lp) + 1 / (2 * PFd$n_gsa)))), na.rm=TRUE),
              mean_signed_diff_lp_minus_gsa=mean(dpf), median_abs_diff=median(abs(dpf)), p99_abs_diff=unname(quantile(abs(dpf), 0.99))),
   "agg_platform_frequency_pca_variants.tsv")
A_(data.table(encoding=ENC, blocks=length(blocks), exported=sum(VQd$exported), used=m_used, dropped_maf=sum(VQd$maf_fail),
              dropped_missing=sum(VQd$miss_fail), missing_values_set_to_mean=nmiss_mean, used_hardcall_version=m_used_hc), "agg_variant_qc.tsv")
fm <- row_miss / n_var_seen
A_(rbindlist(lapply(c("low-pass", "GSA"), function(p){ v <- fm[plat_of[keys] == p]
  data.table(platform=p, records=length(v), median_missing=if(length(v)) median(v) else NA_real_, p95_missing=if(length(v)) unname(quantile(v, 0.95)) else NA_real_,
             max_missing=if(length(v)) max(v) else NA_real_) })), "agg_missingness.tsv")

A_(rbind(data.table(group="reference (in-sample)", n=F$n_ref, kappa_median=0, missing_any_pc1_5=0L),
         data.table(group="related women (projected, corrected)", n=Wm[in_ref == FALSE, .N], kappa_median=if(Wm[in_ref == FALSE, .N]) median(Wm[in_ref == FALSE, kappa_main]) else NA_real_,
                    missing_any_pc1_5=sum(!stats::complete.cases(S[Wm[in_ref == FALSE, key], 1:5, drop=FALSE])))), "agg_projection.tsv")

if(!is.null(ACC) && file.exists(ACC)){
  AP <- if(grepl("\\.rds$", ACC)) as.data.table(readRDS(ACC)) else fread(ACC, colClasses=list(character="IID"))
  AP[, IID := normid(IID)]; ap <- intersect(paste0("PC", 1:10), names(AP))
  M <- merge(PCW[, c("IID", pcs), with=FALSE], unique(AP[, c("IID", ap), with=FALSE], by="IID"), by="IID", suffixes=c("", "_acc"))
  A5 <- as.matrix(M[, paste0("PC", 1:5), with=FALSE]); B5 <- as.matrix(M[, paste0("PC", 1:5, "_acc"), with=FALSE])
  ccr <- cancor_all(A5, B5)
  R <- cor(A5, as.matrix(M[, paste0(ap, "_acc"), with=FALSE]), use="pairwise.complete.obs")
  A_(data.table(n=nrow(M), pc=paste0("PC", 1:5), canonical_correlation=ccr, same_axis_abs_r=abs(diag(R[, 1:5])),
                best_previous_axis=ap[apply(abs(R), 1, which.max)], best_abs_r=apply(abs(R), 1, max)), "agg_vs_previous.tsv")
}

if(ENC == "dosage"){
  Fh <- fit_project(Ghc, ref_keys, kap, K=NPC); Sh <- Fh$scores[Wm$key, 1:5, drop=FALSE]
  ok <- stats::complete.cases(S[, 1:5]) & stats::complete.cases(Sh)
  A_(data.table(variants_dosage=m_used, variants_hardcall=m_used_hc, pc=paste0("PC", 1:5), canonical_correlation=cancor_all(S[ok, 1:5], Sh[ok, ]),
                same_axis_abs_r=abs(diag(cor(S[ok, 1:5], Sh[ok, ])))), "agg_encoding_sensitivity.tsv")
}

if(requireNamespace("png", quietly=TRUE) || capabilities("png")){
  png(file.path(AGG, sprintf("diag_%s_joint_pcs.png", COH)), width=1800, height=620, res=150)
  par(mfrow=c(1, 3), mar=c(4, 4, 2.5, 1), cex=0.8)
  cols <- c("low-pass WGS only"="#3b6ea5", "GSA only"="#d1495b", "both"="#e8a33d", "none"="grey50")
  o <- order(PCW$technology != "low-pass WGS only")
  for(pp in list(c(1, 2), c(3, 4), c(5, 1))){
    plot(S[o, pp[1]], S[o, pp[2]], pch=16, cex=0.35, col=adjustcolor(cols[PCW$technology[o]], 0.6),
         xlab=sprintf("PC%d (%.2f%%)", pp[1], 100 * lam[pp[1]] / tot), ylab=sprintf("PC%d (%.2f%%)", pp[2], 100 * lam[pp[2]] / tot),
         main=sprintf("%s: PC%d vs PC%d", COH, pp[1], pp[2]))
  }
  legend("topright", legend=names(cols)[1:3], col=cols[1:3], pch=16, bty="n", cex=0.9)
  dev.off()
}
say(sprintf("done: %d women (%d reference), PC1 %.2f%% of variance, %d informative axes, %.0f s", nrow(Wm), F$n_ref,
            100 * lam[1] / tot, kstar, as.numeric(difftime(Sys.time(), t0, units="secs"))))
