## cross_cohort_pca.R
##
## Descriptive cross-cohort PCA (figure data as 2-D bins) and pairwise Hudson FST (ratio of averages, sample sizes in
## chromosomes, delete-one-block jackknife over 5-Mb windows).
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(here, "helpers", "ids.R")); source(file.path(here, "pca_lib.R"))
args <- commandArgs(TRUE); arg <- function(f, d=NULL){ i <- which(args == f); if(length(i)) args[i + 1] else d }
RUN <- arg("--rundir"); ENC <- arg("--encoding"); COHS <- strsplit(arg("--cohorts"), " ")[[1]]; SEED <- as.integer(arg("--seed", "20260925"))
X <- file.path(RUN, "work", "cross"); AGG <- file.path(RUN, "aggregate", "cross"); dir.create(AGG, FALSE, TRUE)
NBIN <- 400L; K <- 10L; t0 <- Sys.time()
say <- function(...) cat(sprintf("[%s] cross: ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep="")
SHORT <- c("AMANHI-Bangladesh"="Sylhet", "AMANHI-Pakistan"="Karachi", "GAPPS-Bangladesh"="Matlab", "AMANHI-Pemba"="Pemba", "GAPPS-Zambia"="Lusaka")

Wm <- rbindlist(lapply(COHS, function(c){ w <- fread(file.path(RUN, "work", c, "women.tsv"), colClasses=list(character="IID")); w[, cohort := c]; w }), fill=TRUE)
Wm[, IID := normid(IID)]; if(anyDuplicated(Wm$IID)) stop("a woman appears in two cohorts")
Wm[, key := paste(IID, fifelse(record == "low-pass", "lp", "gsa"), sep="|")]
keys <- Wm$key; n <- length(keys); isref <- Wm$in_ref; grp <- paste(Wm$cohort, Wm$record)
ref_keys <- keys[isref]; kap <- setNames(fifelse(Wm$in_ref, 0, Wm$kappa_main), keys)

nb <- length(list.files(file.path(X, "blocks", COHS[1]), pattern="^lp_[0-9]+\\.raw$"))
G <- matrix(0, n, n, dimnames=list(keys, keys)); m_used <- 0L; ACC <- list()
for(k in seq_len(nb)){
  parts <- list()
  for(c in COHS){
    xl <- read_block(file.path(X, "blocks", c, sprintf("lp_%d.raw", k))); rownames(xl) <- paste(normid(rownames(xl)), "lp", sep="|")
    fg <- file.path(X, "blocks", c, sprintf("gsa_%d.raw", k)); xg <- NULL
    if(file.exists(fg)){ xg <- read_block(fg); rownames(xg) <- paste(normid(rownames(xg)), "gsa", sep="|")
      vl <- sub("_[^_]+$", "", colnames(xl)); vg <- sub("_[^_]+$", "", colnames(xg))
      if(!setequal(vl, vg)) stop("block ", k, " ", c, ": technologies carry different variants")
      o <- match(vl, vg); if(!identical(sub("^.*_", "", colnames(xl)), sub("^.*_", "", colnames(xg))[o])) stop("counted alleles differ: ", c)
      xg <- xg[, o, drop=FALSE]; colnames(xg) <- colnames(xl) }
    parts[[c]] <- rbind(xl, xg)
  }
  cn <- colnames(parts[[1]]); for(c in COHS) if(!identical(colnames(parts[[c]]), cn)) stop("block ", k, ": cohorts carry different variants or alleles")
  Xb <- do.call(rbind, parts); if(!all(keys %chin% rownames(Xb))) stop("block ", k, ": records missing"); Xb <- Xb[keys, , drop=FALSE]
  Xh <- hardcall(Xb); Xu <- if(ENC == "hardcall") Xh else Xb
  s <- std_block(Xu, isref, grp)
  if(!is.null(s$Z)){ G <- G + tcrossprod(s$Z); m_used <- m_used + ncol(s$Z) }
  kv <- which(s$keep)
  if(length(kv)){
    vid <- sub("_[^_]+$", "", cn[kv]); pos <- suppressWarnings(as.numeric(sub(".*:", "", vid))); chr <- sub("^chr", "", sub(":.*", "", vid))
    for(c in COHS){ for(sub_ in c("all", "reference")){
      rr <- Wm$cohort == c & (sub_ == "all" | Wm$in_ref)
      h <- Xh[rr, kv, drop=FALSE]; d <- Xb[rr, kv, drop=FALSE]
      ACC[[length(ACC) + 1]] <- data.table(block=k, variant=vid, chr=chr, pos=pos, cohort=c, subset=sub_,
                                            ac_hc=colSums(h, na.rm=TRUE), an_hc=2 * colSums(!is.na(h)),
                                            ac_dos=colSums(d, na.rm=TRUE), an_dos=2 * colSums(!is.na(d))) } }
  }
  rm(parts, Xb, Xh, Xu, s); gc(verbose=FALSE)
}
if(m_used == 0L) stop("no cross-cohort variant survived QC")
G <- G / m_used
say(sprintf("cross-cohort G: %d women x %d variants, %.0f s", n, m_used, as.numeric(difftime(Sys.time(), t0, units="secs"))))

F <- fit_project(G, ref_keys, kap, K=K, Keig=K)
S <- F$scores[keys, , drop=FALSE]
PC <- data.table(IID=Wm$IID, cohort=Wm$cohort, technology=Wm$technology, record=Wm$record, in_reference=Wm$in_ref, as.data.table(S))
fwrite(PC, file.path(X, "cross_pcs.tsv"), sep="\t")
pct <- 100 * F$lambda / F$trace
fwrite(data.table(pc=paste0("PC", seq_len(K)), eigenvalue=F$lambda, pct_variance=pct, n_reference=F$n_ref, n_women=n, n_variants=m_used,
                  trace_reference=F$trace, eig_method=F$eig_method), file.path(AGG, "agg_cross_eigen.tsv"), sep="\t")
GS <- PC[, c(list(n=.N), lapply(.SD, mean), lapply(.SD, sd)), by=.(cohort, technology), .SDcols=c("PC1", "PC2", "PC3")]
setnames(GS, c("cohort", "technology", "n", paste0(c("PC1", "PC2", "PC3"), "_mean"), paste0(c("PC1", "PC2", "PC3"), "_sd")))
fwrite(GS[order(match(cohort, COHS), technology)], file.path(AGG, "agg_cross_group_summary.tsv"), sep="\t")

rx <- range(S[, 1], na.rm=TRUE); ry <- range(S[, 2], na.rm=TRUE); px <- diff(rx) * 0.02; py <- diff(ry) * 0.02
bx <- seq(rx[1] - px, rx[2] + px, length.out=NBIN + 1); by <- seq(ry[1] - py, ry[2] + py, length.out=NBIN + 1)
PC[, `:=`(bin_x=findInterval(PC1, bx, all.inside=TRUE), bin_y=findInterval(PC2, by, all.inside=TRUE))]
fwrite(PC[, .(count=.N), by=.(cohort, technology, bin_x, bin_y)][order(cohort, technology, bin_x, bin_y)], file.path(AGG, "agg_cross_bins.tsv"), sep="\t")
fwrite(data.table(axis=c("PC1", "PC2"), lower=c(bx[1], by[1]), upper=c(bx[NBIN + 1], by[NBIN + 1]), n_bins=NBIN,
                  pct_variance=pct[1:2]), file.path(AGG, "agg_cross_bin_edges.tsv"), sep="\t")

COL <- c("AMANHI-Bangladesh"="#0072B2", "AMANHI-Pakistan"="#56B4E9", "GAPPS-Bangladesh"="#009E73", "AMANHI-Pemba"="#E69F00", "GAPPS-Zambia"="#D55E00")
PCH <- c("low-pass WGS only"=16, "GSA only"=17, "both"=15)
draw <- function(dev_fun, file, ...){
  dev_fun(file, ...); par(mar=c(4.2, 4.4, 0.6, 11), mgp=c(2.6, 0.7, 0), las=1, cex.axis=0.9)
  plot(NA, xlim=range(bx), ylim=range(by), xaxs="i", yaxs="i", xlab=sprintf("PC1 (%.1f%% of variance)", pct[1]),
       ylab=sprintf("PC2 (%.1f%% of variance)", pct[2]), bty="l")
  o <- sample(seq_len(nrow(PC)))
  usr <- par("usr"); pin <- par("pin")
  if(requireNamespace("png", quietly=TRUE) && !identical(dev_fun, grDevices::png)){
    tf <- tempfile(fileext=".png"); grDevices::png(tf, width=round(pin[1] * 600), height=round(pin[2] * 600), res=600, bg="transparent")
    par(mar=c(0, 0, 0, 0)); plot.new(); plot.window(xlim=usr[1:2], ylim=usr[3:4], xaxs="i", yaxs="i")
    points(PC$PC1[o], PC$PC2[o], pch=PCH[PC$technology[o]], cex=0.45, col=adjustcolor(COL[PC$cohort[o]], 0.55)); grDevices::dev.off()
    rasterImage(png::readPNG(tf), usr[1], usr[3], usr[2], usr[4]); unlink(tf)
  } else points(PC$PC1[o], PC$PC2[o], pch=PCH[PC$technology[o]], cex=0.45, col=adjustcolor(COL[PC$cohort[o]], 0.55))
  par(xpd=NA)
  legend(usr[2] + 0.03 * diff(usr[1:2]), usr[4], title="Cohort", title.adj=0, legend=SHORT[COHS], pch=16, col=COL[COHS], bty="n", cex=0.85)
  legend(usr[2] + 0.03 * diff(usr[1:2]), usr[4] - 0.45 * diff(usr[3:4]), title="Genotyping technology", title.adj=0,
         legend=c("Low-pass WGS only", "GSA only", "Both"), pch=PCH, col="grey25", bty="n", cex=0.85)
  grDevices::dev.off() }
set.seed(SEED)
draw(grDevices::pdf, file.path(AGG, "cross_cohort_pca_cluster_render.pdf"), width=7.2, height=4.8, family="Helvetica")
set.seed(SEED)
draw(grDevices::png, file.path(AGG, "cross_cohort_pca_cluster_render.png"), width=7.2 * 300, height=4.8 * 300, res=300)

A <- rbindlist(ACC)
A[, window := paste(chr, floor(pos / 5e6), sep=":")]
FST <- list()
for(i in seq_along(COHS)) for(j in seq_along(COHS)) if(i < j) for(sub_ in c("all", "reference")) for(enc in c("hardcall", if(ENC == "dosage") "dosage")){
  a <- A[cohort == COHS[i] & subset == sub_]; b <- A[cohort == COHS[j] & subset == sub_]
  stopifnot(identical(a$variant, b$variant))
  h <- if(enc == "hardcall") hudson_terms(a$ac_hc, a$an_hc, b$ac_hc, b$an_hc) else hudson_terms(a$ac_dos, a$an_dos, b$ac_dos, b$an_dos)
  jk <- fst_jackknife(h$N, h$D, a$window)
  FST[[length(FST) + 1]] <- data.table(cohort_a=COHS[i], cohort_b=COHS[j], site_a=SHORT[COHS[i]], site_b=SHORT[COHS[j]], women=sub_,
                                       genotypes=enc, fst_hudson=jk$fst, jackknife_se=jk$se, n_variants=jk$n_variants, n_blocks_5mb=jk$n_blocks,
                                       median_women_a=median(a$an_hc) / 2, median_women_b=median(b$an_hc) / 2)
}
FST <- rbindlist(FST)
fwrite(FST, file.path(AGG, "agg_fst_pairwise.tsv"), sep="\t")
pri <- FST[site_a == "Sylhet" & site_b == "Matlab" & women == "all" & genotypes == "hardcall"]
fwrite(pri, file.path(AGG, "agg_fst_sylhet_matlab_primary.tsv"), sep="\t")
say(sprintf("PC1 %.1f%%, PC2 %.1f%% of variance; Sylhet-Matlab Hudson FST %.6f (jackknife SE %.6f, %d variants); %.0f s",
            pct[1], pct[2], pri$fst_hudson, pri$jackknife_se, pri$n_variants, as.numeric(difftime(Sys.time(), t0, units="secs"))))
