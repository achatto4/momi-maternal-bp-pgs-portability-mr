## pca_lib.R
##
## Shared functions of the PCA scripts.
suppressMessages(library(data.table))

DMIN     <- 0.05
HC_THR   <- 0.4999

read_block <- function(f){
  d <- fread(f, showProgress=FALSE, colClasses=list(character=c("FID", "IID")))
  drop <- intersect(c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"), names(d))
  ids <- as.character(d$IID)
  m <- as.matrix(d[, setdiff(names(d), drop), with=FALSE]); storage.mode(m) <- "double"
  rownames(m) <- ids; m
}
hardcall <- function(x){ r <- round(x); r[!is.na(x) & abs(x - r) > HC_THR] <- NA_real_; r }

std_block <- function(X, isref, plat, maf_min=0.05, miss_max=0.05){
  pr  <- colMeans(X[isref, , drop=FALSE], na.rm=TRUE) / 2
  mr  <- colMeans(is.na(X[isref, , drop=FALSE]))
  mpl <- sapply(unique(plat), function(p) colMeans(is.na(X[plat == p, , drop=FALSE])))
  mpl <- if(is.null(dim(mpl))) matrix(mpl, nrow=ncol(X)) else mpl
  maxmiss <- apply(mpl, 1, max)
  keep <- is.finite(pr) & pmin(pr, 1 - pr) >= maf_min & mr <= miss_max & maxmiss <= miss_max
  Z <- NULL
  if(any(keep)){
    p <- pr[keep]; s <- sqrt(2 * p * (1 - p))
    Z <- sweep(sweep(X[, keep, drop=FALSE], 2, 2 * p, "-"), 2, s, "/")
    Z[is.na(Z)] <- 0
  }
  list(Z=Z, keep=keep, p_ref=pr, miss_ref=mr, miss_max_platform=maxmiss,
       n_missing_set_to_mean=if(any(keep)) sum(is.na(X[, keep, drop=FALSE])) else 0)
}

top_eigen <- function(G, k, full_max=7000L){
  n <- nrow(G)
  if(n <= full_max){ e <- eigen(G, symmetric=TRUE); return(list(values=e$values[seq_len(k)], vectors=e$vectors[, seq_len(k), drop=FALSE],
                                                             all_values=e$values, method="LAPACK dsyevr (full)")) }
  if(requireNamespace("RSpectra", quietly=TRUE)){
    e <- RSpectra::eigs_sym(G, k=k, which="LA", opts=list(tol=1e-12, maxitr=5000))
    o <- order(e$values, decreasing=TRUE)
    return(list(values=e$values[o], vectors=e$vectors[, o, drop=FALSE], all_values=NULL, method="RSpectra eigs_sym (Lanczos)"))
  }

  set.seed(1L); b <- 2L * k; V <- qr.Q(qr(matrix(rnorm(n * b), n, b))); old <- rep(Inf, k)
  for(it in 1:5000){ V <- qr.Q(qr(G %*% V)); H <- crossprod(V, G %*% V); h <- eigen((H + t(H)) / 2, symmetric=TRUE)
    val <- h$values[seq_len(k)]; if(max(abs(val - old) / abs(val)) < 1e-12) break; old <- val }
  W <- V %*% h$vectors
  list(values=h$values[seq_len(k)], vectors=W[, seq_len(k), drop=FALSE], all_values=NULL, method=sprintf("subspace iteration (%d iterations)", it))
}

fit_project <- function(G, ref, kappa, K=10L, Keig=20L){
  ref <- intersect(ref, rownames(G)); oth <- setdiff(rownames(G), ref)
  Grr <- G[ref, ref, drop=FALSE]
  Ke <- min(Keig, length(ref) - 1L)
  e  <- top_eigen(Grr, Ke)
  lam <- e$values; U <- e$vectors; gbar <- mean(diag(Grr))
  S_in <- sweep(U[, seq_len(K), drop=FALSE], 2, sqrt(pmax(lam[seq_len(K)], 0)), "*")
  rownames(S_in) <- ref
  P_raw <- NULL; P_cor <- NULL
  if(length(oth)){
    P_raw <- G[oth, ref, drop=FALSE] %*% sweep(U[, seq_len(K), drop=FALSE], 2, 1 / sqrt(lam[seq_len(K)]), "*")
    kx <- kappa[oth]; kx[is.na(kx)] <- 0
    D  <- 1 - outer(gbar * (1 - kx), lam[seq_len(K)], function(a, l) a / l)
    P_cor <- ifelse(D > DMIN, P_raw / D, NA_real_)
    rownames(P_raw) <- rownames(P_cor) <- oth
  }
  colnames(S_in) <- paste0("PC", seq_len(K))
  if(!is.null(P_cor)){ colnames(P_raw) <- colnames(P_cor) <- colnames(S_in) }
  list(scores=rbind(S_in, P_cor), in_sample=S_in, proj_raw=P_raw, proj_cor=P_cor, lambda=lam, U=U, gbar=gbar,
       trace=sum(diag(Grr)), frob2=sum(Grr^2), n_ref=length(ref), all_values=e$all_values, eig_method=e$method)
}

neff <- function(u){ u <- u / sqrt(sum(u^2)); c(n_eff=1 / sum(u^4), top5_share=sum(sort(u^2, decreasing=TRUE)[1:min(5, length(u))])) }
tw_seq <- function(lam, S1, S2, m){
  p <- m - 1; o <- numeric(length(lam)); np <- numeric(length(lam))
  for(k in seq_along(lam)){ n <- p * S1^2 / (p * S2 - S1^2); L <- p * lam[k] / S1
    mu <- (sqrt(n - 1) + sqrt(p))^2 / n; sg <- (sqrt(n - 1) + sqrt(p)) / n * (1 / sqrt(n - 1) + 1 / sqrt(p))^(1/3)
    o[k] <- (L - mu) / sg; np[k] <- n; S1 <- S1 - lam[k]; S2 <- S2 - lam[k]^2; p <- p - 1 }
  data.table(tw_stat=o, tw_p=tw_p(o), n_eff_markers=np) }
min_cancor <- function(A, B){ k <- ncol(A); min(cancor(A, B)$cor[seq_len(k)]) }
cancor_all <- function(A, B){ ok <- stats::complete.cases(A) & stats::complete.cases(B); cancor(A[ok, , drop=FALSE], B[ok, , drop=FALSE])$cor }
agree <- function(x, y, s){ ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  data.table(n=n, r=if(n >= 3) cor(x, y) else NA_real_, slope=if(n >= 3) unname(coef(lm(y ~ x))[2]) else NA_real_,
             sd_ratio=if(n >= 3) sd(y) / sd(x) else NA_real_, mean_diff_sd=if(n) mean(y - x) / s else NA_real_) }

separation <- function(v, g){
  ok <- is.finite(v); v <- v[ok]; g <- g[ok]; s <- sd(v); lv <- sort(unique(g))
  ss_b <- sum(sapply(lv, function(l) sum(g == l) * (mean(v[g == l]) - mean(v))^2)); ss_t <- sum((v - mean(v))^2)
  list(eta2=if(ss_t > 0) ss_b / ss_t else NA_real_,
       by_level=rbindlist(lapply(lv, function(l) data.table(level=l, n=sum(g == l), mean_sd_units=mean(v[g == l]) / s,
                                                            smd_vs_rest=if(sum(g != l) > 0) (mean(v[g == l]) - mean(v[g != l])) / s else NA_real_,
                                                            var_share=sum(v[g == l]^2) / sum(v^2), n_share=mean(g == l)))))
}

hudson_terms <- function(ac1, an1, ac2, an2){
  p1 <- ac1 / an1; p2 <- ac2 / an2
  N <- (p1 - p2)^2 - p1 * (1 - p1) / (an1 - 1) - p2 * (1 - p2) / (an2 - 1)
  D <- p1 * (1 - p2) + p2 * (1 - p1)
  data.table(N=N, D=D)
}
fst_jackknife <- function(N, D, block){
  ok <- is.finite(N) & is.finite(D); N <- N[ok]; D <- D[ok]; block <- block[ok]
  est <- sum(N) / sum(D); bl <- unique(block); g <- length(bl)
  sN <- tapply(N, block, sum)[as.character(bl)]; sD <- tapply(D, block, sum)[as.character(bl)]
  loo <- (sum(N) - sN) / (sum(D) - sD)
  se <- sqrt((g - 1) / g * sum((loo - mean(loo))^2))
  list(fst=est, se=se, n_variants=length(N), n_blocks=g)
}
