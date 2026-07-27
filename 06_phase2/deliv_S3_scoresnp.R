#!/usr/bin/env Rscript
# ============================================================
# deliv_S3_scoresnp.R  [ B06 -> Supp Table S3 ]
# Score-SNP coverage, imputation quality, and cross-cohort allele-frequency divergence.
#
# WHY THIS EXISTS. Three separate audits have now found the African cohorts measurably
# disadvantaged relative to the South-Asian ones:
#   * variant recovery      80.0% (AFR) vs 85.1% (SAS)            [audit_B04 / S1b]
#   * cross-platform PRS reliability  0.853 vs 0.926               [audit_B02 section K]
#   * loss to follow-up     9.8% in Zambia vs 0.2-3.2% elsewhere   [audit_B01]
# Each is individually far too small to explain a two-to-fourfold R2 gap, and that has been
# checked twice. S3 supplies the missing piece -- imputation INFO at score SNPs, and how far
# allele frequencies at those SNPs diverge between cohorts -- and asks whether these are ONE
# phenomenon or three. If AF divergence at score SNPs tracks the R2 gradient, that is a
# specific, measurable mechanism for poor portability rather than "ancestry" in the abstract.
#
# THE INFORMATIVE COMPARISON IS NOT AFR-vs-SAS. That divergence is expected and tells us
# nothing. The diagnostic contrast is AMANHI-Bangladesh vs GAPPS-Bangladesh: same country,
# same ancestry, so AF divergence there should be ~0. If it is not, the problem is
# genotyping, not populations.
#
# Needs --outroot (per-cohort *_clean filesets) and --scores-dir (harmonised PGS files);
# SKIPs cleanly if either is absent so the driver never breaks.
#
# Writes: tables/S3_scoresnp.tsv, tables/S3b_pairwise_fst.tsv
#   Rscript deliv_S3_scoresnp.R --outroot DIR --scores-dir DIR [--plink2 plink2]
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

OUTROOT <- momi_arg("--outroot",    Sys.getenv("MOMI_OUTROOT"))
SCORES  <- momi_arg("--scores-dir", Sys.getenv("MOMI_SCORES"))
PLINK2  <- momi_arg("--plink2", "plink2")
WORK    <- file.path(P$root, "qc", "s3_scoresnp"); dir.create(WORK, showWarnings=FALSE, recursive=TRUE)

## scores to profile: the CHOSEN instruments only. Doing all eight would multiply runtime
## without changing the answer -- these are the ones Table 2 and the MR actually use.
S3_SCORES <- unique(c(MOMI_EUR$SBP, MOMI_EUR$DBP, MOMI_SAS$SBP, MOMI_SAS$DBP))

momi_deliverable("S3_scoresnp", script="06_phase2/deliv_S3_scoresnp.R",
                 inputs="transfer_grid;genotypes;scores", P=P, stop_on_error=FALSE,
                 body=function(ctx){

  if(is.null(OUTROOT) || !nzchar(OUTROOT) || !dir.exists(OUTROOT))
    return(list(skip=TRUE, reason="no --outroot (per-cohort *_clean filesets) — see docs/paths_reference.md"))
  if(is.null(SCORES) || !nzchar(SCORES) || !dir.exists(SCORES))
    return(list(skip=TRUE, reason="no --scores-dir (harmonised PGS files)"))
  if(Sys.which(PLINK2)=="") return(list(skip=TRUE, reason="plink2 not on PATH"))

  ## ---- 1. score variant IDs, per score ----
  ## PGS Catalog harmonised files: hm_chr/hm_pos preferred, else chr_name/chr_position.
  ## ID convention "chr<N>:<pos>" matches the .bim/.pvar built by 00_preprocess.
  read_score_ids <- function(pid){
    f <- file.path(SCORES, paste0(pid, ".txt.gz"))
    if(!file.exists(f)) return(NULL)
    d <- tryCatch(fread(cmd=sprintf("zcat %s | grep -v '^#'", shQuote(f))), error=function(e) NULL)
    if(is.null(d) || !nrow(d)) return(NULL)
    cc <- if("hm_chr" %in% names(d)) "hm_chr" else "chr_name"
    pp <- if("hm_pos" %in% names(d)) "hm_pos" else "chr_position"
    if(!all(c(cc,pp) %in% names(d))) return(NULL)
    ids <- sprintf("chr%s:%s", sub("^chr","",as.character(d[[cc]])), as.character(d[[pp]]))
    unique(ids[!grepl("NA", ids)])
  }
  ## MEMORY (fix 2026-07-19): this used to build SC = list of ID vectors and HOLD IT for the
  ## whole run -- ~23.6M character strings across the four scores (7.35M + 7.35M + 4.81M +
  ## 4.14M), which is what OOM-killed B06 on an interactive allocation. Everything downstream
  ## needs only the COUNT (n_score_snps); the IDs themselves are consumed by plink2 --extract
  ## from the .ids file on disk. So write, count, and free -- one score in memory at a time.
  ## Cache the .ids too: rereading a 7M-line gzipped score file to regenerate an identical
  ## list is pure waste on a rerun.
  n_score <- integer(0)
  for(pid in S3_SCORES){
    idf <- file.path(WORK, paste0(pid, ".ids"))
    if(file.exists(idf) && file.size(idf) > 0){
      n_score[pid] <- as.integer(system(sprintf("wc -l < %s", shQuote(idf)), intern=TRUE))
      next
    }
    ids <- read_score_ids(pid)
    if(is.null(ids) || !length(ids)) next
    writeLines(ids, idf)
    n_score[pid] <- length(ids)
    rm(ids); gc(verbose=FALSE)
  }
  if(!length(n_score)) return(list(skip=TRUE, reason="no readable score files in --scores-dir"))
  SC_IDS <- names(n_score)
  cat(sprintf("score variant lists: %s\n",
              paste(sprintf("%s=%d", SC_IDS, n_score[SC_IDS]), collapse=", ")))

  ## ---- 2. per cohort x platform x score: plink2 --freq on the score-SNP subset ----
  ## cached: an existing .afreq is reused, so reruns are cheap
  ## mean imputation INFO depends on the FILESET ONLY, not on which score we are profiling,
  ## but the previous version recomputed it inside the score loop -- 4 greps of a
  ## multi-million-line .pvar where 1 would do, repeated on every rerun because (unlike the
  ## .afreq) it was never cached. Compute once per cohort x platform, memoised to disk.
  mean_info_for <- local({
    memo <- new.env(parent=emptyenv())
    function(pre, coh, pl){
      key <- paste(coh, pl, sep="|")
      if(!is.null(memo[[key]])) return(memo[[key]])
      cf <- file.path(WORK, sprintf("info__%s__%s.txt", coh, pl))
      if(file.exists(cf)){
        v <- suppressWarnings(as.numeric(readLines(cf, warn=FALSE)[1]))
        memo[[key]] <- if(length(v) && is.finite(v)) v else NA_real_
        return(memo[[key]])
      }
      pv <- paste0(pre, ".pvar"); v <- NA_real_
      if(file.exists(pv)){
        iv <- tryCatch(as.numeric(system(sprintf(
          "grep -v '^#' %s | grep -oE '(INFO|R2)=[0-9.]+' | cut -d= -f2 | awk '{s+=$1;n++} END{if(n>0) printf \"%%.4f\", s/n}'",
          shQuote(pv)), intern=TRUE)), error=function(e) NA_real_)
        if(length(iv) && is.finite(iv)) v <- iv
      }
      writeLines(as.character(v), cf)
      memo[[key]] <- v; v
    }
  })

  rows <- list(); AFL <- list()
  for(coh in MOMI_COH_ALL) for(pl in MOMI_PLATS) for(pid in SC_IDS){
    pre  <- file.path(OUTROOT, coh, paste0(pl, "_clean"))
    isp  <- file.exists(paste0(pre, ".pgen"))
    if(!isp && !file.exists(paste0(pre, ".bed"))) next
    tag  <- sprintf("%s__%s__%s", pid, coh, pl)
    af   <- file.path(WORK, paste0(tag, ".afreq"))
    if(!file.exists(af)){
      cmd <- sprintf("%s %s %s --extract %s --freq --out %s",
                     shQuote(PLINK2), if(isp) "--pfile" else "--bfile", shQuote(pre),
                     shQuote(file.path(WORK, paste0(pid,".ids"))), shQuote(file.path(WORK, tag)))
      system(paste(cmd, ">", shQuote(file.path(WORK, paste0(tag,".log"))), "2>&1"))
    }
    if(!file.exists(af)) next
    a <- tryCatch(fread(af), error=function(e) NULL); if(is.null(a) || !nrow(a)) next
    idc <- intersect(c("ID","#ID"), names(a))[1]
    fcc <- intersect(c("ALT_FREQS","ALT_FREQ"), names(a))[1]
    if(is.na(idc) || is.na(fcc)) next
    setnames(a, c(idc,fcc), c("ID","AF"))

    ## imputation INFO/R2 from the .pvar (dosage sets carry it) — memoised, see above
    meaninfo <- mean_info_for(pre, coh, pl)
    AFL[[tag]] <- data.table(score_id=pid, cohort=coh, platform=pl, ID=a$ID, AF=as.numeric(a$AF))
    rows[[tag]] <- data.table(
      score_id=pid, cohort=coh, cohort_anc=unname(MOMI_ANC[coh]), platform=pl,
      n_score_snps=n_score[[pid]], n_present=nrow(a),
      coverage_pct=round(100*nrow(a)/n_score[[pid]],1),
      mean_af=round(mean(a$AF, na.rm=TRUE),4), mean_info=meaninfo)
  }
  if(!length(rows)) return(list(skip=TRUE, reason="plink2 produced no .afreq — check ID format vs .bim/.pvar"))
  S3 <- rbindlist(rows, fill=TRUE)
  AF <- rbindlist(AFL, fill=TRUE)

  ## ---- 3. cross-cohort AF divergence + pairwise Hudson Fst (dosage platform) ----
  ## Hudson Fst, ratio-of-averages (Bhatia et al. 2013) — the estimator recommended when
  ## sample sizes differ between populations, which they do here (1441 to 3813).
  nsamp <- c(`AMANHI-Bangladesh`=1803, `AMANHI-Pakistan`=2004, `AMANHI-Pemba`=3813,
             `GAPPS-Bangladesh`=3567, `GAPPS-Zambia`=1441)
  fst_rows <- list()
  D <- AF[platform=="lpwgs_dosage" & score_id==MOMI_EUR$SBP]
  if(nrow(D)){
    W <- dcast(D, ID ~ cohort, value.var="AF")
    cohs <- setdiff(names(W), "ID")
    for(i in seq_along(cohs)) for(j in seq_along(cohs)) if(i < j){
      p1 <- W[[cohs[i]]]; p2 <- W[[cohs[j]]]
      n1 <- nsamp[[cohs[i]]]; n2 <- nsamp[[cohs[j]]]
      ok <- is.finite(p1) & is.finite(p2)
      num <- (p1-p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
      den <- p1*(1-p2) + p2*(1-p1)
      fst_rows[[paste(i,j)]] <- data.table(
        cohort_a=cohs[i], cohort_b=cohs[j],
        anc_a=unname(MOMI_ANC[cohs[i]]), anc_b=unname(MOMI_ANC[cohs[j]]),
        n_snps=sum(ok),
        fst_hudson=round(sum(num[ok], na.rm=TRUE)/sum(den[ok], na.rm=TRUE), 5),
        mean_abs_af_diff=round(mean(abs(p1-p2)[ok], na.rm=TRUE), 4))
    }
  }
  FST <- if(length(fst_rows)) rbindlist(fst_rows) else data.table()
  if(nrow(FST)) FST[, same_ancestry := anc_a == anc_b]

  ## ---- 4. join to the transferability grid: does any of this TRACK R2? ----
  if(momi_has_intermediate("transfer_grid", P)){
    G <- momi_read_intermediate("transfer_grid", P)[definition==MOMI_BPDEF_DEFAULT,
                                                    .(score_id, cohort, R2pct, n_transfer=n)]
    S3 <- merge(S3, G, by=c("score_id","cohort"), all.x=TRUE)
  }
  out  <- momi_write_table(S3[order(score_id, platform, -coverage_pct)], "S3_scoresnp", P)
  out2 <- if(nrow(FST)) momi_write_table(FST[order(-fst_hudson)], "S3b_pairwise_fst", P) else NA_character_

  ## ---- 5. console: the two questions this table exists to answer ----
  cat("\n== coverage and INFO, dosage platform, chosen EUR SBP score ==\n")
  print(S3[platform=="lpwgs_dosage" & score_id==MOMI_EUR$SBP,
           .(cohort, cohort_anc, coverage_pct, mean_info, R2pct)][order(-coverage_pct)])
  if(nrow(FST)){
    cat("\n== pairwise Hudson Fst at EUR score SNPs ==\n"); print(FST[order(-fst_hudson)])
    sa <- FST[same_ancestry==TRUE]
    cat("\nTHE DIAGNOSTIC ROW is AMANHI-Bangladesh vs GAPPS-Bangladesh: same country, same\n",
        "ancestry, so Fst should be ~0. A non-trivial value there means the divergence is\n",
        "genotyping, not populations.\n", sep="")
    if(nrow(sa)) print(sa[, .(cohort_a, cohort_b, fst_hudson, mean_abs_af_diff)])
  }
  if("R2pct" %in% names(S3)){
    d <- S3[platform=="lpwgs_dosage" & is.finite(R2pct) & is.finite(coverage_pct)]
    if(nrow(d) > 3) cat(sprintf("\ncorrelation across cohorts x scores:  coverage vs R2 = %.3f",
                                suppressWarnings(cor(d$coverage_pct, d$R2pct))))
    if(nrow(d) > 3 && any(is.finite(d$mean_info)))
      cat(sprintf(" ;  INFO vs R2 = %.3f",
                  suppressWarnings(cor(d$mean_info, d$R2pct, use="complete.obs"))))
    cat("\nA strong positive correlation would mean measurement, not ancestry, drives the\n",
        "portability gradient. A weak one means the deficit is real -- which is what the\n",
        "recovery and reliability analyses have already implied twice.\n", sep="")
  }

  list(n=nrow(S3),
       key=sprintf("cohorts=%d platforms=%d scores=%d; coverage %.0f-%.0f%%; Fst pairs=%d",
                   uniqueN(S3$cohort), uniqueN(S3$platform), uniqueN(S3$score_id),
                   min(S3$coverage_pct,na.rm=TRUE), max(S3$coverage_pct,na.rm=TRUE), nrow(FST)),
       outputs=c(basename(out), if(!is.na(out2)) basename(out2)))
})
