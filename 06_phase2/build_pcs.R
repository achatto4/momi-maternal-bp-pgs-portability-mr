#!/usr/bin/env Rscript
# ============================================================
# build_pcs.R  [ build step B06b ]  -> intermediates/pcs.rds
# Principal components, two kinds, from one shared variant set.
#
# WHY TWO KINDS. They answer different questions and only one of them needs a merge:
#
#   WITHIN-COHORT PCs  -> for ADJUSTMENT (MR arm B22, and S17/B30).
#     The confounding to control is population structure INSIDE a cohort. PCs computed
#     separately per cohort describe exactly that, and no cross-cohort coordinate system is
#     required. The Burgess MR guidelines list "genomic principal components of ancestry"
#     among the covariates that SHOULD be adjusted for (ref/methods_citations.tsv), and
#     without these our MR arm leaves stratification uncontrolled.
#
#   1000G-PROJECTED PCs -> for DESCRIPTION (SF1/B28, and SF2's PC-distance).
#     Each cohort projected onto axes defined by the 1000 Genomes reference. Anchored to
#     known populations, so a cohort's position is interpretable rather than relative to
#     whatever happens to vary most in our own five cohorts.
#
# WHY NOT PCA THE MERGED FILESET. merged_all mixes GSA and lpWGS, and platform is confounded
# with cohort (AMANHI-B is 31.5% GSA-only, GAPPS-B 0.1%), so its leading component would
# track genotyping platform rather than ancestry and be uninterpretable. Projection avoids
# the merge entirely and each cohort is processed independently.
#
# THE SHARED VARIANT SET. Both PC types are computed on ONE LD-pruned common-variant set
# derived from the 1000G reference, so the within-cohort and projected coordinates rest on
# the same markers and cannot disagree for trivial reasons.
#
# KNOWN ARTEFACT: projected PCs shrink toward the origin relative to the reference (the
# projection is a regression prediction, so it has less variance than the thing it predicts).
# Standard and expected. Do NOT read a cohort sitting "inside" the reference cloud as
# admixture without accounting for it -- SF1 must either correct or state this.
#
#   Rscript build_pcs.R --outroot DIR --g1k PREFIX [--npc 10] [--plink2 plink2]
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

OUTROOT <- momi_arg("--outroot", Sys.getenv("MOMI_OUTROOT"))
G1K     <- momi_arg("--g1k",     Sys.getenv("MOMI_G1K"))
NPC     <- as.integer(momi_arg("--npc", "10"))
PLINK2  <- momi_arg("--plink2", "plink2")
WORK    <- file.path(P$root, "qc", "pcs"); dir.create(WORK, showWarnings=FALSE, recursive=TRUE)

sh <- function(...) system(paste(...), intern=FALSE)

momi_deliverable("B06b_build_pcs", script="06_phase2/build_pcs.R",
                 inputs="genotypes;g1k", P=P, stop_on_error=FALSE, body=function(ctx){
  if(is.null(OUTROOT) || !dir.exists(OUTROOT))
    return(list(skip=TRUE, reason="no --outroot; see docs/paths_reference.md"))
  if(is.null(G1K) || !file.exists(paste0(G1K, ".bed")))
    return(list(skip=TRUE, reason=sprintf("no 1000G reference at %s.bed", G1K)))
  if(Sys.which(PLINK2)=="") return(list(skip=TRUE, reason="plink2 not on PATH"))

  ## ---- 1. LD-pruned common variant set from the reference ----
  prune <- file.path(WORK, "g1k_pruned")
  if(!file.exists(paste0(prune, ".prune.in"))){
    cat("LD-pruning the 1000G reference...\n")
    sh(shQuote(PLINK2), "--bfile", shQuote(G1K),
       "--maf 0.05 --geno 0.02 --hwe 1e-6",
       "--indep-pairwise 200 50 0.2",
       "--out", shQuote(prune), ">", shQuote(paste0(prune, ".log")), "2>&1")
  }
  if(!file.exists(paste0(prune, ".prune.in")))
    return(list(skip=TRUE, reason="LD-pruning produced no .prune.in — check the reference"))
  nprune <- length(readLines(paste0(prune, ".prune.in")))
  cat(sprintf("shared variant set: %d LD-pruned common variants\n", nprune))

  ## ---- 2. reference PCA, saving ALLELE WEIGHTS so cohorts can be projected ----
  ref <- file.path(WORK, "g1k_pca")
  if(!file.exists(paste0(ref, ".eigenvec.allele"))){
    cat("computing reference PCs (with allele weights)...\n")
    sh(shQuote(PLINK2), "--bfile", shQuote(G1K),
       "--extract", shQuote(paste0(prune, ".prune.in")),
       "--freq counts", "--pca", NPC, "allele-wts",
       "--out", shQuote(ref), ">", shQuote(paste0(ref, ".log")), "2>&1")
  }
  if(!file.exists(paste0(ref, ".eigenvec.allele")))
    return(list(skip=TRUE, reason="reference PCA produced no allele weights"))

  ## ---- 2b. TRANSLATE reference variant IDs into the cohort convention ----
  ## The two use incompatible schemes and --extract matched 0 of 279,047 on the first run:
  ##   1000G   "1:15447:A:G"   no chr prefix, alleles appended
  ##   cohorts "chr1:15189"    chr prefix, position only
  ## Translate reference -> cohort form for every file the cohort steps consume: the prune
  ## list, the allele weights, and the allele counts. The reference PCA itself keeps native
  ## IDs. NB chr:pos is not guaranteed unique if multiallelics survived, so dedupe.
  cohort_id <- function(x){
    p <- tstrsplit(x, ":", fixed=TRUE)
    sprintf("chr%s:%s", sub("^chr", "", p[[1]]), p[[2]])
  }
  pin  <- readLines(paste0(prune, ".prune.in"))
  pin2 <- unique(cohort_id(pin))
  prune_c <- file.path(WORK, "prune_cohortIDs.txt")
  writeLines(pin2, prune_c)
  cat(sprintf("translated prune list: %d reference IDs -> %d unique cohort-form IDs\n",
              length(pin), length(pin2)))

  aw <- fread(paste0(ref, ".eigenvec.allele"))
  idc <- if("ID" %in% names(aw)) "ID" else names(aw)[2]
  aw[[idc]] <- cohort_id(aw[[idc]])
  aw <- aw[!duplicated(aw[[idc]])]
  aw_c <- file.path(WORK, "eigenvec_allele_cohortIDs.tsv")
  fwrite(aw, aw_c, sep="\t")

  ac_c <- NULL
  if(file.exists(paste0(ref, ".acount"))){
    ac <- fread(paste0(ref, ".acount"))
    idc2 <- if("ID" %in% names(ac)) "ID" else names(ac)[2]
    ac[[idc2]] <- cohort_id(ac[[idc2]])
    ac <- ac[!duplicated(ac[[idc2]])]
    ac_c <- file.path(WORK, "acount_cohortIDs.tsv")
    fwrite(ac, ac_c, sep="\t")
  }

  ## sanity: do the translated IDs actually exist in a cohort fileset? Fail loudly if not,
  ## rather than running five cohorts and reporting an empty result 3 minutes later.
  probe <- file.path(OUTROOT, MOMI_COH_ALL[1], "lpwgs_dosage_clean.pvar")
  if(file.exists(probe)){
    pv <- fread(cmd=sprintf("grep -v '^##' %s | cut -f3", shQuote(probe)), header=TRUE)
    hit <- sum(pin2 %in% pv[[1]])
    cat(sprintf("ID-translation check against %s: %d of %d prune variants present (%.1f%%)\n",
                MOMI_COH_ALL[1], hit, length(pin2), 100*hit/length(pin2)))
    if(hit < 1000)
      return(list(skip=TRUE, reason=sprintf(
        "ID translation still does not match: only %d of %d prune variants found in %s. Reference IDs look like '%s', cohort IDs like '%s'.",
        hit, length(pin2), basename(probe), pin[1], pv[[1]][1])))
  }

  ## reference sample coordinates + population labels, for the SF1 backdrop
  refev <- fread(paste0(ref, ".eigenvec"))
  setnames(refev, 1:2, c("FID","IID"))
  panel <- paste0(G1K, ".samples.panel")
  if(file.exists(panel)){
    pn <- tryCatch(fread(panel), error=function(e) NULL)
    if(!is.null(pn) && ncol(pn) >= 3){
      setnames(pn, 1:3, c("IID","pop","superpop"))
      refev <- merge(refev, pn[, .(IID=as.character(IID), pop, superpop)],
                     by.x="IID", by.y="IID", all.x=TRUE)
    }
  }
  refev[, cohort := "1000G"]

  ## ---- 3. per cohort: within-cohort PCA, and projection onto the reference axes ----
  wcol <- paste0("PC", seq_len(NPC))          # within-cohort
  pcol <- paste0("PC", seq_len(NPC), "_1kg")  # projected
  rows <- list(); prj <- list()
  for(coh in MOMI_COH_ALL){
    pre <- file.path(OUTROOT, coh, "lpwgs_dosage_clean")
    isp <- file.exists(paste0(pre, ".pgen"))
    if(!isp && !file.exists(paste0(pre, ".bed"))){ cat("skip (no fileset):", coh, "\n"); next }
    gin <- if(isp) "--pfile" else "--bfile"

    ## (a) WITHIN-cohort PCs on the shared variant set
    wout <- file.path(WORK, paste0(coh, "_within"))
    if(!file.exists(paste0(wout, ".eigenvec"))){
      cat("within-cohort PCA:", coh, "\n")
      sh(shQuote(PLINK2), gin, shQuote(pre),
         "--extract", shQuote(prune_c),
         "--pca", NPC, "--out", shQuote(wout),
         ">", shQuote(paste0(wout, ".log")), "2>&1")
    }
    if(file.exists(paste0(wout, ".eigenvec"))){
      w <- fread(paste0(wout, ".eigenvec"))
      setnames(w, 1:2, c("FID","IID"))
      cn <- setdiff(names(w), c("FID","IID"))
      setnames(w, cn[seq_len(min(NPC, length(cn)))], wcol[seq_len(min(NPC, length(cn)))])
      rows[[coh]] <- cbind(data.table(IID=as.character(w$IID), cohort=coh),
                           w[, intersect(wcol, names(w)), with=FALSE])
    }

    ## (b) PROJECTION onto the reference axes, via the allele weights.
    ## plink2 --score with variance-standardize is the documented projection recipe; the
    ## allele-weight file has ID, A1, then one weight column per PC starting at column 5.
    pout <- file.path(WORK, paste0(coh, "_proj"))
    if(!file.exists(paste0(pout, ".sscore"))){
      cat("projecting onto 1000G axes:", coh, "\n")
      sh(shQuote(PLINK2), gin, shQuote(pre),
         "--extract", shQuote(prune_c),
         if(!is.null(ac_c)) paste("--read-freq", shQuote(ac_c)) else "",
         "--score", shQuote(aw_c),
         "2 5 header-read no-mean-imputation variance-standardize",
         "--score-col-nums", sprintf("6-%d", 5+NPC-1),
         "--out", shQuote(pout), ">", shQuote(paste0(pout, ".log")), "2>&1")
    }
    if(file.exists(paste0(pout, ".sscore"))){
      s <- fread(paste0(pout, ".sscore"))
      sc <- grep("SCORE.*AVG|_AVG$", names(s), value=TRUE)
      if(length(sc)){
        k <- min(NPC, length(sc))
        pp <- cbind(data.table(IID=as.character(s[[1]]), cohort=coh),
                    setnames(s[, sc[seq_len(k)], with=FALSE], pcol[seq_len(k)]))
        prj[[coh]] <- pp
      }
    }
  }

  W <- if(length(rows)) rbindlist(rows, fill=TRUE) else data.table()
  J <- if(length(prj))  rbindlist(prj,  fill=TRUE) else data.table()
  if(!nrow(W) && !nrow(J))
    return(list(skip=TRUE, reason="neither within-cohort PCA nor projection produced output"))

  PCS <- if(nrow(W) && nrow(J)) merge(W, J, by=c("IID","cohort"), all=TRUE) else
         if(nrow(W)) W else J
  momi_save_intermediate(PCS, "pcs", P)
  momi_save_intermediate(refev, "pcs_reference", P)

  ## ---- console ----
  cat(sprintf("\npcs.rds: %d mothers, %d cohorts; within-cohort PCs=%d, projected PCs=%d\n",
              nrow(PCS), uniqueN(PCS$cohort),
              sum(wcol %in% names(PCS)), sum(pcol %in% names(PCS))))
  cat("\ncoverage by cohort (should match the dosage sample counts):\n")
  print(PCS[, .(n=.N,
                within = sum(is.finite(get(wcol[1]))),
                projected = if(pcol[1] %in% names(PCS)) sum(is.finite(get(pcol[1]))) else 0L),
            by=cohort][order(cohort)])
  if(pcol[1] %in% names(PCS) && pcol[2] %in% names(PCS)){
    cat("\nprojected PC1/PC2 centroid per cohort (vs 1000G superpopulations in pcs_reference):\n")
    print(PCS[, .(PC1=round(mean(get(pcol[1]), na.rm=TRUE),4),
                  PC2=round(mean(get(pcol[2]), na.rm=TRUE),4)), by=cohort][order(PC1)])
    cat("\nThis is what tests whether AMANHI-Pemba is simply 'AFR' as momi_config asserts.\n",
        "Zanzibar has documented Omani and South Asian admixture; if Pemba sits between the\n",
        "AFR and SAS reference clusters rather than within AFR, the label is doing work the\n",
        "data does not support. NB projected PCs shrink toward the origin -- do not read\n",
        "'closer to the middle' as admixture without accounting for that.\n", sep="")
  }

  list(n=nrow(PCS),
       key=sprintf("mothers=%d cohorts=%d prunedSNPs=%d withinPC=%d projPC=%d",
                   nrow(PCS), uniqueN(PCS$cohort), nprune,
                   sum(wcol %in% names(PCS)), sum(pcol %in% names(PCS))),
       outputs=c("pcs.rds","pcs_reference.rds"))
})
