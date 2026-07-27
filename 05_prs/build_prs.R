#!/usr/bin/env Rscript
# ============================================================
# build_prs.R  [ build step B02 ]
# Canonical polygenic-score intermediates. For every score in MOMI_PANEL and every
# genotyped cohort, standardize the score WITHIN cohort x platform (z), then form the
# de-duplicated per-mother instrument = mean z across platforms (merged). Both the
# merged and per-platform tables are frozen so transferability (merged) and the
# platform supplement S5 (per-platform) read identical numbers.
#
# Writes (results/current/intermediates/):
#   prs_z.rds           long: IID, cohort, score_id, trait, anc, z   (merged avg-z)
#   prs_z_platform.rds  long: IID, cohort, score_id, trait, anc, platform, z
#
#   Rscript build_prs.R --sscore-dir DIR [--pipe PIPE]
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R"))
source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)

# SSC from CLI flag, else the environment (driver exports it).
SSC <- momi_arg("--sscore-dir", Sys.getenv("SSC"))
if(is.null(SSC) || !nzchar(SSC)) stop("no SSC: pass --sscore-dir or export SSC")

momi_deliverable("B02_build_prs", script="05_prs/build_prs.R",
                 inputs=paste0("sscore=",basename(SSC)), P=P, body=function(ctx){

  merged <- list(); byplat <- list()
  for(i in seq_len(nrow(MOMI_PANEL))){
    pid <- MOMI_PANEL$id[i]; tr <- MOMI_PANEL$trait[i]; an <- MOMI_PANEL$anc[i]
    for(coh in MOMI_COH_ALL){
      # per-platform standardized z
      for(pl in MOMI_PLATS){
        f <- file.path(SSC, sprintf("%s__%s__%s.sscore", pid, coh, pl))
        if(!file.exists(f)) next
        s <- fread(f)
        z <- momi_zin(s$SCORE1_AVG)
        ## Variants actually used. DENOM is the allele-observation count; audit_B02 (2026-07-19)
        ## confirmed it is CONSTANT within every file (cv = 0 in all 80), so the median is exact.
        ## Carried through because recovery is strongly DIFFERENTIAL by ancestry -- for the EUR
        ## SBP score, 80.7% of published variants in AMANHI-Pakistan vs 67.6% in GAPPS-Zambia --
        ## which attenuates African R2 independently of any ancestry-mismatch effect. S1 needs
        ## this per cohort, and re-reading 80 sscore files later to recover it would be silly.
        nv <- if("DENOM" %in% names(s))
                as.integer(round(median(as.numeric(s$DENOM), na.rm=TRUE)/2)) else NA_integer_
        d <- data.table(IID=as.character(s[[1]]), z=as.numeric(z))[is.finite(z)]
        if(nrow(d)) byplat[[paste(pid,coh,pl)]] <-
          data.table(IID=d$IID, cohort=coh, score_id=pid, trait=tr, anc=an, platform=pl,
                     z=d$z, n_variants=nv)
      }
      # merged de-duplicated avg-z (config helper: mean z across platforms per IID)
      zt <- momi_getz(pid, coh, SSC)
      if(!is.null(zt) && nrow(zt))
        merged[[paste(pid,coh)]] <-
          data.table(IID=zt$IID, cohort=coh, score_id=pid, trait=tr, anc=an, z=zt$z)
    }
  }
  prs_z    <- rbindlist(merged, fill=TRUE)
  prs_zpl  <- rbindlist(byplat, fill=TRUE)
  if(!nrow(prs_z)) stop("no PRS merged — check --sscore-dir and score/cohort naming")

  ## Attach platform count and variant count to the MERGED instrument. n_plat matters
  ## because a dual-platform mother's z averages two measurements and so carries less
  ## measurement error than a single-platform mother's -- and coverage is highly
  ## differential: 30.9% dual in GAPPS-Zambia vs 0.0% in AMANHI-Pakistan (audit_B02).
  ## Downstream can now condition on, or adjust for, instrument precision.
  np <- prs_zpl[, .(n_plat = uniqueN(platform),
                    n_variants = as.integer(round(mean(n_variants, na.rm=TRUE)))),
                by=.(score_id, cohort, IID)]
  prs_z <- merge(prs_z, np, by=c("score_id","cohort","IID"), all.x=TRUE)
  momi_save_intermediate(prs_z,   "prs_z", P)
  momi_save_intermediate(prs_zpl, "prs_z_platform", P)

  list(n=nrow(prs_z),
       key=sprintf("scores=%d cohorts=%d rows_merged=%d rows_platform=%d",
                   uniqueN(prs_z$score_id), uniqueN(prs_z$cohort), nrow(prs_z), nrow(prs_zpl)),
       outputs=c("prs_z.rds","prs_z_platform.rds"))
})
