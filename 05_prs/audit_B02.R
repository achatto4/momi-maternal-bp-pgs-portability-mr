#!/usr/bin/env Rscript
# ============================================================
# audit_B02.R — independent, READ-ONLY audit of the B02 polygenic-score build.
#
# Purpose: re-derive B02's inputs straight from the .sscore files and check every
# assumption from scratch, WITHOUT modifying build_prs.R. It deliberately uses the same
# momi_zin() the build uses, so it audits the real code path rather than a reimplementation.
#
# The questions, in order:
#   A  file inventory: which score x cohort x platform files exist, and how big
#   B  PLATFORM OVERLAP -- how many mothers are on both platforms, and is it differential
#      by cohort? (a cohort with more dual-platform mothers gets a better-measured
#      instrument for purely technical reasons, which would inflate its R2)
#   C  DENOM spread WITHIN a file -- do some mothers get far fewer variants than others?
#   D  raw SCORE1_AVG sanity, degenerate/zero-variance files
#   E  cross-platform agreement for dual-platform mothers (correlation of the two z's)
#   F  merged-z variance: single-platform vs dual-platform mothers
#   G  predicted vs actual row counts in the frozen intermediates
#   H  do score IIDs match the analytic sample?
#
# Run:
#   Rscript 05_prs/audit_B02.R          # SSC from environment
# Writes: results/current/qc/B02_audit_*.tsv   (console output is the main product)
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R"))
source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)
QC <- file.path(P$root, "qc"); dir.create(QC, showWarnings=FALSE, recursive=TRUE)

SSC <- momi_arg("--sscore-dir", Sys.getenv("SSC"))
if(is.null(SSC) || !nzchar(SSC)) stop("no SSC: export SSC")

hdr <- function(x) cat(sprintf("\n\n========== %s ==========\n", x))
sec <- function(x) cat(sprintf("\n-- %s --\n", x))
sav <- function(dt, nm){ fwrite(dt, file.path(QC, paste0("B02_audit_",nm,".tsv")), sep="\t"); invisible(NULL) }
pr  <- function(dt, n=45) print(head(dt, n), nrows=n, class=FALSE)

cat("audit_B02.R —", format(Sys.time(), "%F %T"), "\ngit:", momi_git_sha(PIPE), "\nSSC:", SSC, "\n")

## ============================================================
hdr("A. FILE INVENTORY — read every score x cohort x platform the build would read")
recs <- list()
for(i in seq_len(nrow(MOMI_PANEL))){
  pid <- MOMI_PANEL$id[i]; tr <- MOMI_PANEL$trait[i]; an <- MOMI_PANEL$anc[i]
  for(coh in MOMI_COH_ALL) for(pl in MOMI_PLATS){
    f <- file.path(SSC, sprintf("%s__%s__%s.sscore", pid, coh, pl))
    if(!file.exists(f)) next
    s <- tryCatch(fread(f), error=function(e) NULL)
    if(is.null(s) || !nrow(s)) next
    recs[[paste(pid,coh,pl)]] <- data.table(
      score_id=pid, trait=tr, anc=an, cohort=coh, platform=pl,
      IID=as.character(s[[1]]),
      raw = if("SCORE1_AVG" %in% names(s)) as.numeric(s$SCORE1_AVG) else NA_real_,
      denom = if("DENOM" %in% names(s)) as.numeric(s$DENOM) else NA_real_)
  }
}
ALL <- rbindlist(recs, fill=TRUE)
if(!nrow(ALL)) stop("no sscore files read — check SSC and naming")
## standardize exactly as build_prs.R does: momi_zin WITHIN cohort x platform
ALL[, z := momi_zin(raw), by=.(score_id, cohort, platform)]

cat(sprintf("files read: %d | rows: %d | scores: %d | cohorts: %d\n",
            length(recs), nrow(ALL), uniqueN(ALL$score_id), uniqueN(ALL$cohort)))
sec("scores found on the USABLE platforms (expect 8: EUR/SAS/EAS/MVP x SBP/DBP)")
sf <- ALL[, .(cohorts=uniqueN(cohort), platforms=uniqueN(platform), rows=.N),
          by=.(score_id, trait, anc)][order(anc, trait)]
pr(sf); sav(sf, "scores_found")
miss <- setdiff(MOMI_PANEL$id, unique(ALL$score_id))
cat(sprintf("\npanel scores with NO usable file (expected: DIVERSE/MULTI/BMI): %s\n",
            paste(miss, collapse=", ")))

sec("samples per cohort x platform")
cp <- ALL[score_id==ALL$score_id[1], .N, by=.(cohort, platform)]
cp <- ALL[, .(files=uniqueN(score_id), samples=.N/uniqueN(score_id)), by=.(cohort,platform)][order(cohort,platform)]
pr(cp); sav(cp, "samples_per_cohort_platform")

## ============================================================
hdr("B. PLATFORM OVERLAP — the differential-precision question")
cat("A mother on BOTH platforms gets a merged z that averages two measurements, which\n",
    "removes noise and RAISES her instrument's correlation with BP. If dual-platform\n",
    "mothers are concentrated in one cohort, that cohort's transferability is inflated for\n",
    "a purely technical reason -- structurally the same trap as the postpartum contamination\n",
    "found in B01. This is the key table.\n", sep="")
one <- ALL[score_id==MOMI_EUR$SBP]          # coverage is identical across scores; use one
ovm <- one[, .(nplat=uniqueN(platform)), by=.(cohort, IID)]
ov  <- ovm[, .(mothers=.N,
               on_one_platform=sum(nplat==1), on_both=sum(nplat==2),
               pct_dual=round(100*mean(nplat==2),1)), by=cohort][order(cohort)]
pr(ov); sav(ov, "platform_overlap")
cat(sprintf("\nTOTAL mothers=%d  dual-platform=%d (%.1f%%)\n",
            nrow(ovm), sum(ovm$nplat==2), 100*mean(ovm$nplat==2)))
cat("\nTable 1 currently claims 'GAPPS = same mothers on both platforms'. Check that against\n",
    "pct_dual for GAPPS-Bangladesh -- the file inventory suggests it cannot exceed ~7%.\n", sep="")

sec("which platform does each cohort actually rely on?")
rel <- one[, .N, by=.(cohort, platform)]
rel <- dcast(rel, cohort ~ platform, value.var="N", fill=0)
pr(rel); sav(rel, "platform_reliance")

## ============================================================
hdr("C. DENOM SPREAD WITHIN A FILE — is every mother's score equally well measured?")
cat("SCORE1_AVG normalises by the variants actually scored, so a mother with half the\n",
    "variants is not wrong, but she IS noisier. Large within-file spread means the\n",
    "instrument's reliability varies between mothers invisibly.\n\n", sep="")
dn <- ALL[is.finite(denom), .(
        n=.N,
        variants_med = round(median(denom)/2),
        variants_min = round(min(denom)/2),
        variants_max = round(max(denom)/2),
        ratio_min_over_med = round(min(denom)/median(denom), 3),
        cv = round(sd(denom)/mean(denom), 4),
        pct_below_90pct_of_median = round(100*mean(denom < 0.9*median(denom)),1)),
      by=.(score_id, cohort, platform)][order(ratio_min_over_med)]
cat("worst 20 files by min/median variant ratio:\n")
pr(dn, 20); sav(dn, "denom_spread")
cat("\nratio_min_over_med near 1 = every mother scored on essentially the same variants.\n",
    "Values well below 1, or a non-trivial pct_below_90pct_of_median, mean differential\n",
    "instrument quality WITHIN a cohort.\n", sep="")

## ============================================================
hdr("D. RAW SCORE SANITY")
rs <- ALL[, .(n=.N, n_finite=sum(is.finite(raw)),
              mean_raw=signif(mean(raw, na.rm=TRUE),4),
              sd_raw=signif(sd(raw, na.rm=TRUE),4),
              sd_zero = as.integer(!is.finite(sd(raw,na.rm=TRUE)) || sd(raw,na.rm=TRUE)==0),
              n_z_finite=sum(is.finite(z))),
          by=.(score_id, cohort, platform)]
bad <- rs[sd_zero==1 | n_z_finite < n]
cat(sprintf("files with zero/undefined SD or dropped z values: %d\n", nrow(bad)))
if(nrow(bad)) pr(bad, 20) else cat("(none — every file standardised cleanly)\n")
sav(rs, "raw_score_sanity")

## ============================================================
hdr("E. CROSS-PLATFORM AGREEMENT (dual-platform mothers only)")
cat("How much do GSA and lpWGS-dosage agree on the SAME mother's score? This decides how\n",
    "much averaging actually buys -- and a LOW correlation would itself be alarming, since\n",
    "both platforms claim to measure the same genome.\n\n", sep="")
W <- dcast(ALL, score_id + trait + anc + cohort + IID ~ platform, value.var="z")
pls <- intersect(MOMI_PLATS, names(W))
if(length(pls) == 2){
  setnames(W, pls, c("z_a","z_b"))
  cc <- W[is.finite(z_a) & is.finite(z_b),
          .(n_dual=.N, r=round(suppressWarnings(cor(z_a, z_b)),3),
            mean_abs_diff=round(mean(abs(z_a-z_b)),3)),
          by=.(score_id, trait, anc, cohort)][order(r)]
  cat(sprintf("(platforms compared: %s vs %s)\n", pls[1], pls[2]))
  cat("\nlowest 20 correlations:\n"); pr(cc, 20)
  cat("\nsummary of r across all score x cohort pairs:\n")
  print(round(summary(cc$r), 3))
  sav(cc, "cross_platform_agreement")
  cat("\nIf r is high (>0.9) averaging removes little noise and the section-B concern is\n",
      "small. If r is moderate (<0.8) the platforms genuinely disagree, averaging matters a\n",
      "lot, and dual-platform mothers have materially better instruments than the rest.\n", sep="")
} else cat("only one platform present — cannot compare.\n")

## ============================================================
hdr("F. MERGED-z VARIANCE: single-platform vs dual-platform mothers")
cat("build_prs.R averages the two z's WITHOUT re-standardising. The mean of two imperfectly\n",
    "correlated unit-variance variables has SD sqrt((1+r)/2) < 1, so dual-platform mothers\n",
    "sit on a compressed scale relative to single-platform mothers in the same cohort.\n",
    "This quantifies the mixture actually present in prs_z.\n\n", sep="")
M <- ALL[, .(z=mean(z), nplat=.N), by=.(score_id, trait, anc, cohort, IID)]
vv <- M[, .(mothers=.N,
            n_single=sum(nplat==1), n_dual=sum(nplat==2),
            sd_all    = round(sd(z), 3),
            sd_single = round(sd(z[nplat==1]), 3),
            sd_dual   = round(sd(z[nplat==2]), 3)),
        by=.(score_id, cohort)][order(cohort, score_id)]
pr(vv, 45); sav(vv, "merged_z_variance")
cat("\nCompare sd_single with sd_dual within a cohort. A gap means the merged instrument is\n",
    "on two different scales inside the same regression.\n", sep="")

## ============================================================
hdr("G. PREDICTED vs ACTUAL ROW COUNTS")
cat(sprintf("predicted prs_z rows (one per score x mother)          : %d\n", nrow(M)))
cat(sprintf("predicted prs_z_platform rows (one per score x sample) : %d\n", nrow(ALL)))
cat(sprintf("difference (= scores x dual-platform mothers)          : %d\n", nrow(ALL)-nrow(M)))
for(nm in c("prs_z","prs_z_platform")){
  if(momi_has_intermediate(nm, P)){
    x <- momi_read_intermediate(nm, P)
    cat(sprintf("  actual %-16s rows=%d scores=%d cohorts=%d\n",
                nm, nrow(x), uniqueN(x$score_id), uniqueN(x$cohort)))
  } else cat(sprintf("  %s not built yet\n", nm))
}

## ============================================================
hdr("H. DO SCORE IIDs MATCH THE ANALYTIC SAMPLE?")
if(momi_has_intermediate("analytic_mothers", P)){
  A <- momi_read_intermediate("analytic_mothers", P)
  sid <- unique(M$IID); aid <- unique(A$IID)
  cat(sprintf("scored IIDs=%d | in analytic_mothers=%d | NOT matched=%d\n",
              length(sid), length(intersect(sid,aid)), length(setdiff(sid,aid))))
  mm <- merge(M[score_id==MOMI_EUR$SBP, .(cohort, IID)],
              A[, .(IID, has_bp=!is.na(S_mean), genotyped)], by="IID", all.x=TRUE)
  cs <- mm[, .(scored=.N, matched=sum(!is.na(genotyped)),
               with_antenatal_BP=sum(has_bp %in% TRUE),
               pct_with_BP=round(100*mean(has_bp %in% TRUE),1)), by=cohort][order(cohort)]
  cat("\nper cohort: of the mothers with a PRS, how many have an ANTENATAL mean BP\n")
  cat("(this is the Part I analytic sample after B01's antenatal restriction):\n")
  pr(cs); sav(cs, "prs_vs_analytic")
} else cat("analytic_mothers not built — run B01 first.\n")

## ============================================================
hdr("I. IS THE CROSS-PLATFORM DISAGREEMENT BIAS OR SCATTER?")
cat("r = 0.69-0.99 tells us the two platforms disagree, but not HOW. A systematic shift\n",
    "(one platform reading higher) is a different problem from symmetric scatter, and only\n",
    "scatter behaves like classical measurement error. Also checks whether dual-platform\n",
    "mothers are a biased subset -- each z is standardised on ALL samples in its file, so a\n",
    "non-random dual subset would not be centred at zero.\n\n", sep="")
if(exists("W") && all(c("z_a","z_b") %in% names(W))){
  bs <- W[is.finite(z_a) & is.finite(z_b),
          .(n=.N,
            mean_gsa=round(mean(z_a),3), mean_dosage=round(mean(z_b),3),
            mean_diff=round(mean(z_a-z_b),3), sd_diff=round(sd(z_a-z_b),3),
            sd_gsa=round(sd(z_a),3), sd_dosage=round(sd(z_b),3)),
          by=.(cohort, score_id)][order(cohort, score_id)]
  pr(bs, 45); sav(bs, "bias_vs_scatter")
  cat("\nmean_diff near 0 with large sd_diff = symmetric scatter (classical measurement\n",
      "error, attenuates but does not bias). mean_gsa/mean_dosage far from 0 = the\n",
      "dual-platform mothers are NOT a random subset of their cohort.\n", sep="")
}

## ============================================================
hdr("J. DOES AGREEMENT TRACK VARIANT RECOVERY?")
cat("Hypothesis: scores built from more variants are more robust to per-variant imputation\n",
    "error, so agreement should rise with the number of variants recovered. If it does, the\n",
    "EAS/MVP scores reading ~0 everywhere may be partly a MEASUREMENT failure rather than a\n",
    "portability failure.\n\n", sep="")
if(exists("cc")){
  nv <- ALL[is.finite(denom), .(n_variants=as.integer(round(median(denom)/2))),
            by=.(score_id, cohort, platform)]
  nvm <- nv[, .(n_variants=as.integer(mean(n_variants))), by=.(score_id, cohort)]
  jj <- merge(cc, nvm, by=c("score_id","cohort"))
  jj <- merge(jj, unique(MOMI_PANEL[, c("id","anc")]), by.x="score_id", by.y="id", all.x=TRUE)
  setorder(jj, n_variants)
  cat("agreement vs variants recovered (sorted by variants):\n")
  pr(jj[, .(score_id, anc=anc.x, cohort, n_dual, n_variants, r)], 45)
  cat(sprintf("\ncorrelation between r and n_variants across the %d score x cohort pairs: %.3f\n",
              nrow(jj), suppressWarnings(cor(jj$r, jj$n_variants, use="complete.obs"))))
  cat("\nmean r by score ancestry:\n")
  print(jj[, .(pairs=.N, mean_r=round(mean(r),3), mean_variants=as.integer(mean(n_variants))),
           by=.(anc=anc.x)][order(mean_r)])
  cat("\nmean r by COHORT ancestry (is imputation worse in African samples?):\n")
  jj[, coh_anc := MOMI_ANC[cohort]]
  print(jj[, .(pairs=.N, mean_r=round(mean(r),3)), by=coh_anc][order(mean_r)])
  sav(jj, "agreement_vs_variants")
}

## ============================================================
hdr("K. CROSS-PLATFORM r AS A RELIABILITY ESTIMATE — IMPLIED ATTENUATION")
cat("Two platforms measuring the same genome are PARALLEL MEASURES, so their correlation\n",
    "estimates the reliability of a single-platform PRS -- exactly as the ICC estimates the\n",
    "reliability of a single BP reading. The paper already disattenuates transferability for\n",
    "BP measurement error; it does NOT disattenuate for PRS measurement error, and that\n",
    "error is worst in the African cohorts.\n",
    "  single-platform reliability  ~ r\n",
    "  dual-platform reliability    ~ 2r/(1+r)   (Spearman-Brown)\n",
    "  R2_true                      ~ R2_observed / reliability\n\n", sep="")
if(exists("cc")){
  kk <- cc[, .(score_id, cohort, n_dual, r)]
  kk[, `:=`(rel_single = round(r,3),
            rel_dual   = round(2*r/(1+r),3),
            R2_inflation_single = round(1/r,3))]
  kk[, chosen := score_id == vapply(cohort, function(c) momi_instrument(c,"SBP"), "")]
  cat("CHOSEN SBP instrument per cohort -- the rows that actually matter for Table 2:\n")
  pr(kk[chosen==TRUE][order(cohort)], 10)
  cat("\nA cohort with rel_single = 0.85 has its observed R2 attenuated by ~15%. If that\n",
      "reliability is systematically lower in AFR cohorts, part of the AFR transferability\n",
      "deficit is genotyping error, not ancestry mismatch.\n", sep="")
  sav(kk, "reliability_from_platforms")
}

## ============================================================
hdr("L. DECISIVE TEST — TRANSFERABILITY BY PLATFORM IN THE *SAME* WOMEN")
cat("Same mothers, same antenatal BP, three instruments: GSA-only z, dosage-only z, and the\n",
    "merged average. Any difference is pure measurement, since genetics and phenotype are\n",
    "held fixed. If merged > both singles, averaging is removing real noise and single-\n",
    "platform cohorts are being penalised.\n\n", sep="")
source(file.path(PIPE,"lib/momi_estimators.R"))
if(momi_has_intermediate("analytic_mothers", P) && exists("W")){
  A2 <- momi_read_intermediate("analytic_mothers", P)
  D  <- merge(W[is.finite(z_a) & is.finite(z_b)],
              A2[, .(IID, AGE, S_mean, D_mean)], by="IID")
  rows <- list()
  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    sid <- momi_instrument(coh, tr)
    bpc <- if(tr=="SBP") "S_mean" else "D_mean"
    d <- D[cohort==coh & score_id==sid & is.finite(get(bpc)) & is.finite(AGE)]
    if(nrow(d) < 60) next
    ta <- momi_transfer(d$z_a,               d[[bpc]], cov=d[, .(AGE)])
    tb <- momi_transfer(d$z_b,               d[[bpc]], cov=d[, .(AGE)])
    tm <- momi_transfer((d$z_a + d$z_b)/2,   d[[bpc]], cov=d[, .(AGE)])
    rows[[paste(coh,tr)]] <- data.table(
      cohort=coh, trait=tr, score=sid, n_dual=nrow(d),
      R2_gsa    = round(100*ta$incR2, 3),
      R2_dosage = round(100*tb$incR2, 3),
      R2_merged = round(100*tm$incR2, 3),
      F_gsa=round(ta$F,1), F_dosage=round(tb$F,1), F_merged=round(tm$F,1))
  }
  if(length(rows)){
    L <- rbindlist(rows)
    pr(L, 20); sav(L, "transfer_by_platform_same_women")
    cat("\nOnly cohorts with >=60 dual-platform mothers appear (Zambia 448, GAPPS-B 249,\n",
        "Pemba 111; AMANHI-B has 8 and Pakistan 0, so neither can be tested).\n", sep="")
  } else cat("no cohort has enough dual-platform mothers with BP to test.\n")
} else cat("analytic_mothers not available.\n")

## ============================================================
hdr("M. PLATFORM-OF-ORIGIN — DOES THE GSA/DOSAGE MIX CONFOUND F3?")
cat("Section L could only use dual-platform mothers, and AMANHI-Bangladesh has 8 of them.\n",
    "But the platform MIX is wildly asymmetric in exactly the two cohorts F3 compares:\n",
    "  AMANHI-Bangladesh   831 GSA-only (31.5%),  1795 dosage-only,   8 dual\n",
    "  GAPPS-Bangladesh      5 GSA-only ( 0.1%),  3318 dosage-only, 249 dual\n",
    "If GSA yields a worse instrument, a third of AMANHI-Sylhet carries a degraded PRS while\n",
    "essentially none of PreSSMat does -- and part of the 3.45%% vs 6.26%% gap that F3\n",
    "attributes to measurement TIMING would actually be genotyping PLATFORM.\n",
    "Caveat throughout: which mothers went onto which platform was not randomised, so these\n",
    "comparisons are confounded by selection. M1 shows how badly.\n", sep="")

if(momi_has_intermediate("analytic_mothers", P)){
  A3 <- momi_read_intermediate("analytic_mothers", P)
  PM <- ALL[, .(z = mean(z), nplat = .N,
                plat = paste(sort(unique(platform)), collapse="+")),
            by=.(score_id, trait, anc, cohort, IID)]
  PM[, plat_group := fifelse(nplat == 2L, "dual", plat)]
  D2 <- merge(PM, A3[, .(IID, AGE, BMI, S_mean, D_mean)], by="IID")

  sec("M1. are the platform subgroups comparable? (selection check)")
  ch <- D2[score_id == MOMI_EUR$SBP,
           .(mothers=.N,
             pct_with_BP = round(100*mean(is.finite(S_mean)),1),
             age  = round(mean(AGE, na.rm=TRUE),1),
             bmi  = round(mean(BMI, na.rm=TRUE),1),
             mean_SBP = round(mean(S_mean, na.rm=TRUE),1),
             sd_SBP   = round(sd(S_mean,   na.rm=TRUE),1)),
           by=.(cohort, plat_group)][order(cohort, plat_group)]
  pr(ch, 20); sav(ch, "platform_subgroup_characteristics")
  cat("\nA difference in sd_SBP between platform groups matters directly: R2 depends on the\n",
      "outcome's variance, so unequal BP spread alone would produce unequal R2.\n", sep="")

  sec("M2. transferability by platform-of-origin, chosen instrument per cohort")
  rows <- list()
  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    sid <- momi_instrument(coh, tr); bpc <- if(tr=="SBP") "S_mean" else "D_mean"
    for(pg in c("gsa","lpwgs_dosage","dual")){
      d <- D2[cohort==coh & score_id==sid & plat_group==pg &
              is.finite(get(bpc)) & is.finite(AGE)]
      if(nrow(d) < 60) next
      tt <- momi_transfer(d$z, d[[bpc]], cov=d[, .(AGE)])
      rows[[paste(coh,tr,pg)]] <- data.table(
        cohort=coh, trait=tr, score=sid, platform=pg, n=tt$n,
        R2pct=round(100*tt$incR2,3), beta=round(tt$beta,3), F=round(tt$F,1))
    }
  }
  if(length(rows)){ M2 <- rbindlist(rows); pr(M2, 40); sav(M2, "transfer_by_platform_origin") }
  cat("\nCompare gsa vs lpwgs_dosage WITHIN a cohort. A large, consistent gap means the\n",
      "platform mix is a real confounder for every cross-cohort comparison in Part I.\n", sep="")

  sec("M3. THE F3 RE-TEST — AMANHI-B vs GAPPS-B on dosage-only mothers")
  cat("F3 uses one score (SAS SBP) across both Bangladeshi cohorts to isolate timing from\n",
      "genetics. Restricting BOTH cohorts to dosage-only mothers removes the platform mix as\n",
      "well. If the gap survives, F3 is safe; if it shrinks materially, F3's explanation and\n",
      "its headline numbers need revising.\n\n", sep="")
  sid <- MOMI_SAS$SBP
  f3 <- list()
  for(pg in c("ALL mothers","dosage-only")){
    for(coh in c("AMANHI-Bangladesh","GAPPS-Bangladesh")){
      d <- D2[cohort==coh & score_id==sid & is.finite(S_mean) & is.finite(AGE)]
      if(pg=="dosage-only") d <- d[plat_group=="lpwgs_dosage"]
      if(nrow(d) < 60) next
      tt <- momi_transfer(d$z, d$S_mean, cov=d[, .(AGE)])
      f3[[paste(pg,coh)]] <- data.table(restriction=pg, cohort=coh, n=tt$n,
                                        R2pct=round(100*tt$incR2,3), F=round(tt$F,1))
    }
  }
  if(length(f3)){
    F3 <- rbindlist(f3)
    pr(F3, 10); sav(F3, "f3_retest_dosage_only")
    for(pg in unique(F3$restriction)){
      x <- F3[restriction==pg]
      if(nrow(x)==2) cat(sprintf("  %-12s gap: GAPPS-B %.2f%% vs AMANHI-B %.2f%%  ->  ratio %.2fx\n",
          pg, x[cohort=="GAPPS-Bangladesh"]$R2pct, x[cohort=="AMANHI-Bangladesh"]$R2pct,
          x[cohort=="GAPPS-Bangladesh"]$R2pct / max(x[cohort=="AMANHI-Bangladesh"]$R2pct, 1e-9)))
    }
    cat("\nIf the ratio is materially smaller under 'dosage-only', platform was inflating it.\n")
  }
} else cat("analytic_mothers not available.\n")

hdr("AUDIT COMPLETE")
cat("Per-section TSVs written to:", QC, "\n")
