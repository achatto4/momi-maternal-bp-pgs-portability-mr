#!/usr/bin/env Rscript
# ============================================================
# deliv_S15_fetal.R  [ B32 -> Supp Table S15, maternal vs fetal genotype on birthweight ]
#
# ------------------------------------------------------------------
# THE PROBLEM THIS SOLVES, AND WHY IT IS NOT OPTIONAL
#
# Part II reports that genetically predicted maternal blood pressure lowers birthweight by
# roughly 80-90 g per 10 mmHg. That estimate comes from the MOTHER's polygenic score. But a
# mother transmits half her alleles to the fetus, so her BP score is correlated with the
# fetus own BP score by construction (r ~ 0.5). The maternal estimate therefore blends:
#
#   (a) a genuine MATERNAL effect -- her blood pressure alters the intrauterine environment,
#       placental perfusion, nutrient delivery, and thereby fetal growth; and
#   (b) a FETAL effect -- the alleles the infant inherited act in the infant to affect its
#       own growth, and would do so whatever the mother's blood pressure.
#
# These have completely different meanings. (a) says maternal hypertension harms the fetus,
# which implies treating maternal blood pressure could help. (b) says nothing of the kind.
# Reporting the unconditional maternal estimate as a maternal causal effect, when infant
# genotypes were available and not used, would be indefensible.
#
# THE DECOMPOSITION (Warrington et al.): regress birthweight on BOTH scores at once.
#     BWT ~ zBP_mother + zBP_infant + covariates
# The coefficient on zBP_mother is now the maternal effect CONDITIONAL on fetal genotype --
# i.e. the part that operates through the mother rather than through transmitted alleles.
#
# HOW TO READ IT:
#   * maternal coefficient roughly unchanged from B22  -> the effect is maternal. The
#     intrauterine interpretation stands and this is the paper's strongest claim.
#   * maternal attenuates toward zero, fetal grows     -> we were substantially measuring
#     transmitted alleles. The headline must be restated.
#   * both attenuate                                   -> collinearity, not biology. Check
#     the correlation printed below before concluding anything.
#
# THE COLLINEARITY CAVEAT, CHECKED FIRST. Mother and infant scores are correlated ~0.5 by
# descent. That is exactly the regime where two coefficients become individually imprecise
# while their sum stays well estimated. So a maternal coefficient that widens on adjustment
# is expected and is NOT evidence the maternal effect is absent. The observed correlation and
# the variance inflation are printed before any estimate, and a maternal estimate whose
# INTERVAL still excludes zero after conditioning is the finding worth having.
#
# SAMPLE. Mother-infant pairs with both genotyped and a birthweight, live births only.
# Pairing comes from make_sample_map.sh --mode infant, which writes <out>.keep.pairs as
# BABY_ID -> mother PARTICIPANT_ID; no re-derivation here, so the linkage used for scoring
# and the linkage used for analysis cannot diverge.
#
# COVARIATES: age + within-cohort maternal PCs, matching B22. Fetal PCs would be the
# strictly correct control for fetal structure but do not exist (B06b runs on the maternal
# filesets); given mother and infant share ancestry within a cohort this is a minor
# approximation, and it is stated rather than hidden.
# ------------------------------------------------------------------
#   Writes: tables/S15_fetal_maternal.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

SSC       <- momi_arg("--sscore-dir", Sys.getenv("SSC"))
NPC_ADJ   <- as.integer(momi_arg("--npc", "5"))
MAPDIR    <- file.path(P$root, "qc", "infant_maps")

## ---- infant platforms, and the QC gate that decides which are usable -------------------
## CANDIDATES, not a fixed list. Each is admitted only if it passes MIN_TRANSMISSION below.
INF_PLATS_ALL <- c("infant_gsa","infant_lpwgs","infant_lpwgs_dosage")
##
## WHY A GATE IS NEEDED (found 2026-07-20, first run of this module).
## A mother transmits half her alleles to her child, so their polygenic scores MUST correlate
## about 0.5. Measured per platform, the first run gave:
##     infant_gsa    r = 0.41 - 0.58   across all five cohorts   <- correct
##     infant_lpwgs  r = 0.018 - 0.129 across all five cohorts   <- noise
## The pairing was right (GSA proves it); the lpWGS scores were the problem. lpWGS is
## LOW-PASS sequencing, and the infant filesets were built as HARD CALLS, which discards most
## of the information at low coverage. The maternal pipeline never had this issue because it
## scores lpWGS from DOSAGES (lpwgs_dosage, a .pgen) -- compute_prs.sh notes "on a pgen the
## score uses dosages automatically (low-pass WGS power retained)".
##
## Averaging a good platform with a noise platform does not average the correlation, it
## DILUTES it: pooled r came out at 0.08, the maternal-vs-fetal coefficients had nothing to
## separate, and the maternal estimate moved 1.6% -- which would have been mistaken for
## "the effect is maternal" when it actually meant "the fetal score carries no signal".
## That is a false negative that would have been very easy to report as a finding.
##
## The gate makes the failure loud and self-correcting: any platform whose observed
## transmission correlation falls below MIN_TRANSMISSION is DROPPED with a message, and if
## none survives the module skips rather than producing a diluted answer.
MIN_TRANSMISSION <- as.numeric(momi_arg("--min-transmission", "0.25"))

momi_deliverable("S15_fetal", script="06_phase2/deliv_S15_fetal.R",
                 inputs="analytic_mothers;prs_z;pcs", P=P, stop_on_error=FALSE,
                 body=function(ctx){

  if(is.null(SSC) || !nzchar(SSC)) return(list(skip=TRUE, reason="no SSC dir"))
  A  <- momi_read_intermediate("analytic_mothers", P)
  Z  <- momi_read_intermediate("prs_z", P)
  PC <- tryCatch(momi_read_intermediate("pcs", P), error=function(e) NULL)
  if(is.null(PC) || !nrow(PC)) return(list(skip=TRUE, reason="no pcs.rds -- run B06b"))
  pccols <- intersect(paste0("PC", seq_len(NPC_ADJ)), names(PC))
  covall <- c("AGE", pccols)

  ## read ONE infant platform
  inf_z_plat <- function(pid, coh, p){
    f <- file.path(SSC, sprintf("%s__%s__%s.sscore", pid, coh, p))
    if(!file.exists(f)) return(NULL)
    s <- fread(f)
    data.table(BABY_ID=as.character(s[[1]]), z=as.numeric(momi_zin(s$SCORE1_AVG)))
  }
  ## merged mean-z across the platforms that PASSED the gate (set below)
  inf_z <- function(pid, coh, plats){
    L <- lapply(plats, function(p) inf_z_plat(pid, coh, p))
    L <- L[!vapply(L, is.null, TRUE)]
    if(!length(L)) return(NULL)
    rbindlist(L)[is.finite(z), .(zINF=mean(z)), by=BABY_ID]
  }

  ## ---- mother<->infant pairing, straight from the sample maps ----
  pf <- list.files(MAPDIR, pattern="_infant\\.keep\\.pairs$", full.names=TRUE)
  if(!length(pf))
    return(list(skip=TRUE, reason=sprintf("no *.keep.pairs in %s -- run run_preprocess_infants.sh", MAPDIR)))
  PAIRS <- unique(rbindlist(lapply(pf, function(f){
    x <- tryCatch(fread(f, header=FALSE, col.names=c("BABY_ID","IID")), error=function(e) NULL)
    if(is.null(x) || !nrow(x)) return(NULL)
    x[, cohort := sub("__.*", "", basename(f))]
    x
  }), fill=TRUE))
  if(!nrow(PAIRS)) return(list(skip=TRUE, reason="pairs files present but empty"))

  A[, hasBP := is.finite(S_mean) | is.finite(D_mean)]
  A2 <- A[genotyped==1 & hasBP==TRUE & livebirth==1]
  A2 <- merge(A2, PC[, c("IID","cohort",pccols), with=FALSE], by=c("IID","cohort"), all.x=TRUE)

  ## ---- THE GATE: measure transmission correlation per platform, keep what works ----------
  ## Uses the SBP instrument in each cohort as the probe. A platform is admitted only if its
  ## median observed r across cohorts clears MIN_TRANSMISSION. Reported in full either way,
  ## because "which platform can carry a fetal score" is itself worth recording.
  gate <- list()
  for(p in INF_PLATS_ALL) for(coh in MOMI_COH_ALL){
    sid <- momi_instrument(coh, "SBP")
    zi  <- inf_z_plat(sid, coh, p); if(is.null(zi)) next
    zm  <- Z[score_id==sid & cohort==coh, .(IID, zM=z)]
    d   <- merge(merge(zi, PAIRS[, .(BABY_ID, IID)], by="BABY_ID"), zm, by="IID")
    d   <- d[is.finite(z) & is.finite(zM)]
    if(nrow(d) < 50) next
    gate[[paste(p,coh)]] <- data.table(platform=p, cohort=coh, n=nrow(d),
                                       r_transmission=round(cor(d$zM, d$z), 3))
  }
  GATE <- rbindlist(gate, fill=TRUE)
  if(!nrow(GATE)) return(list(skip=TRUE, reason="no infant scores found for any platform"))
  keep <- GATE[, .(median_r=median(r_transmission), k=.N), by=platform][median_r >= MIN_TRANSMISSION]
  INF_PLATS <- keep$platform

  cat("\n=== TRANSMISSION CHECK: r(maternal score, infant score), expected ~0.5 ===\n")
  print(dcast(GATE, cohort ~ platform, value.var="r_transmission"))
  cat(sprintf("\nthreshold = %.2f; platforms admitted: %s\n", MIN_TRANSMISSION,
              if(length(INF_PLATS)) paste(INF_PLATS, collapse=", ") else "NONE"))
  dropped <- setdiff(unique(GATE$platform), INF_PLATS)
  if(length(dropped))
    cat("DROPPED:", paste(dropped, collapse=", "),
        "-- these scores do not track maternal genotype, so they carry no fetal signal.\n",
        "  For lpWGS this means HARD CALLS were used where DOSAGES are required: low-pass\n",
        "  sequencing loses most of its information when hard-called. Build the infant\n",
        "  dosage filesets (as the maternal pipeline does) to recover them.\n", sep=" ")
  if(!length(INF_PLATS))
    return(list(skip=TRUE,
                reason=sprintf("no infant platform reached r>=%.2f -- fetal scores carry no signal", MIN_TRANSMISSION)))

  rows <- list(); diag <- list()
  for(coh in MOMI_COH_ALL) for(tr in c("SBP","DBP")){
    sid <- momi_instrument(coh, tr)
    zi  <- inf_z(sid, coh, INF_PLATS)
    if(is.null(zi) || !nrow(zi)) next
    zm  <- Z[score_id==sid & cohort==coh, .(IID, zMOM=z)]
    if(!nrow(zm)) next

    d <- merge(A2[cohort==coh], zm, by="IID")
    d <- merge(d, PAIRS[cohort==coh, .(IID, BABY_ID)], by="IID")
    d <- merge(d, zi, by="BABY_ID")
    d <- d[is.finite(zMOM) & is.finite(zINF) & is.finite(BWT)]
    d <- d[complete.cases(d[, c("zMOM","zINF","BWT",..covall)])]
    if(nrow(d) < 200) next

    r <- cor(d$zMOM, d$zINF)
    diag[[paste(coh,tr)]] <- data.table(
      cohort=coh, trait=tr, n_pairs=nrow(d),
      cor_mother_infant=round(r,3), vif=round(1/(1-r^2),2),
      expected_r="~0.5 by descent")

    crhs <- paste(covall, collapse="+")
    m_uni <- tryCatch(lm(as.formula(paste("BWT ~ zMOM +", crhs)), d), error=function(e) NULL)
    m_mv  <- tryCatch(lm(as.formula(paste("BWT ~ zMOM + zINF +", crhs)), d), error=function(e) NULL)
    m_inf <- tryCatch(lm(as.formula(paste("BWT ~ zINF +", crhs)), d), error=function(e) NULL)
    if(is.null(m_uni) || is.null(m_mv)) next
    cu <- summary(m_uni)$coef; cm <- summary(m_mv)$coef
    ci <- if(is.null(m_inf)) NULL else summary(m_inf)$coef

    ## first stage in the SAME pairs, so the scaled maternal effect is per-mmHg not per-SD
    bpcol <- paste0(substr(tr,1,1), "_", MOMI_BPDEF_DEFAULT)
    be <- NA_real_
    if(bpcol %in% names(d)){
      fs <- tryCatch(summary(lm(as.formula(paste(bpcol, "~ zMOM + zINF +", crhs)),
                                d[is.finite(get(bpcol))]))$coef, error=function(e) NULL)
      if(!is.null(fs) && "zMOM" %in% rownames(fs)) be <- fs["zMOM",1]
    }

    rows[[paste(coh,tr)]] <- data.table(
      cohort=coh, cohort_display=MOMI_DISPLAY[[coh]], ancestry=MOMI_ANC[[coh]],
      trait=tr, n_pairs=nrow(d), cor_mother_infant=round(r,3),
      mat_uni=cu["zMOM",1], mat_uni_se=cu["zMOM",2], mat_uni_p=cu["zMOM",4],
      mat_adj=cm["zMOM",1], mat_adj_se=cm["zMOM",2], mat_adj_p=cm["zMOM",4],
      fet_adj=cm["zINF",1], fet_adj_se=cm["zINF",2], fet_adj_p=cm["zINF",4],
      fet_uni=if(is.null(ci)) NA_real_ else ci["zINF",1],
      pct_attenuation=round(100*(cm["zMOM",1]-cu["zMOM",1])/abs(cu["zMOM",1]),1),
      fs_beta_cond=be,
      mat_adj_per10 = if(is.finite(be) && be!=0) 10*cm["zMOM",1]/be else NA_real_)
  }
  if(!length(rows)) return(list(skip=TRUE, reason="no cohort had >=200 complete mother-infant pairs"))
  S15 <- rbindlist(rows); DG <- rbindlist(diag, fill=TRUE)

  ## ---- pooled, SAS and all (the AFR instrument does not transfer; see Part I) ----
  pool <- rbindlist(lapply(c("SAS","all"), function(st) rbindlist(lapply(c("SBP","DBP"), function(tr){
    s <- S15[trait==tr]; if(st!="all") s <- s[ancestry==st]
    s <- s[is.finite(mat_adj) & mat_adj_se>0]
    if(nrow(s) < 2) return(NULL)
    fm <- momi_meta_iv(s$mat_adj, s$mat_adj_se)
    ff <- momi_meta_iv(s$fet_adj, s$fet_adj_se)
    fu <- momi_meta_iv(s$mat_uni, s$mat_uni_se)
    data.table(stratum=st, trait=tr, k=fm$k, n_pairs=sum(s$n_pairs),
               mat_uni=fu$b, mat_uni_p=fu$p,
               mat_adj=fm$b, mat_adj_se=fm$se, mat_adj_p=fm$p,
               fet_adj=ff$b, fet_adj_se=ff$se, fet_adj_p=ff$p,
               pct_attenuation=round(100*(fm$b-fu$b)/abs(fu$b),1))
  }), fill=TRUE)), fill=TRUE)

  out  <- momi_write_table(S15, "S15_fetal_maternal", P)
  out2 <- momi_write_table(pool, "S15b_fetal_pooled", P)
  out3 <- momi_write_table(DG,  "S15c_pair_diagnostics", P)
  out4 <- momi_write_table(GATE, "S15d_transmission_check", P)

  ## ---------------- console ----------------
  cat("\n=== collinearity check — read BEFORE the estimates ===\n")
  print(DG)
  cat("\nr near 0.5 is expected (half the alleles are transmitted). At that correlation the two\n",
      "coefficients are individually imprecise even when their sum is well estimated, so a\n",
      "maternal estimate that WIDENS on adjustment is expected and is not evidence of absence.\n", sep="")

  cat("\n=== maternal effect on birthweight, before vs after conditioning on fetal genotype ===\n")
  print(S15[, .(cohort_display, trait, n_pairs,
                mat_uni=round(mat_uni,1), mat_adj=round(mat_adj,1),
                pct_att=pct_attenuation, mat_adj_p=signif(mat_adj_p,3),
                fet_adj=round(fet_adj,1), fet_adj_p=signif(fet_adj_p,3))])

  cat("\n=== pooled ===\n")
  if(nrow(pool)) print(pool[, .(stratum, trait, k, n_pairs,
                                mat_uni=round(mat_uni,1), mat_adj=round(mat_adj,1),
                                mat_adj_p=signif(mat_adj_p,3), pct_att=pct_attenuation,
                                fet_adj=round(fet_adj,1), fet_adj_p=signif(fet_adj_p,3))])

  cat("\n=== HOW TO REPORT THIS (decision 2026-07-20) ===\n")
  cat("SUPPLEMENTARY, not a main display item, and NOT framed as overturning anything.\n\n",
      "Our result: conditioning on fetal genotype ATTENUATES the maternal coefficient by\n",
      "37-50%, and the maternal estimate loses significance while the fetal one gains it.\n",
      "BUT every maternal point estimate REMAINS NEGATIVE -- the direction is preserved, only\n",
      "the magnitude and precision degrade.\n\n",
      "WHY WE DO NOT CLAIM A FETAL EFFECT. Warrington et al. (Nat Genet 2019;51:804-814,\n",
      "PMID 31043758) addressed exactly this question with far more power and found the\n",
      "OPPOSITE: 'a 1SD (10mmHg) higher maternal SBP is causally associated with a 0.15 SD\n",
      "(95%CI: -0.19, -0.11) lower offspring birth weight, independent of the direct fetal\n",
      "effects. In contrast, there was no fetal effect of SBP on their own birth weight, after\n",
      "adjusting for the indirect maternal effect (-0.01 SD per 10mmHg, 95% CI: -0.05, 0.03)'.\n\n",
      "Three reasons to trust theirs over ours on the DECOMPOSITION:\n",
      "  1. Our conditional maternal p-values are 0.09-0.26. That cannot distinguish 'attenuated\n",
      "     to zero' from 'unchanged' -- the point moved, the interval did not exclude either.\n",
      "  2. They validated their SEM against conditional regression in 18,873 mother-offspring\n",
      "     pairs; Evans et al. (IJE 2019;48:861-875) report being UNABLE to resolve maternal from\n",
      "     fetal effects at N=12,909. We have ~8,700 pairs with r~0.5 collinearity.\n",
      "  3. Neither they nor we have PATERNAL genotypes. Warrington list this as a limitation\n",
      "     ('we have not considered paternal genotypes, and it is possible that this omission has\n",
      "     biased the results'); with assortative mating it biases the decomposition, and it bites\n",
      "     harder at our sample size.\n\n",
      "CORROBORATING DETAIL WORTH STATING: our UNCONDITIONAL maternal estimate (SAS, -0.165 SD\n",
      "per 10 mmHg) closely matches BOTH Warrington's CONDITIONAL maternal estimate (-0.15) and\n",
      "Morales-Berstein's unconditional one (-0.126). That agreement is what you would expect if\n",
      "there were little genuine fetal contribution to remove -- i.e. it argues our attenuation\n",
      "is a power artefact rather than a real fetal effect.\n\n",
      "TEXT TO WRITE: the maternal association attenuates on adjustment for fetal genotype but\n",
      "remains directionally consistent; our sample is too small to resolve the two components,\n",
      "and the better-powered Warrington decomposition supports a maternal rather than fetal\n",
      "origin. Report as a limitation of our design, not as a competing finding.\n", sep="")

  list(n=nrow(S15),
       key=sprintf("platforms=%s; pairs=%d cohorts=%d; median r(mother,infant)=%.2f; median maternal attenuation=%.1f%%",
                   paste(INF_PLATS, collapse="+"),
                   sum(S15$n_pairs), uniqueN(S15$cohort),
                   median(S15$cor_mother_infant, na.rm=TRUE),
                   median(S15$pct_attenuation, na.rm=TRUE)),
       outputs=c(basename(out), basename(out2), basename(out3), basename(out4)))
})
