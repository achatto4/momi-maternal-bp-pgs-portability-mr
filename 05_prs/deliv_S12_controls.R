#!/usr/bin/env Rscript
# ============================================================
# deliv_S12_controls.R  [ B24 -> S12, does the instrument behave like an instrument? ]
#
# ------------------------------------------------------------------
# WHY THIS MODULE LOOKS NOTHING LIKE THE ORIGINAL PLAN
#
# The obvious validation was to reproduce the external panel's control structure
# (ref/external_mr_estimates.tsv): HBW 0.76 as the mirror of LBW 1.33, POSTTERM 0.94 as the
# mirror of PTB 1.12, and STILL/MISC at exactly 1.00 as explicit nulls. Reproducing a set of
# mirror images and a pair of nulls would have been strong evidence the instrument is doing
# what we claim.
#
# That route is CLOSED: those four outcomes were removed from MOMI on 2026-07-20 (too few
# events; MOMI_OUTCOMES_DROPPED in momi_config.R). So we cannot validate on the outcome side
# at all, and pretending otherwise with 20-event cells would be worse than not trying.
#
# THE SUBSTITUTE, WHICH IS ARGUABLY THE BETTER TEST ANYWAY.
#
# An instrument earns its name by satisfying the exclusion restriction: the PRS may affect
# the outcome ONLY through blood pressure. The outcome-side controls test that indirectly.
# The tests below attack it directly, on the exposure side, and they need no extra outcomes:
#
#   POSITIVE CONTROLS -- the PRS MUST be associated with things blood pressure genuinely
#     causes. If a BP score does not predict preeclampsia or chronic hypertension, the
#     instrument is not measuring BP liability and nothing downstream is interpretable.
#     PE has 238 cases across cohorts (thin per-cohort -- pooled only; AMANHI-Bangladesh has
#     7 and is reported but must not be read alone). CHRON_HTN is the cleaner of the two.
#
#   NEGATIVE CONTROLS -- the PRS must NOT be associated with things blood pressure cannot
#     plausibly cause and that instead mark population structure or social position:
#     maternal EDUCATION, WEALTH, GRAVIDITY, AGE. A genotype is fixed at conception, so it
#     cannot be caused by any of these; an association therefore signals CONFOUNDING BY
#     ANCESTRY (the PRS tracking structure that also tracks social position) or selection.
#     This is the standard stratification check and it is the one that would actually catch
#     the failure mode most likely to be lurking here -- five cohorts, three ancestral
#     backgrounds, scores trained in Europeans.
#
#     NB AGE and the PCs are ADJUSTED FOR in the MR (Burgess), so the negative-control test
#     for those is run UNADJUSTED on purpose: we want to know whether the raw association
#     exists, not whether our adjustment absorbed it. BMI is deliberately EXCLUDED from the
#     negative controls: BP and BMI share genuine biology, so an association there is not
#     evidence of stratification and would be misread as a failure.
#
# HOW TO READ A FAILED NEGATIVE CONTROL. An association between the PRS and education within
# a cohort does not automatically invalidate the MR -- it means the within-cohort PC
# adjustment is incomplete. The relevant comparison is the association BEFORE and AFTER PC
# adjustment: if PCs remove it, the adjustment is working and the MR stands; if it survives,
# the exclusion restriction is in doubt and the causal claim must be softened. Both columns
# are therefore reported for every negative control.
# ------------------------------------------------------------------
#   Writes: tables/S12_controls.tsv
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

NPC_ADJ <- as.integer(momi_arg("--npc", "5"))

## association of z with a variable, with and without PC adjustment
assoc <- function(d, yvar, pccols, type=c("bin","lin")){
  type <- match.arg(type)
  keep <- c("z", yvar, "AGE", pccols)
  x <- d[, intersect(keep, names(d)), with=FALSE]
  setnames(x, yvar, "y")
  x <- x[is.finite(z) & is.finite(y)]
  if(nrow(x) < 50) return(NULL)
  fit <- function(rhs, dd){
    f <- as.formula(paste("y ~", rhs))
    m <- tryCatch(if(type=="bin") glm(f, binomial, dd) else lm(f, dd), error=function(e) NULL)
    if(is.null(m)) return(c(NA,NA,NA))
    cc <- summary(m)$coef
    if(!("z" %in% rownames(cc))) return(c(NA,NA,NA))
    c(cc["z",1], cc["z",2], cc["z",4])
  }
  crude <- fit("z", x)
  xa <- x[complete.cases(x)]
  adj <- if(nrow(xa) >= 50 && all(pccols %in% names(xa)))
           fit(paste(c("z", pccols), collapse="+"), xa) else c(NA,NA,NA)
  data.table(n=nrow(x), n_adj=nrow(xa),
             beta_crude=crude[1], se_crude=crude[2], p_crude=crude[3],
             beta_pcadj=adj[1],  se_pcadj=adj[2],  p_pcadj=adj[3])
}

momi_deliverable("S12_controls", script="05_prs/deliv_S12_controls.R",
                 inputs="analytic_mothers;prs_z;pcs", P=P, stop_on_error=FALSE,
                 body=function(ctx){

  A  <- momi_read_intermediate("analytic_mothers", P)
  Z  <- momi_read_intermediate("prs_z", P)
  PC <- tryCatch(momi_read_intermediate("pcs", P), error=function(e) NULL)
  if(is.null(PC) || !nrow(PC)) return(list(skip=TRUE, reason="no pcs.rds -- run B06b first"))
  pccols <- intersect(paste0("PC", seq_len(NPC_ADJ)), names(PC))

  A[, hasBP := is.finite(S_mean) | is.finite(D_mean)]
  A2 <- A[genotyped==1 & hasBP==TRUE]
  A2 <- merge(A2, PC[, c("IID","cohort",pccols), with=FALSE], by=c("IID","cohort"), all.x=TRUE)

  ## control panel: name -> (column, type, role)
  CTRL <- list(
    list(v="PE",        t="bin", role="positive", why="BP raises preeclampsia risk — the PRS MUST predict this"),
    list(v="CHRON_HTN", t="bin", role="positive", why="chronic hypertension is the exposure itself, near-tautological"),
    list(v="EDU",       t="lin", role="negative", why="genotype cannot be caused by education; assoc => structure"),
    list(v="WEALTH",    t="lin", role="negative", why="as education; also flags selection into cohort"),
    list(v="GRAV",      t="lin", role="negative", why="prior pregnancies cannot be caused by BP genotype"),
    list(v="AGE",       t="lin", role="negative", why="maternal age cannot be caused by genotype"))

  rows <- list()
  for(coh in c(MOMI_COH_ALL, "POOLED")) for(tr in c("SBP","DBP")){
    if(coh=="POOLED"){
      ## pooled = stack all cohorts, adjust for cohort. PCs are within-cohort and therefore
      ## not commensurable across cohorts (see B22 header), so the pooled row is reported
      ## CRUDE ONLY and its PC-adjusted columns are left NA rather than computed wrongly.
      dl <- lapply(MOMI_COH_ALL, function(cc){
        sid <- momi_instrument(cc, tr)
        zc  <- Z[score_id==sid & cohort==cc, .(IID, z)]
        merge(A2[cohort==cc], zc, by="IID")
      })
      d <- rbindlist(dl, fill=TRUE)
      pcuse <- character(0)
    } else {
      sid <- momi_instrument(coh, tr)
      zc  <- Z[score_id==sid & cohort==coh, .(IID, z)]
      d   <- merge(A2[cohort==coh], zc, by="IID")
      pcuse <- pccols
    }
    if(!nrow(d)) next
    for(ct in CTRL){
      if(!(ct$v %in% names(d))) next
      r <- assoc(d, ct$v, pcuse, ct$t)
      if(is.null(r)) next
      rows[[paste(coh,tr,ct$v)]] <- cbind(
        data.table(cohort=coh, trait=tr, variable=ct$v, role=ct$role,
                   type=ct$t, rationale=ct$why), r)
    }
  }
  if(!length(rows)) return(list(skip=TRUE, reason="no control cells estimable"))
  S12 <- rbindlist(rows, fill=TRUE)

  ## a negative control "fails" if the crude association is significant; the question then is
  ## whether PC adjustment removes it. Flag both states separately -- they mean different things.
  S12[, flag := fifelse(role=="negative" & is.finite(p_crude) & p_crude < 0.05 &
                        is.finite(p_pcadj) & p_pcadj < 0.05, "FAILS_AFTER_PC",
                 fifelse(role=="negative" & is.finite(p_crude) & p_crude < 0.05,
                         "crude_only_PC_fixes_it",
                 fifelse(role=="positive" & is.finite(p_crude) & p_crude >= 0.05,
                         "POSITIVE_CONTROL_NULL", "")))]
  setorder(S12, role, variable, trait, cohort)
  out <- momi_write_table(S12, "S12_controls", P)

  ## ---------------- console ----------------
  cat("\n=== POSITIVE controls: the PRS must predict these ===\n")
  print(S12[role=="positive", .(cohort, trait, variable, n,
                                beta=round(beta_crude,3), p=signif(p_crude,3), flag)])
  cat("\nA null here on CHRON_HTN in a cohort with adequate cases means the instrument is not\n",
      "capturing BP liability in that cohort, and its MR row should not be interpreted.\n", sep="")

  cat("\n=== NEGATIVE controls: the PRS must NOT predict these ===\n")
  print(S12[role=="negative", .(cohort, trait, variable, n,
                                p_crude=signif(p_crude,3), p_pcadj=signif(p_pcadj,3), flag)])

  nf <- S12[flag=="FAILS_AFTER_PC"]
  cat(sprintf("\nnegative controls surviving PC adjustment: %d of %d\n",
              nrow(nf), nrow(S12[role=="negative"])))
  if(nrow(nf)){
    cat("These are the ones that matter — an association with a variable a genotype cannot\n",
        "cause, which within-cohort PCs do NOT remove, is evidence the exclusion restriction\n",
        "is violated by residual structure. Any MR row for that cohort must be softened.\n", sep="")
    print(nf[, .(cohort, trait, variable, p_crude=signif(p_crude,3), p_pcadj=signif(p_pcadj,3))])
  } else cat("None — consistent with the exclusion restriction holding.\n")

  cat("\nNB expect ~5% of negative-control tests to trip at p<0.05 by chance; with ",
      nrow(S12[role=="negative"]), " tests that is ~",
      round(0.05*nrow(S12[role=="negative"]),1), " cells. Count before concluding.\n", sep="")

  list(n=nrow(S12),
       key=sprintf("controls=%d; neg failing after PC=%d/%d; pos-control nulls=%d",
                   nrow(S12), nrow(nf), nrow(S12[role=="negative"]),
                   nrow(S12[flag=="POSITIVE_CONTROL_NULL"])),
       outputs=basename(out))
})
