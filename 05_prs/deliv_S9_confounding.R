#!/usr/bin/env Rscript
# ============================================================
# deliv_S9_confounding.R  [ B18 -> Supp Table S9 ]
# Sequential confounder build-up for the OBSERVATIONAL BP -> outcome association.
# First deliverable of Part II, and the first arm of T4's triangulation.
#
# THE CLAIM IT TESTS. The paper's motivation for doing MR at all is that the adjusted
# observational estimate for LBW and SGA is implausible: higher maternal BP appears
# PROTECTIVE against low birth weight once BMI is adjusted for. S9 is where that sign flip
# is either demonstrated or is not. If it does not appear, T4's framing needs rethinking
# before B19-B21 are built on top of it.
#
# THE TRAP THIS MODULE IS BUILT TO AVOID. Adding a confounder changes the model AND the
# sample, because complete-case analysis drops anyone missing that confounder. BMI is
# computable for 21,118 of 21,685 mothers and BIRTH_WEIGHT is 12-32% missing depending on
# cohort, so an M1-vs-M2 comparison on each model's own complete cases is NOT like-for-like:
# an apparent "flip on adjustment" could be a change of sample. Every model is therefore fit
# TWICE --
#   sample="common"  all four models on the rows complete for EVERY covariate (the honest
#                    comparison: only the adjustment changes)
#   sample="own"     each model on its own complete cases (the conventional presentation)
# If the flip appears in "own" but not "common", it is sample composition, not confounding.
#
# Writes: tables/S9_confounding.tsv
#   Rscript deliv_S9_confounding.R
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

## Sequential build-up. Follows the locked confounder decision (age + gravidity + education
## + BMI, site pooled only): PARITY is 18-34% missing at four of five sites so gravidity is
## used, and SMOK_FREQ is ~99% missing so smoking is absent.
## NB audit_B01 found SNIFF_TOBA and PASSIVE_SMOK are near-complete in five of six cohorts
## and have never been tested. Adding a smoking step would need B01 to extract them first,
## so it is deliberately NOT done here -- flagged in the console output instead.
MODELS <- list(
  M0_unadjusted   = character(0),
  M1_age_grav     = c("AGE","GRAV"),
  M2_plus_BMI     = c("AGE","GRAV","BMI"),
  M3_plus_EDU     = c("AGE","GRAV","BMI","EDU"))

OUTCOMES <- list(PTB=list(col="PTB", type="bin", livebirth=FALSE),
                 LBW=list(col="LBW", type="bin", livebirth=TRUE),
                 SGA=list(col="SGA", type="bin", livebirth=TRUE),
                 BWT=list(col="BWT", type="lin", livebirth=TRUE))

momi_deliverable("S9_confounding", script="05_prs/deliv_S9_confounding.R",
                 inputs="analytic_mothers", P=P, body=function(ctx){
  A <- momi_read_intermediate("analytic_mothers", P)

  ## exposures: the primary definition from config, plus `mean` for continuity with every
  ## previously locked table. Never hardcode a definition (tiering changed 2026-07-19).
  defs <- unique(c(MOMI_BPDEF_DEFAULT, "mean"))
  EXPO <- CJ(trait=c("SBP","DBP"), def=defs, sorted=FALSE)
  EXPO[, col := paste0(substr(trait,1,1), "_", def)]
  EXPO <- EXPO[col %in% names(A)]

  scopes <- c("pooled", MOMI_COH_ALL)
  rows <- list()

  for(s in seq_len(nrow(EXPO))){
    bpcol <- EXPO$col[s]
    for(onm in names(OUTCOMES)){
      O <- OUTCOMES[[onm]]
      if(!O$col %in% names(A)) next
      for(sc in scopes){
        d0 <- if(sc=="pooled") copy(A) else A[cohort==sc]
        ## Part I/II analytic sample: genotyped is NOT required for an observational model,
        ## but the live-birth restriction is, for birthweight-derived outcomes.
        if(O$livebirth) d0 <- d0[livebirth==1]
        conf_full <- MODELS[["M3_plus_EDU"]]
        ## site enters ONLY pooled models -- it is constant within a cohort
        add_site <- (sc=="pooled")

        ## the common-sample mask: rows usable by EVERY model in the sequence
        need <- c(bpcol, O$col, conf_full)
        dcom <- d0[complete.cases(d0[, ..need])]

        for(mnm in names(MODELS)){
          conf <- MODELS[[mnm]]
          cf   <- c(conf, if(add_site && length(conf)) "cohort" else NULL)
          for(samp in c("common","own")){
            dd <- if(samp=="common") dcom else d0
            r <- tryCatch(momi_obs(dd, bpcol, O$col, type=O$type, conf=cf, per=10),
                          error=function(e) NULL)
            if(is.null(r)) next
            rows[[paste(bpcol,onm,sc,mnm,samp)]] <- data.table(
              exposure=bpcol, trait=EXPO$trait[s], definition=EXPO$def[s],
              outcome=onm, scope=sc, model=mnm, sample=samp,
              n=r$n, est=round(r$est,4), se=round(r$se,4),
              OR=if(O$type=="bin") round(r$OR,3) else NA_real_,
              lo=if(O$type=="bin") round(r$OR_lo,3) else round(r$lo,3),
              hi=if(O$type=="bin") round(r$OR_hi,3) else round(r$hi,3),
              p=signif(r$p,3),
              ## N1 (2026-07-26): direction must be outcome-aware. For birth weight a POSITIVE
              ## effect (more grams) is protective, not harmful; for the adverse binary outcomes
              ## a positive log-odds is harmful.
              direction=fifelse(onm=="BWT",
                                fifelse(r$est > 0, "protective", "harmful"),
                                fifelse(r$est > 0, "harmful", "protective")))
          }
        }
      }
    }
  }
  if(!length(rows)) stop("S9: no models fit — check outcome/exposure columns")
  S9 <- rbindlist(rows, fill=TRUE)
  out <- momi_write_table(S9[order(trait, definition, outcome, scope, model, sample)],
                          "S9_confounding", P)

  ## ---- console: the sign flip, and whether it survives a fixed sample ----
  pd <- MOMI_BPDEF_DEFAULT
  cat(sprintf("\n== POOLED, %s definition, SBP -- OR per 10 mmHg (COMMON sample) ==\n", pd))
  w <- dcast(S9[scope=="pooled" & trait=="SBP" & definition==pd & sample=="common" &
                outcome %in% c("PTB","LBW","SGA")],
             outcome + n ~ model, value.var="OR")
  print(w)

  cat("\n== does the flip depend on the SAMPLE or the ADJUSTMENT? ==\n")
  fl <- S9[scope=="pooled" & trait=="SBP" & definition==pd &
           outcome %in% c("LBW","SGA") & model %in% c("M1_age_grav","M2_plus_BMI"),
           .(outcome, model, sample, n, OR, direction)]
  print(dcast(fl, outcome + sample ~ model, value.var=c("n","OR")))
  cat("\nIf LBW/SGA turn protective between M1 and M2 in the COMMON sample, the flip is\n",
      "confounding by BMI and T4's motivation holds. If it only appears in the OWN sample,\n",
      "it is driven by who drops out when BMI is required, and the framing needs rethinking.\n",
      sep="")

  cat("\n== sample loss from requiring every covariate (pooled) ==\n")
  print(S9[scope=="pooled" & trait=="SBP" & definition==pd & model=="M0_unadjusted",
           .(outcome, n_own=n[sample=="own"], n_common=n[sample=="common"],
             pct_lost=round(100*(n[sample=="own"]-n[sample=="common"])/n[sample=="own"],1)),
           by=outcome][, -1])

  cat("\nNOTE smoking is absent by decision (SMOK_FREQ ~99% missing), but audit_B01 found\n",
      "SNIFF_TOBA and PASSIVE_SMOK near-complete in five of six cohorts and never tested.\n",
      "Adding a smoking step needs B01 to extract them first -- not done here.\n", sep="")

  key_lbw <- S9[scope=="pooled" & trait=="SBP" & definition==pd & outcome=="LBW" &
                sample=="common" & model=="M2_plus_BMI"]
  list(n=nrow(S9),
       key=sprintf("rows=%d; pooled %s SBP->LBW M2 OR=%s (%s); flip_in_common_sample=%s",
                   nrow(S9), pd,
                   if(nrow(key_lbw)) as.character(key_lbw$OR) else "NA",
                   if(nrow(key_lbw)) key_lbw$direction else "NA",
                   as.character(nrow(key_lbw) && key_lbw$OR < 1)),
       outputs=basename(out))
})
