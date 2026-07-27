#!/usr/bin/env Rscript
# ============================================================
# build_transfer_grid.R  [ build step B03 ]
# The full transferability grid: for every score x cohort x BP-definition, the
# incremental R2 (adj age), first-stage beta/SE/F, N, reliability ICC, and the
# disattenuated R2. This single grid is the source for:
#   Table 2 (chosen instrument, mean def)  = a filtered view
#   Figure 2 (score x cohort heatmap)      = a pivot of one def
#   S4 (full grid)                         = the whole thing
#   S7 (GA-residual)                       = def=="resid" rows
#   S8 (<20 / >=20 wk)                     = def in {lt20, ge20}
# Reads intermediates prs_z, analytic_mothers, bp_icc (from B01/B02).
#
# ------------------------------------------------------------------
# REVISION 2026-07-19. Two changes:
#
# 1. DEFINITION COLUMNS ARE NAMED, NOT PATTERN-MATCHED. The previous version found its BP
#    columns with grep("^S_", names(A)). When B01 gained postnatal companion columns on
#    2026-07-19 that regex also caught S_mean_pn, S_first_pn, S_median_pn and -- worst --
#    S_n_pn, which is a COUNT OF READINGS, not a blood pressure. Unfixed it would have
#    silently grown the grid from 440 cells to 600 and regressed the BP-PRS on a visit
#    count. Naming the columns means no future addition to analytic_mothers can change the
#    grid's shape by accident.
#
# 2. A PARALLEL POSTNATAL GRID is written to transfer_grid_pn. The antenatal grid's
#    contract is UNCHANGED (11 definitions, 440 cells), so T2/F2/S4/S7/S8 need no edits.
#    Motivation: the BP scores were trained on NON-PREGNANT adults, and removing postpartum
#    readings LOWERED AMANHI-Bangladesh's R2 from 3.45% to 2.82% while leaving GAPPS-B at
#    6.26% -- i.e. those readings were helping, which is what you would expect if the score
#    predicts the non-pregnant state better than the pregnant one. This grid is the data
#    needed to test that directly.
# ------------------------------------------------------------------
#
# Writes: intermediates/transfer_grid.rds, intermediates/transfer_grid_pn.rds
#         tables/S4_transfer_grid.tsv, tables/QC_transfer_grid_pn.tsv
#   Rscript build_transfer_grid.R [--pipe PIPE]
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R"))
source(file.path(PIPE,"lib/momi_config.R"))
source(file.path(PIPE,"lib/momi_estimators.R"))
P <- momi_paths(PIPE)

## ---- definition columns (see header note 1) ----
## Taken from MOMI_BPDEFS in config, which is the single source of truth for the tiering
## decided 2026-07-19: PRIMARY = resid, ge20; AUXILIARY = mean, lt20, tri1/2/3, last.
## early/median/mn2 are still computed by B01 but deliberately excluded from the grid
## (MOMI_BPDEFS_UNUSED) -- `early` is retained in analytic_mothers for trajectory work.
## Label == column suffix for every retained definition, so no mapping is needed.
ANTE_DEFS <- MOMI_BPDEFS
POST_DEFS <- c("mean_pn","first_pn","median_pn")   # NB: n_pn deliberately excluded (a count)

momi_deliverable("B03_build_transfer_grid", script="05_prs/build_transfer_grid.R",
                 inputs="prs_z;analytic_mothers;bp_icc", P=P, body=function(ctx){

  A   <- momi_read_intermediate("analytic_mothers", P)
  Z   <- momi_read_intermediate("prs_z", P)
  ICC <- momi_read_intermediate("bp_icc", P)

  defcols <- function(prefix, defs){
    want <- paste0(prefix, "_", defs)
    miss <- setdiff(want, names(A))
    if(length(miss)) message("B03: definition columns absent from analytic_mothers: ",
                             paste(miss, collapse=", "))
    want[want %in% names(A)]
  }
  relab <- function(col, prefix){
    d <- sub(paste0("^",prefix,"_"), "", col)
    d <- sub("_pn$", "", d)
    if(d=="first") "early" else d
  }

  ## ---- one routine, used for both the antenatal and the postnatal grid ----
  ## use_icc=FALSE for the postnatal grid: bp_icc is computed on ANTENATAL readings only
  ## (B01, 2026-07-19), so it is not a valid reliability for postpartum BP. Reported as NA
  ## rather than silently disattenuating with the wrong number.
  make_grid <- function(defs, use_icc){
    rows <- list()
    combos <- unique(Z[, .(score_id, trait, anc, cohort)])
    for(j in seq_len(nrow(combos))){
      sid <- combos$score_id[j]; tr <- combos$trait[j]
      an  <- combos$anc[j];      coh <- combos$cohort[j]
      if(tr=="BMI") next                       # BMI score belongs to MVMR (S14), not transfer
      pre  <- substr(tr,1,1)
      cols <- defcols(pre, defs)
      if(!length(cols)) next
      zc <- Z[score_id==sid & cohort==coh, .(IID, z)]
      dm <- merge(zc, A[cohort==coh, c("IID","AGE", cols), with=FALSE], by="IID")
      ## reliability inputs. bp_icc carries the single-reading ICC plus the counts needed to
      ## get k = mean antenatal readings per mother, which is what the Spearman-Brown
      ## correction needs. Carrying icc, k and the derived reliability into the grid means a
      ## reader can see exactly what turned the observed R2 into the corrected one, rather
      ## than having to trust a single adjusted number.
      icc_ct <- NA_real_; k_ct <- NA_real_
      if(use_icc){
        r <- ICC[cohort==coh & trait==tr]
        if(nrow(r)){
          icc_ct <- r$icc[1]
          if(all(c("n_readings","n_mothers") %in% names(r)) && r$n_mothers[1] > 0)
            k_ct <- r$n_readings[1] / r$n_mothers[1]
        }
      }
      for(dc in cols){
        tt <- momi_transfer(dm$z, dm[[dc]], cov=dm[, .(AGE)])
        ## ---- how many readings does THIS definition actually average? ----
        ## k_ct above is the cohort's mean antenatal readings per mother, which is the right
        ## k only for definitions that average ALL of them. Applying it to every definition
        ## (as this did until 2026-07-19) under-corrects `last` by ~2x -- the mirror image of
        ## the single-ICC bug. From B01:
        ##   mean  = mnv(), mean of all antenatal readings          -> k_ct
        ##   resid = mean of GA-residuals over all antenatal        -> k_ct (needs finite GA,
        ##           so marginally fewer readings; treated as k_ct, a slight over-estimate)
        ##   last  = lv(), the LAST non-missing reading             -> k = 1
        ##   ge20/lt20/tri1/tri2/tri3 = wm(), mean over a WINDOW    -> k unknown to this
        ##           script; bp_icc carries only the overall reading count, not a per-window
        ##           one. Report NA rather than a number we cannot justify. Recovering these
        ##           means B01 emitting a per-definition reading count, which re-locks B01
        ##           for a supplementary column -- not worth it unless S4 is cited on them.
        ddef <- relab(dc, pre)
        k_def <- if(ddef %in% c("mean","resid")) k_ct
                 else if(ddef %in% c("last","early")) 1
                 else NA_real_
        rel_def  <- if(is.na(k_def)) NA_real_ else momi_rel_mean(icc_ct, k_def)
        r2d_def  <- if(use_icc && !is.na(k_def))
                      round(100*momi_disattenuate(tt$incR2, icc_ct, k_def), 3) else NA_real_
        rows[[paste(sid,coh,dc)]] <- data.table(
          score_id=sid, trait=tr, score_anc=an, cohort=coh,
          cohort_anc=unname(MOMI_ANC[coh]), definition=ddef,
          n=tt$n, incR2=tt$incR2, R2pct=round(100*tt$incR2,3),
          beta=tt$beta, se=tt$se, F=tt$F, p=tt$p,
          ## OBSERVED R2 is R2pct above — the measured quantity, and the ONLY one Table 2
          ## reports (decision 2026-07-19). The four columns here are the measurement-error
          ## correction shown in full, so a reader can reproduce it from the row: the
          ## single-reading ICC, the readings THIS definition averages, the resulting
          ## reliability of that average (Spearman-Brown), and only then the corrected R2.
          ## k_readings is NA -> reliability and disattenuated R2 are NA. That is deliberate.
          icc=icc_ct,
          k_readings=round(k_def,2),
          reliability_mean=round(rel_def,3),
          R2pct_disatt=r2d_def)
      }
    }
    rbindlist(rows, fill=TRUE)
  }

  ## ---- antenatal grid: the paper's grid, contract unchanged ----
  grid <- make_grid(ANTE_DEFS, use_icc=TRUE)
  momi_save_intermediate(grid, "transfer_grid", P)
  out  <- momi_write_table(grid[order(trait,cohort,definition,-R2pct)], "S4_transfer_grid", P)

  ## ---- postnatal grid: separate object, no downstream consumer yet ----
  grid_pn <- make_grid(POST_DEFS, use_icc=FALSE)
  momi_save_intermediate(grid_pn, "transfer_grid_pn", P)
  out_pn  <- momi_write_table(grid_pn[order(trait,cohort,definition,-R2pct)],
                              "QC_transfer_grid_pn", P)

  ## ---- console sanity: the antenatal/postnatal contrast this was built to test ----
  cat("\nantenatal grid: ", nrow(grid), " cells (expect ", 8*5*length(ANTE_DEFS),
      " = 8 scores x 5 cohorts x ", length(ANTE_DEFS), " defs)\n", sep="")
  cat("postnatal grid: ", nrow(grid_pn), " cells\n", sep="")

  cat("\nPRIMARY definitions (", paste(MOMI_BPDEF_PRIMARY, collapse=", "),
      "), chosen SBP instrument, coverage and R2:\n", sep="")
  prim <- grid[definition %in% MOMI_BPDEF_PRIMARY & trait=="SBP"]
  prim <- prim[score_id == vapply(cohort, function(c) momi_instrument(c,"SBP"), "")]
  print(dcast(prim, cohort + score_id ~ definition, value.var=c("n","R2pct")))
  cat("\nn tells you whether a definition is usable in that cohort. resid needs a valid GA;\n",
      "if its n is far below the cohort's analytic N, GA quality is limiting the headline\n",
      "exposure and that has to be reported.\n", sep="")

  cat("\nmean-definition R2 by cohort, CHOSEN SBP instrument — antenatal vs postnatal:\n")
  cmp <- rbindlist(list(
    grid   [definition=="mean" & trait=="SBP", .(cohort, score_id, window="antenatal", n, R2pct)],
    grid_pn[definition=="mean" & trait=="SBP", .(cohort, score_id, window="postnatal", n, R2pct)]),
    fill=TRUE)
  cmp <- cmp[score_id == vapply(cohort, function(c) momi_instrument(c,"SBP"), "")]
  print(dcast(cmp, cohort + score_id ~ window, value.var=c("n","R2pct")))
  cat("\nIf R2 is HIGHER postnatally, the score predicts the non-pregnant state better --\n",
      "which would be direct evidence that pregnancy physiology, not ancestry, drives the\n",
      "portability gap. Postnatal N is ~0 in GAPPS by construction.\n", sep="")

  list(n=nrow(grid),
       key=sprintf("cells=%d scores=%d cohorts=%d defs=%d | postnatal cells=%d",
                   nrow(grid), uniqueN(grid$score_id), uniqueN(grid$cohort),
                   uniqueN(grid$definition), nrow(grid_pn)),
       outputs=c("transfer_grid.rds","transfer_grid_pn.rds",
                 basename(out), basename(out_pn)))
})
