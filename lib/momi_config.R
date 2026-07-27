#!/usr/bin/env Rscript
# ============================================================
# momi_config.R — SINGLE SOURCE OF TRUTH for the BP -> perinatal paper.
# Every 05_prs / 06_phase2 R script should `source()` this so cohort labels,
# ancestry, GA cutoffs, outcome definitions, the CHOSEN INSTRUMENT per cohort,
# and the confounder set are identical everywhere. Change a decision HERE, once.
# Aligned to the MOMI consortium paper (Tang et al., JOGH 2026;16:05002) and the
# July-2026 meeting decisions.
#   usage:  source(file.path(PIPE,"lib/momi_config.R"))
# ============================================================

## ---- missing-value tokens (consistent across all reads) ----
MOMI_NA <- c("-88","-77","-99","NA","na",".","")

## ---- SITE_CODE <-> cohort (internal label) and consortium display name ----
MOMI_SITE <- c("1"="GAPPS-Zambia","2"="GAPPS-Bangladesh","3"="AMANHI-Pakistan",
               "4"="AMANHI-Bangladesh","5"="AMANHI-Pemba","6"="THSTI-India")
# Display names for Table 1 (do NOT use as keys in code):
MOMI_DISPLAY <- c("GAPPS-Zambia"="ZAPPS (Lusaka, Zambia)",
                  "GAPPS-Bangladesh"="PreSSMat (Matlab, Bangladesh)",
                  "AMANHI-Pakistan"="AMANHI (Karachi, Pakistan)",
                  "AMANHI-Bangladesh"="AMANHI (Sylhet, Bangladesh)",
                  "AMANHI-Pemba"="AMANHI (Pemba, Tanzania)",
                  "THSTI-India"="GARBH-INi (Delhi, India)")

## ---- ancestry ----
MOMI_ANC <- c("AMANHI-Bangladesh"="SAS","AMANHI-Pakistan"="SAS","GAPPS-Bangladesh"="SAS",
              "AMANHI-Pemba"="AFR","GAPPS-Zambia"="AFR")
MOMI_COH_ALL <- names(MOMI_ANC)                 # 5 genotyped cohorts
MOMI_COH_SAS <- names(MOMI_ANC)[MOMI_ANC=="SAS"]
MOMI_COH_AFR <- names(MOMI_ANC)[MOMI_ANC=="AFR"]

## ---- genotype platforms (merged per mother = mean z across these) ----
MOMI_PLATS <- c("gsa","lpwgs_dosage")

## ---- gestational-age cutoffs (GA_HDLK_NEW is in DAYS; Hadlock, per consortium) ----
MOMI_GA <- list(
  tri1_end = 98,    # <14 wk
  wk20     = 140,   # 20 wk  -> chronic (<20) vs gestational (>=20) split (Baqui)
  tri2_end = 196,   # 28 wk
  ptb_days = 259    # <37 wk = preterm (liveborn), per consortium
)

## ---- BP exposure definitions, TIERED (decision 2026-07-19) ----
# PRIMARY
#   resid  THE headline exposure. GA-standardised (residual from a per-cohort spline of BP
#          on gestational age), so it removes the between-cohort measurement-timing
#          difference rather than being confounded by it -- which is the very problem
#          Part I exists to characterise.
#   ge20   second primary: the only window near-universally covered (96-98% of mothers in
#          every cohort).
# AUXILIARY
#   mean   retained for continuity with every previously locked table (was the old default)
#   lt20   AUXILIARY, not primary: AMANHI barely measures BP before 20 wk (7.6 / 1.0 / 1.0%
#          of mothers in AMANHI-B / Pakistan / Pemba, vs 99.9% in GAPPS-B). That gap IS the
#          timing contrast, so it cannot carry a primary exposure.
#   tri1/2/3, last  -- NB tri1 is ~0% covered outside GAPPS; read with the same care.
# NOT GRIDDED
#   early, median, mn2 -- B01 still COMPUTES these columns (they cost nothing, and `early`
#   is wanted for later trajectory analysis); they simply do not enter the grid.
MOMI_BPDEF_PRIMARY <- c("resid","ge20")
MOMI_BPDEF_AUX     <- c("mean","lt20","tri1","tri2","tri3","last")
MOMI_BPDEFS        <- c(MOMI_BPDEF_PRIMARY, MOMI_BPDEF_AUX)
MOMI_BPDEF_DEFAULT <- "resid"
MOMI_BPDEFS_UNUSED <- c("early","median","mn2")   # computed by B01, excluded from the grid

## ---- outcomes (main set; extras go to the supplementary panel) ----
# PTB = GAGEBRTH_NEW < 259 d ; LBW = BIRTH_WEIGHT < 2500 g ; SGA = SGA_10_NEW == 1 (INTERGROWTH <10th)
MOMI_OUTCOMES_MAIN <- c("PTB","LBW","SGA")           # + BWT (continuous companion)
## DECISION 2026-07-20 (Anagh's advisor): HBW, POSTTERM, STILL and MISC are REMOVED. Too few
## events in MOMI to estimate anything, and reporting near-empty cells invites over-reading.
## Consequence to be aware of: the external comparator (Morales-Berstein 2026,
## ref/external_mr_estimates.tsv) uses exactly these as its mirror-image and null controls
## (HBW 0.76 mirrors LBW 1.33; POSTTERM 0.94 mirrors PTB; STILL/MISC are explicit nulls), so
## we can no longer validate by reproducing that pattern. S12/B24 uses a different control
## strategy instead -- see deliv_S12_controls.R.
MOMI_OUTCOMES_DROPPED <- c("HBW","POSTTERM","STILL","MISC")
MOMI_OUTCOMES_SUPP <- c("PROVDEL","BWT","GAd")
# live-birth restriction applies to birthweight-derived outcomes:
MOMI_LIVEBIRTH_OUTCOMES <- c("LBW","SGA","BWT")   # HBW removed 2026-07-20 with the others

## ---- POLYGENIC SCORE PANEL (all scores used anywhere) ----
# columns: id, trait, ancestry, source, role
# `role` separates the BP instrument candidates from the BMI score, which is not a
# candidate instrument at all -- it is the secondary exposure for the (currently blocked)
# S14 MVMR supplement. Without this column a reader of Table S1 has no way to tell why a
# BMI score is sitting in a table of blood-pressure instruments.
# NB: catalog metadata (variant count, method, PMID, discovery/training N) is NOT stored
# here -- it is externally sourced in ref/pgs_catalog_metadata.tsv and joined by B04.
MOMI_PANEL <- data.frame(
  id  = c("PGS004603","PGS004830","PGS002376","PGS003968","PGS005008","PGS002238",
          "PGS004604","PGS004758","PGS002362","PGS003964","PGS002239","PGS000027"),
  trait=c("SBP","SBP","SBP","SBP","SBP","SBP","DBP","DBP","DBP","DBP","DBP","BMI"),
  anc = c("EUR","SAS","EAS","DIVERSE","MULTI","MVP","EUR","SAS","EAS","DIVERSE","MVP","EUR"),
  source=c("Keaton2024","Truong2024","Weissbrod2022","Kurniansyah2023","Gunn2024","Breeyear2022",
           "Keaton2024","Truong2024","Weissbrod2022","Kurniansyah2023","Breeyear2022","Khera2019"),
  role = c(rep("BP instrument candidate", 11), "MVMR secondary exposure (BMI)"),
  stringsAsFactors=FALSE)

## ---- CHOSEN INSTRUMENT per cohort (Part-I decision) ----
## ANCESTRY-MATCHED, decided 2026-07-21: each cohort takes the score whose training ancestry
## matches it as closely as an available GWAS allows.
##   South Asian cohorts (AMANHI-Bangladesh, AMANHI-Pakistan, GAPPS-Bangladesh) -> SAS score.
##   African cohorts (AMANHI-Pemba, GAPPS-Zambia) -> EUR score. There is no African BP GWAS
##     of comparable size, so EUR is the closest available proxy; this is stated as a caveat
##     in the transportability section rather than presented as a match.
## Selection remains fully disclosed by S4/F2, which report EVERY score in every cohort, so
## STROBE-MR item 6b holds regardless of which one is designated the instrument.
##
## PREVIOUS decision (archived here, 2026-07-19): EUR primary everywhere except GAPPS-B=SAS.
##   MOMI_INSTR_ANC <- c("AMANHI-Bangladesh"="EUR","AMANHI-Pakistan"="EUR","GAPPS-Bangladesh"="SAS",
##                       "AMANHI-Pemba"="EUR","GAPPS-Zambia"="EUR")
MOMI_INSTR_ANC <- c("AMANHI-Bangladesh"="SAS","AMANHI-Pakistan"="SAS","GAPPS-Bangladesh"="SAS",
                    "AMANHI-Pemba"="EUR","GAPPS-Zambia"="EUR")
# resolve to a PGS id for a cohort+trait using the chosen ancestry:
momi_instrument <- function(cohort, trait){
  a <- MOMI_INSTR_ANC[[cohort]]; if(is.null(a)) a <- "EUR"
  hit <- MOMI_PANEL$id[MOMI_PANEL$trait==trait & MOMI_PANEL$anc==a]
  if(!length(hit)) hit <- MOMI_PANEL$id[MOMI_PANEL$trait==trait & MOMI_PANEL$anc=="EUR"]
  hit[1]
}
# EUR ids for convenience (primary instrument, most scripts):
MOMI_EUR <- list(SBP="PGS004603", DBP="PGS004604")
MOMI_SAS <- list(SBP="PGS004830", DBP="PGS004758")
MOMI_BMI <- "PGS000027"

## ---- observational confounder set (keep PE & GDM as MEDIATORS: never adjust) ----
# Data-driven set (missingness audit, July 2026): PARITY is 18-34% missing at 4/5
# genotyped sites -> use GRAVIDITY (0% missing); SMOK_FREQ ~99% missing -> dropped;
# WEALTH_INDEX 50% missing in Zambia -> use PW_EDUCATION for SES; BMI computed from
# MAT_HEIGHT/MAT_WEIGHT (well-measured except Pemba). Decision: all-core incl. BMI.
#   MOMI_CONF_CORE = WITHIN-cohort covariates (vary within a cohort; used in per-cohort
#     transferability + MR first/second stage).
#   MOMI_CONF_SITE = the site term, added ONLY for POOLED observational models (it is
#     constant within a cohort, so it must NOT enter per-cohort fits).
MOMI_CONF_CORE  <- c("AGE","GRAV","EDU","BMI")
MOMI_CONF_SITE  <- "cohort"                       # factor() in pooled models
MOMI_CONF_EXTRA <- character(0)                   # smoking unusable (~99% missing)
# formula piece for pooled models (adds site); per-cohort models pass MOMI_CONF_CORE.
momi_conf_formula <- function(pooled=TRUE){
  paste(c(MOMI_CONF_CORE, if(pooled) MOMI_CONF_SITE, MOMI_CONF_EXTRA), collapse="+")
}

## ---- shared getz(): merged avg-z PRS per IID across platforms ----
momi_zin <- function(x){s<-sd(x,na.rm=TRUE); if(is.na(s)||s==0) x*NA else (x-mean(x,na.rm=TRUE))/s}
momi_getz <- function(pid, coh, sscore_dir, plats=MOMI_PLATS, infant=FALSE){
  suf <- if(infant) "__infant" else ""
  L <- lapply(plats, function(p){
    f <- file.path(sscore_dir, sprintf("%s__%s__%s%s.sscore", pid, coh, p, suf))
    if(!file.exists(f)) return(NULL)
    s <- data.table::fread(f)
    data.table::data.table(IID=as.character(s[[1]]), z=as.numeric(momi_zin(s$SCORE1_AVG)))})
  L <- L[!vapply(L, is.null, TRUE)]; if(!length(L)) return(NULL)
  ## M2 (2026-07-26): averaging two unit-variance platform scores gives SD<1 for dual-platform
  ## mothers, so the merged z is not unit-variance and its scale depends on each cohort's
  ## dual-platform fraction. Re-standardise the merged score to mean 0, SD 1 within cohort so
  ## that "per SD of the score" is honest and beta is comparable across cohorts. (Scale-
  ## invariant quantities -- Wald ratio, F, incremental R^2 -- are unaffected.)
  data.table::rbindlist(L)[is.finite(z), .(z=mean(z)), by=IID][, z := momi_zin(z)][]
}

## ---- inverse-variance meta ----
momi_metaiv <- function(b,se){ok<-is.finite(b)&is.finite(se)&se>0; b<-b[ok]; se<-se[ok]
  if(!length(b)) return(NULL); w<-1/se^2; list(b=sum(w*b)/sum(w), se=sqrt(1/sum(w)), k=length(b))}

invisible(TRUE)
