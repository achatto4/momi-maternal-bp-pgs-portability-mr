## -------------------------------------------------------------------------------
## 02_prepare_genetic_covariates_and_pgs.R
##
## Attaches the genetic covariates and the polygenic scores to the cleaned sample:
##   * principal components 1-5, one row per woman, from the projection prepared upstream;
##   * the merged polygenic score for each cohort and trait. The score is standardised
##     within cohort and platform, averaged across platforms for a woman scored on both,
##     and re-standardised within cohort, with the genotype records excluded by 01 removed
##     before any standardisation, so that one standard deviation means the same thing
##     here as in the portability analysis.
##
## Writes a per-cohort, per-trait record of how each score was constructed.
## -------------------------------------------------------------------------------
if (!exists("MOMI_ROOT")) MOMI_ROOT <- getwd()
if (!exists("config")) source(file.path(MOMI_ROOT,
  if (file.exists(file.path(MOMI_ROOT, "config.R"))) "config.R" else "config.example.R"))
source(file.path(MOMI_ROOT, "analysis", "00_functions.R"))

AN      <- readRDS(file.path(DERIVED, "analysis_sample.rds"))
CLEAN_N <- nrow(AN)

pct <- sub("\\.rds$", ".tsv", PCSF)
if(file.exists(PCSF)){
  PC_USED <- PCSF
  PL <- if(grepl("\\.rds$", PCSF, ignore.case=TRUE)) as.data.table(readRDS(PCSF)) else
        fread(PCSF, colClasses=list(character="IID"), showProgress=FALSE)
} else if(file.exists(pct)){
  PC_USED <- pct
  PL <- fread(pct, colClasses=list(character="IID"), showProgress=FALSE)
} else stop("corrected PC file missing: ", PCSF)
PL <- as.data.table(PL)
if(!"IID" %in% names(PL)) stop("PC table has no IID column")
PL[, IID := normid(IID)]
miss_pc <- setdiff(PCS5, names(PL)); if(length(miss_pc)) stop("PC table lacks: ", paste(miss_pc, collapse=", "))
PCTAB <- unique(PL[, c("IID", PCS5), with=FALSE], by="IID")
n_pc_dup <- nrow(PL) - nrow(PCTAB)
AN <- merge(AN, PCTAB, by="IID", all.x=TRUE)
AN[, has_pc5 := Reduce(`&`, lapply(PCS5, function(p) is.finite(get(p))))]
n_no_pc <- sum(!AN$has_pc5)

NEED <- unique(c(INSTR[COH5]$SBP, INSTR[COH5]$DBP))
ZL <- list(); ZALT <- list(); SC <- list()
for(coh in COH5) for(sid in NEED){
  z <- getz(sid, coh); if(is.null(z)) stop("no .sscore files for ", sid, " / ", coh)
  ZL[[paste(coh, sid)]]   <- z
  ZALT[[paste(coh, sid)]] <- z_alt(sid, coh)
}
for(coh in COH5) for(tr in TRAITS){
  sid <- INSTR[coh][[tr]]
  z   <- ZL[[paste(coh, sid)]]; fls <- score_files(sid, coh)
  za  <- merge(AN[cohort == coh, .(IID)], z, by="IID")
  SC[[length(SC) + 1]] <- data.table(
    cohort=coh, cohort_display=unname(DISPLAY[coh]), site=unname(SITE_SHORT[coh]),
    ancestry_group=unname(GROUP[coh]), bp_trait=tr, exposure=unname(BPCOL[[tr]]),
    pgs_id=sid, prespecified_pgs_id=INSTR[coh][[tr]],
    n_sscore_files=length(fls), sscore_files=jn(basename(fls)),
    sscore_md5=jn(vapply(fls, md5, character(1))),
    n_scored_records=nrow(z), n_dual_platform_scored=sum(z$n_plat == 2),
    mean_scored=mean(z$z), sd_scored=sd(z$z),
    n_in_cleaned_sample=nrow(za), mean_in_cleaned_sample=mean(za$z), sd_in_cleaned_sample=sd(za$z))
}
SC <- rbindlist(SC)

saveRDS(list(AN = AN, ZL = ZL, ZALT = ZALT),
        file.path(DERIVED, "analysis_sample_with_scores.rds"))
W_full(SC, "score_construction.tsv")
cat(sprintf("02: %d cohort x trait scores prepared; %d of %d women have PC1-PC5\n",
            nrow(SC), CLEAN_N - n_no_pc, CLEAN_N))
