## -------------------------------------------------------------------------------
## 01_construct_cohort_and_phenotypes.R
##
## Builds the genetic-analysis sample and the phenotypes.
##
##   1. the eligible sample: first recorded pregnancy, not a known multiple pregnancy,
##      at least one valid antenatal blood-pressure reading, at least one genotype record;
##   2. the prespecified sample-cleaning rules, applied cohort by cohort:
##        R1  a woman is excluded when every genotype record she has failed genotype QC.
##            A failed record belonging to a woman who also has a passing record on the
##            other platform is dropped from scoring and she is retained;
##        R2  a woman whose two platform records are genetically discordant is excluded,
##            and her two scores are never averaged;
##        R3  connected components of genetically identical records. If any edge of a
##            component cannot be resolved to a single woman, every record of that
##            component is excluded; otherwise exactly one record is retained, chosen by
##            earliest eligible pregnancy, then earliest enrolment, then lower genotype
##            missingness, then the first identifier. No blood pressure, outcome
##            availability or outcome value enters this choice;
##        R4  ordinary first- and second-degree relatives are never excluded;
##   3. maternal age and the outcome definitions: preterm birth; low birth weight
##      (<2,500 g); small for gestational age (<10th centile); and birth weight in grams
##      restricted to a 500-6,500 g window, the last three among live births.
##
## Writes the cleaned-sample identifier list and the excluded-genotype-record list that
## the later scripts read, together with aggregate summaries of the cleaning.
## -------------------------------------------------------------------------------
if (!exists("MOMI_ROOT")) MOMI_ROOT <- getwd()
if (!exists("config")) source(file.path(MOMI_ROOT,
  if (file.exists(file.path(MOMI_ROOT, "config.R"))) "config.R" else "config.example.R"))
source(file.path(MOMI_ROOT, "analysis", "00_functions.R"))

COHS     <- COH5
KEEP_CLS <- c("same woman, repeat enrolment", "same pregnancy recorded twice")

S   <- eligible_sample(EPI, SSC)
AN  <- S$AN[cohort %chin% COHS]
GEN <- S$GEN[cohort %chin% COHS]

## Delivery and enrolment dates. Used only to decide which of two genetically identical
## records is retained; no outcome value is read here.

MOMI_NA <- c("-88","-77","-99","NA","na",".","")
DATE_COLS <- c("PARTICIPANT_ID","PREGNANCY_ID","DEL_DATE","DATE_LMP","GAGEBRTH_NEW","BWT_MEASURE_DATE","VISITDT")
hdr <- names(fread(EPI, nrows=0)); miss <- setdiff(DATE_COLS, hdr)
if(length(miss)) stop("EPI lacks columns: ", paste(miss, collapse=", "))
E <- fread(EPI, na.strings=MOMI_NA, select=DATE_COLS)
dt_ <- function(x){ if(inherits(x, c("Date","IDate"))) return(as.Date(x)); suppressWarnings(as.Date(as.character(x), format="%Y-%m-%d")) }
nm_ <- function(x) suppressWarnings(as.numeric(x))
E[, IID := as.character(PARTICIPANT_ID)]
E[, `:=`(DELd=dt_(DEL_DATE), LMPd=dt_(DATE_LMP), BWTd=dt_(BWT_MEASURE_DATE), VISd=dt_(VISITDT), GAd=nm_(GAGEBRTH_NEW))]
E1 <- E[PREGNANCY_ID == 1]
fn <- function(v){ w <- which(!is.na(v)); if(length(w)) v[w[1]] else v[NA_integer_][1] }
DT <- E1[, .(delivery_date=fn(DELd), lmp=fn(LMPd), bwt=fn(BWTd), ga_birth=fn(GAd),
             enrolment_date=suppressWarnings(min(VISd, na.rm=TRUE))), by=IID]

DT[, dd_source := fifelse(!is.na(delivery_date), "recorded delivery date", NA_character_)]
DT[is.na(delivery_date) & !is.na(lmp) & is.finite(ga_birth),
   `:=`(delivery_date = lmp + ga_birth, dd_source = "imputed: LMP + gestational age at birth")]
DT[is.na(delivery_date) & !is.na(bwt), `:=`(delivery_date = bwt, dd_source = "imputed: birthweight-measurement date")]
DT[is.na(delivery_date), dd_source := "no delivery date"]
DT[!is.finite(enrolment_date), enrolment_date := as.Date(NA)]

## The cleaning rules, cohort by cohort.

FLAGS <- list(); COMP <- list(); CMPSUM <- list(); PAIRS <- list(); DROPREC <- list(); DDP <- list()
for(coh in COHS){
  W <- file.path(F4B, "work", coh)
  need <- c("sample_flags.tsv", "pairs_identical_linkage.tsv", "pairs_person_level.tsv")
  for(f in need) if(!file.exists(file.path(W, f))) stop("genotype quality-control file missing: ", file.path(W, f))
  fl <- fread(file.path(W, "sample_flags.tsv"), colClasses=list(character="ID"))
  fl[, ID := normid(ID)]
  idl <- fread(file.path(W, "pairs_identical_linkage.tsv"), colClasses=list(character=c("a", "b")))
  pp <- fread(file.path(W, "pairs_person_level.tsv"), colClasses=list(character=c("a", "b")))
  if(nrow(idl)) idl[, `:=`(a=normid(a), b=normid(b))]
  if(nrow(pp)) pp[, `:=`(a=normid(a), b=normid(b))]

  sm <- function(pfx, plat){ f <- file.path(W, paste0(pfx, ".smiss"))
    if(!file.exists(f)) return(data.table(ID=character(0), platform=character(0), F_MISS=numeric(0)))
    x <- rd_tab(f, c("#FID", "FID", "#IID", "IID")); data.table(ID=x$IID, platform=plat, F_MISS=as.numeric(x$F_MISS)) }
  MISSR <- rbind(sm("lpqc", "low-pass"), sm("gsaqc", "GSA"))

  A <- AN[cohort == coh]
  F <- merge(data.table(ID=A$IID, cohort=coh, technology=A$technology), fl, by="ID", all.x=TRUE)
  F[is.na(lp), lp := FALSE]; F[is.na(gsa), gsa := FALSE]
  F[is.na(lp_qc_fail), lp_qc_fail := FALSE]; F[is.na(gsa_qc_fail), gsa_qc_fail := FALSE]
  F[is.na(dual_discordant), dual_discordant := FALSE]
  F[, n_records := as.integer(lp) + as.integer(gsa)]
  F[, n_records_ok := as.integer(lp & !lp_qc_fail) + as.integer(gsa & !gsa_qc_fail)]

  F[, excl_qc_sole_record := n_records > 0 & n_records_ok == 0]
  F[, record_dropped_partial := n_records_ok > 0 & n_records_ok < n_records]

  F[, excl_dual_discordant := dual_discordant == TRUE & excl_qc_sole_record == FALSE]

  ids <- unique(c(idl$a, idl$b))
  comp <- data.table(ID=character(0), component=integer(0))
  csum <- data.table()
  if(exists("mem", inherits=FALSE)) rm(mem)
  if(length(ids)){

    root <- setNames(ids, ids)
    findr <- function(x){ while(root[[x]] != x){ root[[x]] <<- root[[root[[x]]]]; x <- root[[x]] }; x }
    for(k in seq_len(nrow(idl))){ ra <- findr(idl$a[k]); rb <- findr(idl$b[k]); if(ra != rb) root[[rb]] <- ra }
    cid <- vapply(ids, findr, character(1))
    comp <- data.table(ID=ids, root=unname(cid))[, component := .GRP, by=root][, .(ID, component)]
    idl[comp, ca := i.component, on=c(a="ID")]

    cl <- idl[, .(n_edges=.N, n_bad=sum(!(linkage %chin% KEEP_CLS)),
                  classes=paste(sort(unique(linkage)), collapse="; "),
                  any_possible_mz=any(possible_mz %in% c(TRUE, "TRUE"))), by=.(component=ca)]
    cl[, resolution := fifelse(n_bad > 0, "excluded (identity not verifiable)", "one woman: retain a single record")]
    mem <- merge(comp, F[, .(ID, eligible=TRUE)], by="ID", all.x=TRUE)[is.na(eligible), eligible := FALSE]
    mem[, size := .N, by=component]
    mem <- merge(mem, DT[, .(ID=IID, delivery_date, enrolment_date, dd_source)], by="ID", all.x=TRUE)
    mem[is.na(dd_source), dd_source := "no delivery date"]
    mm <- MISSR[, .(f_miss=min(F_MISS, na.rm=TRUE)), by=ID]; mm[!is.finite(f_miss), f_miss := NA_real_]
    mem <- merge(mem, mm, by="ID", all.x=TRUE)
    mem <- merge(mem, cl[, .(component, resolution)], by="component", all.x=TRUE)

    ex12 <- F[excl_qc_sole_record == TRUE | excl_dual_discordant == TRUE, ID]
    mem[, candidate := eligible == TRUE & !(ID %chin% ex12)]

    mem[, cand_rank := fifelse(candidate, 0L, fifelse(eligible == FALSE, 1L, 2L))]
    setorderv(mem, c("component", "cand_rank", "delivery_date", "enrolment_date", "f_miss", "ID"),
              c(1L, 1L, 1L, 1L, 1L, 1L), na.last=TRUE)
    mem[, retained_in_component := FALSE]
    mem[resolution == "one woman: retain a single record", retained_in_component := seq_len(.N) == 1L, by=component]
    mem[, retained_for_standardization_only := retained_in_component & !candidate]

    mem[, dd_strict := fifelse(dd_source == "recorded delivery date", delivery_date, as.Date(NA))]
    setorderv(mem, c("component", "cand_rank", "dd_strict", "enrolment_date", "f_miss", "ID"),
              c(1L, 1L, 1L, 1L, 1L, 1L), na.last=TRUE)
    mem[, retained_strict := FALSE]
    mem[resolution == "one woman: retain a single record", retained_strict := seq_len(.N) == 1L, by=component]
    DDP[[coh]] <- mem[, .(cohort=coh, retainable=resolution == "one woman: retain a single record",
                          dd_source, retained_in_component, retained_strict)]
    setorderv(mem, c("component", "cand_rank", "delivery_date", "enrolment_date", "f_miss", "ID"),
              c(1L, 1L, 1L, 1L, 1L, 1L), na.last=TRUE)
    csum <- merge(cl, mem[, .(n_members=.N, n_eligible=sum(eligible), n_candidates=sum(candidate),
                              n_retained=sum(retained_in_component),
                              n_retained_std_only=sum(retained_for_standardization_only)), by=component],
                  by="component")[, cohort := coh][]
    F[mem[resolution != "one woman: retain a single record"], excl_identity_unresolved := TRUE, on="ID"]
    F[mem[resolution == "one woman: retain a single record" & eligible == TRUE & retained_in_component == FALSE],
      excl_duplicate_record := TRUE, on="ID"]
    F[mem[resolution == "one woman: retain a single record" & retained_for_standardization_only == TRUE],
      excl_duplicate_record := TRUE, on="ID"]
    F[mem, component := i.component, on="ID"]
    COMP[[coh]] <- mem[, .(cohort=coh, component, ID, eligible, candidate, size, resolution, retained_in_component,
                           retained_for_standardization_only, delivery_date, dd_source, enrolment_date, f_miss)]
    CMPSUM[[coh]] <- csum
  }
  for(v in c("excl_identity_unresolved", "excl_duplicate_record")) if(!v %in% names(F)) F[, (v) := FALSE]
  F[is.na(excl_identity_unresolved), excl_identity_unresolved := FALSE]
  F[is.na(excl_duplicate_record), excl_duplicate_record := FALSE]

  F[, exclusion_reason := fifelse(excl_qc_sole_record, "genotype QC: only record failed",
                          fifelse(excl_dual_discordant, "discordant dual-platform genotypes",
                          fifelse(excl_identity_unresolved, "identical genotypes, identity not verifiable",
                          fifelse(excl_duplicate_record, "duplicate record of a retained genetic individual", "retained"))))]
  F[, retained := exclusion_reason == "retained"]

  rel <- if(nrow(pp)) pp[class %chin% c("first degree", "second degree")] else pp[0]
  F[, related_to_other_retained := FALSE]
  if(nrow(rel)){ rt0 <- F[retained == TRUE, ID]
    both <- rel[a %chin% rt0 & b %chin% rt0]
    F[ID %chin% unique(c(both$a, both$b)), related_to_other_retained := TRUE] }
  FLAGS[[coh]] <- F

  memx <- if(exists("mem", inherits=FALSE) && is.data.table(mem)) mem else data.table(ID=character(0), resolution=character(0), retained_in_component=logical(0))
  drop_ids <- unique(c(F[excl_qc_sole_record == TRUE | excl_dual_discordant == TRUE, ID],
                       memx[resolution != "one woman: retain a single record", ID],
                       memx[resolution == "one woman: retain a single record" & retained_in_component == FALSE, ID]))
  flr <- copy(fl)[, `:=`(lp=lp %in% c(TRUE, "TRUE"), gsa=gsa %in% c(TRUE, "TRUE"),
                         lp_qc_fail=lp_qc_fail %in% c(TRUE, "TRUE"), gsa_qc_fail=gsa_qc_fail %in% c(TRUE, "TRUE"))]
  recs <- rbind(flr[lp == TRUE, .(ID, platform="lpwgs_dosage", qc_fail=lp_qc_fail)],
                flr[gsa == TRUE, .(ID, platform="gsa", qc_fail=gsa_qc_fail)])
  recs[, drop_reason := fifelse(ID %chin% drop_ids, "participant excluded from the genetic-analysis sample",
                        fifelse(qc_fail, "record failed genotype QC", NA_character_))]
  DROPREC[[coh]] <- recs[!is.na(drop_reason), .(cohort=coh, ID, platform, drop_reason)]
  rm(mem)
  rt <- F[retained == TRUE, ID]
  PAIRS[[coh]] <- if(nrow(pp)) pp[, .(pairs=.N, pairs_both_retained=sum(a %chin% rt & b %chin% rt)), by=class][, cohort := coh][] else
                  data.table(class=character(0), pairs=integer(0), pairs_both_retained=integer(0), cohort=character(0))
  cat(sprintf("%s: eligible %d; excluded QC %d, discordant %d, identity %d, duplicate %d; retained %d\n",
              coh, nrow(F), sum(F$excl_qc_sole_record), sum(F$excl_dual_discordant),
              sum(F$exclusion_reason == "identical genotypes, identity not verifiable"),
              sum(F$exclusion_reason == "duplicate record of a retained genetic individual"), sum(F$retained)))
}
FL <- rbindlist(FLAGS, fill=TRUE)

## Identifier lists for the later scripts. These are participant-level files: they are
## written into the derived-data directory, which is not distributed.

fwrite(FL, file.path(DERIVED, "sample_flags.tsv"), sep="\t")
writeLines(FL[retained == TRUE, ID], file.path(DERIVED, "cleaned_sample_ids.txt"))
writeLines(FL[retained == FALSE, ID], file.path(DERIVED, "excluded_ids.txt"))
REC <- rbindlist(DROPREC)
fwrite(REC[, .(ID, platform, drop_reason, cohort)], file.path(DERIVED, "dropped_genotype_records.tsv"), sep="\t")
fwrite(REC[, .(ID, platform)], file.path(DERIVED, "dropped_genotype_records_ids.tsv"), sep="\t")
if(length(COMP)) fwrite(rbindlist(COMP), file.path(DERIVED, "identity_components.tsv"), sep="\t")

W_(FL[, .(n_eligible = .N, retained = sum(retained),
          excluded_genotype_qc = sum(excl_qc_sole_record),
          excluded_discordant_dual = sum(excl_dual_discordant),
          excluded_identity_unresolved = sum(exclusion_reason ==
            "identical genotypes, identity not verifiable"),
          excluded_duplicate_record = sum(exclusion_reason ==
            "duplicate record of a retained genetic individual")), by = cohort],
   "sample_cleaning_summary.tsv")
W_(rbindlist(PAIRS), "relatedness_pairs.tsv")

## The exposures: each woman's mean of her valid antenatal systolic and diastolic
## readings, together with the genotyping technology available to her.

ALLW <- S$ALLW
AN   <- S$AN[cohort %chin% COH5]
AN   <- merge(AN, ALLW[, .(IID, S_mean, D_mean)], by="IID", all.x=TRUE)
ELIG_N <- nrow(AN)

## Maternal age and the outcomes.

hdr  <- names(fread(EPI, nrows=0, showProgress=FALSE))
miss <- setdiff(OUT_COLS, hdr); if(length(miss)) stop("EPI lacks columns: ", paste(miss, collapse=", "))
raw <- fread(EPI, na.strings=MOMI_NA, select=OUT_COLS, showProgress=FALSE)[PREGNANCY_ID == 1]
getcol <- function(nm) if(nm %in% names(raw)) num(raw[[nm]]) else rep(NA_real_, nrow(raw))
raw[, IID := as.character(PARTICIPANT_ID)]
garead <- getcol("GA_HDLK_NEW")
raw[, GA := ifelse(garead >= 0 & garead <= 315, garead, NA_real_)]
raw[, VIS := getcol("VISITDT")]
raw[, AGEv := getcol("PW_AGE")];       raw[, PTBv  := getcol("PTB_NEW")]
raw[, BWTv := getcol("BIRTH_WEIGHT")]; raw[, SGAv  := getcol("SGA_10_NEW")]
raw[, BOUTv := getcol("BIRTH_OUTCOME")]
setorder(raw, IID, VIS, GA)
PH <- raw[, .(AGE=firstnn(AGEv), PTB=firstnn(PTBv), BWTraw=firstnn(BWTv),
              SGA=as.integer(firstnn(SGAv) == 1), BOUT=firstnn(BOUTv)), by=IID]
PH[, BWT := ifelse(BWTraw >= 500 & BWTraw <= 6500, BWTraw, NA_real_)]
PH[, livebirth := as.integer(BOUT == 1)]
AN <- merge(AN, PH[, .(IID, AGE, PTB, BWT, SGA, livebirth)], by="IID", all.x=TRUE)

## Restrict to the cleaned sample.

keep_raw <- as.character(fread(KEEPF, header=FALSE, colClasses="character", showProgress=FALSE)[[1]])
keep_raw <- keep_raw[nzchar(keep_raw)]
keep_ids <- unique(keep_raw)
n_keep_dup <- length(keep_raw) - length(keep_ids)
n_keep_unmatched <- sum(!(keep_ids %chin% AN$IID))
AN <- AN[IID %chin% keep_ids]
CLEAN_N <- nrow(AN)

saveRDS(AN, file.path(DERIVED, "analysis_sample.rds"))
W_(AN[, .(n = .N,
          with_mean_SBP = sum(is.finite(S_mean)), with_mean_DBP = sum(is.finite(D_mean)),
          with_age = sum(is.finite(AGE)), with_preterm_birth = sum(!is.na(PTB)),
          live_births = sum(livebirth == 1, na.rm = TRUE),
          with_birth_weight = sum(is.finite(BWT)), with_sga = sum(!is.na(SGA))),
       by = cohort], "phenotype_availability.tsv")
cat(sprintf("01: cleaned analysis sample, %d women in %d cohorts\n", nrow(AN), length(COHS)))
