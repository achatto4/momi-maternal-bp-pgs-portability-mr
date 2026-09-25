## assemble_joint_pcs.R
##
## The joint principal-component file used by the analysis (one row per woman) and the quality gates that must pass
## before it is used: variant counts, cross-technology frequency agreement, platform separation, duplicates, coverage.
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(here, "helpers", "ids.R")); source(file.path(here, "helpers", "eligible_sample.R"))
args <- commandArgs(TRUE); arg <- function(f, d=NULL){ i <- which(args == f); if(length(i)) args[i + 1] else d }
RUN <- arg("--rundir"); EPI <- arg("--epi"); SSC <- arg("--sscore-dir"); ENC <- arg("--encoding")
COHS <- strsplit(arg("--cohorts"), " ")[[1]]
KEEPF <- arg("--keep"); KEEP_IDS <- normid(fread(KEEPF, header=FALSE, colClasses="character")[[1]])
EXPN <- table(eligible_sample(EPI, SSC)$ALLW[normid(IID) %chin% KEEP_IDS, cohort])
MINV <- c(within=20000L, cross=10000L)

OV_N <- arg("--expected-n"); OV_V <- arg("--min-variants")
if(!is.null(OV_N)){ kv <- strsplit(strsplit(OV_N, ",")[[1]], "="); EXPN <- setNames(as.integer(vapply(kv, `[`, "", 2)), vapply(kv, `[`, "", 1)) }
if(!is.null(OV_V)){ v <- as.integer(strsplit(OV_V, ",")[[1]]); MINV <- c(within=v[1], cross=v[2]) }
AGG <- file.path(RUN, "aggregate"); PLV <- file.path(RUN, "participant_level"); dir.create(PLV, FALSE, TRUE)
rd_all <- function(f){ L <- lapply(COHS, function(c){ p <- file.path(AGG, "cohort", c, f); if(file.exists(p)) fread(p) else NULL })
  L <- L[!vapply(L, is.null, TRUE)]; if(length(L)) rbindlist(L, fill=TRUE) else data.table() }
D <- list(); rule <- function(r, c, v, thr, st, note="") D[[length(D) + 1]] <<- data.table(rule=r, cohort=c, value=v, threshold=thr, status=st, note=note)
if(!is.null(OV_N) || !is.null(OV_V))
  rule("synthetic-test threshold override", "all", sprintf("expected N %s; minimum variants %d / %d", paste(EXPN, collapse="/"), MINV[["within"]], MINV[["cross"]]),
       "no override on MOMI data", "WARN", "--expected-n / --min-variants were given: valid for synthetic tests only")

PC <- rbindlist(lapply(COHS, function(c) fread(file.path(RUN, "work", c, "pcs_joint.tsv"), colClasses=list(character="IID"))))
PC[, IID := normid(IID)]
cov5 <- PC[, .(n=.N, finite_pc1_5=sum(stats::complete.cases(.SD))), by=cohort, .SDcols=paste0("PC", 1:5)]
for(c in COHS){ v <- cov5[cohort == c]
  rule("S5 coverage", c, sprintf("%d women, %d with PC1-PC5", v$n, v$finite_pc1_5), sprintf("%d", EXPN[[c]]),
       if(v$n == EXPN[[c]] && v$finite_pc1_5 == EXPN[[c]]) "PASS" else "STOP") }
rule("S4 duplicates", "all", sprintf("%d duplicated IDs", anyDuplicated(PC$IID) > 0), "0", if(anyDuplicated(PC$IID) == 0) "PASS" else "STOP")
saveRDS(as.data.frame(PC), file.path(PLV, "joint_pcs.rds")); fwrite(PC, file.path(PLV, "joint_pcs.tsv"), sep="\t")

S <- eligible_sample(EPI, SSC)
tm <- merge(PC[, .(IID, technology)], S$ALLW[, .(IID, technology_acc=technology)], by="IID", all.x=TRUE)
nmis <- tm[is.na(technology_acc) | technology != technology_acc, .N]
rule("S7 technology covariate", "all", sprintf("%d mismatches in %d women", nmis, nrow(tm)), "0", if(nmis == 0) "PASS" else "STOP")

VQ <- rd_all("agg_variant_qc.tsv")
for(c in COHS){ v <- VQ[cohort == c]$used; rule("S1 variants (within-cohort PCA)", c, v, sprintf(">= %d", MINV[["within"]]), if(length(v) && v >= MINV[["within"]]) "PASS" else "STOP") }
ce <- file.path(AGG, "cross", "agg_cross_eigen.tsv")
if(file.exists(ce)){ v <- fread(ce)$n_variants[1]; rule("S1 variants (cross-cohort set)", "all", v, sprintf(">= %d", MINV[["cross"]]), if(v >= MINV[["cross"]]) "PASS" else "STOP") } else
  rule("S1 variants (cross-cohort set)", "all", NA, sprintf(">= %d", MINV[["cross"]]), "STOP", "cross-cohort stage missing")

FC <- rd_all("agg_freq_concordance.tsv")
for(c in COHS){ v <- FC[cohort == c]
  rule("S2 frequency correlation", c, sprintf("%.5f", v$r_frequency), ">= 0.98", if(v$r_frequency >= 0.98) "PASS" else "STOP")
  rule("S2 variants beyond |dp| 0.10", c, sprintf("%.3f%%", v$pct_excluded), "<= 1%", if(v$pct_excluded <= 1) "PASS" else "STOP") }
PFq <- rd_all("agg_platform_frequency_pca_variants.tsv")
for(c in COHS){ v <- PFq[cohort == c]
  rule("S2 frequency correlation, genotypes entering the PCA", c, if(nrow(v)) sprintf("%.5f (%d variants)", v$r_frequency, v$variants) else "missing", ">= 0.98",
       if(nrow(v) && v$r_frequency >= 0.98) "PASS" else "STOP")
  rule("S2 variants beyond |dp| 0.10, genotypes entering the PCA", c, if(nrow(v)) sprintf("%.3f%% (sampling alone: %.3f%%)", v$pct_beyond_0.10, v$pct_expected_by_sampling) else "missing",
       "<= 1%", if(nrow(v) && v$pct_beyond_0.10 <= 1) "PASS" else "STOP") }
if(ENC == "dosage"){ DC <- rd_all("agg_array_dosage_vs_hardcall.tsv")
  for(c in COHS){ v <- DC[cohort == c]
    rule("S2 array dosage vs gsa_clean hard calls", c, if(nrow(v)) sprintf("%.5f over %d genotypes", v$concordance, v$genotype_pairs) else "missing", ">= 0.98",
         if(nrow(v) && v$genotype_pairs >= 1000 && v$concordance >= 0.98) "PASS" else "STOP") } }
SP <- rd_all("agg_same_person_grm.tsv")
for(c in COHS){ v <- SP[cohort == c]
  if(nrow(v) && v$n_dual >= 30) rule("S2 same-woman cross-technology relationship", c, sprintf("median %.3f (n = %d)", v$median_same_person_grm, v$n_dual), ">= 0.80",
                                     if(v$median_same_person_grm >= 0.80) "PASS" else "STOP")
  else rule("S2 same-woman cross-technology relationship", c, sprintf("n = %d dual women", if(nrow(v)) v$n_dual else 0L), "applies with >= 30", "NOT APPLICABLE") }

SEPA <- rd_all("agg_separation.tsv")
for(c in COHS){
  s <- SEPA[cohort == c & scope == "all women" & grouping == "record platform" & pc %in% paste0("PC", 1:5)]
  ng <- s[level == "GSA", n][1]; nl <- s[level == "low-pass", n][1]; if(is.na(ng)) ng <- 0L
  e2 <- unique(s[, .(pc, eta2)])
  if(ng >= 5){ m <- max(e2$eta2); rule("S3 eta^2 of record platform, PC1-PC5", c, sprintf("max %.4f (%s)", m, e2$pc[which.max(e2$eta2)]), "<= 0.50 (warn > 0.10)",
                                       if(m > 0.50) "STOP" else if(m > 0.10) "WARN" else "PASS") }
  else rule("S3 eta^2 of record platform, PC1-PC5", c, sprintf("%d array-record women", ng), "applies with >= 5", "NOT APPLICABLE")
  if(ng >= 30 && nl >= 30){ sm <- s[level == "GSA", .(pc, smd_vs_rest)]; m <- max(abs(sm$smd_vs_rest))
    rule("S3 |SMD| array vs low-pass records, PC1-PC5", c, sprintf("max %.3f (%s)", m, sm$pc[which.max(abs(sm$smd_vs_rest))]), "<= 1.0 (warn > 0.25)",
         if(m > 1.0) "STOP" else if(m > 0.25) "WARN" else "PASS") }
  else rule("S3 |SMD| array vs low-pass records, PC1-PC5", c, sprintf("%d array-record women", ng), "applies with >= 30 per group", "NOT APPLICABLE")
}
DU <- rd_all("agg_dual_agreement.tsv")
for(c in COHS){ v <- DU[cohort == c & grepl("^leave-dual-out", fit) & pc %in% paste0("PC", 1:5)]
  if(nrow(v) && min(v$n) >= 30){ m <- max(abs(v$mean_diff_sd)); rule("S3 leave-dual-out |mean difference|, PC1-PC5", c, sprintf("max %.3f SD (%s; n = %d)", m, v$pc[which.max(abs(v$mean_diff_sd))], min(v$n)),
                                                              "<= 0.50 SD (warn > 0.20)", if(m > 0.50) "STOP" else if(m > 0.20) "WARN" else "PASS") }
  else rule("S3 leave-dual-out |mean difference|, PC1-PC5", c, sprintf("n = %d", if(nrow(v)) min(v$n) else 0L), "applies with >= 30", "NOT APPLICABLE") }

DP <- rd_all("agg_duplicates.tsv")
for(c in COHS){ v <- DP[cohort == c]; rule("S4 different women with relationship >= 0.80", c, sprintf("%d pairs; max off-diagonal %.3f", v$pairs_grm_ge_0.80, v$max_offdiag_grm), "0",
                                           if(v$pairs_grm_ge_0.80 == 0 && !isTRUE(v$duplicated_women)) "PASS" else "STOP",
                                           sprintf("liveness: %d of %d same-woman record pairs >= 0.80", v$same_person_pairs_ge_0.80, v$same_person_pairs)) }
RF <- rd_all("agg_reference.tsv")
for(c in COHS){ v <- RF[cohort == c]; rule("S4 identical-genotype pairs between retained women", c, v$identical_pairs_between_retained_women, "0",
                                           if(v$identical_pairs_between_retained_women == 0) "PASS" else "STOP",
                                           sprintf("liveness: %d identical pairs in the relationship file for this cohort", v$identical_pairs_in_file_all_women)) }
EG <- rd_all("agg_eigen.tsv")
for(c in COHS){ v <- EG[cohort == c & pc %in% paste0("PC", 1:5)]; m <- min(v$n_eff)
  rule("n_eff of PC1-PC5 (localisation)", c, sprintf("min %.1f (%s)", m, v$pc[which.min(v$n_eff)]), ">= 30 (warn)", if(m >= 30) "PASS" else "WARN") }

DEC <- rbindlist(D)
overall <- if(any(DEC$status == "STOP")) "STOP" else "PASS"
DEC <- rbind(DEC, data.table(rule="OVERALL", cohort="all", value=overall, threshold="", status=overall,
                             note=sprintf("%d STOP, %d WARN", sum(DEC$status == "STOP"), sum(DEC$status == "WARN"))))
fwrite(DEC, file.path(AGG, "joint_pca_decision.tsv"), sep="\t")
for(f in c("agg_eigen.tsv", "agg_separation.tsv", "agg_dual_agreement.tsv", "agg_vs_previous.tsv", "agg_encoding_sensitivity.tsv", "agg_selection.tsv",
           "agg_reference.tsv", "agg_freq_concordance.tsv", "agg_missingness.tsv", "agg_duplicates.tsv", "agg_projection.tsv", "agg_variant_qc.tsv",
           "agg_same_person_grm.tsv", "agg_cv_folds.tsv", "agg_tracy_widom_1_20.tsv", "agg_array_dosage_vs_hardcall.tsv",
           "agg_platform_frequency_pca_variants.tsv")){
  x <- rd_all(f); if(nrow(x)) fwrite(x, file.path(AGG, sub("^agg_", "joint_pca_", f)), sep="\t") }
VF <- rbindlist(lapply(COHS, function(c){ p <- file.path(AGG, "cohort", c, "variant_funnel.tsv"); if(file.exists(p)) fread(p)[, cohort := c] else NULL }))
if(nrow(VF)) fwrite(VF, file.path(AGG, "joint_pca_variant_funnel.tsv"), sep="\t")
print(DEC[status != "PASS" | rule == "OVERALL"])
cat(sprintf("\njoint-platform PCA decision: %s\n", overall))
if(overall == "STOP") quit(save="no", status=3)
