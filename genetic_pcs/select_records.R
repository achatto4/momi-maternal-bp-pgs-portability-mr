## select_records.R
##
## One genotype record per woman of the cleaned sample (her low-pass sequencing record if usable, otherwise her array
## record), the unrelated reference set, the relatedness weights, and the genotyping-technology check.
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(here, "helpers", "ids.R"))
args <- commandArgs(TRUE); arg <- function(f, d=NULL){ i <- which(args == f); if(length(i)) args[i + 1] else d }
COH <- arg("--cohort"); KEEPF <- arg("--keep"); DROPF <- arg("--drops"); LPS <- arg("--lp-psam"); GSS <- arg("--gsa-samples")
SSC <- arg("--sscore-dir"); PAIRS <- arg("--pairs"); DSELF <- arg("--dual-self"); OUT <- arg("--out")
SEED <- as.integer(arg("--seed", "20260925"))
for(v in c("COH", "KEEPF", "LPS", "GSS", "SSC", "PAIRS", "OUT")) if(is.null(get(v))) stop("missing argument for ", v)
dir.create(OUT, FALSE, TRUE)
T_ID <- 2^-1.5; T_2 <- 2^-3.5
say <- function(...) cat(sprintf("[%s] %s: ", format(Sys.time(), "%H:%M:%S"), COH), ..., "\n", sep="")

rd_samples <- function(f){
  if(grepl("\\.fam$", f)){ x <- fread(f, header=FALSE, colClasses="character"); d <- data.table(FID=x[[1]], IID=x[[2]]) }
  else { d <- rd_psam(f); d <- d[, intersect(c("FID", "IID"), names(d)), with=FALSE] }
  d[, nid := normid(IID)]; d }
LP <- rd_samples(LPS); GS <- rd_samples(GSS)
if(anyDuplicated(LP$nid) || anyDuplicated(GS$nid)) stop("a sample appears twice within one fileset")

KEEP <- normid(fread(KEEPF, header=FALSE, colClasses="character")[[1]])
DROP <- if(!is.null(DROPF) && file.exists(DROPF)) fread(DROPF, colClasses="character") else data.table(ID=character(0), platform=character(0))
if(nrow(DROP)){ setnames(DROP, 1:2, c("ID", "platform")); DROP[, ID := normid(ID)]
  DROP[, platform := fifelse(grepl("gsa", tolower(platform)), "gsa", fifelse(grepl("lp|low", tolower(platform)), "lpwgs_dosage", NA_character_))]
  if(anyNA(DROP$platform)) stop("unrecognised platform label in the drop list") }
dropped <- function(id, pl) paste(id, pl) %chin% DROP[, paste(ID, platform)]

W <- data.table(IID=KEEP[KEEP %chin% c(LP$nid, GS$nid)])
W[, `:=`(in_lp=IID %chin% LP$nid, in_gsa=IID %chin% GS$nid)]
W[, `:=`(has_lp=in_lp & !dropped(IID, "lpwgs_dosage"), has_gsa=in_gsa & !dropped(IID, "gsa"))]
W[, record := fifelse(has_lp, "low-pass", fifelse(has_gsa, "GSA", NA_character_))]
W[, extra_gsa := has_lp & has_gsa]
n_norecord <- W[is.na(record), .N]

probe <- function(pl){ f <- file.path(SSC, sprintf("PGS004603__%s__%s.sscore", COH, pl))
  if(file.exists(f)) normid(as.character(fread(f)[[1]])) else character(0) }
g_ids <- probe("gsa"); l_ids <- probe("lpwgs_dosage")
W[, technology := fifelse(IID %chin% g_ids & IID %chin% l_ids, "both", fifelse(IID %chin% g_ids, "GSA only",
                   fifelse(IID %chin% l_ids, "low-pass WGS only", "none")))]

P <- fread(PAIRS, colClasses=list(character=c("ID1", "ID2", "a", "b")))
P[, `:=`(a=normid(a), b=normid(b), KINSHIP=as.numeric(KINSHIP), class=as.character(class))]
n_identical_file <- P[class == "identical", .N]
P <- P[a %chin% W$IID & b %chin% W$IID & a != b & KINSHIP >= T_2, .(a, b, KINSHIP, class)]
n_identical_kept <- P[KINSHIP >= T_ID, .N]
cand <- W[!is.na(record), IID]
Ed <- P[, .(x=a, y=b)]
set.seed(SEED); rnd <- setNames(runif(length(cand)), cand); removed <- character(0)
while(nrow(Ed)){ dg <- table(c(Ed$x, Ed$y)); tied <- names(dg)[dg == max(dg)]
  pick <- tied[order(rnd[tied])][1]; removed <- c(removed, pick); Ed <- Ed[x != pick & y != pick] }
ref <- setdiff(cand, removed)
stopifnot(!nrow(P[a %chin% ref & b %chin% ref]))
dual <- W[extra_gsa == TRUE, IID]
ref_nd <- setdiff(ref, dual)
kap <- function(ids, refset){ z <- rbind(P[a %chin% ids & b %chin% refset, .(ID=a, k=KINSHIP)], P[b %chin% ids & a %chin% refset, .(ID=b, k=KINSHIP)])
  o <- data.table(ID=ids); o[z[, .(s=sum(2 * k)), by=ID], s := i.s, on="ID"]; o[is.na(s), s := 0]; pmin(1, o$s) }
W[, in_ref := IID %chin% ref]; W[, removed_related := IID %chin% removed]
W[, n_relatives := { d <- table(c(P$a, P$b)); v <- as.integer(d[IID]); v[is.na(v)] <- 0L; v }]
W[, kappa_main := fifelse(in_ref, 0, kap(IID, ref))]

selfk <- numeric(0)
if(!is.null(DSELF) && file.exists(DSELF)){ d <- fread(DSELF, colClasses=list(character="ID"))
  if(nrow(d)) selfk <- setNames(as.numeric(d$KINSHIP), normid(d$ID)) }
W[, self_kinship := 0.5]; W[IID %chin% names(selfk), self_kinship := as.numeric(selfk[IID])]
W[extra_gsa == TRUE, kappa_extra_main := pmin(1, kap(IID, ref) + fifelse(in_ref, 2 * self_kinship, 0))]
W[, kappa_ldo := kap(IID, ref_nd)]
W[, in_ref_nodual := IID %chin% ref_nd]

wr <- function(S, ids, f){ y <- S[nid %chin% ids, intersect(c("FID", "IID"), names(S)), with=FALSE]
  if("FID" %in% names(y)) setnames(y, "FID", "#FID") else setnames(y, "IID", "#IID"); fwrite(y, f, sep="\t"); nrow(y) }
n_sel_lp <- wr(LP, W[record == "low-pass", IID], file.path(OUT, "sel_lp.keep"))
n_sel_gs <- wr(GS, W[record == "GSA", IID], file.path(OUT, "sel_gsa.keep"))
n_ext_gs <- wr(GS, W[extra_gsa == TRUE, IID], file.path(OUT, "extra_gsa.keep"))
n_all_gs <- wr(GS, W[record == "GSA" | extra_gsa == TRUE, IID], file.path(OUT, "gsa_all.keep"))
stopifnot(n_sel_lp == W[record == "low-pass", .N], n_sel_gs == W[record == "GSA", .N])
fwrite(W, file.path(OUT, "women.tsv"), sep="\t")

S <- rbind(W[, .(n=.N), by=.(technology, record)][, level := "technology x record used"],
           W[, .(technology="all", record="all", n=.N, level="women")], fill=TRUE)
fwrite(S[order(level, technology, record)], file.path(OUT, "agg_selection.tsv"), sep="\t")
R <- data.table(cohort=COH, women=nrow(W), women_in_keep_list_total=length(KEEP), no_record=n_norecord,
                low_pass_record=n_sel_lp, gsa_record=n_sel_gs, dual_extra_gsa_records=n_ext_gs,
                lp_record_dropped_uses_gsa=W[in_lp & !has_lp & record == "GSA", .N],
                relationship_pairs_kept=nrow(P), first_degree_pairs=P[class == "first degree", .N],
                second_degree_pairs=P[class == "second degree", .N],
                identical_pairs_between_retained_women=n_identical_kept, identical_pairs_in_file_all_women=n_identical_file,
                reference=length(ref), reference_gsa_records=W[in_ref & record == "GSA", .N], removed_for_relatedness=length(removed),
                nonreference_with_reference_relative=W[!in_ref & kappa_main > 0, .N], reference_without_dual=length(ref_nd),
                technology_none=W[technology == "none", .N], seed=SEED)
fwrite(R, file.path(OUT, "agg_reference.tsv"), sep="\t")
say(sprintf("%d women: %d low-pass, %d GSA records, %d dual extras; reference %d (%d GSA); %d removed for relatedness; identical retained pairs %d",
            nrow(W), n_sel_lp, n_sel_gs, n_ext_gs, length(ref), R$reference_gsa_records, length(removed), n_identical_kept))
if(n_norecord > 0) stop(n_norecord, " women of the cleaned sample have no usable genotype record")
if(n_identical_kept > 0) stop(n_identical_kept, " identical-genotype pairs between retained women: the sample cleaning not honoured")
