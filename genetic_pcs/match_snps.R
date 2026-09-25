## match_snps.R
##
## Position-matched, allele-aligned variants shared by the two technologies (strand-ambiguous and mismatched removed).
suppressMessages(library(data.table))
a <- commandArgs(TRUE); LPV <- a[1]; GSV <- a[2]; OUT <- a[3]
rd <- function(f){ x <- fread(cmd=sprintf("grep -v '^##' '%s'", f), colClasses="character")
  setnames(x, 1, "CHROM"); x[, CHROM := sub("^chr", "", CHROM, ignore.case=TRUE)]
  x[, key := paste(CHROM, POS, sep=":")]; x[, .(key, ID, REF, ALT)] }
lp <- rd(LPV); gs <- rd(GSV); n_lp_raw <- nrow(lp); n_gs_raw <- nrow(gs)
okal <- function(d) d[nchar(REF) == 1 & nchar(ALT) == 1 & REF %chin% c("A","C","G","T") & ALT %chin% c("A","C","G","T") & REF != ALT]
lp <- okal(lp); gs <- okal(gs)
lp <- lp[!(duplicated(key) | duplicated(key, fromLast=TRUE))]; gs <- gs[!(duplicated(key) | duplicated(key, fromLast=TRUE))]
m <- merge(lp, gs, by="key", suffixes=c("_lp", "_gsa"))
amb <- function(r, x) paste0(pmin(r, x), pmax(r, x)) %chin% c("AT", "CG")
comp <- c(A="T", C="G", G="C", T="A")
m[, ambiguous := amb(REF_lp, ALT_lp)]
m[, same := REF_lp == REF_gsa & ALT_lp == ALT_gsa]; m[, swapped := REF_lp == ALT_gsa & ALT_lp == REF_gsa]
m[, mismatch := !(same | swapped)]
m[, complement_matchable := mismatch & ((comp[REF_lp] == REF_gsa & comp[ALT_lp] == ALT_gsa) | (comp[REF_lp] == ALT_gsa & comp[ALT_lp] == REF_gsa))]
k <- m[!ambiguous & !mismatch]
fwrite(k[, .(ID_lp)], file.path(OUT, "shared_lp.ids"), col.names=FALSE)
fwrite(k[, .(ID_gsa)], file.path(OUT, "shared_gsa.ids"), col.names=FALSE)
fwrite(k[, .(ID_gsa, ID_lp)], file.path(OUT, "gsa_rename.txt"), col.names=FALSE, sep="\t")
fwrite(k[, .(ID_lp, REF_lp)], file.path(OUT, "lp_ref.txt"), col.names=FALSE, sep="\t")
cnt <- data.table(step=c("lp_qc_snps", "gsa_qc_snps", "lp_qc_snps_unique_position", "gsa_qc_snps_unique_position", "position_match", "strand_ambiguous_excluded", "allele_mismatch_excluded",
                         "of_which_complement_matchable", "shared_same_ref_alt", "shared_swapped_ref_alt", "shared_aligned"),
                  n=c(n_lp_raw, n_gs_raw, nrow(lp), nrow(gs), nrow(m), sum(m$ambiguous), sum(m$mismatch & !m$ambiguous), sum(m$complement_matchable & !m$ambiguous),
                      sum(k$same), sum(k$swapped), nrow(k)))
fwrite(cnt, file.path(OUT, "match_counts.tsv"), sep="\t")
print(cnt)
