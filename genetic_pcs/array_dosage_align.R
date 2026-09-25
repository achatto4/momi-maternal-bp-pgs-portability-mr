## array_dosage_align.R
##
## Alignment of the array dosages to the low-pass alleles, checked against the array hard calls.
suppressMessages(library(data.table))
a <- commandArgs(TRUE); AV <- a[1]; LV <- a[2]; OUT <- a[3]
rd <- function(f){ x <- fread(cmd=sprintf("grep -v '^##' '%s'", f), colClasses="character"); setnames(x, 1, "CHROM")
  x[, key := paste(sub("^chr", "", CHROM, ignore.case=TRUE), POS, sep=":")]; x[, .(key, ID, REF, ALT)] }
av <- rd(AV); lv <- rd(LV)
av <- av[!(duplicated(key) | duplicated(key, fromLast=TRUE))]
m <- merge(av, lv, by="key", suffixes=c("_arr", "_lp"))
m <- m[(REF_arr == REF_lp & ALT_arr == ALT_lp) | (REF_arr == ALT_lp & ALT_arr == REF_lp)]
fwrite(m[, .(ID_arr, ID_lp)], file.path(OUT, "arr_rename.txt"), col.names=FALSE, sep="\t")
fwrite(m[, .(ID_lp)], file.path(OUT, "arr_extract_lp.ids"), col.names=FALSE)
fwrite(m[, .(ID_lp, REF_lp)], file.path(OUT, "arr_lp_ref.txt"), col.names=FALSE, sep="\t")
fwrite(data.table(pruned_lowpass=nrow(lv), imported_array=nrow(av), matched_alleles=nrow(m), not_recovered=nrow(lv) - nrow(m)),
       file.path(OUT, "agg_array_dosage_import.tsv"), sep="\t")
cat(sprintf("array dosage import: %d of %d pruned variants recovered from the raw VCF\n", nrow(m), nrow(lv)))
