## frequency_concordance.R
##
## Cross-technology allele-frequency concordance of the shared variants.
suppressMessages(library(data.table))
a <- commandArgs(TRUE); LPF <- a[1]; GSF <- a[2]; REN <- a[3]; OUT <- a[4]; COH <- a[5]
MAXD <- if(length(a) >= 6) as.numeric(a[6]) else 0.10
rd <- function(f){ x <- fread(f, colClasses=list(character=c("ID", "REF", "ALT"))); setnames(x, 1, sub("^#", "", names(x)[1]))
  fc <- intersect(c("ALT_FREQS", "ALT_FREQ"), names(x))[1]; x[, .(ID, REF, ALT, p=as.numeric(get(fc)), obs=as.numeric(OBS_CT))] }
lp <- rd(LPF); gs <- rd(GSF)
rn <- fread(REN, header=FALSE, colClasses="character"); setnames(rn, c("ID_gsa", "ID_lp"))
gs <- merge(gs, rn, by.x="ID", by.y="ID_gsa")
m <- merge(lp, gs, by.x="ID", by.y="ID_lp", suffixes=c("_lp", "_gsa"))
m[, same := REF_lp == REF_gsa & ALT_lp == ALT_gsa][, swap := REF_lp == ALT_gsa & ALT_lp == REF_gsa]
if(m[!(same | swap), .N]) stop("alleles not aligned for ", m[!(same | swap), .N], " matched variants")
m[, p_gsa_aligned := fifelse(same, p_gsa, 1 - p_gsa)][, d := p_lp - p_gsa_aligned]
keep <- m[abs(d) <= MAXD & is.finite(d)]

pb <- (m$p_lp * m$obs_lp + m$p_gsa_aligned * m$obs_gsa) / (m$obs_lp + m$obs_gsa)
exp_beyond <- 100 * mean(2 * pnorm(-MAXD / sqrt(pb * (1 - pb) * (1 / m$obs_lp + 1 / m$obs_gsa))), na.rm=TRUE)
fwrite(keep[, .(ID)], file.path(OUT, "concordant_lp.ids"), col.names=FALSE)
fwrite(data.table(cohort=COH, shared_aligned=nrow(m), excluded_abs_diff_gt=nrow(m) - nrow(keep), max_abs_diff_rule=MAXD,
                  pct_excluded=100 * (nrow(m) - nrow(keep)) / nrow(m), pct_expected_by_sampling=exp_beyond, r_frequency=cor(m$p_lp, m$p_gsa_aligned),
                  median_abs_diff=median(abs(m$d)), p99_abs_diff=unname(quantile(abs(m$d), 0.99)), mean_signed_diff_lp_minus_gsa=mean(m$d),
                  median_obs_ct_lp=median(m$obs_lp), median_obs_ct_gsa=median(m$obs_gsa)),
       file.path(OUT, "agg_freq_concordance.tsv"), sep="\t")
cat(sprintf("%s: %d shared; %d excluded (|dp| > %.2f); r = %.5f; median |dp| = %.4f\n", COH, nrow(m), nrow(m) - nrow(keep), MAXD,
            cor(m$p_lp, m$p_gsa_aligned), median(abs(m$d))))
