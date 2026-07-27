#!/usr/bin/env Rscript
# ============================================================
# assemble_dosage_pcs.R  — turn plink2 eigenvec output into the pipeline's pcs.rds
# (within-cohort PCs, for MR adjustment) and SF1_pca_points.tsv (pooled PCA, for the
# colored-dot Figure S1). Called by build_dosage_pca.sh. Sample-wide, no 1000G projection.
#
#   Rscript assemble_dosage_pcs.R WORK RESULTS NPC "coh1,coh2,..."
# ============================================================
suppressMessages(library(data.table))
a <- commandArgs(TRUE)
WORK <- a[1]; RESULTS <- a[2]; NPC <- as.integer(a[3]); cohs <- strsplit(a[4], ",")[[1]]
pcn <- paste0("PC", seq_len(NPC))
dir.create(file.path(RESULTS, "intermediates"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RESULTS, "tables"), showWarnings = FALSE, recursive = TRUE)

read_ev <- function(f) {
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  setnames(d, 1, sub("^#", "", names(d)[1]))          # '#FID'/'#IID' -> 'FID'/'IID'
  if (!"IID" %in% names(d) && names(d)[1] == "FID" && ncol(d) >= 2) setnames(d, 2, "IID")
  pcs <- grep("^PC[0-9]+$", names(d), value = TRUE)
  d <- d[, c("IID", pcs), with = FALSE]
  d[, IID := as.character(IID)]
  d
}

## ---- within-cohort PCs -> pcs.rds (what the MR/S17/S12 steps read) ----
rows <- list()
for (coh in cohs) {
  d <- read_ev(file.path(WORK, paste0(coh, "_within.eigenvec")))
  if (is.null(d)) { cat("  (no within eigenvec for", coh, ")\n"); next }
  d[, cohort := coh]
  rows[[coh]] <- d
}
if (!length(rows)) stop("no within-cohort eigenvecs found in ", WORK)
PCS <- rbindlist(rows, fill = TRUE)
setcolorder(PCS, c("IID", "cohort", intersect(pcn, names(PCS))))
saveRDS(PCS, file.path(RESULTS, "intermediates", "pcs.rds"))
cat(sprintf("pcs.rds written: %d mothers, %d cohorts, within-cohort PC1-%d\n",
            nrow(PCS), uniqueN(PCS$cohort), sum(pcn %in% names(PCS))))
print(PCS[, .(n = .N), by = cohort][order(cohort)])

## ---- pooled dosage PCA -> SF1_pca_points.tsv (colored-dot figure; sample-wide) ----
pool <- read_ev(file.path(WORK, "pooled.eigenvec"))
if (!is.null(pool)) {
  map <- unique(PCS[, .(IID, cohort)])
  pool <- merge(pool, map, by = "IID", all.x = TRUE)
  evf <- file.path(WORK, "pooled.eigenval")
  pct <- if (file.exists(evf)) { v <- scan(evf, quiet = TRUE); round(100 * v / sum(v), 2) } else rep(NA_real_, NPC)
  out <- pool[, .(cohort, PC1, PC2)]
  out[, `:=`(PC1_pct = pct[1], PC2_pct = pct[2])]
  fwrite(out, file.path(RESULTS, "tables", "SF1_pca_points.tsv"), sep = "\t")
  cat(sprintf("SF1_pca_points.tsv written: %d individuals (pooled dosage PCA; PC1 %.1f%%, PC2 %.1f%%)\n",
              nrow(out), pct[1], pct[2]))
} else {
  cat("no pooled PCA eigenvec — SF1 keeps its centroid fallback\n")
}
