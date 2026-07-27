#!/usr/bin/env Rscript
# ============================================================
# check_stale.R — which deliverables are out of date with their inputs?
#
# The manifest records, for every item that ran, an md5 fingerprint of each intermediate it
# READ (in_md5). This script recomputes those fingerprints against the files on disk now and
# reports any mismatch. That is the mechanism by which an upstream change ROLLS DOWN: rebuild
# B01 and analytic_mothers gets a new digest, so every deliverable that reads it is flagged
# STALE until it is rerun.
#
#   Rscript bin/check_stale.R              # report
#   Rscript bin/check_stale.R --quiet      # exit 1 if anything is stale, no table
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R"))
P <- momi_paths(PIPE)
QUIET <- "--quiet" %in% commandArgs(TRUE)

if(!file.exists(P$manifest)) stop("no manifest — run the build first: ", P$manifest)
## fill=TRUE as a safety net: SKIP rows are written by run_build.sh in bash while OK rows are
## written by R, so a schema change in one and not the other yields ragged rows. That has
## happened once (2026-07-19); the schema is now derived from a single MANIFEST_COLS variable,
## but tolerate it here rather than failing the whole check on one malformed line.
M <- suppressWarnings(fread(P$manifest, fill=TRUE))
if(!"in_md5" %in% names(M)){
  cat("manifest predates input fingerprinting — rerun the full build once to populate it.\n")
  quit(status=0)
}

## keep the LAST run of each id (a rerun supersedes an earlier attempt)
M[, ord := .I]
last <- M[order(id, ord)][, .SD[.N], by=id]

## DEFECT FIXED 2026-07-19: the staleness loop below skips any row whose status is not OK,
## which silently DROPPED such items from the report — the script would print "Everything is
## current with its inputs" while a deliverable was broken, and the only trace was the
## "items OK" count quietly falling. Absence of a stale flag must never be read as health.
##
## But FAIL and SKIP are NOT the same condition and must not be shown the same way:
##   FAIL — the script ran and errored. Something is wrong NOW. Loud, itemised.
##   SKIP — the deliverable is not written yet (B19-B32 etc). Expected, and stable for weeks.
## Listing 13 unwritten modules as "BROKEN" on every check trains you to ignore the section,
## which is precisely the failure this fix exists to prevent. SKIP gets one quiet line.
broke <- last[status == "FAIL"]
skipd <- last[status == "SKIP"]

## DEFECT FIXED 2026-07-20: a module is logged under TWO ids -- the bash driver writes SKIP
## rows keyed by BUILD id ("B24") while R writes OK rows keyed by DELIVERABLE id
## ("S12_controls"). Keeping the last row per id therefore preserved the old SKIP row
## forever, so B22/B23/B24 still read "not yet written" after they had run successfully.
## The two rows share a `script` path, which is the reliable join. Drop any SKIP whose
## script has since produced an OK row.
if(nrow(skipd) && "script" %in% names(M)){
  ran_ok <- unique(M[status=="OK" & nzchar(script), script])
  skipd  <- skipd[!(script %in% ran_ok)]
}

rows <- lapply(seq_len(nrow(last)), function(i){
  r <- last[i]
  if(r$status != "OK") return(NULL)
  now <- momi_fingerprint(r$inputs, P)
  data.table(id=r$id, script=r$script, ran=r$ts,
             stale = !identical(now, r$in_md5),
             recorded=r$in_md5, current=now)
})
S <- rbindlist(rows[!vapply(rows, is.null, TRUE)], fill=TRUE)

nst <- sum(S$stale)
## --quiet is used as a gate in scripts, so a BROKEN item must fail it just as a stale one does
if(QUIET){ quit(status=if(nst>0 || nrow(broke)>0) 1 else 0) }

cat(sprintf("\nmanifest: %s\nitems OK: %d | STALE: %d | FAILED: %d | not yet written: %d\n\n",
            P$manifest, nrow(S), nst, nrow(broke), nrow(skipd)))
if(nrow(broke)){
  cat("FAILED — the last attempt ERRORED. These items are absent from the staleness check\n",
      "below, so do not read a clean STALE report as meaning they are fine:\n", sep="")
  print(broke[, .(id, script, ran=ts, message=substr(message,1,60))], class=FALSE)
  cat("\n")
}
if(nrow(skipd))
  cat(sprintf("not yet written (SKIP): %s\n\n", paste(sort(skipd$id), collapse=" ")))
if(nst){
  cat("STALE — an input changed since these last ran:\n")
  print(S[stale==TRUE, .(id, script, ran)], class=FALSE)
  cat("\nfingerprint detail (recorded -> current):\n")
  for(i in which(S$stale)) cat(sprintf("  %-22s %s\n                       -> %s\n",
                                       S$id[i], S$recorded[i], S$current[i]))
  cat("\nRerun the build (or the affected ids) before citing any of the above.\n")
} else if(!nrow(broke)) cat("Everything is current with its inputs.\n") else
  cat("No STALE items — but see FAILED above before citing anything.\n")

## also warn about outputs that have gone missing from disk
outs <- unlist(strsplit(paste(last[status=="OK"]$outputs, collapse=";"), ";"))
outs <- trimws(outs[nzchar(outs)])
gone <- outs[!vapply(outs, function(o)
  any(file.exists(file.path(c(P$tables,P$figures,P$intermediates), o))), TRUE)]
if(length(gone)) cat(sprintf("\nWARNING: %d recorded output(s) not found on disk: %s\n",
                             length(gone), paste(utils::head(gone,8), collapse=", ")))
