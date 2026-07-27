#!/usr/bin/env bash
# ============================================================
# score_panel.sh
# Score every PGS in the transferability panel on every cohort x platform fileset
# (gsa/lpwgs hard-call beds and the lpwgs_dosage pgen, dosages used automatically).
# Output: <out-dir>/<PGS_id>__<COHORT>__<PLATFORM>.sscore
# (used by transferability_analysis.R to build the PGS-ancestry x cohort-ancestry
#  prediction-accuracy matrix).
#
# Usage:
#   score_panel.sh --panel transfer_panel.tsv --scores-dir DIR --outroot OUTROOT --out-dir DIR
#                  [--plink2 plink2]
#   (--scores-dir holds <PGS_id>.txt.gz downloaded by run_transferability.sh)
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PANEL=""; SCORES=""; OUTROOT=""; OUTDIR=""; PLINK2=plink2; PLATFORMS_ARG=""; COHORTS_ARG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --panel) PANEL="$2"; shift 2;;
    --scores-dir) SCORES="$2"; shift 2;;
    --outroot) OUTROOT="$2"; shift 2;;
    --out-dir) OUTDIR="$2"; shift 2;;
    --platforms) PLATFORMS_ARG="$2"; shift 2;;   # default: gsa lpwgs lpwgs_dosage
    --cohorts) COHORTS_ARG="$2"; shift 2;;       # ADDED 2026-07-20 -- see note below
    --plink2) PLINK2="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
for v in PANEL SCORES OUTROOT OUTDIR; do [[ -n "${!v}" ]] || { echo "missing --$v"; exit 1; }; done
mkdir -p "$OUTDIR"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
# COHORTS: overridable as of 2026-07-20 so the work can be split across SLURM ARRAY TASKS.
# Without an axis to parallelise over, filling the 39 missing score x cohort x platform
# combos ran strictly serially -- several of them against 14M-variant lpWGS dosage filesets.
# One array task per (score, cohort) turns that into 15 concurrent jobs. Safe to run
# concurrently because every task writes a DISJOINT set of files (<pid>__<cohort>__<plat>)
# and the "have this .sscore already" check below makes each task idempotent.
COHORTS=(AMANHI-Bangladesh AMANHI-Pakistan AMANHI-Pemba GAPPS-Bangladesh GAPPS-Zambia)
if [[ -n "$COHORTS_ARG" ]]; then read -r -a COHORTS <<< "$COHORTS_ARG"; fi
# platforms: gsa (hard call), lpwgs (hard call), lpwgs_dosage (dosages). compute_prs
# auto-detects bed-vs-pgen; any platform whose fileset is absent for a cohort is skipped.
if [[ -n "$PLATFORMS_ARG" ]]; then read -r -a PLATFORMS <<< "$PLATFORMS_ARG"; else PLATFORMS=(gsa lpwgs lpwgs_dosage); fi

# read panel (skip header)
tail -n +2 "$PANEL" | while IFS=$'\t' read -r pid trait anc src; do
  [[ -n "$pid" ]] || continue
  pgsfile="$SCORES/${pid}.txt.gz"
  [[ -f "$pgsfile" ]] || { log "MISSING score file $pgsfile — skipping $pid"; continue; }
  for c in "${COHORTS[@]}"; do
    for plat in "${PLATFORMS[@]}"; do
      bed="$OUTROOT/$c/${plat}_clean"
      [[ -f "$bed.bed" || -f "$bed.pgen" ]] || continue
      out="$OUTDIR/${pid}__${c}__${plat}"
      [[ -f "$out.sscore" ]] && { log "have $out.sscore"; continue; }
      log "score $pid ($anc/$trait) on $c/$plat"
      bash "$HERE/compute_prs.sh" --pgs "$pgsfile" --bfile "$bed" --out "$out" \
        --plink2 "$PLINK2" --threads 4 --mem-mb 16000 \
        || log "WARN: failed $pid on $c/$plat (see $out.score.log)"
    done
  done
done
log "DONE. sscores: $(ls -1 "$OUTDIR"/*.sscore 2>/dev/null | wc -l)"
