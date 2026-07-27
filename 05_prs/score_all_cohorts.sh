#!/usr/bin/env bash
# ============================================================
# score_all_cohorts.sh
# Compute the 4 BP PGS (SBP/DBP x EUR/SAS) on EACH cohort x platform's own bed.
# No genotype merging — each cohort is scored on its own variants. Output:
#   <out-dir>/<COHORT>__<PLATFORM>__<TAG>.sscore   (TAG in SBP_EUR,DBP_EUR,SBP_SAS,DBP_SAS)
# These per-sample scores are pooled later (pooled_prs_analysis.R).
#
# Usage:
#   score_all_cohorts.sh --outroot DIR --out-dir DIR [--sumdir DIR] [--plink2 plink2]
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUTROOT=""; OUTDIR=""; PLINK2=plink2
SUM="/dcs04/nilanjan/data/Anagh/MOMI/sumstat/BP/PRS_BP"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outroot) OUTROOT="$2"; shift 2;;
    --out-dir) OUTDIR="$2"; shift 2;;
    --sumdir) SUM="$2"; shift 2;;
    --plink2) PLINK2="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTROOT" && -n "$OUTDIR" ]] || { echo "need --outroot --out-dir"; exit 1; }
mkdir -p "$OUTDIR"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

COHORTS=(AMANHI-Bangladesh AMANHI-Pakistan AMANHI-Pemba GAPPS-Bangladesh GAPPS-Zambia)
PLATFORMS=(gsa lpwgs)
declare -A PGS=(
  [SBP_EUR]="$SUM/SBP/PGS004603_hmPOS_GRCh38.txt.gz"
  [DBP_EUR]="$SUM/DBP/PGS004604_hmPOS_GRCh38.txt.gz"
  [SBP_SAS]="$SUM/SBP/PGS004830_hmPOS_GRCh38.txt.gz"
  [DBP_SAS]="$SUM/DBP/PGS004758_hmPOS_GRCh38.txt.gz"
)

for c in "${COHORTS[@]}"; do
  for plat in "${PLATFORMS[@]}"; do
    bed="$OUTROOT/$c/${plat}_clean"
    [[ -f "$bed.bed" ]] || { log "skip $c/$plat (no bed)"; continue; }
    # NOTE: no --keep — cleaned beds are already mothers-only, renamed to PARTICIPANT_ID
    for tag in SBP_EUR DBP_EUR SBP_SAS DBP_SAS; do
      out="$OUTDIR/${c}__${plat}__${tag}"
      [[ -f "$out.sscore" ]] && { log "have $out.sscore"; continue; }
      log "score $c/$plat $tag"
      bash "$HERE/compute_prs.sh" --pgs "${PGS[$tag]}" --bfile "$bed" \
        --out "$out" --plink2 "$PLINK2" --threads 4 --mem-mb 16000 \
        || log "WARN: scoring failed for $c/$plat/$tag (see $out.score.log)"
    done
  done
done
log "DONE scoring. sscores in $OUTDIR"
ls -1 "$OUTDIR"/*.sscore 2>/dev/null | wc -l
