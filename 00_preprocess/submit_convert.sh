#!/usr/bin/env bash
# ============================================================
# submit_convert.sh
# One idempotent command per cohort x platform:
#   1) build sample map (make_sample_map.sh) on the login node (fast, awk)
#   2) submit vcf_to_plink.sh as a SLURM job (heavy conversion)
#
# Guards against accidental double-submission (checks for a running job with the
# same name, and for an existing finished output unless --force).
#
# File conventions inside a cohort dir:
#   gsa   : gsa_reimputed.vcf.gz   + gsa_idmap.csv
#   lpwgs : lpwgs_imputed.vcf.gz   + lpwgs_idmap.csv
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

usage() {
  cat <<USAGE
Usage: $0 --cohort-dir DIR --platform gsa|lpwgs --pid-prefix STR --epi FILE --out DIR
          [--plink-module plink/2.00a4.6] [--bcftools-module bcftools]
          [--r2-min 0.3] [--maf 0.01]
          [--cpus 8] [--mem 72G] [--time 8:00:00] [--force]
USAGE
  exit 1
}

COHORT_DIR=""; PLAT=""; PREFIX=""; EPI=""; OUT=""; EXTRACT=""
PLINK_MOD="plink/2.00a4.6"; BCF_MOD="bcftools"
R2MIN=0.8; MAF=0; CPUS=8; MEM="72G"; TIME="8:00:00"; FORCE=0   # MAF off; quality via R2; union via --extract
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cohort-dir) COHORT_DIR="$2"; shift 2;;
    --platform) PLAT="$2"; shift 2;;
    --pid-prefix) PREFIX="$2"; shift 2;;
    --epi) EPI="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --extract) EXTRACT="$2"; shift 2;;
    --plink-module) PLINK_MOD="$2"; shift 2;;
    --bcftools-module) BCF_MOD="$2"; shift 2;;
    --r2-min) R2MIN="$2"; shift 2;;
    --maf) MAF="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    --mem) MEM="$2"; shift 2;;
    --time) TIME="$2"; shift 2;;
    --force) FORCE=1; shift;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done
[[ -n "$COHORT_DIR" && -n "$PLAT" && -n "$PREFIX" && -n "$EPI" && -n "$OUT" ]] || usage

case "$PLAT" in
  gsa)   VCF="$COHORT_DIR/gsa_reimputed.vcf.gz";  IDMAP="$COHORT_DIR/gsa_idmap.csv";;
  lpwgs) VCF="$COHORT_DIR/lpwgs_imputed.vcf.gz";  IDMAP="$COHORT_DIR/lpwgs_idmap.csv";;
  *) echo "platform must be gsa|lpwgs" >&2; exit 1;;
esac
[[ -f "$VCF"   ]] || { echo "VCF not found:   $VCF"   >&2; exit 2; }
[[ -f "$IDMAP" ]] || { echo "idmap not found: $IDMAP" >&2; exit 2; }
mkdir -p "$OUT"

COH="$(basename "$COHORT_DIR")"
JOB="cv_${COH}_${PLAT}"
MAP="$OUT/${PLAT}_map"
FINAL="$OUT/${PLAT}_clean"

# ---- guard: already finished? ----
if [[ -f "$FINAL.bed" && $FORCE -eq 0 ]]; then
  echo "[skip] $FINAL.bed already exists (use --force to redo). Samples=$(wc -l < "$FINAL.fam" 2>/dev/null || echo '?')"
  exit 0
fi
# ---- guard: already queued/running? ----
if squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -qx "$JOB"; then
  echo "[skip] a job named '$JOB' is already in the queue. Cancel it first or wait."
  exit 0
fi

echo "[map] building sample map for $COH / $PLAT (prefix=$PREFIX)"
bash "$HERE/make_sample_map.sh" --idmap "$IDMAP" --epi "$EPI" --pid-prefix "$PREFIX" --out "$MAP"
NKEEP=$(wc -l < "$MAP.keep")
[[ "$NKEEP" -gt 0 ]] || { echo "ERROR: 0 mothers mapped — check --pid-prefix '$PREFIX'." >&2; exit 3; }

MEMMB=$(( ${MEM%G} * 1000 - 4000 ))   # leave headroom under the SLURM --mem cap
echo "[submit] $JOB  ($NKEEP mothers, cpus=$CPUS mem=$MEM time=$TIME)"
sbatch --job-name="$JOB" --cpus-per-task="$CPUS" --mem="$MEM" --time="$TIME" \
  --output="$OUT/convert_${PLAT}_%j.out" \
  --wrap "module load $PLINK_MOD; module load $BCF_MOD 2>/dev/null || true; \
bash '$HERE/vcf_to_plink.sh' \
  --vcf '$VCF' --platform '$PLAT' \
  --keep '$MAP.keep' --rename '$MAP.rename' \
  --out '$FINAL' ${EXTRACT:+--extract '$EXTRACT'} \
  --r2-min '$R2MIN' --maf '$MAF' --threads '$CPUS' --mem-mb '$MEMMB'"

echo "[ok] submitted. Track: squeue -u $USER  |  log: $OUT/convert_${PLAT}_<jobid>.out"
