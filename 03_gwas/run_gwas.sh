#!/usr/bin/env bash
# ============================================================
# run_gwas.sh  (Stage 03)
# Internal GWAS on the analytic mothers: SBP/DBP (linear) + PTB (logistic),
# adjusting for PW_AGE + genetic PCs + site dummies.
#
#   1) make_covar.R -> covar file (PW_AGE + PCs + SITE dummies)
#   2) plink2 --glm per trait
#
# Outputs: <out-dir>/results_<TRAIT>/gwas_<TRAIT>.PHENO.glm.{linear,logistic.hybrid}
#
# Usage:
#   run_gwas.sh --merged PREFIX --keep mothers.keep --pheno-dir DIR \
#               --covar COVTABLE --eigenvec EVEC --out-dir DIR \
#               [--bp-window latest] [--n-pcs 10] [--plink2 plink2] [--rscript Rscript]
#               [--threads 8] [--mem-mb 60000]
# ============================================================
set -euo pipefail
MERGED=""; KEEP=""; PHDIR=""; COV=""; EVEC=""; OUTDIR=""
WIN=latest; NPC=10; PLINK2=plink2; RSCRIPT=Rscript; THREADS=8; MEMMB=60000
HERE="$(cd "$(dirname "$0")" && pwd)"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --merged) MERGED="$2"; shift 2;;
    --keep) KEEP="$2"; shift 2;;
    --pheno-dir) PHDIR="$2"; shift 2;;
    --covar) COV="$2"; shift 2;;
    --eigenvec) EVEC="$2"; shift 2;;
    --out-dir) OUTDIR="$2"; shift 2;;
    --bp-window) WIN="$2"; shift 2;;
    --n-pcs) NPC="$2"; shift 2;;
    --plink2) PLINK2="$2"; shift 2;;
    --rscript) RSCRIPT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mem-mb) MEMMB="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
for v in MERGED KEEP PHDIR COV EVEC OUTDIR; do [[ -n "${!v}" ]] || { echo "missing --$v"; exit 1; }; done
command -v "$PLINK2" >/dev/null || { echo "plink2 not on PATH (module load?)" >&2; exit 3; }
mkdir -p "$OUTDIR"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

# accept a bed (hard calls) or a pgen (dosages); --glm reads dosages from a pgen directly
GIN=(--bfile "$MERGED"); [[ -f "$MERGED.pgen" ]] && GIN=(--pfile "$MERGED")
[[ -f "$MERGED.bed" || -f "$MERGED.pgen" ]] || { echo "no .bed or .pgen at prefix: $MERGED" >&2; exit 2; }
log "genotype input: ${GIN[*]}"

COVFILE="$OUTDIR/covar4plink.txt"
log "building covariate file (PW_AGE + $NPC PCs + site dummies)"
"$RSCRIPT" "$HERE/make_covar.R" --covar "$COV" --eigenvec "$EVEC" --out "$COVFILE" --n-pcs "$NPC"

run_trait(){  # $1=trait label  $2=pheno file  $3=extra glm flags
  local trait="$1" pheno="$2" extra="$3"
  local rdir="$OUTDIR/results_${trait}"; mkdir -p "$rdir"
  [[ -f "$pheno" ]] || { log "SKIP $trait: pheno not found ($pheno)"; return; }
  log "GWAS $trait"
  "$PLINK2" "${GIN[@]}" --keep "$KEEP" \
    --pheno "$pheno" --pheno-name PHENO \
    --covar "$COVFILE" --covar-variance-standardize \
    --glm hide-covar omit-ref $extra \
    --threads "$THREADS" --memory "$MEMMB" \
    --out "$rdir/gwas_${trait}" > "$rdir/gwas_${trait}.run.log" 2>&1 \
    || { log "WARN: $trait glm returned non-zero; see $rdir/gwas_${trait}.run.log"; }
  ls "$rdir"/gwas_${trait}.PHENO.glm.* 2>/dev/null | while read f; do log "  -> $f ($(wc -l < "$f") lines)"; done
}

# SBP / DBP : quantitative -> linear ;  PTB : 1/2 case-control -> logistic (firth fallback)
run_trait "SBP"     "$PHDIR/mothers_SBP_${WIN}_final.pheno"  ""
run_trait "DBP"     "$PHDIR/mothers_DBP_${WIN}_final.pheno"  ""
run_trait "PTB_NEW" "$PHDIR/mothers_PTB_NEW_final.pheno"     "firth-fallback"

log "GWAS stage done -> $OUTDIR/results_*"
