#!/usr/bin/env bash
# ============================================================
# run_pca.sh
# Stage 02: genetic PCs on the analytic mothers.
#   1) LD-prune (plink1.9) on a MAF/geno/HWE-clean SNP set
#   2) PCA (plink2 --pca) on the pruned set
# Output: <out-dir>/pca_analysis/pca_results.{eigenvec,eigenval} + pruned lists
#
# These PCs are used as GWAS/MR/PRS covariates to control population structure
# (our merged set spans SAS + AFR cohorts, so the top PCs capture ancestry).
#
# Usage:
#   run_pca.sh --merged PREFIX --keep mothers.keep --out-dir DIR
#              [--n-pcs 10] [--maf 0.05] [--geno 0.02] [--hwe 1e-6]
#              [--prune "50 5 0.2"] [--plink1 plink] [--plink2 plink2]
#              [--threads 4] [--mem-mb 28000]
# ============================================================
set -euo pipefail

MERGED=""; KEEP=""; OUTDIR=""; NPC=10
MAF=0.05; GENO=0.02; HWE=1e-6; PRUNE="50 5 0.2"
PLINK1=plink; PLINK2=plink2; THREADS=4; MEMMB=28000
while [[ $# -gt 0 ]]; do
  case "$1" in
    --merged) MERGED="$2"; shift 2;;
    --keep) KEEP="$2"; shift 2;;
    --out-dir) OUTDIR="$2"; shift 2;;
    --n-pcs) NPC="$2"; shift 2;;
    --maf) MAF="$2"; shift 2;;
    --geno) GENO="$2"; shift 2;;
    --hwe) HWE="$2"; shift 2;;
    --prune) PRUNE="$2"; shift 2;;
    --plink1) PLINK1="$2"; shift 2;;
    --plink2) PLINK2="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mem-mb) MEMMB="$2"; shift 2;;
    -h|--help) echo "see header"; exit 1;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$MERGED" && -n "$KEEP" && -n "$OUTDIR" ]] || { echo "need --merged --keep --out-dir" >&2; exit 1; }
[[ -f "$MERGED.bed" || -f "$MERGED.pgen" ]] || { echo "missing $MERGED.bed / $MERGED.pgen" >&2; exit 2; }
[[ -f "$KEEP" ]] || { echo "missing $KEEP" >&2; exit 2; }
command -v "$PLINK2" >/dev/null || { echo "plink2 not on PATH (module load?)" >&2; exit 3; }
# accept a bed (hard calls) or pgen (dosages) prefix
GIN=(--bfile "$MERGED"); [[ -f "$MERGED.pgen" ]] && GIN=(--pfile "$MERGED")
# Note: plink1.9 and plink/2.x share the 'plink' module name on JHPCE and can't
# co-load, so this stage uses plink2 for BOTH pruning and PCA.

PCADIR="$OUTDIR/pca_analysis"; mkdir -p "$PCADIR"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
log "PCA on $(wc -l < "$KEEP") mothers from $MERGED  (n_pcs=$NPC)"

# 1) LD-prune (plink2)
log "LD-pruning: --maf $MAF --geno $GENO --hwe $HWE --indep-pairwise $PRUNE"
"$PLINK2" "${GIN[@]}" --keep "$KEEP" \
  --maf "$MAF" --geno "$GENO" --hwe "$HWE" \
  --indep-pairwise $PRUNE \
  --threads "$THREADS" --memory "$MEMMB" \
  --out "$PCADIR/pruned_for_pca" > "$PCADIR/prune.log" 2>&1
NPRUNE=$(wc -l < "$PCADIR/pruned_for_pca.prune.in")
log "pruned to $NPRUNE independent SNPs"

# 2) PCA
log "computing $NPC PCs"
"$PLINK2" "${GIN[@]}" --keep "$KEEP" \
  --extract "$PCADIR/pruned_for_pca.prune.in" \
  --pca "$NPC" \
  --threads "$THREADS" --memory "$MEMMB" \
  --out "$PCADIR/pca_results" > "$PCADIR/pca.log" 2>&1

[[ -f "$PCADIR/pca_results.eigenvec" ]] || { echo "ERROR: PCA failed; see $PCADIR/pca.log" >&2; exit 4; }
log "DONE: $(wc -l < "$PCADIR/pca_results.eigenvec") rows -> $PCADIR/pca_results.eigenvec"
log "eigenvalues (top): $(head -5 "$PCADIR/pca_results.eigenval" | tr '\n' ' ')"
