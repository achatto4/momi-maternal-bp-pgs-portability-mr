#!/usr/bin/env bash
# ============================================================================
# build_dosage_pca.sh — build pcs.rds from the low-pass dosage filesets, sample-wide
# (no 1000G projection). Produces:
#   * within-cohort PCs  -> intermediates/pcs.rds   (used by MR: B22, and S17, S12, T5, S6)
#   * pooled dosage PCA  -> tables/SF1_pca_points.tsv (colored-dot Figure S1)
#
# This is the "sound fix" from docs/paths_reference.md: single-platform (dosage-only) so no
# GSA/lpWGS batch effect, and within-cohort PCAs so each cohort's structure is described in
# its own coordinates (what the MR adjusts for).
#
#   module load plink/2.00a4.6 conda_R/4.3
#   export MOMI_OUTROOT=/dcs10/chatterj/data/achattop/MOMI/Genomics/pipeline_out
#   bash 06_phase2/build_dosage_pca.sh --results results/momi_output
#
# Options:  --outroot DIR   (or export MOMI_OUTROOT)   --results DIR   --npc N (default 10)
# ============================================================================
set -uo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"     # bp_ptb_pipeline
OUTROOT="${MOMI_OUTROOT:-}"
RESULTS="results/momi_output"
NPC=10
PLINK2="${PLINK2:-plink2}"
COHORTS=(AMANHI-Bangladesh AMANHI-Pakistan AMANHI-Pemba GAPPS-Bangladesh GAPPS-Zambia)

while [[ $# -gt 0 ]]; do case "$1" in
  --outroot) OUTROOT="$2"; shift 2;;
  --results) RESULTS="$2"; shift 2;;
  --npc)     NPC="$2"; shift 2;;
  *) echo "unknown arg: $1"; exit 2;;
esac; done

[[ -n "$OUTROOT" && -d "$OUTROOT" ]] || { echo "Set --outroot (or MOMI_OUTROOT) to the per-cohort clean filesets dir."; exit 1; }
command -v "$PLINK2" >/dev/null 2>&1 || { echo "plink2 not on PATH — run: module load plink/2.00a4.6"; exit 1; }

WORK="$RESULTS/qc/dosage_pca"; mkdir -p "$WORK" "$RESULTS/intermediates" "$RESULTS/tables"
say(){ printf '\n== %s\n' "$*"; }

# per-cohort within-cohort PCA -------------------------------------------------
for coh in "${COHORTS[@]}"; do
  pre="$OUTROOT/$coh/lpwgs_dosage_clean"
  if   [[ -f "$pre.pgen" ]]; then gin="--pfile"
  elif [[ -f "$pre.bed"  ]]; then gin="--bfile"
  else echo "  skip $coh (no lpwgs_dosage_clean fileset at $pre)"; continue; fi
  say "within-cohort PCA: $coh"
  "$PLINK2" $gin "$pre" --maf 0.01 --geno 0.05 --snps-only \
      --indep-pairwise 200 50 0.2 --out "$WORK/${coh}_prune" >"$WORK/${coh}_prune.log" 2>&1 || true
  if [[ -f "$WORK/${coh}_prune.prune.in" ]]; then
    "$PLINK2" $gin "$pre" --extract "$WORK/${coh}_prune.prune.in" \
        --pca "$NPC" --out "$WORK/${coh}_within" >"$WORK/${coh}_within.log" 2>&1 || true
  fi
  [[ -f "$WORK/${coh}_within.eigenvec" ]] && echo "  ok: $WORK/${coh}_within.eigenvec" \
      || echo "  !! $coh within-cohort PCA failed (see log)"
done

# pooled dosage PCA (for the SF1 figure) — best-effort ------------------------
say "pooled dosage PCA (merge of all cohorts)"
ML="$WORK/mergelist.txt"; : > "$ML"
for coh in "${COHORTS[@]}"; do
  pre="$OUTROOT/$coh/lpwgs_dosage_clean"
  [[ -f "$pre.pgen" ]] && echo "$pre" >> "$ML"
done
if [[ -s "$ML" ]]; then
  "$PLINK2" --pmerge-list "$ML" pfile --make-pgen --out "$WORK/merged_dosage" >"$WORK/merge.log" 2>&1 \
    && "$PLINK2" --pfile "$WORK/merged_dosage" --maf 0.01 --geno 0.05 --snps-only \
         --indep-pairwise 200 50 0.2 --out "$WORK/pool_prune" >"$WORK/pool_prune.log" 2>&1 \
    && "$PLINK2" --pfile "$WORK/merged_dosage" --extract "$WORK/pool_prune.prune.in" \
         --pca "$NPC" --out "$WORK/pooled" >"$WORK/pooled.log" 2>&1
  [[ -f "$WORK/pooled.eigenvec" ]] && echo "  ok: pooled PCA" \
      || echo "  !! pooled PCA failed (SF1 will use centroid fallback) — see $WORK/*.log"
else
  echo "  no dosage pgens to merge; skipping pooled PCA"
fi

# assemble into pcs.rds + SF1_pca_points.tsv ----------------------------------
say "assembling pcs.rds and SF1_pca_points.tsv"
CSV=$(IFS=,; echo "${COHORTS[*]}")
Rscript "$PIPE/06_phase2/assemble_dosage_pcs.R" "$WORK" "$RESULTS" "$NPC" "$CSV"

say "done — now rerun the MR/PCA deliverables, e.g.:  bash run_build.sh --results $RESULTS --epi \"\$EPI\" --sscore-dir \"\$SSC\""
