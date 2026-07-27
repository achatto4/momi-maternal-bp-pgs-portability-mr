#!/usr/bin/env bash
# ============================================================================
# build_pooled_pca.sh — robust pooled dosage PCA for the colored-dot Figure S1.
# The within-cohort PCs for the MR are already built by build_dosage_pca.sh; this only
# (re)builds the POOLED PCA, which the plain --pmerge-list can trip over when per-cohort
# filesets disagree on variant IDs/alleles. Strategy: standardise variant IDs to
# chr:pos:ref:alt, drop duplicates, intersect to the common variant set, then merge.
#
#   module load plink/2.00a4.6 conda_R/4.3
#   export MOMI_OUTROOT=/dcs10/chatterj/data/achattop/MOMI/Genomics/pipeline_out
#   bash 06_phase2/build_pooled_pca.sh --results results/momi_output
#   # then re-assemble so SF1_pca_points.tsv gets written:
#   Rscript 06_phase2/assemble_dosage_pcs.R results/momi_output/qc/dosage_pca results/momi_output 10 \
#       "AMANHI-Bangladesh,AMANHI-Pakistan,AMANHI-Pemba,GAPPS-Bangladesh,GAPPS-Zambia"
# ============================================================================
set -uo pipefail
OUTROOT="${MOMI_OUTROOT:-}"; RESULTS="results/momi_output"; NPC=10
PLINK2="${PLINK2:-plink2}"; PLINK1="${PLINK1:-plink}"   # plink2 for pgen; plink 1.9 for sample-merge
COHORTS=(AMANHI-Bangladesh AMANHI-Pakistan AMANHI-Pemba GAPPS-Bangladesh GAPPS-Zambia)
while [[ $# -gt 0 ]]; do case "$1" in
  --outroot) OUTROOT="$2"; shift 2;; --results) RESULTS="$2"; shift 2;; --npc) NPC="$2"; shift 2;;
  *) echo "unknown arg: $1"; exit 2;; esac; done
[[ -n "$OUTROOT" && -d "$OUTROOT" ]] || { echo "set --outroot / MOMI_OUTROOT"; exit 1; }
command -v "$PLINK2" >/dev/null || { echo "plink2 not on PATH (module load plink/2.00a4.6)"; exit 1; }
WORK="$RESULTS/qc/dosage_pca"; mkdir -p "$WORK"
say(){ printf '\n== %s\n' "$*"; }

# 1. standardise IDs + drop duplicates, per cohort ----------------------------
present=()
for coh in "${COHORTS[@]}"; do
  pre="$OUTROOT/$coh/lpwgs_dosage_clean"
  [[ -f "$pre.pgen" ]] || { echo "  skip $coh (no dosage pgen)"; continue; }
  say "standardise IDs: $coh"
  "$PLINK2" --pfile "$pre" \
      --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 200 missing \
      --rm-dup exclude-all --snps-only just-acgt \
      --make-pgen --out "$WORK/std_$coh" >"$WORK/std_$coh.log" 2>&1 \
    && present+=("$coh") || echo "  !! standardise failed for $coh (see log)"
done
[[ ${#present[@]} -ge 2 ]] || { echo "fewer than two cohorts standardised; aborting pooled PCA"; exit 1; }

# 2. common-variant intersection ---------------------------------------------
say "intersecting variant IDs across ${#present[@]} cohorts"
common="$WORK/common_ids.txt"; tmp="$WORK/_ids"
first="${present[0]}"
grep -v '^#' "$WORK/std_$first.pvar" | cut -f3 | sort -u > "$common"
for coh in "${present[@]:1}"; do
  grep -v '^#' "$WORK/std_$coh.pvar" | cut -f3 | sort -u > "$tmp"
  comm -12 "$common" "$tmp" > "$common.next" && mv "$common.next" "$common"
done
echo "  common variants: $(wc -l < "$common")"

# 3. LD-prune the common set on one cohort, so the merge is small ------------
first="${present[0]}"
say "LD-pruning common variants (on $first) to shrink the merge"
"$PLINK2" --pfile "$WORK/std_$first" --extract "$common" --maf 0.01 --geno 0.05 \
    --indep-pairwise 200 50 0.2 --out "$WORK/pool_prune" >"$WORK/pool_prune.log" 2>&1
[[ -f "$WORK/pool_prune.prune.in" ]] || { echo "  !! pruning produced no prune.in"; tail -15 "$WORK/pool_prune.log"; exit 1; }
echo "  pruned to $(wc -l < "$WORK/pool_prune.prune.in") independent variants"

# 4. convert each cohort's pruned subset to .bed (hardcalls) -----------------
allbeds="$WORK/allbeds.txt"; : > "$allbeds"
for coh in "${present[@]}"; do
  "$PLINK2" --pfile "$WORK/std_$coh" --extract "$WORK/pool_prune.prune.in" \
      --make-bed --out "$WORK/bed_$coh" >"$WORK/bed_$coh.log" 2>&1 \
    && echo "$WORK/bed_$coh" >> "$allbeds" || echo "  !! bed conversion failed for $coh"
done

# 5. sample-merge with plink 1.9, then PCA with plink2 -----------------------
command -v "$PLINK1" >/dev/null || { echo "plink 1.9 not on PATH (module load plink/1.90b)"; exit 1; }
say "merging cohorts with plink 1.9"
"$PLINK1" --merge-list "$allbeds" --make-bed --out "$WORK/merged_pooled" >"$WORK/merge.log" 2>&1
if [[ ! -f "$WORK/merged_pooled.bed" ]]; then
  # a 3+-way merge can fail on strand/allele mismatches; plink writes a .missnp — drop and retry
  if [[ -f "$WORK/merged_pooled-merge.missnp" ]]; then
    echo "  merge flagged $(wc -l < "$WORK/merged_pooled-merge.missnp") problem variants; excluding and retrying"
    : > "$allbeds"
    for coh in "${present[@]}"; do
      "$PLINK1" --bfile "$WORK/bed_$coh" --exclude "$WORK/merged_pooled-merge.missnp" \
          --make-bed --out "$WORK/bed2_$coh" >>"$WORK/merge.log" 2>&1 && echo "$WORK/bed2_$coh" >> "$allbeds"
    done
    "$PLINK1" --merge-list "$allbeds" --make-bed --out "$WORK/merged_pooled" >>"$WORK/merge.log" 2>&1
  fi
fi
[[ -f "$WORK/merged_pooled.bed" ]] || { echo "  !! merge failed — tail:"; tail -20 "$WORK/merge.log"; exit 1; }

say "pooled PCA"
"$PLINK2" --bfile "$WORK/merged_pooled" --pca "$NPC" --out "$WORK/pooled" >"$WORK/pooled.log" 2>&1
[[ -f "$WORK/pooled.eigenvec" ]] && echo "  ok: $WORK/pooled.eigenvec" \
    || { echo "  !! PCA failed — tail:"; tail -15 "$WORK/pooled.log"; exit 1; }
say "done — now re-run assemble_dosage_pcs.R to write SF1_pca_points.tsv"
