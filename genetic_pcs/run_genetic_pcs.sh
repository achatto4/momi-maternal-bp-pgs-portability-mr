#!/usr/bin/env bash
# run_genetic_pcs.sh -- the within-cohort joint-platform principal-component analysis, then the descriptive
# cross-cohort PCA with the pairwise FST, then the joint principal-component file and its quality gates.
#   bash genetic_pcs/run_genetic_pcs.sh all                 everything, in order
#   bash genetic_pcs/run_genetic_pcs.sh cohort <cohort>     one cohort's joint PCA
#   bash genetic_pcs/run_genetic_pcs.sh cross               the cross-cohort PCA and FST (after all five cohorts)
#   bash genetic_pcs/run_genetic_pcs.sh assemble            the joint-PC file and the quality gates
# The joint-PC file is written to $JPCA_OUT/participant_level/joint_pcs.rds (set pc_file in config.R to it) and
# the gates to $JPCA_OUT/aggregate/joint_pca_decision.tsv; the analysis must not use the file unless every gate passes.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/config_pcs.sh"
mkdir -p "$JPCA_OUT/aggregate"
abspath(){ printf '%s/%s' "$(cd "$(dirname "$1")" && pwd)" "$(basename "$1")"; }   # the stage scripts change directory
export JPCA_OUT="$(cd "$JPCA_OUT" && pwd)" JPCA_KEEP="$(abspath "$JPCA_KEEP")" JPCA_DROPS="$(abspath "$JPCA_DROPS")"
case "${1:-all}" in
  cohort)   bash "$HERE/cohort_genotypes.sh" "$2" "$JPCA_ENCODING" ;;
  cross)    bash "$HERE/cross_cohort_genotypes.sh" "$JPCA_ENCODING" ;;
  assemble) Rscript "$HERE/assemble_joint_pcs.R" --rundir "$JPCA_OUT" --epi "$JPCA_EPI" --sscore-dir "$JPCA_SSC" --keep "$JPCA_KEEP" --encoding "$JPCA_ENCODING" \
              --cohorts "${JPCA_COHORTS[*]}" ${JPCA_TEST_EXPECTED_N:+--expected-n "$JPCA_TEST_EXPECTED_N"} \
              ${JPCA_TEST_MIN_VARIANTS:+--min-variants "$JPCA_TEST_MIN_VARIANTS"} ;;
  all)      for c in "${JPCA_COHORTS[@]}"; do bash "$0" cohort "$c"; done
            bash "$0" cross
            bash "$0" assemble ;;
  *)        sed -n '2,9p' "$0"; exit 2 ;;
esac
