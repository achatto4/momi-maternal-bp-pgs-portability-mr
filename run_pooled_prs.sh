#!/usr/bin/env bash
# ============================================================
# run_pooled_prs.sh — run on JHPCE.
# Combined PRS-based MR across ALL cohorts x platforms (~14k mothers), WITHOUT
# pooling genotypes (which collapses to ~280 SNPs). Each cohort is scored on its
# own variants; per-sample scores are pooled and standardized within group.
#
# Chains (SLURM --dependency):
#   J1 phenotypes (build_phenotypes on the union of all cohort mothers)
#   J2 scoring    (4 PGS x 10 cohort-platform beds)        [parallel to J1]
#   J3 pooled PRS-MR analysis                               [after J1 & J2]
#
# Usage:  bash run_pooled_prs.sh
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUTROOT="/dcs10/chatterj/data/achattop/MOMI/Genomics/pipeline_out"
EPI="/dcs04/nilanjan/data/Anagh/MOMI/Epi Data/MOMI_Selected_Variables_All_Sites.txt"
SUM="/dcs04/nilanjan/data/Anagh/MOMI/sumstat/BP/PRS_BP"
PLINK2_MOD="plink/2.00a4.6"; R_MOD="R"
ANALYSIS="$OUTROOT/analysis_pooled"
COHORTS=(AMANHI-Bangladesh AMANHI-Pakistan AMANHI-Pemba GAPPS-Bangladesh GAPPS-Zambia)
mkdir -p "$ANALYSIS" "$ANALYSIS/prs_scores"

# ---- union of all cohort x platform mothers -> normalized 2-col TAB FID IID, dedup ----
# (gsa beds are plink2 [tab], lpwgs are plink1.9 [space]; awk default FS handles both)
FAM="$ANALYSIS/all_mothers.fam"; RAW="$ANALYSIS/.all_fams.raw"; : > "$RAW"
for c in "${COHORTS[@]}"; do
  for plat in gsa lpwgs; do
    f="$OUTROOT/$c/${plat}_clean.fam"; [[ -f "$f" ]] && cat "$f" >> "$RAW"
  done
done
awk '!seen[$2]++{print $1"\t"$2}' "$RAW" > "$FAM" && rm -f "$RAW"
echo "[fam] union mothers (pre-analytic-filter): $(wc -l < "$FAM")"

SB(){ sbatch --parsable "$@"; }

J1=$(SB --job-name=pph --cpus-per-task=2 --mem=32G --time=1:00:00 --output="$ANALYSIS/01_pheno_%j.out" \
  --wrap "module load $R_MOD; Rscript '$HERE/01_phenotypes/build_phenotypes.R' \
    --merged-fam '$FAM' --epi '$EPI' --out-dir '$ANALYSIS' --ptb-case-code 2")

J2=$(SB --job-name=pscore --cpus-per-task=4 --mem=24G --time=4:00:00 --output="$ANALYSIS/prs_scores/02_score_%j.out" \
  --wrap "module load $PLINK2_MOD; bash '$HERE/05_prs/score_all_cohorts.sh' \
    --outroot '$OUTROOT' --out-dir '$ANALYSIS/prs_scores' --sumdir '$SUM' --plink2 plink2")

J3=$(SB --job-name=pprsmr --dependency=afterok:$J1:$J2 --cpus-per-task=2 --mem=24G --time=1:00:00 \
  --output="$ANALYSIS/prs_scores/03_pooled_mr_%j.out" \
  --wrap "module load $R_MOD; Rscript '$HERE/05_prs/pooled_prs_analysis.R' \
    --prs-dir '$ANALYSIS/prs_scores' --pheno-dir '$ANALYSIS' --out-dir '$ANALYSIS/prs_scores' \
    --bp-window latest --ptb-case-code 2")

cat <<EOF

==================  SUBMITTED  ==================
 01 phenotypes (all cohorts) : $J1
 02 score all cohorts (4x10) : $J2   (parallel to 01)
 03 pooled PRS-MR            : $J3   (after 01 & 02)
=================================================
 monitor: squeue -u $USER
 result:  $ANALYSIS/prs_scores/pooled_prs_mr_results.txt
          tail -20 $ANALYSIS/prs_scores/03_pooled_mr_*.out
EOF
