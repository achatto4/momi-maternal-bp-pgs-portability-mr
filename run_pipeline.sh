#!/usr/bin/env bash
# ============================================================
# run_pipeline.sh  — run on JHPCE. Submits the WHOLE analytic chain on an
# already-merged genotype set as chained SLURM jobs (each waits for the prior):
#
#   01 phenotypes -> 02 PCA -> 05 PRS(compute->analysis)
#                         \--> 03 GWAS -> 04 instruments -> 04 MR
#
# You run ONE command; SLURM sequences everything via --dependency=afterok.
#
# Usage:
#   bash run_pipeline.sh --merged /path/merged_gsa --analysis-dir /path/analysis_gsa
#        [--bp-window latest] [--n-pcs 10] [--ptb-case-code 2]
#
# Prereqs: merged.{bed,bim,fam} exists (stage 00 done); PGS GRCh38 files present
# (incl. PGS004830_hmPOS_GRCh38 for SAS-SBP); TwoSampleMR installed in R.
# Edit the PATHS block below if your locations differ.
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

# ---------- PATHS (edit if needed) ----------
EPI="/dcs04/nilanjan/data/Anagh/MOMI/Epi Data/MOMI_Selected_Variables_All_Sites.txt"
SUM="/dcs04/nilanjan/data/Anagh/MOMI/sumstat/BP"
PGS_SBP_EUR="$SUM/PRS_BP/SBP/PGS004603_hmPOS_GRCh38.txt.gz"
PGS_DBP_EUR="$SUM/PRS_BP/DBP/PGS004604_hmPOS_GRCh38.txt.gz"
PGS_SBP_SAS="$SUM/PRS_BP/SBP/PGS004830_hmPOS_GRCh38.txt.gz"
PGS_DBP_SAS="$SUM/PRS_BP/DBP/PGS004758_hmPOS_GRCh38.txt.gz"
EXT_SBP="$SUM/mill_EUR/SBP/GCST90310294.h.tsv.gz"
EXT_DBP="$SUM/mill_EUR/DBP/GCST90310295.h.tsv.gz"
PLINK2_MOD="plink/2.00a4.6"; PLINK1_MOD="plink/1.90b"; R_MOD="R"
# --------------------------------------------

MERGED=""; ANALYSIS=""; WIN=latest; NPC=10; PTBCASE=2
while [[ $# -gt 0 ]]; do
  case "$1" in
    --merged) MERGED="$2"; shift 2;;
    --analysis-dir) ANALYSIS="$2"; shift 2;;
    --bp-window) WIN="$2"; shift 2;;
    --n-pcs) NPC="$2"; shift 2;;
    --ptb-case-code) PTBCASE="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$MERGED" && -n "$ANALYSIS" ]] || { echo "usage: run_pipeline.sh --merged PREFIX --analysis-dir DIR" >&2; exit 1; }
[[ -f "$MERGED.bed" ]] || { echo "merged not found: $MERGED.bed" >&2; exit 2; }
mkdir -p "$ANALYSIS" "$ANALYSIS/prs" "$ANALYSIS/mr"

SB(){ sbatch --parsable "$@"; }   # echoes the job id

# 01 — phenotypes / covariates / mothers.keep
J1=$(SB --job-name=ph --cpus-per-task=2 --mem=32G --time=1:00:00 --output="$ANALYSIS/01_pheno_%j.out" \
  --wrap "module load $R_MOD; Rscript '$HERE/01_phenotypes/build_phenotypes.R' \
    --merged-fam '$MERGED.fam' --epi '$EPI' --out-dir '$ANALYSIS' --ptb-case-code $PTBCASE")

# 02 — PCA (+ plots)
J2=$(SB --job-name=pca --dependency=afterok:$J1 --cpus-per-task=4 --mem=32G --time=2:00:00 --output="$ANALYSIS/02_pca_%j.out" \
  --wrap "module load $PLINK2_MOD; \
bash '$HERE/02_pca/run_pca.sh' --merged '$MERGED' --keep '$ANALYSIS/mothers.keep' --out-dir '$ANALYSIS' --n-pcs $NPC --threads 4 --mem-mb 28000 && \
{ module load $R_MOD; Rscript '$HERE/02_pca/plot_pca.R' --eigenvec '$ANALYSIS/pca_analysis/pca_results.eigenvec' --eigenval '$ANALYSIS/pca_analysis/pca_results.eigenval' --covar '$ANALYSIS/covariate_table_analytic_mothers.txt' --out-dir '$ANALYSIS/pca_analysis'; }")

# 05a — PRS compute (4 scores)
J3=$(SB --job-name=prsc --dependency=afterok:$J2 --cpus-per-task=4 --mem=24G --time=2:00:00 --output="$ANALYSIS/prs/05_compute_%j.out" \
  --wrap "module load $PLINK2_MOD; \
bash '$HERE/05_prs/compute_prs.sh' --pgs '$PGS_SBP_EUR' --bfile '$MERGED' --keep '$ANALYSIS/mothers.keep' --out '$ANALYSIS/prs/SBP_EUR'; \
bash '$HERE/05_prs/compute_prs.sh' --pgs '$PGS_DBP_EUR' --bfile '$MERGED' --keep '$ANALYSIS/mothers.keep' --out '$ANALYSIS/prs/DBP_EUR'; \
bash '$HERE/05_prs/compute_prs.sh' --pgs '$PGS_SBP_SAS' --bfile '$MERGED' --keep '$ANALYSIS/mothers.keep' --out '$ANALYSIS/prs/SBP_SAS'; \
bash '$HERE/05_prs/compute_prs.sh' --pgs '$PGS_DBP_SAS' --bfile '$MERGED' --keep '$ANALYSIS/mothers.keep' --out '$ANALYSIS/prs/DBP_SAS'")

# 05b — PRS analysis + PRS-based MR
J4=$(SB --job-name=prsa --dependency=afterok:$J3 --cpus-per-task=2 --mem=16G --time=0:30:00 --output="$ANALYSIS/prs/05_analysis_%j.out" \
  --wrap "module load $R_MOD; Rscript '$HERE/05_prs/prs_analysis.R' \
    --prs-dir '$ANALYSIS/prs' --pheno-dir '$ANALYSIS' --covar '$ANALYSIS/covariate_table_analytic_mothers.txt' \
    --eigenvec '$ANALYSIS/pca_analysis/pca_results.eigenvec' --out-dir '$ANALYSIS/prs' --bp-window $WIN --n-pcs $NPC --ptb-case-code $PTBCASE")

# 03 — internal GWAS (SBP/DBP/PTB)
J5=$(SB --job-name=gwas --dependency=afterok:$J2 --cpus-per-task=8 --mem=64G --time=4:00:00 --output="$ANALYSIS/03_gwas_%j.out" \
  --wrap "module load $PLINK2_MOD; module load $R_MOD; \
bash '$HERE/03_gwas/run_gwas.sh' --merged '$MERGED' --keep '$ANALYSIS/mothers.keep' --pheno-dir '$ANALYSIS' \
  --covar '$ANALYSIS/covariate_table_analytic_mothers.txt' --eigenvec '$ANALYSIS/pca_analysis/pca_results.eigenvec' \
  --out-dir '$ANALYSIS' --bp-window $WIN --n-pcs $NPC --threads 8 --mem-mb 60000")

# 04a — MR instruments from external EUR BP GWAS
J6=$(SB --job-name=instr --dependency=afterok:$J2 --cpus-per-task=4 --mem=24G --time=2:00:00 --output="$ANALYSIS/mr/04_instr_%j.out" \
  --wrap "module load $PLINK1_MOD; \
bash '$HERE/04_mr/prep_instruments.sh' --ext '$EXT_SBP' --trait SBP --merged '$MERGED' --out-dir '$ANALYSIS/mr'; \
bash '$HERE/04_mr/prep_instruments.sh' --ext '$EXT_DBP' --trait DBP --merged '$MERGED' --out-dir '$ANALYSIS/mr'")

# 04b — two-sample MR (needs GWAS PTB outcome + instruments)
J7=$(SB --job-name=mr --dependency=afterok:$J5:$J6 --cpus-per-task=2 --mem=16G --time=1:00:00 --output="$ANALYSIS/mr/04_mr_%j.out" \
  --wrap "module load $R_MOD; \
Rscript '$HERE/04_mr/run_mr.R' --instruments '$ANALYSIS/mr/instruments_SBP.txt' --ptb-gwas '$ANALYSIS/results_PTB_NEW/gwas_PTB_NEW.PHENO.glm.logistic.hybrid' --trait SBP --out-dir '$ANALYSIS/mr'; \
Rscript '$HERE/04_mr/run_mr.R' --instruments '$ANALYSIS/mr/instruments_DBP.txt' --ptb-gwas '$ANALYSIS/results_PTB_NEW/gwas_PTB_NEW.PHENO.glm.logistic.hybrid' --trait DBP --out-dir '$ANALYSIS/mr'")

cat <<EOF

==================  SUBMITTED  ==================
 01 phenotypes   : $J1
 02 PCA          : $J2   (after 01)
 05a PRS compute : $J3   (after 02)
 05b PRS + PRS-MR: $J4   (after 05a)
 03 GWAS         : $J5   (after 02)
 04a instruments : $J6   (after 02)
 04b SNP-MR      : $J7   (after 03 & 04a)
=================================================
 monitor:  squeue -u $USER
 results land in: $ANALYSIS  (prs/, mr/, results_*/)
 a failed job auto-cancels its dependents (afterok).
EOF
