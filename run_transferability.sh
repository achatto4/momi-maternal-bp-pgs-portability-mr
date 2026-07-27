#!/usr/bin/env bash
# ============================================================
# run_transferability.sh   — PRS TRANSFERABILITY STUDY (one command)
# "How well does each ancestry's BP PGS predict measured BP in each of our
#  cohorts/ancestries?"  Downloads the panel scores, scores every PGS on every
#  cohort x platform, then builds the PGS-ancestry x cohort-ancestry R2 matrix.
#
# Steps:
#   1. download each panel PGS (GRCh38 harmonized) from the PGS Catalog FTP   [login node, has internet]
#   2. sbatch score_panel.sh        (PGS x cohort x platform -> .sscore)
#   3. sbatch --dependency=afterok  transferability_analysis.R (-> matrices)
#
# Usage (from repo root on JHPCE):
#   bash bp_ptb_pipeline/run_transferability.sh \
#        --outroot   /dcs10/chatterj/data/achattop/MOMI/Genomics/pipeline_out \
#        --pheno-dir /dcs10/chatterj/data/achattop/MOMI/Genomics/pipeline_out/phenotypes \
#        --workdir   /dcs10/chatterj/data/achattop/MOMI/Genomics/transferability
# Notes:
#   --outroot must contain <COHORT>/{gsa,lpwgs}_clean.bed  (the per-cohort beds the
#     GSA + lpWGS converters already produce — same ones fed to the all-10 merge)
#   --pheno-dir must contain mothers_{SBP,DBP}_<win>_final.pheno + covariate_table_analytic_mothers.txt
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; PRS="$HERE/05_prs"
OUTROOT=""; PHENO=""; WORK=""; WIN="latest"; ACCOUNT=""; PART=""; NOSUB=0
while [[ $# -gt 0 ]]; do case "$1" in
  --outroot) OUTROOT="$2"; shift 2;;
  --pheno-dir) PHENO="$2"; shift 2;;
  --workdir) WORK="$2"; shift 2;;
  --bp-window) WIN="$2"; shift 2;;
  --account) ACCOUNT="$2"; shift 2;;
  --partition) PART="$2"; shift 2;;
  --no-submit) NOSUB=1; shift;;          # just download + write scripts, don't sbatch
  *) echo "Unknown arg: $1" >&2; exit 1;;
esac; done
for v in OUTROOT PHENO WORK; do [[ -n "${!v}" ]] || { echo "missing --$v"; exit 1; }; done
PANEL="$PRS/transfer_panel.tsv"; SCORES="$WORK/scores"; SSC="$WORK/sscore"; RES="$WORK/results"; LOG="$WORK/logs"
mkdir -p "$SCORES" "$SSC" "$RES" "$LOG"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
FTP="https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/%s/ScoringFiles/Harmonized/%s_hmPOS_GRCh38.txt.gz"

# ---- 1. download panel scores (login node) ----
log "downloading panel scores -> $SCORES"
tail -n +2 "$PANEL" | while IFS=$'\t' read -r pid trait anc src; do
  [[ -n "$pid" ]] || continue
  out="$SCORES/${pid}.txt.gz"
  if [[ -s "$out" ]]; then log "have $pid"; continue; fi
  url="$(printf "$FTP" "$pid" "$pid")"
  log "wget $pid ($anc/$trait) $url"
  wget -q -O "$out" "$url" || { log "FAILED download $pid"; rm -f "$out"; }
done
got=$(ls -1 "$SCORES"/*.txt.gz 2>/dev/null | wc -l)
log "have $got/$(($(wc -l < "$PANEL")-1)) panel score files"
[[ "$got" -gt 0 ]] || { echo "no score files downloaded; aborting"; exit 1; }

SB=(sbatch --parsable); [[ -n "$ACCOUNT" ]] && SB+=(--account "$ACCOUNT"); [[ -n "$PART" ]] && SB+=(--partition "$PART")

if [[ "$NOSUB" -eq 1 ]]; then
  log "DRY: scores in $SCORES. To run scoring locally:"
  echo "  bash $PRS/score_panel.sh --panel $PANEL --scores-dir $SCORES --outroot $OUTROOT --out-dir $SSC"
  echo "  Rscript $PRS/transferability_analysis.R --panel $PANEL --sscore-dir $SSC --pheno-dir $PHENO --out-dir $RES --bp-window $WIN"
  exit 0
fi

# module names match the rest of the pipeline (run_pipeline.sh / run_pooled_prs.sh)
PLINK2_MOD="plink/2.00a4.6"; R_MOD="R"
# ---- 2. score panel (one SLURM job; compute_prs handles plink2 module) ----
J1=$("${SB[@]}" -J prs_score -t 8:00:00 --mem=20G -c 4 -o "$LOG/score_%j.out" -e "$LOG/score_%j.err" \
   --wrap "module load $PLINK2_MOD; \
           bash '$PRS/score_panel.sh' --panel '$PANEL' --scores-dir '$SCORES' --outroot '$OUTROOT' --out-dir '$SSC'")
log "submitted scoring job $J1"

# ---- 3. transferability matrix (after scoring) ----
J2=$("${SB[@]}" --dependency=afterok:"$J1" -J prs_xfer -t 1:00:00 --mem=8G -c 1 \
   -o "$LOG/xfer_%j.out" -e "$LOG/xfer_%j.err" \
   --wrap "module load $R_MOD; \
           Rscript '$PRS/transferability_analysis.R' --panel '$PANEL' --sscore-dir '$SSC' \
              --pheno-dir '$PHENO' --out-dir '$RES' --bp-window '$WIN'")
log "submitted analysis job $J2 (afterok:$J1)"
echo
log "watch:   squeue -u \$USER"
log "results: $RES/transferability_matrix_*.txt  and  transferability_long.txt"
