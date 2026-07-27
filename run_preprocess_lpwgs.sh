#!/usr/bin/env bash
# ============================================================
# run_preprocess_lpwgs.sh  — run on JHPCE.
# Submits lpWGS conversion for the 5 cohorts, then (when all succeed) the FINAL
# all-10 merge (5 GSA + 5 lpWGS) -> merged_all, via SLURM --dependency=afterok.
#
# Per-cohort: make_sample_map (login, fast) + sbatch vcf_to_plink --platform lpwgs.
# lpWGS jobs are HEAVY (100s of GB; bcftools -f PASS | plink2). Give them lots of time.
#
# SAFE first run (test one cohort, no merge):
#   bash run_preprocess_lpwgs.sh --cohorts GAPPS-Zambia --no-final-merge
# Full run (all 5 + final merge):
#   bash run_preprocess_lpwgs.sh
#
# Options:
#   --cohorts "A B ..."     default: all 5
#   --no-final-merge        just convert; skip the all-10 merge
#   --merged-out PREFIX      default: <OUTROOT>/merged_all
#   --force                  re-convert even if lpwgs_clean.bed exists
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

# ---------- PATHS (edit if needed) ----------
OUTROOT="/dcs10/chatterj/data/achattop/MOMI/Genomics/pipeline_out"
NEW="/dcs10/chatterj/data/achattop/MOMI/Genomics/newdata_Jiong"
EPI="/dcs04/nilanjan/data/Anagh/MOMI/Epi Data/MOMI_Selected_Variables_All_Sites.txt"
PLINK2_MOD="plink/2.00a4.6"; PLINK1_MOD="plink/1.90b"; BCFTOOLS_MOD="bcftools"
R2MIN=0.8
declare -A PFX=( [AMANHI-Bangladesh]=AMANHIB [AMANHI-Pakistan]=AMANHIP \
                 [AMANHI-Pemba]=AMANHIT [GAPPS-Bangladesh]=GAPPSB [GAPPS-Zambia]=ZAPPS )
ALL_COHORTS=(AMANHI-Bangladesh AMANHI-Pakistan AMANHI-Pemba GAPPS-Bangladesh GAPPS-Zambia)
# --------------------------------------------

COHORTS=("${ALL_COHORTS[@]}"); DO_MERGE=1; MERGED_OUT="$OUTROOT/merged_all"; FORCE=0; EXTRACT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cohorts) read -r -a COHORTS <<< "$2"; shift 2;;
    --no-final-merge) DO_MERGE=0; shift;;
    --merged-out) MERGED_OUT="$2"; shift 2;;
    --extract) EXTRACT="$2"; shift 2;;   # union-KEEP variant list (chr:pos) for mega build
    --force) FORCE=1; shift;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

CONV_IDS=()
for c in "${COHORTS[@]}"; do
  pfx="${PFX[$c]:-}"; [[ -n "$pfx" ]] || { echo "no PID prefix for cohort '$c'"; exit 1; }
  cdir="$NEW/$c"; out="$OUTROOT/$c"; mkdir -p "$out"
  vcf="$cdir/lpwgs_imputed.vcf.gz"; idmap="$cdir/lpwgs_idmap.csv"
  [[ -f "$vcf" && -f "$idmap" ]] || { echo "missing lpWGS files for $c (vcf/idmap)"; exit 2; }
  if [[ -f "$out/lpwgs_clean.bed" && $FORCE -eq 0 ]]; then
    log "[skip] $c lpwgs_clean.bed exists ($(wc -l < "$out/lpwgs_clean.fam") mothers) — use --force to redo"
    continue
  fi
  log "[map] $c lpwgs (prefix=$pfx)"
  bash "$HERE/00_preprocess/make_sample_map.sh" --idmap "$idmap" --epi "$EPI" --pid-prefix "$pfx" --out "$out/lpwgs_map" >/dev/null
  nkeep=$(wc -l < "$out/lpwgs_map.keep")
  log "      $nkeep mothers; submitting 22-chromosome array + merge"
  # 22-way parallel per-chromosome conversion (uses the VCF tabix index)
  jarr=$(sbatch --parsable --job-name="lpc_${c}" --array=1-22 --cpus-per-task=4 --mem=32G --time=8:00:00 \
        --output="$out/lpwgs_chr_%a.out" \
        --wrap "module load $PLINK2_MOD; module load $BCFTOOLS_MOD 2>/dev/null || true; \
bash '$HERE/00_preprocess/lpwgs_chr_to_plink.sh' --vcf '$vcf' --keep '$out/lpwgs_map.keep' \
  --out-prefix '$out/lpwgs_chrtmp' --maf 0 ${EXTRACT:+--extract '$EXTRACT'} --threads 4 --mem-mb 28000")
  # concatenate the 22 chromosomes + rename to PARTICIPANT_ID (after all array tasks succeed)
  jmrg=$(sbatch --parsable --job-name="lpm_${c}" --dependency=afterok:$jarr --cpus-per-task=4 --mem=48G --time=4:00:00 \
        --output="$out/lpwgs_merge_%j.out" \
        --wrap "module load $PLINK1_MOD; \
bash '$HERE/00_preprocess/merge_chr_and_rename.sh' --prefix '$out/lpwgs_chrtmp' \
  --rename '$out/lpwgs_map.rename' --out '$out/lpwgs_clean' --plink1 plink --threads 4 --mem-mb 44000")
  CONV_IDS+=("$jmrg")
  log "      submitted array $jarr -> merge+rename $jmrg"
done

if [[ $DO_MERGE -eq 1 ]]; then
  # only do the all-10 merge if we're running the full 5 cohorts
  if [[ ${#COHORTS[@]} -ne 5 ]]; then
    log "[note] --cohorts is a subset; skipping final all-10 merge. Re-run without --cohorts (and with all lpwgs_clean present) to merge."
  else
    beds=()
    for c in "${ALL_COHORTS[@]}"; do beds+=("$OUTROOT/$c/gsa_clean" "$OUTROOT/$c/lpwgs_clean"); done
    dep=""
    if [[ ${#CONV_IDS[@]} -gt 0 ]]; then dep="--dependency=afterok:$(IFS=:; echo "${CONV_IDS[*]}")"; fi
    log "[merge] submitting all-10 merge -> $MERGED_OUT ${dep:+(after: ${CONV_IDS[*]})}"
    mjid=$(sbatch --parsable --job-name=merge_all $dep --cpus-per-task=4 --mem=48G --time=6:00:00 \
          --output="$OUTROOT/merge_all_%j.out" \
          --wrap "module load $PLINK1_MOD; \
bash '$HERE/00_preprocess/merge_filesets.sh' --out '$MERGED_OUT' --plink1 plink --geno 0.05 --mind 0.05 --maf 0.01 \
  $(printf "'%s' " "${beds[@]}")")
    log "[merge] job $mjid"
  fi
fi

cat <<EOF

==================  SUBMITTED  ==================
 lpWGS conversions: ${CONV_IDS[*]:-（none / all skipped)}
 final all-10 merge: ${mjid:-(not submitted)}  -> $MERGED_OUT
=================================================
 monitor: squeue -u $USER
 after merge_all completes, run the analytic chain on it:
   bash $HERE/run_pipeline.sh --merged '$MERGED_OUT' --analysis-dir '$OUTROOT/analysis_all'
EOF
