#!/usr/bin/env bash
# ============================================================
# run_preprocess_lpwgs_dosage.sh  — run on JHPCE.  DOSAGE track (alongside hard-call).
# Builds per-cohort lpWGS DOSAGE filesets (lpwgs_dosage_clean.{pgen,pvar,psam}) that
# keep the imputed DS dosages, so GWAS/PRS can use plink2's dosage support directly.
# It does NOT touch the hard-call lpwgs_clean.* outputs.
#
# Per-cohort: make_sample_map (login, fast) + sbatch 22-chr dosage array + merge/rename.
# Optional: a pooled cross-cohort dosage fileset (--pool). To keep that merge tractable,
#           restrict variants AT CONVERSION with --extract KEEP_concordant.txt (the same
#           concordant KEEP list as the mega build); the pool is then a plain pmerge of
#           already-small per-cohort filesets. (plink2 --pmerge-list does not filter
#           variants mid-merge, so restricting up front is both correct and far cheaper.)
#
# SAFE first run (one cohort, no pool):
#   bash run_preprocess_lpwgs_dosage.sh --cohorts GAPPS-Zambia
# Full run (all 5; per-cohort only):
#   bash run_preprocess_lpwgs_dosage.sh
# Full run + pooled dosage merge on the concordant KEEP set:
#   bash run_preprocess_lpwgs_dosage.sh --extract /path/KEEP_concordant.txt --pool
#
# Options:
#   --cohorts "A B ..."     default: all 5
#   --extract LIST           per-cohort variant restriction at conversion (chr:pos list)
#   --pool                   after all cohorts convert, pmerge them -> one dosage fileset
#   --pool-out PREFIX        default: <OUTROOT>/merged_lpwgs_dosage
#   --force                  re-convert even if lpwgs_dosage_clean.pgen exists
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

# ---------- PATHS (match run_preprocess_lpwgs.sh) ----------
OUTROOT="/dcs10/chatterj/data/achattop/MOMI/Genomics/pipeline_out"
NEW="/dcs10/chatterj/data/achattop/MOMI/Genomics/newdata_Jiong"
EPI="/dcs04/nilanjan/data/Anagh/MOMI/Epi Data/MOMI_Selected_Variables_All_Sites.txt"
PLINK2_MOD="plink/2.00a4.6"; BCFTOOLS_MOD="bcftools"
declare -A PFX=( [AMANHI-Bangladesh]=AMANHIB [AMANHI-Pakistan]=AMANHIP \
                 [AMANHI-Pemba]=AMANHIT [GAPPS-Bangladesh]=GAPPSB [GAPPS-Zambia]=ZAPPS )
ALL_COHORTS=(AMANHI-Bangladesh AMANHI-Pakistan AMANHI-Pemba GAPPS-Bangladesh GAPPS-Zambia)
# -----------------------------------------------------------

COHORTS=("${ALL_COHORTS[@]}"); FORCE=0; EXTRACT=""; DO_POOL=0; POOL_OUT="$OUTROOT/merged_lpwgs_dosage"
NOPASS=1; MAF=0.01   # DEFAULT: keep all variants (Gencove FILTER=PASS is genotype-certainty-biased
                     # against common variants in low-pass data) + dosage MAF>=1%. Use --pass to restore PASS.
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cohorts) read -r -a COHORTS <<< "$2"; shift 2;;
    --extract) EXTRACT="$2"; shift 2;;
    --pool) DO_POOL=1; shift;;
    --pool-out) POOL_OUT="$2"; shift 2;;
    --maf) MAF="$2"; shift 2;;
    --pass) NOPASS=0; shift;;          # restore old (wrong) FILTER=PASS behaviour
    --no-pass) NOPASS=1; shift;;
    --force) FORCE=1; shift;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$EXTRACT" && ! -f "$EXTRACT" ]] && { echo "ERROR: --extract list not found: $EXTRACT" >&2; exit 2; }
NPOPT=""; [[ $NOPASS -eq 1 ]] && NPOPT="--no-pass"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
log "variant policy: $([[ $NOPASS -eq 1 ]] && echo 'NO PASS filter (correct)' || echo 'FILTER=PASS (legacy/wrong)'), dosage MAF>=$MAF"

CONV_IDS=()
for c in "${COHORTS[@]}"; do
  pfx="${PFX[$c]:-}"; [[ -n "$pfx" ]] || { echo "no PID prefix for cohort '$c'"; exit 1; }
  cdir="$NEW/$c"; out="$OUTROOT/$c"; mkdir -p "$out"
  vcf="$cdir/lpwgs_imputed.vcf.gz"; idmap="$cdir/lpwgs_idmap.csv"
  [[ -f "$vcf" && -f "$idmap" ]] || { echo "missing lpWGS files for $c (vcf/idmap)"; exit 2; }
  if [[ -f "$out/lpwgs_dosage_clean.pgen" && $FORCE -eq 0 ]]; then
    log "[skip] $c lpwgs_dosage_clean.pgen exists — use --force to redo"; continue
  fi
  log "[map] $c lpwgs dosage (prefix=$pfx)"
  bash "$HERE/00_preprocess/make_sample_map.sh" --idmap "$idmap" --epi "$EPI" --pid-prefix "$pfx" --out "$out/lpwgs_map" >/dev/null
  log "      $(wc -l < "$out/lpwgs_map.keep") mothers; submitting 22-chr DOSAGE array + merge"
  jarr=$(sbatch --parsable --job-name="lpcd_${c}" --array=1-22 --cpus-per-task=4 --mem=32G --time=8:00:00 \
        --output="$out/lpwgs_dosage_chr_%a.out" \
        --wrap "module load $PLINK2_MOD; module load $BCFTOOLS_MOD 2>/dev/null || true; \
bash '$HERE/00_preprocess/lpwgs_chr_to_plink_dosage.sh' --vcf '$vcf' --keep '$out/lpwgs_map.keep' \
  --out-prefix '$out/lpwgs_dosage_chrtmp' $NPOPT --maf '$MAF' ${EXTRACT:+--extract '$EXTRACT'} --threads 4 --mem-mb 28000")
  jmrg=$(sbatch --parsable --job-name="lpmd_${c}" --dependency=afterok:$jarr --cpus-per-task=4 --mem=48G --time=4:00:00 \
        --output="$out/lpwgs_dosage_merge_%j.out" \
        --wrap "module load $PLINK2_MOD; \
bash '$HERE/00_preprocess/merge_chr_dosage_and_rename.sh' --prefix '$out/lpwgs_dosage_chrtmp' \
  --rename '$out/lpwgs_map.rename' --out '$out/lpwgs_dosage_clean' --plink2 plink2 --threads 4 --mem-mb 44000")
  CONV_IDS+=("$jmrg")
  log "      submitted array $jarr -> merge+rename $jmrg"
done

POOL_JID=""
if [[ $DO_POOL -eq 1 ]]; then
  if [[ ${#COHORTS[@]} -ne 5 ]]; then
    log "[note] --pool given but --cohorts is a subset; skipping pooled merge."
  else
    [[ -n "$EXTRACT" ]] || log "[warn] --pool without --extract: the dosage union may be very large (memory/disk). Consider --extract KEEP_concordant.txt."
    # pmerge the 5 per-cohort dosage filesets (already KEEP-restricted if --extract was used)
    plist="$OUTROOT/merged_lpwgs_dosage.pmerge.list"; : > "$plist"
    for c in "${ALL_COHORTS[@]}"; do echo "$OUTROOT/$c/lpwgs_dosage_clean" >> "$plist"; done
    dep=""; [[ ${#CONV_IDS[@]} -gt 0 ]] && dep="--dependency=afterok:$(IFS=:; echo "${CONV_IDS[*]}")"
    log "[pool] submitting cross-cohort dosage pmerge -> $POOL_OUT ${EXTRACT:+(cohorts restricted to $(basename "$EXTRACT") at conversion)}"
    POOL_JID=$(sbatch --parsable --job-name=merge_lpwgs_dosage $dep --cpus-per-task=4 --mem=64G --time=8:00:00 \
          --output="$OUTROOT/merge_lpwgs_dosage_%j.out" \
          --wrap "module load $PLINK2_MOD; \
plink2 --pmerge-list '$plist' --threads 4 --memory 60000 --out '$POOL_OUT'")
    log "[pool] job $POOL_JID"
  fi
fi

cat <<EOF

==================  SUBMITTED (lpWGS DOSAGE)  ==================
 per-cohort dosage conversions: ${CONV_IDS[*]:-(none / all skipped)}
 pooled dosage merge:           ${POOL_JID:-(not requested)}  $([[ $DO_POOL -eq 1 ]] && echo "-> $POOL_OUT")
===============================================================
 monitor: squeue -u $USER
 per-cohort outputs: <OUTROOT>/<COHORT>/lpwgs_dosage_clean.{pgen,pvar,psam}
 downstream plink2 stages (run_pca/run_gwas/compute_prs) auto-detect pgen vs bed,
 so point --merged / --bfile at a 'lpwgs_dosage_clean' prefix to run on dosages.
EOF
