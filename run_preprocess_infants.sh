#!/usr/bin/env bash
# ============================================================
# run_preprocess_infants.sh — build INFANT filesets for S15/B32.
#
# WHY THIS IS A THIN WRAPPER AND NOT A NEW PIPELINE. The fetal-vs-maternal decomposition
# regresses birthweight on the maternal AND fetal polygenic score at once. If those two
# scores come from filesets built under different QC, the comparison between them is
# meaningless -- any difference could be a QC artefact rather than biology. So this script
# calls EXACTLY the same conversion scripts as the maternal path
# (lpwgs_chr_to_plink.sh / merge_chr_and_rename.sh / vcf_to_plink.sh), with identical
# parameters, changing only the sample list. Nothing about the variant QC differs.
#
# The only new input is make_sample_map.sh --mode infant, which selects children instead of
# mothers and renames them to their Epi BABY_ID. It also emits <out>.keep.pairs mapping each
# BABY_ID to its mother PARTICIPANT_ID, which B32 reads to form the trios.
#
# OUTPUT NAMING. <OUTROOT>/<COHORT>/infant_gsa_clean and infant_lpwgs_clean, deliberately
# alongside the maternal gsa_clean / lpwgs_clean rather than in a separate tree, so the two
# are obviously parallel and a later reader cannot mistake one for the other.
#
#   bash run_preprocess_infants.sh --outroot DIR [--cohorts "A B"] [--platforms "gsa lpwgs"] [--force]
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
: "${EPI:?export EPI first}"

NEW=/dcs10/chatterj/data/achattop/MOMI/Genomics/newdata_Jiong
# Module names must match run_preprocess_lpwgs.sh EXACTLY. `plink/1.90b6.21` was written here
# on 2026-07-20 by guessing a version suffix; that module does not exist on JHPCE and every
# merge job died in 1 second with a Lmod error. The maternal runner uses `plink/1.90b`.
PLINK2_MOD=plink/2.00a4.6; PLINK1_MOD=plink/1.90b; BCFTOOLS_MOD=bcftools

# ---- VARIANT FILTERS: must be IDENTICAL to the maternal runners ---------------------------
# These are not free parameters. The whole point of the infant filesets is that maternal and
# fetal polygenic scores be computed on the same variants; if the two differ, the scores are
# not comparable and the fetal-vs-maternal decomposition is meaningless. Copied verbatim from
# run_preprocess_lpwgs_dosage.sh (MAF=0.01, NOPASS=1) and run_preprocess_lpwgs.sh (--maf 0).
#
# Getting this wrong on 2026-07-20 cost a full conversion cycle: `--maf 0` with no
# `--no-pass` produced 33.5M infant variants vs 8.1M maternal, through a different bcftools
# filter, and the resulting infant scores did not track maternal genotype at all
# (r = 0.02-0.13 vs the ~0.5 that descent requires). GSA was unaffected because its filters
# happened to match. If either maternal runner changes its filters, change these too.
DOSAGE_MAF=0.01           # run_preprocess_lpwgs_dosage.sh: MAF=0.01
DOSAGE_NPOPT="--no-pass"  # run_preprocess_lpwgs_dosage.sh: NOPASS=1 -> keeps LOWCONF variants
HARDCALL_MAF=0            # run_preprocess_lpwgs.sh passes --maf 0
declare -A PFX=( [AMANHI-Bangladesh]=AMANHIB [AMANHI-Pakistan]=AMANHIP \
                 [AMANHI-Pemba]=AMANHIT [GAPPS-Bangladesh]=GAPPSB [GAPPS-Zambia]=ZAPPS )
ALL_COHORTS=(AMANHI-Bangladesh AMANHI-Pakistan AMANHI-Pemba GAPPS-Bangladesh GAPPS-Zambia)

OUTROOT="${MOMI_OUTROOT:-}"; FORCE=0; MERGE_ONLY=0
COHORTS=("${ALL_COHORTS[@]}"); PLATFORMS=(gsa lpwgs)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outroot) OUTROOT="$2"; shift 2;;
    --cohorts) read -r -a COHORTS <<< "$2"; shift 2;;
    --platforms) read -r -a PLATFORMS <<< "$2"; shift 2;;
    --merge-only) MERGE_ONLY=1; shift;;
    --force) FORCE=1; shift;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

# ---- --merge-only : resubmit just the lpWGS chromosome merges ----------------------------
# WHY THIS MODE EXISTS. The lpWGS path submits a 22-task chromosome array plus a merge job
# gated on --dependency=afterok. If the dependency is not satisfied -- which happened on
# 2026-07-20 when an interactive session dropped mid-array -- the merge silently sits in
# DependencyNeverSatisfied and then fails, while all 22 chromosome filesets sit on disk
# perfectly intact. Redoing the array costs ~3 hours per cohort (bcftools dominates); redoing
# the merge costs minutes. Without this mode the recovery is a hand-pasted sbatch line, which
# went wrong three times in one afternoon: twice from an unset $MOMI_OUTROOT expanding to
# nonsense paths, and once from a module name that does not exist on this cluster.
#
# It refuses to submit unless all 22 chromosome filesets and the .rename are actually present,
# so a missing input fails HERE with a clear message rather than 1 second into a batch job.
if [[ $MERGE_ONLY -eq 1 ]]; then
  [[ -n "$OUTROOT" ]] || { echo "need --outroot (or export MOMI_OUTROOT)" >&2; exit 1; }
  sub=0
  for c in "${COHORTS[@]}"; do
    out="$OUTROOT/$c"; tag=infant_lpwgs
    n=$(ls "$out/${tag}_chrtmp_chr"*.bed 2>/dev/null | wc -l)
    if [[ "$n" -ne 22 ]]; then
      echo "[skip] $c: $n/22 chromosome filesets — rerun the array, not the merge"; continue
    fi
    if [[ ! -f "$out/${tag}_map.rename" ]]; then
      echo "[skip] $c: no ${tag}_map.rename"; continue
    fi
    if [[ -f "$out/${tag}_clean.bed" && $FORCE -eq 0 ]]; then
      echo "[skip] $c: ${tag}_clean.bed exists ($(wc -l < "$out/${tag}_clean.fam") infants) — --force to redo"; continue
    fi
    jid=$(sbatch --parsable --job-name="inf_lpm_$c" --cpus-per-task=4 --mem=48G --time=4:00:00 \
          --output="$out/${tag}_merge_%j.out" \
          --wrap "module load $PLINK1_MOD; \
bash '$HERE/00_preprocess/merge_chr_and_rename.sh' --prefix '$out/${tag}_chrtmp' \
  --rename '$out/${tag}_map.rename' --out '$out/${tag}_clean' \
  --plink1 plink --threads 4 --mem-mb 44000")
    echo "[merge] $c -> job $jid  (22/22 filesets present)"
    sub=$((sub+1))
  done
  echo "submitted $sub merge job(s). watch: squeue -u \$USER"
  exit 0
fi
[[ -n "$OUTROOT" ]] || { echo "need --outroot (or export MOMI_OUTROOT)" >&2; exit 1; }
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

for c in "${COHORTS[@]}"; do
  pfx="${PFX[$c]:-}"; [[ -n "$pfx" ]] || { echo "no PID prefix for '$c'" >&2; exit 1; }
  cdir="$NEW/$c"; out="$OUTROOT/$c"; mkdir -p "$out"

  for plat in "${PLATFORMS[@]}"; do
    ## both lpWGS tracks (hard call and dosage) read the SAME idmap and the SAME VCF --
    ## they differ only in whether plink2 keeps DS dosages. Resolve that here, before the
    ## existence check, or `lpwgs_dosage_idmap.csv` is looked for and the platform silently
    ## skips.
    src_plat="$plat"; [[ "$plat" == "lpwgs_dosage" ]] && src_plat="lpwgs"
    idmap="$cdir/${src_plat}_idmap.csv"
    [[ -f "$idmap" ]] || { log "[skip] $c $plat: no idmap ($idmap)"; continue; }
    case "$plat" in
      gsa)                vcf="$cdir/gsa_reimputed.vcf.gz";;
      lpwgs|lpwgs_dosage) vcf="$cdir/lpwgs_imputed.vcf.gz";;
      *)                  log "[skip] unknown platform $plat"; continue;;
    esac
    [[ -f "$vcf" ]] || { log "[skip] $c $plat: no vcf"; continue; }

    tag="infant_${plat}"
    ## dosage output is a .pgen, hard-call output a .bed
    done_marker="$out/${tag}_clean.bed"
    [[ "$plat" == "lpwgs_dosage" ]] && done_marker="$out/${tag}_clean.pgen"
    if [[ -f "$done_marker" && $FORCE -eq 0 ]]; then
      log "[skip] $c $(basename "$done_marker") exists — --force to redo"
      continue
    fi
    ## sample map first (fast, awk, runs here on the login node)
    log "[map] $c $plat infants (prefix=$pfx)"
    bash "$HERE/00_preprocess/make_sample_map.sh" --idmap "$idmap" --epi "$EPI" \
         --pid-prefix "$pfx" --out "$out/${tag}_map" --mode infant >/dev/null
    nkeep=$(wc -l < "$out/${tag}_map.keep" 2>/dev/null || echo 0)
    if [[ "$nkeep" -eq 0 ]]; then
      log "      0 infants mapped — skipping (see $out/${tag}_map.report)"
      continue
    fi
    log "      $nkeep infants mapped"

    if [[ "$plat" == "lpwgs_dosage" ]]; then
      ## ---- DOSAGE track, added 2026-07-20 -------------------------------------------
      ## The hard-call infant lpWGS filesets produced polygenic scores that did NOT track
      ## maternal genotype: r(mother, infant) came out at 0.02-0.13 where descent requires
      ## ~0.5, while the GSA filesets gave 0.41-0.58 in the same mothers. lpWGS is LOW-PASS
      ## sequencing, so hard-calling discards most of the information; the maternal pipeline
      ## never hit this because it scores lpWGS from DOSAGES. B32 now gates on that
      ## correlation and drops any platform that fails, which is why the first real run used
      ## GSA alone (1,300 pairs across 2 cohorts instead of ~9,500 across 5).
      ## This track rebuilds the infant lpWGS filesets as dosages, mirroring
      ## run_preprocess_lpwgs_dosage.sh exactly so maternal and fetal scores stay comparable.
      ## FILTERS MUST MATCH run_preprocess_lpwgs_dosage.sh EXACTLY -- see DOSAGE_MAF/NOPASS
      ## above. The first attempt (2026-07-20) passed `--maf 0` and omitted `--no-pass`,
      ## producing 33.5M infant variants against 8.1M maternal, built through a DIFFERENT
      ## bcftools filter as well. Maternal and fetal scores computed on different variant
      ## sets are not comparable, and the symptom was r(mother, infant) = 0.02-0.13 where
      ## descent requires ~0.5 -- while GSA, whose filters DID match, gave 0.41-0.58.
      jarr=$(sbatch --parsable --job-name="infd_lpc_${c}" --array=1-22 --cpus-per-task=4 \
            --mem=32G --time=8:00:00 --output="$out/${tag}_chr_%a.out" \
            --wrap "module load $PLINK2_MOD; module load $BCFTOOLS_MOD 2>/dev/null || true; \
bash '$HERE/00_preprocess/lpwgs_chr_to_plink_dosage.sh' --vcf '$vcf' --keep '$out/${tag}_map.keep' \
  --out-prefix '$out/${tag}_chrtmp' $DOSAGE_NPOPT --maf '$DOSAGE_MAF' --threads 4 --mem-mb 28000")
      jmrg=$(sbatch --parsable --job-name="infd_lpm_${c}" --dependency=afterok:$jarr \
            --cpus-per-task=4 --mem=48G --time=4:00:00 --output="$out/${tag}_merge_%j.out" \
            --wrap "module load $PLINK2_MOD; \
bash '$HERE/00_preprocess/merge_chr_dosage_and_rename.sh' --prefix '$out/${tag}_chrtmp' \
  --rename '$out/${tag}_map.rename' --out '$out/${tag}_clean' --plink2 plink2 --threads 4 --mem-mb 44000")
      log "      submitted DOSAGE array $jarr -> merge+rename $jmrg"
    elif [[ "$plat" == "lpwgs" ]]; then
      ## same 22-way per-chromosome array + merge as the maternal lpWGS path
      jarr=$(sbatch --parsable --job-name="inf_lpc_${c}" --array=1-22 --cpus-per-task=4 \
            --mem=32G --time=8:00:00 --output="$out/${tag}_chr_%a.out" \
            --wrap "module load $PLINK2_MOD; module load $BCFTOOLS_MOD 2>/dev/null || true; \
bash '$HERE/00_preprocess/lpwgs_chr_to_plink.sh' --vcf '$vcf' --keep '$out/${tag}_map.keep' \
  --out-prefix '$out/${tag}_chrtmp' --maf '$HARDCALL_MAF' --threads 4 --mem-mb 28000")
      jmrg=$(sbatch --parsable --job-name="inf_lpm_${c}" --dependency=afterok:$jarr \
            --cpus-per-task=4 --mem=48G --time=4:00:00 --output="$out/${tag}_merge_%j.out" \
            --wrap "module load $PLINK1_MOD; \
bash '$HERE/00_preprocess/merge_chr_and_rename.sh' --prefix '$out/${tag}_chrtmp' \
  --rename '$out/${tag}_map.rename' --out '$out/${tag}_clean' --plink1 plink --threads 4 --mem-mb 44000")
      log "      submitted array $jarr -> merge+rename $jmrg"
    else
      ## GSA is a single VCF, one job, same script the maternal path used
      jid=$(sbatch --parsable --job-name="inf_gsa_${c}" --cpus-per-task=8 \
           --mem=60G --time=8:00:00 --output="$out/${tag}_%j.out" \
           --wrap "module load $PLINK2_MOD; module load $BCFTOOLS_MOD 2>/dev/null || true; \
bash '$HERE/00_preprocess/vcf_to_plink.sh' --vcf '$vcf' --platform gsa \
  --keep '$out/${tag}_map.keep' --rename '$out/${tag}_map.rename' --out '$out/${tag}_clean' \
  --threads 8 --mem-mb 56000")
      log "      submitted $jid"
    fi
  done
done

cat <<'EOF'

submitted. watch with:  squeue -u $USER
when it clears:
  1) score the infants  — the sscore filenames must not collide with the maternal ones,
     so pass a distinct platform label:
       bash $MOMI_PIPE/05_prs/score_panel.sh --panel $MOMI_PIPE/05_prs/transfer_panel.tsv \
         --scores-dir $MOMI_SCORES --outroot $MOMI_OUTROOT --out-dir $SSC \
         --platforms "infant_gsa infant_lpwgs"
  2) then B32:  bash $MOMI_PIPE/run_build.sh --only B32
EOF
