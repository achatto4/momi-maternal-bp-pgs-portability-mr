#!/usr/bin/env bash
# ============================================================
# submit_score_array.sh — fill missing PGS scorings as a SLURM JOB ARRAY.
#
# WHY AN ARRAY. score_panel.sh loops score x cohort x platform serially. Filling the 39
# combos left unscored for the three multi-ancestry scores (PGS003968, PGS005008,
# PGS003964) that way means ~39 sequential plink2 runs, several against 14M-variant lpWGS
# dosage filesets. Splitting on (score, cohort) gives 15 independent tasks.
#
# SAFE TO RUN CONCURRENTLY. Each task writes only <pid>__<cohort>__<plat>.sscore, so the
# tasks touch disjoint files. score_panel.sh already skips any .sscore that exists, which
# makes every task idempotent -- a failed or requeued task can simply be rerun.
#
# THROTTLE. The array is submitted with %5, capping concurrent tasks at five. JHPCE
# documents this "step size" specifically to avoid saturating a shared file server, and
# these jobs are I/O bound on the same filesystem, so running all 15 at once would be
# antisocial and probably slower.
#   https://jhpce.jhu.edu/slurm/crafting-jobs/#job-arrays
#
#   bash bin/submit_score_array.sh [--panel FILE] [--throttle 5] [--dry-run]
# ============================================================
set -euo pipefail
PIPE="${MOMI_PIPE:?export MOMI_PIPE first}"
PANEL="$PIPE/05_prs/transfer_panel.tsv"
THROTTLE=5
DRY=0
## Platforms default to the MATERNAL filesets. Pass --platforms "infant_gsa infant_lpwgs"
## to score the infant filesets for S15/B32. Added 2026-07-20: the infant scoring was being
## run serially through score_panel.sh, which is ~60 plink2 invocations against filesets up
## to 16 GB -- an hour of wall clock for work that is embarrassingly parallel.
PLATFORMS=(gsa lpwgs lpwgs_dosage)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --panel) PANEL="$2"; shift 2;;
    --throttle) THROTTLE="$2"; shift 2;;
    --platforms) read -r -a PLATFORMS <<< "$2"; shift 2;;
    --dry-run) DRY=1; shift;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
: "${MOMI_SCORES:?export MOMI_SCORES}"; : "${MOMI_OUTROOT:?export MOMI_OUTROOT}"; : "${SSC:?export SSC}"

COHORTS=(AMANHI-Bangladesh AMANHI-Pakistan AMANHI-Pemba GAPPS-Bangladesh GAPPS-Zambia)
WORK="$PIPE/results/current/qc/score_array"; mkdir -p "$WORK"
TASKS="$WORK/tasks.tsv"

## ---- build the task list: one line per (score, cohort) with work left to do ----
## Scores already complete are skipped entirely rather than submitted as no-op tasks, so the
## array size reflects real work and a glance at squeue tells you what is outstanding.
: > "$TASKS"
tail -n +2 "$PANEL" | while IFS=$'\t' read -r pid trait anc src; do
  [[ -n "${pid:-}" ]] || continue
  [[ -f "$MOMI_SCORES/${pid}.txt.gz" ]] || { echo "skip $pid: no weights" >&2; continue; }
  for c in "${COHORTS[@]}"; do
    need=0
    for plat in "${PLATFORMS[@]}"; do
      bed="$MOMI_OUTROOT/$c/${plat}_clean"
      [[ -f "$bed.bed" || -f "$bed.pgen" ]] || continue
      [[ -f "$SSC/${pid}__${c}__${plat}.sscore" ]] || need=1
    done
    ## if/fi rather than `[[ ... ]] && printf` -- clearer, and immune to set -e edge cases
    ## in the && form regardless of whether they bite in this particular bash version.
    if [[ $need -eq 1 ]]; then
      printf '%s\t%s\n' "$pid" "$c" >> "$TASKS"
    fi
  done
done

N=$(wc -l < "$TASKS")
if [[ "$N" -eq 0 ]]; then echo "nothing to do — every panel score is fully scored."; exit 0; fi
echo "tasks to run ($N):"; cat "$TASKS"

RUNNER="$WORK/run_task.sh"
cat > "$RUNNER" <<'EOF'
#!/bin/bash
set -uo pipefail
set +u; module load plink/2.00a4.6; set -u
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKS")
PID=$(cut -f1 <<< "$LINE"); COH=$(cut -f2 <<< "$LINE")
echo "[task ${SLURM_ARRAY_TASK_ID}] $PID on $COH"
# a one-row panel so score_panel.sh does exactly this score
TMP=$(mktemp); printf 'PGS_id\ttrait\tancestry\tsource\n%s\tNA\tNA\tNA\n' "$PID" > "$TMP"
bash "$MOMI_PIPE/05_prs/score_panel.sh" \
  --panel "$TMP" --scores-dir "$MOMI_SCORES" --outroot "$MOMI_OUTROOT" \
  --out-dir "$SSC" --cohorts "$COH" --platforms "$PLATS"
rc=$?
rm -f "$TMP"
exit $rc
EOF
chmod +x "$RUNNER"

if [[ $DRY -eq 1 ]]; then echo "(dry run — not submitting)"; exit 0; fi

sbatch --job-name=momi_score \
       --array="1-${N}%${THROTTLE}" \
       --mem=24G --cpus-per-task=4 --time=4:00:00 \
       --output="$WORK/score_%A_%a.out" \
       --export=ALL,TASKS="$TASKS",PLATS="${PLATFORMS[*]}" \
       "$RUNNER"

echo
echo "submitted. watch with:  squeue -u \$USER"
echo "when it clears:         bash \$MOMI_PIPE/run_build.sh --from B02 --to B17"
