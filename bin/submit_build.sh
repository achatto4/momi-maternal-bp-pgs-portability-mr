#!/usr/bin/env bash
# ============================================================
# submit_build.sh — run the pipeline (or one step) as a SLURM batch job.
#
# Long steps should not run in an interactive shell: they die with the session, and a wall
# limit can kill them halfway. This generates a job script, bakes in the environment the
# modules need, and submits it. Everything after the wrapper's own options is passed
# straight through to run_build.sh.
#
#   bash bin/submit_build.sh --only B06
#   bash bin/submit_build.sh --only B06 --time 8:00:00 --mem 64G --cpus 8
#   bash bin/submit_build.sh                       # full build
#   bash bin/submit_build.sh --from B07 --to B17
#
# Wrapper options (consumed here, NOT passed on): --time --mem --cpus --partition --name
# Everything else goes to run_build.sh verbatim.
#
# NB the job script writes the environment with explicit quoting rather than relying on
# --export=ALL, because the EPI path contains a space ("Epi Data") and has bitten this
# pipeline before.
# ============================================================
set -euo pipefail

PIPE="${MOMI_PIPE:-$(cd "$(dirname "$0")/.." && pwd)}"
RESULTS="${MOMI_RESULTS:-$PIPE/results/current}"
TIME="4:00:00"; MEM="32G"; CPUS="4"; PART="shared"; NAME=""
PASS=()
while [[ $# -gt 0 ]]; do case "$1" in
  --time) TIME="$2"; shift 2;;
  --mem)  MEM="$2";  shift 2;;
  --cpus) CPUS="$2"; shift 2;;
  --partition) PART="$2"; shift 2;;
  --name) NAME="$2"; shift 2;;
  *) PASS+=("$1"); shift;;
esac; done

# name the job after the step if --only was given
if [[ -z "$NAME" ]]; then
  NAME="momi_build"
  for ((i=0; i<${#PASS[@]}; i++)); do
    [[ "${PASS[$i]}" == "--only" ]] && NAME="momi_${PASS[$((i+1))]}"
  done
fi

SLURMDIR="$RESULTS/logs/slurm"; mkdir -p "$SLURMDIR"
JOB="$SLURMDIR/${NAME}_$(date +%Y%m%d-%H%M).sbatch"

{
  echo "#!/usr/bin/env bash"
  echo "#SBATCH --job-name=$NAME"
  echo "#SBATCH --partition=$PART"
  echo "#SBATCH --time=$TIME"
  echo "#SBATCH --mem=$MEM"
  echo "#SBATCH --cpus-per-task=$CPUS"
  echo "#SBATCH --output=$SLURMDIR/${NAME}_%j.out"
  echo "#SBATCH --error=$SLURMDIR/${NAME}_%j.err"
  echo
  echo "set -o pipefail"
  echo "echo \"host=\$(hostname) job=\$SLURM_JOB_ID started=\$(date '+%F %T')\""
  echo
  echo "# the batch shell is fresh: run_build.sh deliberately does NOT load modules itself"
  echo "# (lmod's shell function aborts non-interactive scripts), so load them here."
  echo "#"
  echo "# NOTE: NO 'set -u' anywhere near this. conda_R/4.3's activation script references"
  echo "# an unbound SYS_SYSROOT, so under 'set -u' the module load dies instantly and takes"
  echo "# the whole job with it (job 34258691 failed in 2 s exactly this way). Same family as"
  echo "# the documented lmod-vs-non-interactive-shell problem. Belt and braces below."
  echo "set +u"
  echo "module load conda_R/4.3 plink/2.00a4.6"
  echo "set +u   # keep it off: some module activation runs lazily on first use"
  echo
  echo "# explicit quoting, not --export=ALL: the EPI path contains a space."
  printf 'export MOMI_PIPE=%q\n'    "$PIPE"
  printf 'export MOMI_RESULTS=%q\n' "$RESULTS"
  [[ -n "${EPI:-}"          ]] && printf 'export EPI=%q\n'          "$EPI"
  [[ -n "${SSC:-}"          ]] && printf 'export SSC=%q\n'          "$SSC"
  [[ -n "${MOMI_OUTROOT:-}" ]] && printf 'export MOMI_OUTROOT=%q\n' "$MOMI_OUTROOT"
  [[ -n "${MOMI_SCORES:-}"  ]] && printf 'export MOMI_SCORES=%q\n'  "$MOMI_SCORES"
  echo
  printf 'bash %q' "$PIPE/run_build.sh"
  for a in ${PASS[@]+"${PASS[@]}"}; do printf ' %q' "$a"; done
  echo
  echo "echo \"finished=\$(date '+%F %T') exit=\$?\""
} > "$JOB"

chmod +x "$JOB"
echo "job script: $JOB"
JID=$(sbatch --parsable "$JOB")
echo "submitted:  $JID"
echo
echo "  watch:    tail -f $SLURMDIR/${NAME}_${JID}.out"
echo "  status:   squeue -j $JID"
echo "  cancel:   scancel $JID"
echo "  module log: $RESULTS/logs/<ID>.log"
