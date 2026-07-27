#!/usr/bin/env bash
# ============================================================
# run_build.sh — single synchronized driver for the MOMI BP->perinatal paper.
# Runs every deliverable module in DEPENDENCY ORDER (build order B01..B32, which is
# NOT manuscript order), each into its own log, all reporting to one manifest.tsv.
# A module reads only the frozen intermediates from the build steps, so no two
# deliverables can disagree on the data. Planned-but-unwritten modules are recorded
# as SKIP so the manifest always shows the full intended progression.
#
# Usage:
#   bash run_build.sh --epi EPI --sscore-dir SSC [options]
#   options:
#     --only  ID        run just one step (e.g. --only B03)
#     --from  ID        run from this step to the end
#     --to    ID        run up to and including this step
#     --list            print the ordered registry and exit
#     --dry-run         show what would run, don't execute
#     --pipe  DIR       pipeline root (default: script dir)
#     --results DIR     results root (default: <pipe>/results/current)
#     --no-modules      skip `module load` (already loaded)
# ============================================================
set -uo pipefail

PIPE="$(cd "$(dirname "$0")" && pwd)"
# inherit EPI/SSC from the environment if already exported; a --epi/--sscore-dir flag
# overrides. (Do NOT blank them here — that would clobber an exported EPI to empty.)
EPI="${EPI:-}"; SSC="${SSC:-}"
ONLY=""; FROM=""; TO=""; LIST=0; DRY=0; NOMOD=0
RESULTS="$PIPE/results/current"
while [[ $# -gt 0 ]]; do case "$1" in
  --epi) EPI="$2"; shift 2;;
  --sscore-dir) SSC="$2"; shift 2;;
  --only) ONLY="$2"; shift 2;;
  --from) FROM="$2"; shift 2;;
  --to) TO="$2"; shift 2;;
  --list) LIST=1; shift;;
  --dry-run) DRY=1; shift;;
  --pipe) PIPE="$2"; shift 2;;
  --results) RESULTS="$2"; shift 2;;
  --no-modules) NOMOD=1; shift;;
  *) echo "unknown option: $1"; exit 2;;
esac; done

export MOMI_PIPE="$PIPE"
export MOMI_RESULTS="$RESULTS"
# Pass EPI/SSC to R modules via the ENVIRONMENT (robust to spaces in paths, e.g.
# ".../Epi Data/..."). Only export if provided on the CLI, else inherit the shell's.
[[ -n "$EPI" ]] && export EPI
[[ -n "$SSC" ]] && export SSC
mkdir -p "$RESULTS/logs" "$RESULTS/intermediates" "$RESULTS/tables" "$RESULTS/figures"
MANIFEST="$RESULTS/manifest.tsv"

# ---- ordered registry: ID | script (relative to PIPE) | extra-args | policy ----
# policy: stop = abort run on failure; cont = log FAIL/SKIP and continue (Stage 7 etc.)
# NB: EPI/SSC reach modules via the environment (exported above), NOT as CLI args, so
# paths with spaces are safe. Modules read momi_arg("--epi", Sys.getenv("EPI")).
REGISTRY=(
  # ---- Stage 0: foundations (frozen intermediates) ----
  "B01|01_phenotypes/build_analytic.R||stop"
  "B02|05_prs/build_prs.R||stop"
  "B03|05_prs/build_transfer_grid.R||stop"
  # ---- Stage 1: static metadata ----
  "B04|05_prs/deliv_S1_prs_manifest.R||cont"
  # ---- Stage 2: genotype QC (needs genotype paths; stubbed to SKIP) ----
  "B05|06_phase2/deliv_S2_variant_qc.R||cont"
  "B06|06_phase2/deliv_S3_scoresnp.R||cont"
  # B06b: added 2026-07-19. Build step, not a deliverable — produces intermediates/pcs.rds.
  # WITHIN-cohort PCs are needed to adjust the MR arm for population stratification (Burgess
  # guidelines, ref/methods_citations.tsv); 1000G-PROJECTED PCs are needed for SF1's ancestry
  # figure and SF2's PC-distance. Must run before B22.
  "B06b|06_phase2/build_pcs.R||cont"
  # ---- Stage 3: sample description ----
  "B07|01_phenotypes/deliv_F1_flow.R||stop"
  "B08|01_phenotypes/deliv_T1_cohorts.R||stop"
  # ---- Stage 4: Part I transferability (views of the grid) ----
  "B09|05_prs/deliv_S4_grid_table.R||stop"
  "B10|05_prs/deliv_T2_transfer.R||stop"
  "B11|05_prs/deliv_F2_transfer.R||stop"
  "B12|05_prs/deliv_S5_platform.R||cont"
  "B13|05_prs/deliv_S7_residual.R||cont"
  "B14|05_prs/deliv_S8_wk20.R||cont"
  "B15|05_prs/deliv_S16_bpdist.R||cont"
  "B16|05_prs/deliv_F3_diagnostic.R||cont"
  "B17|05_prs/deliv_T3_instrument.R||stop"
  # ---- Stage 5: Part II observational ----
  "B18|05_prs/deliv_S9_confounding.R||stop"
  # B18b: added 2026-07-19. S9 produced per-cohort ORs spanning 0.94-1.70 with two cohorts
  # null on all six outcome x trait estimates. T4/F4 and the MR panel all assume a single
  # pooled effect, so whether pooling is defensible has to be settled BEFORE they are built.
  "B18b|05_prs/deliv_S9b_heterogeneity.R||cont"
  "B19|05_prs/deliv_T5_subtype.R||stop"
  "B20|05_prs/deliv_S6_downstream.R||cont"
  # ---- Stage 6: Part II MR ----
  # ORDERING FIX 2026-07-23: T4 (triangulation) READS mr_panel/mr_pooled/power_grid, which are
  # produced by S10 (B22) and S13 (B25). It was previously listed as B21, i.e. it RAN BEFORE
  # its own inputs, so in any single build sweep it consumed the PREVIOUS run's MR panel. In
  # steady state the two matched and this was invisible; the moment the instrument changed
  # (ancestry-matched scores), T4 lagged one build behind and reported the old numbers while
  # S10b reported the new ones. T4 is therefore run AFTER B25 below (kept its id "B21" so
  # --only B21 still resolves; only its position in the sweep changed).
  "B22|05_prs/deliv_S10_mrpanel.R||cont"
  # B23 (leave-one-out sensitivity) REMOVED FROM THE BUILD 2026-07-20 (Anagh).
  # The per-cohort forest (B27/SF4) shows the same thing more directly: a reader can see that
  # PreSSMat's interval is the only one excluding 1 for indicated PTB, or that four of five
  # cohorts sit below zero for birthweight, without needing a sensitivity framework explained.
  # The script is retained (not deleted) because it is the RECORD OF WHY certain results are
  # described the way they are -- SBP->LBW is reported as Karachi's rather than pooled, and
  # indicated PTB as Matlab's rather than a five-cohort effect. Those descriptions must stand
  # whether or not the analysis behind them appears. Re-run manually if ever queried:
  #   bash run_build.sh --only B23     (still resolvable by id; just not in the default list)
  # "B23|05_prs/deliv_S11_sensitivity.R||cont"
  "B24|05_prs/deliv_S12_controls.R||cont"
  "B25|05_prs/deliv_S13_power.R||cont"
  # T4 runs here, after its inputs (S10=B22 mr_panel/mr_pooled, S13=B25 power_grid) exist.
  "B21|05_prs/deliv_T4_triangulation.R||stop"
  # B26 (non-linear / doubly-ranked MR) REMOVED 2026-07-20 (Anagh) -- no power at this N.
  # "B26|05_prs/deliv_F5_nonlinear.R||cont"
  # SF4 became load-bearing on 2026-07-20 when S11 (leave-one-out) was cut from the paper:
  # the per-cohort intervals show the same concentration directly. Runs after B19/B22.
  "B27|05_prs/deliv_SF4_persite.R||cont"
  # ---- Stage 7: Phase-2 genetics (blocked on --merged/--analysis-dir; SKIP-safe) ----
  "B28|06_phase2/deliv_SF1_pca.R||cont"
  "B29|06_phase2/deliv_SF2_distance.R||cont"
  "B30|06_phase2/deliv_S17_pcadj.R||cont"
  # B33 (S20, instrument-choice sensitivity) was written 2026-07-23 and then dropped by Anagh
  # before it ran -- not in the build. Script retained in 06_phase2/deliv_S20_instrument_sens.R
  # if the EUR-vs-SAS instrument comparison is ever wanted.
  # "B33|06_phase2/deliv_S20_instrument_sens.R||cont"
  # B31 (S14, MVMR with BMI) and B32 (S15, fetal-vs-maternal decomposition) REMOVED FROM THE
  # BUILD 2026-07-21 (Anagh). The scripts are retained -- they run cleanly and are the record
  # of the analyses -- but they are no longer part of the paper:
  #   S14: the MVMR-BMI adjustment is dropped from the analysis plan.
  #   S15: the fetal decomposition conflicted with the far-better-powered Warrington 2019 and
  #        was already being reported only as a limitation; removed at the analysis stage.
  # Run manually if ever queried:  bash run_build.sh --only B31   (still resolvable by id)
  # "B31|06_phase2/deliv_S14_mvmr.R||cont"
  # "B32|06_phase2/deliv_S15_fetal.R||cont"
)

if [[ $LIST -eq 1 ]]; then
  printf "%-5s %-45s %-8s %s\n" ID SCRIPT POLICY STATUS
  for e in "${REGISTRY[@]}"; do IFS='|' read -r id sc ar po <<<"$e"
    st="planned"; [[ -f "$PIPE/$sc" ]] && st="implemented"
    printf "%-5s %-45s %-8s %s\n" "$id" "$sc" "$po" "$st"; done
  exit 0
fi

# ---- selection window over the ordered list ----
ids=(); for e in "${REGISTRY[@]}"; do ids+=("${e%%|*}"); done
pos() { local t="$1" i=0; for x in "${ids[@]}"; do [[ "$x" == "$t" ]] && { echo $i; return; }; i=$((i+1)); done; echo -1; }
FROM_I=0; TO_I=$(( ${#ids[@]} - 1 ))
[[ -n "$FROM" ]] && FROM_I=$(pos "$FROM")
[[ -n "$TO"   ]] && TO_I=$(pos "$TO")

# ---- fresh manifest + provenance only for a full run ----
# ---- manifest schema, declared ONCE ----
# This MUST match the column order built in momi_manifest_append() (lib/momi_io.R). They are
# written by different languages, so they drifted the first time columns were added: R wrote
# 13 fields while this header and skip_row() still wrote 11, and fread() then failed to parse
# the file. Keeping the list in one variable and deriving skip_row's padding from it means a
# future column can only be added in one place here.
MANIFEST_COLS="ts id script status n key inputs outputs in_md5 out_md5 seconds git message"
NCOL=$(echo $MANIFEST_COLS | wc -w)

if [[ -z "$ONLY" && -z "$FROM" && -z "$TO" ]]; then
  echo "$MANIFEST_COLS" | tr ' ' '\t' > "$MANIFEST"
  { echo "MOMI paper build — $(date '+%F %T')";
    echo "pipe=$PIPE"; echo "results=$RESULTS";
    echo "epi=$EPI"; echo "sscore=$SSC";
    echo "git=$(git -C "$PIPE" rev-parse --short HEAD 2>/dev/null)"; } > "$RESULTS/PROVENANCE.txt"
fi

# NOTE: load environment modules in your SHELL before running:
#   module load conda_R/4.3 plink/2.00a4.6
# The driver does NOT call `module load` itself — the lmod `module` function can abort a
# non-interactive script. --no-modules is accepted for backward compat and ignored.
: "${NOMOD:=0}"
command -v Rscript >/dev/null 2>&1 || { echo "Rscript not on PATH — run 'module load conda_R/4.3' first"; exit 3; }

skip_row() { # id script message
  # Row layout: ts, id, script, status, <NCOL-6 empty fields>, git, message.
  # The padding is derived from MANIFEST_COLS so it cannot drift out of step with the header.
  local blanks="" i
  for ((i=0; i<NCOL-6; i++)); do blanks="${blanks}\t"; done
  printf "%s\t%s\t%s\tSKIP${blanks}\t%s\t%s\n" \
    "$(date '+%F %T')" "${1:-}" "${2:-}" "$(git -C "$PIPE" rev-parse --short HEAD 2>/dev/null)" "${3:-}" >> "$MANIFEST"
}

# ---- concurrency lock ----
# Two builds sharing results/current/ will interleave. A full run TRUNCATES the manifest and
# rewrites it as it goes, so a concurrent --only appends into a half-written file; worse, it
# may read an intermediate while the other run is rewriting it and silently produce results
# from a partial .rds. Observed 2026-07-19: a --only B25 during a full rebuild reported
# "items OK: 6" because the full run had only reached its sixth module.
LOCK="$RESULTS/.build.lock"
if [[ -f "$LOCK" ]]; then
  OTHER=$(head -1 "$LOCK" 2>/dev/null)
  if [[ -n "$OTHER" ]] && kill -0 "$OTHER" 2>/dev/null; then
    echo "REFUSING TO START: another build is running (pid $OTHER, started $(sed -n 2p "$LOCK"))."
    echo "  wait for it, or if you are certain it is dead:  rm $LOCK"
    exit 4
  fi
  echo "note: stale lock from dead pid $OTHER — removing"; rm -f "$LOCK"
fi
{ echo $$; date '+%F %T'; echo "${ONLY:-${FROM:-full}}"; } > "$LOCK"
trap 'rm -f "$LOCK"' EXIT

echo "== MOMI build == pipe=$PIPE results=$RESULTS"
i=-1
for e in "${REGISTRY[@]}"; do
  i=$((i+1))
  IFS='|' read -r id sc ar po <<<"$e"
  [[ -n "$ONLY" && "$id" != "$ONLY" ]] && continue
  if [[ -z "$ONLY" ]]; then (( i < FROM_I || i > TO_I )) && continue; fi
  args="${ar//@EPI@/$EPI}"; args="${args//@SSC@/$SSC}"
  log="$RESULTS/logs/$id.log"
  if [[ ! -f "$PIPE/$sc" ]]; then
    echo "  [$id] SKIP (planned, not yet implemented): $sc"
    skip_row "$id" "$sc" "planned: module not yet written"; continue
  fi
  echo "  [$id] $sc $args  -> $log"
  [[ $DRY -eq 1 ]] && continue
  if Rscript "$PIPE/$sc" $args >"$log" 2>&1; then
    echo "     done ($id) — tail:"; tail -n 2 "$log" | sed 's/^/       /'
  else
    echo "     FAIL ($id) — see $log"; tail -n 8 "$log" | sed 's/^/       /'
    [[ "$po" == "stop" ]] && { echo "aborting: $id is stop-on-error"; exit 1; }
  fi
done
echo "== manifest: $MANIFEST"

# ---- after a FULL run: snapshot into results/archive/ and report staleness ----
# results/current/ is scratch and gets overwritten; the archive is the record. A full run is
# the only thing worth freezing, since a partial run leaves the tree internally inconsistent.
if [[ -z "$ONLY" && -z "$FROM" && -z "$TO" && $DRY -eq 0 ]]; then
  bash "$PIPE/bin/archive_results.sh" || echo "(archive step failed — results/current is intact)"
fi
if [[ $DRY -eq 0 ]] && command -v Rscript >/dev/null 2>&1; then
  echo; Rscript "$PIPE/bin/check_stale.R" || true
fi
