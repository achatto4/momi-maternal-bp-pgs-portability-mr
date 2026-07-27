#!/usr/bin/env bash
# ============================================================
# attic_move.sh — move superseded / exploratory scripts into attic/.
#
# WHY MOVE RATHER THAN DELETE. Git history preserves both equally, so deletion buys nothing.
# A move keeps the files visible: if a number in the manuscript later turns out to trace back
# to one of these, it is still on disk rather than requiring an archaeology session in the
# reflog. The working tree gets clean; nothing becomes unrecoverable.
#
# HOW THE LIST WAS BUILT. Transitive closure of file references starting from the pipeline
# entry points (run_build.sh, the preprocessing drivers, bin/*). Anything not reachable is a
# candidate. Then four categories were protected BY HAND, because reachability cannot judge
# them:
#
#   push.sh                       -- the sync helper this workflow depends on
#   build_common_union.sh         -- produced the KEEP sets the real filesets were built from
#   build_keep_sets.sh            -- ditto; provenance for variant selection
#   pgs_catalog_manifest.py       -- generated ref/pgs_catalog_metadata.tsv, which we ship
#   audit_B01.R / audit_B02.R     -- the from-scratch verification of B01/B02
#   deliv_S11_sensitivity.R       -- cut from the paper, but the record of WHY results are
#                                    described as they are (SBP->LBW is Karachi's, indicated
#                                    PTB is Matlab's). Those descriptions stand regardless.
#   deliv_F5_nonlinear.R          -- cut for lack of power; keep the reasoning discoverable
#
# GUARD. The script refuses to move anything referenced by the build registry in
# run_build.sh. If that guard ever trips, the list is wrong -- investigate, do not override.
#
#   bash bin/attic_move.sh            # from the repo root or anywhere
#   bash bin/attic_move.sh --dry-run
# ============================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PIPE="$(dirname "$HERE")"                 # .../bp_ptb_pipeline
ROOT="$(cd "$PIPE/.." && pwd)"            # repo root (contains .git)
REL="$(basename "$PIPE")"                 # bp_ptb_pipeline
DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1

cd "$ROOT" || exit 1
[[ -d .git ]] || { echo "no .git at $ROOT" >&2; exit 1; }

FILES=(
00_preprocess/_hardcall_lpwgs.sh
00_preprocess/_mega_merge_one.sh
00_preprocess/check_instrument_recovery.sh
00_preprocess/check_score_in_raw_vcf.sh
00_preprocess/diagnose_lpwgs_build.sh
00_preprocess/diagnose_lpwgs_coverage.sh
00_preprocess/diagnose_lpwgs_imputation.sh
00_preprocess/lpwgs_evidence_report.sh
00_preprocess/mega_shared_count.sh
00_preprocess/report_nopass.sh
01_phenotypes/check_sites.R
01_phenotypes/epi_eda_ptb.R
01_phenotypes/epi_extra_outcomes.R
01_phenotypes/epi_ptb_subtype.R
01_phenotypes/explore_bd_cohorts.R
01_phenotypes/extra_figs_data.R
01_phenotypes/obs_confounders.R
01_phenotypes/ptb_by_bp_bins.R
01_phenotypes/ptb_ushape_siteadj.R
01_phenotypes/sample_flow.R
01_phenotypes/sensitivity_codebook.R
01_phenotypes/site_outcomes.R
05_prs/bp_def_downstream.R
05_prs/bp_definitions.R
05_prs/bp_dist_bd.R
05_prs/bp_var_traj.R
05_prs/confirm_nopass_prs.sh
05_prs/cv_power.R
05_prs/exposure_primary.R
05_prs/make_figures.R
05_prs/mr_panel_extras.R
05_prs/prs_bestinstrument_mr.R
05_prs/prs_ebbp.R
05_prs/prs_extra_outcomes.R
05_prs/prs_meanbp_wald.R
05_prs/prs_mr.R
05_prs/prs_mr_meta.R
05_prs/prs_nonlinear_ptb.R
05_prs/prs_outcome_meta.R
05_prs/prs_site_outcomes.R
05_prs/prs_to_pe_meta.R
05_prs/transfer_diag_bd.R
05_prs/transferability_merged.R
05_prs/transferability_table.R
05_prs/verify_pipeline.R
06_phase2/a1_fst_pcdist.R
06_phase2/a1_fst_pcdist.sh
06_phase2/a2_scoresnp_freq.sh
06_phase2/a2_scoresnp_summary.R
06_phase2/a3_king.sh
06_phase2/a3_pcadj_r2.R
06_phase2/c1_bmi_mvmr.R
06_phase2/c1_bmi_score.sh
06_phase2/c2_fetal_maternal.R
06_phase2/c2_fetal_score.sh
06_phase2/run_phase2.sh
run_build_keep.sh
run_lpwgs_nopass_test.sh
run_mega_build.sh
run_paper.sh
run_prs_mr_gsa.sh
)

## ---- guard: nothing in the build registry may appear above ----
BUILD=$(grep -oE '^[[:space:]]*"B[0-9]+[a-z]?\|[^|]+\.R' "$PIPE/run_build.sh" | sed 's/.*|//' | sort -u)
bad=0
for f in "${FILES[@]}"; do
  if grep -qxF "$f" <<< "$BUILD"; then echo "ABORT: $f is a live build module" >&2; bad=1; fi
done
[[ $bad -eq 1 ]] && { echo "list is wrong — investigate, do not override" >&2; exit 1; }
echo "guard passed: none of the ${#FILES[@]} files is a live build module"

if [[ -f .git/index.lock ]]; then
  echo "NOTE: .git/index.lock exists. If no git process is running, remove it:"
  echo "      rm -f '$ROOT/.git/index.lock'"
  exit 1
fi

moved=0; skipped=0
for f in "${FILES[@]}"; do
  src="$REL/$f"; dst="$REL/attic/$f"
  if [[ ! -f "$src" ]]; then echo "  absent: $f"; skipped=$((skipped+1)); continue; fi
  if [[ $DRY -eq 1 ]]; then echo "  would move: $f"; moved=$((moved+1)); continue; fi
  mkdir -p "$(dirname "$dst")"
  if git mv "$src" "$dst"; then moved=$((moved+1)); else echo "  FAILED: $f"; fi
done

if [[ $DRY -eq 1 ]]; then echo "(dry run) $moved would move, $skipped absent"; exit 0; fi

cat > "$REL/attic/README.md" <<'EOF'
# attic

Superseded and exploratory code, kept rather than deleted.

Nothing here runs. Every one of these files is either **replaced by a module in the current
build** or was a one-off exploration. They are retained because a number in the manuscript
may later need tracing to its origin, and that is easier with the file on disk than in the
reflog.

Direct replacements, for anyone looking for the current version:

| attic file | replaced by |
|---|---|
| `06_phase2/c1_bmi_mvmr.R`, `c1_bmi_score.sh` | `06_phase2/deliv_S14_mvmr.R` (B31) |
| `06_phase2/c2_fetal_maternal.R`, `c2_fetal_score.sh` | `06_phase2/deliv_S15_fetal.R` (B32) |
| `06_phase2/a1_fst_pcdist.*` | `06_phase2/deliv_S3_scoresnp.R` (B06), `deliv_SF2_distance.R` (B29) |
| `06_phase2/a2_scoresnp_*` | `06_phase2/deliv_S3_scoresnp.R` (B06) |
| `06_phase2/a3_pcadj_r2.R` | `06_phase2/deliv_S17_pcadj.R` (B30) |
| `05_prs/prs_meanbp_wald.R`, `prs_mr.R`, `prs_site_outcomes.R`, `mr_panel_extras.R` | `05_prs/deliv_S10_mrpanel.R` (B22) |
| `05_prs/transferability_table.R`, `transferability_merged.R` | `05_prs/build_transfer_grid.R` (B03), `deliv_T2_transfer.R` (B10) |
| `05_prs/cv_power.R` | `05_prs/deliv_S13_power.R` (B25) |
| `05_prs/bp_def_downstream.R` | `05_prs/deliv_S6_downstream.R` (B20) |
| `01_phenotypes/sample_flow.R` | `01_phenotypes/deliv_F1_flow.R` (B07) |
| `01_phenotypes/obs_confounders.R` | `05_prs/deliv_S9_confounding.R` (B18) |

The authoritative index of what produces each figure and table is
`docs/figure_provenance.md`; the build order is the registry in `run_build.sh`.
EOF

git add "$REL/attic/README.md"
echo
echo "moved $moved files (skipped $skipped absent)"
echo "review with:  git status --short"
echo "then:         git commit -m 'move superseded and exploratory scripts to attic/'"
