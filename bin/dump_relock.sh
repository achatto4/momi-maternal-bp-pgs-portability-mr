#!/usr/bin/env bash
# ============================================================
# dump_relock.sh — print every table needed to re-verify B09-B17 in ONE paste.
#
# The re-lock is a text problem, not a compute problem: B09-B17 were all rebuilt after B01's
# antenatal restriction and are marked "rebuilt, not re-verified". Three of them carry PROSE
# that now contradicts their own numbers (S16's conclusion reversed, S7's framing inverted,
# F3's panel C caption). Verifying that requires reading the tables, not rerunning anything.
#
#   bash bin/dump_relock.sh | tee relock_dump.txt
# ============================================================
set -uo pipefail
R="${MOMI_RESULTS:?set MOMI_RESULTS}"
T="$R/tables"

show(){  # show <file> <n_lines> <label>
  local f="$T/$1"; shift; local n="$1"; shift
  echo ""
  echo "=============================================================="
  echo "== $* "
  echo "== file: $(basename "$f")"
  echo "=============================================================="
  if [[ -f "$f" ]]; then
    head -n "$n" "$f" | column -t -s $'\t'
    local tot; tot=$(( $(wc -l < "$f") - 1 ))
    echo "... [$tot data rows total]"
  else
    echo "MISSING: $f"
  fi
}

echo "MOMI re-lock dump — $(date '+%Y-%m-%d %H:%M:%S')"
echo "results root: $R"
echo "git: $(git -C "${MOMI_PIPE:-.}" rev-parse --short HEAD 2>/dev/null || echo NA)"

show S4_transfer_grid.tsv 25 "B09 / S4 — full transferability grid (8 defs x 5 cohorts x 8 scores)"
show T2_transfer.tsv      12 "B10 / T2 — chosen instrument, OBSERVED R2 only"
show S5_platform.tsv      12 "B12 / S5 — GSA vs lpWGS dosage"
show S7_residual.tsv      12 "B13 / S7 — GA-residual vs mean  [FRAMING INVERTED: resid is now PRIMARY]"
show S8_wk20.tsv          12 "B14 / S8 — <20wk vs >=20wk"
show S16_bpdist.tsv       20 "B15 / S16 — BP distributions  [CONCLUSION REVERSED: gap persists at every def]"
show T3_instrument.tsv     6 "B17 / T3 — instrument decision"

echo ""
echo "=============================================================="
echo "== B11 / F2 and B16 / F3 are FIGURES — key numbers from the manifest"
echo "=============================================================="
awk -F'\t' 'NR==1 || $2 ~ /^(F2_transfer|F3_diagnostic|S4_transferability|S16_bpdist|S7_residual)$/' \
    "$R/manifest.tsv" | cut -f1,2,5,6 | column -t -s $'\t' | tail -20

echo ""
echo "=============================================================="
echo "== manifest: current status of every item"
echo "=============================================================="
Rscript "${MOMI_PIPE:?}/bin/check_stale.R" 2>/dev/null | head -20
