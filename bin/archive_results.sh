#!/usr/bin/env bash
# ============================================================
# archive_results.sh — snapshot a completed build into results/archive/.
#
# Every full run is frozen into its own directory, named by timestamp and git SHA, so that
# a number quoted in the manuscript can always be traced back to the exact build that
# produced it. Nothing in results/current/ is ever the record of anything -- it is scratch
# that the next run overwrites.
#
#   bash bin/archive_results.sh                 # snapshot tables/figures/logs/manifest
#   bash bin/archive_results.sh --with-intermediates   # also freeze the .rds intermediates
#   bash bin/archive_results.sh --label "part-I-locked"
#
# The archive is APPEND-ONLY by convention: never edit or delete a snapshot. If numbers
# change, take a new one -- the diff between snapshots is the audit trail.
# ============================================================
set -euo pipefail

PIPE="${MOMI_PIPE:-$(cd "$(dirname "$0")/.." && pwd)}"
RESULTS="${MOMI_RESULTS:-$PIPE/results/current}"
ARCHIVE="$PIPE/results/archive"
WITH_INT=0; LABEL=""
while [[ $# -gt 0 ]]; do case "$1" in
  --with-intermediates) WITH_INT=1; shift;;
  --label) LABEL="$2"; shift 2;;
  *) echo "unknown option: $1"; exit 2;;
esac; done

[[ -f "$RESULTS/manifest.tsv" ]] || { echo "no manifest at $RESULTS — run the build first"; exit 1; }

SHA="$(git -C "$PIPE" rev-parse --short HEAD 2>/dev/null || echo nogit)"
DIRTY=""; git -C "$PIPE" diff --quiet 2>/dev/null || DIRTY="-dirty"
STAMP="$(date +%Y%m%d-%H%M)"
NAME="${STAMP}_${SHA}${DIRTY}${LABEL:+_$LABEL}"
DEST="$ARCHIVE/$NAME"

mkdir -p "$DEST"
for d in tables figures logs qc; do
  [[ -d "$RESULTS/$d" ]] && cp -R "$RESULTS/$d" "$DEST/" || true
done
cp "$RESULTS/manifest.tsv" "$DEST/" 2>/dev/null || true
cp "$RESULTS/PROVENANCE.txt" "$DEST/" 2>/dev/null || true
[[ $WITH_INT -eq 1 && -d "$RESULTS/intermediates" ]] && cp -R "$RESULTS/intermediates" "$DEST/"

# a README so a snapshot is self-describing years later
{
  echo "# MOMI results snapshot: $NAME"
  echo
  echo "created:      $(date '+%F %T')"
  echo "git SHA:      $SHA${DIRTY:+  (WORKING TREE DIRTY — not reproducible from git alone)}"
  echo "pipeline:     $PIPE"
  echo "source:       $RESULTS"
  echo "intermediates included: $([[ $WITH_INT -eq 1 ]] && echo yes || echo no)"
  [[ -n "$LABEL" ]] && echo "label:        $LABEL"
  echo
  echo "## items in this build"
  echo
  awk -F'\t' 'NR==1{next} {printf "- %-24s %-6s %s\n", $2, $4, $6}' "$RESULTS/manifest.tsv"
} > "$DEST/README.md"

echo "archived -> $DEST"
[[ -n "$DIRTY" ]] && echo "WARNING: working tree was dirty; this snapshot is not reproducible from git alone."
ls -1 "$ARCHIVE" | tail -5 | sed 's/^/  recent: /'
