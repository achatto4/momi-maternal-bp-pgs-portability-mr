#!/usr/bin/env bash
# ============================================================
# merge_chr_dosage_and_rename.sh   — DOSAGE variant of merge_chr_and_rename.sh
# Concatenate the 22 per-chromosome PGEN filesets from lpwgs_chr_to_plink_dosage.sh
# (same samples, disjoint variants) into one pgen, then rename samples -> Epi
# PARTICIPANT_ID. All plink2 (pgen merge + --update-ids), dosages preserved end-to-end.
#
# Usage:
#   merge_chr_dosage_and_rename.sh --prefix <out-prefix> --rename RENAMEFILE --out FINAL
#       [--plink2 plink2] [--threads 4] [--mem-mb 28000]
# Expects <prefix>_chr1 .. <prefix>_chr22 (pgen) ; writes FINAL.{pgen,pvar,psam}
# ============================================================
set -euo pipefail
PREFIX=""; RENAME=""; OUT=""; PLINK2=plink2; THREADS=4; MEMMB=28000
while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) PREFIX="$2"; shift 2;;
    --rename) RENAME="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --plink2) PLINK2="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mem-mb) MEMMB="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
for v in PREFIX RENAME OUT; do [[ -n "${!v}" ]] || { echo "missing --$v"; exit 1; }; done
command -v "$PLINK2" >/dev/null || { echo "plink2 not on PATH (module load?)" >&2; exit 3; }
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
mkdir -p "$(dirname "$OUT")"

# collect the per-chr pgen filesets that exist
LIST="$OUT.chrmerge.list"; : > "$LIST"; n=0
for chr in $(seq 1 22); do
  if [[ -f "${PREFIX}_chr${chr}.pgen" ]]; then echo "${PREFIX}_chr${chr}" >> "$LIST"; n=$((n+1)); fi
done
log "found $n per-chromosome PGEN filesets"
[[ "$n" -ge 1 ]] || { echo "ERROR: no per-chr pgen at ${PREFIX}_chr*" >&2; exit 4; }

TMP="$OUT.allchr"
if [[ "$n" -eq 1 ]]; then
  pfx="$(head -1 "$LIST")"; for e in pgen pvar psam; do cp "$pfx.$e" "$TMP.$e"; done
else
  log "pmerging $n chromosomes (same samples, union of variants; dosages kept)"
  # default --pmerge-list mode is pfile (bare prefixes -> pgen+pvar+psam); --pmerge-list
  # writes a pgen to the --out prefix by default, so no separate --make-pgen.
  "$PLINK2" --pmerge-list "$LIST" \
    --threads "$THREADS" --memory "$MEMMB" \
    --out "$TMP" > "$OUT.chrmerge.log" 2>&1 \
    || { echo "chr pmerge failed; see $OUT.chrmerge.log"; tail -10 "$OUT.chrmerge.log"; exit 5; }
fi
log "merged (pre-rename): $(grep -cv '^#' "$TMP.pvar") variants x $(grep -cv '^#' "$TMP.psam") samples"

log "renaming samples -> PARTICIPANT_ID (plink2 --update-ids)"
"$PLINK2" --pfile "$TMP" --update-ids "$RENAME" --make-pgen \
  --threads "$THREADS" --memory "$MEMMB" --out "$OUT" > "$OUT.rename.log" 2>&1 \
  || { echo "rename failed; see $OUT.rename.log"; tail -10 "$OUT.rename.log"; exit 6; }

[[ -f "$OUT.pgen" ]] || { echo "ERROR: final pgen not produced" >&2; exit 7; }
log "DONE: $(grep -cv '^#' "$OUT.pvar") variants x $(grep -cv '^#' "$OUT.psam") samples -> $OUT.{pgen,pvar,psam}"
rm -f "$TMP".{pgen,pvar,psam} "$LIST"
log "sample IIDs (first 3):"; awk '!/^#/{print $1; c++} c>=3{exit}' "$OUT.psam"
