#!/usr/bin/env bash
# ============================================================
# merge_chr_and_rename.sh
# Concatenate the 22 per-chromosome beds from lpwgs_chr_to_plink.sh (same samples,
# disjoint variants) into one fileset, then rename samples -> Epi PARTICIPANT_ID.
# All plink1.9 (merge + --update-ids), so no plink module conflict.
#
# Usage:
#   merge_chr_and_rename.sh --prefix <out-prefix> --rename RENAMEFILE --out FINAL
#       [--plink1 plink] [--threads 4] [--mem-mb 28000]
# Expects <prefix>_chr1 .. <prefix>_chr22 ; writes FINAL.{bed,bim,fam}
# ============================================================
set -euo pipefail
PREFIX=""; RENAME=""; OUT=""; PLINK1=plink; THREADS=4; MEMMB=28000
while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) PREFIX="$2"; shift 2;;
    --rename) RENAME="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --plink1) PLINK1="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mem-mb) MEMMB="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
for v in PREFIX RENAME OUT; do [[ -n "${!v}" ]] || { echo "missing --$v"; exit 1; }; done
command -v "$PLINK1" >/dev/null || { echo "plink (1.9) not on PATH; module load plink/1.90b" >&2; exit 3; }
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
WORK="$(dirname "$OUT")"; mkdir -p "$WORK"

# collect the per-chr beds that exist
LIST="$OUT.chrmerge.list"; : > "$LIST"; n=0
for chr in $(seq 1 22); do
  [[ -f "${PREFIX}_chr${chr}.bed" ]] && { echo "${PREFIX}_chr${chr}"; n=$((n+1)); }
done > "$LIST"
log "found $n per-chromosome filesets"
[[ "$n" -ge 1 ]] || { echo "ERROR: no per-chr beds at ${PREFIX}_chr*" >&2; exit 4; }

TMP="$OUT.allchr"
if [[ "$n" -eq 1 ]]; then
  cp "$(head -1 "$LIST").bed" "$TMP.bed"; cp "$(head -1 "$LIST").bim" "$TMP.bim"; cp "$(head -1 "$LIST").fam" "$TMP.fam"
else
  log "merging $n chromosomes (same samples, union of variants)"
  "$PLINK1" --merge-list "$LIST" --make-bed --threads "$THREADS" --memory "$MEMMB" \
    --out "$TMP" > "$OUT.chrmerge.log" 2>&1 \
    || { echo "chr-merge failed; see $OUT.chrmerge.log"; tail -8 "$OUT.chrmerge.log"; exit 5; }
fi
log "merged (pre-rename): $(wc -l < "$TMP.bim") variants x $(wc -l < "$TMP.fam") samples"

log "renaming samples -> PARTICIPANT_ID"
"$PLINK1" --bfile "$TMP" --update-ids "$RENAME" --make-bed \
  --threads "$THREADS" --memory "$MEMMB" --out "$OUT" > "$OUT.rename.log" 2>&1 \
  || { echo "rename failed; see $OUT.rename.log"; tail -8 "$OUT.rename.log"; exit 6; }

[[ -f "$OUT.bed" ]] || { echo "ERROR: final bed not produced" >&2; exit 7; }
log "DONE: $(wc -l < "$OUT.bim") variants x $(wc -l < "$OUT.fam") samples -> $OUT.{bed,bim,fam}"
# cleanup intermediates (keep per-chr beds for now; remove if you want space)
rm -f "$TMP".{bed,bim,fam} "$LIST"
log "sample IIDs:"; awk 'NR<=3{print $2}' "$OUT.fam"   # no | head (avoids SIGPIPE under pipefail)
