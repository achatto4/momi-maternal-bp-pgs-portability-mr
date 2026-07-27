#!/usr/bin/env bash
# ============================================================
# fetch_pgs.sh — download a harmonised PGS Catalog scoring file into $MOMI_SCORES.
#
# WHICH FILE AND WHY. We fetch the *_hmPOS_GRCh38.txt.gz variant, not the raw scoring file.
# compute_prs.sh builds its variant IDs from hm_chr/hm_pos, which exist only in the
# harmonised files, and our genotypes are GRCh38. Downloading the raw file would parse to
# empty or, worse, to GRCh37 coordinates that silently fail to match the bim -- producing a
# score with near-zero overlap rather than an error.
#
#   https://www.pgscatalog.org/downloads/
#   ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/<ID>/ScoringFiles/Harmonized/<ID>_hmPOS_GRCh38.txt.gz
#
# VERIFY BEFORE USE. The script checks the md5 published alongside the file, and then
# confirms the required columns are present, because a truncated or redirected download can
# still be a valid gzip. A score file that parses to the wrong build is the failure mode that
# would be hardest to notice downstream -- it looks like poor transferability, not an error.
#
#   bash bin/fetch_pgs.sh PGS000027 [PGS000123 ...]
# ============================================================
set -uo pipefail
: "${MOMI_SCORES:?export MOMI_SCORES first}"
BASE="https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores"
[[ $# -gt 0 ]] || { echo "usage: fetch_pgs.sh PGSID [PGSID ...]" >&2; exit 1; }
mkdir -p "$MOMI_SCORES"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

rc_all=0
for ID in "$@"; do
  OUT="$MOMI_SCORES/${ID}.txt.gz"
  if [[ -s "$OUT" ]]; then log "$ID: already present at $OUT — skipping"; continue; fi
  URL="$BASE/${ID}/ScoringFiles/Harmonized/${ID}_hmPOS_GRCh38.txt.gz"
  TMP="$OUT.part"
  log "$ID: fetching $URL"
  if ! curl -fsSL --retry 3 --retry-delay 5 -o "$TMP" "$URL"; then
    log "$ID: DOWNLOAD FAILED. Check the ID exists and that this node has outbound https."
    rm -f "$TMP"; rc_all=1; continue
  fi

  ## md5, if the catalog publishes one next to the file
  if curl -fsSL --retry 2 -o "$TMP.md5" "$URL.md5" 2>/dev/null; then
    want=$(awk '{print $1}' "$TMP.md5")
    got=$(md5sum "$TMP" | awk '{print $1}')
    if [[ -n "$want" && "$want" != "$got" ]]; then
      log "$ID: MD5 MISMATCH (want $want got $got) — refusing to install"
      rm -f "$TMP" "$TMP.md5"; rc_all=1; continue
    fi
    log "$ID: md5 ok"
    rm -f "$TMP.md5"
  else
    log "$ID: no published md5 — continuing without checksum verification"
  fi

  ## structural check: gunzip cleanly, and carry the columns compute_prs.sh needs
  if ! gzip -t "$TMP" 2>/dev/null; then
    log "$ID: not a valid gzip — refusing to install"; rm -f "$TMP"; rc_all=1; continue
  fi
  hdr=$(zcat "$TMP" | grep -v '^#' | head -1)
  miss=""
  for c in hm_chr hm_pos effect_allele effect_weight; do
    grep -qw -- "$c" <<< "$hdr" || miss="$miss $c"
  done
  if [[ -n "$miss" ]]; then
    log "$ID: MISSING REQUIRED COLUMNS:$miss"
    log "     header was: $hdr"
    log "     compute_prs.sh needs hm_chr/hm_pos (GRCh38). Refusing to install."
    rm -f "$TMP"; rc_all=1; continue
  fi

  mv "$TMP" "$OUT"
  n=$(zcat "$OUT" | grep -vc '^#')
  log "$ID: installed $OUT ($((n-1)) variant rows)"
done

echo
echo "next: add the id to 05_prs/transfer_panel.tsv if it is not there, then"
echo "      bash \$MOMI_PIPE/bin/submit_score_array.sh"
exit $rc_all
