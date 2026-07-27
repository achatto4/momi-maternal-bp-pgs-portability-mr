#!/usr/bin/env bash
# ============================================================
# build_keep_sets.sh
# From the per-cohort x platform common-variant lists (<COHORT>__<PLATFORM>.common
# produced by build_common_union.sh), construct two candidate KEEP sets:
#
#   KEEP_anyunion.txt   = common (MAF>=0.01) in >=1 cohort  (liberal)
#   KEEP_concordant.txt = common on BOTH GSA AND lpWGS within >=1 ancestry
#                         (recommended for MEGA: drops platform-discordant artifacts
#                          while keeping ancestry-specific common variants)
#
# Ancestry groups (edit if cohorts change):
#   SAS = AMANHI-Bangladesh, AMANHI-Pakistan, GAPPS-Bangladesh
#   AFR = AMANHI-Pemba, GAPPS-Zambia
#
# Usage: build_keep_sets.sh --keepdir DIR
# ============================================================
set -euo pipefail
KEEPDIR=""
while [[ $# -gt 0 ]]; do case "$1" in --keepdir) KEEPDIR="$2"; shift 2;; *) echo "Unknown arg: $1" >&2; exit 1;; esac; done
[[ -n "$KEEPDIR" ]] || { echo "need --keepdir"; exit 1; }
cd "$KEEPDIR"
SAS=(AMANHI-Bangladesh AMANHI-Pakistan GAPPS-Bangladesh)
AFR=(AMANHI-Pemba GAPPS-Zambia)

# union of a set of .common files for one platform within an ancestry
plat_union(){ local plat="$1"; shift; local fs=(); for c in "$@"; do [[ -f "${c}__${plat}.common" ]] && fs+=("${c}__${plat}.common"); done
  [[ ${#fs[@]} -gt 0 ]] && sort -u "${fs[@]}" || true; }

echo "building per-ancestry per-platform common sets..."
plat_union gsa   "${SAS[@]}" > .sas_gsa.tmp
plat_union lpwgs "${SAS[@]}" > .sas_lp.tmp
plat_union gsa   "${AFR[@]}" > .afr_gsa.tmp
plat_union lpwgs "${AFR[@]}" > .afr_lp.tmp

# concordant within ancestry = common on BOTH platforms
comm -12 .sas_gsa.tmp .sas_lp.tmp > .sas_concord.tmp
comm -12 .afr_gsa.tmp .afr_lp.tmp > .afr_concord.tmp
sort -u .sas_concord.tmp .afr_concord.tmp > KEEP_concordant.txt

# liberal = common in any cohort/platform
sort -u ./*.common > KEEP_anyunion.txt

echo "=== KEEP set sizes ==="
printf "  SAS GSA-common      : %s\n" "$(wc -l < .sas_gsa.tmp)"
printf "  SAS lpWGS-common    : %s\n" "$(wc -l < .sas_lp.tmp)"
printf "  SAS concordant      : %s\n" "$(wc -l < .sas_concord.tmp)"
printf "  AFR GSA-common      : %s\n" "$(wc -l < .afr_gsa.tmp)"
printf "  AFR lpWGS-common    : %s\n" "$(wc -l < .afr_lp.tmp)"
printf "  AFR concordant      : %s\n" "$(wc -l < .afr_concord.tmp)"
printf "  KEEP_concordant.txt : %s  (recommended for mega)\n" "$(wc -l < KEEP_concordant.txt)"
printf "  KEEP_anyunion.txt   : %s  (liberal)\n" "$(wc -l < KEEP_anyunion.txt)"
rm -f .sas_gsa.tmp .sas_lp.tmp .afr_gsa.tmp .afr_lp.tmp .sas_concord.tmp .afr_concord.tmp
