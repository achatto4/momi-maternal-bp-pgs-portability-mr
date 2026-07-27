#!/usr/bin/env bash
# ============================================================
# build_common_union.sh
# Emit chr:pos of variants that are COMMON (MAF>=maf) and quality-pass in ONE
# cohort x platform, read straight from the VCF allele frequencies (no genotypes).
#   GSA  : INFO/MAF >= maf  AND  INFO/R2 >= r2
#   lpWGS: FILTER == PASS   AND  maf <= AF <= 1-maf
# Run once per cohort x platform; the union of all outputs = "common in >=1 ancestry".
#
# Usage:
#   build_common_union.sh --vcf F.vcf.gz --platform gsa|lpwgs --out FILE
#       [--maf 0.01] [--r2 0.8] [--bcftools bcftools] [--threads 4]
# Output: FILE with one "chrN:pos" per line (matches our '@:#' bim IDs), sorted-unique.
# ============================================================
set -euo pipefail
VCF=""; PLAT=""; OUT=""; MAF=0.01; R2=0.8; BCFTOOLS=bcftools; THREADS=4
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2;;
    --platform) PLAT="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --maf) MAF="$2"; shift 2;;
    --r2) R2="$2"; shift 2;;
    --bcftools) BCFTOOLS="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
for v in VCF PLAT OUT; do [[ -n "${!v}" ]] || { echo "missing --$v"; exit 1; }; done
[[ -f "$VCF" ]] || { echo "VCF not found: $VCF" >&2; exit 2; }
command -v "$BCFTOOLS" >/dev/null || { echo "bcftools not on PATH" >&2; exit 3; }
mkdir -p "$(dirname "$OUT")"
HI=$(awk "BEGIN{print 1-$MAF}")
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

case "$PLAT" in
  gsa)   EXPR="INFO/MAF>=$MAF && INFO/R2>=$R2";;
  lpwgs) EXPR="FILTER=\"PASS\" && AF>=$MAF && AF<=$HI";;
  *) echo "platform must be gsa|lpwgs" >&2; exit 1;;
esac
log "querying common variants ($PLAT): $EXPR"
# NOTE: bcftools 'query' does not accept --threads (only 'view' does)
"$BCFTOOLS" query -f '%CHROM:%POS\n' -i "$EXPR" "$VCF" \
  | sort -u > "$OUT"
log "DONE: $(wc -l < "$OUT") common variants -> $OUT"
