#!/usr/bin/env bash
# ============================================================
# lpwgs_chr_to_plink.sh
# Convert ONE chromosome of a lpWGS VCF -> plink bed (mothers only, QC, chr:pos IDs).
# Uses the VCF's tabix index to read just that chromosome (fast, parallelizable).
# No sample rename here (done once after the 22 chromosomes are merged).
#
# Usage (typically as a SLURM array task, CHR from $SLURM_ARRAY_TASK_ID):
#   lpwgs_chr_to_plink.sh --vcf F.vcf.gz --keep KEEP --out-prefix PFX --chr 7
#       [--maf 0.01] [--hard-call 0.4999] [--bcftools bcftools] [--plink2 plink2]
#       [--threads 4] [--mem-mb 28000] [--chr-name-prefix chr]
# Produces: <out-prefix>_chr<chr>.{bed,bim,fam}
# ============================================================
set -euo pipefail
VCF=""; KEEP=""; OUTPFX=""; CHR=""; MAF=0; HCT=0.4999; INFOMIN=0.8; EXTRACT=""
BCFTOOLS=bcftools; PLINK2=plink2; THREADS=4; MEMMB=28000; CHRNP=chr
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2;;
    --keep) KEEP="$2"; shift 2;;
    --out-prefix) OUTPFX="$2"; shift 2;;
    --chr) CHR="$2"; shift 2;;
    --maf) MAF="$2"; shift 2;;
    --info-min) INFOMIN="$2"; shift 2;;
    --extract) EXTRACT="$2"; shift 2;;   # restrict to union-KEEP variant list (chr:pos)
    --hard-call) HCT="$2"; shift 2;;
    --bcftools) BCFTOOLS="$2"; shift 2;;
    --plink2) PLINK2="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mem-mb) MEMMB="$2"; shift 2;;
    --chr-name-prefix) CHRNP="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -z "$CHR" && -n "${SLURM_ARRAY_TASK_ID:-}" ]] && CHR="$SLURM_ARRAY_TASK_ID"
for v in VCF KEEP OUTPFX CHR; do [[ -n "${!v}" ]] || { echo "missing --$v"; exit 1; }; done
[[ -f "$VCF" ]] || { echo "VCF not found: $VCF" >&2; exit 2; }
command -v "$BCFTOOLS" >/dev/null || { echo "bcftools not on PATH" >&2; exit 3; }
command -v "$PLINK2"   >/dev/null || { echo "plink2 not on PATH" >&2; exit 3; }
OUT="${OUTPFX}_chr${CHR}"
log(){ printf '[%s chr%s] %s\n' "$(date '+%F %T')" "$CHR" "$*"; }

REGION="${CHRNP}${CHR}"
TMPB="${OUT}.pass.bcf"   # real file (NOT a pipe): plink2 cannot finish a streamed BCF import
# Quality = FILTER==PASS (Gencove max(GP)>0.9). The header declares an INFO/INFO score
# but it is NOT populated in the records, so we cannot filter on it. MAF is left to the
# pooled merge (--maf "$MAF" with MAF=0 here = no per-cohort MAF).
log "step 1: bcftools view -f PASS -r $REGION -> $TMPB"
"$BCFTOOLS" view -f PASS -r "$REGION" --threads "$THREADS" -Ob -o "$TMPB" "$VCF" 2> "$OUT.bcftools.log" \
  || { echo "chr$CHR bcftools failed; tail:"; tail -6 "$OUT.bcftools.log"; rm -f "$TMPB"; exit 4; }
log "step 2: plink2 --bcf $TMPB dosage=DS  (mothers only${EXTRACT:+ ; --extract})"
EXOPT=(); [[ -n "$EXTRACT" ]] && EXOPT=(--extract "$EXTRACT")
MAFOPT=(); [[ -n "$MAF" && "$MAF" != "0" ]] && MAFOPT=(--maf "$MAF")
"$PLINK2" --bcf "$TMPB" dosage=DS \
    --double-id --keep "$KEEP" \
    --hard-call-threshold "$HCT" \
    --snps-only just-acgt --max-alleles 2 \
    --set-all-var-ids '@:#' --rm-dup force-first \
    "${MAFOPT[@]}" "${EXOPT[@]}" --output-chr chrM \
    --make-bed --threads "$THREADS" --memory "$MEMMB" \
    --out "$OUT" > "$OUT.log" 2>&1 \
  || { echo "chr$CHR plink2 failed; tail of log:"; tail -8 "$OUT.log"; rm -f "$TMPB"; exit 5; }
rm -f "$TMPB"   # free the temp BCF once the bed is written

if [[ -f "$OUT.bim" ]]; then
  log "DONE: $(wc -l < "$OUT.bim") variants x $(wc -l < "$OUT.fam") samples -> $OUT.{bed,bim,fam}"
else
  echo "ERROR: no bed for chr$CHR (maybe 0 PASS variants?); see $OUT.log" >&2; exit 6
fi
