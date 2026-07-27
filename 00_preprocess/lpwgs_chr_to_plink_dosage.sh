#!/usr/bin/env bash
# ============================================================
# lpwgs_chr_to_plink_dosage.sh   — DOSAGE variant of lpwgs_chr_to_plink.sh
# Convert ONE chromosome of a lpWGS VCF -> plink2 PGEN, KEEPING the imputed
# dosages (DS) instead of hard-calling. Low-pass WGS's power lives in the
# dosages, so downstream --glm/--score read them directly (plink2 pgen stores
# dosages natively). This runs ALONGSIDE the hard-call converter; it does not
# touch lpwgs_clean.* — it writes lpwgs_dosage_* instead.
#
# Usage (SLURM array task, CHR from $SLURM_ARRAY_TASK_ID):
#   lpwgs_chr_to_plink_dosage.sh --vcf F.vcf.gz --keep KEEP --out-prefix PFX --chr 7
#       [--maf 0] [--extract LIST] [--bcftools bcftools] [--plink2 plink2]
#       [--threads 4] [--mem-mb 28000] [--chr-name-prefix chr]
# Produces: <out-prefix>_chr<chr>.{pgen,pvar,psam}   (dosages retained)
# ============================================================
set -euo pipefail
VCF=""; KEEP=""; OUTPFX=""; CHR=""; MAF=0; EXTRACT=""; NOPASS=0
BCFTOOLS=bcftools; PLINK2=plink2; THREADS=4; MEMMB=28000; CHRNP=chr
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2;;
    --keep) KEEP="$2"; shift 2;;
    --out-prefix) OUTPFX="$2"; shift 2;;
    --chr) CHR="$2"; shift 2;;
    --maf) MAF="$2"; shift 2;;
    --no-pass) NOPASS=1; shift;;          # keep ALL variants (incl. LOWCONF); dosage analysis uses --maf instead
    --extract) EXTRACT="$2"; shift 2;;   # restrict to union-KEEP variant list (chr:pos)
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
FOPT=(); [[ "$NOPASS" -eq 0 ]] && FOPT=(-f PASS)   # PASS is genotype-certainty-biased; --no-pass keeps LOWCONF too
log "step 1: bcftools view ${FOPT[*]:-(all variants)} -r $REGION -> $TMPB"
"$BCFTOOLS" view "${FOPT[@]}" -r "$REGION" --threads "$THREADS" -Ob -o "$TMPB" "$VCF" 2> "$OUT.bcftools.log" \
  || { echo "chr$CHR bcftools failed; tail:"; tail -6 "$OUT.bcftools.log"; rm -f "$TMPB"; exit 4; }

# KEY DIFFERENCE vs hard-call converter: --make-pgen (NOT --make-bed) and NO
# --hard-call-threshold, so plink2 keeps the DS dosages in the pgen.
log "step 2: plink2 --bcf $TMPB dosage=DS  --make-pgen (dosages kept${EXTRACT:+ ; --extract})"
EXOPT=(); [[ -n "$EXTRACT" ]] && EXOPT=(--extract "$EXTRACT")
MAFOPT=(); [[ -n "$MAF" && "$MAF" != "0" ]] && MAFOPT=(--maf "$MAF")
"$PLINK2" --bcf "$TMPB" dosage=DS \
    --double-id --keep "$KEEP" \
    --snps-only just-acgt --max-alleles 2 \
    --set-all-var-ids '@:#' --rm-dup force-first \
    "${MAFOPT[@]}" "${EXOPT[@]}" --output-chr chrM \
    --make-pgen --threads "$THREADS" --memory "$MEMMB" \
    --out "$OUT" > "$OUT.log" 2>&1 \
  || { echo "chr$CHR plink2 failed; tail of log:"; tail -8 "$OUT.log"; rm -f "$TMPB"; exit 5; }
rm -f "$TMPB"   # free the temp BCF once the pgen is written

if [[ -f "$OUT.pvar" ]]; then
  log "DONE: $(grep -cv '^#' "$OUT.pvar") variants x $(grep -cv '^#' "$OUT.psam") samples -> $OUT.{pgen,pvar,psam}"
else
  echo "ERROR: no pgen for chr$CHR (maybe 0 PASS variants?); see $OUT.log" >&2; exit 6
fi
