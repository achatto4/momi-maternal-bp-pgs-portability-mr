#!/usr/bin/env bash
# ============================================================
# vcf_to_plink.sh
# Convert ONE cohort x platform imputed VCF -> QC'd PLINK1 bed, mothers only,
# renamed to Epi PARTICIPANT_ID, variant IDs as chr:pos (autosomes), GRCh38.
#
# Two-pass design (avoids plink2 --keep/--update-ids ordering ambiguity):
#   PASS 1 (heavy): read VCF, QC, keep mothers by VCF id (--double-id), set chr:pos IDs
#   PASS 2 (fast) : rename samples to Epi PARTICIPANT_ID on the small subset
#
# Platform differences:
#   gsa   : dosage=HDS ; imputation quality = INFO/R2  (keep R2 >= r2_min)
#   lpwgs : GT:DS:GP    ; no R2 -> keep FILTER==PASS (bcftools pre-filter)
#
# Inputs: the .keep and .rename produced by make_sample_map.sh
# Output: <out>.{bed,bim,fam} + logs
# Designed for SLURM (one job per cohort x platform). Validate on GSA first.
# ============================================================
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $0 --vcf FILE --platform gsa|lpwgs --keep FILE --rename FILE --out PREFIX
          [--r2-min 0.3] [--maf 0.01] [--hard-call 0.1]
          [--plink2 plink2] [--bcftools bcftools] [--threads 8] [--mem-mb 60000]
USAGE
  exit 1
}

VCF=""; PLAT=""; KEEP=""; RENAME=""; OUT=""; EXTRACT=""
R2MIN=0.8; MAF=0; HCT=0.4999   # with --extract (union KEEP), MAF off; quality via R2
PLINK2=plink2; BCFTOOLS=bcftools; THREADS=8; MEMMB=60000
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2;;
    --platform) PLAT="$2"; shift 2;;
    --keep) KEEP="$2"; shift 2;;
    --rename) RENAME="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --extract) EXTRACT="$2"; shift 2;;   # restrict to a union-KEEP variant list (chr:pos)
    --r2-min) R2MIN="$2"; shift 2;;
    --maf) MAF="$2"; shift 2;;
    --hard-call) HCT="$2"; shift 2;;
    --plink2) PLINK2="$2"; shift 2;;
    --bcftools) BCFTOOLS="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mem-mb) MEMMB="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done
[[ -n "$VCF" && -n "$PLAT" && -n "$KEEP" && -n "$RENAME" && -n "$OUT" ]] || usage
[[ -f "$VCF" ]] || { echo "VCF not found: $VCF" >&2; exit 2; }
command -v "$PLINK2" >/dev/null || { echo "plink2 not on PATH (load module?)" >&2; exit 3; }
mkdir -p "$(dirname "$OUT")"

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
log "platform=$PLAT  vcf=$VCF"
log "keep=$(wc -l < "$KEEP") mothers   r2_min=$R2MIN maf=$MAF hardcall=$HCT"

TMP="${OUT}.tmp"

# ---- shared QC / formatting options (applied in PASS 1) ----
#   --double-id            FID=IID=VCF sample name (matches make_sample_map .keep)
#   --keep                 mothers only (BE), by VCF id
#   --chr 1-22             autosomes
#   --snps-only just-acgt  drop indels / non-ACGT for clean positional merge
#   --max-alleles 2        biallelic only
#   --set-all-var-ids @:#  chr:pos IDs (positional merge key, mirrors prior pipeline)
#   --rm-dup force-first   one record per chr:pos (keep first)
#   --maf                  minor-allele-freq filter (post-subset)
#   --output-chr chrM      keep 'chr' prefix consistently (GRCh38)
QC=( --double-id
     --keep "$KEEP"
     --chr 1-22
     --snps-only just-acgt
     --max-alleles 2
     --set-all-var-ids '@:#'
     --rm-dup force-first
     --maf "$MAF"
     --output-chr chrM
     --make-bed
     --threads "$THREADS" --memory "$MEMMB"
     --out "$TMP" )
# restrict to the union-KEEP variant list, if given (applied after --set-all-var-ids)
[[ -n "$EXTRACT" ]] && QC+=( --extract "$EXTRACT" )

log "PASS 1: VCF -> QC'd bed (mothers, by VCF id)${EXTRACT:+ ; --extract $EXTRACT}"
case "$PLAT" in
  gsa)
    "$PLINK2" --vcf "$VCF" dosage=HDS \
      --extract-if-info "R2 >= $R2MIN" \
      --hard-call-threshold "$HCT" \
      "${QC[@]}" 2>&1 | tee "$OUT.pass1.log"
    ;;
  lpwgs)
    command -v "$BCFTOOLS" >/dev/null || { echo "bcftools not on PATH" >&2; exit 3; }
    log "(lpWGS) bcftools view -f PASS | plink2 --bcf  [validate throughput on GSA first]"
    "$BCFTOOLS" view -f PASS -Ou "$VCF" \
      | "$PLINK2" --bcf /dev/stdin dosage=DS \
          --hard-call-threshold "$HCT" \
          "${QC[@]}" 2>&1 | tee "$OUT.pass1.log"
    ;;
  *) echo "platform must be gsa|lpwgs" >&2; exit 1;;
esac

[[ -f "$TMP.bed" ]] || { log "ERROR: PASS 1 produced no bed; see $OUT.pass1.log"; exit 4; }
log "PASS 1 done: $(wc -l < "$TMP.bim") variants x $(wc -l < "$TMP.fam") samples"

log "PASS 2: rename samples -> Epi PARTICIPANT_ID"
"$PLINK2" --bfile "$TMP" \
  --update-ids "$RENAME" \
  --make-bed --threads "$THREADS" --memory "$MEMMB" \
  --out "$OUT" 2>&1 | tee "$OUT.pass2.log"

[[ -f "$OUT.bed" ]] || { log "ERROR: PASS 2 failed; see $OUT.pass2.log"; exit 5; }
log "cleaning tmp"; rm -f "$TMP".{bed,bim,fam,log}

log "DONE: $(wc -l < "$OUT.bim") variants x $(wc -l < "$OUT.fam") samples -> $OUT.{bed,bim,fam}"
log "sample of final IIDs:"; awk 'NR<=3{print $2}' "$OUT.fam"   # no | head (SIGPIPE under pipefail)
