#!/usr/bin/env bash
# ============================================================
# prep_instruments.sh  (Stage 04a)
# Select MR instruments from an EXTERNAL BP GWAS (exposure): genome-wide
# significant SNPs, LD-clumped to independence using the merged genotypes.
# Variant IDs are rewritten to chr:pos to match our bim and the internal GWAS.
#
# Output: <out-dir>/instruments_<trait>.txt  with header:
#   SNP EA OA BETA SE EAF P N         (SNP = chr<chr>:<pos>, GRCh38)
#
# Usage:
#   prep_instruments.sh --ext FILE.tsv.gz --trait SBP --merged PREFIX --out-dir DIR
#       [--p1 5e-8] [--clump-r2 0.001] [--clump-kb 10000]
#       [--chr-prefix chr] [--plink2 plink2] [--threads 4] [--mem-mb 24000]
# Expects external GWAS-Catalog harmonized columns:
#   chromosome base_pair_location effect_allele other_allele beta
#   standard_error effect_allele_frequency p_value rsid n
# ============================================================
set -euo pipefail
EXT=""; TRAIT=""; MERGED=""; OUTDIR=""; P1=5e-8; R2=0.001; KB=10000
CHRPFX=chr; PLINK1=plink; THREADS=4; MEMMB=24000   # clumping uses plink1.9 (plink2 has no --clump)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ext) EXT="$2"; shift 2;;
    --trait) TRAIT="$2"; shift 2;;
    --merged) MERGED="$2"; shift 2;;
    --out-dir) OUTDIR="$2"; shift 2;;
    --p1) P1="$2"; shift 2;;
    --clump-r2) R2="$2"; shift 2;;
    --clump-kb) KB="$2"; shift 2;;
    --chr-prefix) CHRPFX="$2"; shift 2;;
    --plink1) PLINK1="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mem-mb) MEMMB="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
for v in EXT TRAIT MERGED OUTDIR; do [[ -n "${!v}" ]] || { echo "missing --$v"; exit 1; }; done
[[ -f "$EXT" ]] || { echo "external GWAS not found: $EXT" >&2; exit 2; }
command -v "$PLINK1" >/dev/null || { echo "plink (1.9) not on PATH; module load plink/1.90b" >&2; exit 3; }
mkdir -p "$OUTDIR"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

SIG="$OUTDIR/exposure_${TRAIT}_sig.txt"      # genome-wide-sig exposure rows (ID-keyed)
CLIN="$OUTDIR/clump_in_${TRAIT}.txt"         # plink2 --clump input (ID P)
log "extracting genome-wide-significant exposure SNPs (p<=$P1) from $EXT"
zcat "$EXT" | awk -v OFS='\t' -v p1="$P1" -v chrpfx="$CHRPFX" -v sig="$SIG" -v clin="$CLIN" '
  BEGIN{ FS="\t" }
  NR==1{
    for(i=1;i<=NF;i++) c[$i]=i
    chr=c["chromosome"]; pos=c["base_pair_location"]; ea=c["effect_allele"]; oa=c["other_allele"]
    b=c["beta"]; se=c["standard_error"]; eaf=c["effect_allele_frequency"]; p=c["p_value"]; n=c["n"]
    if(!chr||!pos||!ea||!oa||!b||!se||!p){ print "ERROR: missing required external columns" > "/dev/stderr"; exit 3 }
    print "SNP","EA","OA","BETA","SE","EAF","P","N" > sig
    print "ID","P" > clin
    next
  }
  {
    if($p=="" || $p=="NA") next
    if($p+0 > p1+0) next
    id = chrpfx $chr ":" $pos
    print id, $ea, $oa, $b, $se, (eaf?$eaf:"NA"), $p, (n?$n:"NA") >> sig
    print id, $p >> clin
  }'
NSIG=$(( $(wc -l < "$SIG") - 1 ))
log "genome-wide-sig SNPs: $NSIG"
[[ "$NSIG" -gt 0 ]] || { echo "ERROR: no SNPs passed p<=$P1" >&2; exit 4; }

log "LD-clumping against $MERGED (plink1.9; r2<$R2, ${KB}kb)"
"$PLINK1" --bfile "$MERGED" \
  --clump "$CLIN" --clump-p1 "$P1" --clump-r2 "$R2" --clump-kb "$KB" \
  --clump-snp-field ID --clump-field P \
  --threads "$THREADS" --memory "$MEMMB" \
  --out "$OUTDIR/clump_${TRAIT}" > "$OUTDIR/clump_${TRAIT}.run.log" 2>&1 \
  || { echo "plink --clump failed; see $OUTDIR/clump_${TRAIT}.run.log"; tail -5 "$OUTDIR/clump_${TRAIT}.run.log"; exit 5; }

CLUMPS="$OUTDIR/clump_${TRAIT}.clumped"   # plink1.9 writes .clumped
[[ -f "$CLUMPS" ]] || { echo "no .clumped produced (0 clumps?); see log" >&2; exit 6; }
# lead-SNP IDs = the SNP column of .clumped (skip header + trailing blank lines)
awk 'NR==1{for(i=1;i<=NF;i++) if($i=="SNP") idc=i; next} NF>0 && $idc!=""{print $idc}' "$CLUMPS" | sort -u > "$OUTDIR/lead_${TRAIT}.ids"
NLEAD=$(wc -l < "$OUTDIR/lead_${TRAIT}.ids")
log "independent instruments: $NLEAD"

# join lead IDs back to exposure stats
awk -v OFS='\t' 'NR==FNR{keep[$1]=1; next} FNR==1{print; next} ($1 in keep)' \
  "$OUTDIR/lead_${TRAIT}.ids" "$SIG" > "$OUTDIR/instruments_${TRAIT}.txt"
log "DONE: $OUTDIR/instruments_${TRAIT}.txt ($(( $(wc -l < "$OUTDIR/instruments_${TRAIT}.txt") - 1 )) instruments)"
