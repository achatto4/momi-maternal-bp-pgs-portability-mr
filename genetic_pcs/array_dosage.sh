#!/usr/bin/env bash
# array_dosage.sh
#
# Array dosages at the pruned variants from the re-imputed array VCF.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/config_pcs.sh"
COH="$1"; IDS="$2"; KP="$3"; OP="$4"; T="${OP}_tmp"; mkdir -p "$T"
P2=( "$JPCA_PLINK2" --threads "$JPCA_THREADS" --memory "$JPCA_MEMMB" )
VCF="$JPCA_ARRAY_VCF_DIR/$COH/gsa_reimputed.vcf.gz"; MAP="$JPCA_FILESETS/$COH/gsa_map.rename"; LP="$JPCA_FILESETS/$COH/lpwgs_dosage_clean"
[[ -f "$VCF" && -f "$MAP" ]] || { echo "dosage mode needs $VCF and $MAP" >&2; exit 3; }
HAVE_BCF=0; command -v "$JPCA_BCFTOOLS" >/dev/null 2>&1 && HAVE_BCF=1
hdr_vcf(){ if [[ $HAVE_BCF -eq 1 ]]; then "$JPCA_BCFTOOLS" view -h "$1"; else { gzip -dc "$1" 2>/dev/null || true; } | awk '/^#CHROM/{print; exit} /^##/{print}'; fi; }   # awk stops at #CHROM; gzip's SIGPIPE is expected
hdr_vcf "$VCF" > "$T/vcf.header"
FIELDS=$(grep '^##FORMAT' "$T/vcf.header" | sed 's/.*ID=\([^,]*\),.*/\1/' | paste -sd, -)
if [[ ",$FIELDS," == *",HDS,"* ]]; then DFIELD=HDS; elif [[ ",$FIELDS," == *",DS,"* ]]; then DFIELD=DS; else echo "no HDS/DS in $VCF" >&2; exit 3; fi
"${P2[@]}" --pfile "$LP" --extract "$IDS" --make-just-pvar --out "$T/lp_sub" > "$T/lp_sub.log" 2>&1
PFX=$(grep -m1 '^##contig=<ID=' "$T/vcf.header" | sed 's/^##contig=<ID=\([^,>]*\).*/\1/' | grep -q '^chr' && echo chr || echo "")
grep -v '^#' "$T/lp_sub.pvar" | awk -v p="$PFX" 'BEGIN{OFS="\t"}{c=$1; sub(/^chr/,"",c); print p c, $2}' | sort -k1,1V -k2,2n -u > "$T/regions.tsv"
awk 'NR==FNR{if($1 !~ /^#/) want[($2!="")?$2:$1]=1; next} ($4 in want){print $2}' "$KP" "$MAP" | sort -u > "$T/vcf_samples.txt"
NS=$(wc -l < "$T/vcf_samples.txt"); NW=$(grep -vc '^#' "$KP" || true)
[[ "$NS" -eq "$NW" ]] || { echo "only $NS of $NW array records map to a VCF sample" >&2; exit 4; }
if [[ $HAVE_BCF -eq 1 ]]; then                          # indexed region and sample subset, then plink2 on the small BCF
  if [[ -f "$VCF.tbi" || -f "$VCF.csi" ]]; then RSEL=(-R "$T/regions.tsv"); else RSEL=(-T "$T/regions.tsv"); fi
  "$JPCA_BCFTOOLS" view "${RSEL[@]}" -S "$T/vcf_samples.txt" --force-samples -Ob -o "$T/arr.bcf" "$VCF" 2> "$T/bcftools.log"
  "${P2[@]}" --bcf "$T/arr.bcf" dosage="$DFIELD" --double-id --snps-only just-acgt --max-alleles 2 --set-all-var-ids '@:#' \
    --rm-dup force-first --make-pgen --out "$T/arr_raw" > "$T/arr_raw.log" 2>&1
else                                                    # no bcftools: plink2 reads the VCF itself, keeping the same positions and samples
  awk 'BEGIN{OFS="\t"}{print $1, $2, $2, "r" NR}' "$T/regions.tsv" > "$T/regions.range"
  awk '{print $1"\t"$1}' "$T/vcf_samples.txt" > "$T/vcf_samples.keep"
  "${P2[@]}" --vcf "$VCF" dosage="$DFIELD" --double-id --keep "$T/vcf_samples.keep" --extract range "$T/regions.range" --snps-only just-acgt \
    --max-alleles 2 --set-all-var-ids '@:#' --rm-dup force-first --make-pgen --out "$T/arr_raw" > "$T/arr_raw.log" 2>&1
fi
"${P2[@]}" --pfile "$T/arr_raw" --update-ids "$MAP" --make-pgen --out "$T/arr_ren" > "$T/arr_ren.log" 2>&1
Rscript "$HERE/array_dosage_align.R" "$T/arr_ren.pvar" "$T/lp_sub.pvar" "$T" > "$T/align.log" 2>&1
"${P2[@]}" --pfile "$T/arr_ren" --update-name "$T/arr_rename.txt" 2 1 --make-pgen --out "$T/arr_ren2" > "$T/arr_ren2.log" 2>&1
"${P2[@]}" --pfile "$T/arr_ren2" --extract "$T/arr_extract_lp.ids" --ref-allele force "$T/arr_lp_ref.txt" 2 1 \
  --make-pgen --out "$OP" > "$OP.log" 2>&1
cp "$T/agg_array_dosage_import.tsv" "${OP}_import.tsv"; echo "$DFIELD $([[ $HAVE_BCF -eq 1 ]] && echo bcftools || echo plink2-vcf)" > "${OP}_field.txt"
rm -rf "$T"
