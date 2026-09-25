#!/usr/bin/env bash
# cohort_genotypes.sh
#
# Per cohort: variant quality control on each technology, the shared allele-aligned variant set, the cross-technology
# frequency check, LD pruning, the array dosages or hard calls, the blockwise relationship matrix and the PCA.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/config_pcs.sh"
COH="$1"; ENC="$2"
[[ "$ENC" == dosage || "$ENC" == hardcall ]] || { echo "encoding must be dosage or hardcall" >&2; exit 2; }
W="$JPCA_OUT/work/$COH"; AGG="$JPCA_OUT/aggregate/cohort/$COH"; mkdir -p "$W/blocks" "$AGG"
echo "started" > "$W/status.txt"; STEP="setup"; trap 'echo "failed at: $STEP" > "$W/status.txt"' ERR
P2=( "$JPCA_PLINK2" --threads "$JPCA_THREADS" --memory "$JPCA_MEMMB" )
say(){ printf '[%s] %s: %s\n' "$(date '+%F %T')" "$COH" "$*" | tee -a "$W/run_log.txt"; }
LP="$JPCA_FILESETS/$COH/lpwgs_dosage_clean"; GS="$JPCA_FILESETS/$COH/gsa_clean"
if   [[ -f "$GS.pgen" ]]; then GSIN=(--pfile "$GS"); GSSAMP="$GS.psam"
elif [[ -f "$GS.bed"  ]]; then GSIN=(--bfile "$GS"); GSSAMP="$GS.fam"
else echo "no gsa_clean fileset for $COH" >&2; exit 3; fi
HL="$HERE/highld_regions_hg38.txt"
QC=(--autosome --snps-only just-acgt --max-alleles 2 --exclude range "$HL" --geno 0.05 --maf 0.05)
nrow(){ local n; n=$(grep -vc '^#' "$1" 2>/dev/null || true); echo "${n:-0}"; }

STEP="1 record selection and reference"; say "$STEP"
Rscript "$HERE/select_records.R" --cohort "$COH" --keep "$JPCA_KEEP" --drops "$JPCA_DROPS" --lp-psam "$LP.psam" \
  --gsa-samples "$GSSAMP" --sscore-dir "$JPCA_SSC" --pairs "$JPCA_RELATEDNESS_DIR/$COH/pairs_person_level.tsv" \
  --dual-self "$JPCA_RELATEDNESS_DIR/$COH/xking_dual_self.tsv" --out "$W" --seed "$JPCA_SEED" 2>&1 | tee -a "$W/run_log.txt"
NGS=$(nrow "$W/gsa_all.keep")
[[ "$NGS" -gt 0 ]] || { echo "no array records in $COH: the joint-platform design needs both technologies" >&2; exit 3; }

STEP="2 variant QC on each technology"; say "$STEP"
"${P2[@]}" --pfile "$LP" --keep "$W/sel_lp.keep" "${QC[@]}" --mach-r2-filter 0.8 --make-just-pvar --out "$W/lp_qc" > "$W/lp_qc.log" 2>&1
"${P2[@]}" "${GSIN[@]}" --keep "$W/gsa_all.keep" "${QC[@]}" --make-just-pvar --out "$W/gsa_qc" > "$W/gsa_qc.log" 2>&1
{ printf 'step\tn\n'; printf 'lp_variants_input\t%s\n' "$(nrow "$LP.pvar")"; printf 'lp_after_qc\t%s\n' "$(nrow "$W/lp_qc.pvar")"
  if [[ -f "$GS.pvar" ]]; then printf 'gsa_variants_input\t%s\n' "$(nrow "$GS.pvar")"; else printf 'gsa_variants_input\t%s\n' "$(wc -l < "$GS.bim")"; fi
  printf 'gsa_after_qc\t%s\n' "$(nrow "$W/gsa_qc.pvar")"; } > "$W/funnel_qc.tsv"

STEP="3 allele-aligned shared SNPs"; say "$STEP"
Rscript "$HERE/match_snps.R" "$W/lp_qc.pvar" "$W/gsa_qc.pvar" "$W" > "$W/match.log" 2>&1

STEP="4 frequency concordance between technologies"; say "$STEP"
"${P2[@]}" --pfile "$LP" --keep "$W/sel_lp.keep" --extract "$W/shared_lp.ids" --freq --out "$W/lp" > "$W/lp_freq.log" 2>&1
"${P2[@]}" "${GSIN[@]}" --keep "$W/gsa_all.keep" --extract "$W/shared_gsa.ids" --freq --out "$W/gsa" > "$W/gsa_freq.log" 2>&1
Rscript "$HERE/frequency_concordance.R" "$W/lp.afreq" "$W/gsa.afreq" "$W/gsa_rename.txt" "$W" "$COH" 0.10 2>&1 | tee -a "$W/run_log.txt"

STEP="5 LD pruning (low-pass records)"; say "$STEP"
"${P2[@]}" --pfile "$LP" --keep "$W/sel_lp.keep" --extract "$W/concordant_lp.ids" --indep-pairwise 200 50 0.2 --out "$W/prune" > "$W/prune.log" 2>&1
NSNP=$(wc -l < "$W/prune.prune.in"); say "pruned to $NSNP shared SNPs"

STEP="6 array genotypes aligned to the low-pass IDs and REF alleles ($ENC)"; say "$STEP"
if [[ "$NGS" -gt 0 ]]; then
"${P2[@]}" "${GSIN[@]}" --keep "$W/gsa_all.keep" --extract "$W/shared_gsa.ids" --make-pgen --out "$W/gsa_sub" > "$W/gsa_sub.log" 2>&1
"${P2[@]}" --pfile "$W/gsa_sub" --update-name "$W/gsa_rename.txt" 2 1 --make-pgen --out "$W/gsa_ren" > "$W/gsa_ren.log" 2>&1
rm -f "$W"/gsa_sub.{pgen,pvar,psam}
"${P2[@]}" --pfile "$W/gsa_ren" --extract "$W/prune.prune.in" --ref-allele force "$W/lp_ref.txt" 2 1 --make-pgen --out "$W/gsa_hc_aligned" > "$W/gsa_hc_aligned.log" 2>&1
if [[ "$ENC" == dosage ]]; then
  bash "$HERE/array_dosage.sh" "$COH" "$W/prune.prune.in" "$W/gsa_all.keep" "$W/gsa_aligned"
  head -2000 "$W/prune.prune.in" > "$W/dchk.ids"
  "${P2[@]}" --pfile "$W/gsa_aligned" --extract "$W/dchk.ids" --export A --out "$W/dchk_dos" > /dev/null 2>&1
  "${P2[@]}" --pfile "$W/gsa_hc_aligned" --extract "$W/dchk.ids" --export A --out "$W/dchk_hc" > /dev/null 2>&1
  Rscript -e 'suppressMessages(library(data.table)); a <- commandArgs(TRUE); rd <- function(f){ d <- fread(f); m <- as.matrix(d[, -(1:6)]); rownames(m) <- d$IID; m }
    x <- rd(a[1]); y <- rd(a[2]); cn <- intersect(colnames(x), colnames(y)); rn <- intersect(rownames(x), rownames(y)); x <- x[rn, cn]; y <- y[rn, cn]
    r <- round(x); r[abs(x - r) > 0.4999] <- NA; ok <- !is.na(r) & !is.na(y)
    fwrite(data.table(cohort=a[3], records=length(rn), variants=length(cn), genotype_pairs=sum(ok), concordant=sum(r[ok] == y[ok]),
                      concordance=mean(r[ok] == y[ok]), mean_abs_dosage_minus_hardcall=mean(abs(x[ok] - y[ok]))), a[4], sep="\t")' \
    "$W/dchk_dos.raw" "$W/dchk_hc.raw" "$COH" "$W/agg_array_dosage_vs_hardcall.tsv"
  rm -f "$W"/dchk_*.raw
else
  for e in pgen pvar psam; do cp "$W/gsa_hc_aligned.$e" "$W/gsa_aligned.$e"; done
fi
NAL=$(nrow "$W/gsa_aligned.pvar"); say "array genotypes aligned on $NAL of $NSNP pruned variants"
else
  NAL=0; say "no array records in this cohort"
fi

STEP="7 block export"; say "$STEP"
rm -f "$W"/blocks/*.raw "$W"/blk_*
if [[ "$NGS" -gt 0 ]]; then grep -Fxf <(grep -v '^#' "$W/gsa_aligned.pvar" | cut -f3) "$W/prune.prune.in" > "$W/pca_variants.ids" || true
else cp "$W/prune.prune.in" "$W/pca_variants.ids"; fi
"${P2[@]}" --pfile "$LP" --keep "$W/sel_lp.keep" --extract "$W/pca_variants.ids" --make-pgen --out "$W/lp_pca" > "$W/lp_pca.log" 2>&1
"${P2[@]}" --pfile "$W/gsa_aligned" --keep "$W/gsa_all.keep" --extract "$W/pca_variants.ids" --make-pgen --out "$W/gsa_pca" > "$W/gsa_pca.log" 2>&1
split -l "$JPCA_BLOCK" -d -a 4 "$W/pca_variants.ids" "$W/blk_"
k=0
for f in "$W"/blk_*; do
  k=$((k+1))
  "${P2[@]}" --pfile "$W/lp_pca" --extract "$f" --export A --out "$W/blocks/lp_$k" > "$W/blocks/lp_$k.log" 2>&1
  "${P2[@]}" --pfile "$W/gsa_pca" --extract "$f" --export A --out "$W/blocks/gsa_$k" > "$W/blocks/gsa_$k.log" 2>&1
done
rm -f "$W"/blk_* "$W"/lp_pca.{pgen,pvar,psam} "$W"/gsa_pca.{pgen,pvar,psam}; say "exported $k blocks of up to $JPCA_BLOCK variants"
{ cat "$W/funnel_qc.tsv"; tail -n +2 "$W/match_counts.tsv" | sed 's/^/match:/'; printf 'concordant_after_freq_check\t%s\npruned\t%s\narray_aligned_pruned\t%s\nexported_for_pca\t%s\n' \
  "$(wc -l < "$W/concordant_lp.ids")" "$NSNP" "$NAL" "$(wc -l < "$W/pca_variants.ids")"; } > "$W/variant_funnel.tsv"

STEP="8 joint PCA"; say "$STEP"
Rscript "$HERE/joint_pca.R" --work "$W" --cohort "$COH" --encoding "$ENC" --previous-pcs "$JPCA_PREVIOUS_PCS" --agg "$AGG" \
  --seed "$JPCA_SEED" 2>&1 | tee -a "$W/run_log.txt"
for f in agg_selection.tsv agg_reference.tsv agg_freq_concordance.tsv variant_funnel.tsv agg_array_dosage_vs_hardcall.tsv gsa_aligned_import.tsv; do
  [[ -f "$W/$f" ]] && cp "$W/$f" "$AGG/$f"; done
[[ "${JPCA_KEEP_BLOCKS:-0}" == 1 ]] || rm -f "$W"/blocks/*.raw         # participant-level; deleted once used
echo "ok" > "$W/status.txt"; say "done"
