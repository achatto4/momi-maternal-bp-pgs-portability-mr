#!/usr/bin/env bash
# cross_cohort_genotypes.sh
#
# The variants passing the within-cohort filters in all five cohorts, pruned in each and thinned, for the descriptive
# cross-cohort PCA and the pairwise Hudson FST.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/config_pcs.sh"
ENC="$1"
X="$JPCA_OUT/work/cross"; AGG="$JPCA_OUT/aggregate/cross"; mkdir -p "$X/blocks" "$AGG"
echo "started" > "$X/status.txt"; STEP="setup"; trap 'echo "failed at: $STEP" > "$X/status.txt"' ERR
P2=( "$JPCA_PLINK2" --threads "$JPCA_THREADS" --memory "$JPCA_MEMMB" )
say(){ printf '[%s] cross: %s\n' "$(date '+%F %T')" "$*" | tee -a "$X/run_log.txt"; }
for coh in "${JPCA_COHORTS[@]}"; do [[ "$(cat "$JPCA_OUT/work/$coh/status.txt" 2>/dev/null)" == ok ]] || { echo "cohort $coh not finished" >&2; exit 3; }; done

STEP="1 variants concordant in all five cohorts with identical alleles"; say "$STEP"
for coh in "${JPCA_COHORTS[@]}"; do
  "${P2[@]}" --pfile "$JPCA_FILESETS/$coh/lpwgs_dosage_clean" --extract "$JPCA_OUT/work/$coh/concordant_lp.ids" --make-just-pvar --out "$X/lp_$coh" > "$X/lp_$coh.log" 2>&1
done
Rscript -e 'suppressMessages(library(data.table)); a <- commandArgs(TRUE); out <- a[1]; fs <- a[-1]
  L <- lapply(fs, function(f){ x <- fread(cmd=sprintf("grep -v \"^##\" %s", shQuote(f)), colClasses="character"); setnames(x, 1, "CHROM"); x[, .(ID, CHROM, POS, al=paste(REF, ALT))] })
  ids <- Reduce(intersect, lapply(L, `[[`, "ID")); A <- rbindlist(lapply(L, function(x) x[ID %chin% ids]))
  ok <- A[, .(n_al=uniqueN(al)), by=ID][n_al == 1, ID]
  fwrite(data.table(ID=ok), file.path(out, "cross_candidates.ids"), col.names=FALSE)
  fwrite(data.table(in_all_five=length(ids), identical_alleles=length(ok)), file.path(out, "cross_candidates_counts.tsv"), sep="\t")' \
  "$X" $(for coh in "${JPCA_COHORTS[@]}"; do echo "$X/lp_$coh.pvar"; done)

STEP="2 sequential LD pruning"; say "$STEP"
cp "$X/cross_candidates.ids" "$X/cur.ids"
for coh in "${JPCA_COHORTS[@]}"; do
  "${P2[@]}" --pfile "$JPCA_FILESETS/$coh/lpwgs_dosage_clean" --keep "$JPCA_OUT/work/$coh/sel_lp.keep" --extract "$X/cur.ids" \
    --indep-pairwise 200 50 0.2 --out "$X/prune_$coh" > "$X/prune_$coh.log" 2>&1
  cp "$X/prune_$coh.prune.in" "$X/cur.ids"; say "after $coh: $(wc -l < "$X/cur.ids")"
done
Rscript -e 'suppressMessages(library(data.table)); a <- commandArgs(TRUE); ids <- fread(a[1], header=FALSE)[[1]]; set.seed(as.integer(a[3]))
  mx <- as.integer(a[2]); s <- if(length(ids) > mx) sort(sample(seq_along(ids), mx)) else seq_along(ids)
  fwrite(data.table(ids[s]), a[4], col.names=FALSE)' "$X/cur.ids" "$JPCA_CROSS_MAXSNP" "$JPCA_SEED" "$X/cross.ids"
{ printf 'step\tn\n'; tail -n +2 "$X/cross_candidates_counts.tsv" | awk -F'\t' '{print "in_all_five_cohorts\t"$1"\nidentical_alleles\t"$2}'
  printf 'after_sequential_pruning\t%s\nused_after_thinning\t%s\nthinning_cap\t%s\n' "$(wc -l < "$X/cur.ids")" "$(wc -l < "$X/cross.ids")" "$JPCA_CROSS_MAXSNP"; } > "$AGG/agg_cross_variant_funnel.tsv"

STEP="3 aligned array genotypes and block export per cohort"; say "$STEP"
split -l "$JPCA_BLOCK" -d -a 4 "$X/cross.ids" "$X/blk_"
for coh in "${JPCA_COHORTS[@]}"; do
  Wc="$JPCA_OUT/work/$coh"; mkdir -p "$X/blocks/$coh"
  NG=$(grep -vc '^#' "$Wc/sel_gsa.keep" || true)
  if [[ "${NG:-0}" -gt 0 ]]; then
    if [[ "$ENC" == dosage ]]; then bash "$HERE/array_dosage.sh" "$coh" "$X/cross.ids" "$Wc/sel_gsa.keep" "$X/gsa_$coh"
    else "${P2[@]}" --pfile "$Wc/gsa_ren" --keep "$Wc/sel_gsa.keep" --extract "$X/cross.ids" --ref-allele force "$Wc/lp_ref.txt" 2 1 \
           --make-pgen --out "$X/gsa_$coh" > "$X/gsa_$coh.log" 2>&1; fi
  fi
  "${P2[@]}" --pfile "$JPCA_FILESETS/$coh/lpwgs_dosage_clean" --keep "$Wc/sel_lp.keep" --extract "$X/cross.ids" --make-pgen --out "$X/lpx_$coh" > "$X/lpx_$coh.log" 2>&1
  k=0; for f in "$X"/blk_*; do k=$((k+1))
    "${P2[@]}" --pfile "$X/lpx_$coh" --extract "$f" --export A --out "$X/blocks/$coh/lp_$k" > /dev/null 2>&1
    if [[ "${NG:-0}" -gt 0 ]]; then "${P2[@]}" --pfile "$X/gsa_$coh" --keep "$Wc/sel_gsa.keep" --extract "$f" --export A --out "$X/blocks/$coh/gsa_$k" > /dev/null 2>&1; fi
  done
  rm -f "$X"/lpx_"$coh".{pgen,pvar,psam} "$X"/gsa_"$coh".{pgen,pvar,psam}
done
rm -f "$X"/blk_*

STEP="4 cross-cohort PCA, figure and FST"; say "$STEP"
Rscript "$HERE/cross_cohort_pca.R" --rundir "$JPCA_OUT" --encoding "$ENC" --cohorts "${JPCA_COHORTS[*]}" --seed "$JPCA_SEED" 2>&1 | tee -a "$X/run_log.txt"
[[ "${JPCA_KEEP_BLOCKS:-0}" == 1 ]] || rm -rf "$X"/blocks/*/
echo "ok" > "$X/status.txt"; say "done"
