#!/usr/bin/env bash
# ============================================================
# merge_filesets.sh
# Merge N QC'd PLINK1 filesets (cohort x platform 'clean' beds) into one.
# Samples are NON-OVERLAPPING across inputs -> this is a SNP-intersection +
# sample-stack. Variant IDs are uniform chr:pos (set in stage-00 conversion).
#
# Steps:
#   1) common SNP set = variant IDs present in ALL inputs
#   2) extract common SNPs from each input
#   3) plink1.9 --merge-list ; on allele/strand mismatch:
#        retry after --flip of the .missnp in every input,
#        then (if still failing) --exclude the .missnp from every input
#   4) final QC: --geno (variant missingness), --mind (sample missingness), --maf
#
# Usage:
#   merge_filesets.sh --out PREFIX [--plink1 plink] [--workdir DIR]
#                     [--geno 0.05] [--mind 0.05] [--maf 0.01]
#                     BED1 BED2 [BED3 ...]
# (BEDi = plink1 fileset prefix, i.e. path without .bed/.bim/.fam)
# ============================================================
set -euo pipefail

PLINK1=plink; OUT=""; WORKDIR=""; GENO=0.05; MIND=0.05; MAF=0.01
BEDS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out) OUT="$2"; shift 2;;
    --plink1) PLINK1="$2"; shift 2;;
    --workdir) WORKDIR="$2"; shift 2;;
    --geno) GENO="$2"; shift 2;;
    --mind) MIND="$2"; shift 2;;
    --maf) MAF="$2"; shift 2;;
    -h|--help) echo "see header"; exit 1;;
    --*) echo "Unknown arg: $1" >&2; exit 1;;
    *) BEDS+=("$1"); shift;;
  esac
done
[[ -n "$OUT" ]] || { echo "need --out" >&2; exit 1; }
[[ ${#BEDS[@]} -ge 2 ]] || { echo "need >=2 input bed prefixes" >&2; exit 1; }
command -v "$PLINK1" >/dev/null || { echo "plink1 not on PATH (module load?)" >&2; exit 3; }
for b in "${BEDS[@]}"; do [[ -f "$b.bed" && -f "$b.bim" && -f "$b.fam" ]] || { echo "missing fileset: $b" >&2; exit 2; }; done

mkdir -p "$(dirname "$OUT")"
[[ -n "$WORKDIR" ]] || WORKDIR="$(dirname "$OUT")/merge_work"
mkdir -p "$WORKDIR"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
N=${#BEDS[@]}
log "merging $N filesets -> $OUT  (geno=$GENO mind=$MIND maf=$MAF)"

# ---- 1) common SNPs (present in all N) ----
: > "$WORKDIR/all_ids.txt"
for b in "${BEDS[@]}"; do awk '{print $2}' "$b.bim" | sort -u >> "$WORKDIR/all_ids.txt"; done
sort "$WORKDIR/all_ids.txt" | uniq -c | awk -v n="$N" '$1==n{print $2}' > "$WORKDIR/common.snps"
NCOM=$(wc -l < "$WORKDIR/common.snps")
log "common SNPs across all $N inputs: $NCOM"
[[ "$NCOM" -gt 0 ]] || { echo "ERROR: no common SNPs — check variant-ID consistency" >&2; exit 4; }

# ---- 2) extract common SNPs from each input ----
extract_all() {  # $1 = suffix tag, $2... = extra plink args (e.g. --flip file / --exclude file)
  local tag="$1"; shift
  : > "$WORKDIR/mergelist_$tag.txt"
  local i=0
  for b in "${BEDS[@]}"; do
    "$PLINK1" --bfile "$b" --extract "$WORKDIR/common.snps" "$@" \
      --make-bed --out "$WORKDIR/part_${tag}_$i" >/dev/null 2>&1
    [[ $i -gt 0 ]] && echo "$WORKDIR/part_${tag}_$i" >> "$WORKDIR/mergelist_$tag.txt"
    i=$((i+1))
  done
}

try_merge() {  # $1 = tag ; merges part_${tag}_0 + mergelist_${tag}
  local tag="$1"
  "$PLINK1" --bfile "$WORKDIR/part_${tag}_0" \
    --merge-list "$WORKDIR/mergelist_$tag.txt" \
    --make-bed --out "$WORKDIR/merged_$tag" > "$WORKDIR/merge_$tag.log" 2>&1 || true
}

# ---- 3) merge with flip/exclude retries ----
log "attempt 1: straight merge"
extract_all a
try_merge a
MERGED="$WORKDIR/merged_a"

if [[ -f "$WORKDIR/merged_a-merge.missnp" ]]; then
  NM=$(wc -l < "$WORKDIR/merged_a-merge.missnp")
  log "attempt 2: $NM mismatching variants -> flip them in all inputs and retry"
  extract_all b --flip "$WORKDIR/merged_a-merge.missnp"
  try_merge b
  MERGED="$WORKDIR/merged_b"
  if [[ -f "$WORKDIR/merged_b-merge.missnp" ]]; then
    NM2=$(wc -l < "$WORKDIR/merged_b-merge.missnp")
    log "attempt 3: $NM2 still mismatching -> exclude them from all inputs"
    cat "$WORKDIR/merged_a-merge.missnp" "$WORKDIR/merged_b-merge.missnp" | sort -u > "$WORKDIR/drop.missnp"
    extract_all c --exclude "$WORKDIR/drop.missnp"
    try_merge c
    MERGED="$WORKDIR/merged_c"
  fi
fi

[[ -f "$MERGED.bed" ]] || { echo "ERROR: merge failed; see $WORKDIR/merge_*.log" >&2; exit 5; }
log "merged (pre-QC): $(wc -l < "$MERGED.bim") variants x $(wc -l < "$MERGED.fam") samples"

# ---- 4) final QC ----
log "final QC: --geno $GENO --mind $MIND --maf $MAF"
"$PLINK1" --bfile "$MERGED" --geno "$GENO" --mind "$MIND" --maf "$MAF" \
  --make-bed --out "$OUT" > "$OUT.qc.log" 2>&1

[[ -f "$OUT.bed" ]] || { echo "ERROR: final QC failed; see $OUT.qc.log" >&2; exit 6; }
log "DONE: $(wc -l < "$OUT.bim") variants x $(wc -l < "$OUT.fam") samples -> $OUT.{bed,bim,fam}"
log "per-cohort sample counts in merged fam (by PARTICIPANT_ID prefix):"
awk '{print $2}' "$OUT.fam" | sed 's/[0-9].*//' | sort | uniq -c
