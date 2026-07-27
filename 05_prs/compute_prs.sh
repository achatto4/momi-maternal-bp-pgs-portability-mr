#!/usr/bin/env bash
# ============================================================
# compute_prs.sh
# Compute ONE polygenic score on the merged genotypes.
#   1) parse a PGS Catalog scoring file -> weights (variant_id, effect_allele, weight)
#      variant_id = chr<hm_chr>:<hm_pos>  to match our 'chrN:pos' bim IDs (GRCh38)
#   2) plink2 --score  -> <out>.sscore  (column SCORE1_AVG)
#
# Column detection (by NAME from the header, robust to column order):
#   chr   : hm_chr (preferred, GRCh38) else chr_name
#   pos   : hm_pos (preferred, GRCh38) else chr_position
#   allele: effect_allele
#   weight: effect_weight
# Override with --chr-col/--pos-col/--ea-col/--w-col if needed.
#
# Usage:
#   compute_prs.sh --pgs FILE.txt.gz --bfile PREFIX --keep FILE --out PREFIX
#                  [--chr-prefix chr] [--plink2 plink2] [--threads 4] [--mem-mb 16000]
# ============================================================
set -euo pipefail
PGS=""; BFILE=""; KEEP=""; OUT=""; CHRPFX="chr"
CHRCOL=""; POSCOL=""; EACOL="effect_allele"; WCOL="effect_weight"
PLINK2=plink2; THREADS=4; MEMMB=16000
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pgs) PGS="$2"; shift 2;;
    --bfile) BFILE="$2"; shift 2;;
    --keep) KEEP="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --chr-prefix) CHRPFX="$2"; shift 2;;
    --chr-col) CHRCOL="$2"; shift 2;;
    --pos-col) POSCOL="$2"; shift 2;;
    --ea-col) EACOL="$2"; shift 2;;
    --w-col) WCOL="$2"; shift 2;;
    --plink2) PLINK2="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mem-mb) MEMMB="$2"; shift 2;;
    -h|--help) echo "see header"; exit 1;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$PGS" && -n "$BFILE" && -n "$OUT" ]] || { echo "need --pgs --bfile --out (--keep optional)" >&2; exit 1; }
[[ -f "$PGS" ]] || { echo "PGS file not found: $PGS" >&2; exit 2; }
command -v "$PLINK2" >/dev/null || { echo "plink2 not on PATH (module load?)" >&2; exit 3; }
mkdir -p "$(dirname "$OUT")"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

# accept either a bed (hard calls) or a pgen (dosages) prefix; --pfile uses dosages
GIN=(--bfile "$BFILE"); [[ -f "$BFILE.pgen" ]] && GIN=(--pfile "$BFILE")
[[ -f "$BFILE.bed" || -f "$BFILE.pgen" ]] || { echo "no .bed or .pgen at prefix: $BFILE" >&2; exit 2; }

WEIGHTS="$OUT.weights.txt"
log "parsing $PGS -> weights ($CHRPFX{chr}:{pos})"
zcat "$PGS" | awk -v chrpfx="$CHRPFX" -v chrcol="$CHRCOL" -v poscol="$POSCOL" -v eacol="$EACOL" -v wcol="$WCOL" '
  BEGIN { FS="\t" }    # PGS Catalog files are TAB-delimited; preserve empty fields
  /^#/ { next }
  hdr==0 {
    hdr=1
    for(i=1;i<=NF;i++) col[$i]=i
    # choose chr/pos columns
    c = chrcol!="" ? chrcol : ("hm_chr" in col ? "hm_chr" : "chr_name")
    p = poscol!="" ? poscol : ("hm_pos" in col ? "hm_pos" : "chr_position")
    if(!(c in col) || !(p in col) || !(eacol in col) || !(wcol in col)){
      print "ERROR: required columns not found. Have: " > "/dev/stderr"
      for(k in col) printf " %s", k > "/dev/stderr"; print "" > "/dev/stderr"
      exit 3
    }
    ci=col[c]; pi=col[p]; ai=col[eacol]; wi=col[wcol]
    next
  }
  {
    chr=$ci; pos=$pi; ea=$ai; w=$wi
    if(chr=="" || pos=="" || ea=="" || w=="" || w=="NA") next
    sub(/^chr/,"",chr)             # normalize, then re-add prefix
    print chrpfx chr ":" pos "\t" ea "\t" w
  }
' > "$WEIGHTS"

NW=$(wc -l < "$WEIGHTS")
log "weights: $NW variants"
[[ "$NW" -gt 0 ]] || { echo "ERROR: 0 weights parsed" >&2; exit 4; }

# ---- DROP VARIANTS WHOSE chr:pos APPEARS MORE THAN ONCE (added 2026-07-20) ----
# Our variant IDs are chr:pos with NO alleles, because that is what the bim files carry.
# A multi-allelic site, or the same position listed twice in the scoring file, therefore
# collapses to one ID with two different effect alleles, and plink2 refuses outright:
#   "Error: REF allele for variant 'chr1:2293437' appears multiple times in --score file."
# This is why PGS003964/PGS003968/PGS005008 (Kurniansyah 2023, Gunn 2024) were never scored:
# they are ~1.3M-variant scores that include duplicated positions. The eight scores already
# in the paper have none, so this step is a NO-OP for them and cannot change any locked result.
#
# WE DROP RATHER THAN PICK. Keeping the first occurrence would assign a weight to whichever
# allele happened to come first in the file, and since the ID carries no allele we cannot
# check that against the genotypes -- a silently wrong sign on those variants. Dropping loses
# a little signal; guessing risks corrupting it. The count is logged so the loss is visible,
# and it flows into S1's per-score variant counts.
DUPS=$(cut -f1 "$WEIGHTS" | sort | uniq -d | wc -l)
if [[ "$DUPS" -gt 0 ]]; then
  awk 'NR==FNR{c[$1]++; next} c[$1]==1' "$WEIGHTS" "$WEIGHTS" > "$WEIGHTS.dedup"
  mv "$WEIGHTS.dedup" "$WEIGHTS"
  ND=$(wc -l < "$WEIGHTS")
  log "duplicate chr:pos IDs: $DUPS distinct -> dropped ALL their rows; weights $NW -> $ND"
  NW="$ND"
  [[ "$NW" -gt 0 ]] || { echo "ERROR: every weight was a duplicate" >&2; exit 4; }
fi

# plink2 --score: col1=ID, col2=allele, col3=weight. --keep optional (filesets are already mothers-only).
# On a pgen the score uses dosages automatically (low-pass WGS power retained).
log "plink2 --score ${GIN[0]}${KEEP:+ (keep $KEEP)}"
"$PLINK2" "${GIN[@]}" ${KEEP:+--keep "$KEEP"} \
  --score "$WEIGHTS" 1 2 3 cols=nallele,denom,scoresums,scoreavgs \
  --threads "$THREADS" --memory "$MEMMB" \
  --out "$OUT" > "$OUT.score.log" 2>&1 || { echo "plink2 --score failed; see $OUT.score.log"; tail -5 "$OUT.score.log"; exit 5; }

# overlap = variants actually used
USED=$(grep -oE 'valid predictors loaded|--score: [0-9]+ variant' "$OUT.score.log" | tail -1 || true)
log "DONE: $OUT.sscore  ($(awk 'END{print NR-1}' "$OUT.sscore") samples scored). plink note: $USED"
head -2 "$OUT.sscore"
