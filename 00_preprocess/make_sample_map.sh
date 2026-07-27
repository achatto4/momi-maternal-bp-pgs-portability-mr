#!/usr/bin/env bash
# ============================================================
# make_sample_map.sh  (v2)
# Build, for ONE cohort x platform, the lists stage-00 conversion needs:
#   <out>.keep     PLINK --keep (FID IID) of MOTHER samples, as named in the VCF
#   <out>.rename   PLINK --update-ids: oldFID oldIID newFID newIID
#                  (rename each mother VCF sample -> Epi PARTICIPANT_ID)
#   <out>.report   counts (mothers matched / dropped dup / suspected-but-unmatched)
#
# MOTHER RULE (validated across all 5 cohorts x both platforms):
#   A) Subject_ID is itself an Epi PARTICIPANT_ID         -> mother; PID = Subject_ID
#      (GSA x5 ; lpWGS for AMANHI-* and GAPPS-Bangladesh)
#   B) else Subject_ID ends in "-M" (Zambia lpWGS style)  -> mother; strip "-M":
#        - if stripped id is an Epi ORIG_ID  -> PID = orig2pid[stripped]   (primary)
#        - elif "<prefix>-<stripped>" is an Epi PARTICIPANT_ID -> that PID  (fallback)
#   Everything else (children, non-Epi samples) is excluded.
#
# Linkage is to Epi PARTICIPANT_ID, restricted to this cohort's --pid-prefix so
# ORIG_IDs can't collide across cohorts. Canonical analysis IID = PARTICIPANT_ID.
#
# idmap.csv columns: id,sid,sex,ORIG_ID,Subject_ID   (sid/ORIG_ID may be unused/empty)
# Epi (tab): col1=SITE_CODE col2=ORIG_ID col3=PARTICIPANT_ID ...
# Pair the output with plink2 --double-id (FID=IID=VCF sample name).
# ============================================================
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $0 --idmap FILE --epi FILE --pid-prefix STR --out PREFIX
          [--idcol 1] [--origcol 4] [--subjcol 5]
USAGE
  exit 1
}

IDCOL=1; ORIGCOL=4; SUBJCOL=5
IDMAP=""; EPI=""; PREFIX=""; OUT=""; MODE="mother"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --idmap) IDMAP="$2"; shift 2;;
    --epi) EPI="$2"; shift 2;;
    --pid-prefix) PREFIX="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --idcol) IDCOL="$2"; shift 2;;
    --origcol) ORIGCOL="$2"; shift 2;;
    --subjcol) SUBJCOL="$2"; shift 2;;
    --mode) MODE="$2"; shift 2;;      # mother (default) | infant  -- ADDED 2026-07-20
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done
[[ "$MODE" == "mother" || "$MODE" == "infant" ]] || { echo "--mode must be mother|infant" >&2; exit 1; }

# ---- INFANT MODE (added 2026-07-20 for S15/B32, fetal-vs-maternal decomposition) ----
# The mother rule below was validated across all 5 cohorts x both platforms and is left
# BYTE-IDENTICAL in behaviour when --mode mother (the default), so nothing already built
# can shift. Infant mode mirrors it exactly, swapping two things:
#   target set : Epi BABY_ID instead of PARTICIPANT_ID
#   suffix     : "-C" instead of "-M"   (the idmap encodes child samples as ...-C)
#
# WHY THIS IS WORTH THE EDIT RATHER THAN A NEW SCRIPT. The infant filesets must go through
# EXACTLY the QC the maternal ones did, or the fetal and maternal scores are not comparable
# and the decomposition is meaningless. Sharing the sample-map logic is the cheapest way to
# guarantee that; a parallel implementation would drift.
#
# BABY_ID column is located BY NAME from the Epi header rather than hardcoded, because the
# mother path's positional assumption (col3=PARTICIPANT_ID) was validated and BABY_ID's
# position was not. Fails loudly if absent.
BABYCOL=0
if [[ "$MODE" == "infant" ]]; then
  BABYCOL=$(head -1 "$EPI" | tr '\t' '\n' | grep -nx "BABY_ID" | cut -d: -f1 || true)
  [[ -n "$BABYCOL" && "$BABYCOL" -gt 0 ]] || { echo "ERROR: no BABY_ID column in $EPI" >&2; exit 2; }
  echo "[make_sample_map] infant mode: BABY_ID is Epi column $BABYCOL"
fi
[[ -n "$IDMAP" && -n "$EPI" && -n "$PREFIX" && -n "$OUT" ]] || usage
[[ -f "$IDMAP" ]] || { echo "idmap not found: $IDMAP" >&2; exit 2; }
[[ -f "$EPI"   ]] || { echo "epi not found: $EPI" >&2; exit 2; }
mkdir -p "$(dirname "$OUT")"
: > "$OUT.report"   # truncate (awk appends with >>)

awk -v prefix="$PREFIX" -v idc="$IDCOL" -v origc="$ORIGCOL" -v subjc="$SUBJCOL" \
    -v keepf="$OUT.keep" -v renamef="$OUT.rename" -v repf="$OUT.report" \
    -v mode="$MODE" -v babyc="$BABYCOL" \
    '
  # suffix is -M for mothers, -C for children; strip() is used by rule B
  function strip(s){ sub(sfx "$","",s); return s }
  BEGIN { sfx = (mode=="infant") ? "-C" : "-M"
          who = (mode=="infant") ? "infant" : "mother" }
  # ---- pass 1: Epi (tab) for this cohort ----
  FNR==NR {
    if (FNR==1) next
    pid=$3; orig=$2
    if (pid ~ ("^" prefix)) {
      pidset[pid]=1
      if (orig!="" && orig!="NA") {
        if ((orig in orig2pid) && orig2pid[orig]!=pid) conflict[orig]=1
        orig2pid[orig]=pid
      }
      ## infant mode: the target identifier is BABY_ID, and we remember which mother it
      ## belongs to so the decomposition can pair them without a second lookup.
      if (mode=="infant" && babyc>0) {
        bid=$(babyc)
        if (bid!="" && bid!="NA") { babyset[bid]=1; baby2mom[bid]=pid; mom2baby[pid]=bid }
      }
    }
    next
  }
  # ---- pass 2: idmap (csv) ----
  {
    if (FNR==1) next
    n=split($0, f, ",")
    id=f[idc]; orig=f[origc]; subj=f[subjc]
    gsub(/^[ \t"]+|[ \t"]+$/,"",id);  gsub(/^[ \t"]+|[ \t"]+$/,"",orig); gsub(/^[ \t"]+|[ \t"]+$/,"",subj)
    ntot++
    mother=0; pid=""
    if (mode=="mother") {
      if (subj in pidset) { mother=1; pid=subj }                 # rule A
      else if (subj ~ /-M$/) {                                   # rule B (Zambia lpWGS)
        key=strip(subj)
        if (key in orig2pid)                { mother=1; pid=orig2pid[key] }
        else if ((prefix "-" key) in pidset) { mother=1; pid=prefix "-" key }
      }
    } else {
      ## INFANT: mirror of the above. Rule A is a direct BABY_ID hit; rule B handles the
      ## Zambia-lpWGS "-C" style. Note there is no ORIG_ID equivalent for babies, so rule B
      ## can only try the prefixed form -- any child sample not resolvable that way is
      ## reported as UNMATCHED rather than guessed at.
      if (subj in babyset) { mother=1; pid=subj }                # rule A
      else if (subj ~ /-C$/) {                                   # rule B
        key=strip(subj)
        if (key in babyset)                 { mother=1; pid=key }
        else if ((prefix "-" key) in babyset){ mother=1; pid=prefix "-" key }
        ## RULE C, added 2026-07-20 for GAPPS-Zambia lpWGS.
        ## NB no apostrophes anywhere in this comment: the whole awk program is a single
        ## quoted string, so one apostrophe terminates it and bash parses the remainder.
        ## Zambia encodes samples as <mother-ORIG_ID>-M / <mother-ORIG_ID>-C, i.e. a CHILD
        ## sample carries the MATERNAL identifier rather than its own. Rules A and B both
        ## look for a BABY_ID belonging to the child, and so matched 0 of 439 -C samples.
        ## The link is: strip -C -> maternal ORIG_ID -> maternal PARTICIPANT_ID -> BABY_ID.
        ## Composes the orig2pid map the mother rule already uses with the mom2baby map
        ## built in pass 1, so it adds no new assumption about the data.
        else {
          if (key in orig2pid && (orig2pid[key] in mom2baby)) {
            mother=1; pid=mom2baby[orig2pid[key]]
          } else if (((prefix "-" key) in mom2baby)) {
            mother=1; pid=mom2baby[prefix "-" key]
          }
        }
      }
    }
    if (!mother) {
      # looked like the target (carried the suffix) but had no Epi link -> report
      if (subj ~ (sfx "$") || id ~ (sfx "$")) { susp_unmatched++; print "UNMATCHED suspected-" who ": id=" id " Subject_ID=" subj > repf }
      else nonmother++
      next
    }
    nmother++
    if (pid in seenpid) {
      dropped_dup++
      print "DROPPED dup-PID: id=" id " PID=" pid " (kept " seenpid[pid] ")" > repf
      next
    }
    seenpid[pid]=id
    print id, id           > keepf
    print id, id, pid, pid > renamef
    ## infant mode also emits the pairing, so B32 never has to re-derive who belongs to whom
    if (mode=="infant") print pid, baby2mom[pid] > (keepf ".pairs")
    kept++
  }
  END {
    print "=== make_sample_map (v2) report ===" >> repf
    print "mode                       : " who >> repf
    print "cohort PID prefix          : " prefix >> repf
    print "Epi PARTICIPANT_IDs        : " length(pidset) >> repf
    if (mode=="infant") print "Epi BABY_IDs               : " length(babyset) >> repf
    print "idmap rows                 : " (ntot+0) >> repf
    print who "s identified            : " (nmother+0) >> repf
    print "  unique " who "s kept      : " (kept+0) >> repf
    print "  dup-ID dropped           : " (dropped_dup+0) >> repf
    print "suspected-" who " UNMATCHED : " (susp_unmatched+0) " (carried " sfx " but no Epi link; usually no phenotype)" >> repf
    print "other rows (excluded)      : " (nonmother+0) >> repf
    print "ORIG_ID->PID conflicts     : " (length(conflict)+0) >> repf
    print "kept .keep/.rename rows    : " (kept+0) >> repf
  }
  ' "$EPI" "$IDMAP"

echo "[make_sample_map] $OUT  -> kept $(wc -l < "$OUT.keep" 2>/dev/null || echo 0) ${MODE}s"
tail -n 13 "$OUT.report"
