# ============================================================
# momi_env.sh — one source line to make a fresh JHPCE shell usable.
#
#   source /dcs10/chatterj/data/achattop/MOMI/code/MOMI/bp_ptb_pipeline/env/momi_env.sh
#
# Every variable below is required by some part of the pipeline, and a MISSING one fails in
# a way that looks like something else: an unset MOMI_PIPE turns `bash "$MOMI_PIPE/run_build.sh"`
# into `bash /run_build.sh` ("No such file or directory"), and an unset MOMI_OUTROOT makes the
# scoring array report "nothing to do" because it finds no filesets. Neither error names the
# real cause, which is why this file exists rather than a note in the docs.
#
# Add to ~/.bashrc if you want it automatic:
#   source /dcs10/chatterj/data/achattop/MOMI/code/MOMI/bp_ptb_pipeline/env/momi_env.sh
# ============================================================

module load conda_R/4.3 plink/2.00a4.6 2>/dev/null

export MOMI_PIPE=/dcs10/chatterj/data/achattop/MOMI/code/MOMI/bp_ptb_pipeline
export MOMI_RESULTS="$MOMI_PIPE/results/current"

# raw phenotype file (tab-delimited EPI extract)
export EPI="/dcs04/nilanjan/data/Anagh/MOMI/Epi Data/MOMI_Selected_Variables_All_Sites.txt"

# polygenic scores: weights in, per-sample scores out
export MOMI_SCORES=/dcs10/chatterj/data/achattop/MOMI/Genomics/transferability/scores
export SSC=/dcs10/chatterj/data/achattop/MOMI/Genomics/transferability/sscore

# QC'd per-cohort filesets: <OUTROOT>/<COHORT>/{gsa,lpwgs,lpwgs_dosage,infant_*}_clean
# NB docs/paths_reference.md at one point listed newdata_Jiong here; that is the RAW VCF
# location, not the fileset root. Confirmed 2026-07-20: only pipeline_out has *_clean files.
export MOMI_OUTROOT=/dcs10/chatterj/data/achattop/MOMI/Genomics/pipeline_out

# 1000 Genomes hg38 reference for PC projection (B06b/SF1/SF2)
export MOMI_G1K=/dcs04/nilanjan/data/Anagh/tools/g1k.hg38/plink_common.id

# quick self-check: report anything that does not resolve
_momi_check(){
  local bad=0
  for v in MOMI_PIPE MOMI_RESULTS MOMI_SCORES SSC MOMI_OUTROOT; do
    [[ -d "${!v}" ]] || { echo "  !! $v does not exist: ${!v}"; bad=1; }
  done
  [[ -f "$EPI" ]] || { echo "  !! EPI not found: $EPI"; bad=1; }
  [[ -f "$MOMI_G1K.bed" ]] || echo "  ,, MOMI_G1K.bed absent (only B06b/SF1/SF2 need it)"
  [[ $bad -eq 0 ]] && echo "momi env ok — pipe=$MOMI_PIPE"
}
_momi_check
