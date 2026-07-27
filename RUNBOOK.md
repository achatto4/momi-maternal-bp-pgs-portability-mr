# RUNBOOK — run order on JHPCE

All commands assume you've `git pull`ed the latest in the repo. Set these once per shell:

```bash
cd /dcs10/chatterj/data/achattop/MOMI/code/MOMI && git pull && cd bp_ptb_pipeline
REPO=$(pwd)
OUTROOT=/dcs10/chatterj/data/achattop/MOMI/Genomics/pipeline_out
NEW=/dcs10/chatterj/data/achattop/MOMI/Genomics/newdata_Jiong
EPI="/dcs04/nilanjan/data/Anagh/MOMI/Epi Data/MOMI_Selected_Variables_All_Sites.txt"
SUM=/dcs04/nilanjan/data/Anagh/MOMI/sumstat/BP
ANALYSIS=$OUTROOT/analysis_gsa
declare -A PFX=( [AMANHI-Bangladesh]=AMANHIB [AMANHI-Pakistan]=AMANHIP \
                 [AMANHI-Pemba]=AMANHIT [GAPPS-Bangladesh]=GAPPSB [GAPPS-Zambia]=ZAPPS )
```

Status: stage 00 (GSA merge), 01 (phenotypes), 02 (PCA), 05 (PRS-MR) are DONE on `analysis_gsa`.
This runbook covers (1) GWAS + SNP-MR on GSA, then (2) lpWGS + final all-10 merge + full rerun.

---

## PART 1 — GWAS → SNP-based MR (GSA)

### 1.1 Internal GWAS (SBP/DBP linear, PTB logistic; PW_AGE + 10 PCs + site dummies)
```bash
sbatch --job-name=gwas_gsa --cpus-per-task=8 --mem=64G --time=4:00:00 \
  --output="$ANALYSIS/gwas_%j.out" \
  --wrap "module load plink/2.00a4.6; module load R; \
bash 03_gwas/run_gwas.sh --merged '$OUTROOT/merged_gsa' --keep '$ANALYSIS/mothers.keep' \
  --pheno-dir '$ANALYSIS' --covar '$ANALYSIS/covariate_table_analytic_mothers.txt' \
  --eigenvec '$ANALYSIS/pca_analysis/pca_results.eigenvec' --out-dir '$ANALYSIS' \
  --bp-window latest --n-pcs 10 --threads 8 --mem-mb 60000"
```
Check when done:
```bash
cat "$ANALYSIS"/gwas_*.out | tail -20
wc -l "$ANALYSIS"/results_*/gwas_*.glm.* 2>/dev/null
```
Expect: results_SBP/DBP `.glm.linear` and results_PTB_NEW `.glm.logistic.hybrid`.

### 1.2 Build MR instruments from external EUR BP GWAS (clump vs merged)
```bash
mkdir -p "$ANALYSIS/mr"
sbatch --job-name=instr_gsa --cpus-per-task=4 --mem=24G --time=2:00:00 \
  --output="$ANALYSIS/mr/instr_%j.out" \
  --wrap "module load plink/2.00a4.6; \
bash 04_mr/prep_instruments.sh --ext '$SUM/mill_EUR/SBP/GCST90310294.h.tsv.gz' --trait SBP --merged '$OUTROOT/merged_gsa' --out-dir '$ANALYSIS/mr'; \
bash 04_mr/prep_instruments.sh --ext '$SUM/mill_EUR/DBP/GCST90310295.h.tsv.gz' --trait DBP --merged '$OUTROOT/merged_gsa' --out-dir '$ANALYSIS/mr'"
```
Check: `wc -l "$ANALYSIS"/mr/instruments_*.txt` (expect tens–hundreds of independent instruments).

### 1.3 Two-sample MR (BP exposure -> PTB outcome)
> Needs the R package **TwoSampleMR**. First check it's available:
> `module load R; Rscript -e 'library(TwoSampleMR)'` — if it errors, install once:
> `Rscript -e 'install.packages("remotes"); remotes::install_github("MRCIEU/TwoSampleMR")'`
```bash
sbatch --job-name=mr_gsa --cpus-per-task=2 --mem=16G --time=1:00:00 \
  --output="$ANALYSIS/mr/mr_%j.out" \
  --wrap "module load R; \
Rscript 04_mr/run_mr.R --instruments '$ANALYSIS/mr/instruments_SBP.txt' --ptb-gwas '$ANALYSIS/results_PTB_NEW/gwas_PTB_NEW.PHENO.glm.logistic.hybrid' --trait SBP --out-dir '$ANALYSIS/mr'; \
Rscript 04_mr/run_mr.R --instruments '$ANALYSIS/mr/instruments_DBP.txt' --ptb-gwas '$ANALYSIS/results_PTB_NEW/gwas_PTB_NEW.PHENO.glm.logistic.hybrid' --trait DBP --out-dir '$ANALYSIS/mr'"
```
Check: `cat "$ANALYSIS"/mr/mr_*.out` and `column -t "$ANALYSIS"/mr/mr_results_*_PTB.txt`.

---

## PART 2 — lpWGS, then final all-10 merge, then full rerun

### 2.0 Verify the big lpWGS files (Pemba & GAPPS-Bangladesh were split & concatenated)
```bash
module load bcftools 2>/dev/null || true
for c in AMANHI-Pemba GAPPS-Bangladesh AMANHI-Bangladesh AMANHI-Pakistan GAPPS-Zambia; do
  f=$NEW/$c/lpwgs_imputed.vcf.gz
  echo "=== $c ==="; ls -lh "$f"*
  bcftools view -h "$f" >/dev/null 2>&1 && echo "header OK" || echo "HEADER READ FAILED (possible truncated cat)"
done
```
If any "HEADER READ FAILED", re-concatenate that site's parts: `cat lpwgs_imputed.vcf.gz* > lpwgs_imputed.vcf.gz` (parts named gz0..gz4), then re-index `bcftools index -t`.

### 2.1 lpWGS conversion — TEST ONE COHORT FIRST (smallest = GAPPS-Zambia, ~122 GB)
> The lpWGS path (`bcftools view -f PASS | plink2 --bcf`) is NEW/untested. Validate on the
> smallest cohort before launching all five. lpWGS jobs are heavy: give them more time/mem.
```bash
bash 00_preprocess/submit_convert.sh \
  --cohort-dir "$NEW/GAPPS-Zambia" --platform lpwgs \
  --pid-prefix ZAPPS --epi "$EPI" --out "$OUTROOT/GAPPS-Zambia" \
  --r2-min 0.8 --cpus 8 --mem 96G --time 24:00:00 --force
squeue -u achattop
```
When done, sanity-check missingness (should be ~0 like GSA):
```bash
sbatch --job-name=lpchk --cpus-per-task=2 --mem=16G --time=0:30:00 --output=/tmp/lpchk_%j.out \
  --wrap "module load plink/2.00a4.6; plink2 --bfile '$OUTROOT/GAPPS-Zambia/lpwgs_clean' --missing --out /tmp/lpz"
# then: awk 'NR>1{s+=$NF;n++}END{print \"mean F_miss\",s/n}' /tmp/lpz.smiss   (plink2 .smiss: last col is F_MISS)
```
If the bcftools|plink2 pipe fails or is too slow, ping me — we may pre-filter to a temp BCF per chromosome instead.

### 2.2 Convert the remaining 4 lpWGS cohorts
```bash
for c in AMANHI-Bangladesh AMANHI-Pakistan AMANHI-Pemba GAPPS-Bangladesh; do
  bash 00_preprocess/submit_convert.sh \
    --cohort-dir "$NEW/$c" --platform lpwgs \
    --pid-prefix "${PFX[$c]}" --epi "$EPI" --out "$OUTROOT/$c" \
    --r2-min 0.8 --cpus 8 --mem 96G --time 24:00:00 --force
done
squeue -u achattop
# check: for c in <all 5>; do echo "$c: $(wc -l < $OUTROOT/$c/lpwgs_clean.fam) mothers"; done
```

### 2.3 Build the FINAL all-10 merge (5 GSA + 5 lpWGS) -> merged_all
```bash
sbatch --job-name=merge_all --cpus-per-task=4 --mem=48G --time=6:00:00 \
  --output="$OUTROOT/merge_all_%j.out" \
  --wrap "module load plink/1.90b; \
bash 00_preprocess/merge_filesets.sh --out '$OUTROOT/merged_all' --plink1 plink \
  --geno 0.05 --mind 0.05 --maf 0.01 \
  '$OUTROOT/AMANHI-Bangladesh/gsa_clean'  '$OUTROOT/AMANHI-Pakistan/gsa_clean' \
  '$OUTROOT/AMANHI-Pemba/gsa_clean'       '$OUTROOT/GAPPS-Bangladesh/gsa_clean' \
  '$OUTROOT/GAPPS-Zambia/gsa_clean' \
  '$OUTROOT/AMANHI-Bangladesh/lpwgs_clean' '$OUTROOT/AMANHI-Pakistan/lpwgs_clean' \
  '$OUTROOT/AMANHI-Pemba/lpwgs_clean'      '$OUTROOT/GAPPS-Bangladesh/lpwgs_clean' \
  '$OUTROOT/GAPPS-Zambia/lpwgs_clean'"
```
Check the tail: common-SNP count (will be smaller than GSA-only since GSA∩lpWGS intersect) and total samples (~all GSA + lpWGS mothers). If the common-SNP count is too thin, ping me — we'd switch to per-platform analysis + meta-analysis instead of one pooled set.

### 2.4 Full rerun on the combined set (new analysis dir)
```bash
ALL=$OUTROOT/analysis_all; mkdir -p "$ALL"
# Stage 01 phenotypes
sbatch --job-name=ph_all --mem=32G --time=1:00:00 --output=$ALL/pheno_%j.out \
  --wrap "module load R; Rscript 01_phenotypes/build_phenotypes.R --merged-fam '$OUTROOT/merged_all.fam' --epi '$EPI' --out-dir '$ALL' --ptb-case-code 2"
# Stage 02 PCA  (after pheno done)
sbatch --job-name=pca_all --cpus-per-task=4 --mem=32G --time=2:00:00 --output=$ALL/pca_%j.out \
  --wrap "module load plink/2.00a4.6; bash 02_pca/run_pca.sh --merged '$OUTROOT/merged_all' --keep '$ALL/mothers.keep' --out-dir '$ALL' --n-pcs 10 --threads 4 --mem-mb 28000 && { module load R; Rscript 02_pca/plot_pca.R --eigenvec '$ALL/pca_analysis/pca_results.eigenvec' --eigenval '$ALL/pca_analysis/pca_results.eigenval' --covar '$ALL/covariate_table_analytic_mothers.txt' --out-dir '$ALL/pca_analysis'; }"
# Stage 05 PRS  (compute 4 scores, then analysis) — same commands as analysis_gsa but with $ALL and merged_all
# Stage 03 GWAS + Stage 04 MR — same as PART 1 but with $ALL and merged_all
```
(For 2.4 PRS/GWAS/MR, reuse the PART-1 / stage-05 commands, swapping `$ANALYSIS`->`$ALL` and `merged_gsa`->`merged_all`.)

---

## Notes / known risks to watch
- **TwoSampleMR** may need a one-time install (see 1.3).
- **lpWGS conversion** path is untested — validate on GAPPS-Zambia first (2.1).
- **Combined common-SNP count** (2.3): GSA and lpWGS overlap may be modest; if too small, we pivot to per-platform GWAS/PRS + meta-analysis (tell me the number).
- plink2 `.smiss` per-sample missingness column is `F_MISS` (last col); plink1.9 used `.imiss`.
