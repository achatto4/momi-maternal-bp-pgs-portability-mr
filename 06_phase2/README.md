# Stage 06 — Phase-2 genotype / scoring analyses

Genotype-level analyses that go beyond the epi + sscore files used in Phase 1.
They need `plink2`, the merged genotypes, the stage-02 PCA eigenvectors, and the
per-cohort clean filesets. Run everything with `run_phase2.sh` (edit its PATHS
block), or call the scripts individually.

| ID | Question | Scripts | Key output |
|----|----------|---------|-----------|
| **A1** | Does transferability R² track genetic distance (Fst / PC) between cohorts? | `a1_fst_pcdist.sh` → `a1_fst_pcdist.R` | `a1_distance_vs_r2.{tsv,png}` |
| **A2** | Is weak transfer explained by missing score SNPs / low imputation INFO / divergent allele frequencies? | `a2_scoresnp_freq.sh` → `a2_scoresnp_summary.R` | `a2_scoresnp_summary.tsv` |
| **A3** | Is the R² robust to population structure (PC adjustment) and cryptic relatedness (KING)? | `a3_king.sh` → `a3_pcadj_r2.R` | `a3_pcadj_r2.tsv` |
| **C1** | Is the BP→fetal-growth effect independent of maternal adiposity? | `c1_bmi_score.sh` → `c1_bmi_mvmr.R` | printed MVMR table |
| **C2** | How much of the "maternal" birth-weight effect is fetal-genotype transmission? | `c2_fetal_score.sh` → `c2_fetal_maternal.R` | printed maternal/fetal decomposition |

## Inputs (produced by earlier stages)
- `--merged` : merged all-cohort `.bed/.bim/.fam` (or `.pgen`); IID = `PARTICIPANT_ID`.
- `--analysis-dir` : contains `pca_analysis/pca_results.eigenvec` (stage 02) and
  `covariate_table_analytic_mothers.txt` (stage 01).
- `--outroot` : holds `<COHORT>/<plat>_clean` filesets for `gsa`, `lpwgs`, `lpwgs_dosage`.
- `--scores-dir` : `<PGS_id>.txt.gz` harmonized PGS-Catalog files (BP panel; BMI/infant reuse it).
- `--sscore-dir` : stage-05 BP `.sscore` files (`<PGS>__<COHORT>__<plat>.sscore`).

## Notes / caveats
- **A1** uses the best-transferring cohort (GAPPS-Bangladesh) as an internal anchor
  because no European reference sample is genotyped here; interpret distances as
  *relative*. Default R² are the Phase-1 mean-BP EUR values; override with `--r2-tsv`.
- **C1** default BMI score is `PGS000027` (Khera 2019, EUR). Swap with `--bmi-pgs`.
- **C2** needs infant genotype filesets (`<COHORT>/<plat>_infant_clean`) **and** a
  mother↔infant ID map (`--mi-link`, two columns: mother_IID, infant_IID). If either
  is absent both scripts print a clear message and exit 0 — no failure.
- KING and Fst run on the merged set; if it is very large, give A3/A1 more memory.

## Quick start
```bash
module load plink/2.00a4.6 R
bash 06_phase2/run_phase2.sh \
  --merged   /path/merged_gsa \
  --analysis-dir /path/analysis_gsa \
  --outroot  /path/cohort_filesets \
  --sscore-dir /path/sscore \
  --scores-dir /path/scores \
  --no-submit          # drop this to actually sbatch
```
