# Portability of blood-pressure polygenic scores and Mendelian randomization of maternal blood pressure on perinatal outcomes

Analysis code for a study of (i) the cross-ancestry portability of blood-pressure
polygenic scores across five South Asian and sub-Saharan African pregnancy cohorts, and
(ii) one-sample Mendelian randomization (MR) of maternal blood pressure on perinatal
outcomes (preterm birth, low birth weight, small-for-gestational-age, and birth weight).

Cohorts are drawn from the Maternal and Offspring Mid-life Investigation (MOMI) Consortium:
AMANHI-Sylhet (Bangladesh), AMANHI-Karachi (Pakistan), GAPPS/PreSSMat-Matlab (Bangladesh),
AMANHI-Pemba (Tanzania), and ZAPPS-Lusaka (Zambia).

This repository contains the **analysis code only**. Individual-level genotype and phenotype
data are not included and cannot be redistributed here (see *Data availability* below).

## What the pipeline does

1. **Phenotype construction** — harmonizes antenatal blood-pressure readings and perinatal
   outcomes across cohorts, builds the analytic samples, and derives the gestational-age–
   standardized blood-pressure exposure.
2. **Polygenic scoring** — applies ancestry-matched published blood-pressure polygenic scores
   and evaluates transferability as the incremental R² of the score for measured blood pressure
   in each cohort.
3. **Mendelian randomization** — one-sample MR (reduced-form test and Wald-ratio magnitude) of
   maternal blood pressure on each perinatal outcome, with observational and external
   two-sample estimates for triangulation.
4. **Reporting** — regenerates all manuscript tables, figures, and the supplement directly from
   the frozen result tables, so every number in the paper is reproducible from code.

## Repository layout

```
lib/            shared R libraries (I/O, config, estimators, plotting theme)
00_preprocess/  genotype/dosage preprocessing
01_phenotypes/  analytic-sample and phenotype construction
02_pca/         principal-component analysis
03_gwas/        association utilities
04_mr/          Mendelian-randomization routines
05_prs/         polygenic scoring, transferability, MR deliverables
06_phase2/      within-cohort PCA, pooled PCA, distance, PC-adjusted analyses
report/         Python/R generators for the paper, supplement, slides, and figures
ref/            public reference metadata (literature MR estimates, citations, PGS-catalog metadata)
config/         example configuration
env/, bin/      cluster environment setup and helper scripts
run_*.sh        stage launchers; run_build.sh drives the reporting build (steps B01–B30)
RUNBOOK.md      step-by-step run instructions
```

## Requirements

- R (≥ 4.3) with `data.table` and standard modeling packages
- Python 3 (matplotlib; `graphviz` for the flow diagram and DAG)
- PLINK 2.00 and PLINK 1.90 (the latter for sample-wise merges)

The pipeline was developed to run on an HPC cluster (SLURM). Paths in the launcher scripts are
cluster-specific and will need to be adapted to your environment.

## Running

The reporting build is registry-driven. With the analytic intermediates in place:

```bash
bash run_build.sh --only B08 --results results/momi_output   # one step
bash run_build.sh --results results/momi_output              # full build
```

See `RUNBOOK.md` for the full stage sequence and dependencies.

## Data availability

Individual-level data from the MOMI cohorts are held under the cohorts' own governance and are
not distributed with this code. Access is subject to the data-sharing policies of the MOMI
Consortium and the contributing studies (AMANHI, GAPPS/PreSSMat, ZAPPS). The published
blood-pressure polygenic scores used as instruments are available from the PGS Catalog; the
relevant score metadata are recorded in `ref/`.

## Citation

If you use this code, please cite the accompanying manuscript (in preparation) and this
repository. Citation details will be added on publication.

## License

Released under the MIT License — see [LICENSE](LICENSE).
