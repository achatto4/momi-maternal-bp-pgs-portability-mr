# Maternal blood pressure in pregnancy: polygenic score portability and one-sample Mendelian randomization

Analysis code, published figures and aggregate source tables for a study of antenatal
blood pressure and pregnancy outcomes in five South Asian and African pregnancy cohorts.

The study asks two questions. First, how well do published polygenic scores for blood
pressure, all developed in adult populations outside these cohorts, transfer to them?
Each score is tested against measured antenatal blood pressure, cohort by cohort, and
its incremental explanatory contribution is reported with bootstrap confidence
intervals. Second, using one prespecified score per cohort and trait as an instrument,
what do one-sample Mendelian randomization estimates say about the effect of higher
genetically predicted antenatal blood pressure on preterm birth, low birth weight,
small for gestational age and birth weight? Cohort estimates are pooled by ancestry
group and overall with a random-effects model, and the sensitivity of the result to how
the antenatal blood-pressure exposure is defined is examined under four definitions.

## Data availability

The individual-level data of the contributing cohorts are controlled and are **not**
distributed here. Access is governed by the data-access arrangements of the
contributing studies. This repository contains the analysis code, which reads the
inputs described below from paths you supply, and the aggregate tables behind the
published displays. No participant identifier, no participant-level record and no
derived participant-level file is committed; the directory those are written to is
ignored by git.

## Repository layout

```
analysis/              the six analysis steps and their shared definitions
  00_functions.R       cohort labels, instruments, outcome and covariate coding, helpers
  01_construct_cohort_and_phenotypes.R
  02_prepare_genetic_covariates_and_pgs.R
  03_pgs_portability.R
  04_mr_cohort_models_and_pooling.R
  05_exposure_sensitivity.R
  06_create_tables_and_figures.R
config.example.R       every path the pipeline uses; copy to config.R and edit
run_all.R              runs the six steps in order
data/aggregate/        the aggregate values behind each published display
figures/               the published figures, as PDF and PNG
metadata/
  pgs_manifest.tsv     the 11 published blood-pressure scores that were evaluated
  software_versions.tsv the versions the reported analysis was run with
LICENSE                MIT
```

A run writes nothing into `figures/` or `data/aggregate/`: the example configuration
sends all generated output to `output/`, which is ignored by git.

## Requirements

R, with `data.table`, `metafor` and the base package `splines`, and a build of R with
`cairo_pdf` available, which the figure scripts use. Polygenic scoring and the
principal-component projection were run with PLINK 2.

`metadata/software_versions.tsv` records the versions the reported analysis was run
with, and which step each was used for:

| component | version | where it was used |
|---|---|---|
| R | 4.3.1 Patched (2023-07-19 r84711) | cohort construction through exposure sensitivity, in the secure environment |
| data.table | 1.17.6 | data handling in those steps |
| splines | base package, distributed with R | `ns(GA_days, 3)` for the gestational-age-standardized exposures |
| PLINK | v2.00a4.6LM 64-bit Intel (29 Aug 2023) | polygenic scoring and the principal-component projection |
| R | 4.3.3 (2024-02-29) | pooling and the display scripts |
| metafor | 4.4.0 | `rma(yi, sei, method = "REML", test = "z")` |

The analysis uses long-stable interfaces of these packages, so a nearby version will
usually behave identically; the table records what was actually used rather than a
lower bound.

## Inputs

All paths come from a single configuration file. Copy `config.example.R` to `config.R`
and edit it; `config.R` is ignored by git, so the paths to controlled data are never
committed. Nothing else in the repository contains a path.

| setting | what it points at |
|---|---|
| `epi_file` | The phenotype extract: one tab-separated file, one row per participant visit. |
| `sscore_dir` | Polygenic-score files, one per score, cohort and genotyping platform. |
| `genotype_qc_dir` | Genotype quality-control tables, per cohort. |
| `pc_file` | Principal components, one row per participant. |

**Phenotype extract.** The pipeline reads these columns and no others:
`SITE_CODE`, `PARTICIPANT_ID`, `PREGNANCY_ID`, `SBP`, `DBP`, `GA_HDLK_NEW`, `VISITDT`,
`PW_AGE`, `SINGLE_TWIN`, `DEL_DATE`, `DATE_LMP`, `GAGEBRTH_NEW`, `BWT_MEASURE_DATE`,
`PTB_NEW`, `BIRTH_WEIGHT`, `SGA_10_NEW`, `BIRTH_OUTCOME`. Missing values may be coded
`-88`, `-77`, `-99`, `NA`, `na`, `.` or empty. `SITE_CODE` identifies the cohort.

**Polygenic-score files.** One file per score, cohort and platform, named

```
<PGS catalogue id>__<cohort>__<platform>.sscore
```

with platform `gsa` or `lpwgs_dosage`, carrying the participant identifier in the first
column and an average score in `SCORE1_AVG`. A woman scored on both platforms
contributes to both files. The scores are those listed in `metadata/pgs_manifest.tsv`.

**Genotype quality-control directory.** For each cohort:

```
work/<cohort>/sample_flags.tsv             ID, lp, gsa, lp_qc_fail, gsa_qc_fail, dual_discordant
work/<cohort>/pairs_identical_linkage.tsv  a, b, linkage, possible_mz
work/<cohort>/pairs_person_level.tsv       a, b, class
work/<cohort>/lpqc.smiss                   per-record call-rate table for the low-pass platform
work/<cohort>/gsaqc.smiss                  per-record call-rate table for the array platform
```

`sample_flags.tsv` says which platforms a woman has a record on, which of those records
failed quality control, and whether her two records are genetically discordant.
`pairs_identical_linkage.tsv` lists pairs of genetically identical records with the
linkage class that was established for each pair; `pairs_person_level.tsv` lists
ordinary relative pairs with their degree.

**Principal components.** Either an `.rds` holding a data frame or a tab-separated file,
with a column `IID` and columns `PC1` to `PC5` or more.

## Running the analysis

```sh
cp config.example.R config.R     # then edit the paths
Rscript run_all.R                # all six steps, in order
Rscript run_all.R 04             # one step, once the earlier ones have been run
```

| script | what it does |
|---|---|
| `analysis/00_functions.R` | Shared definitions: cohort labels, the prespecified score instruments, outcome definitions, covariate coding, and the helper functions. Sourced by the others; fits nothing itself. |
| `analysis/01_construct_cohort_and_phenotypes.R` | Builds the eligible sample, applies the sample-cleaning rules, and derives maternal age and the four outcomes. Writes the cleaned-sample identifier list. |
| `analysis/02_prepare_genetic_covariates_and_pgs.R` | Attaches principal components and builds the merged polygenic score for each cohort and trait. |
| `analysis/03_pgs_portability.R` | Fits the portability models for every cohort, trait and published score, with the bootstrap, and the cross-platform agreement diagnostic. |
| `analysis/04_mr_cohort_models_and_pooling.R` | Fits the cohort-specific one-sample Mendelian randomization models and pools them by ancestry group and overall. |
| `analysis/05_exposure_sensitivity.R` | Repeats the Mendelian randomization under four definitions of the antenatal blood-pressure exposure. |
| `analysis/06_create_tables_and_figures.R` | Builds the forest plot, the supplementary estimate tables and the exposure-definition figure. |

Generated output goes to the four directories named in the configuration, all under
`output/` and all ignored by git:

* `output/derived/` — participant-level intermediates. These stay inside the secure
  environment.
* `output/results/` — aggregate result tables, written at full precision: every number
  is the shortest decimal string that reads back as the identical double, so a table
  reproduces the fit it came from exactly.
* `output/figures/` — figures, each accompanied by a table of the values it plots.
* `output/tables/` — the supplementary estimate tables.

## Published figures

`figures/` holds the figures as published, each as PDF and PNG.

| file | display |
|---|---|
| `figure1_participant_flow` | Figure 1, participant flow |
| `figure2_portability_heatmap` | Figure 2, polygenic score portability |
| `figure3_mr_forest` | Figure 3, Mendelian randomization forest plot |
| `figureS1_exposure_definitions` | Supplementary Figure S1, exposure-definition sensitivity |
| `figureS2_mr_framework` | Supplementary Figure S2, the assumptions of the design (a diagram; it has no underlying data) |

## Aggregate source values

`data/aggregate/` holds the values behind each published display, so a reader can check
a figure or table without access to the controlled data. Estimates are written at full
precision; the published displays round them.

| file | display it underlies |
|---|---|
| `figure1_participant_flow.tsv` | Figure 1: each box and exclusion, and the reasons given in the legend |
| `table1_cohort_characteristics.tsv` | Table 1, as printed |
| `figure2_pgs_portability.tsv` | Figure 2: the 50 plotted cells, with the incremental R² behind each printed label |
| `table_S2_panelA_platform_composition.tsv` | Supplementary Table S2, panel A |
| `table_S2_panelB_cross_platform_agreement.tsv` | Supplementary Table S2, panel B |
| `table_S3_pgs_portability_full.tsv` | Supplementary Table S3, as printed |
| `mr_cohort_estimates.tsv` | the cohort-specific Mendelian randomization results behind Figure 3 and Supplementary Table S4, panel A |
| `mr_pooled_estimates.tsv` | the pooled results behind Figure 3 and Supplementary Table S4, panel B, with the heterogeneity statistics |
| `figure3_plotted_values.tsv` | Figure 3: the 64 plotted estimates, their intervals, and the axis of each panel |
| `table_S4_panelA_cohort_estimates.tsv` | Supplementary Table S4, panel A, as printed |
| `table_S4_panelB_pooled_estimates.tsv` | Supplementary Table S4, panel B, as printed |
| `figure_S1_exposure_definitions.tsv` | Supplementary Figure S1: the pooled estimate under each of the four exposure definitions |

Odds ratios and birth-weight differences are per 10 mmHg. In the Mendelian
randomization tables `yi` and `sei` are the Wald-ratio estimate and its delta-method
standard error on the analysis scale (log odds ratio, or grams), which are what the
random-effects model is fitted to.

## Polygenic scores

`metadata/pgs_manifest.tsv` lists the 11 published blood-pressure scores that were
evaluated: the ten primary scores, five paired systolic and diastolic families, and one
supplementary African-ancestry systolic score. For each score it records the trait,
score family, PGS Catalog identifier and URL, development ancestry, publication, method,
source sample size and type, the number of published variants, the role the score had in
the analysis, and the percentage of published variants recovered on each genotyping
platform. The file matches Supplementary Table S1.

## Licence

MIT. See `LICENSE`.
