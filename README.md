# Maternal blood pressure in pregnancy: polygenic score portability and one-sample Mendelian randomization

Analysis code, published figures and aggregate source tables for a study of antenatal
blood pressure and pregnancy outcomes in five South Asian and African pregnancy cohorts.

The study asks two questions. First, how well do published polygenic scores for blood
pressure, all developed in adult populations outside these cohorts, transfer to them?
Each of eight scores, four paired systolic and diastolic families, is tested against
measured antenatal blood pressure, cohort by cohort, over a model of maternal age,
genotyping technology and the first five within-cohort genetic principal components, and
its incremental explanatory contribution is reported with bootstrap confidence intervals.
Second, using one prespecified score per cohort and trait as an instrument, what do
one-sample Mendelian randomization estimates say about the effect of higher genetically
predicted antenatal blood pressure on preterm birth, low birth weight, small for
gestational age and birth weight? Cohort estimates are pooled by ancestry group and overall
with an inverse-variance fixed-effect model, with between-cohort heterogeneity summarized
by Cochran's Q test and I², and the sensitivity of the result to how the antenatal
blood-pressure exposure is defined is examined under four definitions.

Within each ancestry group and blood-pressure trait, the instrument was the score family
with the largest mean incremental R² across the contributing cohorts, provided that it
was highest in at least half of those cohorts and had a first-stage F statistic of at
least 10 in every cohort.

The genetic principal components are computed within each cohort jointly from the two
genotyping technologies, low-pass whole-genome sequencing and the Global Screening Array,
with one genotype record per woman (`genetic_pcs/`).

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
genetic_pcs/           the within-cohort joint-platform principal-component analysis, the
                       descriptive cross-cohort analysis and the pairwise FST (see its README)
config.example.R       every path the R steps use; copy to config.R and edit
run_all.R              runs the six steps in order
data/aggregate/        the aggregate values behind each published display
figures/               the published figures, as PDF and PNG
metadata/
  pgs_manifest.tsv     the eight published blood-pressure scores that were evaluated
  software_versions.tsv the versions the reported analysis was run with
LICENSE                MIT
```

A run writes nothing into `figures/` or `data/aggregate/`: the example configurations
send all generated output to `output/`, which is ignored by git.

## Requirements

R, with `data.table`, `metafor` and the base package `splines`, and a build of R with
`cairo_pdf` available, which the figure scripts use. Polygenic scoring and the
principal-component analysis use PLINK 2; `bcftools` is used by `genetic_pcs/` when present.

`metadata/software_versions.tsv` records the versions the reported analysis was run
with, and which step each was used for:

| component | version | where it was used |
|---|---|---|
| R | 4.3.1 Patched (2023-07-19 r84711) | cohort construction, phenotypes, genetic covariates and scores, the joint-platform principal components, polygenic-score portability, cohort Mendelian randomization, exposure sensitivity |
| data.table | 1.17.6 | data handling throughout the steps above |
| splines | distributed with R (base package, no separate version) | natural cubic spline of gestational age, ns(GA_days, 3), for the gestational-age-standardized exposure definitions |
| PLINK | v2.00a4.6LM 64-bit Intel (29 Aug 2023) | polygenic scoring from the harmonized scoring files; variant quality control, LD pruning and genotype export for the principal-component analysis |
| bcftools | 1.18 | reading the re-imputed array VCF for the array dosages |
| R | 4.3.3 (2024-02-29) | pooling of the cohort estimates and the display scripts |
| metafor | 4.4.0 | inverse-variance fixed-effect meta-analysis, rma(yi, sei, method = "FE", test = "z") |

The analysis uses long-stable interfaces of these packages, so a nearby version will
usually behave identically; the table records what was actually used rather than a
lower bound.

## Inputs

All paths come from two configuration files. Copy `config.example.R` to `config.R` and
`genetic_pcs/config_pcs.example.sh` to `genetic_pcs/config_pcs.sh`, and edit them; both
copies are ignored by git, so the paths to controlled data are never committed. Nothing
else in the repository contains a path.

| setting | what it points at |
|---|---|
| `epi_file` | The phenotype extract: one tab-separated file, one row per participant visit. |
| `sscore_dir` | Polygenic-score files, one per score, cohort and genotyping platform. |
| `genotype_qc_dir` | Genotype quality-control tables, per cohort. |
| `pca_dir`, `pc_file` | The output directory of `genetic_pcs/` and the principal-component file it writes. |

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

**Genotypes for the principal components.** `genetic_pcs/config_pcs.example.sh` lists
them: the cleaned low-pass and array filesets of each cohort, the re-imputed array VCF,
and the relative pairs.

## Running the analysis

```sh
cp config.example.R config.R                                   # then edit the paths
cp genetic_pcs/config_pcs.example.sh genetic_pcs/config_pcs.sh # then edit the paths
Rscript run_all.R 01                     # the cleaned sample and its identifier lists
bash genetic_pcs/run_genetic_pcs.sh all  # the joint-platform principal components and their checks
Rscript run_all.R 02 03 04 05 06         # everything that uses them
Rscript run_all.R 04                     # one step, once the earlier ones have been run
```

| script | what it does |
|---|---|
| `analysis/00_functions.R` | Shared definitions: cohort labels, the prespecified score instruments, outcome definitions, covariate coding, and the helper functions. Sourced by the others; fits nothing itself. |
| `analysis/01_construct_cohort_and_phenotypes.R` | Builds the eligible sample, applies the sample-cleaning rules, and derives maternal age and the four outcomes. Writes the cleaned-sample identifier list. |
| `genetic_pcs/run_genetic_pcs.sh` | Computes the within-cohort joint-platform principal components from one genotype record per woman, checks them against prespecified quality gates, and runs the descriptive cross-cohort analysis and the pairwise FST. |
| `analysis/02_prepare_genetic_covariates_and_pgs.R` | Attaches principal components 1–5 and builds the merged polygenic score for each cohort and trait. |
| `analysis/03_pgs_portability.R` | Fits the portability models for every cohort, trait and score of the four families, with maternal age, genotyping technology and principal components 1–5 in both the base and the expanded model, with the bootstrap, and the cross-platform agreement diagnostic. |
| `analysis/04_mr_cohort_models_and_pooling.R` | Fits the cohort-specific one-sample Mendelian randomization models and pools them by ancestry group and overall with the inverse-variance fixed-effect model, with Cochran's Q and I². |
| `analysis/05_exposure_sensitivity.R` | Repeats the Mendelian randomization under four definitions of the antenatal blood-pressure exposure, with the same pooling. |
| `analysis/06_create_tables_and_figures.R` | Builds the figures and the supplementary estimate tables. |

Generated output goes to the directories named in the configurations, all under
`output/` and all ignored by git:

* `output/derived/` — participant-level intermediates. These stay inside the secure
  environment.
* `output/genetic_pcs/` — the principal-component file (participant-level, stays inside
  the secure environment) and its aggregate diagnostics.
* `output/results/` — aggregate result tables, written at full precision: every number
  is the shortest decimal string that reads back as the identical double, so a table
  reproduces the fit it came from exactly.
* `output/figures/` — figures, each accompanied by a table of the values it plots.
* `output/tables/` — the supplementary estimate tables.

## Published figures

`figures/` holds the figures as published, each as PDF and PNG.

| file | display |
|---|---|
| `figure1_portability` | Figure 1, polygenic score portability |
| `figure2_mr_forest` | Figure 2, Mendelian randomization forest plot |
| `figureS1_mr_framework` | Supplementary Figure S1, the assumptions of the design (a diagram; it has no underlying data) |
| `figureS2_participant_flow_by_site` | Supplementary Figure S2, participant flow by study site |
| `figureS3_cross_cohort_pca` | Supplementary Figure S3, descriptive cross-cohort principal components |
| `figureS4_exposure_definitions` | Supplementary Figure S4, exposure-definition sensitivity |

## Aggregate source values

`data/aggregate/` holds the values behind each published display, so a reader can check
a figure or table without access to the controlled data. Estimates are written at full
precision; the published displays round them.

| file | display it underlies |
|---|---|
| `table1_cohort_characteristics.tsv` | Table 1, as printed |
| `figure1_pgs_portability.tsv` | Figure 1: the 40 plotted cells, with the incremental R² behind each printed label |
| `table_S2_panelA_platform_composition.tsv` | Supplementary Table S2, panel A |
| `table_S2_panelB_cross_platform_agreement.tsv` | Supplementary Table S2, panel B |
| `table_S3_pgs_portability_full.tsv` | Supplementary Table S3, as printed |
| `pgs_portability_estimates.tsv` | the 40 portability estimates behind Figure 1 and Supplementary Table S3, at full precision |
| `mr_cohort_estimates.tsv` | the cohort-specific Mendelian randomization results behind Figure 2 and Supplementary Table S4, panel A |
| `mr_pooled_estimates.tsv` | the fixed-effect pooled results behind Figure 2 and Supplementary Table S4, panel B, with Cochran's Q and I² |
| `figure2_plotted_values.tsv` | Figure 2: the 64 plotted estimates, their intervals, the axis of each panel, and the heterogeneity text of the overall rows |
| `table_S4_panelA_cohort_estimates.tsv` | Supplementary Table S4, panel A, as printed |
| `table_S4_panelB_pooled_estimates.tsv` | Supplementary Table S4, panel B, as printed |
| `figureS2_participant_flow_by_site.tsv` | Supplementary Figure S2: each count, by study site and in total (the audited participant-flow counts) |
| `figureS3_cross_cohort_pca_counts.tsv` | Supplementary Figure S3: the number of women plotted per cohort and genotyping technology |
| `figureS3_cross_cohort_pca_axes.tsv` | Supplementary Figure S3: the percentage of variance of each plotted axis, and the numbers of variants and women |
| `figureS4_exposure_definitions.tsv` | Supplementary Figure S4: the pooled estimate under each of the four exposure definitions |
| `joint_pca_summary.tsv` | the within-cohort joint-platform principal components: women, reference women, variants and the largest share of the variance of PC1–PC5 explained by the record's platform |
| `fst_sylhet_matlab.tsv` | the Hudson FST between Sylhet and Matlab quoted in the Discussion, with its jackknife standard error |

Odds ratios and birth-weight differences are per 10 mmHg. In the Mendelian
randomization tables `yi` and `sei` are the Wald-ratio estimate and its delta-method
standard error on the analysis scale (log odds ratio, or grams), which are what the
fixed-effect model is fitted to.

## Polygenic scores

`metadata/pgs_manifest.tsv` lists the eight published blood-pressure scores that were
evaluated: four paired systolic and diastolic families, named for the ancestry of the
samples in which they were developed (European, South Asian, East Asian and diverse
ancestry). For each score it records the trait, score family, PGS Catalog identifier and
URL, development ancestry, publication, method, source sample size and type, the number of
published variants, the role the score had in the analysis, and the percentage of published
variants recovered on each genotyping platform. The file matches Supplementary Table S1.

## Licence

MIT. See `LICENSE`.
