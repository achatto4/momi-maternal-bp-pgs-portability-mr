# Stage 01 — Analytic mothers, covariates, phenotypes

**Generalizes:** `step1.0_merged_filter_mothers_AS.sh`, `step1.2_covasriates_AS.R`,
`step3_merged_prep_sample_phenotypes_ext_AS.sh`

**Input:** `merged.fam` (stage 00) + Epi pheno file.

**Output:** `analytic_mothers.txt`, `mothers.keep`, `covariate_table_analytic_mothers.txt`,
per-window `mothers_{SBP,DBP}_<window>_final.pheno`, `mothers_PTB_NEW_final.pheno`.

**Generalization:** site-ID handling is config-driven (idmap + site_code), replacing the
hardcoded ZAPPS/`17-|Z2-|30-`/`AMANHIB-` regexes. Canonical PTB coding (2=case).
