# Stage 02 — PCA

**Generalizes:** `step1.3_compute_and_plot_pca.sh`, `step1.4_compute_and_plot_pca.R`

**Input:** `merged.{bed,bim,fam}` + `mothers.keep` + covariate table (for site labels).

**Output:** `pca_analysis/pca_results.{eigenvec,eigenval}` + plots colored by site.

**Generalization:** site labels/colors come from `config.cohorts` (no hardcoded 5-site map,
no defaulting missing sites to "Zambia"). PCs are exported and **used as GWAS covariates**
(fixes the old gap where PCs were computed but ignored).
