# Stage 03 — Internal GWAS

**Generalizes:** `step4.0_merged_run_plink_gwas.sh`

**Input:** `merged` + `mothers.keep` + phenotypes + covariate table + PCs.

**Output:** `results_<TRAIT>/gwas_<TRAIT>.PHENO.glm.{linear,logistic.hybrid}`.

**Generalization:** trait list, BP window, and covariate set (incl. N PCs) from config.
PLINK2 binary/module from config.
