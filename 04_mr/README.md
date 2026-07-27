# Stage 04 — Mendelian Randomization

**Generalizes:** `step2_subset_external_sumstats.sh`, `step4.1/4.2`, `step5.x`, `step6_run_mr.R`,
and the comprehensive `step_MR_Analysis_pipeline.R`.

**Input:** internal GWAS (outcome) + external BP GWAS (exposure, GRCh38) + LD ref.

**Output:** instruments, harmonized tables, MR results (IVW/Egger/median/mode) + sensitivity
(heterogeneity, pleiotropy, single-SNP, leave-one-out) + plots.

**Generalization:** external GWAS source, build, ancestry, column names, site labels, and all
paths from config. Build-aware harmonization (positional, GRCh38).
