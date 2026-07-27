# Stage 05 — Polygenic scores → PTB

**Generalizes:** `step7_UNIFIED_PRS_AS.R` + `step7_UNIFIED_PRS_plot_AS.R`

**Input:** `merged` + mothers + phenotypes + covariates + PCs + PGS manifest (config.prs.scores).

**Output:** per-score `.sscore`, incremental-R2 results `.rds`, PTB association tables, figures.

**Generalization:** loops over a PGS manifest (trait × ancestry) instead of 4 hardcoded PGS IDs.
**Ancestry strategy:** SAS scores for Bangladesh/Pakistan, AFR for Pemba/Zambia (old pipeline was
SAS/EUR only — AFR scores still need sourcing; see config TODO). Site list / ancestries are
data-driven, and plotting no longer assumes exactly SAS-vs-EUR or 5 fixed sites.
