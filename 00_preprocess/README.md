# Stage 00 — Preprocessing (VCF → QC'd merged PLINK)

**Generalizes:** `new_data/{vcf2plink,dataprep,merge_data_LPS,mergedata,merge_all_sites_lps_gsa}.sh`

**Input:** per-cohort `gsa_reimputed.vcf.gz` + `lpwgs_imputed.vcf.gz` (GRCh38) and
`*_idmap.csv`, for the 5 cohorts in `config.cohorts`.

**Steps (planned):**
1. Per cohort × platform: `plink2 --vcf ... dosage=<field> --hard-call-threshold` → bed.
2. Variant-quality QC: imputation INFO/R2, MAF, geno, HWE; autosomes only.
3. Mother-filter: keep `-M` samples; link to Epi `PARTICIPANT_ID` via idmap.
4. Standardize variant IDs to `chr:pos`; dedup by position.
5. Merge all cohorts × platforms (non-overlapping samples → stack on common SNP set)
   with flip/exclude retry on allele mismatches.
6. Sample-missingness QC (`mind`).

**Output (the downstream contract):**
`<work_dir>/merged.{bed,bim,fam}` + `mothers.keep` + a QC report.

**Notes:** non-overlapping samples across platforms (different people) — so merging is a
SNP-intersection + sample-stack, not a per-person reconciliation.
