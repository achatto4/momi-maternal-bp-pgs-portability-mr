# config_pcs.example.sh -- copy to config_pcs.sh (ignored by git) and edit the paths. Every path the genetic_pcs/
# scripts use comes from here, and each value can also be overridden by exporting the variable before a run.
: "${JPCA_FILESETS:=/path/to/genotype_filesets}"        # <cohort>/lpwgs_dosage_clean.{pgen,pvar,psam}; <cohort>/gsa_clean fileset and gsa_map.rename
: "${JPCA_ARRAY_VCF_DIR:=/path/to/array_imputation}"     # <cohort>/gsa_reimputed.vcf.gz, the re-imputed array VCF (dosage encoding)
: "${JPCA_RELATEDNESS_DIR:=/path/to/relatedness}"        # <cohort>/pairs_person_level.tsv and <cohort>/xking_dual_self.tsv
: "${JPCA_EPI:=/path/to/phenotype_extract.txt}"          # the phenotype extract of config.R (for the genotyping-technology variable)
: "${JPCA_SSC:=/path/to/pgs_scores}"                     # the score files of config.R (for the genotyping-technology variable)
: "${JPCA_KEEP:=output/derived/cleaned_sample_ids.txt}"               # written by analysis/01
: "${JPCA_DROPS:=output/derived/dropped_genotype_records_ids.tsv}"    # written by analysis/01
: "${JPCA_OUT:=output/genetic_pcs}"                      # everything genetic_pcs/ writes goes below here (ignored by git)
: "${JPCA_PREVIOUS_PCS:=}"                               # optional: an earlier principal-component file to compare with
: "${JPCA_ENCODING:=dosage}"                             # dosage (used for the reported analysis) or hardcall
: "${JPCA_PLINK2:=plink2}"
: "${JPCA_BCFTOOLS:=bcftools}"                           # optional; without it the VCF is read by plink2 --vcf
: "${JPCA_THREADS:=8}"
: "${JPCA_MEMMB:=32000}"
: "${JPCA_BLOCK:=5000}"                                  # variants per exported block
: "${JPCA_CROSS_MAXSNP:=50000}"                          # cap on the cross-cohort variant set (seeded thinning)
: "${JPCA_SEED:=20260925}"
JPCA_COHORTS=(AMANHI-Bangladesh AMANHI-Pakistan GAPPS-Bangladesh AMANHI-Pemba GAPPS-Zambia)
