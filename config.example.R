## -------------------------------------------------------------------------------
## config.example.R
##
## Copy this file to config.R and edit the paths for your installation.  Every path
## used anywhere in the R pipeline comes from here; no script contains a path of its own.
## config.R is listed in .gitignore, so a local copy is never committed.  The
## principal-component analysis in genetic_pcs/ has its own configuration file,
## genetic_pcs/config_pcs.example.sh.
## -------------------------------------------------------------------------------

config <- list(

  ## ---- inputs ----------------------------------------------------------------
  ## Phenotype extract: one tab-separated file, one row per participant visit.  The
  ## columns the pipeline reads are listed in the README.
  epi_file = "/path/to/phenotype_extract.txt",

  ## Directory of polygenic-score files, one per score, cohort and genotyping
  ## platform, named  <PGS id>__<cohort>__<platform>.sscore  and carrying the
  ## columns the scoring tool writes (the identifier in the first column and
  ## SCORE1_AVG).  Platforms are "gsa" and "lpwgs_dosage".
  sscore_dir = "/path/to/pgs_scores",

  ## Genotype quality-control directory.  It must contain, for each cohort,
  ##   work/<cohort>/sample_flags.tsv            ID, lp, gsa, lp_qc_fail,
  ##                                             gsa_qc_fail, dual_discordant
  ##   work/<cohort>/pairs_identical_linkage.tsv a, b, linkage, possible_mz
  ##   work/<cohort>/pairs_person_level.tsv      a, b, class
  ##   work/<cohort>/lpqc.smiss, gsaqc.smiss     per-record missingness
  genotype_qc_dir = "/path/to/genotype_qc",

  ## Output directory of genetic_pcs/ (its JPCA_OUT), and the principal-component file
  ## it writes there: one row per participant, columns IID and PC1..PC5 or more.
  pca_dir = "output/genetic_pcs",
  pc_file = "output/genetic_pcs/participant_level/joint_pcs.rds",

  ## ---- outputs ---------------------------------------------------------------
  ## derived_dir holds participant-level intermediates and never leaves the secure
  ## environment.  The other three hold aggregate results only.  All four sit under
  ## output/, which is ignored by git, so a run never writes over the published
  ## figures in figures/ or the aggregate source tables in data/aggregate/.
  derived_dir = "output/derived",
  results_dir = "output/results",
  figures_dir = "output/figures",
  tables_dir  = "output/tables",

  ## ---- analysis settings -----------------------------------------------------
  cohorts = c("AMANHI-Bangladesh", "AMANHI-Pakistan", "GAPPS-Bangladesh",
              "AMANHI-Pemba", "GAPPS-Zambia"),

  ## Bootstrap for the portability confidence intervals.  The seed rule is
  ## base + 10 * cohort index + trait index, so a run is reproducible.
  bootstrap_replicates = 1000,
  bootstrap_seed       = 20260921
)
