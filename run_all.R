## -------------------------------------------------------------------------------
## run_all.R -- runs the analysis end to end, in order.
##
##   Rscript run_all.R 01                             the cleaned sample
##   bash genetic_pcs/run_genetic_pcs.sh all          the joint-platform principal components
##   Rscript run_all.R 02 03 04 05 06                 everything that uses them
##
## Copy config.example.R to config.R and genetic_pcs/config_pcs.example.sh to
## genetic_pcs/config_pcs.sh and edit the paths first.  Step 01 writes the cleaned-sample
## identifier lists that genetic_pcs/ reads; steps 02 to 05 read the principal components
## it writes.  Each script can also be run on its own, in this order.
## -------------------------------------------------------------------------------
MOMI_ROOT <- getwd()
cfg <- file.path(MOMI_ROOT, if (file.exists(file.path(MOMI_ROOT, "config.R")))
                              "config.R" else "config.example.R")
message("configuration: ", basename(cfg))
source(cfg)

STEPS <- c("01_construct_cohort_and_phenotypes.R",
           "02_prepare_genetic_covariates_and_pgs.R",
           "03_pgs_portability.R",
           "04_mr_cohort_models_and_pooling.R",
           "05_exposure_sensitivity.R",
           "06_create_tables_and_figures.R")

only <- commandArgs(TRUE)
if (length(only)) STEPS <- STEPS[substr(STEPS, 1, 2) %in% only]
if (any(substr(STEPS, 1, 2) %in% c("02", "03", "04", "05")) && !file.exists(config$pc_file))
  stop("the principal-component file is missing (", config$pc_file, "): run genetic_pcs/run_genetic_pcs.sh after step 01")

for (s in STEPS) {
  message("\n=== ", s, " ===")
  t <- system.time(source(file.path(MOMI_ROOT, "analysis", s), local = new.env()))
  message(sprintf("    %.1f s", t[["elapsed"]]))
}
message("\ndone: results in ", config$results_dir, ", figures in ", config$figures_dir,
        ", tables in ", config$tables_dir)
