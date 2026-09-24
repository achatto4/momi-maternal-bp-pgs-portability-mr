## -------------------------------------------------------------------------------
## run_all.R -- runs the analysis end to end, in order.
##
##   Rscript run_all.R
##
## Copy config.example.R to config.R and edit the paths first.  Each script can also
## be run on its own, in the same order, once the earlier ones have been run.
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

for (s in STEPS) {
  message("\n=== ", s, " ===")
  t <- system.time(source(file.path(MOMI_ROOT, "analysis", s), local = new.env()))
  message(sprintf("    %.1f s", t[["elapsed"]]))
}
message("\ndone: results in ", config$results_dir, ", figures in ", config$figures_dir,
        ", tables in ", config$tables_dir)
