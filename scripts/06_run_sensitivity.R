# 06 | Run the BRAN sensitivity analysis
# - Run the 37 configurations with five paired seeds per case.
# - Vary grid, depth treatment, horizon, gradient floor, boundary policy and tally thresholds.
# - Reuse the matched particle engine and compare classes within the reporting polygon.
# - Estimate attributable switching, its within-case Monte Carlo interval and the seed-noise flag.
#
# Particle stepping and classification follow the same functions as script 05.
# The case plan changes the specified settings. Only the reporting summary clips to Gascoyne.

if (!exists("demo_dir", inherits = FALSE)) stop("Run 00_config.R first.")
for (file in c("raw_data.R", "region.R", "velocity.R", "checkpoints.R", "matched_cases.R")) {
  source(file.path(demo_dir, "support", file))
}
check_packages(c("terra", "geosphere"))
engine_files <- file.path(method_dir, c("00_src/revision_config.R", "01_scripts/610_functions_3D_fast.R"))
for (file in engine_files) source(file)
recipe_files <- c(engine_files, file.path(demo_dir, "support", c("region.R", "velocity.R", "matched_cases.R")))
if (!file.exists(file.path(run_dir, "velocity", "corrected_field_plan.csv"))) {
  stop("Run 03_calculate_velocity.R first.")
}
sensitivity_cases <- subset(matched_case_plan(revision_connectivity_case_plan()), phase == "sensitivity")
if (pilot_only) sensitivity_cases <- sensitivity_cases[match(sensitivity_pilot_ids, sensitivity_cases$case_id), ]
if (anyNA(sensitivity_cases$case_id)) stop("Unknown sensitivity pilot case.")

# Centre sampling retains selected levels. The aggregation cases use fields
# built from layer-averaged temperatures. Both use five matched pairs and
# compare two independent pairs of 2D runs to estimate seed noise.
completed <- vector("list", nrow(sensitivity_cases))
for (index in seq_len(nrow(sensitivity_cases))) {
  completed[index] <- list(run_regional_case(sensitivity_cases[index, ], run_dir, fixture_dir,
                                             recipe_files, pilot_only, run_matched_case))
}
if (pilot_only) {
  message("Section complete: three sensitivity timing seeds. Set pilot_only to FALSE to finish.")
} else {
  publish_ledger(completed, file.path(run_dir, "cases", "sensitivity_ledger.csv"))
  message("Section complete: 37 sensitivity cases and 185 paired seeds.")
}
