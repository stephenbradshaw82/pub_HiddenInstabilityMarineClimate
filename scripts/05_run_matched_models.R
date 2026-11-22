# 05 | Run the matched 2D control and full 3D particle model
# - Use 20 paired seeds per case, with 250 particles per voxel for BRAN T1
#   and 750 for each CMIP6 SSP2-4.5 period.
# - Share a jittered seed plan, horizontal fields, valid-voxel mask and time step in both modes.
# - Step from the current voxel every 0.5 years for the 10-year baseline horizon.
# - Tally starts, ends and through-flow, then assign source, sink, corridor and movement classes.
# - Estimate vertical-interior class switching and separate it from 2D seed noise.
#
# Model code in method/01_scripts/610_functions_3D_fast.R
#   run_matched_case() prepares the common seed plan and runs both modes.
#   simulate_particle_chunk() updates longitude, latitude and depth in the time loop.
#   classify_voxels() turns the start, end and through-flow tallies into classes.

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
headline <- subset(matched_case_plan(revision_connectivity_case_plan()), phase == "matched")
if (pilot_only) headline <- headline[match(matched_pilot_ids, headline$case_id), ]
if (anyNA(headline$case_id)) stop("Unknown headline pilot case.")

# Each call uses the paper density and the same starting particles in both
# modes. The first mode fixes vertical velocity at zero. The second keeps it.
completed <- vector("list", nrow(headline))
for (index in seq_len(nrow(headline))) {
  completed[index] <- list(run_regional_case(headline[index, ], run_dir, fixture_dir,
                                             recipe_files, pilot_only, run_matched_case))
}
if (pilot_only) {
  message("Section complete: two headline timing seeds. Set pilot_only to FALSE to finish.")
} else {
  publish_ledger(completed, file.path(run_dir, "cases", "matched_ledger.csv"))
  message("Section complete: three headline cases and 60 paired seeds.")
}
