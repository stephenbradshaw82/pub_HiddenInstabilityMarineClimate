# 03 | Calculate horizontal and vertical climate velocity
# If you chose "prepared", run this after script 01.
# If you chose "download", run this after script 02 has prepared the NetCDFs.
# Both routes provide the same checked GeoTIFF layout in the run's prepared folder.
# - Read the monthly temperatures and form annual means for T1 and T2.
# - Estimate the temperature trend by linear regression and use local spatial gradients.
# - Save HVoCC in km per year with a compass bearing and signed VVoCC in m per year.
# - Build the sensitivity fields with the specified grid, gradient floor and layer thickness.
# - Average temperatures over each layer before recalculating velocity for aggregation cases.

if (!exists("demo_dir", inherits = FALSE)) stop("Run 00_config.R first.")
source(file.path(demo_dir, "support", "raw_data.R"))
source(file.path(demo_dir, "support", "region.R"))
source(file.path(demo_dir, "support", "velocity.R"))
source(file.path(demo_dir, "support", "checkpoints.R"))
check_packages("terra")
production_file <- file.path(method_dir, "00_src", "revision_velocity_functions.R")
source(production_file)
prepared_dir <- file.path(run_dir, "prepared")
prepared_receipt <- file.path(prepared_dir, "temperature_manifest.csv")
if (!file.exists(prepared_receipt)) {
  stop('Complete script 01 in "prepared" mode, or scripts 01 and 02 in "download" mode, before running 03.')
}
check_file_manifest(prepared_receipt)
velocity_dir <- file.path(run_dir, "velocity")
dir.create(velocity_dir, recursive = TRUE, showWarnings = FALSE)
script_path <- file.path(demo_dir, "scripts", "03_calculate_velocity.R")
check_stage_recipe(file.path(velocity_dir, "recipe.rds"),
                   list(periods = time_periods, depths = analysis_depths_m,
                        md5 = unname(tools::md5sum(c(prepared_receipt, production_file, script_path,
                                                   file.path(demo_dir, "support", "velocity.R"))))))

corrected_field_plan <- function() {
  bran <- expand.grid(grid_deg = c(0.25, 0.50), gradient_floor = c(0, 0.01, 0.05))
  bran$aggregation_m <- NA_real_
  aggregate <- data.frame(grid_deg = 0.25, gradient_floor = c(0, 0.01, 0.05, 0.01),
                          aggregation_m = c(10, 10, 10, 20))
  bran <- rbind(bran, aggregate)
  bran$dataset <- "BRAN"
  bran$period <- "T1"
  cmip <- data.frame(grid_deg = 0.25, gradient_floor = 0.01, aggregation_m = NA_real_,
                     dataset = "CMIP6_ssp245", period = c("T1", "T2"))
  rbind(bran, cmip)
}

velocity_plan <- corrected_field_plan()
velocity_plan$period_id <- unname(time_periods[velocity_plan$period])
velocity_plan$variant <- mapply(velocity_variant_name, velocity_plan$grid_deg,
                                velocity_plan$gradient_floor, velocity_plan$aggregation_m)
velocity_plan$bundle <- file.path(velocity_plan$dataset, velocity_plan$period,
                                  velocity_plan$variant, "velocity_bundle.rds")
velocity_plan$elapsed_s <- NA_real_
velocity_plan$valid_voxels <- NA_integer_
for (index in seq_len(nrow(velocity_plan))) {
  case <- velocity_plan[index, ]
  destination <- dirname(file.path(velocity_dir, case$bundle))
  receipt <- file.path(destination, "files.csv")
  if (!dir.exists(destination)) {
    paths <- file.path(prepared_dir, sprintf("SpatRaster_%s_depth_%05.1fm.tif",
                                             case$dataset, analysis_depths_m))
    for (path in paths) verify_velocity_window(path, case$period_id)
    names(paths) <- as.character(analysis_depths_m)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    staged <- tempfile(".velocity-", tmpdir = dirname(destination))
    # The same annual trends and gradient rules are used for both matched
    # models. True layer aggregation averages temperatures before velocities.
    elapsed <- system.time({
      fields <- build_revision_velocity_bundle(
        input_paths = paths, dataset_id = case$dataset, period_id = case$period_id,
        grid_deg = case$grid_deg, depth_centres_m = analysis_depths_m,
        gradient_floor_quantile = case$gradient_floor, output_dir = staged,
        aggregation_thickness_m = case$aggregation_m
      )
    })[["elapsed"]]
    verify_grid_origin(bundle_grid(fields), case$grid_deg)
    if (length(fields$depths_m) < 3L || !any(fields$valid_mask)) stop("Empty or shallow field: ", case$variant)
    publish_csv(data.frame(elapsed_s = elapsed, valid_voxels = sum(fields$valid_mask)),
                file.path(staged, "timing.csv"))
    files <- list.files(staged)
    publish_csv(data.frame(file = files, md5 = unname(tools::md5sum(file.path(staged, files)))),
                file.path(staged, "files.csv"))
    if (!file.rename(staged, destination)) stop("Could not publish ", destination)
  }
  check_file_manifest(receipt)
  timing <- utils::read.csv(file.path(destination, "timing.csv"))
  velocity_plan$elapsed_s[index] <- timing$elapsed_s
  velocity_plan$valid_voxels[index] <- timing$valid_voxels
  message("Velocity ", index, "/", nrow(velocity_plan), ": ", case$dataset, " ", case$period, " ", case$variant)
}
velocity_plan$production_md5 <- unname(tools::md5sum(production_file))
velocity_plan$temperatures_manifest_md5 <- unname(tools::md5sum(prepared_receipt))
plan_path <- file.path(velocity_dir, "corrected_field_plan.csv")
if (!file.exists(plan_path)) publish_csv(velocity_plan, plan_path)
message("Section complete: ", nrow(velocity_plan), " corrected velocity variants in ", velocity_dir)
