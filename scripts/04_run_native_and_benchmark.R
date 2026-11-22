# 04 | Classify native VoCC trajectories and the velocity benchmark
# - Use the VoCC R package for native horizontal velocity and two-dimensional trajectories.
# - Summarise native trajectory starts, ends and through-flow at each depth.
# - Apply the direct HVoCC and absolute VVoCC thresholds of 20 km/year and 1 m/year.
# - Save both class stacks and the BRAN surface speed and bearing fields.
#
# The version pinned in script 00 supplies tempTrend(), spatGrad(), gVoCC()
# and voccTraj() for the native method.

if (!exists("demo_dir", inherits = FALSE)) stop("Run 00_config.R first.")
for (file in c("raw_data.R", "velocity.R", "checkpoints.R")) {
  source(file.path(demo_dir, "support", file), local = TRUE)
}

native_horizontal_velocity <- function(annual) {
  layers <- raster::brick(annual)
  names(layers) <- format(terra::time(annual), "%Y")
  trend <- VoCC::tempTrend(layers, th = raster::nlayers(layers))
  gradient <- VoCC::spatGrad(layers, th = 0.001, projected = FALSE)
  velocity <- VoCC::gVoCC(trend, gradient)
  list(magnitude = abs(velocity[[1]]), bearing = velocity[[2]],
       mean = raster::calc(layers, mean, na.rm = TRUE))
}

native_trajectory_classes <- function(horizontal, derive_starts, classify_trajectories) {
  vel <- horizontal$magnitude
  ang <- horizontal$bearing
  mn <- horizontal$mean
  valid <- is.finite(raster::values(vel)) & is.finite(raster::values(ang)) &
    is.finite(raster::values(mn))
  valid_raster <- raster::setValues(raster::raster(ang), ifelse(valid, 1, NA_real_))
  cells <- which(valid)
  if (!length(cells)) {
    return(list(classes = terra::rast(valid_raster), starts = 0L, paths = NULL,
                outside_grid_records = 0L))
  }

  starts <- as.data.frame(raster::xyFromCell(vel, cells))
  starts$vel <- raster::values(vel)[cells]
  starts$ang <- raster::values(ang)[cells]
  starts$mn <- raster::values(mn)[cells]
  trajectories <- VoCC::voccTraj(starts, vel, ang, mn, tyr = 30,
                                 trajID = seq_len(nrow(starts)), correct = TRUE)
  counts <- derive_starts(trajectories, ang, valid_raster)
  xy <- as.matrix(trajectories[, 1:2])
  if (!all(is.finite(xy))) stop("The native trajectories contain non-finite coordinates.")
  outside <- is.na(raster::cellFromXY(vel, xy))
  classifier_inputs <- raster::stack(vel, ang, mn)
  if (any(outside)) {
    # Assign finite cell IDs to off-grid records. The added cells have no
    # velocity or temperature and are cropped away after classification.
    spacing <- raster::res(vel)
    classifier_extent <- raster::extent(
      min(raster::xmin(vel), min(xy[outside, 1]) - spacing[1] / 2),
      max(raster::xmax(vel), max(xy[outside, 1]) + spacing[1] / 2),
      min(raster::ymin(vel), min(xy[outside, 2]) - spacing[2] / 2),
      max(raster::ymax(vel), max(xy[outside, 2]) + spacing[2] / 2)
    )
    classifier_inputs <- raster::extend(classifier_inputs, classifier_extent)
    if (anyNA(raster::cellFromXY(classifier_inputs, xy))) {
      stop("Could not assign all trajectory records to classifier cells.")
    }
  }
  classes <- withCallingHandlers(
    classify_trajectories(trajectories, classifier_inputs[[1]], classifier_inputs[[2]],
                          classifier_inputs[[3]], trajSt = counts$trajSt,
                          tyr = 10, nmL = 5, smL = 20, Nend = 10, Nst = 10, NFT = 10),
    warning = function(w) {
      if (grepl("number of items to replace|replacement length", conditionMessage(w))) {
        stop("Native VoCC could not assign trajectory counts safely: ", conditionMessage(w))
      }
    }
  )
  if (any(outside)) classes <- raster::crop(classes, vel)
  list(classes = terra::rast(classes[["TrajClas"]]), starts = counts$n_trajectories,
       paths = trajectories, outside_grid_records = sum(outside))
}

benchmark_classes <- function(horizontal, vertical) {
  # Benchmark class IDs are 0 stable, 1 upward, 2 downward, 3 horizontal,
  # 4 horizontal/upward and 5 horizontal/downward.
  high_h <- horizontal >= 20
  upward <- vertical >= 1
  downward <- vertical <= -1
  terra::ifel(high_h & downward, 5,
    terra::ifel(high_h & upward, 4,
      terra::ifel(high_h, 3, terra::ifel(downward, 2, terra::ifel(upward, 1, 0)))))
}

check_vocc_version <- function() {
  installed_sha <- utils::packageDescription("VoCC", fields = "RemoteSha")
  if (!identical(installed_sha, vocc_sha)) stop("Install the pinned VoCC revision: ", vocc_sha)
}

check_packages(c("terra", "raster", "VoCC", "data.table", "geosphere"))
check_vocc_version()
# The pinned VoCC revision looks up these functions on the search path.
suppressPackageStartupMessages({
  library(raster)
  library(data.table)
  library(geosphere)
})
production_files <- file.path(method_dir, "00_src", c("functions.R", "revision_velocity_functions.R"))
for (path in production_files) source(path, local = TRUE)
prepared_dir <- file.path(run_dir, "prepared")
prepared_receipt <- file.path(prepared_dir, "temperature_manifest.csv")
check_file_manifest(prepared_receipt)
region_path <- file.path(fixture_dir, "regions.gpkg")
simulation_region <- terra::vect(region_path, layer = "simulation")
output_dir <- file.path(run_dir, "comparators")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
check_stage_recipe(file.path(output_dir, "recipe.rds"),
                   list(periods = time_periods, depths = prepared_depths_m, vocc = vocc_sha,
                        md5 = unname(tools::md5sum(c(prepared_receipt, region_path, production_files,
                          file.path(demo_dir, "scripts", "04_run_native_and_benchmark.R"))))))

comparator_plan <- data.frame(dataset = c("BRAN", "CMIP6_ssp245", "CMIP6_ssp245"),
                              period = c("T1", "T1", "T2"))
comparator_plan$elapsed_s <- NA_real_
for (index in seq_len(nrow(comparator_plan))) {
  dataset <- comparator_plan$dataset[index]
  period_id <- comparator_plan$period[index]
  period <- time_periods[[period_id]]
  destination <- file.path(output_dir, dataset, period_id)
  if (!dir.exists(destination)) {
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    staged <- tempfile(".comparator-", tmpdir = dirname(destination))
    dir.create(staged)
    started <- proc.time()[["elapsed"]]
    annual <- lapply(prepared_depths_m, function(depth) {
      path <- file.path(prepared_dir, sprintf("SpatRaster_%s_depth_%05.1fm.tif", dataset, depth))
      verify_velocity_window(path, period)
      revision_annual_temperature(terra::rast(path), period)
    })
    temperature_mean <- lapply(annual, revision_temperature_mean)
    template <- temperature_mean[[1]]
    means <- lapply(temperature_mean, revision_raster_matrix)
    gradient <- revision_vertical_gradient(means, prepared_depths_m)
    calculation_mask <- terra::rasterize(simulation_region, template, field = 1, touches = FALSE)
    native <- benchmark <- vector("list", length(prepared_depths_m))
    starts <- outside <- integer(length(prepared_depths_m))
    for (k in seq_along(prepared_depths_m)) {
      # Estimate velocity before applying the buffer, so gradient neighbours
      # are available. The native model then follows horizontal trajectories.
      horizontal <- native_horizontal_velocity(annual[[k]])
      horizontal$magnitude <- raster::mask(horizontal$magnitude, raster::raster(calculation_mask))
      horizontal$bearing <- raster::mask(horizontal$bearing, raster::raster(calculation_mask))
      horizontal$mean <- raster::mask(horizontal$mean, raster::raster(calculation_mask))
      classified <- native_trajectory_classes(horizontal, derive_trajSt_from_trajectory_starts,
                                                reference_trajClas)
      native[[k]] <- classified$classes
      starts[k] <- classified$starts
      outside[k] <- classified$outside_grid_records

      # The benchmark assigns classes from speed and direction thresholds.
      # It does not follow trajectories or count where they start and end.
      dz <- revision_matrix_raster(gradient[[k]], template, "vertical_gradient")
      trend <- revision_ols_trend(annual[[k]])
      vertical <- terra::clamp(trend / (abs(dz) + 1e-6) * terra::ifel(dz < 0, -1, 1), -20, 20)
      benchmark[[k]] <- benchmark_classes(terra::rast(horizontal$magnitude),
                                          terra::mask(vertical, calculation_mask))
      if (dataset == "BRAN" && k == 1L) {
        terra::writeRaster(terra::rast(horizontal$magnitude), file.path(staged, "surface_speed.tif"))
        terra::writeRaster(terra::rast(horizontal$bearing), file.path(staged, "surface_bearing.tif"))
      }
    }
    for (view in c("native", "benchmark")) {
      layers <- if (view == "native") native else benchmark
      stack <- terra::rast(layers)
      names(stack) <- paste0("depth_", prepared_depths_m)
      terra::writeRaster(stack, file.path(staged, paste0(view, "_class_ids.tif")))
    }
    publish_csv(data.frame(depth_m = prepared_depths_m, starts = starts,
                           outside_grid_records = outside), file.path(staged, "native_starts.csv"))
    publish_csv(data.frame(elapsed_s = proc.time()[["elapsed"]] - started), file.path(staged, "timing.csv"))
    files <- list.files(staged)
    publish_csv(data.frame(file = files, md5 = unname(tools::md5sum(file.path(staged, files)))),
                file.path(staged, "files.csv"))
    if (!file.rename(staged, destination)) stop("Could not publish ", destination)
  }
  check_file_manifest(file.path(destination, "files.csv"))
  comparator_plan$elapsed_s[index] <- utils::read.csv(file.path(destination, "timing.csv"))$elapsed_s
  message(dataset, " ", period_id, " comparator fields (", comparator_plan$elapsed_s[index], " s)")
}
comparator_plan$vocc_sha <- vocc_sha
receipt <- file.path(output_dir, "comparator_manifest.csv")
if (!file.exists(receipt)) publish_csv(comparator_plan, receipt)
message("Section complete: native trajectories and direct threshold benchmark.")
