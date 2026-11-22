# 02 | Prepare monthly temperatures on the common analysis grid
# Run this script only when script 01 is set to "download" and all raw NetCDFs are ready.
# If script 01 is set to "prepared", skip 02 and run 03.
# - Read the native monthly sea temperatures and check their dates and depth coordinates.
# - Interpolate BRAN and CMIP6 to 42 levels at 5 m intervals and a separate 2.5 m layer.
# - Align the fields to the NOAA OISST 0.25-degree grid. Join CMIP6 historical and SSP2-4.5.
# - At each depth, keep cells with finite values across all input months and datasets.
# - Save the dated GeoTIFFs and intermediate chunks so preparation can restart.

if (!exists("demo_dir", inherits = FALSE)) stop("Run 00_config.R first.")
source(file.path(demo_dir, "support", "raw_data.R"))
source(file.path(demo_dir, "support", "prepared_inputs.R"))

regular_raster <- function(x, y) {
  dx <- diff(x)[1]
  dy <- diff(y)[1]
  if (length(x) < 2L || length(y) < 2L || dx <= 0 || dy <= 0 ||
      any(abs(diff(x) - dx) > 1e-4) || any(abs(diff(y) - dy) > 1e-4)) {
    stop("The regular grid must have increasing, equally spaced coordinates.")
  }
  # Use the mean spacing to avoid preserving NetCDF coordinate round-off.
  dx <- mean(diff(x))
  dy <- mean(diff(y))
  terra::rast(ncols = length(x), nrows = length(y),
              xmin = min(x) - dx / 2, xmax = max(x) + dx / 2,
              ymin = min(y) - dy / 2, ymax = max(y) + dy / 2, crs = "EPSG:4326")
}

profile_at_depths <- function(native_depths, values, target_depths, rule) {
  valid <- is.finite(values)
  if (!any(valid)) return(rep(NA_real_, length(target_depths)))
  if (sum(valid) == 1L) {
    if (rule == 2L) return(rep(values[valid], length(target_depths)))
    return(ifelse(abs(target_depths - native_depths[valid]) < 1e-8, values[valid], NA_real_))
  }
  stats::approx(native_depths[valid], values[valid], xout = target_depths, rule = rule)$y
}

cmip_regrid <- function(plane, lon, lat, lon_jittered, lat_jittered, target_x, target_y) {
  valid <- is.finite(plane) & is.finite(lon) & is.finite(lat)
  if (sum(valid) < 3L) return(matrix(NA_real_, length(target_x), length(target_y)))
  values <- akima::interp(lon_jittered[valid], lat_jittered[valid], plane[valid],
                          xo = target_x, yo = target_y, linear = TRUE, extrap = FALSE)$z
  if (all(valid)) return(values)
  # Interpolate finite-water coverage separately, then mask cells below 50%.
  # akima reports when it resolves collinearity with jitter on a local point cloud.
  coords <- is.finite(lon) & is.finite(lat)
  mask <- withCallingHandlers(
    akima::interp(lon[coords], lat[coords], as.numeric(valid[coords]),
                   xo = target_x, yo = target_y, linear = FALSE, extrap = FALSE)$z,
    warning = function(w) {
      if (grepl("collinear points|collinearities reduced", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  values[!is.finite(mask) | mask < 0.5] <- NA_real_
  values
}

prepare_source_months <- function(path, dataset, target) {
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  axes <- read_raw_axes(nc, dataset)
  variable <- c(BRAN = "temp", CMIP6 = "thetao", NOAA = "sst")[[dataset]]
  time_axis <- tail(axes, 1L)[[1]]
  calendar <- ncdf4::ncatt_get(nc, time_axis$name, "calendar")$value
  months <- read_cf_months(time_axis$vals, time_axis$units, calendar)
  temperature <- ncdf4::ncvar_get(nc, variable, collapse_degen = FALSE)
  if (dataset == "NOAA") {
    raw_grid <- regular_raster(axes[[1]]$vals, axes[[2]]$vals)
    n_depth <- 1L
  } else {
    native_depths <- axes[[3]]$vals
    n_depth <- length(prepared_depths_m)
    if (dataset == "BRAN") {
      raw_grid <- regular_raster(axes[[1]]$vals, axes[[2]]$vals)
      depth_rule <- 1L
    } else {
      lon <- as.vector(ncdf4::ncvar_get(nc, "lon"))
      lat <- as.vector(ncdf4::ncvar_get(nc, "lat"))
      set.seed(42)
      lon_jittered <- lon + runif(length(lon), -1e-6, 1e-6)
      lat_jittered <- lat + runif(length(lat), -1e-6, 1e-6)
      wa_x <- seq(108.75, 130, by = 0.25)
      wa_y <- seq(-39.25, -9, by = 0.25)
      x <- which(wa_x >= terra::xmin(target) - 0.5 & wa_x <= terra::xmax(target) + 0.5)
      y <- which(wa_y >= terra::ymin(target) - 0.5 & wa_y <= terra::ymax(target) + 0.5)
      # Build the full WA coordinate geometry before slicing so cell positions
      # remain consistent with the established 0.25-degree grid.
      x_step <- diff(range(wa_x)) / length(wa_x)
      y_step <- diff(range(wa_y)) / length(wa_y)
      raw_grid <- regular_raster(wa_x[1] + (x - 0.5) * x_step,
                                 wa_y[1] + (y - 0.5) * y_step)
      depth_rule <- 2L
    }
  }
  aligned <- array(NA_real_, c(terra::ncell(target), n_depth, length(months)))
  same_grid <- terra::compareGeom(raw_grid, target, stopOnError = FALSE)
  for (i in seq_along(months)) {
    if (dataset == "NOAA") {
      profiles <- matrix(as.vector(temperature[, rev(seq_len(dim(temperature)[2])), i]), ncol = 1L)
    } else {
      if (dataset == "BRAN") {
        native_layers <- temperature[, , , i]
      } else {
        native_layers <- array(NA_real_, c(length(x), length(y), length(native_depths)))
        for (k in seq_along(native_depths)) {
          native_layers[, , k] <- cmip_regrid(as.vector(temperature[, , k, i]), lon, lat,
                                             lon_jittered, lat_jittered, wa_x[x], wa_y[y])
        }
      }
      # NetCDF y increases northward. Raster rows start at the north edge.
      native_layers <- native_layers[, rev(seq_len(dim(native_layers)[2])), , drop = FALSE]
      raw_profiles <- matrix(native_layers, ncol = length(native_depths))
      profiles <- t(apply(raw_profiles, 1L, function(row) {
        profile_at_depths(native_depths, row, prepared_depths_m, depth_rule)
      }))
    }
    r <- terra::setValues(terra::rast(raw_grid, nlyrs = n_depth), profiles)
    if (!same_grid) {
      # NOAA is already aligned and only needs an exact geographic crop.
      r <- if (dataset == "NOAA") terra::crop(r, target) else
        terra::resample(r, target, method = "bilinear")
    }
    aligned[, , i] <- terra::values(r, mat = TRUE)
  }
  list(values = aligned, months = months)
}

prepare_fixture <- function(fixture_dir, prepared_dir) {
  check_packages(c("ncdf4", "terra", "akima"))
  receipt <- file.path(fixture_dir, "sources.csv")
  if (!file.exists(receipt)) stop("Run 01_load_data.R first. Missing ", receipt)
  sources <- utils::read.csv(receipt)
  if (!nrow(sources) || anyDuplicated(sources$file) || length(unique(sources$pilot)) > 1L) {
    stop("Invalid raw fixture inventory: ", receipt)
  }
  pilot <- if ("pilot" %in% names(sources)) isTRUE(sources$pilot[1]) else
    all(sources$first_month >= "1994-01-01" & sources$last_month <= "1994-12-01")
  inputs <- file.path(fixture_dir, sources$file)
  if (!identical(unname(tools::md5sum(inputs)), sources$md5)) {
    stop("Missing or modified raw temperature file; check the fixture receipts.")
  }

  # Use the NOAA centres as the target, including the cells needed for gradients.
  noaa_path <- inputs[sources$dataset == "NOAA"]
  if (length(noaa_path) != 1L) stop("A fixture needs one NOAA grid source.")
  noaa <- ncdf4::nc_open(noaa_path)
  target <- tryCatch({
    grid <- regular_raster(noaa$dim$lon$vals, noaa$dim$lat$vals)
    terra::crop(grid, terra::ext(unname(analysis_bounds)), snap = "near")
  }, finally = ncdf4::nc_close(noaa))
  if (!isTRUE(all.equal(unname(as.vector(terra::ext(target))), unname(analysis_bounds)))) {
    stop("The observed NOAA grid does not cover the declared preparation bounds.")
  }

  dir.create(prepared_dir, recursive = TRUE, showWarnings = FALSE)
  chunks <- file.path(dirname(prepared_dir), "preparation_chunks")
  dir.create(chunks, showWarnings = FALSE)
  recipe_file <- file.path(chunks, "recipe.rds")
  recipe <- list(raw_md5 = sources$md5, bounds = analysis_bounds, depths = prepared_depths_m,
                 code = unname(tools::md5sum(c(file.path(demo_dir, "scripts", "02_prepare_data.R"),
                                                         file.path(demo_dir, "support", "raw_data.R")))),
                 terra = as.character(utils::packageVersion("terra")),
                 akima = as.character(utils::packageVersion("akima")))
  if (file.exists(recipe_file)) {
    if (!identical(recipe, readRDS(recipe_file))) stop("Preparation inputs or code changed. Choose a new run ID.")
  } else {
    if (length(list.files(prepared_dir))) stop("Choose a fresh prepared directory: ", prepared_dir)
    saveRDS(recipe, recipe_file)
  }

  completed <- file.path(prepared_dir, "temperature_manifest.csv")
  if (file.exists(completed)) {
    files <- utils::read.csv(completed)
    if (nrow(files) != 87L || !identical(unname(tools::md5sum(file.path(prepared_dir, files$file))), files$md5)) {
      stop("Prepared outputs differ from the completion receipt.")
    }
    return(invisible(files))
  }

  common_mask <- matrix(TRUE, terra::ncell(target), length(prepared_depths_m))
  group_files <- list(BRAN = character(), NOAA = character(), CMIP6_ssp245 = character())
  group_months <- lapply(group_files, function(x) as.Date(character()))
  elapsed <- numeric(nrow(sources))
  for (i in seq_len(nrow(sources))) {
    chunk <- file.path(chunks, sprintf("input_%03d.rds", i))
    dataset <- sources$dataset[i]
    if (!file.exists(chunk)) {
      elapsed[i] <- system.time(monthly <- prepare_source_months(inputs[i], dataset, target))[["elapsed"]]
      coverage <- apply(is.finite(monthly$values), c(1, 2), all)
      group <- if (dataset != "CMIP6") dataset else {
        if (sources$experiment[i] %in% c("historical", projection_scenario)) "CMIP6_ssp245" else NA_character_
      }
      selected <- seq_along(monthly$months)
      if (identical(group, "CMIP6_ssp245")) {
        selected <- which(monthly$months >= as.Date("1994-01-01") &
                            monthly$months <= as.Date("2053-12-01"))
      }
      # All retained files contribute to the common mask. Select only the
      # T1 and T2 months when writing the CMIP6 output stack.
      values <- if (!is.na(group) && length(selected)) monthly$values[, , selected, drop = FALSE] else NULL
      entry <- list(group = group, months = monthly$months[selected], coverage = coverage,
                    values = values, elapsed_s = elapsed[i])
      staged <- tempfile(".chunk-", tmpdir = chunks)
      saveRDS(entry, staged, compress = FALSE)
      if (!file.rename(staged, chunk)) stop("Could not publish prepared source chunk: ", chunk)
      message("Prepared ", i, "/", nrow(sources), " ", sources$file[i], " (", round(elapsed[i], 2), " s)")
    } else {
      entry <- readRDS(chunk)
      elapsed[i] <- entry$elapsed_s
    }
    if (dataset == "NOAA") {
      common_mask[, 1] <- common_mask[, 1] & entry$coverage[, 1]
    } else {
      common_mask <- common_mask & entry$coverage
    }
    if (!is.null(entry$values)) {
      group_files[[entry$group]] <- c(group_files[[entry$group]], chunk)
      group_months[[entry$group]] <- c(group_months[[entry$group]], entry$months)
    }
  }
  if (!any(common_mask[, 2:4])) stop("No common finite coverage at the first three analysis depths.")
  for (dataset in names(group_files)) {
    months <- group_months[[dataset]]
    window <- if (pilot) c("1994-01", "1994-12") else switch(
      dataset, BRAN = c("1993-01", "2023-12"), NOAA = c("1981-09", "2025-06"),
      CMIP6_ssp245 = c("1994-01", "2053-12")
    )
    check_months(months, month_sequence(window[1], window[2]), dataset)
    depths <- if (dataset == "NOAA") 2.5 else prepared_depths_m
    # Assemble only one depth at a time, after the final common mask is known.
    for (k in seq_along(depths)) {
      name <- sprintf("SpatRaster_%s_depth_%05.1fm.tif", dataset, depths[k])
      output <- file.path(prepared_dir, name)
      if (file.exists(output)) next
      values <- do.call(cbind, lapply(group_files[[dataset]], function(chunk) {
        matrix(readRDS(chunk)$values[, k, , drop = FALSE], nrow = terra::ncell(target))
      }))
      values[!common_mask[, k], ] <- NA_real_
      stack <- terra::setValues(terra::rast(target, nlyrs = length(months)), values)
      terra::time(stack) <- months
      names(stack) <- format(months, "%Y-%m")
      staged <- tempfile(".temperature-", tmpdir = prepared_dir, fileext = ".tif")
      terra::writeRaster(stack, staged, overwrite = FALSE)
      if (!file.rename(staged, output)) stop("Could not publish ", name)
    }
  }
  utils::write.csv(data.frame(input = sources$file, elapsed_s = elapsed),
                   file.path(prepared_dir, "preparation_times.csv"), row.names = FALSE)
  files <- list.files(prepared_dir, "[.]tif$", full.names = TRUE)
  validate_prepared_temperature(prepared_dir, pilot)
  if (length(files) != 87L) stop("Expected 87 prepared temperature stacks.")
  staged <- tempfile(".temperature-manifest-", tmpdir = prepared_dir)
  utils::write.csv(data.frame(file = basename(files), md5 = unname(tools::md5sum(files))),
                   staged, row.names = FALSE)
  if (!file.rename(staged, completed)) stop("Could not publish the temperature manifest.")
}

if (identical(get0("input_source", ifnotfound = ""), "prepared")) {
  message("Skip script 02 in prepared mode. Continue with script 03.")
} else if (identical(get0("input_source", ifnotfound = ""), "download")) {
  prepare_fixture(fixture_dir, file.path(run_dir, "prepared"))
  message("Section complete: 87 checked temperature stacks in ", file.path(run_dir, "prepared"))
} else {
  stop("Run script 01 and select the input source before script 02.")
}
