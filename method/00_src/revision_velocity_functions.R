# Revision velocity-bundle helpers
#
# These are plain functions for the corrected connectivity analysis. They do
# not read global variables, install packages, download data, or write outside
# the explicitly supplied output directory.

#' Require the packages used by velocity-bundle helpers
#'
#' @returns Invisibly `TRUE`. Stops with installation instructions when `terra`
#'   is unavailable; it never attempts package installation.
revision_require_velocity_packages <- function() {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Missing required package `terra`. Install it outside the analysis scripts, then rerun.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Parse a depth centre from a processed temperature filename
#'
#' @param path Character file path containing `_depth_XXX.Xm`.
#' @returns Numeric depth centre in metres, or `NA_real_` when the filename
#'   does not use the required convention.
revision_depth_from_temperature_path <- function(path) {
  hit <- regexec("_depth_([0-9]+(?:\\.[0-9]+)?)m", basename(path))
  value <- regmatches(basename(path), hit)[[1]]
  if (length(value) < 2L) NA_real_ else as.numeric(value[2])
}

#' Select an exact prepared set of temperature-depth files
#'
#' @param input_paths Named character vector of existing local temperature
#'   paths. Names must be numeric depth centres in metres, or each filename
#'   must contain the `_depth_XXX.Xm` convention.
#' @param depth_centres_m Increasing numeric source-depth centres, in metres,
#'   that are authorised for the requested analysis.
#' @returns Named character vector ordered exactly as `depth_centres_m`. Files
#'   at other available depths are excluded rather than silently included.
#'
#' @details This small boundary helper is used where a directory can contain
#' more processed depths than the configured prepared common grid. It selects
#' only the explicitly requested source files; it does not aggregate, read, or
#' modify a raster.
#'
#' @preconditions Every path exists locally and every requested depth has one
#' unique matching input file.
#'
#' @failure_modes Stops for non-character/missing paths, unparseable or
#' duplicate depth labels, malformed requested depths, or missing requested
#' coverage. It writes no files and makes no network calls.
revision_select_prepared_temperature_paths <- function(input_paths, depth_centres_m) {
  if (!is.character(input_paths) || !length(input_paths) || any(!file.exists(input_paths))) {
    stop("input_paths must contain existing local temperature paths.", call. = FALSE)
  }
  requested <- as.numeric(depth_centres_m)
  if (!length(requested) || any(!is.finite(requested)) || any(diff(requested) <= 0)) {
    stop("depth_centres_m must contain increasing finite numeric metres.", call. = FALSE)
  }
  source_depths <- suppressWarnings(as.numeric(names(input_paths)))
  if (is.null(names(input_paths)) || anyNA(source_depths)) {
    source_depths <- vapply(input_paths, revision_depth_from_temperature_path, numeric(1))
  }
  if (anyNA(source_depths) || anyDuplicated(source_depths)) {
    stop("input_paths must resolve to one unique numeric depth centre per file.", call. = FALSE)
  }
  index <- match(requested, source_depths)
  if (anyNA(index)) {
    stop("Temperature inputs are missing prepared depth centres: ",
         paste(requested[is.na(index)], collapse = ", "), call. = FALSE)
  }
  out <- input_paths[index]
  names(out) <- as.character(requested)
  out
}

#' Make a regular target raster for a velocity bundle
#'
#' @param template Single-layer `terra::SpatRaster` defining the source extent
#'   and coordinate reference system.
#' @param grid_deg Positive target longitude/latitude resolution in degrees.
#' @returns Empty `SpatRaster` at `grid_deg`. Stops when the source is not in a
#'   longitude/latitude coordinate reference system or `grid_deg` is invalid.
revision_target_grid <- function(template, grid_deg) {
  if (!is.numeric(grid_deg) || length(grid_deg) != 1L || !is.finite(grid_deg) || grid_deg <= 0) {
    stop("grid_deg must be one positive finite number of degrees.", call. = FALSE)
  }
  if (!terra::is.lonlat(template)) {
    stop("Revision velocity inputs must be longitude/latitude rasters.", call. = FALSE)
  }
  terra::rast(ext = terra::ext(template), resolution = grid_deg, crs = terra::crs(template))
}

#' Convert monthly or irregular temperature layers to annual means
#'
#' @param temperature Time-indexed `terra::SpatRaster` of temperature.
#' @param period_id Inclusive `YYYYMMDD-YYYYMMDD` analysis period.
#' @param target_grid Optional target `SpatRaster`; every layer is bilinearly
#'   resampled when its geometry differs.
#' @returns Time-indexed annual-mean `SpatRaster` with dates set to 1 July.
#'   Stops for absent time metadata, fewer than three annual means, empty
#'   periods, or incompatible geometry.
revision_annual_temperature <- function(temperature, period_id, target_grid = NULL) {
  if (!inherits(temperature, "SpatRaster")) {
    stop("temperature must be a terra::SpatRaster.", call. = FALSE)
  }
  parts <- strsplit(period_id, "-", fixed = TRUE)[[1]]
  if (length(parts) != 2L || any(!grepl("^[0-9]{8}$", parts))) {
    stop("period_id must use YYYYMMDD-YYYYMMDD.", call. = FALSE)
  }
  dates <- terra::time(temperature)
  if (is.null(dates) || length(dates) != terra::nlyr(temperature) || anyNA(dates)) {
    stop("Temperature raster must have one non-missing date per layer.", call. = FALSE)
  }
  from <- as.Date(parts[1], format = "%Y%m%d")
  to <- as.Date(parts[2], format = "%Y%m%d")
  keep <- dates >= from & dates <= to
  if (!any(keep)) {
    stop("No temperature layers occur within period: ", period_id, call. = FALSE)
  }
  annual_input <- temperature[[which(keep)]]
  years <- format(dates[keep], "%Y")
  annual <- terra::tapp(annual_input, index = years, fun = mean, na.rm = TRUE)
  annual_years <- sort(unique(years))
  if (terra::nlyr(annual) < 3L) {
    stop("At least three annual means are required for an OLS trend.", call. = FALSE)
  }
  terra::time(annual) <- as.Date(paste0(annual_years, "-07-01"))
  if (!is.null(target_grid) && !isTRUE(terra::compareGeom(annual[[1]], target_grid, stopOnError = FALSE))) {
    annual <- terra::resample(annual, target_grid, method = "bilinear")
    terra::time(annual) <- as.Date(paste0(annual_years, "-07-01"))
  }
  annual
}

#' Calculate a per-cell annual-mean OLS temperature slope
#'
#' @param annual_temperature Time-indexed annual-mean `SpatRaster`.
#' @returns Single-layer `SpatRaster` in temperature units per year. Cells with
#'   fewer than three finite annual values are `NA`. The slope uses ordinary
#'   least squares against elapsed decimal years.
revision_ols_trend <- function(annual_temperature) {
  if (!inherits(annual_temperature, "SpatRaster") || terra::nlyr(annual_temperature) < 3L) {
    stop("annual_temperature must contain at least three annual raster layers.", call. = FALSE)
  }
  dates <- terra::time(annual_temperature)
  years <- as.numeric(dates - min(dates)) / 365.25
  terra::app(annual_temperature, fun = function(values) {
    keep <- is.finite(values)
    if (sum(keep) < 3L) return(NA_real_)
    stats::coef(stats::lm(values[keep] ~ years[keep]))[2]
  })
}

#' Calculate a mean temperature field
#'
#' @param annual_temperature Time-indexed annual-mean `SpatRaster`.
#' @returns Single-layer temporal mean `SpatRaster`; cells with no finite value
#'   are `NA`.
revision_temperature_mean <- function(annual_temperature) {
  if (!inherits(annual_temperature, "SpatRaster")) {
    stop("annual_temperature must be a terra::SpatRaster.", call. = FALSE)
  }
  terra::app(annual_temperature, mean, na.rm = TRUE)
}

#' Convert a single-layer raster to a row-major numeric matrix
#'
#' @param raster Single-layer `terra::SpatRaster`.
#' @returns Numeric matrix in raster row/column order.
revision_raster_matrix <- function(raster) {
  if (!inherits(raster, "SpatRaster") || terra::nlyr(raster) != 1L) {
    stop("raster must be a single-layer terra::SpatRaster.", call. = FALSE)
  }
  matrix(terra::values(raster, mat = FALSE), nrow = terra::nrow(raster), ncol = terra::ncol(raster), byrow = TRUE)
}

#' Convert a row-major numeric matrix to a single-layer raster
#'
#' @param values Numeric matrix matching `template` rows and columns.
#' @param template Single-layer `terra::SpatRaster` supplying geometry.
#' @param name Layer name.
#' @returns Single-layer `SpatRaster`. Stops for incompatible dimensions.
revision_matrix_raster <- function(values, template, name) {
  if (!is.matrix(values) || nrow(values) != terra::nrow(template) || ncol(values) != terra::ncol(template)) {
    stop("values dimensions must match template geometry.", call. = FALSE)
  }
  out <- template
  terra::values(out) <- as.vector(t(values))
  names(out) <- name
  out
}

#' Calculate eastward and northward horizontal temperature gradients
#'
#' @param mean_temperature Single-layer temporal mean temperature raster.
#' @returns List of eastward and northward gradients in temperature per metre.
#'   Interior cells use centred differences; edges use one-sided differences.
#'   Gradients adjacent to missing values are `NA` rather than extrapolated.
revision_horizontal_gradient <- function(mean_temperature) {
  x <- revision_raster_matrix(mean_temperature)
  nr <- nrow(x)
  nc <- ncol(x)
  if (nr < 2L || nc < 2L) stop("Temperature grid must have at least two rows and columns.", call. = FALSE)
  rs <- terra::res(mean_temperature)
  ys <- terra::yFromRow(mean_temperature, seq_len(nr))
  dx <- 111320 * cos(ys * pi / 180) * rs[1]
  dy <- 111320 * rs[2]
  east <- matrix(NA_real_, nr, nc)
  north <- matrix(NA_real_, nr, nc)
  for (row in seq_len(nr)) {
    for (col in seq_len(nc)) {
      if (!is.finite(x[row, col])) next
      if (col == 1L) {
        left <- x[row, col]
        right <- x[row, col + 1L]
        divisor_x <- dx[row]
      } else if (col == nc) {
        left <- x[row, col - 1L]
        right <- x[row, col]
        divisor_x <- dx[row]
      } else {
        left <- x[row, col - 1L]
        right <- x[row, col + 1L]
        divisor_x <- 2 * dx[row]
      }
      if (is.finite(left) && is.finite(right) && is.finite(divisor_x) && divisor_x > 0) {
        east[row, col] <- (right - left) / divisor_x
      }
      if (row == 1L) {
        south <- x[row + 1L, col]
        north_value <- x[row, col]
        divisor_y <- dy
      } else if (row == nr) {
        south <- x[row, col]
        north_value <- x[row - 1L, col]
        divisor_y <- dy
      } else {
        south <- x[row + 1L, col]
        north_value <- x[row - 1L, col]
        divisor_y <- 2 * dy
      }
      if (is.finite(south) && is.finite(north_value)) {
        north[row, col] <- (north_value - south) / divisor_y
      }
    }
  }
  list(east = east, north = north)
}

#' Calculate vertical temperature gradients across depth centres
#'
#' @param mean_matrices List of row-major mean-temperature matrices ordered by
#'   increasing positive-downward `depth_centres_m`.
#' @param depth_centres_m Numeric increasing depth centres in metres.
#' @returns List of vertical gradients in temperature per metre. Boundary
#'   layers use one-sided differences and interior layers centred differences.
revision_vertical_gradient <- function(mean_matrices, depth_centres_m) {
  if (length(mean_matrices) != length(depth_centres_m) || length(mean_matrices) < 2L || any(diff(depth_centres_m) <= 0)) {
    stop("mean_matrices and strictly increasing depth_centres_m must have matching length of at least two.", call. = FALSE)
  }
  lapply(seq_along(mean_matrices), function(k) {
    if (k == 1L) {
      (mean_matrices[[2L]] - mean_matrices[[1L]]) / (depth_centres_m[2L] - depth_centres_m[1L])
    } else if (k == length(mean_matrices)) {
      (mean_matrices[[k]] - mean_matrices[[k - 1L]]) / (depth_centres_m[k] - depth_centres_m[k - 1L])
    } else {
      (mean_matrices[[k + 1L]] - mean_matrices[[k - 1L]]) / (depth_centres_m[k + 1L] - depth_centres_m[k - 1L])
    }
  })
}

#' Compute a depth-specific absolute-gradient floor
#'
#' @param gradient Numeric matrix of a horizontal-gradient magnitude or a
#'   vertical gradient at one depth.
#' @param quantile Probability in `[0, 1)`. Zero retains non-zero gradients
#'   without a floor.
#' @returns One non-negative gradient floor. Stops for invalid probabilities or
#'   no finite non-zero gradients.
revision_gradient_floor <- function(gradient, quantile) {
  if (!is.numeric(quantile) || length(quantile) != 1L || !is.finite(quantile) || quantile < 0 || quantile >= 1) {
    stop("gradient floor quantile must be one number in [0, 1).", call. = FALSE)
  }
  values <- abs(as.numeric(gradient))
  values <- values[is.finite(values) & values > 0]
  if (!length(values)) stop("Cannot calculate a gradient floor from no finite non-zero gradients.", call. = FALSE)
  if (identical(quantile, 0)) min(values) else as.numeric(stats::quantile(values, probs = quantile, names = FALSE, type = 8))
}

#' Build a regularised horizontal and vertical velocity bundle
#'
#' @param input_paths Named character vector of existing temperature GeoTIFF
#'   paths. Names must be numeric depth centres in metres, or paths must follow
#'   the `_depth_XXX.Xm` convention.
#' @param dataset_id Dataset label recorded in the bundle and manifest.
#' @param period_id Inclusive `YYYYMMDD-YYYYMMDD` period used for annual means
#'   and OLS trends.
#' @param grid_deg Positive target resolution in degrees.
#' @param depth_centres_m Increasing depth centres to retain.
#' @param gradient_floor_quantile Depth-specific regularisation quantile in
#'   `[0, 1)`; primary runs use `0.01`.
#' @param aggregation_thickness_m `NA` for the established source-centre
#'   workflow, or a positive multiple of the common source-depth spacing in
#'   metres. When supplied, adjacent annual temperature layers are averaged
#'   within depth bins before trends and gradients are calculated. A terminal
#'   incomplete bin is retained and recorded explicitly.
#' @param output_dir New empty output directory for the RDS bundle and CSV
#'   manifest. No existing output is overwritten.
#' @returns Named list containing temperature trends, horizontal and vertical
#'   gradients, HVoCC/VVoCC arrays, valid mask, geometry, source manifest, and
#'   regularisation floors. It also writes `velocity_bundle.rds`,
#'   `velocity_manifest.csv`, and `velocity_distributions.csv`.
#'
#' All input fields must be local time-indexed temperature GeoTIFFs. Missing
#' files, depth coverage, dates, geometry, non-empty output directories, and
#' non-finite field construction stop with an error. This function makes no
#' network calls, never alters its input files, and writes only to a new output
#' directory.
build_revision_velocity_bundle <- function(input_paths, dataset_id, period_id,
                                           grid_deg, depth_centres_m,
                                           gradient_floor_quantile, output_dir,
                                           aggregation_thickness_m = NA_real_) {
  revision_require_velocity_packages()
  if (!is.character(input_paths) || !length(input_paths) || any(!file.exists(input_paths))) {
    stop("input_paths must contain existing local temperature GeoTIFF paths.", call. = FALSE)
  }
  if (!nzchar(dataset_id)) stop("dataset_id must be non-empty.", call. = FALSE)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Velocity output directory is non-empty: ", output_dir, call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Could not create velocity output directory: ", output_dir, call. = FALSE)

  parsed_depths <- suppressWarnings(as.numeric(names(input_paths)))
  if (is.null(names(input_paths)) || anyNA(parsed_depths)) {
    parsed_depths <- vapply(input_paths, revision_depth_from_temperature_path, numeric(1))
  }
  if (anyNA(parsed_depths) || anyDuplicated(parsed_depths)) {
    stop("Every input temperature path must resolve to one unique depth centre.", call. = FALSE)
  }
  source_depths <- sort(unique(parsed_depths))
  requested_source <- sort(unique(as.numeric(depth_centres_m)))
  if (!length(requested_source) || any(diff(requested_source) <= 0)) {
    stop("depth_centres_m must contain increasing finite depth centres.", call. = FALSE)
  }
  lookup <- setNames(input_paths, sprintf("%.4f", parsed_depths))
  source_keys <- sprintf("%.4f", source_depths)
  if (!all(source_keys %in% names(lookup))) {
    stop("Input paths do not cover every parsed source depth.", call. = FALSE)
  }
  aggregate_layers <- !is.na(aggregation_thickness_m)
  if (aggregate_layers && (!is.numeric(aggregation_thickness_m) || length(aggregation_thickness_m) != 1L ||
                           !is.finite(aggregation_thickness_m) || aggregation_thickness_m <= 0)) {
    stop("aggregation_thickness_m must be NA or one positive finite number.", call. = FALSE)
  }
  if (!aggregate_layers) {
    requested <- requested_source
    keys <- sprintf("%.4f", requested)
    if (!all(keys %in% names(lookup))) {
      stop("Input paths do not cover every requested depth: ", paste(requested[!keys %in% names(lookup)], collapse = ", "), call. = FALSE)
    }
    groups <- lapply(seq_along(requested), function(i) match(requested[i], source_depths))
  } else {
    if (length(source_depths) < 2L) stop("Depth aggregation requires at least two source layers.", call. = FALSE)
    if (!isTRUE(all.equal(requested_source, source_depths, tolerance = 1e-10))) {
      stop(
        "Depth aggregation must use exactly the declared prepared source-depth centres; it cannot silently add, omit, or infer layers.",
        call. = FALSE
      )
    }
    spacing <- unique(round(diff(source_depths), 10L))
    if (length(spacing) != 1L || aggregation_thickness_m %% spacing != 0) {
      stop("aggregation_thickness_m must be a positive multiple of the regular source-depth spacing.", call. = FALSE)
    }
    group_size <- as.integer(round(aggregation_thickness_m / spacing))
    if (group_size < 2L) stop("aggregation_thickness_m must combine at least two source layers.", call. = FALSE)
    groups <- unname(split(seq_along(source_depths), ceiling(seq_along(source_depths) / group_size)))
    requested <- unname(vapply(groups, function(index) mean(source_depths[index]), numeric(1L)))
    keys <- source_keys
  }

  first_temperature <- terra::rast(lookup[[source_keys[1L]]])
  target_grid <- revision_target_grid(first_temperature[[1]], grid_deg)
  annual <- vector("list", length(groups))
  means <- vector("list", length(requested))
  trends <- vector("list", length(requested))
  if (!aggregate_layers) {
    for (k in seq_along(requested)) {
      annual[[k]] <- revision_annual_temperature(terra::rast(lookup[[sprintf("%.4f", requested[k])]]), period_id, target_grid)
    }
  } else {
    annual_source <- lapply(source_keys, function(key) {
      revision_annual_temperature(terra::rast(lookup[[key]]), period_id, target_grid)
    })
    annual <- lapply(groups, function(index) {
      reference_time <- as.character(terra::time(annual_source[[index[1L]]]))
      if (any(vapply(index, function(i) !identical(as.character(terra::time(annual_source[[i]])), reference_time), logical(1L)))) {
        stop("Source layers in one aggregate depth bin have incompatible annual time coordinates.", call. = FALSE)
      }
      out <- annual_source[[index[1L]]]
      if (length(index) > 1L) {
        for (i in index[-1L]) out <- out + annual_source[[i]]
      }
      out / length(index)
    })
  }
  for (k in seq_along(requested)) {
    means[[k]] <- revision_temperature_mean(annual[[k]])
    trends[[k]] <- revision_ols_trend(annual[[k]])
  }

  mean_matrices <- lapply(means, revision_raster_matrix)
  trend_matrices <- lapply(trends, revision_raster_matrix)
  vertical_gradient <- revision_vertical_gradient(mean_matrices, requested)
  horizontal_gradient <- lapply(mean_matrices, function(x) {
    revision_horizontal_gradient(revision_matrix_raster(x, target_grid, "mean_temperature"))
  })
  nr <- terra::nrow(target_grid)
  nc <- terra::ncol(target_grid)
  nd <- length(requested)
  hcv_km_yr <- array(NA_real_, dim = c(nr, nc, nd))
  hcv_angle_deg <- array(NA_real_, dim = c(nr, nc, nd))
  vcv_m_yr <- array(NA_real_, dim = c(nr, nc, nd))
  h_floor <- numeric(nd)
  v_floor <- numeric(nd)
  h_gradient_magnitude <- vector("list", nd)

  for (k in seq_len(nd)) {
    grad_e <- horizontal_gradient[[k]]$east
    grad_n <- horizontal_gradient[[k]]$north
    grad_h <- sqrt(grad_e^2 + grad_n^2)
    h_gradient_magnitude[[k]] <- grad_h
    h_floor[k] <- revision_gradient_floor(grad_h, gradient_floor_quantile)
    v_floor[k] <- revision_gradient_floor(vertical_gradient[[k]], gradient_floor_quantile)
    h_denom <- pmax(grad_h, h_floor[k])^2
    h_east <- -trend_matrices[[k]] * grad_e / h_denom
    h_north <- -trend_matrices[[k]] * grad_n / h_denom
    v_denom <- sign(vertical_gradient[[k]]) * pmax(abs(vertical_gradient[[k]]), v_floor[k])
    vcv_m_yr[, , k] <- trend_matrices[[k]] / v_denom
    hcv_km_yr[, , k] <- sqrt(h_east^2 + h_north^2) / 1000
    hcv_angle_deg[, , k] <- (atan2(h_east, h_north) * 180 / pi) %% 360
  }
  valid_mask <- is.finite(hcv_km_yr) & is.finite(hcv_angle_deg) & is.finite(vcv_m_yr)
  hcv_km_yr[!valid_mask] <- NA_real_
  hcv_angle_deg[!valid_mask] <- NA_real_
  vcv_m_yr[!valid_mask] <- NA_real_

  if (!aggregate_layers) {
    paths <- unname(lookup[sprintf("%.4f", requested)])
    source_info <- file.info(paths)
    source_manifest <- data.frame(
      depth_m = requested, source_depths_m = as.character(requested),
      layer_thickness_m = rep(NA_real_, length(requested)),
      path = normalizePath(paths, winslash = "/", mustWork = TRUE),
      bytes = unname(source_info$size), modified_time = as.character(source_info$mtime),
      md5 = unname(tools::md5sum(paths)), stringsAsFactors = FALSE
    )
  } else {
    source_manifest <- do.call(rbind, lapply(seq_along(groups), function(k) {
      paths <- unname(lookup[source_keys[groups[[k]]]])
      source_info <- file.info(paths)
      data.frame(
        depth_m = requested[k], source_depths_m = paste(source_depths[groups[[k]]], collapse = ";"),
        layer_thickness_m = length(groups[[k]]) * spacing,
        path = paste(normalizePath(paths, winslash = "/", mustWork = TRUE), collapse = ";"),
        bytes = sum(source_info$size),
        modified_time = paste(as.character(source_info$mtime), collapse = ";"),
        md5 = paste(unname(tools::md5sum(paths)), collapse = ";"), stringsAsFactors = FALSE
      )
    }))
  }
  distribution <- do.call(rbind, lapply(seq_len(nd), function(k) {
    valid <- valid_mask[, , k]
    data.frame(
      depth_m = requested[k], valid_voxels = sum(valid),
      h_gradient_floor = h_floor[k], v_gradient_floor = v_floor[k],
      hcv_q50_km_yr = stats::median(hcv_km_yr[, , k][valid], na.rm = TRUE),
      hcv_q95_km_yr = as.numeric(stats::quantile(hcv_km_yr[, , k][valid], 0.95, na.rm = TRUE)),
      abs_vcv_q50_m_yr = stats::median(abs(vcv_m_yr[, , k][valid]), na.rm = TRUE),
      abs_vcv_q95_m_yr = as.numeric(stats::quantile(abs(vcv_m_yr[, , k][valid]), 0.95, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
  # Do not persist a SpatRaster or SpatExtent object in the bundle. terra uses
  # external pointers for these objects; an RDS can load successfully while the
  # pointer is no longer valid in a later R process. The arrays below are the
  # complete scientific fields, and this plain geometry is sufficient to
  # reconstruct coordinates without that persistence hazard.
  target_extent <- terra::ext(target_grid)
  target_geometry <- list(
    xmin = as.numeric(target_extent$xmin), xmax = as.numeric(target_extent$xmax),
    ymin = as.numeric(target_extent$ymin), ymax = as.numeric(target_extent$ymax)
  )
  bundle <- list(
    method_version = "connectivity_revision_1.2.0",
    dataset_id = dataset_id, period_id = period_id, grid_deg = grid_deg,
    depths_m = requested, ext = target_geometry,
    res = as.numeric(terra::res(target_grid)), crs_wkt = terra::crs(target_grid),
    x_centers = terra::xFromCol(target_grid, seq_len(nc)),
    y_centers = terra::yFromRow(target_grid, seq_len(nr)),
    annual_mean_ols = TRUE, gradient_floor_quantile = gradient_floor_quantile,
    vertical_method = if (aggregate_layers) "temperature_layer_aggregate" else "source_depth_centres",
    aggregation_thickness_m = if (aggregate_layers) aggregation_thickness_m else NA_real_,
    aggregation_groups = lapply(groups, function(index) source_depths[index]),
    software = list(
      R = as.character(getRversion()),
      terra = as.character(utils::packageVersion("terra"))
    ),
    h_gradient_floor = h_floor, v_gradient_floor = v_floor,
    temperature_trend_c_yr = array(unlist(trend_matrices), dim = c(nr, nc, nd)),
    h_gradient_east_c_m = array(unlist(lapply(horizontal_gradient, `[[`, "east")), dim = c(nr, nc, nd)),
    h_gradient_north_c_m = array(unlist(lapply(horizontal_gradient, `[[`, "north")), dim = c(nr, nc, nd)),
    v_gradient_c_m = array(unlist(vertical_gradient), dim = c(nr, nc, nd)),
    hcv_km_yr = hcv_km_yr, hcv_angle_deg = hcv_angle_deg, vcv_m_yr = vcv_m_yr,
    valid_mask = valid_mask, source_manifest = source_manifest, distribution = distribution
  )
  saveRDS(bundle, file.path(output_dir, "velocity_bundle.rds"))
  utils::write.csv(source_manifest, file.path(output_dir, "velocity_manifest.csv"), row.names = FALSE)
  utils::write.csv(distribution, file.path(output_dir, "velocity_distributions.csv"), row.names = FALSE)
  bundle
}

#' Summarise velocity distributions by depth and named bioregion
#'
#' @param velocity_bundle Bundle returned by `build_revision_velocity_bundle()`.
#' @param bioregion_path Existing vector file readable by `sf::st_read()` and
#'   containing a character `bioregion` attribute.
#' @param output_path New CSV path for the summary. Its parent must already
#'   exist; an existing file is never overwritten.
#' @returns Data frame with valid-cell counts and median/95th percentile HVoCC
#'   and absolute VVoCC by bioregion and depth. HVoCC is km/year and VVoCC is
#'   metres/year. Stops for missing `sf`, invalid bundle geometry, an absent
#'   attribute, points outside all bioregions, or an occupied output path.
#'   This helper performs only local file reads and a local CSV write.
revision_velocity_bioregion_distribution <- function(velocity_bundle,
                                                      bioregion_path,
                                                      output_path) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Missing required package `sf`. Install it outside the analysis scripts, then rerun.", call. = FALSE)
  }
  if (!file.exists(bioregion_path)) stop("Bioregion vector file does not exist: ", bioregion_path, call. = FALSE)
  if (file.exists(output_path)) stop("Refusing to overwrite velocity bioregion summary: ", output_path, call. = FALSE)
  required <- c("depths_m", "x_centers", "y_centers", "crs_wkt", "hcv_km_yr", "vcv_m_yr", "valid_mask")
  if (!is.list(velocity_bundle) || length(setdiff(required, names(velocity_bundle)))) {
    stop("velocity_bundle is missing required geometry or velocity arrays.", call. = FALSE)
  }
  shape <- dim(velocity_bundle$valid_mask)
  if (length(shape) != 3L || !identical(dim(velocity_bundle$hcv_km_yr), shape) ||
      !identical(dim(velocity_bundle$vcv_m_yr), shape)) {
    stop("Velocity bundle has incompatible 3D arrays.", call. = FALSE)
  }
  bioregions <- sf::st_read(bioregion_path, quiet = TRUE)
  if (!"bioregion" %in% names(bioregions)) stop("Bioregion vector lacks a `bioregion` attribute.", call. = FALSE)
  points <- expand.grid(row = seq_len(shape[1L]), col = seq_len(shape[2L]))
  if (length(velocity_bundle$x_centers) != shape[2L] || length(velocity_bundle$y_centers) != shape[1L] ||
      !is.character(velocity_bundle$crs_wkt) || length(velocity_bundle$crs_wkt) != 1L || !nzchar(velocity_bundle$crs_wkt)) {
    stop("Velocity bundle has invalid plain coordinate geometry.", call. = FALSE)
  }
  points$longitude <- velocity_bundle$x_centers[points$col]
  points$latitude <- velocity_bundle$y_centers[points$row]
  points <- sf::st_as_sf(points, coords = c("longitude", "latitude"), crs = velocity_bundle$crs_wkt)
  points <- sf::st_transform(points, sf::st_crs(bioregions))
  points <- sf::st_join(points, bioregions[, "bioregion", drop = FALSE], left = TRUE)
  point_regions <- as.character(sf::st_drop_geometry(points)$bioregion)
  out <- do.call(rbind, lapply(seq_len(shape[3L]), function(k) {
    valid <- velocity_bundle$valid_mask[, , k]
    idx <- which(valid, arr.ind = TRUE)
    if (!nrow(idx)) return(NULL)
    values <- data.frame(
      bioregion = point_regions[idx[, 1L] + (idx[, 2L] - 1L) * shape[1L]],
      depth_m = velocity_bundle$depths_m[k],
      hcv_km_yr = velocity_bundle$hcv_km_yr[, , k][valid],
      abs_vcv_m_yr = abs(velocity_bundle$vcv_m_yr[, , k][valid]),
      stringsAsFactors = FALSE
    )
    values <- values[!is.na(values$bioregion), , drop = FALSE]
    if (!nrow(values)) return(NULL)
    do.call(rbind, lapply(split(values, values$bioregion), function(group) {
      data.frame(
        bioregion = group$bioregion[1L], depth_m = group$depth_m[1L], valid_voxels = nrow(group),
        hcv_q50_km_yr = stats::median(group$hcv_km_yr),
        hcv_q95_km_yr = as.numeric(stats::quantile(group$hcv_km_yr, 0.95)),
        abs_vcv_q50_m_yr = stats::median(group$abs_vcv_m_yr),
        abs_vcv_q95_m_yr = as.numeric(stats::quantile(group$abs_vcv_m_yr, 0.95)),
        stringsAsFactors = FALSE
      )
    }))
  }))
  if (is.null(out) || !nrow(out)) stop("No valid velocity cells intersect the supplied bioregions.", call. = FALSE)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out, output_path, row.names = FALSE)
  out
}
