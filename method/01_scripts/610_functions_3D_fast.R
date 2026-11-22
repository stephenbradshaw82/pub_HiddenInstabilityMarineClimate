# language: R
################################################################################
# 01_scripts/610_functions_3D_fast.R
#
# Legacy helper assumptions:
# - HCV RDS files contain gVoCC_mag (km/yr) and gVoCC_ang (degrees).
# - VCV RDS files contain VCV (m/yr) for one matching depth layer.
# - Depths are positive meters below surface and ordered shallow to deep.
#
# Legacy helper preconditions:
# - HCV and VCV rasters are aligned or can be resampled to a common template.
# - Time period key exists in both HCV and VCV objects, or the first common key
#   is acceptable.
#
# Legacy helper limits:
# - Horizontal movement is sampled from the current voxel only. This is a
#   memory-efficient cell-transition engine, not a full continuous 3D PDE solver.
# - Dateline wrapping is optional and only relevant for global domains.
#
# Legacy helper failure modes:
# - Missing files or missing time-period keys.
# - PackedSpatRaster inputs that need unwrapping.
# - Empty rasters, geometry mismatches, invalid seeds, and particles leaving
#   the domain or entering NA cells.
#
# Current workflow:
# - Scripts 05 and 06 source this file and use velocity bundles created by
#   script 03.
#
# Expected output:
# - A list with velocity arrays, seed table, count arrays, diagnostics, and
#   integer classification rasters ready for export.
################################################################################

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Extract a numeric depth from a standard velocity filename
#'
#' @param fname Input filename or path containing `depth_<metres>m`.
#' @returns Numeric depth in metres, or `NA_real_` when the pattern is absent.
extract_depth_3d_fast <- function(fname) {
  as.numeric(stringr::str_match(basename(fname), "depth_([0-9\\.]+)m")[, 2])
}

#' Convert a supported raster representation to a terra raster
#'
#' @param x A `PackedSpatRaster`, `SpatRaster`, or legacy raster object.
#' @returns A `terra::SpatRaster`. Stops for unsupported classes.
unwrap_to_spatraster_3d_fast <- function(x) {
  if (inherits(x, "PackedSpatRaster")) {
    return(terra::unwrap(x))
  }
  if (inherits(x, "SpatRaster")) {
    return(terra::rast(x))
  }
  if (inherits(x, c("RasterLayer", "RasterBrick", "RasterStack"))) {
    if (!raster::hasValues(x)) {
      x <- raster::readAll(x)
    }
    return(terra::rast(x))
  }
  stop("Unsupported raster object class: ", paste(class(x), collapse = "/"), call. = FALSE)
}

#' Locate a raster inside a supported object or nested list
#'
#' @param x Raster object or non-empty nested list containing a supported
#'   raster under common HVoCC/VVoCC field names.
#' @returns A `terra::SpatRaster`. Stops when no supported raster is found.
find_spatraster_3d_fast <- function(x) {
  if (inherits(x, c("PackedSpatRaster", "SpatRaster", "RasterLayer", "RasterBrick", "RasterStack"))) {
    return(unwrap_to_spatraster_3d_fast(x))
  }
  if (is.list(x) && length(x) > 0) {
    keys <- c("gVoCC_mag", "gVoCC_ang", "VCV", "mag", "ang", "rast", "raster")
    for (nm in keys) {
      if (!is.null(x[[nm]])) {
        out <- try(find_spatraster_3d_fast(x[[nm]]), silent = TRUE)
        if (!inherits(out, "try-error")) {
          return(out)
        }
      }
    }
    for (el in x) {
      out <- try(find_spatraster_3d_fast(el), silent = TRUE)
      if (!inherits(out, "try-error")) {
        return(out)
      }
    }
  }
  stop("Could not locate a raster within the supplied object.", call. = FALSE)
}

#' Resolve a shared HVoCC/VVoCC time-period key
#'
#' @param hcv_obj Deserialised HVoCC object with named period entries.
#' @param vcv_obj Deserialised VVoCC object with named period entries.
#' @param time_period Optional requested key; `NULL` selects the first shared
#'   key.
#' @returns One shared time-period key. Stops when a requested key is absent or
#'   no common key exists.
resolve_time_period_3d_fast <- function(hcv_obj, vcv_obj, time_period = NULL) {
  h_keys <- names(hcv_obj)
  v_keys <- names(vcv_obj)
  common_keys <- intersect(h_keys, v_keys)

  if (!is.null(time_period)) {
    if (!(time_period %in% common_keys)) {
      stop("Requested time_period not found in both HCV and VCV objects: ", time_period, call. = FALSE)
    }
    return(time_period)
  }

  if (length(common_keys) == 0L) {
    stop("No common time period key found between HCV and VCV objects.", call. = FALSE)
  }

  common_keys[1]
}

#' Unwrap unnamed singleton lists
#'
#' @param x Object that may be nested in unnamed length-one lists.
#' @returns The innermost non-singleton or named object.
unwrap_singleton_list_3d_fast <- function(x) {
  out <- x
  while (is.list(out) && length(out) == 1L && is.null(names(out))) {
    out <- out[[1]]
  }
  out
}

#' Read one aligned HVoCC/VVoCC depth pair
#'
#' @param hcv_file HVoCC RDS path containing magnitude and angle rasters.
#' @param vcv_file VVoCC RDS path containing a vertical-velocity raster.
#' @param time_period Optional shared time key; `NULL` selects the first common
#'   key.
#' @param template Optional target `SpatRaster`; nonmatching fields are bilinear
#'   resampled to this geometry.
#' @returns List containing the selected time period plus `mag`, `ang`, and
#'   `vcv` `SpatRaster` objects. Stops for unreadable files, missing keys, or
#'   unsupported raster structures.
extract_hcv_vcv_pair_3d_fast <- function(hcv_file, vcv_file, time_period = NULL, template = NULL) {
  hcv_obj <- readRDS(hcv_file)
  vcv_obj <- readRDS(vcv_file)
  tp <- resolve_time_period_3d_fast(hcv_obj, vcv_obj, time_period = time_period)

  h_entry <- unwrap_singleton_list_3d_fast(hcv_obj[[tp]])
  v_entry <- unwrap_singleton_list_3d_fast(vcv_obj[[tp]])

  mag <- find_spatraster_3d_fast(h_entry$gVoCC_mag %||% h_entry$mag %||% h_entry)
  ang <- find_spatraster_3d_fast(h_entry$gVoCC_ang %||% h_entry$ang)
  vcv <- find_spatraster_3d_fast(v_entry$VCV %||% v_entry)

  if (!is.null(template)) {
    same_grid <- function(a, b) {
      identical(dim(a), dim(b)) &&
        isTRUE(all.equal(as.vector(terra::ext(a)), as.vector(terra::ext(b)), tolerance = 1e-10)) &&
        isTRUE(all.equal(terra::res(a), terra::res(b), tolerance = 1e-10))
    }

    if (!same_grid(mag, template)) {
      mag <- terra::resample(mag, template, method = "bilinear")
    }
    if (!same_grid(ang, template)) {
      ang <- terra::resample(ang, template, method = "bilinear")
    }
    if (!same_grid(vcv, template)) {
      vcv <- terra::resample(vcv, template, method = "bilinear")
    }
  }

  list(
    time_period = tp,
    mag = mag,
    ang = ang,
    vcv = vcv
  )
}

#' Convert a single raster layer to a row-major matrix
#'
#' @param r Single-layer `terra::SpatRaster`.
#' @returns Numeric matrix ordered by raster row then column.
raster_to_matrix_3d_fast <- function(r) {
  vals <- terra::values(r, mat = FALSE)
  matrix(vals, nrow = terra::nrow(r), ncol = terra::ncol(r), byrow = TRUE)
}

#' Prepare matched HVoCC/VVoCC arrays for particle simulation
#'
#' @param hcv_files HVoCC RDS paths, one per depth.
#' @param vcv_files VVoCC RDS paths, one per depth.
#' @param depths Optional numeric depth centres in metres; `NULL` derives them
#'   from filenames.
#' @param time_period Optional shared time key.
#' @returns List containing aligned geometry, depth edges, HVoCC magnitude and
#'   degree-bearing arrays, east/north components in metres/year, VVoCC in
#'   metres/year, and the valid-voxel mask. Stops for missing depth coverage,
#'   invalid filenames, or absent shared periods.
prepare_velocity_fields_3d_fast <- function(hcv_files, vcv_files, depths = NULL, time_period = NULL) {
  if (length(hcv_files) == 0L || length(vcv_files) == 0L) {
    stop("At least one HCV file and one VCV file are required.", call. = FALSE)
  }

  if (is.null(depths)) {
    depths <- vapply(hcv_files, extract_depth_3d_fast, numeric(1))
  }

  if (anyNA(depths)) {
    stop("Could not extract all depths from filenames. Supply depths explicitly.", call. = FALSE)
  }

  hcv_depths <- vapply(hcv_files, extract_depth_3d_fast, numeric(1))
  vcv_depths <- vapply(vcv_files, extract_depth_3d_fast, numeric(1))
  common_depths <- sort(intersect(unique(hcv_depths), unique(vcv_depths)))

  if (length(common_depths) == 0L) {
    stop("No common depths between HCV and VCV files.", call. = FALSE)
  }

  if (!all(depths %in% common_depths)) {
    stop("Requested depths are not all available in both HCV and VCV files.", call. = FALSE)
  }

  h_lookup <- setNames(hcv_files, sprintf("%.4f", hcv_depths))
  v_lookup <- setNames(vcv_files, sprintf("%.4f", vcv_depths))
  depths <- sort(unique(as.numeric(depths)))

  template_pair <- extract_hcv_vcv_pair_3d_fast(
    hcv_file = h_lookup[[sprintf("%.4f", depths[1])]],
    vcv_file = v_lookup[[sprintf("%.4f", depths[1])]],
    time_period = time_period
  )
  template <- template_pair$mag

  nrow_r <- terra::nrow(template)
  ncol_r <- terra::ncol(template)
  ndepth <- length(depths)

  hcv_mag_km_yr <- array(NA_real_, dim = c(nrow_r, ncol_r, ndepth))
  hcv_ang_deg <- array(NA_real_, dim = c(nrow_r, ncol_r, ndepth))
  east_m_yr <- array(NA_real_, dim = c(nrow_r, ncol_r, ndepth))
  north_m_yr <- array(NA_real_, dim = c(nrow_r, ncol_r, ndepth))
  vcv_m_yr <- array(NA_real_, dim = c(nrow_r, ncol_r, ndepth))
  valid_mask <- array(FALSE, dim = c(nrow_r, ncol_r, ndepth))

  ex <- terra::ext(template)
  rs <- terra::res(template)
  x_centers <- seq(ex$xmin + rs[1] / 2, ex$xmax - rs[1] / 2, by = rs[1])
  y_centers <- seq(ex$ymax - rs[2] / 2, ex$ymin + rs[2] / 2, by = -rs[2])
  lon_m_per_deg <- 111320 * cos(y_centers * pi / 180)
  lat_m_per_deg <- rep(111320, length(y_centers))
  cell_width_m <- rs[1] * lon_m_per_deg
  cell_height_m <- rep(rs[2] * 111320, length(y_centers))

  period_used <- NULL

  for (k in seq_along(depths)) {
    key <- sprintf("%.4f", depths[k])
    pair <- extract_hcv_vcv_pair_3d_fast(
      hcv_file = h_lookup[[key]],
      vcv_file = v_lookup[[key]],
      time_period = time_period,
      template = template
    )
    period_used <- period_used %||% pair$time_period

    mag_mat <- raster_to_matrix_3d_fast(pair$mag)
    ang_mat <- raster_to_matrix_3d_fast(pair$ang)
    vcv_mat <- raster_to_matrix_3d_fast(pair$vcv)

    theta <- ang_mat * pi / 180
    east_mat <- mag_mat * 1000 * sin(theta)
    north_mat <- mag_mat * 1000 * cos(theta)

    hcv_mag_km_yr[, , k] <- mag_mat
    hcv_ang_deg[, , k] <- ang_mat
    east_m_yr[, , k] <- east_mat
    north_m_yr[, , k] <- north_mat
    vcv_m_yr[, , k] <- vcv_mat
    valid_mask[, , k] <- !is.na(mag_mat) & !is.na(ang_mat) & !is.na(vcv_mat)
  }

  depth_edges <- numeric(ndepth + 1)
  if (ndepth == 1L) {
    depth_edges[] <- c(depths[1] - 0.5, depths[1] + 0.5)
  } else {
    mids <- (depths[-1] + depths[-ndepth]) / 2
    depth_edges[1] <- depths[1] - (mids[1] - depths[1])
    depth_edges[ndepth + 1] <- depths[ndepth] + (depths[ndepth] - mids[length(mids)])
    depth_edges[2:ndepth] <- mids
  }

  list(
    template = template,
    time_period = period_used,
    depths = depths,
    depth_edges = depth_edges,
    x_centers = x_centers,
    y_centers = y_centers,
    ext = ex,
    res = rs,
    nrow = nrow_r,
    ncol = ncol_r,
    ndepth = ndepth,
    lon_m_per_deg = lon_m_per_deg,
    lat_m_per_deg = lat_m_per_deg,
    cell_width_m = cell_width_m,
    cell_height_m = cell_height_m,
    hcv_mag_km_yr = hcv_mag_km_yr,
    hcv_ang_deg = hcv_ang_deg,
    east_m_yr = east_m_yr,
    north_m_yr = north_m_yr,
    vcv_m_yr = vcv_m_yr,
    valid_mask = valid_mask
  )
}

#' Apply observed or zero vertical velocity to prepared fields
#'
#' @param fields Prepared list from `prepare_velocity_fields_3d_fast()`.
#' @param vertical_mode `"observed"` retains VVoCC; `"zero"` sets VVoCC to zero
#'   only in valid voxels for the matched control.
#' @returns Copy of `fields` with `vertical_mode` recorded and, where requested,
#'   zeroed valid VVoCC values.
apply_vertical_mode_3d_fast <- function(fields, vertical_mode = c("observed", "zero")) {
  vertical_mode <- match.arg(vertical_mode)
  out <- fields
  out$vertical_mode <- vertical_mode

  if (identical(vertical_mode, "zero")) {
    vcv_zero <- out$vcv_m_yr
    vcv_zero[out$valid_mask] <- 0
    out$vcv_m_yr <- vcv_zero
  }

  out
}

#' Create reproducible particle starts from valid voxels
#'
#' @param fields Prepared velocity-field list with a valid-voxel mask.
#' @param seed_mode `"all"` seeds every valid voxel; `"sample"` samples valid
#'   voxels independently at each depth.
#' @param trajSt Positive particles per valid voxel in `"all"` mode.
#' @param n_seeds_per_depth Positive sample size required in `"sample"` mode.
#' @param seed Integer random seed.
#' @param jitter Logical; horizontally jitter starts within cell boundaries.
#' @returns Data frame of seed IDs, indices, and x/y/depth starting positions.
#'   Stops when no valid starts are produced or sampling settings are invalid.
seed_table_3d_fast <- function(fields, seed_mode = c("all", "sample"), trajSt = 1L, n_seeds_per_depth = NULL, seed = 1L, jitter = TRUE) {
  seed_mode <- match.arg(seed_mode)
  set.seed(seed)

  out <- vector("list", fields$ndepth)
  seed_id <- 1L

  for (k in seq_len(fields$ndepth)) {
    valid_k <- fields$valid_mask[, , k, drop = FALSE]
    dim(valid_k) <- dim(fields$valid_mask)[1:2]
    valid_idx <- which(valid_k, arr.ind = TRUE)
    if (nrow(valid_idx) == 0L) {
      next
    }

    if (seed_mode == "sample") {
      if (is.null(n_seeds_per_depth) || n_seeds_per_depth < 1L) {
        stop("n_seeds_per_depth must be provided for seed_mode='sample'.", call. = FALSE)
      }
      picks <- sample.int(nrow(valid_idx), size = n_seeds_per_depth, replace = nrow(valid_idx) < n_seeds_per_depth)
      valid_idx <- valid_idx[picks, , drop = FALSE]
      reps <- rep.int(1L, nrow(valid_idx))
    } else {
      reps <- rep.int(as.integer(trajSt), nrow(valid_idx))
    }

    row_rep <- rep(valid_idx[, 1], reps)
    col_rep <- rep(valid_idx[, 2], reps)
    n_seed <- length(row_rep)
    if (n_seed == 0L) {
      next
    }

    x <- fields$x_centers[col_rep]
    y <- fields$y_centers[row_rep]

    if (isTRUE(jitter)) {
      x <- x + stats::runif(n_seed, -0.5, 0.5) * fields$res[1]
      y <- y + stats::runif(n_seed, -0.5, 0.5) * fields$res[2]
      x <- pmin(pmax(x, fields$ext$xmin + .Machine$double.eps), fields$ext$xmax - .Machine$double.eps)
      y <- pmin(pmax(y, fields$ext$ymin + .Machine$double.eps), fields$ext$ymax - .Machine$double.eps)
    }

    out[[k]] <- data.frame(
      seed_id = seq.int(seed_id, length.out = n_seed),
      row = row_rep,
      col = col_rep,
      k = rep.int(k, n_seed),
      x = x,
      y = y,
      z = rep.int(fields$depths[k], n_seed),
      stringsAsFactors = FALSE
    )
    seed_id <- seed_id + n_seed
  }

  seeds <- dplyr::bind_rows(out)
  if (nrow(seeds) == 0L) {
    stop("No valid seeds generated from the available rasters.", call. = FALSE)
  }

  seeds
}

#' Summarise valid cells and generated starts by depth
#'
#' @param fields Prepared velocity-field list.
#' @param seeds Seed table from `seed_table_3d_fast()`.
#' @returns Data frame of valid cells, seeded particles, unique seeded cells,
#'   and valid-start checks by depth.
seed_summary_3d_fast <- function(fields, seeds) {
  valid_cells_per_depth <- vapply(seq_len(fields$ndepth), function(k) {
    sum(fields$valid_mask[, , k], na.rm = TRUE)
  }, numeric(1))

  seeded_per_depth <- as.integer(table(factor(seeds$k, levels = seq_len(fields$ndepth))))

  data.frame(
    depth_index = seq_len(fields$ndepth),
    depth_m = fields$depths,
    valid_cells = valid_cells_per_depth,
    seeded_particles = seeded_per_depth,
    seeded_unique_cells = vapply(seq_len(fields$ndepth), function(k) {
      tmp <- seeds[seeds$k == k, c("row", "col"), drop = FALSE]
      if (nrow(tmp) == 0L) {
        return(0)
      }
      nrow(unique(tmp))
    }, numeric(1)),
    starts_only_on_valid = vapply(seq_len(fields$ndepth), function(k) {
      tmp <- seeds[seeds$k == k, , drop = FALSE]
      if (nrow(tmp) == 0L) {
        return(TRUE)
      }
      all(fields$valid_mask[cbind(tmp$row, tmp$col, tmp$k)])
    }, logical(1)),
    stringsAsFactors = FALSE
  )
}

#' Convert continuous longitude/latitude to raster row and column
#'
#' @param x Longitude coordinate or vector.
#' @param y Latitude coordinate or vector.
#' @param fields Prepared velocity-field list defining extent and resolution.
#' @returns List of integer-like `row` and `col` vectors; range validation is
#'   performed by the caller.
row_col_from_xy_3d_fast <- function(x, y, fields) {
  col <- floor((x - fields$ext$xmin) / fields$res[1]) + 1L
  row <- floor((fields$ext$ymax - y) / fields$res[2]) + 1L
  list(row = row, col = col)
}

#' Convert positive-downward depth to a depth-layer index
#'
#' @param z Depth coordinate or vector in metres.
#' @param fields Prepared velocity-field list defining depth edges.
#' @returns Integer depth index, with out-of-range values handled by the caller.
depth_index_from_z_3d_fast <- function(z, fields) {
  findInterval(z, fields$depth_edges, rightmost.closed = TRUE)
}

#' Encode row, column, and depth indices as a linear array index
#'
#' @param row Raster-row index or vector.
#' @param col Raster-column index or vector.
#' @param k Depth-layer index or vector.
#' @param fields Prepared velocity-field list defining dimensions.
#' @returns Integer-like linear array index. Inputs must already be in range.
linear_index_3d_fast <- function(row, col, k, fields) {
  row + (col - 1L) * fields$nrow + (k - 1L) * fields$nrow * fields$ncol
}

#' Accumulate visits into a 3D voxel-count array
#'
#' @param row Raster-row indices.
#' @param col Raster-column indices.
#' @param k Depth-layer indices.
#' @param fields Prepared velocity-field list defining dimensions.
#' @returns Integer start/end/visit count array with field dimensions.
accumulate_counts_3d_fast <- function(row, col, k, fields) {
  idx <- linear_index_3d_fast(row, col, k, fields)
  tab <- tabulate(idx, nbins = fields$nrow * fields$ncol * fields$ndepth)
  array(as.integer(tab), dim = c(fields$nrow, fields$ncol, fields$ndepth))
}

#' Step and tally one production-style particle chunk
#'
#' @param seed_chunk Seed-table subset from `seed_table_3d_fast()`.
#' @param fields Prepared velocity-field list, optionally modified by
#'   `apply_vertical_mode_3d_fast()`.
#' @param dt_years Positive Euler time step in years.
#' @param tyr Positive projection horizon in years.
#' @param sample_path_ids Seed IDs for paths retained as diagnostics.
#' @param record_every Positive integration-step interval for retained paths.
#' @param correct Logical date-line correction switch.
#' @param DateLine Logical; whether horizontal positions wrap at the date line.
#' @returns List with chunk-local start, end, and through arrays; per-seed
#'   diagnostics; and sampled paths. Latitude is updated before longitude;
#'   rejected proposals remain in continuous endpoint diagnostics while the end
#'   voxel remains the last valid index. Through counts are unique visits within
#'   this chunk after start/end subtraction, so aggregate values can depend on
#'   `chunk_size`.
simulate_chunk_3d_fast <- function(seed_chunk, fields, dt_years, tyr, sample_path_ids = integer(0), record_every = 1L, correct = TRUE, DateLine = FALSE) {
  steps <- as.integer(round(tyr / dt_years))
  if (steps < 1L) {
    stop("tyr / dt_years must be at least 1.", call. = FALSE)
  }

  n <- nrow(seed_chunk)
  row <- as.integer(seed_chunk$row)
  col <- as.integer(seed_chunk$col)
  k <- as.integer(seed_chunk$k)
  x <- as.numeric(seed_chunk$x)
  y <- as.numeric(seed_chunk$y)
  z <- as.numeric(seed_chunk$z)
  start_row <- row
  start_col <- col
  start_k <- k
  start_x <- x
  start_y <- y
  start_z <- z
  active <- rep(TRUE, n)
  moved_from_start <- rep(FALSE, n)
  seed_idx <- seq_len(n)

  sample_path_ids <- intersect(sample_path_ids, seed_chunk$seed_id)
  sample_paths <- vector("list", length(sample_path_ids))
  names(sample_paths) <- as.character(sample_path_ids)
  visit_seed_parts <- vector("list", steps + 1L)
  visit_lin_parts <- vector("list", steps + 1L)

  start_lin_idx <- linear_index_3d_fast(start_row, start_col, start_k, fields)
  visit_seed_parts[[1]] <- seed_idx
  visit_lin_parts[[1]] <- start_lin_idx

  if (length(sample_path_ids) > 0L) {
    for (sid in sample_path_ids) {
      ii <- which(seed_chunk$seed_id == sid)[1]
      sample_paths[[as.character(sid)]] <- data.frame(
        seed_id = sid,
        step = 0L,
        x = start_x[ii],
        y = start_y[ii],
        z = start_z[ii],
        row = start_row[ii],
        col = start_col[ii],
        k = start_k[ii]
      )
    }
  }

  for (step in seq_len(steps)) {
    idx_active <- which(active)
    if (length(idx_active) == 0L) {
      break
    }

    lin_idx <- linear_index_3d_fast(row[idx_active], col[idx_active], k[idx_active], fields)
    east <- fields$east_m_yr[lin_idx]
    north <- fields$north_m_yr[lin_idx]
    vert <- fields$vcv_m_yr[lin_idx]

    valid_vel <- !is.na(east) & !is.na(north) & !is.na(vert)
    if (!all(valid_vel)) {
      active[idx_active[!valid_vel]] <- FALSE
      idx_active <- idx_active[valid_vel]
      east <- east[valid_vel]
      north <- north[valid_vel]
      vert <- vert[valid_vel]
    }

    if (length(idx_active) == 0L) {
      next
    }

    y[idx_active] <- y[idx_active] + (north * dt_years) / 111320
    lon_scale <- 111320 * cos(y[idx_active] * pi / 180)
    lon_scale[abs(lon_scale) < 1e-8] <- NA_real_
    x[idx_active] <- x[idx_active] + (east * dt_years) / lon_scale
    # Depth is positive downward: negative VCV deepens, positive VCV shoals.
    z[idx_active] <- z[idx_active] - vert * dt_years

    if (isTRUE(DateLine) && isTRUE(correct)) {
      width <- fields$ext$xmax - fields$ext$xmin
      x[idx_active] <- ((x[idx_active] - fields$ext$xmin) %% width) + fields$ext$xmin
    }

    rc <- row_col_from_xy_3d_fast(x[idx_active], y[idx_active], fields)
    new_row <- as.integer(rc$row)
    new_col <- as.integer(rc$col)
    new_k <- as.integer(depth_index_from_z_3d_fast(z[idx_active], fields))

    inside <- new_row >= 1L & new_row <= fields$nrow &
      new_col >= 1L & new_col <= fields$ncol &
      new_k >= 1L & new_k <= fields$ndepth

    valid_cell <- rep(FALSE, length(idx_active))
    if (any(inside)) {
      new_lin <- linear_index_3d_fast(new_row[inside], new_col[inside], new_k[inside], fields)
      valid_cell[inside] <- fields$valid_mask[new_lin]
    }

    good <- inside & valid_cell
    bad_global <- idx_active[!good]
    if (length(bad_global) > 0L) {
      active[bad_global] <- FALSE
    }

    if (any(good)) {
      good_global <- idx_active[good]
      new_lin_good <- linear_index_3d_fast(new_row[good], new_col[good], new_k[good], fields)
      row[good_global] <- new_row[good]
      col[good_global] <- new_col[good]
      k[good_global] <- new_k[good]
      moved_from_start[good_global] <- moved_from_start[good_global] |
        row[good_global] != start_row[good_global] |
        col[good_global] != start_col[good_global] |
        k[good_global] != start_k[good_global]
      visit_seed_parts[[step + 1L]] <- good_global
      visit_lin_parts[[step + 1L]] <- new_lin_good
    }

    if (length(sample_path_ids) > 0L && (step %% record_every) == 0L) {
      for (sid in sample_path_ids) {
        ii <- which(seed_chunk$seed_id == sid)[1]
        sample_paths[[as.character(sid)]] <- rbind(
          sample_paths[[as.character(sid)]],
          data.frame(
            seed_id = sid,
            step = step,
            x = x[ii],
            y = y[ii],
            z = z[ii],
            row = row[ii],
            col = col[ii],
            k = k[ii]
          )
        )
      }
    }
  }

  start_counts <- accumulate_counts_3d_fast(start_row, start_col, start_k, fields)
  end_counts <- accumulate_counts_3d_fast(row, col, k, fields)

  visit_seed <- unlist(visit_seed_parts, use.names = FALSE)
  visit_lin <- unlist(visit_lin_parts, use.names = FALSE)
  unique_visits <- unique(data.frame(
    seed_idx = visit_seed,
    lin_idx = visit_lin
  ))
  visit_tab <- tabulate(
    unique_visits$lin_idx,
    nbins = fields$nrow * fields$ncol * fields$ndepth
  )
  total_visit_counts <- array(
    as.integer(visit_tab),
    dim = c(fields$nrow, fields$ncol, fields$ndepth)
  )
  through_counts <- total_visit_counts - start_counts - end_counts
  through_counts[through_counts < 0L] <- 0L

  list(
    start = start_counts,
    end = end_counts,
    through = through_counts,
    diagnostics = data.frame(
      seed_id = seed_chunk$seed_id,
      start_row = start_row,
      start_col = start_col,
      start_k = start_k,
      end_row = row,
      end_col = col,
      end_k = k,
      start_x = start_x,
      start_y = start_y,
      start_z = start_z,
      end_x = x,
      end_y = y,
      end_z = z,
      moved_from_start = moved_from_start,
      horizontal_distance_km = geosphere::distHaversine(
        cbind(start_x, start_y),
        cbind(x, y)
      ) / 1000,
      vertical_distance_m = abs(z - start_z),
      active_at_end = active,
      stringsAsFactors = FALSE
    ),
    sample_paths = sample_paths
  )
}

#' Run the chunked transition engine across all seeds
#'
#' @param fields Prepared velocity-field list.
#' @param seeds Complete seed table.
#' @param dt_years Positive Euler time step in years.
#' @param tyr Positive projection horizon in years.
#' @param cores Requested worker count; Windows uses serial processing.
#' @param chunk_size Positive number of particles per tally chunk.
#' @param sample_paths_n Number of seed paths retained for diagnostics.
#' @param record_every Positive integration-step interval for retained paths.
#' @param correct Logical date-line correction switch.
#' @param DateLine Logical; whether horizontal positions wrap at the date line.
#' @returns List of aggregated start/end/through arrays, per-seed diagnostics,
#'   and sampled paths. The tally preserves the production chunk convention.
run_transition_engine_3d_fast <- function(fields,
                                          seeds,
                                          dt_years = 0.5,
                                          tyr = 30,
                                          cores = max(1L, parallel::detectCores() - 2L),
                                          chunk_size = 20000L,
                                          sample_paths_n = 250L,
                                          record_every = 1L,
                                          correct = TRUE,
                                          DateLine = FALSE) {
  split_idx <- split(seq_len(nrow(seeds)), ceiling(seq_len(nrow(seeds)) / chunk_size))
  sample_path_ids <- if (sample_paths_n > 0L) {
    sample(seeds$seed_id, size = min(sample_paths_n, nrow(seeds)))
  } else {
    integer(0)
  }

  worker_fun <- function(idx) {
    simulate_chunk_3d_fast(
      seed_chunk = seeds[idx, , drop = FALSE],
      fields = fields,
      dt_years = dt_years,
      tyr = tyr,
      sample_path_ids = sample_path_ids,
      record_every = record_every,
      correct = correct,
      DateLine = DateLine
    )
  }

  if (.Platform$OS.type == "windows" || cores <= 1L) {
    parts <- lapply(split_idx, worker_fun)
  } else {
    parts <- parallel::mclapply(split_idx, worker_fun, mc.cores = cores)
  }

  start <- array(0L, dim = c(fields$nrow, fields$ncol, fields$ndepth))
  end <- array(0L, dim = c(fields$nrow, fields$ncol, fields$ndepth))
  through <- array(0L, dim = c(fields$nrow, fields$ncol, fields$ndepth))
  diag_df <- vector("list", length(parts))
  sample_paths <- list()

  for (i in seq_along(parts)) {
    start <- start + parts[[i]]$start
    end <- end + parts[[i]]$end
    through <- through + parts[[i]]$through
    diag_df[[i]] <- parts[[i]]$diagnostics
    sample_paths <- c(sample_paths, parts[[i]]$sample_paths)
  }

  sample_paths <- sample_paths[!duplicated(names(sample_paths))]

  list(
    start = start,
    end = end,
    through = through,
    diagnostics = dplyr::bind_rows(diag_df),
    sample_paths = sample_paths
  )
}

#' Classify matched-particle voxel tallies
#'
#' @param fields Prepared velocity-field list.
#' @param counts Start/end/through arrays from `run_transition_engine_3d_fast()`.
#' @param seeds Optional complete seed table used to construct diagnostic depth
#'   percentages when `seeds_per_depth` is not supplied.
#' @param seeds_per_depth Optional non-negative numeric vector giving the number
#'   of seeded particles at each depth. This is required when starts are
#'   streamed from a memory-bounded seed plan rather than materialised.
#' @param tyr Projection horizon in years, retained for run provenance and the
#'   default `class_tyr` value.
#' @param class_tyr Horizon in years for local velocity-derived movement
#'   distance.
#' @param nmL Non-moving distance limit in kilometres over `class_tyr`.
#' @param smL Slow-moving upper distance limit in kilometres over `class_tyr`.
#' @param Nend Local ending-percentage threshold for sink eligibility.
#' @param Nst Local low-starting-percentage threshold for Relative Sinks.
#' @param NFT Local through-flow-percentage threshold for Corridors.
#' @param sink_balance_tol Minimum local end-minus-start difference in
#'   percentage points for Internal/Boundary Sinks.
#' @param corridor_balance_tol Maximum absolute local start/end difference in
#'   percentage points for Corridors.
#' @param min_count Minimum ending count for Internal/Boundary Sink eligibility
#'   only.
#' @param norm_mode Diagnostic normalisation retained as total or depth shares;
#'   it does not set the class denominator.
#' @returns List of diagnostic total/depth percentages, local class percentages,
#'   class IDs and labels, movement distance, boundary mask, and class lookup.
#'
#' Class thresholds use local `start + end + through` proportions. `D` is
#' derived from local HVoCC/VVoCC and `class_tyr`, not particle endpoint
#' displacement. Non-moving (`D <= nmL`) and Slow-moving (`nmL < D <= smL`)
#' labels are assigned first and protected from later sink assignment.
classify_3d_fast <- function(fields,
                             counts,
                             seeds = NULL,
                             seeds_per_depth = NULL,
                             tyr,
                             class_tyr = tyr,
                             nmL = 5,
                             smL = 20,
                             Nend = 10,
                             Nst = 10,
                             NFT = 10,
                             sink_balance_tol = 5,
                             corridor_balance_tol = 5,
                             min_count = 1L,
                             norm_mode = c("total", "depth")) {
  norm_mode <- match.arg(norm_mode)
  total_traj <- sum(counts$start)
  if (is.null(seeds_per_depth)) {
    if (!is.data.frame(seeds) || !"k" %in% names(seeds)) {
      stop("Supply a seed table with k or seeds_per_depth for classification.", call. = FALSE)
    }
    seeds_per_depth <- as.integer(table(factor(seeds$k, levels = seq_len(fields$ndepth))))
  } else {
    if (!is.numeric(seeds_per_depth) || length(seeds_per_depth) != fields$ndepth ||
        any(!is.finite(seeds_per_depth)) || any(seeds_per_depth < 0)) {
      stop("seeds_per_depth must be a non-negative numeric vector for every depth.", call. = FALSE)
    }
    seeds_per_depth <- as.integer(seeds_per_depth)
  }
  local_total <- counts$start + counts$end + counts$through

  pct_total_start <- if (total_traj > 0) 100 * counts$start / total_traj else counts$start * 0
  pct_total_end <- if (total_traj > 0) 100 * counts$end / total_traj else counts$end * 0
  pct_total_through <- if (total_traj > 0) 100 * counts$through / total_traj else counts$through * 0

  pct_depth_start <- array(0, dim = dim(counts$start))
  pct_depth_end <- array(0, dim = dim(counts$start))
  pct_depth_through <- array(0, dim = dim(counts$start))

  for (k in seq_len(fields$ndepth)) {
    if (seeds_per_depth[k] > 0L) {
      pct_depth_start[, , k] <- 100 * counts$start[, , k] / seeds_per_depth[k]
      pct_depth_end[, , k] <- 100 * counts$end[, , k] / seeds_per_depth[k]
      pct_depth_through[, , k] <- 100 * counts$through[, , k] / seeds_per_depth[k]
    }
  }

  # `trajClas` thresholds are based on local voxel proportions:
  # Nend / Nst / NFT are percentages of the voxel's own
  # start + through + end tally, not percentages of all seeded particles.
  pct_local_start <- array(0, dim = dim(counts$start))
  pct_local_end <- array(0, dim = dim(counts$start))
  pct_local_through <- array(0, dim = dim(counts$start))
  has_local_total <- local_total > 0
  pct_local_start[has_local_total] <- 100 * counts$start[has_local_total] / local_total[has_local_total]
  pct_local_end[has_local_total] <- 100 * counts$end[has_local_total] / local_total[has_local_total]
  pct_local_through[has_local_total] <- 100 * counts$through[has_local_total] / local_total[has_local_total]

  # Keep the broader normalizations as diagnostics, but use the local
  # proportions for trajectory-style class assignment.
  pst <- pct_local_start
  pend <- pct_local_end
  pth <- pct_local_through

  horiz_dist_km <- fields$hcv_mag_km_yr * class_tyr
  vert_dist_km <- abs(fields$vcv_m_yr) * class_tyr / 1000
  dist_3d_km <- sqrt(horiz_dist_km^2 + vert_dist_km^2)

  boundary_mask <- array(FALSE, dim = dim(fields$valid_mask))
  for (k in seq_len(fields$ndepth)) {
    valid_k <- fields$valid_mask[, , k]
    nr <- nrow(valid_k)
    nc <- ncol(valid_k)
    bnd_k <- matrix(FALSE, nrow = nr, ncol = nc)
    valid_idx <- which(valid_k, arr.ind = TRUE)
    if (nrow(valid_idx) > 0L) {
      for (ii in seq_len(nrow(valid_idx))) {
        r <- valid_idx[ii, 1]
        c <- valid_idx[ii, 2]
        if (r == 1L || r == nr || c == 1L || c == nc) {
          bnd_k[r, c] <- TRUE
          next
        }
        nbrs <- c(
          valid_k[r - 1L, c],
          valid_k[r + 1L, c],
          valid_k[r, c - 1L],
          valid_k[r, c + 1L]
        )
        if (any(!nbrs)) {
          bnd_k[r, c] <- TRUE
        }
      }
    }
    boundary_mask[, , k] <- bnd_k
  }

  class_id <- array(NA_integer_, dim = c(fields$nrow, fields$ncol, fields$ndepth))
  class_label <- array(NA_character_, dim = c(fields$nrow, fields$ncol, fields$ndepth))

  movement_non <- dist_3d_km <= nmL
  movement_slow <- dist_3d_km > nmL & dist_3d_km <= smL
  sink_dominant <- (pend - pst) > sink_balance_tol
  is_sink <- pend >= Nend & counts$end >= min_count & sink_dominant

  # Tally logic:
  # - `counts$start` is how many seeded trajectories began in a voxel
  # - `counts$end` is how many finished in a voxel
  # - `counts$through` is how many passed through a voxel during stepping
  # The classification rules are then applied to these tallies in sequence,
  # rather than to the raw path geometry directly.
  class_id[movement_non & fields$valid_mask] <- 1L
  class_label[movement_non & fields$valid_mask] <- "Non-moving"

  class_id[is.na(class_id) & movement_slow & fields$valid_mask] <- 2L
  class_label[is.na(class_label) & movement_slow & fields$valid_mask] <- "Slow-moving"

  sink_idx <- is_sink & fields$valid_mask
  internal_sink_idx <- sink_idx & !boundary_mask
  boundary_sink_idx <- sink_idx & boundary_mask
  class_id[is.na(class_id) & internal_sink_idx] <- 3L
  class_label[is.na(class_label) & internal_sink_idx] <- "Internal sinks"
  class_id[is.na(class_id) & boundary_sink_idx] <- 4L
  class_label[is.na(class_label) & boundary_sink_idx] <- "Boundary sinks"

  remaining_idx <- which(fields$valid_mask & is.na(class_id), arr.ind = TRUE)
  if (nrow(remaining_idx) > 0L) {
    rem_start <- pst[remaining_idx]
    rem_end <- pend[remaining_idx]
    rem_through <- pth[remaining_idx]
    rem_balanced <- abs(rem_start - rem_end) <= corridor_balance_tol
    rem_class <- ifelse(
      rem_end == 0, 5L,
      ifelse(
        rem_end > Nend & rem_start < Nst, 6L,
        ifelse(
          rem_through > NFT & rem_balanced, 7L,
          ifelse(rem_end < rem_start, 8L, 9L)
        )
      )
    )
    rem_label <- c(
      "5" = "Sources",
      "6" = "Relative sinks",
      "7" = "Corridors",
      "8" = "Divergence",
      "9" = "Convergence"
    )[as.character(rem_class)]
    class_id[remaining_idx] <- rem_class
    class_label[remaining_idx] <- unname(rem_label)
  }

  class_id[fields$valid_mask & is.na(class_id)] <- 0L
  class_label[fields$valid_mask & is.na(class_label)] <- "Unclassified"

  lookup <- data.frame(
    class_id = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
    class_label = c(
      "Unclassified",
      "Non-moving",
      "Slow-moving",
      "Internal sinks",
      "Boundary sinks",
      "Sources",
      "Relative sinks",
      "Corridors",
      "Divergence",
      "Convergence"
    ),
    stringsAsFactors = FALSE
  )

  list(
    pct_total = list(start = pct_total_start, end = pct_total_end, through = pct_total_through),
    pct_depth = list(start = pct_depth_start, end = pct_depth_end, through = pct_depth_through),
    pct_local = list(start = pct_local_start, end = pct_local_end, through = pct_local_through),
    class_id = class_id,
    class_label = class_label,
    class_lookup = lookup,
    movement_distance_km = dist_3d_km,
    boundary_mask = boundary_mask
  )
}

#' Aggregate matched 3D classes by depth for profiles
#'
#' @param result Completed matched 2D/3D result containing runtime fields and
#'   3D classification arrays.
#' @returns Data frame of depth, class ID/label, valid-cell count, and percent
#'   of classified cells at each depth. Invalid and unclassified cells are not
#'   counted in the percentage denominator.
aggregate_class_profile_3d_fast <- function(result) {
  class_lookup <- class_lookup_named_3d_fast(result$classification$class_lookup)
  valid_mask <- result$fields_runtime$valid_mask
  class_id <- result$classification$class_id

  out <- dplyr::bind_rows(lapply(seq_along(result$fields_runtime$depths), function(k) {
    vals <- class_id[, , k]
    vals <- vals[valid_mask[, , k] & !is.na(vals)]
    tb <- table(factor(vals, levels = as.integer(names(class_lookup))))
    total <- sum(tb)
    data.frame(
      depth_m = result$fields_runtime$depths[k],
      class_id = as.integer(names(tb)),
      class_label = unname(class_lookup[names(tb)]),
      n_cells = as.integer(tb),
      pct_cells = if (total > 0) 100 * as.integer(tb) / total else 0,
      stringsAsFactors = FALSE
    )
  }))

  out
}

#' Exclude specified depths from matched-workflow tables
#'
#' @param df Data frame that may contain a `depth_m` column.
#' @param exclude_depths Numeric depths in metres to omit; by default the
#'   non-manuscript 2.5 m layer.
#' @returns Input data frame without excluded depth rows, or unchanged when
#'   `depth_m` is absent.
exclude_depths_3d_fast <- function(df, exclude_depths = c(2.5)) {
  if (is.null(df$depth_m)) {
    return(df)
  }
  df[!(as.numeric(df$depth_m) %in% exclude_depths), , drop = FALSE]
}

array_to_raster_stack_3d_fast <- function(arr, template, depths, prefix) {
  layers <- vector("list", length(depths))
  for (k in seq_along(depths)) {
    r <- template
    # Arrays are stored as [row, col, depth]; transpose on write-back so the
    # raster footprint matches the original terra cell ordering.
    terra::values(r) <- as.vector(t(arr[, , k]))
    names(r) <- paste0(prefix, "_", sprintf("%05.1f", depths[k]), "m")
    layers[[k]] <- r
  }
  terra::rast(layers)
}

matrix_to_raster_3d_fast <- function(mat, template, name) {
  if (!all(dim(mat) == c(terra::nrow(template), terra::ncol(template)))) {
    stop("Matrix dimensions do not match template raster rows/columns.", call. = FALSE)
  }
  r <- template
  # Matrices are stored as [row, col]; transpose on write-back so terra cells
  # receive row-major values matching the map footprint.
  terra::values(r) <- as.vector(t(mat))
  names(r) <- name
  r
}

project_class_2d_3d_fast <- function(class_id, valid_mask) {
  nr <- dim(class_id)[1]
  nc <- dim(class_id)[2]
  out <- matrix(NA_integer_, nrow = nr, ncol = nc)
  priority <- c(9L, 8L, 7L, 6L, 5L, 4L, 3L, 2L, 1L, 0L)

  # Depth collapse:
  # after each voxel in the 3D stack has its own class, a single surface
  # summary can be created by choosing the highest-priority class found in each
  # vertical water column.
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      vals <- class_id[i, j, ]
      vals <- vals[!is.na(vals) & valid_mask[i, j, ]]
      if (length(vals) == 0L) {
        next
      }
      hit <- priority[priority %in% vals]
      out[i, j] <- if (length(hit) > 0L) hit[1] else vals[1]
    }
  }

  out
}

diagnostics_summary_3d_fast <- function(run_obj) {
  diag_df <- run_obj$diagnostics
  same_end_voxel <- with(
    diag_df,
    start_row == end_row & start_col == end_col & start_k == end_k
  )
  returned_to_start_voxel <- with(
    diag_df,
    moved_from_start & same_end_voxel
  )
  data.frame(
    n_seeds = nrow(diag_df),
    pct_moved_from_start = mean(diag_df$moved_from_start) * 100,
    pct_same_end_voxel = mean(same_end_voxel) * 100,
    pct_returned_to_start_voxel = mean(returned_to_start_voxel) * 100,
    mean_horizontal_distance_km = mean(diag_df$horizontal_distance_km, na.rm = TRUE),
    median_horizontal_distance_km = stats::median(diag_df$horizontal_distance_km, na.rm = TRUE),
    mean_vertical_distance_m = mean(diag_df$vertical_distance_m, na.rm = TRUE),
    median_vertical_distance_m = stats::median(diag_df$vertical_distance_m, na.rm = TRUE),
    pct_active_at_end = mean(diag_df$active_at_end) * 100
  )
}

voxel_class_diagnostics_3d_fast <- function(result) {
  fields <- result$fields_runtime
  class_lookup <- class_lookup_named_3d_fast(result$classification$class_lookup)
  grid <- expand.grid(
    row = seq_len(fields$nrow),
    col = seq_len(fields$ncol),
    k = seq_len(fields$ndepth)
  )
  keep <- fields$valid_mask[cbind(grid$row, grid$col, grid$k)]
  grid <- grid[keep, , drop = FALSE]
  grid$x <- fields$x_centers[grid$col]
  grid$y <- fields$y_centers[grid$row]
  grid$depth_m <- fields$depths[grid$k]
  grid$class_id <- result$classification$class_id[cbind(grid$row, grid$col, grid$k)]
  grid$class_label <- unname(class_lookup[as.character(grid$class_id)])
  grid$start_count <- result$counts$start[cbind(grid$row, grid$col, grid$k)]
  grid$end_count <- result$counts$end[cbind(grid$row, grid$col, grid$k)]
  grid$through_count <- result$counts$through[cbind(grid$row, grid$col, grid$k)]
  grid$pct_start_local <- result$classification$pct_local$start[cbind(grid$row, grid$col, grid$k)]
  grid$pct_end_local <- result$classification$pct_local$end[cbind(grid$row, grid$col, grid$k)]
  grid$pct_through_local <- result$classification$pct_local$through[cbind(grid$row, grid$col, grid$k)]
  grid$start_minus_end_count <- grid$start_count - grid$end_count
  grid$end_minus_start_count <- grid$end_count - grid$start_count
  grid$start_minus_end_pct <- grid$pct_start_local - grid$pct_end_local
  grid$abs_start_end_pct <- abs(grid$start_minus_end_pct)
  grid$terminal_pct <- grid$pct_start_local + grid$pct_end_local
  grid$through_to_terminal_ratio <- ifelse(
    grid$terminal_pct > 0,
    grid$pct_through_local / grid$terminal_pct,
    NA_real_
  )
  grid$weak_imbalance_le_1pct <- grid$abs_start_end_pct <= 1
  grid$weak_imbalance_le_5pct <- grid$abs_start_end_pct <= 5
  grid$weak_imbalance_le_10pct <- grid$abs_start_end_pct <= 10
  rownames(grid) <- NULL
  grid
}

corridor_divergence_summary_3d_fast <- function(result) {
  voxel_df <- voxel_class_diagnostics_3d_fast(result)
  focus <- voxel_df[voxel_df$class_label %in% c("Corridors", "Divergence"), , drop = FALSE]
  if (nrow(focus) == 0L) {
    return(data.frame())
  }
  dplyr::group_by(focus, depth_m, class_id, class_label) %>%
    dplyr::summarise(
      n_voxels = dplyr::n(),
      mean_start_count = mean(start_count, na.rm = TRUE),
      mean_end_count = mean(end_count, na.rm = TRUE),
      mean_through_count = mean(through_count, na.rm = TRUE),
      mean_pct_start_local = mean(pct_start_local, na.rm = TRUE),
      mean_pct_end_local = mean(pct_end_local, na.rm = TRUE),
      mean_pct_through_local = mean(pct_through_local, na.rm = TRUE),
      mean_start_minus_end_count = mean(start_minus_end_count, na.rm = TRUE),
      mean_start_minus_end_pct = mean(start_minus_end_pct, na.rm = TRUE),
      mean_abs_start_end_pct = mean(abs_start_end_pct, na.rm = TRUE),
      median_abs_start_end_pct = stats::median(abs_start_end_pct, na.rm = TRUE),
      prop_weak_imbalance_le_1pct = mean(weak_imbalance_le_1pct, na.rm = TRUE),
      prop_weak_imbalance_le_5pct = mean(weak_imbalance_le_5pct, na.rm = TRUE),
      prop_weak_imbalance_le_10pct = mean(weak_imbalance_le_10pct, na.rm = TRUE),
      mean_through_to_terminal_ratio = mean(through_to_terminal_ratio, na.rm = TRUE),
      .groups = "drop"
    )
}

class_balance_summary_3d_fast <- function(result) {
  voxel_df <- voxel_class_diagnostics_3d_fast(result)
  dplyr::group_by(voxel_df, depth_m, class_id, class_label) %>%
    dplyr::summarise(
      n_voxels = dplyr::n(),
      mean_start_count = mean(start_count, na.rm = TRUE),
      mean_end_count = mean(end_count, na.rm = TRUE),
      mean_through_count = mean(through_count, na.rm = TRUE),
      mean_pct_start_local = mean(pct_start_local, na.rm = TRUE),
      mean_pct_end_local = mean(pct_end_local, na.rm = TRUE),
      mean_pct_through_local = mean(pct_through_local, na.rm = TRUE),
      mean_start_minus_end_count = mean(start_minus_end_count, na.rm = TRUE),
      mean_start_minus_end_pct = mean(start_minus_end_pct, na.rm = TRUE),
      mean_abs_start_end_pct = mean(abs_start_end_pct, na.rm = TRUE),
      prop_weak_imbalance_le_5pct = mean(weak_imbalance_le_5pct, na.rm = TRUE),
      prop_weak_imbalance_le_10pct = mean(weak_imbalance_le_10pct, na.rm = TRUE),
      mean_through_to_terminal_ratio = mean(through_to_terminal_ratio, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Create a named matched-class lookup vector
#'
#' @param class_lookup Data frame with integer `class_id` and character
#'   `class_label` columns.
#' @returns Character vector of labels named by class ID, for plot/table
#'   labelling. Missing required columns fail naturally.
class_lookup_named_3d_fast <- function(class_lookup) {
  out <- class_lookup$class_label
  names(out) <- as.character(class_lookup$class_id)
  out
}

#' Return the ordered matched-class labels
#'
#' @param include_unclassified Logical; include the zero-valued
#'   `Unclassified` display level when `TRUE`.
#' @returns Ordered character vector used by the matched-workflow plots and
#'   tables.
class_levels_3d_fast <- function(include_unclassified = TRUE) {
  levs <- c(
    "Non-moving",
    "Slow-moving",
    "Internal sinks",
    "Boundary sinks",
    "Sources",
    "Relative sinks",
    "Corridors",
    "Divergence",
    "Convergence"
  )
  if (isTRUE(include_unclassified)) {
    c(levs, "Unclassified")
  } else {
    levs
  }
}

#' Return the matched-class display palette
#'
#' @returns Named character vector of hexadecimal colours keyed by the ordered
#'   matched-class labels, including `Unclassified`.
class_colors_3d_fast <- function() {
  c(
    "Non-moving" = "#000000",
    "Slow-moving" = "#E69F00",
    "Internal sinks" = "#56B4E9",
    "Boundary sinks" = "#009E73",
    "Sources" = "#F0E442",
    "Relative sinks" = "#0072B2",
    "Corridors" = "#D55E00",
    "Divergence" = "#CC79A7",
    "Convergence" = "#999999",
    "Unclassified" = "#d9d9d9"
  )
}

#' Select colours for class labels present in an interactive plot
#'
#' @param labels Character or factor vector of class labels represented in a
#'   plot.
#' @param include_unclassified Logical; allow the `Unclassified` colour.
#' @returns Unnamed hexadecimal colour vector in the fixed class order for
#'   labels present in `labels`.
class_color_scale_3d_fast <- function(labels, include_unclassified = TRUE) {
  palette <- class_colors_3d_fast()
  levels_all <- class_levels_3d_fast(include_unclassified = include_unclassified)
  present <- levels_all[levels_all %in% unique(as.character(stats::na.omit(labels)))]
  unname(palette[present])
}

#' Select the fixed palette for static matched-workflow plots
#'
#' @param include_unclassified Logical; include the `Unclassified` entry.
#' @returns Named hexadecimal colour vector in the fixed class order.
class_palette_static_3d_fast <- function(include_unclassified = TRUE) {
  palette <- class_colors_3d_fast()
  levels_all <- class_levels_3d_fast(include_unclassified = include_unclassified)
  palette[levels_all]
}

#' Attach matched class labels and colours to a raster
#'
#' @param r Raster of matched integer class IDs.
#' @param class_lookup Data frame with matched `class_id` and `class_label`
#'   columns, normally from `result$classification$class_lookup`.
#' @returns Factor `SpatRaster` with levels and a colour table attached. Stops
#'   when the raster or lookup is incompatible with terra.
decorate_class_raster_3d_fast <- function(r, class_lookup) {
  r <- terra::as.factor(r)
  levels(r) <- data.frame(
    value = class_lookup$class_id,
    class_label = class_lookup$class_label,
    stringsAsFactors = FALSE
  )
  terra::coltab(r) <- data.frame(
    value = class_lookup$class_id,
    col = unname(class_colors_3d_fast()[class_lookup$class_label]),
    stringsAsFactors = FALSE
  )
  r
}

#' Extract one matched-workflow class profile through depth
#'
#' @param result Completed matched 2D/3D result with classification and tally
#'   arrays.
#' @param x Longitude in the result raster coordinate reference system.
#' @param y Latitude in the result raster coordinate reference system.
#' @returns Data frame of class, start/end/through counts, and local/depth
#'   percentages for one water column. Stops when x/y lies outside the grid.
class_profile_for_pixel_3d_fast <- function(result, x, y) {
  rc <- row_col_from_xy_3d_fast(x, y, result$fields_runtime)
  if (rc$row < 1L || rc$row > result$fields_runtime$nrow || rc$col < 1L || rc$col > result$fields_runtime$ncol) {
    stop("Requested x/y falls outside the raster extent.", call. = FALSE)
  }

  ids <- result$classification$class_id[rc$row, rc$col, ]
  labels <- class_lookup_named_3d_fast(result$classification$class_lookup)[as.character(ids)]

  data.frame(
    depth_m = result$fields_runtime$depths,
    class_id = ids,
    class_label = unname(labels),
    start_count = result$counts$start[rc$row, rc$col, ],
    end_count = result$counts$end[rc$row, rc$col, ],
    through_count = result$counts$through[rc$row, rc$col, ],
    pct_start_local = result$classification$pct_local$start[rc$row, rc$col, ],
    pct_end_local = result$classification$pct_local$end[rc$row, rc$col, ],
    pct_through_local = result$classification$pct_local$through[rc$row, rc$col, ],
    pct_start = result$classification$pct_depth$start[rc$row, rc$col, ],
    pct_end = result$classification$pct_depth$end[rc$row, rc$col, ],
    pct_through = result$classification$pct_depth$through[rc$row, rc$col, ],
    stringsAsFactors = FALSE
  )
}

#' Draw interactive matched particle seeds and sampled paths
#'
#' @param result Completed matched 2D/3D result with seeds, diagnostics, and
#'   optional sampled paths.
#' @param color_by Either `"depth"` or `"moved"` for seed-marker colour.
#' @param max_paths Maximum sampled paths/seeds displayed.
#' @param path_line_width Positive Plotly line width for retained paths.
#' @returns Interactive Plotly object. Requires `plotly`; stops when no seeds
#'   are available or the package is missing. It does not write a file.
plot_seeded_paths_plotly_3d_fast <- function(result,
                                             color_by = c("depth", "moved"),
                                             max_paths = 200L,
                                             path_line_width = 5) {
  color_by <- match.arg(color_by)
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive plots.", call. = FALSE)
  }

  seeds_df <- result$seeds
  if (nrow(seeds_df) == 0L) {
    stop("No seeds available in result.", call. = FALSE)
  }

  depth_cols <- grDevices::rgb(
    t(grDevices::col2rgb(viridisLite::viridis(max(2L, length(unique(seeds_df$z))))))
      / 255
  )

  if (length(result$sample_paths) > 0L) {
    show_ids <- as.integer(head(names(result$sample_paths), max_paths))
  } else {
    show_ids <- head(seeds_df$seed_id, max_paths)
  }
  seeds_df <- seeds_df[seeds_df$seed_id %in% show_ids, , drop = FALSE]
  diag_df <- result$diagnostics$per_seed
  diag_df <- diag_df[diag_df$seed_id %in% show_ids, , drop = FALSE]
  seeds_df <- dplyr::left_join(seeds_df, diag_df[, c("seed_id", "moved_from_start")], by = "seed_id")

  p <- plotly::plot_ly()

  p <- plotly::add_markers(
    p,
    data = seeds_df,
    x = ~x,
    y = ~y,
    z = ~z,
    color = if (color_by == "depth") ~z else ~moved_from_start,
    colors = if (color_by == "depth") depth_cols else c("#9e9e9e", "#d95f02"),
    marker = list(size = 6, opacity = 0.85),
    hovertemplate = paste(
      "seed_id=%{customdata[1]}",
      "<br>depth=%{customdata[2]} m",
      "<br>row=%{customdata[3]} col=%{customdata[4]}",
      "<br>moved=%{customdata[5]}<extra></extra>"
    ),
    customdata = cbind(
      seeds_df$seed_id,
      seeds_df$z,
      seeds_df$row,
      seeds_df$col,
      seeds_df$moved_from_start
    ),
    name = "Seeds"
  )

  if (length(result$sample_paths) > 0L) {
    for (sid in names(result$sample_paths)) {
      path_df <- result$sample_paths[[sid]]
      if (nrow(path_df) < 2L || !(as.integer(sid) %in% show_ids)) {
        next
      }
      p <- plotly::add_trace(
        p,
        type = "scatter3d",
        mode = "lines",
        x = path_df$x,
        y = path_df$y,
        z = path_df$z,
        line = list(width = path_line_width, color = "rgba(33,113,181,0.55)"),
        hoverinfo = "skip",
        showlegend = FALSE
      )
    }
  }

  plotly::layout(
    p,
    legend = list(
      font = list(size = 16),
      title = list(font = list(size = 18)),
      itemsizing = "constant"
    ),
    scene = list(
      xaxis = list(title = "Longitude"),
      yaxis = list(title = "Latitude"),
      zaxis = list(title = "Depth (m)", autorange = "reversed")
    )
  )
}

#' Draw an interactive matched 3D class or velocity view
#'
#' @param result Completed matched 2D/3D result with runtime fields and
#'   classifications.
#' @param metric One of `"class"`, `"hcv"`, `"vcv"`, or `"movement"`.
#' @param max_points Maximum valid voxels plotted; excess points are randomly
#'   sampled for display only.
#' @returns Interactive Plotly object. Requires `plotly`; it does not alter the
#'   result or write a file.
plot_cv_layers_plotly_3d_fast <- function(result,
                                          metric = c("class", "hcv", "vcv", "movement"),
                                          max_points = 25000L) {
  metric <- match.arg(metric)
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive plots.", call. = FALSE)
  }

  fields <- result$fields_runtime
  tyr_years <- if (!is.null(result$parameters) && !is.null(result$parameters$tyr)) result$parameters$tyr else "projection"
  grid <- expand.grid(
    row = seq_len(fields$nrow),
    col = seq_len(fields$ncol),
    k = seq_len(fields$ndepth)
  )
  keep <- fields$valid_mask[cbind(grid$row, grid$col, grid$k)]
  grid <- grid[keep, , drop = FALSE]
  if (nrow(grid) > max_points) {
    grid <- grid[sample.int(nrow(grid), max_points), , drop = FALSE]
  }

  grid$x <- fields$x_centers[grid$col]
  grid$y <- fields$y_centers[grid$row]
  grid$depth_m <- fields$depths[grid$k]

  if (metric == "class") {
    grid$value <- result$classification$class_id[cbind(grid$row, grid$col, grid$k)]
    grid$label <- class_lookup_named_3d_fast(result$classification$class_lookup)[as.character(grid$value)]
    grid$label <- factor(grid$label, levels = class_levels_3d_fast())
    color_arg <- ~label
    color_scale <- class_color_scale_3d_fast(grid$label)
  } else if (metric == "hcv") {
    grid$value <- fields$hcv_mag_km_yr[cbind(grid$row, grid$col, grid$k)]
    grid$label <- sprintf("%.3f km/yr", grid$value)
    color_arg <- ~value
    color_scale <- viridisLite::viridis(64)
  } else if (metric == "vcv") {
    grid$value <- fields$vcv_m_yr[cbind(grid$row, grid$col, grid$k)]
    grid$label <- sprintf("%.3f m/yr", grid$value)
    color_arg <- ~value
    color_scale <- viridisLite::viridis(64)
  } else {
    grid$value <- result$classification$movement_distance_km[cbind(grid$row, grid$col, grid$k)]
    grid$label <- sprintf("%.3f km over %s years", grid$value, tyr_years)
    color_arg <- ~value
    color_scale <- viridisLite::viridis(64)
  }

  grid$hover_text <- paste0(
    "lon=", signif(grid$x, 6),
    "<br>lat=", signif(grid$y, 6),
    "<br>depth=", signif(grid$depth_m, 6), " m",
    "<br>value=", grid$label
  )

  plotly::plot_ly(
    data = grid,
    x = ~x,
    y = ~y,
    z = ~depth_m,
    type = "scatter3d",
    mode = "markers",
    color = color_arg,
    colors = color_scale,
    marker = list(size = 5, opacity = 0.85),
    text = ~hover_text,
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      legend = list(
        font = list(size = 16),
        title = list(font = list(size = 18)),
        itemsizing = "constant"
      ),
      scene = list(
        xaxis = list(title = "Longitude"),
        yaxis = list(title = "Latitude"),
        zaxis = list(title = "Depth (m)", autorange = "reversed")
      )
    )
}

#' Draw an interactive matched class strip for one water column
#'
#' @param result Completed matched 2D/3D result.
#' @param x Longitude in the result raster coordinate reference system.
#' @param y Latitude in the result raster coordinate reference system.
#' @returns Interactive Plotly bar chart of classes through depth. Requires
#'   `plotly` and stops when x/y falls outside the raster.
plot_class_profile_plotly_3d_fast <- function(result, x, y) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive plots.", call. = FALSE)
  }

  profile_df <- class_profile_for_pixel_3d_fast(result, x = x, y = y)
  profile_df$class_label <- factor(profile_df$class_label, levels = class_levels_3d_fast())
  profile_df$bar_value <- 1
  profile_df$hover_text <- paste0(
    "depth=", profile_df$depth_m, " m",
    "<br>class=", profile_df$class_label,
    "<br>start=", signif(profile_df$pct_start, 6), "%",
    "<br>end=", signif(profile_df$pct_end, 6), "%",
    "<br>through=", signif(profile_df$pct_through, 6), "%",
    "<br>start_count=", profile_df$start_count,
    "<br>end_count=", profile_df$end_count,
    "<br>through_count=", profile_df$through_count
  )

  plotly::plot_ly(
    profile_df,
    x = ~bar_value,
    y = ~factor(depth_m, levels = rev(unique(depth_m))),
    color = ~class_label,
    type = "bar",
    orientation = "h",
    colors = class_color_scale_3d_fast(profile_df$class_label),
    text = ~hover_text,
    hoverinfo = "text",
    marker = list(line = list(width = 0))
  ) %>%
    plotly::layout(
      barmode = "overlay",
      legend = list(
        font = list(size = 16),
        title = list(font = list(size = 18)),
        itemsizing = "constant"
      ),
      xaxis = list(title = "Class strip", showticklabels = FALSE),
      yaxis = list(title = "Depth (m)", autorange = "reversed")
    )
}

#' Draw an interactive all-domain matched class profile
#'
#' @param result Completed matched 2D/3D result.
#' @returns Interactive stacked Plotly profile after excluding configured
#'   display depths. Requires `plotly`; it does not write a file.
plot_class_profile_all_plotly_3d_fast <- function(result) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive plots.", call. = FALSE)
  }

  prof <- aggregate_class_profile_3d_fast(result)
  prof <- exclude_depths_3d_fast(prof)
  prof <- prof[prof$n_cells > 0, , drop = FALSE]
  prof$class_label <- factor(prof$class_label, levels = class_levels_3d_fast(FALSE))
  prof$hover_text <- paste0(
    "depth=", prof$depth_m, " m",
    "<br>class=", prof$class_label,
    "<br>cells=", prof$n_cells,
    "<br>percent=", signif(prof$pct_cells, 6), "%"
  )

  plotly::plot_ly(
    prof,
    x = ~pct_cells,
    y = ~factor(depth_m, levels = rev(unique(depth_m))),
    color = ~class_label,
    type = "bar",
    orientation = "h",
    colors = class_color_scale_3d_fast(prof$class_label, include_unclassified = FALSE),
    text = ~hover_text,
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      barmode = "stack",
      legend = list(
        font = list(size = 16),
        title = list(font = list(size = 18)),
        itemsizing = "constant"
      ),
      xaxis = list(title = "Percent of classified cells at depth"),
      yaxis = list(title = "Depth (m)", autorange = "reversed")
    )
}

#' Draw a static all-domain matched class profile
#'
#' @param result Completed matched 2D/3D result.
#' @returns Stacked horizontal-bar `ggplot` object after excluding configured
#'   display depths. Requires `ggplot2`; it does not write a file.
plot_class_profile_all_ggplot_3d_fast <- function(result) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for static plots.", call. = FALSE)
  }

  prof <- aggregate_class_profile_3d_fast(result)
  prof <- exclude_depths_3d_fast(prof)
  prof <- prof[prof$n_cells > 0, , drop = FALSE]
  prof$class_label <- factor(prof$class_label, levels = class_levels_3d_fast(FALSE))

  ggplot2::ggplot(
    prof,
    ggplot2::aes(
      x = pct_cells,
      y = factor(depth_m, levels = rev(sort(unique(depth_m)))),
      fill = class_label
    )
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_manual(
      values = class_palette_static_3d_fast(include_unclassified = FALSE),
      drop = TRUE
    ) +
    ggplot2::labs(
      title = "3D trajectory classes by depth",
      x = "Percent of classified cells at depth",
      y = "Depth (m)",
      fill = "Class"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

#' Tabulate matched 3D class counts by depth and overall
#'
#' @param result Completed matched 2D/3D result with class arrays and valid
#'   mask.
#' @returns List with `by_depth` and `overall` class-count data frames. Only
#'   valid, classified cells are counted.
class_tables_3d_fast <- function(result) {
  class_lookup <- class_lookup_named_3d_fast(result$classification$class_lookup)
  class_id <- result$classification$class_id
  valid_mask <- result$fields_runtime$valid_mask

  depth_table <- dplyr::bind_rows(lapply(seq_along(result$fields_runtime$depths), function(k) {
    vals <- class_id[, , k]
    vals <- vals[valid_mask[, , k] & !is.na(vals)]
    tb <- table(factor(vals, levels = as.integer(names(class_lookup))))
    data.frame(
      depth_m = result$fields_runtime$depths[k],
      class_id = as.integer(names(tb)),
      class_label = unname(class_lookup[names(tb)]),
      n_cells = as.integer(tb),
      stringsAsFactors = FALSE
    )
  }))

  overall_vals <- class_id[valid_mask & !is.na(class_id)]
  overall_tb <- table(factor(overall_vals, levels = as.integer(names(class_lookup))))
  overall_table <- data.frame(
    class_id = as.integer(names(overall_tb)),
    class_label = unname(class_lookup[names(overall_tb)]),
    n_cells = as.integer(overall_tb),
    stringsAsFactors = FALSE
  )

  list(
    by_depth = depth_table,
    overall = overall_table
  )
}

#' Tabulate the projected matched 2D class counts
#'
#' @param result Completed matched 2D/3D result with `classification$class_2d`.
#' @returns Data frame of class IDs, labels, and non-missing projected-2D cell
#'   counts.
class_table_2d_fast <- function(result) {
  class_lookup <- class_lookup_named_3d_fast(result$classification$class_lookup)
  vals <- as.vector(result$classification$class_2d)
  vals <- vals[!is.na(vals)]
  tb <- table(factor(vals, levels = as.integer(names(class_lookup))))
  data.frame(
    class_id = as.integer(names(tb)),
    class_label = unname(class_lookup[names(tb)]),
    n_cells = as.integer(tb),
    stringsAsFactors = FALSE
  )
}

#' Draw an interactive projected matched-2D class map
#'
#' @param result Completed matched 2D/3D result with a projected 2D class map.
#' @returns Interactive Plotly map. Requires `plotly`; it omits missing class
#'   cells and does not write a file.
plot_class_map_plotly_2d_fast <- function(result) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive plots.", call. = FALSE)
  }

  fields <- result$fields_runtime
  grid <- expand.grid(
    row = seq_len(fields$nrow),
    col = seq_len(fields$ncol)
  )
  grid$x <- fields$x_centers[grid$col]
  grid$y <- fields$y_centers[grid$row]
  grid$class_id <- as.vector(result$classification$class_2d)
  lut <- class_lookup_named_3d_fast(result$classification$class_lookup)
  grid$class_label <- unname(lut[as.character(grid$class_id)])
  grid <- grid[!is.na(grid$class_id), , drop = FALSE]
  grid$class_label <- factor(grid$class_label, levels = class_levels_3d_fast())
  grid$hover_text <- paste0(
    "lon=", signif(grid$x, 6),
    "<br>lat=", signif(grid$y, 6),
    "<br>class=", grid$class_label
  )

  plotly::plot_ly(
    data = grid,
    x = ~x,
    y = ~y,
    type = "scatter",
    mode = "markers",
    color = ~class_label,
    colors = class_color_scale_3d_fast(grid$class_label),
    marker = list(size = 8, opacity = 0.9, symbol = "square"),
    text = ~hover_text,
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      legend = list(
        font = list(size = 16),
        title = list(font = list(size = 18)),
        itemsizing = "constant"
      ),
      xaxis = list(title = "Longitude"),
      yaxis = list(title = "Latitude", scaleanchor = "x", scaleratio = 1)
    )
}

#' Draw a static projected matched-2D class map
#'
#' @param result Completed matched 2D/3D result with a projected 2D class map.
#' @returns Static `ggplot` raster map. Requires `ggplot2`; it omits missing
#'   class cells and does not write a file.
plot_class_map_ggplot_2d_fast <- function(result) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for static plots.", call. = FALSE)
  }

  fields <- result$fields_runtime
  grid <- expand.grid(
    row = seq_len(fields$nrow),
    col = seq_len(fields$ncol)
  )
  grid$x <- fields$x_centers[grid$col]
  grid$y <- fields$y_centers[grid$row]
  grid$class_id <- as.vector(result$classification$class_2d)
  lut <- class_lookup_named_3d_fast(result$classification$class_lookup)
  grid$class_label <- unname(lut[as.character(grid$class_id)])
  grid <- grid[!is.na(grid$class_id), , drop = FALSE]
  grid$class_label <- factor(grid$class_label, levels = class_levels_3d_fast())

  ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y, fill = class_label)) +
    ggplot2::geom_raster() +
    ggplot2::coord_quickmap() +
    ggplot2::scale_fill_manual(
      values = class_palette_static_3d_fast(),
      drop = TRUE,
      na.value = "white"
    ) +
    ggplot2::labs(
      title = "Projected 2D class map",
      x = "Longitude",
      y = "Latitude",
      fill = "Class"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

#' Partition matched depths for faceted map exports
#'
#' @param depths Numeric vector of depth values in metres.
#' @param n_per_group Positive integer number of unique depths per group.
#' @returns Named list of sorted, unique depth vectors. Empty input returns an
#'   empty list; callers should validate any intended non-empty export.
split_depth_groups_3d_fast <- function(depths, n_per_group = 9L) {
  depth_vals <- sort(unique(as.numeric(depths)))
  split(depth_vals, ceiling(seq_along(depth_vals) / n_per_group))
}

#' Draw static matched 3D class maps by depth
#'
#' @param result Completed matched 2D/3D result with runtime fields and 3D
#'   class IDs.
#' @param depths_subset Optional numeric depths in metres to display; `NULL`
#'   displays all available depths.
#' @returns Faceted `ggplot` raster map. Requires `ggplot2`; an unmatched depth
#'   selection yields an empty plot rather than generating an output file.
plot_class_depth_maps_ggplot_3d_fast <- function(result, depths_subset = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for static plots.", call. = FALSE)
  }

  fields <- result$fields_runtime
  depth_vals <- fields$depths
  if (!is.null(depths_subset)) {
    keep_k <- which(depth_vals %in% depths_subset)
  } else {
    keep_k <- seq_along(depth_vals)
  }

  grid <- expand.grid(
    row = seq_len(fields$nrow),
    col = seq_len(fields$ncol),
    k = keep_k
  )
  keep <- fields$valid_mask[cbind(grid$row, grid$col, grid$k)]
  grid <- grid[keep, , drop = FALSE]
  grid$x <- fields$x_centers[grid$col]
  grid$y <- fields$y_centers[grid$row]
  grid$depth_m <- depth_vals[grid$k]
  grid$depth_f <- factor(
    paste0(signif(grid$depth_m, 6), " m"),
    levels = paste0(signif(sort(unique(as.numeric(grid$depth_m))), 6), " m")
  )
  grid$class_id <- result$classification$class_id[cbind(grid$row, grid$col, grid$k)]
  lut <- class_lookup_named_3d_fast(result$classification$class_lookup)
  grid$class_label <- unname(lut[as.character(grid$class_id)])
  grid$class_label <- factor(grid$class_label, levels = class_levels_3d_fast())

  ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y, fill = class_label)) +
    ggplot2::geom_raster() +
    ggplot2::coord_quickmap() +
    ggplot2::facet_wrap(~ depth_f, ncol = 3) +
    ggplot2::scale_fill_manual(
      values = class_palette_static_3d_fast(),
      drop = TRUE,
      na.value = "white"
    ) +
    ggplot2::labs(
      title = "3D trajectory classes by depth",
      x = "Longitude",
      y = "Latitude",
      fill = "Class"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

#' Draw static matched movement-distance maps by depth
#'
#' @param result Completed matched 2D/3D result with runtime fields and local
#'   velocity-derived movement distance.
#' @param depths_subset Optional numeric depths in metres to display; `NULL`
#'   displays all available depths.
#' @returns Faceted `ggplot` raster map of movement kilometres over the
#'   classification horizon. Requires `ggplot2`; it does not write a file.
plot_movement_depth_maps_ggplot_3d_fast <- function(result, depths_subset = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for static plots.", call. = FALSE)
  }

  fields <- result$fields_runtime
  depth_vals <- fields$depths
  if (!is.null(depths_subset)) {
    keep_k <- which(depth_vals %in% depths_subset)
  } else {
    keep_k <- seq_along(depth_vals)
  }

  grid <- expand.grid(
    row = seq_len(fields$nrow),
    col = seq_len(fields$ncol),
    k = keep_k
  )
  keep <- fields$valid_mask[cbind(grid$row, grid$col, grid$k)]
  grid <- grid[keep, , drop = FALSE]
  grid$x <- fields$x_centers[grid$col]
  grid$y <- fields$y_centers[grid$row]
  grid$depth_m <- depth_vals[grid$k]
  grid$depth_f <- factor(
    paste0(signif(grid$depth_m, 6), " m"),
    levels = paste0(signif(sort(unique(as.numeric(grid$depth_m))), 6), " m")
  )
  grid$movement_km <- result$classification$movement_distance_km[cbind(grid$row, grid$col, grid$k)]

  ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y, fill = movement_km)) +
    ggplot2::geom_raster() +
    ggplot2::coord_quickmap() +
    ggplot2::facet_wrap(~ depth_f, ncol = 3) +
    ggplot2::scale_fill_viridis_c(
      option = "magma",
      na.value = "white",
      name = "Movement\n(km over tyr)"
    ) +
    ggplot2::labs(
      title = "3D movement distance by depth",
      x = "Longitude",
      y = "Latitude"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

comparison_table_3d_fast <- function(result,
                                     method = "3D_HCV_VCV_trajectories",
                                     dataset = "BRAN",
                                     scenario = NA_character_) {
  prof <- aggregate_class_profile_3d_fast(result)
  prof <- exclude_depths_3d_fast(prof)
  prof$method <- method
  prof$dataset <- dataset
  prof$scenario <- scenario
  prof$time_period <- result$fields_runtime$time_period
  prof[, c("method", "dataset", "scenario", "time_period", "depth_m", "class_id", "class_label", "n_cells", "pct_cells")]
}

metadata_table_3d_fast <- function(result,
                                   method = "3D_HCV_VCV_trajectories",
                                   script_path = "01_scripts/610_functions_3D_fast.R",
                                   dataset = "BRAN",
                                   scenario = NA_character_) {
  data.frame(
    method = method,
    script_path = script_path,
    dataset = dataset,
    scenario = scenario,
    time_period = result$fields_runtime$time_period,
    movement_distance_units = result$assumptions$movement_distance_units,
    horizontal_distance_units = result$assumptions$horizontal_distance_units,
    vertical_distance_units = result$assumptions$vertical_distance_units,
    nmL_definition = result$assumptions$nmL_definition,
    smL_definition = result$assumptions$smL_definition,
    sink_balance_tol = result$parameters$sink_balance_tol %||% NA_real_,
    corridor_balance_tol = result$parameters$corridor_balance_tol %||% NA_real_,
    jitter_behavior = result$assumptions$jitter_behavior,
    classification_basis = result$assumptions$classification_basis,
    stringsAsFactors = FALSE
  )
}

###############################################################################
################ CONNECTIVITY REVISION PARTICLE ENGINE ########################
###############################################################################

# The functions below form the self-contained particle engine used by stages
# 05 and 06 of this demo.

#' Convert a revision velocity bundle to a particle-field list
#'
#' @param velocity_bundle List created by `build_revision_velocity_bundle()`.
#' @param vertical_mode Either `"observed"` to retain VVoCC or `"zero"` for
#'   the matched horizontal-only control.
#' @returns Plain list of arrays, geometry, depth edges, and valid mask used by
#'   the corrected particle engine. Stops for malformed arrays or geometry.
connectivity_fields_from_bundle <- function(velocity_bundle,
                                            vertical_mode = c("observed", "zero")) {
  vertical_mode <- match.arg(vertical_mode)
  required <- c("ext", "res", "x_centers", "y_centers", "depths_m", "hcv_km_yr", "hcv_angle_deg", "vcv_m_yr", "valid_mask")
  missing <- setdiff(required, names(velocity_bundle))
  if (length(missing)) stop("Velocity bundle is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  depths <- as.numeric(velocity_bundle$depths_m)
  dims <- dim(velocity_bundle$valid_mask)
  if (length(dims) != 3L || !identical(dim(velocity_bundle$hcv_km_yr), dims) ||
      !identical(dim(velocity_bundle$hcv_angle_deg), dims) || !identical(dim(velocity_bundle$vcv_m_yr), dims) ||
      length(depths) != dims[3L] || any(diff(depths) <= 0)) {
    stop("Velocity-bundle arrays and increasing depths must have matching 3D geometry.", call. = FALSE)
  }
  mids <- (depths[-1L] + depths[-length(depths)]) / 2
  depth_edges <- c(depths[1L] - (mids[1L] - depths[1L]), mids, depths[length(depths)] + (depths[length(depths)] - mids[length(mids)]))
  angle_rad <- velocity_bundle$hcv_angle_deg * pi / 180
  vcv <- velocity_bundle$vcv_m_yr
  if (identical(vertical_mode, "zero")) vcv[velocity_bundle$valid_mask] <- 0
  ex <- velocity_bundle$ext
  rs <- as.numeric(velocity_bundle$res)
  x_centers <- as.numeric(velocity_bundle$x_centers)
  y_centers <- as.numeric(velocity_bundle$y_centers)
  ext_names <- c("xmin", "xmax", "ymin", "ymax")
  if (!is.list(ex) || !identical(names(ex), ext_names) ||
      any(!is.finite(unlist(ex, use.names = FALSE))) || ex$xmin >= ex$xmax || ex$ymin >= ex$ymax ||
      length(rs) != 2L || any(!is.finite(rs)) || any(rs <= 0) ||
      length(x_centers) != dims[2L] || length(y_centers) != dims[1L] ||
      any(!is.finite(x_centers)) || any(!is.finite(y_centers))) {
    stop("Velocity bundle has invalid plain geometry; rebuild it with the revision velocity workflow.", call. = FALSE)
  }
  list(
    ext = ex, res = rs, depths = depths, depth_edges = depth_edges,
    nrow = dims[1L], ncol = dims[2L], ndepth = dims[3L],
    x_centers = x_centers,
    y_centers = y_centers,
    hcv_mag_km_yr = velocity_bundle$hcv_km_yr,
    hcv_ang_deg = velocity_bundle$hcv_angle_deg,
    east_m_yr = velocity_bundle$hcv_km_yr * 1000 * sin(angle_rad),
    north_m_yr = velocity_bundle$hcv_km_yr * 1000 * cos(angle_rad),
    vcv_m_yr = vcv, valid_mask = velocity_bundle$valid_mask,
    vertical_mode = vertical_mode,
    time_period = velocity_bundle$period_id,
    dataset_id = velocity_bundle$dataset_id,
    cell_width_m = velocity_bundle$res[1L] * 111320 * cos(y_centers * pi / 180),
    cell_height_m = rep(velocity_bundle$res[2L] * 111320, dims[1L])
  )
}

#' Select a depth policy from a revision velocity bundle
#'
#' @param velocity_bundle Bundle from `build_revision_velocity_bundle()`.
#' @param depth_mode One of `"5m"`, `"10m"`, `"20m"`, `"all"`, or
#'   `"native"`. `"all"` retains every layer in a separately constructed
#'   temperature-aggregation bundle.
#' @param analysis_floor_m Deepest requested analysis centre in metres.
#' @returns A copy of the bundle containing the selected depth array slices.
#'   `native` fails explicitly until source-native BRAN levels are provided;
#'   it is never silently substituted with interpolated layers.
connectivity_select_depth_policy <- function(velocity_bundle, depth_mode,
                                             analysis_floor_m) {
  if (!depth_mode %in% c("5m", "10m", "20m", "all", "native")) {
    stop("Unsupported depth_mode: ", depth_mode, call. = FALSE)
  }
  if (identical(depth_mode, "native")) {
    stop("The revision bundle contains interpolated common depths, not native BRAN levels. Supply native-level temperature inputs before running depth_mode='native'.", call. = FALSE)
  }
  depths <- velocity_bundle$depths_m
  if (analysis_floor_m > max(depths)) {
    stop("Requested analysis floor ", analysis_floor_m, " m exceeds available velocity-bundle depth ", max(depths), " m.", call. = FALSE)
  }
  keep <- depths <= analysis_floor_m
  if (identical(depth_mode, "10m")) keep <- keep & abs((depths - min(depths)) %% 10) < 1e-8
  if (identical(depth_mode, "20m")) keep <- keep & abs((depths - min(depths)) %% 20) < 1e-8
  if (sum(keep) < 3L) stop("Selected depth policy retains fewer than three layers.", call. = FALSE)
  out <- velocity_bundle
  out$depths_m <- depths[keep]
  for (name in c("temperature_trend_c_yr", "h_gradient_east_c_m", "h_gradient_north_c_m", "v_gradient_c_m",
                 "hcv_km_yr", "hcv_angle_deg", "vcv_m_yr", "valid_mask")) {
    if (!is.null(out[[name]])) out[[name]] <- out[[name]][, , keep, drop = FALSE]
  }
  if (!is.null(out$h_gradient_floor)) out$h_gradient_floor <- out$h_gradient_floor[keep]
  if (!is.null(out$v_gradient_floor)) out$v_gradient_floor <- out$v_gradient_floor[keep]
  if (!is.null(out$distribution)) out$distribution <- out$distribution[out$distribution$depth_m %in% out$depths_m, , drop = FALSE]
  out
}

#' Generate matched longitude, latitude, and depth-jittered particle starts
#'
#' @param fields Field list returned by `connectivity_fields_from_bundle()`.
#' @param particles_per_voxel Positive integer starts in every valid voxel.
#' @param random_seed Integer seed used only for reproducible jitter.
#' @returns Data frame of particle IDs, start indices, and coordinates. The
#'   depth draw is uniform within the start layer's edges. Stops for invalid
#'   particle counts or an empty valid mask.
connectivity_seed_table <- function(fields, particles_per_voxel, random_seed) {
  if (!is.numeric(particles_per_voxel) || length(particles_per_voxel) != 1L ||
      is.na(particles_per_voxel) || particles_per_voxel < 1L || particles_per_voxel != as.integer(particles_per_voxel)) {
    stop("particles_per_voxel must be one positive integer.", call. = FALSE)
  }
  set.seed(as.integer(random_seed))
  idx <- which(fields$valid_mask, arr.ind = TRUE)
  if (!nrow(idx)) stop("Cannot seed an empty valid mask.", call. = FALSE)
  idx <- idx[rep(seq_len(nrow(idx)), each = as.integer(particles_per_voxel)), , drop = FALSE]
  n <- nrow(idx)
  z_low <- fields$depth_edges[idx[, 3L]]
  z_high <- fields$depth_edges[idx[, 3L] + 1L]
  data.frame(
    seed_id = seq_len(n), row = idx[, 1L], col = idx[, 2L], k = idx[, 3L],
    x = fields$x_centers[idx[, 2L]] + stats::runif(n, -0.5, 0.5) * fields$res[1L],
    y = fields$y_centers[idx[, 1L]] + stats::runif(n, -0.5, 0.5) * fields$res[2L],
    z = stats::runif(n, z_low, z_high), stringsAsFactors = FALSE
  )
}

#' Build a deterministic, memory-bounded plan for corrected particle starts
#'
#' @param fields Field list returned by `connectivity_fields_from_bundle()`.
#' @param particles_per_voxel Positive integer number of starts in every valid
#'   voxel.
#' @param random_seed Integer master seed for the L'Ecuyer-CMRG chunk streams.
#' @param generation_chunk_size Positive integer number of starts generated at
#'   once. It changes peak memory only: the stored stream state fixes starts
#'   independently of later simulation sub-chunking.
#' @returns A `connectivity_seed_plan` list containing valid voxel indices,
#'   particle-ID bounds, one reproducible RNG state per generation chunk, and
#'   particles seeded by depth. It never materialises the complete start table.
#'
#' @details The function restores the caller's RNG kind and `.Random.seed` on
#' exit. It stops for an empty mask, invalid arguments, or a particle count
#' that cannot be represented by R integer seed IDs.
connectivity_seed_chunk_plan <- function(fields, particles_per_voxel, random_seed,
                                         generation_chunk_size = 15000L) {
  if (!is.numeric(particles_per_voxel) || length(particles_per_voxel) != 1L ||
      is.na(particles_per_voxel) || particles_per_voxel < 1L ||
      particles_per_voxel != as.integer(particles_per_voxel)) {
    stop("particles_per_voxel must be one positive integer.", call. = FALSE)
  }
  if (!is.numeric(generation_chunk_size) || length(generation_chunk_size) != 1L ||
      is.na(generation_chunk_size) || generation_chunk_size < 1L ||
      generation_chunk_size != as.integer(generation_chunk_size)) {
    stop("generation_chunk_size must be one positive integer.", call. = FALSE)
  }
  if (!is.numeric(random_seed) || length(random_seed) != 1L || is.na(random_seed) ||
      random_seed < 0 || random_seed > .Machine$integer.max) {
    stop("random_seed must be one non-negative R integer seed.", call. = FALSE)
  }
  valid_voxels <- which(fields$valid_mask, arr.ind = TRUE)
  if (!nrow(valid_voxels)) stop("Cannot seed an empty valid mask.", call. = FALSE)
  n_particles <- nrow(valid_voxels) * as.double(particles_per_voxel)
  if (!is.finite(n_particles) || n_particles > .Machine$integer.max) {
    stop("Requested particle count exceeds R integer seed-ID capacity.", call. = FALSE)
  }
  n_particles <- as.integer(n_particles)
  generation_chunk_size <- as.integer(generation_chunk_size)
  n_chunks <- as.integer(ceiling(n_particles / generation_chunk_size))
  chunk_start <- as.integer(seq.int(1L, n_particles, by = generation_chunk_size))
  chunk_end <- pmin.int(chunk_start + generation_chunk_size - 1L, n_particles)

  prior_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  prior_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(prior_kind))
    if (had_seed) {
      assign(".Random.seed", prior_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(random_seed))
  rng_states <- matrix(NA_integer_, nrow = n_chunks, ncol = length(.Random.seed))
  for (chunk_i in seq_len(n_chunks)) {
    rng_states[chunk_i, ] <- .Random.seed
    .Random.seed <- parallel::nextRNGStream(.Random.seed)
  }
  colnames(rng_states) <- paste0("rng_state_", seq_len(ncol(rng_states)))
  structure(list(
    valid_voxels = valid_voxels,
    particles_per_voxel = as.integer(particles_per_voxel),
    random_seed = as.integer(random_seed),
    generation_chunk_size = generation_chunk_size,
    n_particles = n_particles, chunk_start = chunk_start, chunk_end = chunk_end,
    rng_states = rng_states,
    seeds_per_depth = as.integer(tabulate(valid_voxels[, 3L], nbins = fields$ndepth) * particles_per_voxel)
  ), class = "connectivity_seed_plan")
}

#' Generate one deterministic corrected-particle start chunk
#'
#' @param seed_plan Plan from `connectivity_seed_chunk_plan()`.
#' @param fields Field list with the same valid-mask geometry used to create
#'   `seed_plan`.
#' @param chunk_index One-based generation-chunk index in `seed_plan`.
#' @returns Data frame with particle IDs, start voxel indices, and jittered
#'   longitude, latitude, and depth. The caller's RNG state is restored.
#'
#' @details Each generation chunk has an independent recorded L'Ecuyer-CMRG
#' state. Thus starts are identical for matched modes and reproducible without
#' allocating starts from other chunks. Stops when the plan and field mask are
#' incompatible or the requested chunk does not exist.
connectivity_seed_chunk <- function(seed_plan, fields, chunk_index) {
  if (!inherits(seed_plan, "connectivity_seed_plan")) {
    stop("seed_plan must be created by connectivity_seed_chunk_plan().", call. = FALSE)
  }
  if (!is.numeric(chunk_index) || length(chunk_index) != 1L || is.na(chunk_index) ||
      chunk_index < 1L || chunk_index > nrow(seed_plan$rng_states) ||
      chunk_index != as.integer(chunk_index)) {
    stop("chunk_index is outside the seed-plan range.", call. = FALSE)
  }
  if (nrow(seed_plan$valid_voxels) != sum(fields$valid_mask) ||
      any(!fields$valid_mask[cbind(seed_plan$valid_voxels[, 1L], seed_plan$valid_voxels[, 2L], seed_plan$valid_voxels[, 3L])])) {
    stop("seed_plan valid voxels do not match the supplied fields.", call. = FALSE)
  }
  chunk_index <- as.integer(chunk_index)
  seed_id <- seq.int(seed_plan$chunk_start[chunk_index], seed_plan$chunk_end[chunk_index])
  voxel_index <- ((seed_id - 1L) %/% seed_plan$particles_per_voxel) + 1L
  idx <- seed_plan$valid_voxels[voxel_index, , drop = FALSE]
  n <- nrow(idx)

  prior_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  prior_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(prior_kind))
    if (had_seed) {
      assign(".Random.seed", prior_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", as.integer(seed_plan$rng_states[chunk_index, ]), envir = .GlobalEnv)
  x_jitter <- stats::runif(n, -0.5, 0.5)
  y_jitter <- stats::runif(n, -0.5, 0.5)
  z_jitter <- stats::runif(n)
  z_low <- fields$depth_edges[idx[, 3L]]
  z_high <- fields$depth_edges[idx[, 3L] + 1L]
  data.frame(
    seed_id = seed_id, row = idx[, 1L], col = idx[, 2L], k = idx[, 3L],
    x = fields$x_centers[idx[, 2L]] + x_jitter * fields$res[1L],
    y = fields$y_centers[idx[, 1L]] + y_jitter * fields$res[2L],
    z = z_low + z_jitter * (z_high - z_low), stringsAsFactors = FALSE
  )
}

#' Identify vertically interior voxels for the primary estimand
#'
#' @param valid_mask Three-dimensional logical valid-water mask ordered row,
#'   column, depth.
#' @returns Logical array. `TRUE` requires a valid voxel that is not the
#'   shallowest/deepest analysis layer and has valid layers immediately above
#'   and below it. Seabed-adjacent cells are therefore excluded.
connectivity_vertical_interior_mask <- function(valid_mask) {
  if (length(dim(valid_mask)) != 3L) stop("valid_mask must be a 3D array.", call. = FALSE)
  out <- array(FALSE, dim = dim(valid_mask))
  if (dim(valid_mask)[3L] < 3L) return(out)
  for (k in 2:(dim(valid_mask)[3L] - 1L)) {
    out[, , k] <- valid_mask[, , k - 1L] & valid_mask[, , k] & valid_mask[, , k + 1L]
  }
  out
}

#' Draw one nonparametric bootstrap mean by index
#'
#' @param values Non-empty finite numeric vector to resample.
#' @returns One numeric mean calculated after drawing `length(values)` indices
#'   with replacement. A scalar input is returned exactly, rather than being
#'   passed to `sample()`'s special numeric-population interpretation.
#'
#' This helper uses R's active random-number state. Callers must set and record
#' a seed when reproducible bootstrap intervals are required. It stops for an
#' empty, non-numeric, or non-finite vector.
connectivity_bootstrap_mean_600 <- function(values) {
  if (!is.numeric(values) || !length(values) || any(!is.finite(values))) {
    stop("values must be a non-empty finite numeric vector.", call. = FALSE)
  }
  mean(values[sample.int(length(values), size = length(values), replace = TRUE)])
}

#' Step one corrected particle chunk and tally exact through visits
#'
#' @param seed_chunk Data frame subset from `connectivity_seed_table()`.
#' @param fields Particle-field list from `connectivity_fields_from_bundle()`.
#' @param dt_years Positive Euler time step in years.
#' @param horizon_years Positive total integration horizon in years divisible by
#'   `dt_years`.
#' @param boundary_policy `"stop"` terminates vertical boundary attempts;
#'   `"slide"` holds depth at the surface, analysis floor, or seabed and
#'   continues the horizontal step when valid.
#' @param sample_path_ids Optional seed IDs retained for diagnostics.
#' @returns Start/end/through arrays, per-particle termination flags, and
#'   sampled paths. Through tallying is per particle and excludes its final
#'   start/end assignments, so results do not depend on chunk boundaries.
simulate_particle_chunk <- function(seed_chunk, fields, dt_years, horizon_years,
                                    boundary_policy = c("stop", "slide"),
                                    sample_path_ids = integer(0)) {
  boundary_policy <- match.arg(boundary_policy)
  steps <- horizon_years / dt_years
  if (!is.finite(steps) || steps < 1L || abs(steps - round(steps)) > 1e-10) {
    stop("horizon_years must be an integer multiple of positive dt_years.", call. = FALSE)
  }
  steps <- as.integer(round(steps))
  n <- nrow(seed_chunk)
  row <- as.integer(seed_chunk$row); col <- as.integer(seed_chunk$col); k <- as.integer(seed_chunk$k)
  x <- as.numeric(seed_chunk$x); y <- as.numeric(seed_chunk$y); z <- as.numeric(seed_chunk$z)
  start_row <- row; start_col <- col; start_k <- k; start_x <- x; start_y <- y; start_z <- z
  active <- rep(TRUE, n)
  hit_surface <- hit_floor <- hit_seabed <- hit_lateral <- hit_land <- hit_missing <- rep(FALSE, n)
  visits <- list(data.frame(seed_index = seq_len(n), lin = linear_index_3d_fast(row, col, k, fields)))
  sample_path_ids <- intersect(as.integer(sample_path_ids), as.integer(seed_chunk$seed_id))
  sample_paths <- vector("list", length(sample_path_ids))
  names(sample_paths) <- as.character(sample_path_ids)
  if (length(sample_path_ids)) {
    for (sample_id in sample_path_ids) {
      sample_index <- match(sample_id, seed_chunk$seed_id)
      sample_paths[[as.character(sample_id)]] <- data.frame(
        seed_id = sample_id, step = 0L, year = 0,
        x = x[sample_index], y = y[sample_index], z = z[sample_index], status = "start"
      )
    }
  }

  for (step in seq_len(steps)) {
    active_i <- which(active)
    if (!length(active_i)) break
    linear <- linear_index_3d_fast(row[active_i], col[active_i], k[active_i], fields)
    east <- fields$east_m_yr[linear]; north <- fields$north_m_yr[linear]; vert <- fields$vcv_m_yr[linear]
    usable <- is.finite(east) & is.finite(north) & is.finite(vert)
    if (any(!usable)) {
      bad <- active_i[!usable]; active[bad] <- FALSE; hit_missing[bad] <- TRUE
      active_i <- active_i[usable]; east <- east[usable]; north <- north[usable]; vert <- vert[usable]
    }
    if (!length(active_i)) next
    prop_y <- y[active_i] + north * dt_years / 111320
    longitude_scale <- 111320 * cos(prop_y * pi / 180)
    prop_x <- x[active_i] + east * dt_years / longitude_scale
    prop_z <- z[active_i] - vert * dt_years
    rc <- row_col_from_xy_3d_fast(prop_x, prop_y, fields)
    new_row <- as.integer(rc$row); new_col <- as.integer(rc$col)
    above <- prop_z < fields$depth_edges[1L]
    below <- prop_z > fields$depth_edges[length(fields$depth_edges)]
    if (any(above)) hit_surface[active_i[above]] <- TRUE
    if (any(below)) hit_floor[active_i[below]] <- TRUE
    if (identical(boundary_policy, "slide")) {
      prop_z <- pmin(pmax(prop_z, fields$depth_edges[1L] + .Machine$double.eps),
                     fields$depth_edges[length(fields$depth_edges)] - .Machine$double.eps)
    }
    new_k <- as.integer(depth_index_from_z_3d_fast(prop_z, fields))
    lateral <- new_row < 1L | new_row > fields$nrow | new_col < 1L | new_col > fields$ncol
    if (any(lateral)) hit_lateral[active_i[lateral]] <- TRUE
    valid <- rep(FALSE, length(active_i))
    in_domain <- !lateral & new_k >= 1L & new_k <= fields$ndepth
    if (any(in_domain)) {
      pos <- which(in_domain)
      target <- linear_index_3d_fast(new_row[pos], new_col[pos], new_k[pos], fields)
      valid[pos] <- fields$valid_mask[target]
      invalid_pos <- pos[!valid[pos]]
      for (j in invalid_pos) {
        valid_depths <- which(fields$valid_mask[new_row[j], new_col[j], ])
        if (!length(valid_depths)) {
          hit_land[active_i[j]] <- TRUE
        } else {
          hit_seabed[active_i[j]] <- TRUE
          if (identical(boundary_policy, "slide")) {
            new_k[j] <- if (prop_z[j] > fields$depth_edges[min(valid_depths)]) max(valid_depths) else min(valid_depths)
            prop_z[j] <- fields$depths[new_k[j]]
            target_j <- linear_index_3d_fast(new_row[j], new_col[j], new_k[j], fields)
            valid[j] <- fields$valid_mask[target_j]
          }
        }
      }
    }
    vertical_stop <- (above | below) & identical(boundary_policy, "stop")
    good <- in_domain & valid & !vertical_stop
    if (any(!good)) active[active_i[!good]] <- FALSE
    if (any(good)) {
      keep <- active_i[good]
      row[keep] <- new_row[good]; col[keep] <- new_col[good]; k[keep] <- new_k[good]
      x[keep] <- prop_x[good]; y[keep] <- prop_y[good]; z[keep] <- prop_z[good]
      visits[[length(visits) + 1L]] <- data.frame(
        seed_index = keep,
        lin = linear_index_3d_fast(row[keep], col[keep], k[keep], fields)
      )
    }
    if (length(sample_path_ids)) {
      for (sample_id in sample_path_ids) {
        sample_index <- match(sample_id, seed_chunk$seed_id)
        sample_paths[[as.character(sample_id)]] <- rbind(
          sample_paths[[as.character(sample_id)]],
          data.frame(
            seed_id = sample_id, step = step, year = step * dt_years,
            x = x[sample_index], y = y[sample_index], z = z[sample_index],
            status = if (active[sample_index] && step < steps) "active" else if (active[sample_index]) "horizon" else "terminated"
          )
        )
      }
    }
  }

  start_counts <- accumulate_counts_3d_fast(start_row, start_col, start_k, fields)
  end_counts <- accumulate_counts_3d_fast(row, col, k, fields)
  all_visits <- unique(do.call(rbind, visits))
  start_lin <- linear_index_3d_fast(start_row, start_col, start_k, fields)
  end_lin <- linear_index_3d_fast(row, col, k, fields)
  through_visits <- all_visits[
    all_visits$lin != start_lin[all_visits$seed_index] & all_visits$lin != end_lin[all_visits$seed_index],
    , drop = FALSE
  ]
  through_counts <- array(
    as.integer(tabulate(through_visits$lin, nbins = fields$nrow * fields$ncol * fields$ndepth)),
    dim = c(fields$nrow, fields$ncol, fields$ndepth)
  )
  reason <- ifelse(active, "horizon", "")
  flags <- list(surface = hit_surface, analysis_floor = hit_floor, seabed = hit_seabed,
                lateral_domain = hit_lateral, land_invalid = hit_land, missing_field = hit_missing)
  for (label in names(flags)) {
    flagged <- flags[[label]]
    reason[flagged & nzchar(reason)] <- paste(reason[flagged & nzchar(reason)], label, sep = ";")
    reason[flagged & !nzchar(reason)] <- label
  }
  diagnostics <- data.frame(
    seed_id = seed_chunk$seed_id, start_row, start_col, start_k, end_row = row, end_col = col, end_k = k,
    start_x, start_y, start_z, end_x = x, end_y = y, end_z = z, active_at_end = active,
    termination_reason = reason, termination_surface = hit_surface,
    termination_analysis_floor = hit_floor, termination_seabed = hit_seabed,
    termination_lateral_domain = hit_lateral, termination_land_invalid = hit_land,
    termination_missing_field = hit_missing,
    horizontal_distance_km = geosphere::distHaversine(cbind(start_x, start_y), cbind(x, y)) / 1000,
    vertical_distance_m = abs(z - start_z), stringsAsFactors = FALSE
  )
  list(start = start_counts, end = end_counts, through = through_counts, diagnostics = diagnostics,
       sample_paths = sample_paths)
}

#' Run the corrected matched particle engine
#'
#' @param fields Particle fields from `connectivity_fields_from_bundle()`.
#' @param seeds Starts from `connectivity_seed_table()`.
#' @param dt_years Positive Euler time step in years.
#' @param horizon_years Positive integration horizon in years.
#' @param boundary_policy Stop or slide policy passed to `simulate_particle_chunk()`.
#' @param chunk_size Positive number of particles per independent chunk.
#' @param cores Requested workers; serial execution is used on Windows.
#' @param sample_path_ids Optional particle IDs whose exact accepted-step paths
#'   are retained for diagnostics. They do not change simulation or tallying.
#' @returns Aggregated exact tally arrays and all particle diagnostics. Results
#'   are invariant to chunk size because chunks contain disjoint particle IDs.
run_connectivity_transition_engine_table_600 <- function(fields, seeds, dt_years, horizon_years,
                                                         boundary_policy = c("stop", "slide"),
                                                         chunk_size = 15000L, cores = 1L,
                                                         sample_path_ids = integer(0)) {
  boundary_policy <- match.arg(boundary_policy)
  if (!is.numeric(chunk_size) || chunk_size < 1L) stop("chunk_size must be positive.", call. = FALSE)
  indexes <- split(seq_len(nrow(seeds)), ceiling(seq_len(nrow(seeds)) / as.integer(chunk_size)))
  worker <- function(index) simulate_particle_chunk(
    seeds[index, , drop = FALSE], fields, dt_years, horizon_years, boundary_policy,
    sample_path_ids = sample_path_ids
  )
  parts <- if (.Platform$OS.type == "windows" || cores <= 1L) lapply(indexes, worker) else parallel::mclapply(indexes, worker, mc.cores = cores)
  combine <- function(name) Reduce(`+`, lapply(parts, `[[`, name))
  list(
    start = combine("start"), end = combine("end"), through = combine("through"),
    diagnostics = do.call(rbind, lapply(parts, `[[`, "diagnostics")),
    sample_paths = unlist(lapply(parts, `[[`, "sample_paths"), recursive = FALSE),
    boundary_policy = boundary_policy
  )
}

#' Incrementally add retention counts from one diagnostic chunk
#'
#' @param retention List of two-dimensional integer arrays named `started`,
#'   `same_column`, and `completed`.
#' @param diagnostics Per-particle diagnostic data frame from
#'   `simulate_particle_chunk()`.
#' @param fields Field list defining the horizontal grid.
#' @returns Updated retention count arrays. The arrays are sufficient to
#'   reproduce water-column retention means without retaining all endpoints.
connectivity_add_retention_counts_600 <- function(retention, diagnostics, fields) {
  start_linear <- (diagnostics$start_row - 1L) * fields$ncol + diagnostics$start_col
  end_linear <- (diagnostics$end_row - 1L) * fields$ncol + diagnostics$end_col
  n_cells <- fields$nrow * fields$ncol
  add_matrix <- function(current, values) {
    current + matrix(as.integer(tabulate(values, nbins = n_cells)), nrow = fields$nrow, ncol = fields$ncol)
  }
  retention$started <- add_matrix(retention$started, start_linear)
  retention$same_column <- add_matrix(retention$same_column, start_linear[diagnostics$start_row == diagnostics$end_row & diagnostics$start_col == diagnostics$end_col])
  retention$completed <- add_matrix(retention$completed, start_linear[diagnostics$active_at_end])
  retention
}

#' Convert compact retention counts to the historical water-column table
#'
#' @param retention List from `connectivity_add_retention_counts_600()`.
#' @returns Data frame with start row/column and the mean retained/completed
#'   indicators for every horizontal cell with at least one particle start.
connectivity_retention_table_600 <- function(retention) {
  if (!all(c("started", "same_column", "completed") %in% names(retention))) {
    stop("retention must contain started, same_column, and completed arrays.", call. = FALSE)
  }
  index <- which(retention$started > 0L, arr.ind = TRUE)
  data.frame(
    start_row = index[, 1L], start_col = index[, 2L],
    retained = retention$same_column[index] / retention$started[index],
    completed = retention$completed[index] / retention$started[index],
    stringsAsFactors = FALSE
  )
}

#' Produce a termination summary from either diagnostic representation
#'
#' @param tallies Output from `run_connectivity_transition_engine()`.
#' @returns Data frame with one row per termination flag, count, particle
#'   denominator, and percentage. Stops when neither complete diagnostics nor
#'   streamed termination counts are present.
connectivity_termination_summary_600 <- function(tallies) {
  flag_names <- c(
    "termination_surface", "termination_analysis_floor", "termination_seabed",
    "termination_lateral_domain", "termination_land_invalid", "termination_missing_field"
  )
  if (!is.null(tallies$diagnostics)) {
    counts <- colSums(tallies$diagnostics[, flag_names, drop = FALSE])
    n_particles <- nrow(tallies$diagnostics)
  } else if (!is.null(tallies$termination_counts) && !is.null(tallies$n_particles)) {
    counts <- tallies$termination_counts[flag_names]
    n_particles <- tallies$n_particles
  } else {
    stop("Tallies do not contain termination diagnostics or streamed counts.", call. = FALSE)
  }
  data.frame(
    cause = sub("^termination_", "", flag_names),
    n_particles = as.integer(n_particles), n_terminated = as.integer(counts),
    percent_particles = 100 * as.numeric(counts) / n_particles,
    stringsAsFactors = FALSE
  )
}

#' Write compact per-particle termination flags for one seed-plan chunk
#'
#' @param diagnostics Diagnostic data frame from one generated seed chunk.
#' @param diagnostic_dir Existing or new directory for immutable chunk RDS files.
#' @param chunk_index One-based seed-plan generation chunk index.
#' @returns One-row manifest data frame. It persists only particle ID, horizon
#'   completion, and all termination flags, which are the required per-particle
#'   termination record; endpoint positions are represented in aggregate
#'   retention counts to bound memory and disk use.
connectivity_write_termination_part_600 <- function(diagnostics, diagnostic_dir, chunk_index) {
  dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
  file_name <- sprintf("termination_flags_chunk_%05d.rds", as.integer(chunk_index))
  file_path <- file.path(diagnostic_dir, file_name)
  if (file.exists(file_path)) stop("Refusing to overwrite termination part: ", file_path, call. = FALSE)
  keep <- c("seed_id", "active_at_end", grep("^termination_", names(diagnostics), value = TRUE))
  saveRDS(diagnostics[, keep, drop = FALSE], file_path, compress = "gzip")
  data.frame(chunk_index = as.integer(chunk_index), file = file_name,
             n_particles = nrow(diagnostics), stringsAsFactors = FALSE)
}

#' Run corrected particles from a deterministic seed plan with bounded memory
#'
#' @param fields Particle fields from `connectivity_fields_from_bundle()`.
#' @param seed_plan `connectivity_seed_plan` from
#'   `connectivity_seed_chunk_plan()`.
#' @param dt_years Positive Euler step in years.
#' @param horizon_years Positive horizon in years.
#' @param boundary_policy `"stop"` or `"slide"`.
#' @param chunk_size Positive simulation sub-chunk size. It does not alter
#'   start positions because generation chunks and their RNG states are fixed
#'   in `seed_plan`.
#' @param diagnostic_dir Optional new directory for compact per-particle
#'   termination-flag chunks. `NULL` retains only aggregate diagnostics.
#' @param sample_path_ids Optional particle IDs for path diagnostics.
#' @returns Exact aggregate tallies, depth-specific seed counts, termination
#'   totals, water-column retention counts, optional termination-part manifest,
#'   and retained paths. It processes one generated seed chunk at a time and
#'   never materialises the full seed table or full diagnostics table.
#'
#' @details Streaming is serial by design. Parallel workers would multiply the
#' peak field/chunk memory and undermine this recovery path. It stops before
#' processing when `diagnostic_dir` is occupied or inputs are invalid.
run_connectivity_transition_engine_plan_600 <- function(fields, seed_plan, dt_years, horizon_years,
                                                        boundary_policy = c("stop", "slide"),
                                                        chunk_size = 15000L,
                                                        diagnostic_dir = NULL,
                                                        sample_path_ids = integer(0)) {
  boundary_policy <- match.arg(boundary_policy)
  if (!inherits(seed_plan, "connectivity_seed_plan")) stop("seed_plan must be a connectivity_seed_plan.", call. = FALSE)
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L || is.na(chunk_size) ||
      chunk_size < 1L || chunk_size != as.integer(chunk_size)) {
    stop("chunk_size must be one positive integer.", call. = FALSE)
  }
  if (!is.null(diagnostic_dir) && dir.exists(diagnostic_dir) &&
      length(list.files(diagnostic_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Refusing to overwrite non-empty termination directory: ", diagnostic_dir, call. = FALSE)
  }
  shape <- c(fields$nrow, fields$ncol, fields$ndepth)
  start <- array(0L, dim = shape)
  end <- array(0L, dim = shape)
  through <- array(0L, dim = shape)
  retention <- list(
    started = matrix(0L, nrow = fields$nrow, ncol = fields$ncol),
    same_column = matrix(0L, nrow = fields$nrow, ncol = fields$ncol),
    completed = matrix(0L, nrow = fields$nrow, ncol = fields$ncol)
  )
  flag_names <- c(
    "termination_surface", "termination_analysis_floor", "termination_seabed",
    "termination_lateral_domain", "termination_land_invalid", "termination_missing_field"
  )
  termination_counts <- stats::setNames(rep.int(0L, length(flag_names)), flag_names)
  manifest <- vector("list", nrow(seed_plan$rng_states))
  retained_paths <- list()
  chunk_size <- as.integer(chunk_size)

  for (generation_i in seq_len(nrow(seed_plan$rng_states))) {
    generated <- connectivity_seed_chunk(seed_plan, fields, generation_i)
    indexes <- split(seq_len(nrow(generated)), ceiling(seq_len(nrow(generated)) / chunk_size))
    diagnostics <- vector("list", length(indexes))
    for (part_i in seq_along(indexes)) {
      part <- simulate_particle_chunk(
        generated[indexes[[part_i]], , drop = FALSE], fields, dt_years, horizon_years,
        boundary_policy, sample_path_ids = sample_path_ids
      )
      start <- start + part$start
      end <- end + part$end
      through <- through + part$through
      diagnostics[[part_i]] <- part$diagnostics
      if (length(part$sample_paths)) retained_paths <- c(retained_paths, part$sample_paths)
    }
    diagnostics <- do.call(rbind, diagnostics)
    termination_counts <- termination_counts + colSums(diagnostics[, flag_names, drop = FALSE])
    retention <- connectivity_add_retention_counts_600(retention, diagnostics, fields)
    if (!is.null(diagnostic_dir)) {
      manifest[[generation_i]] <- connectivity_write_termination_part_600(diagnostics, diagnostic_dir, generation_i)
    }
    rm(generated, diagnostics)
    if (generation_i %% 25L == 0L) gc(FALSE)
  }
  diagnostic_manifest <- if (is.null(diagnostic_dir)) {
    NULL
  } else {
    do.call(rbind, manifest)
  }
  if (!is.null(diagnostic_manifest)) {
    utils::write.csv(diagnostic_manifest, file.path(diagnostic_dir, "termination_flags_manifest.csv"), row.names = FALSE)
  }
  list(
    start = start, end = end, through = through,
    seeds_per_depth = seed_plan$seeds_per_depth,
    n_particles = seed_plan$n_particles,
    termination_counts = termination_counts,
    retention = retention,
    diagnostic_manifest = diagnostic_manifest,
    sample_paths = retained_paths,
    boundary_policy = boundary_policy
  )
}

#' Run the corrected matched particle engine
#'
#' @param fields Particle fields from `connectivity_fields_from_bundle()`.
#' @param seeds Either a complete start data frame from `connectivity_seed_table()`
#'   (compatibility path) or a `connectivity_seed_plan` (memory-bounded path).
#' @param dt_years Positive Euler time step in years.
#' @param horizon_years Positive integration horizon in years.
#' @param boundary_policy Stop or slide policy passed to `simulate_particle_chunk()`.
#' @param chunk_size Positive number of particles per simulation sub-chunk.
#' @param cores Requested workers. Seed-plan execution requires one worker;
#'   the data-frame compatibility path retains its existing parallel behaviour.
#' @param sample_path_ids Optional particle IDs whose exact accepted-step paths
#'   are retained for diagnostics. They do not change simulation or tallying.
#' @param diagnostic_dir Optional new directory for seed-plan termination parts.
#' @returns Aggregate exact tally arrays and either full compatibility
#'   diagnostics or streamed termination/retention diagnostics. Results are
#'   invariant to simulation chunk size for a fixed start table or seed plan.
run_connectivity_transition_engine <- function(fields, seeds, dt_years, horizon_years,
                                               boundary_policy = c("stop", "slide"),
                                               chunk_size = 15000L, cores = 1L,
                                               sample_path_ids = integer(0),
                                               diagnostic_dir = NULL) {
  if (inherits(seeds, "connectivity_seed_plan")) {
    if (!identical(as.integer(cores), 1L)) {
      stop("Memory-bounded seed-plan execution requires cores = 1.", call. = FALSE)
    }
    return(run_connectivity_transition_engine_plan_600(
      fields, seeds, dt_years, horizon_years, boundary_policy, chunk_size,
      diagnostic_dir = diagnostic_dir, sample_path_ids = sample_path_ids
    ))
  }
  if (!is.data.frame(seeds)) stop("seeds must be a data frame or connectivity_seed_plan.", call. = FALSE)
  if (!is.null(diagnostic_dir)) stop("diagnostic_dir is available only with a connectivity_seed_plan.", call. = FALSE)
  run_connectivity_transition_engine_table_600(
    fields, seeds, dt_years, horizon_years, boundary_policy, chunk_size,
    cores, sample_path_ids
  )
}

#' Classify corrected particle tallies and expose vertical-boundary status
#'
#' @param fields Particle field list.
#' @param tallies Result from `run_connectivity_transition_engine()`.
#' @param case_row One-row case data frame from `revision_connectivity_case_plan()`.
#' @param seeds Optional complete start table for the compatibility engine. It
#'   is required only when `tallies` lacks streamed `seeds_per_depth` metadata.
#' @returns Classification list from `classify_3d_fast()` plus vertical-interior
#'   mask, vertical-boundary mask, and detailed sink labels. The broad class IDs
#'   remain comparable with legacy outputs while detailed labels identify
#'   vertical-boundary-influenced sinks.
classify_voxels <- function(fields, tallies, case_row, seeds = NULL) {
  if (!is.data.frame(case_row) || nrow(case_row) != 1L) stop("case_row must be a one-row data frame.", call. = FALSE)
  required <- c("horizon_years", "Nend", "Nst", "NFT")
  if (length(setdiff(required, names(case_row)))) stop("case_row is missing classifier settings.", call. = FALSE)
  seeds_per_depth <- tallies$seeds_per_depth
  if (is.null(seeds_per_depth) && (is.null(seeds) || !is.data.frame(seeds))) {
    stop("tallies must record seeds_per_depth unless a compatibility seed table is supplied.", call. = FALSE)
  }
  out <- classify_3d_fast(
    fields = fields, counts = tallies, seeds = seeds, seeds_per_depth = seeds_per_depth, tyr = case_row$horizon_years,
    class_tyr = case_row$horizon_years, nmL = 5, smL = 20,
    Nend = case_row$Nend, Nst = case_row$Nst, NFT = case_row$NFT,
    sink_balance_tol = 5, corridor_balance_tol = 5, min_count = 1L, norm_mode = "depth"
  )
  interior <- connectivity_vertical_interior_mask(fields$valid_mask)
  vertical_boundary <- fields$valid_mask & !interior
  detail <- out$class_label
  detail[out$class_id %in% c(3L, 4L) & vertical_boundary] <- "Vertical-boundary-influenced sinks"
  out$vertical_interior_mask <- interior
  out$vertical_boundary_mask <- vertical_boundary
  out$class_detail <- detail
  out
}

#' Summarise a paired corrected 2D/3D run
#'
#' @param result_2d Corrected zero-vertical run list.
#' @param result_3d Corrected full-3D run list.
#' @returns List with all-voxel and vertical-interior switching, Cohen's kappa,
#'   a transition table with row proportions, a depth/boundary breakdown, and
#'   water-column retention. Stops when masks or starts differ.
summarise_matched_pair <- function(result_2d, result_3d) {
  seed_identity_2d <- if (!is.null(result_2d$seed_plan)) result_2d$seed_plan else result_2d$seeds
  seed_identity_3d <- if (!is.null(result_3d$seed_plan)) result_3d$seed_plan else result_3d$seeds
  if (!identical(result_2d$fields$valid_mask, result_3d$fields$valid_mask) ||
      !identical(seed_identity_2d, seed_identity_3d)) {
    stop("Matched pair must have identical valid masks and seed tables.", call. = FALSE)
  }
  valid <- result_2d$fields$valid_mask
  interior <- result_3d$classification$vertical_interior_mask
  c2 <- result_2d$classification$class_id; c3 <- result_3d$classification$class_id
  switching <- function(mask) mean(c2[mask] != c3[mask]) * 100
  classes <- sort(unique(c(result_2d$classification$class_lookup$class_id, result_3d$classification$class_lookup$class_id)))
  transitions <- as.data.frame(table(
    matched_2d_class_id = factor(c2[valid], levels = classes),
    matched_3d_class_id = factor(c3[valid], levels = classes)
  ), stringsAsFactors = FALSE)
  transitions$n_voxels <- as.integer(transitions$Freq); transitions$Freq <- NULL
  transitions$row_pct <- ave(transitions$n_voxels, transitions$matched_2d_class_id,
                             FUN = function(x) if (sum(x)) 100 * x / sum(x) else 0)
  observed <- mean(c2[valid] == c3[valid])
  p2 <- prop.table(table(factor(c2[valid], levels = classes)))
  p3 <- prop.table(table(factor(c3[valid], levels = classes)))
  expected <- sum(p2 * p3)
  kappa <- if (expected < 1) (observed - expected) / (1 - expected) else NA_real_
  by_depth <- do.call(rbind, lapply(seq_along(result_2d$fields$depths), function(k) {
    valid_k <- valid[, , k]; interior_k <- interior[, , k]
    data.frame(
      depth_m = result_2d$fields$depths[k], n_valid = sum(valid_k), n_vertical_interior = sum(interior_k),
      pct_switching_all = if (any(valid_k)) switching(array(valid_k, dim = dim(valid_k))) else NA_real_,
      pct_switching_vertical_interior = if (any(interior_k)) switching(array(interior_k, dim = dim(interior_k))) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  retention <- function(tallies) {
    if (!is.null(tallies$retention)) return(connectivity_retention_table_600(tallies$retention))
    diagnostics <- tallies$diagnostics
    same_column <- diagnostics$start_row == diagnostics$end_row & diagnostics$start_col == diagnostics$end_col
    aggregate(cbind(retained = as.integer(same_column), completed = as.integer(diagnostics$active_at_end)) ~ start_row + start_col,
              data = diagnostics, FUN = mean)
  }
  list(
    summary = data.frame(
      n_valid_voxels = sum(valid), n_vertical_interior_voxels = sum(interior),
      pct_switching_all = switching(valid), pct_switching_vertical_interior = switching(interior),
      cohens_kappa = kappa, stringsAsFactors = FALSE
    ),
    transitions = transitions, by_depth = by_depth,
    retention_2d = retention(result_2d$tallies),
    retention_3d = retention(result_3d$tallies)
  )
}

#' Run and persist one seed of a matched corrected case
#'
#' @param velocity_bundle Bundle from `build_revision_velocity_bundle()`.
#' @param case_row One-row revision case table row.
#' @param random_seed Integer seed shared by zero-vertical and full-3D modes.
#' @param output_dir New empty destination for this case/seed.
#' @param cores Positive worker count for the particle engine.
#' @returns List of matched runs and summary. Writes RDS and CSV products only
#'   below `output_dir`; stops before writing when inputs are invalid or the
#'   output directory is non-empty.
run_matched_case <- function(velocity_bundle, case_row, random_seed, output_dir, cores = 1L) {
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Refusing to overwrite non-empty matched-case directory: ", output_dir, call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Could not create matched-case directory.", call. = FALSE)
  if (is.na(case_row$particles_per_voxel[1L])) stop("case_row particles_per_voxel must be resolved before execution.", call. = FALSE)
  selected_bundle <- connectivity_select_depth_policy(
    velocity_bundle, case_row$depth_mode[1L], case_row$analysis_floor_m[1L]
  )
  fields_observed <- connectivity_fields_from_bundle(selected_bundle, "observed")
  fields_zero <- connectivity_fields_from_bundle(selected_bundle, "zero")
  seed_plan <- connectivity_seed_chunk_plan(
    fields_observed, case_row$particles_per_voxel[1L], random_seed,
    generation_chunk_size = 250000L
  )
  run_mode <- function(fields, mode_id) {
    tallies <- run_connectivity_transition_engine(
      fields, seed_plan, dt_years = 0.5, horizon_years = case_row$horizon_years[1L],
      boundary_policy = case_row$boundary_policy[1L], chunk_size = 15000L, cores = cores,
      diagnostic_dir = file.path(output_dir, paste0(mode_id, "_termination_flags"))
    )
    list(fields = fields, seed_plan = seed_plan, tallies = tallies,
         classification = classify_voxels(fields, tallies, case_row))
  }
  result_2d <- run_mode(fields_zero, "matched_2d")
  result_3d <- run_mode(fields_observed, "matched_3d")
  summary <- summarise_matched_pair(result_2d, result_3d)
  saveRDS(result_2d, file.path(output_dir, "matched_2d_result.rds"))
  saveRDS(result_3d, file.path(output_dir, "matched_3d_result.rds"))
  utils::write.csv(summary$summary, file.path(output_dir, "pair_summary.csv"), row.names = FALSE)
  utils::write.csv(summary$transitions, file.path(output_dir, "class_transition.csv"), row.names = FALSE)
  utils::write.csv(summary$by_depth, file.path(output_dir, "switching_by_depth.csv"), row.names = FALSE)
  utils::write.csv(summary$retention_2d, file.path(output_dir, "retention_2d.csv"), row.names = FALSE)
  utils::write.csv(summary$retention_3d, file.path(output_dir, "retention_3d.csv"), row.names = FALSE)
  termination <- rbind(
    transform(connectivity_termination_summary_600(result_2d$tallies), mode = "matched_2d"),
    transform(connectivity_termination_summary_600(result_3d$tallies), mode = "matched_3d")
  )
  termination <- termination[, c("mode", "cause", "n_particles", "n_terminated", "percent_particles")]
  utils::write.csv(termination, file.path(output_dir, "termination_summary.csv"), row.names = FALSE)
  list(matched_2d = result_2d, matched_3d = result_3d, summary = summary, output_dir = output_dir)
}
