# language: R
#
# Final-package paper-figure helpers
#
# These helpers read paired matched outputs only from immutable case
# directories produced by this demo. They never write files or search
# alternative input locations. The calling workflows own all output creation
# after a successful preflight.

connectivity_final_figure_required_seed_files <- function() {
  c(
    "matched_2d_result.rds", "matched_3d_result.rds", "pair_summary.csv",
    "class_transition.csv", "switching_by_depth.csv", "retention_2d.csv",
    "retention_3d.csv", "termination_summary.csv"
  )
}

connectivity_final_figure_assert_scalar_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1L || value != as.integer(value)) {
    stop(name, " must be one positive integer.", call. = FALSE)
  }
  as.integer(value)
}

connectivity_final_figure_assert_array <- function(value, label, expected_dim = NULL) {
  if (!is.array(value) || length(dim(value)) != 3L) {
    stop(label, " must be a three-dimensional array.", call. = FALSE)
  }
  if (!is.null(expected_dim) && !identical(dim(value), expected_dim)) {
    stop(label, " has incompatible dimensions.", call. = FALSE)
  }
  invisible(TRUE)
}

connectivity_final_figure_validate_result <- function(result, label) {
  if (!is.list(result) || !all(c("fields", "seed_plan", "classification") %in% names(result))) {
    stop(label, " is missing fields, seed_plan, or classification.", call. = FALSE)
  }
  fields <- result$fields
  classification <- result$classification
  seed_plan <- result$seed_plan
  if (!is.list(fields) || !is.list(classification) || !is.list(seed_plan)) {
    stop(label, " has malformed fields, seed_plan, or classification.", call. = FALSE)
  }
  required_fields <- c("valid_mask", "depths", "ext", "res", "nrow", "ncol", "ndepth")
  required_classification <- c("class_id", "vertical_interior_mask", "class_lookup")
  required_seed_plan <- c("valid_voxels", "particles_per_voxel", "seeds_per_depth")
  missing_fields <- setdiff(required_fields, names(fields))
  missing_classification <- setdiff(required_classification, names(classification))
  missing_seed_plan <- setdiff(required_seed_plan, names(seed_plan))
  if (length(missing_fields) || length(missing_classification) || length(missing_seed_plan)) {
    stop(
      label, " is missing required components: ",
      paste(c(missing_fields, missing_classification, missing_seed_plan), collapse = ", "),
      call. = FALSE
    )
  }
  expected_dim <- c(as.integer(fields$nrow), as.integer(fields$ncol), as.integer(fields$ndepth))
  if (any(!is.finite(expected_dim)) || any(expected_dim < 1L)) {
    stop(label, " has invalid field dimensions.", call. = FALSE)
  }
  connectivity_final_figure_assert_array(fields$valid_mask, paste0(label, "$fields$valid_mask"), expected_dim)
  connectivity_final_figure_assert_array(classification$class_id, paste0(label, "$classification$class_id"), expected_dim)
  connectivity_final_figure_assert_array(classification$vertical_interior_mask, paste0(label, "$classification$vertical_interior_mask"), expected_dim)
  if (!is.logical(fields$valid_mask) || anyNA(fields$valid_mask) ||
      !is.logical(classification$vertical_interior_mask) || anyNA(classification$vertical_interior_mask)) {
    stop(label, " must use non-missing logical valid and vertical-interior masks.", call. = FALSE)
  }
  if (length(fields$depths) != expected_dim[3L] || any(!is.finite(fields$depths))) {
    stop(label, " has invalid configured depths.", call. = FALSE)
  }
  lookup <- classification$class_lookup
  if (!is.data.frame(lookup) || !all(c("class_id", "class_label") %in% names(lookup)) ||
      !nrow(lookup) || anyDuplicated(lookup$class_id) || anyNA(lookup$class_id) ||
      anyNA(lookup$class_label)) {
    stop(label, " has an invalid class lookup.", call. = FALSE)
  }
  if (any(classification$vertical_interior_mask & !fields$valid_mask)) {
    stop(label, " has vertical-interior voxels outside its valid mask.", call. = FALSE)
  }
  if (!is.numeric(seed_plan$particles_per_voxel) || length(seed_plan$particles_per_voxel) != 1L ||
      !is.finite(seed_plan$particles_per_voxel) || seed_plan$particles_per_voxel < 1L) {
    stop(label, " has invalid particles_per_voxel.", call. = FALSE)
  }
  invisible(TRUE)
}

connectivity_final_figure_same_fields <- function(reference, current, label) {
  checks <- c(
    valid_mask = identical(reference$fields$valid_mask, current$fields$valid_mask),
    depths = identical(reference$fields$depths, current$fields$depths),
    ext = identical(reference$fields$ext, current$fields$ext),
    res = identical(reference$fields$res, current$fields$res),
    nrow = identical(reference$fields$nrow, current$fields$nrow),
    ncol = identical(reference$fields$ncol, current$fields$ncol),
    ndepth = identical(reference$fields$ndepth, current$fields$ndepth),
    vertical_interior_mask = identical(
      reference$classification$vertical_interior_mask,
      current$classification$vertical_interior_mask
    ),
    class_lookup = identical(reference$classification$class_lookup, current$classification$class_lookup)
  )
  if (!all(checks)) {
    stop(label, " differs from the reference paired geometry, mask, depths, or class lookup: ",
         paste(names(checks)[!checks], collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

connectivity_final_figure_read_seed_pair <- function(seed_dir,
                                                      expected_particles_per_voxel = NULL,
                                                      reference = NULL) {
  if (!is.character(seed_dir) || length(seed_dir) != 1L || is.na(seed_dir) || !dir.exists(seed_dir)) {
    stop("seed_dir must be one existing seed directory.", call. = FALSE)
  }
  required_files <- file.path(seed_dir, connectivity_final_figure_required_seed_files())
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files)) {
    stop("Final matched seed is missing required products: ", paste(missing_files, collapse = "; "), call. = FALSE)
  }
  pair_summary <- utils::read.csv(file.path(seed_dir, "pair_summary.csv"), stringsAsFactors = FALSE)
  required_pair_columns <- c(
    "n_valid_voxels", "n_vertical_interior_voxels", "pct_switching_all",
    "pct_switching_vertical_interior", "cohens_kappa"
  )
  if (nrow(pair_summary) != 1L || !all(required_pair_columns %in% names(pair_summary)) ||
      any(!is.finite(unlist(pair_summary[1L, required_pair_columns], use.names = FALSE)))) {
    stop("Final matched seed has a malformed pair_summary.csv: ", seed_dir, call. = FALSE)
  }
  matched_2d <- readRDS(file.path(seed_dir, "matched_2d_result.rds"))
  matched_3d <- readRDS(file.path(seed_dir, "matched_3d_result.rds"))
  connectivity_final_figure_validate_result(matched_2d, paste0(seed_dir, "/matched_2d_result.rds"))
  connectivity_final_figure_validate_result(matched_3d, paste0(seed_dir, "/matched_3d_result.rds"))
  connectivity_final_figure_same_fields(matched_2d, matched_3d, paste0("Paired modes in ", seed_dir))
  if (!identical(matched_2d$seed_plan$valid_voxels, matched_3d$seed_plan$valid_voxels) ||
      !identical(matched_2d$seed_plan$seeds_per_depth, matched_3d$seed_plan$seeds_per_depth)) {
    stop("Paired modes have divergent valid-voxel plans or seed plans: ", seed_dir, call. = FALSE)
  }
  if (!is.null(expected_particles_per_voxel) &&
      (!identical(as.numeric(matched_2d$seed_plan$particles_per_voxel), as.numeric(expected_particles_per_voxel)) ||
       !identical(as.numeric(matched_3d$seed_plan$particles_per_voxel), as.numeric(expected_particles_per_voxel)))) {
    stop("Final matched seed has unexpected particles_per_voxel: ", seed_dir, call. = FALSE)
  }
  if (!is.null(reference)) {
    connectivity_final_figure_same_fields(reference$matched_2d, matched_2d, paste0("2D replicate ", seed_dir))
    connectivity_final_figure_same_fields(reference$matched_3d, matched_3d, paste0("3D replicate ", seed_dir))
  }
  list(
    seed_dir = normalizePath(seed_dir, winslash = "/", mustWork = TRUE),
    matched_2d = matched_2d,
    matched_3d = matched_3d,
    pair_summary = pair_summary
  )
}

#' Validate the final paired-seed input contract without writing output.
#'
#' @param case_dir Existing final-package case directory.
#' @param expected_seeds Required paired seed count; the approved paper cases use 20.
#' @param expected_particles_per_voxel Required particle density for the case.
#' @returns A contract object holding seed paths, the shared template, and
#'   source metadata. It stops on missing products, malformed summaries,
#'   mismatched geometries/masks/depths, divergent paired seed plans, or an
#'   unexpected seed set.
connectivity_final_figure_validate_case <- function(case_dir,
                                                     expected_seeds = 20L,
                                                     expected_particles_per_voxel = NULL) {
  expected_seeds <- connectivity_final_figure_assert_scalar_integer(expected_seeds, "expected_seeds")
  if (!is.character(case_dir) || length(case_dir) != 1L || is.na(case_dir) || !dir.exists(case_dir)) {
    stop("case_dir must be one existing final-package case directory.", call. = FALSE)
  }
  if (!is.null(expected_particles_per_voxel)) {
    expected_particles_per_voxel <- connectivity_final_figure_assert_scalar_integer(
      expected_particles_per_voxel, "expected_particles_per_voxel"
    )
  }
  expected_names <- sprintf("seed_%03d", seq_len(expected_seeds))
  seed_dirs <- list.dirs(case_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("^seed_[0-9]+$", basename(seed_dirs))]
  seed_dirs <- seed_dirs[order(basename(seed_dirs))]
  if (!identical(basename(seed_dirs), expected_names)) {
    stop("Final case has incomplete or unexpected paired seed directories: ", case_dir, call. = FALSE)
  }
  reference <- NULL
  pair_summaries <- vector("list", length(seed_dirs))
  for (i in seq_along(seed_dirs)) {
    current <- connectivity_final_figure_read_seed_pair(
      seed_dir = seed_dirs[i],
      expected_particles_per_voxel = expected_particles_per_voxel,
      reference = reference
    )
    if (is.null(reference)) reference <- current
    pair_summary <- current$pair_summary
    pair_summary$seed <- basename(seed_dirs[i])
    pair_summaries[[i]] <- pair_summary[, c("seed", setdiff(names(pair_summary), "seed")), drop = FALSE]
  }
  template <- list(
    fields = reference$matched_2d$fields,
    classification = list(
      vertical_interior_mask = reference$matched_2d$classification$vertical_interior_mask,
      class_lookup = reference$matched_2d$classification$class_lookup
    )
  )
  list(
    case_dir = normalizePath(case_dir, winslash = "/", mustWork = TRUE),
    seed_dirs = normalizePath(seed_dirs, winslash = "/", mustWork = TRUE),
    expected_seeds = expected_seeds,
    expected_particles_per_voxel = expected_particles_per_voxel,
    template = template,
    pair_summaries = do.call(rbind, pair_summaries)
  )
}

connectivity_final_figure_normalise_zones <- function(zones, nrow, ncol) {
  if (!is.list(zones) || !length(zones) || is.null(names(zones)) || any(!nzchar(names(zones)))) {
    stop("zones must be a named non-empty list.", call. = FALSE)
  }
  lapply(seq_along(zones), function(i) {
    zone <- zones[[i]]
    if (!is.list(zone) || !all(c("mask", "zone_type") %in% names(zone))) {
      stop("Each zone must provide mask and zone_type.", call. = FALSE)
    }
    if (!is.matrix(zone$mask) || !identical(dim(zone$mask), c(nrow, ncol))) {
      stop("Zone mask has incompatible dimensions: ", names(zones)[i], call. = FALSE)
    }
    if (!is.logical(zone$mask) || anyNA(zone$mask) || !is.character(zone$zone_type) ||
        length(zone$zone_type) != 1L || is.na(zone$zone_type) || !nzchar(zone$zone_type)) {
      stop("Zone mask or zone_type is invalid: ", names(zones)[i], call. = FALSE)
    }
    zone
  }) |>
    stats::setNames(names(zones))
}

connectivity_final_figure_class_rows <- function(class_ids, eligible, lookup,
                                                  zone_name, zone_type, depth_m,
                                                  method_label, column_label, seed_name) {
  values <- class_ids[eligible]
  finite_values <- values[is.finite(values)]
  classified <- finite_values[finite_values %in% lookup$class_id]
  out <- data.frame(
    zone_name = zone_name,
    zone_type = zone_type,
    depth_m = depth_m,
    method_label = method_label,
    column_label = column_label,
    seed = seed_name,
    class_id = lookup$class_id,
    class_label = as.character(lookup$class_label),
    n_cells = vapply(lookup$class_id, function(id) sum(finite_values == id), numeric(1)),
    n_eligible = length(values),
    n_classified = length(classified),
    stringsAsFactors = FALSE
  )
  out
}

connectivity_final_figure_summarise_profile_rows <- function(raw) {
  required_columns <- c(
    "zone_name", "zone_type", "depth_m", "method_label", "column_label",
    "seed", "class_id", "class_label", "n_cells", "n_eligible", "n_classified"
  )
  if (!is.data.frame(raw) || !nrow(raw) || !all(required_columns %in% names(raw))) {
    stop("raw profile rows are empty or malformed.", call. = FALSE)
  }
  keys <- c("zone_name", "zone_type", "depth_m", "method_label", "column_label", "class_id", "class_label")
  grouping <- interaction(raw[keys], drop = TRUE, lex.order = TRUE)
  summed <- lapply(split(raw, grouping), function(x) {
    out <- x[1L, keys, drop = FALSE]
    out$n_cells <- sum(x$n_cells)
    out$n_eligible <- sum(x$n_eligible)
    out$n_classified <- sum(x$n_classified)
    out$n_seed_replicates <- length(unique(x$seed))
    out
  })
  out <- do.call(rbind, summed)
  out$pct_eligible <- ifelse(out$n_eligible > 0, 100 * out$n_cells / out$n_eligible, NA_real_)
  out$pct_classified <- ifelse(out$n_classified > 0, 100 * out$n_cells / out$n_classified, NA_real_)
  rownames(out) <- NULL
  out[order(out$zone_name, out$method_label, out$depth_m, out$class_id), , drop = FALSE]
}

#' Aggregate final matched classifications across paired seeds.
#'
#' The numerator is the count of classified voxel-seed observations in each
#' class; the primary denominator is the vertical-interior eligible
#' voxel-seed count. The output also retains a classified-only percentage so
#' a missing native/benchmark class can never be hidden by renormalisation.
connectivity_final_figure_aggregate_profiles <- function(contract, zones, column_label) {
  if (!is.list(contract) || !all(c("seed_dirs", "template", "expected_seeds") %in% names(contract))) {
    stop("contract must come from connectivity_final_figure_validate_case().", call. = FALSE)
  }
  fields <- contract$template$fields
  lookup <- contract$template$classification$class_lookup
  zones <- connectivity_final_figure_normalise_zones(zones, fields$nrow, fields$ncol)
  rows <- list()
  for (seed_dir in contract$seed_dirs) {
    pair <- connectivity_final_figure_read_seed_pair(
      seed_dir = seed_dir,
      expected_particles_per_voxel = contract$expected_particles_per_voxel
    )
    connectivity_final_figure_same_fields(
      list(fields = fields, classification = contract$template$classification),
      pair$matched_2d,
      paste0("Matched 2D aggregation input ", seed_dir)
    )
    eligible <- pair$matched_2d$fields$valid_mask &
      pair$matched_2d$classification$vertical_interior_mask
    for (method in c("Matched 2D", "Matched 3D")) {
      class_ids <- if (identical(method, "Matched 2D")) {
        pair$matched_2d$classification$class_id
      } else {
        pair$matched_3d$classification$class_id
      }
      for (zone_name in names(zones)) {
        zone <- zones[[zone_name]]
        for (k in seq_along(fields$depths)) {
          rows[[length(rows) + 1L]] <- connectivity_final_figure_class_rows(
            class_ids = class_ids[, , k],
            eligible = eligible[, , k] & zone$mask,
            lookup = lookup,
            zone_name = zone_name,
            zone_type = zone$zone_type,
            depth_m = fields$depths[k],
            method_label = method,
            column_label = column_label,
            seed_name = basename(seed_dir)
          )
        }
      }
    }
  }
  raw <- do.call(rbind, rows)
  connectivity_final_figure_summarise_profile_rows(raw)
}

#' Aggregate paired 2D--3D disagreement using the primary eligible denominator.
#'
#' @returns Voxel and horizontal-column numerator/denominator arrays plus a
#' data frame of per-column values. Every denominator counts paired
#' seed--depth observations satisfying the final vertical-interior eligibility
#' mask; it is calculated from the paired seed results rather than inferred
#' from a class stack.
connectivity_final_figure_aggregate_disagreement <- function(contract) {
  if (!is.list(contract) || !all(c("seed_dirs", "template", "expected_seeds") %in% names(contract))) {
    stop("contract must come from connectivity_final_figure_validate_case().", call. = FALSE)
  }
  fields <- contract$template$fields
  dims <- dim(fields$valid_mask)
  eligible_count <- array(0L, dim = dims)
  switch_count <- array(0L, dim = dims)
  for (seed_dir in contract$seed_dirs) {
    pair <- connectivity_final_figure_read_seed_pair(
      seed_dir = seed_dir,
      expected_particles_per_voxel = contract$expected_particles_per_voxel
    )
    connectivity_final_figure_same_fields(
      list(fields = fields, classification = contract$template$classification),
      pair$matched_2d,
      paste0("Disagreement aggregation input ", seed_dir)
    )
    eligible <- pair$matched_2d$fields$valid_mask &
      pair$matched_2d$classification$vertical_interior_mask
    switched <- eligible &
      is.finite(pair$matched_2d$classification$class_id) &
      is.finite(pair$matched_3d$classification$class_id) &
      pair$matched_2d$classification$class_id != pair$matched_3d$classification$class_id
    eligible_count <- eligible_count + as.integer(eligible)
    switch_count <- switch_count + as.integer(switched)
  }
  switch_proportion <- ifelse(eligible_count > 0L, switch_count / eligible_count, NA_real_)
  column_eligible_count <- apply(eligible_count, c(1L, 2L), sum)
  column_switch_count <- apply(switch_count, c(1L, 2L), sum)
  column_switch_proportion <- ifelse(
    column_eligible_count > 0L,
    column_switch_count / column_eligible_count,
    NA_real_
  )
  index <- which(column_eligible_count > 0L, arr.ind = TRUE)
  columns <- data.frame(
    row = index[, 1L],
    col = index[, 2L],
    longitude = fields$x_centers[index[, 2L]],
    latitude = fields$y_centers[index[, 1L]],
    eligible_seed_depth_pairs = column_eligible_count[index],
    switched_seed_depth_pairs = column_switch_count[index],
    switch_proportion = column_switch_proportion[index],
    stringsAsFactors = FALSE
  )
  list(
    eligible_seed_depth_pairs = eligible_count,
    switched_seed_depth_pairs = switch_count,
    switch_proportion = switch_proportion,
    column_eligible_seed_depth_pairs = column_eligible_count,
    column_switched_seed_depth_pairs = column_switch_count,
    column_switch_proportion = column_switch_proportion,
    columns = columns,
    n_seed_replicates = length(contract$seed_dirs)
  )
}
