# Helpers for native VoCC 2D classification used by script 04.
# The active classifier is provided by the pinned VoCC package.

#' Derive native-2D trajectory starts per valid analysis cell
#'
#' @param traj Trajectory data frame from `VoCC::voccTraj()` with `trajIDs` and
#'   x/y coordinates in its first two columns.
#' @param template_raster Raster defining the analysis-cell geometry used to
#'   locate trajectory starts.
#' @param valid_raster Raster whose non-missing cells define the valid analysis
#'   mask; defaults to `template_raster`.
#' @returns List containing the scalar `trajSt`, total trajectory and valid-cell
#'   counts, and a per-cell start-count table. Stops if starts lie outside the
#'   mask or if valid cells do not have a common number of starts.
#'
#' `trajSt` populates native local start counts. It is not the denominator for
#' native class percentages, which use local `start + end + through` totals.
derive_trajSt_from_trajectory_starts <- function(traj, template_raster, valid_raster = template_raster) {
  require(data.table)
  require(raster)

  if (!"trajIDs" %in% names(traj)) {
    stop("Trajectory table must contain a `trajIDs` column.", call. = FALSE)
  }
  if (ncol(traj) < 2) {
    stop("Trajectory table must contain x/y coordinate columns in columns 1 and 2.", call. = FALSE)
  }

  traj_dt <- data.table::as.data.table(traj)
  coord_cols <- names(traj_dt)[1:2]
  start_dt <- traj_dt[, .SD[1], by = trajIDs]
  start_xy <- as.data.frame(start_dt[, ..coord_cols])
  start_cells <- raster::cellFromXY(template_raster, start_xy)
  valid_start <- !is.na(start_cells) &
    start_cells > 0L &
    start_cells <= raster::ncell(template_raster)

  if (!all(valid_start)) {
    n_invalid <- sum(!valid_start)
    stop("Found ", n_invalid, " trajectory start(s) outside the analysis raster.", call. = FALSE)
  }

  start_counts <- data.table::data.table(cid = start_cells)[, .N, by = cid]
  valid_cells <- which(!is.na(raster::values(valid_raster)))
  missing_start_cells <- setdiff(valid_cells, start_counts$cid)

  if (length(missing_start_cells) > 0L) {
    stop(
      "Found ", length(missing_start_cells),
      " valid analysis cell(s) without a trajectory start.",
      call. = FALSE
    )
  }

  extra_start_cells <- setdiff(start_counts$cid, valid_cells)
  if (length(extra_start_cells) > 0L) {
    stop(
      "Found ", length(extra_start_cells),
      " trajectory start cell(s) outside the valid analysis mask.",
      call. = FALSE
    )
  }

  count_values <- sort(unique(start_counts$N))
  if (length(count_values) != 1L) {
    dist <- paste(
      paste0(names(table(start_counts$N)), " starts: ", as.integer(table(start_counts$N))),
      collapse = "; "
    )
    stop(
      "Trajectory starts per valid analysis cell are not constant: ",
      dist,
      call. = FALSE
    )
  }

  list(
    trajSt = as.integer(count_values),
    n_trajectories = nrow(start_dt),
    n_valid_cells = length(valid_cells),
    starts_per_cell = as.data.frame(start_counts[order(cid)])
  )
}

#' Run the version-pinned native VoCC 2D classifier
#'
#' @param traj Trajectory data frame returned by `VoCC::voccTraj()`.
#' @param vel Raster of HVoCC magnitude in kilometres per year.
#' @param ang Raster of HVoCC compass bearings in degrees.
#' @param mn Raster of mean climate values supplied to `VoCC::trajClas()`.
#' @param trajSt Integer starts per valid analysis cell, derived with
#'   `derive_trajSt_from_trajectory_starts()`.
#' @param tyr Classification horizon in years.
#' @param nmL Non-moving distance threshold in kilometres over `tyr`.
#' @param smL Slow-moving distance threshold in kilometres over `tyr`.
#' @param Nend Local ending-percentage threshold for sink rules.
#' @param Nst Local low-starting-percentage threshold for Relative Sink rules.
#' @param NFT Local through-flow-percentage threshold for Corridor rules.
#' @param DateLine Logical; whether the raster extent crosses the date line.
#' @returns The `VoCC::trajClas()` result. Stops with installation guidance when
#'   VoCC is unavailable; this wrapper does not substitute another classifier.
reference_trajClas <- function(traj, vel, ang, mn, trajSt, tyr, nmL, smL,
                               Nend, Nst, NFT, DateLine = FALSE) {
  if (!requireNamespace("VoCC", quietly = TRUE)) {
    stop(
      "VoCC is required for the reference 2D trajectory classifier. ",
      "Install/version-pin VoCC before running the native/reference 2D workflow; ",
      "no substitute classifier is used.",
      call. = FALSE
    )
  }

  VoCC::trajClas(
    traj = traj,
    vel = vel,
    ang = ang,
    mn = mn,
    trajSt = trajSt,
    tyr = tyr,
    nmL = nmL,
    smL = smL,
    Nend = Nend,
    Nst = Nst,
    NFT = NFT,
    DateLine = DateLine
  )
}
