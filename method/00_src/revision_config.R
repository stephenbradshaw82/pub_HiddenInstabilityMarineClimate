# Time windows, settings and case definitions from the study.

#' Return the configured historical and projection time periods
#'
#' @returns Named character vector with `T1` and `T2` dates in
#'   `YYYYMMDD-YYYYMMDD` format.
revision_time_periods <- function() {
  c(
    T1 = "19940101-20231231",
    T2 = "20240101-20531231"
  )
}

#' Return the fixed configuration for the connectivity revision
#'
#' @returns Named list of revision decisions. The velocity fields use annual
#'   means and an OLS trend at every voxel. The primary analysis retains 5--210
#'   m centres at 0.25 degrees, applies a depth-specific first-percentile
#'   absolute-gradient floor, and reports vertical-interior voxels as the
#'   primary estimand. This function has no side effects.
revision_connectivity_defaults <- function() {
  list(
    method_version = "connectivity_revision_1.2.0",
    primary_grid_deg = 0.25,
    # The supplied, harmonised analysis cube has these 42 centres only. Every
    # active calculation must be a subset of this common grid; do not infer
    # native or deeper layers.
    primary_depths_m = seq(5, 210, by = 5),
    prepared_depth_policy_id = "2026-11-04_common_grid_5_210m_only",
    trend_method = "annual_mean_ols",
    gradient_floor_quantile = 0.01,
    gradient_floor_sensitivity = c(0, 0.05),
    boundary_policy = "stop",
    boundary_policy_sensitivity = "slide",
    primary_analysis_floor_m = 210,
    truncation_sensitivity_floor_m = 180,
    dt_years = 0.5,
    primary_horizon_years = 10,
    horizon_sensitivity_years = c(5, 20),
    precision_protocol_id = "2026-10-30_extend_to_200_250_keep_noise_gate",
    particles_precision = c(5L, 10L, 25L, 50L, 100L, 200L, 250L),
    precision_seed_count = 2L,
    final_seed_count = 20L,
    cmip6_precision_protocol_id = "2026-11-01_common_t1_t2_density_gate",
    cmip6_precision_particles = c(500L, 750L, 1000L),
    cmip6_precision_seed_count = 2L,
    sensitivity_seed_count = 5L,
    seed_start = 20261029L,
    seed_scheduler = list(
      scheduler_version = "seed_scheduler_1.0.0",
      max_seed_workers = 2L,
      estimated_worker_memory_gib = 3,
      memory_reserve_gib = 4,
      estimated_worker_disk_gib = 2,
      disk_reserve_gib = 20
    ),
    chunk_size = 15000L,
    nmL = 5,
    smL = 20,
    Nend = 10,
    Nst = 10,
    NFT = 10,
    sink_balance_tol = 5,
    corridor_balance_tol = 5,
    min_count = 1L
  )
}

#' Build the prespecified connectivity revision case table
#'
#' @param defaults Output of `revision_connectivity_defaults()`.
#' @returns Data frame with one row per requested precision, primary,
#'   sensitivity, factorial, CMIP6 baseline, CMIP6 precision, CMIP6 refined,
#'   or reconciliation case. Every case is retained; this function never
#'   selects a favourable result.
revision_connectivity_case_plan <- function(defaults = revision_connectivity_defaults()) {
  make_case <- function(stage, case_id, dataset, scenario = NA_character_,
                        time_period, particles_per_voxel, n_seeds,
                        grid_deg = defaults$primary_grid_deg,
                        depth_mode = "5m", horizon_years = defaults$primary_horizon_years,
                        gradient_floor_quantile = defaults$gradient_floor_quantile,
                        boundary_policy = defaults$boundary_policy,
                        analysis_floor_m = defaults$primary_analysis_floor_m,
                        Nend = defaults$Nend, Nst = defaults$Nst, NFT = defaults$NFT,
                        vertical_method = "centre_subsample",
                        aggregation_thickness_m = NA_real_) {
    data.frame(
      stage = stage, case_id = case_id, dataset = dataset, scenario = scenario,
      time_period = time_period, particles_per_voxel = as.integer(particles_per_voxel),
      n_seeds = as.integer(n_seeds), grid_deg = grid_deg, depth_mode = depth_mode,
      horizon_years = horizon_years,
      gradient_floor_quantile = gradient_floor_quantile,
      boundary_policy = boundary_policy, analysis_floor_m = analysis_floor_m,
      Nend = Nend, Nst = Nst, NFT = NFT,
      vertical_method = vertical_method,
      aggregation_thickness_m = aggregation_thickness_m,
      stringsAsFactors = FALSE
    )
  }

  precision <- do.call(rbind, lapply(defaults$particles_precision, function(particles) {
    make_case("precision", paste0("precision_ppv_", particles), "BRAN",
              time_period = revision_time_periods()[["T1"]], particles_per_voxel = particles,
              n_seeds = defaults$precision_seed_count)
  }))
  primary <- make_case("bran_primary", "bran_primary", "BRAN",
                       time_period = revision_time_periods()[["T1"]],
                       particles_per_voxel = NA_integer_, n_seeds = defaults$final_seed_count)
  one_at_a_time <- rbind(
    make_case("sensitivity", "grid_050", "BRAN", time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_, n_seeds = defaults$sensitivity_seed_count, grid_deg = 0.5),
    make_case("sensitivity", "depth_010m", "BRAN", time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_, n_seeds = defaults$sensitivity_seed_count, depth_mode = "10m"),
    make_case("sensitivity", "depth_020m", "BRAN", time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_, n_seeds = defaults$sensitivity_seed_count, depth_mode = "20m"),
    make_case("sensitivity", "horizon_005yr", "BRAN", time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_, n_seeds = defaults$sensitivity_seed_count, horizon_years = 5),
    make_case("sensitivity", "horizon_020yr", "BRAN", time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_, n_seeds = defaults$sensitivity_seed_count, horizon_years = 20),
    make_case("sensitivity", "gradient_floor_none", "BRAN", time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_, n_seeds = defaults$sensitivity_seed_count, gradient_floor_quantile = 0),
    make_case("sensitivity", "gradient_floor_q05", "BRAN", time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_, n_seeds = defaults$sensitivity_seed_count, gradient_floor_quantile = 0.05),
    make_case("sensitivity", "boundary_slide", "BRAN", time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_, n_seeds = defaults$sensitivity_seed_count, boundary_policy = "slide"),
    # This retains 5--180 m from the prepared grid and never requests deeper data.
    make_case("sensitivity", "analysis_floor_180m", "BRAN", time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_, n_seeds = defaults$sensitivity_seed_count, analysis_floor_m = defaults$truncation_sensitivity_floor_m),
    make_case("sensitivity", "tallies_15pct", "BRAN", time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_, n_seeds = defaults$sensitivity_seed_count, Nend = 15, Nst = 15, NFT = 15),
    make_case("sensitivity", "tallies_20pct", "BRAN", time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_, n_seeds = defaults$sensitivity_seed_count, Nend = 20, Nst = 20, NFT = 20)
  )
  factorial <- expand.grid(
    depth_mode = c("5m", "10m"), horizon_years = c(5, 20),
    gradient_floor_quantile = c(0, 0.05), boundary_policy = c("stop", "slide"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  factorial <- do.call(rbind, lapply(seq_len(nrow(factorial)), function(i) {
    x <- factorial[i, ]
    make_case("factorial", paste0("factorial_", sprintf("%02d", i)), "BRAN",
              time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_,
              n_seeds = defaults$sensitivity_seed_count, depth_mode = x$depth_mode,
              horizon_years = x$horizon_years,
              gradient_floor_quantile = x$gradient_floor_quantile,
              boundary_policy = x$boundary_policy)
  }))
  # These cases implement the separately labelled physical-layer aggregation
  # sensitivity. They are not replacements for the completed centre-subsample
  # cases above. Aggregate centres are midpoints of 5 m source-layer groups;
  # the 20 m design retains a documented 10 m terminal group (205--210 m).
  aggregate_one_at_a_time <- rbind(
    make_case("sensitivity", "depth_010m_aggregate", "BRAN",
              time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_,
              n_seeds = defaults$sensitivity_seed_count, depth_mode = "all",
              analysis_floor_m = 207.5, vertical_method = "temperature_layer_aggregate",
              aggregation_thickness_m = 10),
    make_case("sensitivity", "depth_020m_aggregate", "BRAN",
              time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_,
              n_seeds = defaults$sensitivity_seed_count, depth_mode = "all",
              analysis_floor_m = 207.5, vertical_method = "temperature_layer_aggregate",
              aggregation_thickness_m = 20)
  )
  aggregate_factorial <- factorial[factorial$depth_mode == "10m", , drop = FALSE]
  aggregate_factorial <- do.call(rbind, lapply(seq_len(nrow(aggregate_factorial)), function(i) {
    x <- aggregate_factorial[i, , drop = FALSE]
    make_case("factorial", paste0(x$case_id, "_aggregate"), "BRAN",
              time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_,
              n_seeds = defaults$sensitivity_seed_count, depth_mode = "all",
              horizon_years = x$horizon_years,
              gradient_floor_quantile = x$gradient_floor_quantile,
              boundary_policy = x$boundary_policy, analysis_floor_m = 207.5,
              vertical_method = "temperature_layer_aggregate",
              aggregation_thickness_m = 10)
  }))
  cmip6 <- rbind(
    make_case("cmip6", "cmip6_ssp245_t1", "CMIP6", "ssp245", "19940101-20231231", NA_integer_, defaults$final_seed_count),
    make_case("cmip6", "cmip6_ssp245_t2", "CMIP6", "ssp245", "20240101-20531231", NA_integer_, defaults$final_seed_count)
  )
  cmip6_precision <- do.call(rbind, lapply(defaults$cmip6_precision_particles, function(particles) {
    rbind(
      make_case("cmip6_precision", paste0("cmip6_precision_t1_ppv_", particles), "CMIP6", "ssp245",
                "19940101-20231231", particles, defaults$cmip6_precision_seed_count),
      make_case("cmip6_precision", paste0("cmip6_precision_t2_ppv_", particles), "CMIP6", "ssp245",
                "20240101-20531231", particles, defaults$cmip6_precision_seed_count)
    )
  }))
  cmip6_refined <- rbind(
    make_case("cmip6_refined", "cmip6_ssp245_t1_refined", "CMIP6", "ssp245", "19940101-20231231", NA_integer_, defaults$final_seed_count),
    make_case("cmip6_refined", "cmip6_ssp245_t2_refined", "CMIP6", "ssp245", "20240101-20531231", NA_integer_, defaults$final_seed_count)
  )
  reconcile <- make_case("reconciliation", "native_engine_reconciliation", "BRAN",
                          time_period = revision_time_periods()[["T1"]], particles_per_voxel = NA_integer_,
                          n_seeds = defaults$sensitivity_seed_count)
  out <- rbind(precision, primary, one_at_a_time, factorial,
               aggregate_one_at_a_time, aggregate_factorial,
               cmip6, cmip6_precision, cmip6_refined, reconcile)
  rownames(out) <- NULL
  out
}
