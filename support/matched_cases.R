# Keep the fixed case settings and seeds. The selected particle densities are
# inputs to this demo.

matched_case_plan <- function(cases) {
  cases <- cases[cases$stage %in% c("bran_primary", "cmip6_refined", "sensitivity", "factorial"), ]
  cases$particles_per_voxel[cases$dataset == "BRAN"] <- 250L
  cases$particles_per_voxel[cases$dataset == "CMIP6"] <- 750L
  cases$phase <- ifelse(cases$stage %in% c("sensitivity", "factorial"), "sensitivity", "matched")
  cases$output_group <- ifelse(cases$phase == "sensitivity", "04_sensitivity",
                               ifelse(cases$dataset == "CMIP6", "05_cmip6_ssp245", "03_bran_primary"))
  if (sum(cases$phase == "matched" & cases$n_seeds == 20L) != 3L ||
      sum(cases$phase == "sensitivity" & cases$n_seeds == 5L) != 37L ||
      nrow(cases) != 40L || anyDuplicated(cases$case_id)) {
    stop("Expected the three paper headline cases and all 37 amended sensitivity cases.")
  }
  cases
}

validate_matched_pair <- function(two, three, particles_per_voxel) {
  valid <- two$fields$valid_mask
  interior <- two$classification$vertical_interior_mask
  if (!identical(valid, three$fields$valid_mask) ||
      !identical(interior, three$classification$vertical_interior_mask) ||
      !identical(two$seed_plan, three$seed_plan) ||
      !identical(two$fields$east_m_yr, three$fields$east_m_yr) ||
      !identical(two$fields$north_m_yr, three$fields$north_m_yr) ||
      !all(two$fields$vcv_m_yr[valid] == 0) || !any(interior)) {
    stop("Invalid matched fields, starting positions or eligibility mask.")
  }
  total <- sum(valid) * particles_per_voxel
  for (mode in list(two, three)) {
    if (sum(mode$tallies$start) != total || sum(mode$tallies$end) != total ||
        anyNA(mode$classification$class_id[valid])) {
      stop("Particle counts or class IDs do not reconcile.")
    }
  }
}

regional_replicate_summary <- function(seed_dirs, reporting_mask, report_seed = 20261029L) {
  # Read one seed pair at a time. Only cells in the reporting region contribute
  # to the ledger, even though particles can start and travel in the buffer.
  raw <- all_depths <- kappa <- numeric(length(seed_dirs))
  noise <- matrix(NA_real_, length(seed_dirs) %/% 2L, 2)
  for (i in seq_along(seed_dirs)) {
    two <- readRDS(file.path(seed_dirs[i], "matched_2d_result.rds"))
    three <- readRDS(file.path(seed_dirs[i], "matched_3d_result.rds"))
    valid <- two$fields$valid_mask
    if (!identical(dim(reporting_mask), dim(valid)[1:2])) stop("Reporting and model grids differ.")
    region <- array(rep(reporting_mask, dim(valid)[3]), dim(valid))
    valid <- valid & region
    interior <- two$classification$vertical_interior_mask & region
    if (!any(interior)) stop("No eligible vertical-interior cells in the reporting region.")
    if (i == 1L) {
      reference_valid <- valid
      reference_interior <- interior
    } else if (!identical(valid, reference_valid) || !identical(interior, reference_interior)) {
      stop("Reporting eligibility changed between seeds.")
    }
    two_class <- two$classification$class_id
    three_class <- three$classification$class_id
    raw[i] <- 100 * mean(two_class[interior] != three_class[interior])
    all_depths[i] <- 100 * mean(two_class[valid] != three_class[valid])
    classes <- sort(unique(c(two$classification$class_lookup$class_id,
                             three$classification$class_lookup$class_id)))
    p2 <- prop.table(table(factor(two_class[valid], levels = classes)))
    p3 <- prop.table(table(factor(three_class[valid], levels = classes)))
    expected <- sum(p2 * p3)
    kappa[i] <- if (expected < 1) (1 - all_depths[i] / 100 - expected) / (1 - expected) else NA_real_
    # Independent noise estimates use adjacent, non-overlapping seed pairs.
    if (i %% 2L) {
      previous_two <- two_class[interior]
      previous_three <- three_class[interior]
    } else {
      noise[i %/% 2L, ] <- 100 * c(mean(previous_two != two_class[interior]),
                                   mean(previous_three != three_class[interior]))
    }
  }
  if (nrow(noise) < 1L) stop("A replicate summary needs at least two seeds.")
  interval <- c(NA_real_, NA_real_)
  if (nrow(noise) >= 2L) {
    set.seed(report_seed)
    draws <- replicate(10000L, mean(sample(raw, length(raw), replace = TRUE)) -
                                mean(sample(noise[, 1], nrow(noise), replace = TRUE)))
    interval <- unname(stats::quantile(draws, c(0.025, 0.975)))
  }
  data.frame(n_mode_pairs = length(raw), n_noise_pairs = nrow(noise), n_eligible_per_seed = sum(interior),
             mean_switching_vertical_interior = mean(raw), mean_switching_all = mean(all_depths),
             mean_noise_2d = mean(noise[, 1]), mean_noise_3d = mean(noise[, 2]),
             attributable_switching = mean(raw) - mean(noise[, 1]),
             bootstrap_ci_available = nrow(noise) >= 2L,
             attributable_ci_low = interval[1], attributable_ci_high = interval[2],
             mean_cohens_kappa = mean(kappa))
}

run_paired_seed <- function(bundle, case, random_seed, seed_dir, run_pair) {
  # Retain incomplete attempts. Publish one complete pair only after validation.
  attempt <- tempfile(paste0(".", basename(seed_dir), "-attempt-"), tmpdir = dirname(seed_dir))
  elapsed <- system.time({
    result <- run_pair(bundle, case, random_seed, attempt, cores = 1L)
    validate_matched_pair(result$matched_2d, result$matched_3d, case$particles_per_voxel)
  })[["elapsed"]]
  termination <- utils::read.csv(file.path(attempt, "termination_summary.csv"))
  expected_counts <- sum(result$matched_2d$fields$valid_mask) * case$particles_per_voxel
  if (nrow(termination) != 12L || length(unique(termination$mode)) != 2L ||
      any(termination$n_particles != expected_counts) ||
      any(termination$n_terminated < 0L | termination$n_terminated > expected_counts)) {
    stop("Invalid termination counts: ", attempt)
  }
  outputs <- list.files(attempt, recursive = TRUE, full.names = TRUE)
  file_names <- substring(outputs, nchar(attempt) + 2L)
  publish_csv(data.frame(seed = random_seed, elapsed_s = elapsed,
                         output_bytes = sum(file.info(outputs)$size)),
              file.path(attempt, "timing.csv"))
  outputs <- c(outputs, file.path(attempt, "timing.csv"))
  publish_csv(data.frame(file = c(file_names, "timing.csv"), md5 = unname(tools::md5sum(outputs))),
              file.path(attempt, "output_checksums.csv"))
  if (file.exists(seed_dir) || !file.rename(attempt, seed_dir)) stop("Could not publish ", seed_dir)
  message(case$case_id, " ", basename(seed_dir), " (", round(elapsed, 2), " s)")
}

run_regional_case <- function(case, run_dir, fixture_dir, production_files, pilot_only, run_pair) {
  dataset <- if (case$dataset == "CMIP6") "CMIP6_ssp245" else "BRAN"
  period <- names(time_periods)[match(case$time_period, time_periods)]
  variant <- velocity_variant_name(case$grid_deg, case$gradient_floor_quantile,
                                   case$aggregation_thickness_m)
  field_dir <- file.path(run_dir, "velocity", dataset, period, variant)
  check_file_manifest(file.path(field_dir, "files.csv"))
  field_path <- file.path(field_dir, "velocity_bundle.rds")
  bundle <- readRDS(field_path)
  grid <- bundle_grid(bundle)
  verify_grid_origin(grid, case$grid_deg)
  region_path <- file.path(fixture_dir, "regions.gpkg")
  reporting_mask <- region_mask(terra::vect(region_path, layer = "reporting"), grid)
  simulation_mask <- region_mask(terra::vect(region_path, layer = "simulation"), grid)
  if (any(reporting_mask & !simulation_mask)) stop("Reporting cells fall outside the calculation mask.")
  # The engine uses this field mask for both seeding and movement.
  bundle$valid_mask <- bundle$valid_mask &
    array(rep(simulation_mask, length(bundle$depths_m)), dim(bundle$valid_mask))
  selected <- connectivity_select_depth_policy(bundle, case$depth_mode, case$analysis_floor_m)
  reporting_layers <- array(rep(reporting_mask, length(selected$depths_m)), dim(selected$valid_mask))
  if (!any(connectivity_vertical_interior_mask(selected$valid_mask) & reporting_layers)) {
    stop("No eligible Gascoyne cells in ", case$case_id)
  }

  root <- file.path(run_dir, "cases", case$output_group, paste0("case_", case$case_id))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(root, ".running")
  if (!dir.create(lock, showWarnings = FALSE)) stop("This case is already running or has a stale lock: ", lock)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  seed_start <- revision_connectivity_defaults()$seed_start
  case_recipe <- list(case = case, seed_start = seed_start,
                      files = unname(tools::md5sum(c(field_path, region_path, production_files))),
                      r_version = as.character(getRversion()),
                      geosphere = as.character(utils::packageVersion("geosphere")))
  check_stage_recipe(file.path(root, "recipe.rds"), case_recipe)

  seed_dirs <- file.path(root, sprintf("seed_%03d", seq_len(if (pilot_only) 1L else case$n_seeds)))
  for (i in seq_along(seed_dirs)) {
    if (dir.exists(seed_dirs[i])) {
      check_file_manifest(file.path(seed_dirs[i], "output_checksums.csv"))
    } else {
      # The production call runs both modes from one seed plan. Only the
      # zero-vertical control has its vertical velocity set to zero.
      run_paired_seed(bundle, case, seed_start + i - 1L, seed_dirs[i], run_pair)
    }
  }
  if (pilot_only) return(invisible(NULL))

  summary <- regional_replicate_summary(seed_dirs, reporting_mask, seed_start)
  summary_path <- file.path(root, "replicate_summary.csv")
  if (!file.exists(summary_path)) {
    publish_csv(summary, summary_path)
  } else if (!isTRUE(all.equal(utils::read.csv(summary_path, colClasses = vapply(summary, class, "")),
                              summary, tolerance = 1e-12))) {
    stop("Stored summary differs from its validated seeds: ", summary_path)
  }
  message("Section case complete: ", case$case_id)
  cbind(case, summary, raw_switching_percent = summary$mean_switching_vertical_interior,
        seed_noise_percent = summary$mean_noise_2d, attributable_pp = summary$attributable_switching,
        interval_low_pp = summary$attributable_ci_low, interval_high_pp = summary$attributable_ci_high,
        kappa = summary$mean_cohens_kappa, noise_pairs = summary$n_noise_pairs,
        precision_limited = summary$mean_noise_2d >= 5)
}

publish_ledger <- function(completed, path) {
  ledger <- do.call(rbind, completed)
  if (!file.exists(path)) {
    publish_csv(ledger, path)
  } else if (!isTRUE(all.equal(utils::read.csv(path, colClasses = vapply(ledger, class, "")),
                              ledger, tolerance = 1e-12, check.attributes = FALSE))) {
    stop("Existing ledger differs from the verified case summaries: ", path)
  }
}
