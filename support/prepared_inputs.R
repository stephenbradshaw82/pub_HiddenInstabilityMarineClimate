# Checks shared by the prepared and raw NetCDF routes.

prepared_stack_names <- function() {
  c(sprintf("SpatRaster_BRAN_depth_%05.1fm.tif", prepared_depths_m),
    sprintf("SpatRaster_CMIP6_ssp245_depth_%05.1fm.tif", prepared_depths_m),
    "SpatRaster_NOAA_depth_002.5m.tif")
}

prepared_archive_groups <- function() {
  stacks <- prepared_stack_names()
  cmip6 <- startsWith(stacks, "SpatRaster_CMIP6_")
  list(
    gascoyne_ssp245_wa_bran_noaa.zip = c(stacks[!cmip6], "temperature_manifest.csv",
                                        "regions.gpkg", "bundle_manifest.csv"),
    gascoyne_ssp245_wa_cmip6.zip = stacks[cmip6]
  )
}

install_prepared_inputs <- function(archives, prepared_dir) {
  stacks <- prepared_stack_names()
  inputs <- c(stacks, "temperature_manifest.csv", "regions.gpkg")
  groups <- prepared_archive_groups()
  for (i in seq_along(groups)) {
    entries <- utils::unzip(archives[i], list = TRUE)$Name
    if (anyDuplicated(entries) || !setequal(entries, groups[[i]])) {
      stop("Unexpected file list in the prepared ZIP: ", archives[i])
    }
  }

  bundle <- utils::read.csv(unz(archives[1L], "bundle_manifest.csv"))
  folder <- prepared_dir
  if (!dir.exists(prepared_dir)) {
    folder <- tempfile(".prepared-", tmpdir = dirname(prepared_dir))
    dir.create(folder)
    # Only the temporary extraction directory is removed on exit.
    on.exit(unlink(folder, recursive = TRUE), add = TRUE)
    for (archive in archives) utils::unzip(archive, exdir = folder)
  }

  if (!identical(utils::read.csv(file.path(folder, "bundle_manifest.csv")), bundle)) {
    stop("The prepared inputs changed. Use a new run ID.")
  }
  manifests <- c("temperature_manifest.csv", "bundle_manifest.csv")
  expected <- list(stacks, inputs)
  for (i in seq_along(manifests)) {
    path <- file.path(folder, manifests[i])
    if (!setequal(utils::read.csv(path)$file, expected[[i]])) stop("Unexpected files in ", path)
    check_file_manifest(path)
  }
  if (!identical(folder, prepared_dir)) {
    validate_prepared_temperature(folder, pilot = FALSE)
    regions <- file.path(folder, "regions.gpkg")
    if (!setequal(sf::st_layers(regions)$name, c("bioregions", "reporting", "simulation", "land")) ||
        !identical(terra::vect(regions, layer = "reporting")$bioregion, region_name) ||
        !nrow(terra::vect(regions, layer = "simulation"))) {
      stop("The prepared ZIP has unexpected region layers.")
    }
    if (!file.rename(folder, prepared_dir)) stop("Could not publish the checked prepared inputs.")
  }
  invisible(prepared_dir)
}

validate_prepared_temperature <- function(prepared_dir, pilot) {
  for (depth in prepared_depths_m) {
    suffix <- sprintf("_depth_%05.1fm.tif", depth)
    datasets <- c("BRAN", "CMIP6_ssp245")
    if (depth == 2.5) datasets <- c(datasets, "NOAA")
    masks <- vector("list", length(datasets))
    availability <- c(BRAN = "1993-01/2023-12", CMIP6_ssp245 = "1994-01/2053-12",
                      NOAA = "1981-09/2025-06")
    for (i in seq_along(datasets)) {
      path <- file.path(prepared_dir, paste0("SpatRaster_", datasets[i], suffix))
      r <- terra::rast(path)
      if (!isTRUE(all.equal(unname(as.vector(terra::ext(r))), unname(analysis_bounds))) ||
          !identical(terra::res(r), c(0.25, 0.25)) || !terra::is.lonlat(r)) {
        stop("Unexpected analysis grid: ", path)
      }
      window <- if (pilot) c("1994-01", "1994-12") else
        strsplit(availability[[datasets[i]]], "/", fixed = TRUE)[[1]]
      check_months(as.Date(terra::time(r)), month_sequence(window[1], window[2]), basename(path))
      finite <- is.finite(terra::values(r, mat = TRUE))
      masks[[i]] <- finite[, 1L]
      if (any(finite != masks[[i]])) stop("Coverage changes between months: ", path)
      if (!identical(masks[[1]], masks[[i]])) stop("The common masks differ at ", depth, " m.")
    }
  }
}

