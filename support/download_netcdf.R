# Raw NetCDF sourcing for the download option in script 01.

source_inventory <- function() {
  bran_months <- month_sequence("1993-01", "2023-12")
  bran_names <- paste0("ocean_temp_mth_", format(bran_months, "%Y_%m"), ".nc")
  bran <- data.frame(dataset = "BRAN", experiment = "historical", first_month = bran_months,
                     last_month = bran_months, file = bran_names,
                     url = paste0("https://dapds00.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/", bran_names),
                     mirror = NA_character_)
  noaa <- data.frame(dataset = "NOAA", experiment = "observed",
                     first_month = as.Date("1981-09-01"), last_month = as.Date("2025-06-01"),
                     file = "sst.mon.mean.nc",
                     url = "https://psl.noaa.gov/thredds/dodsC/Datasets/noaa.oisst.v2.highres/sst.mon.mean.nc",
                     mirror = NA_character_)

  years <- seq(1990L, 2050L, by = 5L)
  experiment <- rep(c("historical", projection_scenario), c(5, 8))
  mip <- rep(c("CMIP", "ScenarioMIP"), c(5, 8))
  version <- rep(c("v20191021", "v20191202"), c(5, 8))
  cmip_names <- sprintf("thetao_Omon_CNRM-CM6-1-HR_%s_r1i1p1f2_gn_%d01-%d12.nc",
                        experiment, years, years + 4L)
  path <- paste(mip, "CNRM-CERFACS/CNRM-CM6-1-HR", experiment,
                "r1i1p1f2/Omon/thetao/gn", version, cmip_names, sep = "/")
  cmip <- data.frame(dataset = "CMIP6", experiment = experiment,
                     first_month = as.Date(sprintf("%d-01-01", years)),
                     last_month = as.Date(sprintf("%d-12-01", years + 4L)),
                     file = cmip_names,
                     url = paste0("https://esgf.nci.org.au/thredds/dodsC/esgcet/replica/CMIP6/", path),
                     mirror = paste0("https://esgf-node.ornl.gov/thredds/dodsC/css03_data/CMIP6/", path))
  plan <- rbind(bran, noaa, cmip)
  # ESGF's file records for the two pinned versions. Values are in file order.
  plan$archive_bytes <- NA_real_
  plan$archive_sha256 <- NA_character_
  cmip_rows <- plan$dataset == "CMIP6"
  plan$archive_bytes[cmip_rows] <- c(
    8638232236, 8640242123, 8638264540, 8638472310, 8635175300,
    8638565927, 8626892816, 8631749570, 8628630011, 8628841670,
    8631790369, 8632080838, 8627629292
  )
  plan$archive_sha256[cmip_rows] <- c(
    "4152a47d150162c1bf96351c29c3f00cfb18e1992522c4af7621a3fbe2872a84",
    "77ced573219ef9624ceb2a76249335d1256826e324320757ba7b9a8992cb1c3e",
    "46d7634677aa9f256fd48d87f62aec572467eaedb493f3714188db0a259ac202",
    "f5672596f36b0af4bf1221c866751cbb0036c22a0377f9e933cc51b429258b20",
    "fa73dd45fd5bf96a27e7df2bef8d1d23772c27db146ab90065bc7541c8d87138",
    "42ddde093990a4c0c1b4279d124561f59602fde6644d36f7e5f47d47711571df",
    "f87ad3703e910bc287064e848a85c3d840c13a783ccb7f33bec196a6ca588c40",
    "43c84a6f28c2e635dbc5cc4865be609565dcb1507d42f40d86610ca9679af7e9",
    "bef8eca8cf5d86a503d8316e4eee52a4c041a0d94284341a883a1c7de10cdbe8",
    "3d742296167cd6ed4c3e5a5672563c4aafbbe4661b95d9488c119c804c5ed4a7",
    "125d4ea8cf59e6552187daf2309cbd1df3dc78a9b1917b420dafadefac3ea801",
    "6bcda920b5251561fa8134fc35651263bc120d4f518ae3ce9bb9abafb2833f55",
    "e1b2a4ba7ce368dc62e3e91fb2d6579a9567282634d6f545e9754745f02006c5"
  )
  plan
}

save_source_slice <- function(row, input_path, output, checkpoint_dir = NULL, url = input_path) {
  input <- tryCatch(ncdf4::nc_open(input_path), error = identity)
  if (inherits(input, "error")) {
    message("Could not open ", url, ". Retrying in 15 seconds.")
    Sys.sleep(15)
    input <- ncdf4::nc_open(input_path)
  }
  on.exit(ncdf4::nc_close(input), add = TRUE)
  axes <- read_raw_axes(input, row$dataset)
  variable <- c(BRAN = "temp", NOAA = "sst", CMIP6 = "thetao")[[row$dataset]]
  # These one-based indices span the reference Gascoyne crop on each native grid.
  # The depth levels bracket 210 m and the surrounding grid cells support regridding.
  indices <- switch(row$dataset,
                    BRAN = list(1090:1179, 452:609, 1:26),
                    NOAA = list(436:472, 241:304),
                    CMIP6 = list(140:182, 403:471, 1:32))
  time <- axes[[length(axes)]]
  calendar <- ncdf4::ncatt_get(input, time$name, "calendar")
  calendar <- if (calendar$hasatt) tolower(calendar$value) else "standard"
  months <- read_cf_months(time$vals, time$units, calendar)
  requested <- seq(row$first_month, row$last_month, by = "month")
  time_indices <- which(months >= min(requested) & months <= max(requested))
  check_months(months[time_indices], requested, row$file)
  indices <- c(indices, list(time_indices))
  start <- vapply(indices, min, integer(1))
  count <- lengths(indices)
  if (any(!vapply(indices, function(i) all(diff(i) == 1L), logical(1)))) {
    stop("The expected slice is no longer contiguous.")
  }
  selected <- Map(function(axis, index) as.numeric(axis$vals[index]), axes, indices)
  # Catch a shifted source grid before it can change the area or the layer depths.
  limits <- switch(row$dataset,
                   BRAN = c(108.95, 117.85, -29.85, -14.15, 2.5, 217.054489),
                   NOAA = c(108.875, 117.875, -29.875, -14.125),
                   CMIP6 = c(140, 182, 403, 471, 0.505760, 221.141180))
  edges <- unlist(lapply(head(selected, -1), function(axis) axis[c(1L, length(axis))]))
  if (any(!is.finite(edges)) || any(abs(edges - limits) > 1e-4)) stop("The native grid changed.")

  dims <- Map(function(axis, values) ncdf4::ncdim_def(axis$name, axis$units, values), axes, selected)
  meta <- input$var[[variable]]
  definitions <- list(ncdf4::ncvar_def(variable, meta$units, dims, missval = meta$missval,
                                      longname = meta$longname, prec = "float", compression = 4))
  tracking <- ncdf4::ncatt_get(input, 0, "tracking_id")
  if (row$dataset == "CMIP6") {
    if (!tracking$hasatt || !nzchar(tracking$value)) stop("The CMIP6 source has no tracking ID.")
    definitions <- c(definitions, list(
      ncdf4::ncvar_def("lon", "degrees_east", dims[1:2], prec = "double"),
      ncdf4::ncvar_def("lat", "degrees_north", dims[1:2], prec = "double")
    ))
  }
  source_urls <- url
  local <- ncdf4::nc_create(output, definitions, force_v4 = TRUE)
  tryCatch({
    # THREDDS can reject a long CMIP6 request while still serving smaller slices.
    time_dim <- length(count)
    batch_months <- if (row$dataset != "CMIP6") count[time_dim] else
      if (endsWith(url, "#mode=bytes")) 1L else 3L
    first <- 1L
    has_water <- FALSE
    retried <- FALSE
    while (first <= count[time_dim]) {
      read_start <- start
      read_start[time_dim] <- start[time_dim] + first - 1L
      batch <- count
      batch[time_dim] <- min(batch_months, count[time_dim] - first + 1L)
      checkpoint <- if (!is.null(checkpoint_dir)) {
        file.path(checkpoint_dir, paste0(format(requested[first], "%Y-%m"), ".rds"))
      }
      if (length(checkpoint) && file.exists(checkpoint)) {
        saved_batch <- readRDS(checkpoint)
        batch <- dim(saved_batch$values)
        if (!identical(saved_batch$tracking_id, tracking$value) ||
            !identical(saved_batch$start, read_start) ||
            length(batch) != time_dim || !identical(head(batch, -1L), head(unname(count), -1L)) ||
            batch[time_dim] > count[time_dim] - first + 1L) {
          stop("The saved CMIP6 block differs from the pinned source: ", checkpoint)
        }
        values <- saved_batch$values
        source_urls <- union(source_urls, saved_batch$url)
      } else {
        if (row$dataset == "CMIP6" && (first - 1L) %% 12L == 0L) {
          message("  Fetching from ", requested[first], " (", first, "/", count[time_dim], " months)")
        }
        values <- tryCatch(ncdf4::ncvar_get(input, variable, read_start, batch, collapse_degen = FALSE),
                           error = identity)
        if (inherits(values, "error")) {
          if (row$dataset != "CMIP6") stop(values)
          if (batch[time_dim] > 1L) {
            batch_months <- 1L
            message("Retrying ", row$file, " one month at a time from ", requested[first], ".")
          } else {
            if (retried) stop(values)
            retried <- TRUE
            message("Retrying ", requested[first], " in 10 seconds.")
            Sys.sleep(10)
          }
          next
        }
        if (!identical(dim(values), unname(batch))) stop("The remote slice has the wrong dimensions.")
        if (length(checkpoint)) {
          # Publish each block before asking the server for the next one.
          staged_batch <- file.path(checkpoint_dir, ".pending.rds")
          saveRDS(list(tracking_id = tracking$value, start = read_start, url = url, values = values), staged_batch)
          if (!file.rename(staged_batch, checkpoint)) stop("Could not save the CMIP6 block: ", checkpoint)
        }
      }
      has_water <- has_water || any(is.finite(values))
      write_start <- c(rep(1L, time_dim - 1L), first)
      ncdf4::ncvar_put(local, variable, values, write_start, batch)
      first <- first + batch[time_dim]
      retried <- FALSE
    }
    if (!has_water) stop("The remote temperature slice has no finite values.")
    ncdf4::ncatt_put(local, time$name, "calendar", calendar)
    ncdf4::ncatt_put(local, 0, "source_url", paste(source_urls, collapse = "; "))
    ncdf4::ncatt_put(local, 0, "source_start", start)
    ncdf4::ncatt_put(local, 0, "source_count", count)
    for (name in c("tracking_id", "license", "institution", "title")) {
      attribute <- ncdf4::ncatt_get(input, 0, name)
      if (attribute$hasatt) ncdf4::ncatt_put(local, 0, name, attribute$value)
    }
    if (row$dataset == "CMIP6") {
      for (coordinate in c("lon", "lat")) {
        values <- ncdf4::ncvar_get(input, coordinate, start[1:2], count[1:2])
        if (!all(is.finite(values))) stop("Missing native CMIP6 coordinates.")
        ncdf4::ncvar_put(local, coordinate, values)
      }
      ncdf4::ncatt_put(local, variable, "coordinates", "lon lat")
    }
  }, finally = ncdf4::nc_close(local))

  saved <- ncdf4::nc_open(output)
  on.exit(ncdf4::nc_close(saved), add = TRUE)
  saved_axes <- read_raw_axes(saved, row$dataset)
  saved_coordinates <- lapply(saved_axes, function(axis) as.numeric(axis$vals))
  if (!identical(unname(saved_coordinates), unname(selected))) {
    stop("The saved NetCDF axes differ from the source slice.")
  }
  invisible(source_urls)
}

download_cmip6_original <- function(row, urls, path) {
  check_packages(c("curl", "digest"))
  if (isTRUE(file.info(path)$size > row$archive_bytes)) {
    stop("The partial original exceeds the ESGF file size. Move it aside and rerun: ", path)
  }
  url <- urls[1L]
  if (!isTRUE(file.info(path)$size == row$archive_bytes)) {
    message("Downloading the original ", row$file, " (",
            sprintf("%.2f", row$archive_bytes / 1e9), " GB) for local cropping.")
    for (url in urls) {
      transfer <- curl::multi_download(url, path, resume = TRUE, failonerror = TRUE,
                                       connecttimeout = 20L, low_speed_limit = 32768L, low_speed_time = 120L)
      if (is.na(transfer$success)) stop("Transfer interrupted. Rerun script 01 to resume the original NetCDF.")
      if (isTRUE(file.info(path)$size == row$archive_bytes)) break
      if (!transfer$success) {
        message("HTTPS transfer failed at ", url, ": ", transfer$error, ". Keeping the partial file.")
      }
    }
  }
  if (!isTRUE(file.info(path)$size == row$archive_bytes)) {
    stop("The original NetCDF is incomplete: ", path, ". Rerun script 01 to resume.")
  }
  if (!identical(digest::digest(path, algo = "sha256", file = TRUE), row$archive_sha256)) {
    stop("The original NetCDF differs from the ESGF checksum. Move it aside and rerun: ", path)
  }
  url
}

download_regions <- function(destination, work) {
  dpird_url <- paste0(
    "https://public-services.slip.wa.gov.au/public/rest/services/SLIP_Public_Services/Boundaries/MapServer/51/query",
    "?where=1%3D1&outFields=bioregion&returnGeometry=true&outSR=4326&f=geojson"
  )
  coast_url <- "https://naturalearth.s3.amazonaws.com/50m_cultural/ne_50m_admin_0_countries.zip"
  bioregion_file <- file.path(work, "bioregions.geojson")
  coast_zip <- file.path(work, "countries.zip")
  utils::download.file(dpird_url, bioregion_file, mode = "wb")
  utils::download.file(coast_url, coast_zip, mode = "wb")
  bioregions <- sf::st_read(bioregion_file, quiet = TRUE)
  if (nrow(bioregions) != 6L || !"bioregion" %in% names(bioregions)) stop("Unexpected DPIRD bioregions.")
  report <- bioregions[bioregions$bioregion == region_name, ]
  if (nrow(report) != 1L) stop("Expected one DPIRD reporting polygon for ", region_name)
  report_albers <- sf::st_make_valid(sf::st_transform(report, 3577))
  buffered <- sf::st_transform(sf::st_buffer(sf::st_union(report_albers), buffer_m), 4326)
  wa <- sf::st_bbox(setNames(paper_bounds, c("xmin", "xmax", "ymin", "ymax")), crs = sf::st_crs(4326))
  old_s2 <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  calculation <- suppressWarnings(sf::st_crop(sf::st_sf(geometry = buffered), wa))
  utils::unzip(coast_zip, exdir = work)
  countries <- sf::st_read(file.path(work, "ne_50m_admin_0_countries.shp"), quiet = TRUE)
  land <- countries[countries$ADMIN == "Australia", ]
  if (!nrow(calculation) || nrow(land) != 1L) stop("Missing calculation polygon or coastline.")
  staged <- file.path(work, "regions.gpkg")
  sf::st_write(bioregions, staged, layer = "bioregions", quiet = TRUE)
  sf::st_write(report, staged, layer = "reporting", quiet = TRUE)
  sf::st_write(calculation, staged, layer = "simulation", quiet = TRUE)
  sf::st_write(land, staged, layer = "land", quiet = TRUE)
  if (!file.rename(staged, destination)) stop("Could not save the region layers.")
}

download_inputs <- function(plan, fixture_dir) {
  dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(fixture_dir, ".downloading")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Input download already running or locked: ", lock)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  receipts <- file.path(fixture_dir, ".receipts")
  dir.create(receipts, showWarnings = FALSE)

  records <- vector("list", nrow(plan))
  for (i in seq_len(nrow(plan))) {
    row <- plan[i, ]
    output <- file.path(fixture_dir, row$file)
    receipt <- file.path(receipts, paste0(row$file, ".rds"))
    if (file.exists(output) && file.exists(receipt)) {
      record <- readRDS(receipt)
      if (identical(unname(tools::md5sum(output)), record$md5)) {
        records[[i]] <- record
        next
      }
    }
    if (file.exists(output)) stop("Move the unverified file aside and rerun script 01: ", output)

    message("Input ", i, "/", nrow(plan), "  ", row$file)
    staged <- file.path(lock, "slice.nc")
    started <- proc.time()[["elapsed"]]
    urls <- na.omit(c(row$url, row$mirror))
    checkpoint_dir <- NULL
    if (row$dataset == "CMIP6") {
      checkpoint_dir <- file.path(fixture_dir, ".cmip6_chunks", row$file)
      dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
      # These HTTPS services expose the same pinned NetCDF without the DAP server.
      nci_archive <- sub("/dodsC/", "/fileServer/", row$url, fixed = TRUE)
      dkrz_archive <- sub("esgf.nci.org.au/thredds/fileServer/esgcet/replica/CMIP6/",
                          "esgf3.dkrz.de/thredds/fileServer/cmip6/", nci_archive, fixed = TRUE)
      archive_urls <- c(nci_archive, dkrz_archive, sub("/dodsC/", "/fileServer/", row$mirror, fixed = TRUE))
      urls <- c(urls, paste0(c(nci_archive, dkrz_archive), "#mode=bytes"))
      original <- file.path(checkpoint_dir, "original.nc")
      # Resume a started whole-file transfer before opening the subset services again.
      if (file.exists(original) && cmip6_full_file_fallback) urls <- character()
    }
    failure <- "Incomplete CMIP6 input"
    for (url in urls) {
      failure <- tryCatch({
        source_urls <- save_source_slice(row, url, staged, checkpoint_dir)
        NULL
      }, error = conditionMessage)
      if (is.null(failure)) break
      message("Source failed at ", url, "\n", failure)
      unlink(staged)
    }
    if (!is.null(failure) && row$dataset == "CMIP6" && cmip6_full_file_fallback) {
      url <- download_cmip6_original(row, archive_urls, original)
      source_urls <- save_source_slice(row, original, staged, checkpoint_dir, url = url)
      failure <- NULL
    }
    if (!is.null(failure)) stop("Download stopped at ", row$file, ". Rerun 01 to resume from the saved inputs.")
    record <- data.frame(row[c("dataset", "experiment", "first_month", "last_month", "file")],
                         source_url = paste(unique(c(url, source_urls)), collapse = "; "),
                         bytes = file.info(staged)$size,
                         md5 = unname(tools::md5sum(staged)), pilot = FALSE,
                         elapsed_s = proc.time()[["elapsed"]] - started)
    saveRDS(record, receipt)
    if (!file.rename(staged, output)) stop("Could not save ", output)
    # The completed regional NetCDF and its receipt replace these temporary blocks.
    if (length(checkpoint_dir)) unlink(checkpoint_dir, recursive = TRUE)
    records[[i]] <- record
  }
  region_path <- file.path(fixture_dir, "regions.gpkg")
  if (!file.exists(region_path)) download_regions(region_path, lock)
  inventory <- do.call(rbind, records)
  staged <- file.path(lock, "sources.csv")
  utils::write.csv(inventory, staged, row.names = FALSE)
  if (!file.rename(staged, file.path(fixture_dir, "sources.csv"))) stop("Could not save sources.csv.")
  inventory
}

load_raw_inputs <- function(fixture_dir) {
  check_packages(c("ncdf4", "terra", "sf"))
  started <- proc.time()[["elapsed"]]
  plan <- source_inventory()
  manifest <- file.path(fixture_dir, "sources.csv")
  region_path <- file.path(fixture_dir, "regions.gpkg")
  raw_files <- file.path(fixture_dir, plan$file)
  if (!all(file.exists(c(manifest, region_path, raw_files)))) {
    download_inputs(plan, fixture_dir)
  }
  inventory <- utils::read.csv(manifest)
  columns <- c("dataset", "experiment", "first_month", "last_month", "file")
  if (!all(c(columns, "md5", "pilot", "source_url") %in% names(inventory)) ||
      !identical(lapply(inventory[columns], as.character), lapply(plan[columns], as.character)) ||
      anyNA(inventory$md5) || !all(grepl("^[0-9a-f]{32}$", inventory$md5))) {
    stop("The source inventory differs from the configured download: ", manifest)
  }
  actual <- unname(tools::md5sum(raw_files))
  invalid <- is.na(actual) | actual != inventory$md5
  if (any(invalid)) {
    stop("Input checksum failed. Move the named file aside and rerun script 01: ", raw_files[which(invalid)[1L]])
  }
  reporting_region <- terra::vect(region_path, layer = "reporting")
  simulation_region <- terra::vect(region_path, layer = "simulation")
  if (!identical(reporting_region$bioregion, region_name) || !nrow(simulation_region)) {
    stop("Missing or mismatched Gascoyne reporting and calculation polygons.")
  }
  elapsed <- proc.time()[["elapsed"]] - started
  message("Section complete: validated ", nrow(inventory), " raw NetCDFs in ", fixture_dir,
          " (", sprintf("%.1f", elapsed), " s this session).")
  inventory
}

