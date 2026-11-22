# Shared velocity variant names and time-window checks.

velocity_variant_name <- function(grid_deg, gradient_floor, aggregation_m) {
  name <- sprintf("grid_%.2f_floor_q%03d", grid_deg, round(100 * gradient_floor))
  if (is.na(aggregation_m)) name else paste0(name, "_depthagg_", aggregation_m, "m")
}

verify_velocity_window <- function(path, period) {
  r <- terra::rast(path)
  months <- as.Date(terra::time(r))
  if (length(months) != terra::nlyr(r) || anyNA(months)) stop("Missing GeoTIFF dates: ", path)
  years <- substr(strsplit(period, "-", fixed = TRUE)[[1]], 1, 4)
  required <- month_sequence(paste0(years[1], "-01"), paste0(years[2], "-12"))
  check_months(months[months >= min(required) & months <= max(required)], required, basename(path))
}

