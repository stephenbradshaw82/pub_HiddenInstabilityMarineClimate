# Shared package, NetCDF and monthly calendar checks.

check_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install these R packages first: ", paste(missing, collapse = ", "))
}

month_sequence <- function(first, last) {
  seq(as.Date(paste0(first, "-01")), as.Date(paste0(last, "-01")), by = "month")
}

read_cf_months <- function(values, units, calendar = "standard") {
  if (!calendar %in% c("standard", "gregorian", "proleptic_gregorian")) {
    stop("Unsupported source calendar: ", calendar)
  }
  if (!grepl("^days since [0-9]{4}-[0-9]{2}-[0-9]{2}", units) ||
      !all(is.finite(values))) {
    stop("Invalid CF time axis: ", units)
  }
  origin <- as.Date(substr(units, 12, 21))
  as.Date(format(origin + floor(values), "%Y-%m-01"))
}

check_months <- function(actual, expected, label) {
  if (!identical(as.numeric(actual), as.numeric(expected))) {
    stop(label, " has missing, duplicate or out-of-order months.")
  }
}

read_raw_axes <- function(nc, dataset) {
  variables <- c(BRAN = "temp", NOAA = "sst", CMIP6 = "thetao")
  axis_names <- list(BRAN = c("xt_ocean", "yt_ocean", "st_ocean", "Time"),
                     NOAA = c("lon", "lat", "time"),
                     CMIP6 = c("x", "y", "lev", "time"))
  variable <- variables[[dataset]]
  meta <- nc$var[[variable]]
  if (is.null(meta) ||
      !identical(unname(vapply(meta$dim, function(d) d$name, "")), axis_names[[dataset]])) {
    stop("Unexpected ", dataset, " temperature variable or axis order in ", nc$filename)
  }
  if (!meta$units %in% c("degrees C", "degree_C", "degC", "degree_Celsius", "degrees_Celsius")) {
    stop("Expected Celsius in ", nc$filename, ". Found units: ", meta$units)
  }
  meta$dim
}
