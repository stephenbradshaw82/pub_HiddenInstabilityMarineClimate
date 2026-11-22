# Shared reporting and calculation grid operations.

region_mask <- function(polygon, grid) {
  values <- terra::values(terra::rasterize(polygon, grid, field = 1, touches = FALSE))
  matrix(is.finite(values), nrow = terra::nrow(grid), byrow = TRUE)
}

bundle_grid <- function(bundle) {
  terra::rast(nrows = length(bundle$y_centers), ncols = length(bundle$x_centers),
              extent = unlist(bundle$ext), crs = "EPSG:4326")
}

verify_grid_origin <- function(grid, degrees) {
  # Both resolutions start from the archived WA west and south edges.
  offset <- (c(terra::xmin(grid), terra::ymin(grid)) - paper_bounds[c("west", "south")]) / degrees
  if (any(abs(offset - round(offset)) > 1e-7) || any(abs(terra::res(grid) - degrees) > 1e-7)) {
    stop("Grid centres differ from the paper's ", degrees, "-degree grid.")
  }
}
