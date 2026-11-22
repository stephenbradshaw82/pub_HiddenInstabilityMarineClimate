# 00 | Set up R and configure the Gascoyne Coast analysis
# - Start a fresh R session. The demo was tested with R 4.6.1 on arm64 macOS.
# - This script installs missing R packages from CRAN and the pinned VoCC from GitHub.
# - Use a writable R package library and demo folder. Allow about 15 GiB for the run.
# - The optional raw download in script 01 needs HTTPS and NetCDF-4 with OPeNDAP.
# - Script 07 needs PNG graphics support.
# - For source installs, set up the platform's R compiler toolchain and development
#   libraries for NetCDF-C, libcurl, GDAL, GEOS, PROJ, SQLite, UDUNITS2 and TBB first.
# - Use the full DPIRD reporting polygon and a 100 km calculation buffer.
# - Report BRAN T1 and CNRM-CM6-1-HR SSP2-4.5 T1 and T2 on a common grid.
# - Analyse 42 levels from 5 to 210 m. Prepare 2.5 m separately for the surface map.
# - Choose a run ID here, then continue in the same R session.

if (file.exists("demo/scripts/00_config.R")) {
  demo_dir <- normalizePath("demo")
} else if (file.exists("scripts/00_config.R")) {
  demo_dir <- normalizePath(".")
} else {
  stop("Open the repository or standalone folder in RStudio first.")
}

# Script 01 writes temporary source data here and reuses it on reruns.
fixture_id <- "gascoyne_ssp245"
run_id <- "gascoyne_ssp245"
pilot_only <- FALSE       # TRUE runs the first seed of the named timing cases.
# Use the original ESGF NetCDF if remote subsetting fails. Each original is about 8.6 GB.
cmip6_full_file_fallback <- TRUE
if (!all(grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", c(fixture_id, run_id))) || run_id == "packages") {
  stop("Choose simple fixture and run IDs using letters, numbers, underscores or hyphens.")
}
raw_fixture_dir <- file.path(demo_dir, "data", "tmp", fixture_id)
fixture_dir <- raw_fixture_dir
run_dir <- file.path(demo_dir, "runs", run_id)
method_dir <- file.path(demo_dir, "method")
if (!dir.exists(method_dir)) stop("The bundled method folder is missing: ", method_dir)
if (getRversion() < "4.1.0") stop("Use R 4.1 or newer. The tested version is R 4.6.1.")

vocc_sha <- "0bfd857d837114eeb7004da588e6de1cb8092e3c"
cran <- "https://cloud.r-project.org"
cran_packages <- c(
  "ncdf4", "sf", "terra", "akima", "curl", "digest", # Inputs, transfer checks and interpolation
  "raster", "data.table", "geosphere", "dplyr",    # VoCC and model helpers
  "ggplot2", "cowplot", "remotes"                 # PNG figures and pinned VoCC installation
)
installed_sha <- suppressWarnings(utils::packageDescription("VoCC", fields = "RemoteSha"))
install_vocc <- !identical(installed_sha, vocc_sha) || !requireNamespace("VoCC", quietly = TRUE)
if (install_vocc) {
  if ("VoCC" %in% loadedNamespaces()) stop("Restart R, then rerun script 00 to install the pinned VoCC.")
}
missing <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) || install_vocc) {
  package_library <- .libPaths()[1L]
  if (file.access(package_library, 2L) != 0L) {
    user_libraries <- strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep, fixed = TRUE)[[1L]]
    if (!length(user_libraries) || !nzchar(user_libraries[1L])) {
      stop("Set R_LIBS_USER to a writable personal library and restart R.")
    }
    package_library <- path.expand(user_libraries[1L])
    dir.create(package_library, recursive = TRUE, showWarnings = FALSE)
    if (file.access(package_library, 2L) != 0L) stop("Cannot write to the R library: ", package_library)
    .libPaths(c(package_library, .libPaths()))
  }

  if (length(missing)) {
    message("Installing missing R packages: ", paste(missing, collapse = ", "))
    utils::install.packages(missing, lib = package_library, repos = cran)
    unavailable <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
    if (length(unavailable)) {
      stop("These packages could not be loaded: ", paste(unavailable, collapse = ", "),
           ". Check the installation log and system requirements, then rerun script 00.")
    }
  }
  if (install_vocc) {
    message("Installing VoCC at ", vocc_sha)
    remotes::install_github("JorGarMol/VoCC", ref = vocc_sha, repos = cran, auth_token = NULL,
                           lib = package_library, upgrade = "never", force = TRUE)
  }
}
if (!requireNamespace("VoCC", quietly = TRUE) ||
    !identical(utils::packageDescription("VoCC", fields = "RemoteSha"), vocc_sha)) {
  stop("The pinned VoCC could not be loaded. Check the installation log and rerun script 00.")
}

region_name <- "Gascoyne Coast"
buffer_m <- 100000
paper_bounds <- c(west = 108.75, east = 130, south = -39.25, north = -9)
# This rectangle covers the clipped buffer plus gradient neighbours. Its
# edges are on the paper's 0.50-degree grid so both grid origins are retained.
analysis_bounds <- c(west = 108.75, east = 116.75, south = -28.75, north = -15.25)
analysis_depths_m <- seq(5, 210, by = 5)
prepared_depths_m <- c(2.5, analysis_depths_m)
time_periods <- c(T1 = "19940101-20231231", T2 = "20240101-20531231")
# This fixture includes SSP2-4.5 only.
projection_scenario <- "ssp245"
figure_names <- c("01_context", "02_surface", "03_depth_profiles", "04_disagreement", "05_sensitivity")
matched_pilot_ids <- c("bran_primary", "cmip6_ssp245_t2_refined")
sensitivity_pilot_ids <- c("grid_050", "horizon_020yr", "depth_010m_aggregate")

message("Section complete: R packages ready and ", run_id, " configured.")
