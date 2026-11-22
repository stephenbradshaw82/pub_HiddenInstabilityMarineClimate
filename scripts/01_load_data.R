# 01 | Choose the input source
# "prepared" (default): use the supplied WA-derived GeoTIFFs. Run 03 next.
# "download": download raw NetCDFs from the providers. Run 02 next.
# Pick one option below. Choose a new run ID in 00 if you switch input sources.

input_source <- "prepared"  # Change to "download" to start from raw NetCDFs.

if (!exists("demo_dir", inherits = FALSE)) stop("Run 00_config.R first.")
source(file.path(demo_dir, "support", "raw_data.R"))
source(file.path(demo_dir, "support", "checkpoints.R"))

load_demo_inputs <- function(input_source) {
  if (!input_source %in% c("prepared", "download")) {
    stop('In script 01, choose "prepared" or "download".')
  }

  recipe_path <- file.path(run_dir, "input_source.rds")
  if (input_source == "prepared") {
    source(file.path(demo_dir, "support", "prepared_inputs.R"), local = TRUE)
    archives <- file.path(demo_dir, "data", "prepared", names(prepared_archive_groups()))
    missing <- archives[!file.exists(archives)]
    if (length(missing)) {
      stop("Missing prepared archive in data/prepared/: ", paste(basename(missing), collapse = ", "),
           ". Restore the missing file or choose the download option in script 01.")
    }
    if (!file.exists(recipe_path) && length(list.files(run_dir, all.files = TRUE, no.. = TRUE))) {
      stop("Choose a new run ID before using prepared inputs with this run folder.")
    }
  }
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  # The prepared installer checks the data manifest even if the ZIPs were repacked.
  if (file.exists(recipe_path)) {
    if (!identical(readRDS(recipe_path)$source, input_source)) {
      stop("Choose a new run ID before switching input sources.")
    }
  } else {
    check_stage_recipe(recipe_path, list(source = input_source))
  }

  if (input_source == "prepared") {
    prepared_dir <- file.path(run_dir, "prepared")
    install_prepared_inputs(archives, prepared_dir)
    message("Section complete: prepared WA inputs ready. Skip script 02 and run script 03.")
    return(prepared_dir)
  }

  source(file.path(demo_dir, "support", "download_netcdf.R"), local = TRUE)
  load_raw_inputs(raw_fixture_dir)
  message("Continue with script 02 to prepare the downloaded NetCDFs.")
  raw_fixture_dir
}

# Later stages read the region layers from the selected input folder.
fixture_dir <- load_demo_inputs(input_source)
