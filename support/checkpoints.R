check_stage_recipe <- function(path, recipe) {
  if (file.exists(path)) {
    if (!identical(readRDS(path), recipe)) stop("Inputs or code changed. Use a new run ID: ", path)
  } else {
    staged <- tempfile(".recipe-", tmpdir = dirname(path))
    saveRDS(recipe, staged)
    if (!file.rename(staged, path)) stop("Could not publish ", path)
  }
}

publish_csv <- function(table, path) {
  staged <- tempfile(".receipt-", tmpdir = dirname(path))
  utils::write.csv(table, staged, row.names = FALSE)
  if (file.exists(path) || !file.rename(staged, path)) stop("Could not publish ", path)
}

check_file_manifest <- function(path, directory = dirname(path)) {
  if (!file.exists(path)) stop("Missing stage receipt: ", path)
  receipt <- utils::read.csv(path)
  if (!nrow(receipt) || anyDuplicated(receipt$file) || anyNA(receipt$md5) ||
      !identical(unname(tools::md5sum(file.path(directory, receipt$file))), receipt$md5)) {
    stop("Files differ from the stage receipt: ", path)
  }
  invisible(receipt)
}
