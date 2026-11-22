# 07 | Draw the Gascoyne Coast figures
# - Map the four WA marine bioregions, the calculation buffer and BRAN surface HVoCC.
# - Compare native VoCC 2D, matched 2D, matched 3D and benchmark depth profiles.
# - Map count-weighted matched class disagreement and plot the 37 sensitivity results.
# - Save five PNGs, the two result tables and the data used in the figures.

if (!exists("demo_dir", inherits = FALSE)) stop("Run 00_config.R first.")
for (file in c("raw_data.R", "region.R", "checkpoints.R")) {
  source(file.path(demo_dir, "support", file))
}

save_demo_png <- function(plot, name, directory, width, height) {
  ggplot2::ggsave(file.path(directory, paste0(name, ".png")), plot = plot,
                  width = width, height = height, dpi = 300, bg = "white")
}

replace_output_folder <- function(staged, destination, history) {
  previous <- NULL
  if (dir.exists(destination)) {
    dir.create(history, recursive = TRUE, showWarnings = FALSE)
    previous <- tempfile(paste0(basename(destination), "-"), tmpdir = history)
    if (!file.rename(destination, previous)) stop("Could not preserve ", destination)
  }
  if (!file.rename(staged, destination)) {
    if (!is.null(previous) && !file.rename(previous, destination)) {
      stop("The previous output is preserved at ", previous)
    }
    stop("The completed staging folder is at ", staged)
  }
}

publish_demo_outputs <- function(figure_dir, run_dir, output_dir) {
  images <- paste0(figure_names, ".png")
  output_files <- c(images, "matched_results.csv", "sensitivity_results.csv")
  inputs <- c(file.path(figure_dir, images),
              file.path(run_dir, "cases", c("matched_ledger.csv", "sensitivity_ledger.csv")))
  if (dir.exists(output_dir)) {
    previous <- check_file_manifest(file.path(output_dir, "output_manifest.csv"))
    existing <- setdiff(list.files(output_dir, all.files = TRUE, no.. = TRUE), ".DS_Store")
    if (!setequal(previous$file, output_files) ||
        !setequal(existing, c(previous$file, "output_manifest.csv"))) {
      stop("Move personal files out of the display folder before publishing: ", output_dir)
    }
  }
  staged <- tempfile(".outputs-", tmpdir = dirname(output_dir))
  if (!dir.create(staged) || !all(file.copy(inputs, file.path(staged, output_files)))) {
    stop("Could not stage the figures and result tables in ", staged)
  }
  md5 <- unname(tools::md5sum(file.path(staged, output_files)))
  if (anyNA(md5) || !identical(md5, unname(tools::md5sum(inputs)))) stop("Staged files changed.")
  publish_csv(data.frame(run_id = basename(run_dir), file = output_files, md5 = md5),
              file.path(staged, "output_manifest.csv"))
  replace_output_folder(staged, output_dir, file.path(run_dir, "figure_history"))
}

surface_velocity_arrows <- function(surface, bearings) {
  if (!terra::compareGeom(surface, bearings, stopOnError = FALSE)) stop("Surface grids differ.")
  names(surface) <- "speed_km_year"
  names(bearings) <- "bearing_deg"
  cells <- terra::as.data.frame(c(surface, bearings), xy = TRUE, na.rm = TRUE)
  names(cells)[1:2] <- c("longitude", "latitude")
  cells <- subset(cells, is.finite(bearing_deg) & is.finite(speed_km_year) & speed_km_year > 0)
  if (!nrow(cells)) stop("No finite positive surface speeds with bearings in the reporting region.")
  # Show one vector per degree square. Use the cell nearest its centre.
  bin_x <- floor(cells$longitude)
  bin_y <- floor(cells$latitude)
  distance <- (cells$longitude - bin_x - 0.5)^2 + (cells$latitude - bin_y - 0.5)^2
  order_cells <- order(distance, cells$longitude, cells$latitude)
  bins <- paste(bin_x[order_cells], bin_y[order_cells])
  selected <- cells[order_cells[!duplicated(bins)], ]
  # Keep all bearings legible. Colour alone encodes speed.
  half_length <- rep(18000, nrow(selected))
  centres <- as.matrix(selected[c("longitude", "latitude")])
  start <- geosphere::destPoint(centres, (selected$bearing_deg + 180) %% 360, half_length)
  end <- geosphere::destPoint(centres, selected$bearing_deg, half_length)
  transform(selected, x = start[, 1], y = start[, 2], xend = end[, 1], yend = end[, 2])
}

report_profile_plot <- function(profiles, trajectory_palette, benchmark_palette, caption) {
  profiles <- subset(profiles, n_eligible > 0 & n_cells > 0)
  method_order <- c("Native VoCC 2D", "Matched 2D", "Matched 3D", "Benchmark")
  column_order <- c("BRAN T1", "CMIP6 SSP2-4.5 T1", "CMIP6 SSP2-4.5 T2")
  profiles$method_label <- factor(profiles$method_label, levels = method_order)
  profiles$column_label <- factor(profiles$column_label, levels = column_order)
  if (!"Unclassified" %in% profiles$class_label) {
    trajectory_palette <- trajectory_palette[names(trajectory_palette) != "Unclassified"]
  }

  profile_block <- function(data, palette, legend_title, show_x) {
    if (!all(data$class_label %in% names(palette))) stop("A profile class has no colour.")
    data$class_label <- factor(data$class_label, levels = names(palette))
    ggplot2::ggplot(data, ggplot2::aes(pct_eligible, depth_m, fill = class_label)) +
      ggplot2::geom_col(orientation = "y", width = 4.2, show.legend = TRUE) +
      ggplot2::facet_grid(method_label ~ column_label) +
      ggplot2::scale_y_reverse(breaks = seq(25, 200, 25), limits = c(212.5, 2.5),
                              expand = ggplot2::expansion(mult = 0)) +
      ggplot2::scale_x_continuous(breaks = seq(0, 100, 25),
                                  expand = ggplot2::expansion(mult = 0, add = c(1, 2))) +
      ggplot2::scale_fill_manual(values = palette, limits = names(palette), drop = FALSE,
                                 name = legend_title) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        strip.text.x = if (show_x) ggplot2::element_blank() else ggplot2::element_text(face = "bold"),
        strip.text.y = ggplot2::element_text(face = "bold"),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "grey65", fill = NA, linewidth = 0.35),
        panel.spacing.x = grid::unit(0.14, "in"),
        axis.text.x = if (show_x) ggplot2::element_text() else ggplot2::element_blank(),
        legend.title = ggplot2::element_text(face = "bold"),
        legend.position = "right",
        legend.key.height = grid::unit(0.19, "in"),
        legend.key.width = grid::unit(0.19, "in"),
        plot.margin = ggplot2::margin(2, 5, 2, 2)
      )
  }

  trajectory <- profile_block(subset(profiles, method_label != "Benchmark"),
                              trajectory_palette, "2D / 3D classes", FALSE)
  benchmark <- profile_block(subset(profiles, method_label == "Benchmark"),
                             benchmark_palette, "Benchmark classes", TRUE)
  combined <- cowplot::plot_grid(trajectory, benchmark, ncol = 1, align = "v", axis = "lr",
                                  rel_heights = c(3, 1.06))
  body <- cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label("Depth (m)", angle = 90),
                             combined, nrow = 1, rel_widths = c(0.035, 0.965))
  title <- cowplot::ggdraw() +
    cowplot::draw_label("Gascoyne Coast trajectory and benchmark profiles", x = 0.035, y = 0.68,
                        hjust = 0, fontface = "bold", size = 18) +
    cowplot::draw_label("Matched panels pool 20 paired seeds", x = 0.035, y = 0.22, hjust = 0, size = 11)
  footer <- cowplot::ggdraw() +
    cowplot::draw_label("Eligible observations (%)", x = 0.43, y = 0.76, size = 12) +
    cowplot::draw_label(caption, x = 0.985, y = 0.23, hjust = 1, size = 9, colour = "grey35")
  cowplot::plot_grid(title, body, footer, ncol = 1, rel_heights = c(0.065, 0.88, 0.055))
}

external_profile_rows <- function(stack_path, fields, eligible, template, lookup, method, period,
                                  count_classes) {
  stack <- terra::crop(terra::rast(stack_path), template)
  if (!terra::compareGeom(stack, template, stopOnError = FALSE)) {
    stop("Comparator and matched grids differ: ", stack_path)
  }
  index <- match(paste0("depth_", fields$depths), names(stack))
  if (anyNA(index)) stop("Comparator is missing required depth layers: ", stack_path)
  rows <- lapply(seq_along(index), function(k) {
    classes <- matrix(terra::values(stack[[index[k]]]), nrow = fields$nrow, byrow = TRUE)
    count_classes(classes, eligible[, , k], lookup, region_name, "Bioregion",
                  fields$depths[k], method, period, "one_field")
  })
  do.call(rbind, rows)
}

check_packages(c("ggplot2", "cowplot", "terra", "sf", "geosphere"))
library(ggplot2)
for (file in c("00_src/revision_final_figure_functions.R", "01_scripts/610_functions_3D_fast.R")) {
  source(file.path(method_dir, file))
}

case_dir <- file.path(run_dir, "cases")
ledgers <- file.path(case_dir, c("matched_ledger.csv", "sensitivity_ledger.csv"))
if (!all(file.exists(ledgers))) stop("Finish 05_run_matched_models.R and 06_run_sensitivity.R first.")
matched <- read.csv(ledgers[1])
sensitivity <- read.csv(ledgers[2])
if (nrow(matched) != 3L || any(matched$n_seeds != 20L | matched$noise_pairs != 10L) ||
    nrow(sensitivity) != 37L || any(sensitivity$n_seeds != 5L | sensitivity$noise_pairs != 2L) ||
    anyDuplicated(c(matched$case_id, sensitivity$case_id))) {
  stop("Expected three finished headline cases and 37 finished sensitivity cases.")
}
out <- tempfile(".figures-", tmpdir = run_dir)
if (!dir.create(out)) stop("Could not create the figure staging folder.")
caption <- paste("Gascoyne Coast replay ·", run_id)
region_path <- file.path(fixture_dir, "regions.gpkg")
reporting_region <- terra::vect(region_path, layer = "reporting")
bioregions <- sf::st_read(region_path, layer = "bioregions", quiet = TRUE)
simulation <- sf::st_read(region_path, layer = "simulation", quiet = TRUE)
land <- sf::st_read(region_path, layer = "land", quiet = TRUE)
report_bbox <- sf::st_bbox(sf::st_as_sf(reporting_region))
map_extent <- coord_sf(xlim = unname(report_bbox[c("xmin", "xmax")]) + c(-0.25, 0.25),
                       ylim = unname(report_bbox[c("ymin", "ymax")]) + c(-0.25, 0.25), expand = FALSE)

# 01. Locate the calculation area on the same WA extent used in the paper.
bioregions$bioregion <- as.character(bioregions$bioregion)
region_palette <- c("Gascoyne Coast" = "#e9ae53", "North Coast" = "#c1d8d4",
                    "West Coast" = "#b7cbdc", "South Coast" = "#d0c5d9")
bioregions <- subset(bioregions, bioregion %in% names(region_palette))
if (!setequal(unique(bioregions$bioregion), names(region_palette))) stop("Unexpected DPIRD bioregions.")
context <- ggplot() +
  geom_sf(data = land, fill = "grey93", colour = "grey55", linewidth = 0.25) +
  geom_sf(data = bioregions, aes(fill = bioregion), colour = "grey48", linewidth = 0.28) +
  geom_sf(data = subset(bioregions, bioregion == region_name), fill = NA,
          colour = "#9b460d", linewidth = 0.65) +
  geom_sf(data = simulation, aes(linetype = "Calculation area"), fill = NA,
          colour = "#463716", linewidth = 0.62) +
  scale_fill_manual(values = region_palette) +
  scale_linetype_manual(values = c("Calculation area" = "dashed"), name = NULL) +
  coord_sf(xlim = unname(paper_bounds[c("west", "east")]),
           ylim = unname(paper_bounds[c("south", "north")]), expand = FALSE) +
  theme_bw(base_size = 11) +
  labs(title = "Gascoyne Coast marine climate demo", fill = "DPIRD bioregion",
       x = "Longitude (°E)", y = "Latitude",
       caption = "The model uses the reporting region and a 100 km buffer clipped to the WA model grid.")
save_demo_png(context, "01_context", out, 8, 7)
write.csv(data.frame(file = "regions.gpkg", md5 = unname(tools::md5sum(region_path)),
                     buffer_m = buffer_m), file.path(out, "context_inputs.csv"), row.names = FALSE)

# 02. Colour shows local BRAN speed and arrows show bearings at selected cells.
comparator_root <- file.path(run_dir, "comparators")
surface_root <- file.path(comparator_root, "BRAN", "T1")
check_file_manifest(file.path(surface_root, "files.csv"))
surface <- terra::rast(file.path(surface_root, "surface_speed.tif"))
bearings <- terra::rast(file.path(surface_root, "surface_bearing.tif"))
surface_mask <- terra::rasterize(reporting_region, surface, field = 1, touches = FALSE)
surface <- terra::mask(surface, surface_mask)
bearings <- terra::mask(bearings, surface_mask)
speed <- terra::as.data.frame(surface, xy = TRUE, na.rm = TRUE)
names(speed) <- c("longitude", "latitude", "speed_km_year")
vectors <- surface_velocity_arrows(surface, bearings)
write.csv(speed, file.path(out, "surface_speed.csv"), row.names = FALSE)
write.csv(vectors, file.path(out, "surface_velocity_arrows.csv"), row.names = FALSE)
surface_plot <- ggplot(speed, aes(longitude, latitude)) +
  geom_raster(aes(fill = speed_km_year)) +
  geom_sf(data = land, fill = "grey93", colour = "grey60", linewidth = 0.25, inherit.aes = FALSE) +
  geom_segment(data = vectors, aes(x, y, xend = xend, yend = yend), inherit.aes = FALSE,
               arrow = grid::arrow(length = grid::unit(1.7, "mm"), type = "closed"),
               colour = "white", linewidth = 0.65) +
  geom_segment(data = vectors, aes(x, y, xend = xend, yend = yend), inherit.aes = FALSE,
               arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "closed"),
               colour = "#222222", linewidth = 0.35) +
  scale_fill_viridis_c(name = "km/year", option = "plasma", begin = 0.18, end = 0.99) +
  map_extent + theme_bw(base_size = 11) +
  labs(title = "BRAN surface climate velocity", subtitle = "1994–2023 · 2.5 m depth",
       x = "Longitude (°E)", y = "Latitude",
       caption = "Arrows show direction at up to one grid cell per degree square. Colour shows speed.")
save_demo_png(surface_plot, "02_surface", out, 7.8, 7.8)

# 03 and 04. Use the paper's profile and switch counts inside the DPIRD polygon.
native_lookup <- data.frame(class_id = 1:9, class_label = class_levels_3d_fast(FALSE))
benchmark_palette <- c(Stable = "#d9f0a3", "Upward-dominant" = "#74add1",
                       "Downward-dominant" = "#4575b4", "Horizontal-dominant" = "#fdae61",
                       "Horizontal + upward" = "#e78ac3", "Horizontal + downward" = "#d73027")
benchmark_lookup <- data.frame(class_id = 0:5, class_label = names(benchmark_palette))
profiles <- vector("list", nrow(matched))
for (i in seq_len(nrow(matched))) {
  case <- matched[i, ]
  dataset <- if (case$dataset == "BRAN") "BRAN" else "CMIP6_ssp245"
  period <- names(time_periods)[match(case$time_period, time_periods)]
  column <- paste(if (dataset == "BRAN") "BRAN" else "CMIP6 SSP2-4.5", period)
  root <- file.path(case_dir, case$output_group, paste0("case_", case$case_id))
  contract <- connectivity_final_figure_validate_case(root, 20L, case$particles_per_voxel)
  fields <- contract$template$fields
  template <- bundle_grid(fields)
  reporting_mask <- region_mask(reporting_region, template)
  eligible <- contract$template$classification$vertical_interior_mask & fields$valid_mask &
    array(rep(reporting_mask, length(fields$depths)), dim(fields$valid_mask))
  zones <- setNames(list(list(mask = reporting_mask, zone_type = "Bioregion")), region_name)
  matched_rows <- connectivity_final_figure_aggregate_profiles(contract, zones, column)
  comparator <- file.path(comparator_root, dataset, period)
  check_file_manifest(file.path(comparator, "files.csv"))
  native_rows <- external_profile_rows(file.path(comparator, "native_class_ids.tif"),
                                       fields, eligible, template, native_lookup,
                                       "Native VoCC 2D", column, connectivity_final_figure_class_rows)
  benchmark_rows <- external_profile_rows(file.path(comparator, "benchmark_class_ids.tif"),
                                          fields, eligible, template, benchmark_lookup,
                                          "Benchmark", column, connectivity_final_figure_class_rows)
  profiles[[i]] <- rbind(matched_rows, connectivity_final_figure_summarise_profile_rows(
    rbind(native_rows, benchmark_rows)
  ))

  if (dataset == "BRAN") {
    columns <- connectivity_final_figure_aggregate_disagreement(contract)$columns
    columns <- columns[reporting_mask[cbind(columns$row, columns$col)], ]
    switched <- sum(columns$switched_seed_depth_pairs)
    observations <- sum(columns$eligible_seed_depth_pairs)
    if (observations <= 0 || observations != case$n_eligible_per_seed * case$n_seeds ||
        abs(100 * switched / observations - case$raw_switching_percent) > 1e-10) {
      stop("The disagreement map counts differ from the BRAN case ledger.")
    }
    write.csv(columns, file.path(out, "disagreement_by_cell.csv"), row.names = FALSE)
    difference <- ggplot(columns, aes(longitude, latitude)) +
      geom_tile(aes(fill = 100 * switch_proportion), width = fields$res[1], height = fields$res[2]) +
      scale_fill_viridis_c(name = "Switched (%)", limits = c(0, 100)) +
      map_extent + theme_bw(base_size = 11) +
      labs(title = "Matched 2D and 3D disagreement",
           subtitle = sprintf("%s of %s eligible depth and seed observations switched class",
                              format(switched, big.mark = ","), format(observations, big.mark = ",")),
           x = "Longitude (°E)", y = "Latitude", caption = caption)
    save_demo_png(difference, "04_disagreement", out, 7.8, 7.8)
  }
}
profiles <- do.call(rbind, profiles)
write.csv(profiles, file.path(out, "profiles_by_depth.csv"), row.names = FALSE)
profile <- report_profile_plot(profiles, class_colors_3d_fast(), benchmark_palette, caption)
save_demo_png(profile, "03_depth_profiles", out, 11.69, 11.02)

# 05. Each point removes the estimated seed noise from raw switching.
sensitivity$group <- ifelse(sensitivity$vertical_method == "temperature_layer_aggregate",
                            "Temperature-layer aggregation", "Centre sampling")
sensitivity$case_id <- factor(sensitivity$case_id, levels = rev(sensitivity$case_id))
write.csv(sensitivity, file.path(out, "sensitivity.csv"), row.names = FALSE)
stress <- ggplot(sensitivity, aes(attributable_pp, case_id)) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.45) +
  geom_errorbar(aes(xmin = interval_low_pp, xmax = interval_high_pp), orientation = "y", width = 0.3) +
  geom_point(aes(shape = precision_limited), size = 2.4) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(width = 19)) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +
  theme_bw(base_size = 11) +
  theme(strip.text.y = element_text(angle = 0),
        plot.caption = element_text(hjust = 0, size = 11, lineheight = 1.15, margin = margin(t = 9))) +
  labs(title = "Gascoyne Coast BRAN sensitivity",
       subtitle = "37 cases · five paired seeds per case",
       x = "Switching above 2D seed noise (percentage points)", y = NULL,
       caption = paste("Points show raw 2D–3D switching minus 2D seed noise.",
                       "The zero line marks equality with the seed noise estimate.",
                       "Bars show 95% within-case Monte Carlo intervals.",
                       "Triangles mark cases with 2D seed noise of 5% or more.", sep = "\n"))
save_demo_png(stress, "05_sensitivity", out, 10.5, 12)

write.csv(data.frame(figure = figure_names,
                     tables = c("context_inputs.csv", "surface_speed.csv, surface_velocity_arrows.csv",
                                "profiles_by_depth.csv", "disagreement_by_cell.csv", "sensitivity.csv")),
          file.path(out, "figure_sources.csv"), row.names = FALSE)
artifacts <- list.files(out, full.names = TRUE)
publish_csv(data.frame(run_id = run_id, file = basename(artifacts), bytes = file.info(artifacts)$size,
                       md5 = unname(tools::md5sum(artifacts))), file.path(out, "figure_manifest.csv"))
check_file_manifest(file.path(out, "figure_manifest.csv"))
figure_dir <- file.path(run_dir, "figures")
replace_output_folder(out, figure_dir, file.path(run_dir, "figure_history"))
publish_demo_outputs(figure_dir, run_dir, file.path(demo_dir, "outputs"))
message("Section complete: five PNGs and two result tables in ", file.path(demo_dir, "outputs"))
