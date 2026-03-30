# Purpose: Filter drone-detected trees and footrpints (both from composite) to Yuba focal area
# (buffered by 10 km). Filter ground plot footprints and trees to same focal area.

library(tidyverse)
library(sf)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"

# Input paths
YUBA_BOUNDARY_FILEPATH    = file.path(DELIVERABLES_DIR, "raw-input/north_yuba_area.kml")
GROUND_TREES_FILEPATH     = file.path(DELIVERABLES_DIR, "ground-reference/ofo_ground-reference_trees_foc.gpkg")
GROUND_PLOTS_FILEPATH     = file.path(DELIVERABLES_DIR, "ground-reference/ofo_ground-reference_plots.gpkg")
FOOTPRINTS_FILEPATH           = file.path(DELIVERABLES_DIR, "composite-missions/overall/forested-mission-footprints_analysis-ready-with-stats.gpkg")
PLOT_SUMMARIES_FILEPATH       = file.path(DELIVERABLES_DIR, "composite-missions/overall/composite-drone-plot-summaries.gpkg")
NADIR_PLOT_SUMMARIES_FILEPATH = file.path(DELIVERABLES_DIR, "individual-missions/overall/individual-drone-plot-summaries.gpkg")
DETECTED_TREES_FILEPATH   = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-trees_analysis-ready.gpkg")
DETECTED_CROWNS_FILEPATH  = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-crowns_analysis-ready.gpkg")

# Output paths
GROUND_OUTPUT_DIR         = file.path(DELIVERABLES_DIR, "ground-reference-delivery")
FOOTPRINTS_OUTPUT_DIR     = file.path(DELIVERABLES_DIR, "composite-missions/footprint-boundaries-delivery")
TREES_OUTPUT_DIR          = file.path(DELIVERABLES_DIR, "composite-missions/detected-trees-delivery")

YUBA_BUFFER_M = 10000  # 10 km buffer around Yuba boundary

dir.create(GROUND_OUTPUT_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(FOOTPRINTS_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TREES_OUTPUT_DIR,     recursive = TRUE, showWarnings = FALSE)


# --- Load Yuba boundary and create buffered focal area ---

yuba_boundary = st_read(YUBA_BOUNDARY_FILEPATH, quiet = TRUE)

# Project to UTM zone 10N for metric buffering, then buffer, then back to WGS84
yuba_utm = st_transform(yuba_boundary, 32610)
yuba_buffered_utm = st_buffer(st_union(yuba_utm), YUBA_BUFFER_M)
yuba_buffered = st_transform(yuba_buffered_utm, 4326)

cat("Yuba focal area (buffered by", YUBA_BUFFER_M / 1000, "km) computed\n")


# --- Load and filter ground data ---

# Note: ofo_ground-reference_trees_foc.gpkg was written by 21_prep-itd-stats.R and already
# incorporates the following typo/erroneous-value fixes from the raw trees file:
#   - Plot 0108: DBH 894.08 corrected to 89.408 (likely a decimal-place typo)
#   - Plot 0042: tree with DBH 82.296 removed (erroneous — implausibly large for that height)
# It also already filters to height > 10 m and live_dead == "L".
# We repeat the fixes below as a safeguard in case the upstream file is ever regenerated without them.

ground_trees = st_read(GROUND_TREES_FILEPATH, quiet = TRUE)
ground_plots = st_read(GROUND_PLOTS_FILEPATH, quiet = TRUE)

# Propagate typo/erroneous-value fixes (no-ops if already corrected in the source file)
ground_trees[ground_trees$plot_id == "0108" & !is.na(ground_trees$dbh) & ground_trees$dbh == 894.08, "dbh"] = 89.408
ground_trees = ground_trees |> filter(!(plot_id == "0042" & !is.na(dbh) & dbh == 82.296))

# Filter to Yuba focal area
ground_trees_yuba = ground_trees |>
  st_transform(4326) |>
  st_filter(yuba_buffered)

ground_plots_yuba = ground_plots |>
  st_transform(4326) |>
  st_filter(yuba_buffered)

cat("Ground trees: ", nrow(ground_trees), "->", nrow(ground_trees_yuba), "within Yuba focal area\n")
cat("Ground plots: ", nrow(ground_plots), "->", nrow(ground_plots_yuba), "within Yuba focal area\n")

st_write(ground_trees_yuba, file.path(GROUND_OUTPUT_DIR, "ofo_ground-reference_trees_yuba.gpkg"), delete_dsn = TRUE)
st_write(ground_plots_yuba, file.path(GROUND_OUTPUT_DIR, "ofo_ground-reference_plots_yuba.gpkg"), delete_dsn = TRUE)

cat("Saved ground data to", GROUND_OUTPUT_DIR, "\n")


# --- Load and filter drone footprints ---

footprints        = st_read(FOOTPRINTS_FILEPATH, quiet = TRUE)
plot_summaries    = st_read(PLOT_SUMMARIES_FILEPATH, quiet = TRUE)
nadir_plot_summaries = st_read(NADIR_PLOT_SUMMARIES_FILEPATH, quiet = TRUE)

footprints_yuba = footprints |>
  st_transform(4326) |>
  st_filter(yuba_buffered)

plot_summaries_yuba = plot_summaries |>
  st_transform(4326) |>
  st_filter(yuba_buffered)

nadir_plot_summaries_yuba = nadir_plot_summaries |>
  st_transform(4326) |>
  st_filter(yuba_buffered)

cat("Footprints: ", nrow(footprints), "->", nrow(footprints_yuba), "within Yuba focal area\n")
cat("Composite plot summaries: ", nrow(plot_summaries), "->", nrow(plot_summaries_yuba), "within Yuba focal area\n")
cat("Nadir plot summaries: ", nrow(nadir_plot_summaries), "->", nrow(nadir_plot_summaries_yuba), "within Yuba focal area\n")

st_write(footprints_yuba,          file.path(FOOTPRINTS_OUTPUT_DIR, "composite-drone-detected-tree-stats_yuba.gpkg"), delete_dsn = TRUE)
st_write(plot_summaries_yuba,      file.path(FOOTPRINTS_OUTPUT_DIR, "composite-drone-plot-summaries_yuba.gpkg"), delete_dsn = TRUE)
st_write(nadir_plot_summaries_yuba, file.path(FOOTPRINTS_OUTPUT_DIR, "individual-drone-plot-summaries_yuba.gpkg"), delete_dsn = TRUE)

cat("Saved footprint data to", FOOTPRINTS_OUTPUT_DIR, "\n")


# --- Load and filter detected trees and crowns ---

detected_trees = st_read(DETECTED_TREES_FILEPATH, quiet = TRUE)
detected_crowns = st_read(DETECTED_CROWNS_FILEPATH, quiet = TRUE)

# Keep only composites whose footprint intersects the Yuba focal area
yuba_composite_ids = footprints_yuba$composite_id

detected_trees_yuba = detected_trees |>
  filter(composite_id %in% yuba_composite_ids)

detected_crowns_yuba = detected_crowns |>
  filter(composite_id %in% yuba_composite_ids)

cat("Detected trees: ", nrow(detected_trees), "->", nrow(detected_trees_yuba), "within Yuba composite footprints\n")
cat("Detected crowns:", nrow(detected_crowns), "->", nrow(detected_crowns_yuba), "within Yuba composite footprints\n")

st_write(detected_trees_yuba,  file.path(TREES_OUTPUT_DIR, "all-detected-trees_analysis-ready_yuba.gpkg"),  delete_dsn = TRUE)
st_write(detected_crowns_yuba, file.path(TREES_OUTPUT_DIR, "all-detected-crowns_analysis-ready_yuba.gpkg"), delete_dsn = TRUE)

cat("Saved detected tree/crown data to", TREES_OUTPUT_DIR, "\n")

cat("\nDone.\n")
