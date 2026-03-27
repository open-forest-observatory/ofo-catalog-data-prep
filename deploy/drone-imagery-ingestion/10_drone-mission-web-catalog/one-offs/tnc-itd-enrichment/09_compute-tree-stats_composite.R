# Purpose: Compute per-composite summary statistics (tree density, basal area, large tree density,
# proportion pine), add as attributes to the composite drone plot summaries GPKG, and produce
# heatmap rasters for pine proportion and large tree density.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)
library(terra)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
ALL_DETECTED_TREES_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-trees_analysis-ready.gpkg")
FOOTPRINTS_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/forested-mission-footprints_analysis-ready.gpkg")
FOOTPRINTS_WITH_STATS_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/forested-mission-footprints_analysis-ready-with-stats.gpkg")

HEATMAPS_PINE_DIR = file.path(DELIVERABLES_DIR, "composite-missions/heatmaps-pine-proportion")
HEATMAPS_LARGE_TREES_DIR = file.path(DELIVERABLES_DIR, "composite-missions/heatmaps-large-trees")

LARGE_TREE_DBH_THRESHOLD = 30  # cm
HEATMAP_RESOLUTION = 10  # meters
HEATMAP_WINDOW_SIZE = 100  # meters (1 ha square)

dir.create(HEATMAPS_PINE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(HEATMAPS_LARGE_TREES_DIR, recursive = TRUE, showWarnings = FALSE)


# --- Load data ---

all_trees = st_read(ALL_DETECTED_TREES_FILEPATH, quiet = TRUE)
footprints = st_read(FOOTPRINTS_FILEPATH, quiet = TRUE)


# --- Derived columns ---

# Compute basal area (sq m) from DBH (cm) for each tree
all_trees$basal_area_sqm = pi * (all_trees$dbh / 100 / 2)^2

# Identify pine trees (species code starting with 'PI')
all_trees$is_pine = str_detect(all_trees$species_prediction, "^PI")
all_trees$is_large = all_trees$dbh > LARGE_TREE_DBH_THRESHOLD


# --- Compute per-composite summary statistics ---

compute_composite_stats = function(composite_id, footprint) {
  trees = all_trees |> filter(composite_id == !!composite_id)

  if (nrow(trees) == 0) {
    return(tibble(
      tree_density_per_ha = NA_real_, tree_density_per_acre = NA_real_,
      basal_area_sqm_per_ha = NA_real_, basal_area_sqft_per_acre = NA_real_,
      large_tree_density_per_ha = NA_real_, large_tree_density_per_acre = NA_real_,
      proportion_pine = NA_real_
    ))
  }

  # Compute forested area from footprint
  area_sqm = as.numeric(st_area(footprint))
  area_ha = area_sqm / 10000
  area_acre = area_ha * 2.47105

  # Tree density
  tree_density_per_ha = nrow(trees) / area_ha
  tree_density_per_acre = nrow(trees) / area_acre

  # Basal area
  total_basal_area_sqm = sum(trees$basal_area_sqm, na.rm = TRUE)
  basal_area_sqm_per_ha = total_basal_area_sqm / area_ha
  basal_area_sqft_per_acre = (total_basal_area_sqm * 10.7639) / area_acre

  # Large tree density
  n_large = sum(trees$is_large, na.rm = TRUE)
  large_tree_density_per_ha = n_large / area_ha
  large_tree_density_per_acre = n_large / area_acre

  # Proportion pine (by basal area)
  pine_basal_area = sum(trees$basal_area_sqm[trees$is_pine & !is.na(trees$is_pine)], na.rm = TRUE)
  proportion_pine = pine_basal_area / total_basal_area_sqm

  tibble(
    tree_density_per_ha, tree_density_per_acre,
    basal_area_sqm_per_ha, basal_area_sqft_per_acre,
    large_tree_density_per_ha, large_tree_density_per_acre,
    proportion_pine
  )
}

cat("Computing per-composite summary statistics...\n")

stats_list = map(seq_len(nrow(footprints)), function(i) {
  cat("  Processing", footprints$composite_id[i], "(", i, "/", nrow(footprints), ")\n")
  tryCatch(
    compute_composite_stats(footprints$composite_id[i], footprints[i, ]),
    error = function(e) {
      cat("    WARNING: Failed:", conditionMessage(e), "\n")
      tibble(
        tree_density_per_ha = NA_real_, tree_density_per_acre = NA_real_,
        basal_area_sqm_per_ha = NA_real_, basal_area_sqft_per_acre = NA_real_,
        large_tree_density_per_ha = NA_real_, large_tree_density_per_acre = NA_real_,
        proportion_pine = NA_real_
      )
    }
  )
})

stats_df = bind_rows(stats_list)
footprints = bind_cols(footprints, stats_df)

st_write(footprints, FOOTPRINTS_WITH_STATS_FILEPATH, delete_dsn = TRUE)
cat("Saved footprints with summary statistics\n")


# --- Heatmap rasters ---

cat("Generating heatmap rasters...\n")

# Add pine_ba column for rasterization
all_trees$pine_ba = ifelse(all_trees$is_pine & !is.na(all_trees$is_pine), all_trees$basal_area_sqm, 0)

# Project trees and footprint to UTM for a given composite
prepare_utm = function(composite_id, footprint) {
  trees = all_trees |> filter(composite_id == !!composite_id)
  if (nrow(trees) == 0) return(NULL)

  utm_zone = floor((st_coordinates(st_centroid(footprint))[1] + 180) / 6) + 1
  utm_crs = st_crs(paste0("+proj=utm +zone=", utm_zone, " +datum=WGS84"))

  list(
    trees = st_transform(trees, utm_crs),
    footprint = st_transform(footprint, utm_crs),
    crs = utm_crs
  )
}

# Create a coarse (1 ha = 100 m) raster template from a UTM footprint
make_coarse_template = function(footprint_utm, utm_crs) {
  bbox = st_bbox(footprint_utm)
  rast(
    xmin = bbox["xmin"], xmax = bbox["xmax"],
    ymin = bbox["ymin"], ymax = bbox["ymax"],
    resolution = HEATMAP_WINDOW_SIZE,
    crs = st_crs(utm_crs)$wkt
  )
}

# Pine proportion heatmap (by basal area): sum pine BA / sum total BA per 1 ha cell,
# then disaggregate to 10 m with bilinear interpolation
generate_pine_heatmap = function(utm_data, output_path) {
  trees_vect = vect(utm_data$trees)
  template = make_coarse_template(utm_data$footprint, utm_data$crs)

  pine_ba = rasterize(trees_vect, template, field = "pine_ba", fun = "sum")
  total_ba = rasterize(trees_vect, template, field = "basal_area_sqm", fun = "sum")
  proportion = pine_ba / total_ba

  # Disaggregate from 100 m to 10 m with bilinear smoothing
  factor = HEATMAP_WINDOW_SIZE / HEATMAP_RESOLUTION
  result = disagg(proportion, fact = factor, method = "bilinear")

  result = mask(result, vect(utm_data$footprint))
  writeRaster(result, output_path, overwrite = TRUE)
}

# Large tree density heatmap: count large trees per 1 ha cell (= trees/ha),
# then disaggregate to 10 m with bilinear interpolation
generate_large_tree_heatmap = function(utm_data, output_path) {
  large_trees = utm_data$trees |> filter(is_large)
  if (nrow(large_trees) == 0) {
    # Write an empty raster
    template = make_coarse_template(utm_data$footprint, utm_data$crs)
    factor = HEATMAP_WINDOW_SIZE / HEATMAP_RESOLUTION
    result = disagg(template, fact = factor)
    result = mask(result, vect(utm_data$footprint))
    writeRaster(result, output_path, overwrite = TRUE)
    return()
  }

  template = make_coarse_template(utm_data$footprint, utm_data$crs)
  large_trees_vect = vect(large_trees)
  count_rast = rasterize(large_trees_vect, template, fun = "length")

  # Each cell is 1 ha (100m x 100m), so count = density per ha
  factor = HEATMAP_WINDOW_SIZE / HEATMAP_RESOLUTION
  result = disagg(count_rast, fact = factor, method = "bilinear")

  result = mask(result, vect(utm_data$footprint))
  writeRaster(result, output_path, overwrite = TRUE)
}

for (i in seq_len(nrow(footprints))) {
  cid = footprints$composite_id[i]
  cat("  Heatmaps for", cid, "(", i, "/", nrow(footprints), ")\n")

  tryCatch({
    utm_data = prepare_utm(cid, footprints[i, ])
    if (is.null(utm_data)) {
      cat("    No trees, skipping\n")
      next
    }

    pine_path = file.path(HEATMAPS_PINE_DIR, paste0(cid, "_pine-proportion.tif"))
    generate_pine_heatmap(utm_data, pine_path)

    large_path = file.path(HEATMAPS_LARGE_TREES_DIR, paste0(cid, "_large-tree-density.tif"))
    generate_large_tree_heatmap(utm_data, large_path)
  }, error = function(e) {
    cat("    WARNING: Failed:", conditionMessage(e), "\n")
  })
}

cat("Done generating heatmaps\n")
