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

LARGE_TREE_DBH_THRESHOLD = 30*2.54  # cm, representing 30 inch DBH
HEATMAP_RESOLUTION = 10  # meters
HEATMAP_WINDOW_SIZE = 50  # meters (1 ha square)

dir.create(HEATMAPS_PINE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(HEATMAPS_LARGE_TREES_DIR, recursive = TRUE, showWarnings = FALSE)


# --- Load data ---

all_trees = st_read(ALL_DETECTED_TREES_FILEPATH, quiet = TRUE)
footprints = st_read(FOOTPRINTS_FILEPATH, quiet = TRUE)


# --- Filter out dead trees ---

all_trees = all_trees |>
  filter(live_dead_prediction == "Live")



# --- Derived columns ---

# Compute basal area (sq m) from DBH (cm) for each tree
all_trees$basal_area_sqm = pi * (all_trees$dbh / 100 / 2)^2

# Identify pine trees (species code starting with 'PI')
all_trees$is_pine = str_detect(all_trees$species_prediction, "^PI")
all_trees$is_large = all_trees$dbh > LARGE_TREE_DBH_THRESHOLD


# --- Compute per-composite summary statistics ---

compute_composite_stats = function(trees, footprint) {
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

trees_split = split(all_trees, all_trees$composite_id)
footprints_split = split(footprints, footprints$composite_id)

# Align trees to footprints order, using empty sf for missing composites
trees_split = trees_split[names(footprints_split)]
empty_tree_sf = st_sf(composite_id = character(0), geometry = st_sfc(), crs = st_crs(footprints))
trees_split[sapply(trees_split, is.null)] = list(empty_tree_sf)

# Ensure they're aligned
all(names(trees_split) == names(footprints_split))

stats_list = map2(trees_split, footprints_split, compute_composite_stats, .progress = TRUE)

stats_df = bind_rows(stats_list)
footprints = bind_cols(footprints, stats_df)

st_write(footprints, FOOTPRINTS_WITH_STATS_FILEPATH, delete_dsn = TRUE)
cat("Saved footprints with summary statistics\n")


# --- Heatmap rasters ---

cat("Generating heatmap rasters...\n")

# Add pine_ba column for rasterization
all_trees$pine_ba = ifelse(all_trees$is_pine & !is.na(all_trees$is_pine), all_trees$basal_area_sqm, 0)

# Project trees and footprint to UTM for a given composite
prepare_utm = function(x) {
  centroid = st_coordinates(st_centroid(x))
  utm_zone = floor((centroid[1] + 180) / 6) + 1
  epsg = 32600 + utm_zone + ifelse(centroid[2] < 0, 100, 0)
  st_transform(x, epsg)
}

# Create a raster template at heatmap resolution from a UTM footprint
make_template = function(footprint) {
  bbox = st_bbox(footprint)
  rast(
    xmin = bbox["xmin"], xmax = bbox["xmax"],
    ymin = bbox["ymin"], ymax = bbox["ymax"],
    resolution = HEATMAP_RESOLUTION,
    crs = st_crs(footprint)$wkt
  )
}

# Pine proportion heatmap (by basal area): compute proportion at fine resolution,
# then smooth with a moving window mean
generate_pine_heatmap = function(trees, footprint, output_path) {
  trees_vect = vect(trees)
  template = make_template(footprint)

  pine_ba = rasterize(trees_vect, template, field = "pine_ba", fun = "sum")
  total_ba = rasterize(trees_vect, template, field = "basal_area_sqm", fun = "sum")
  proportion = pine_ba / total_ba

  # Smooth with a moving window mean
  window_cells = round(HEATMAP_WINDOW_SIZE / HEATMAP_RESOLUTION)
  result = focal(proportion, w = window_cells, fun = "mean", na.rm = TRUE)

  result = mask(result, vect(footprint))
  writeRaster(result, output_path, overwrite = TRUE)
}

# Large tree density heatmap: count large trees at fine resolution,
# then smooth with a moving window mean
generate_large_tree_heatmap = function(trees, footprint, output_path) {
  large_trees = trees |> filter(is_large)

  if (nrow(large_trees) == 0) {
    # If no large trees, create an empty raster with 0 values
    template = make_template(footprint)
    empty_rast = setValues(template, 0)
    writeRaster(empty_rast, output_path, overwrite = TRUE)
    return()
  }

  large_trees_vect = vect(large_trees)
  template = make_template(footprint)

  count_rast = rasterize(large_trees_vect, template, fun = "length")
  # Set cells within footprint but with no trees to 0 (not NA)
  footprint_rast = rasterize(vect(footprint), template)
  count_rast[!is.na(footprint_rast) & is.na(count_rast)] = 0

  # Smooth with a moving window mean
  window_cells = round(HEATMAP_WINDOW_SIZE / HEATMAP_RESOLUTION)
  result = focal(count_rast, w = window_cells, fun = "mean", na.rm = TRUE)

  # Convert from trees per grid cell to trees per hectare
  cell_area_ha = (HEATMAP_RESOLUTION^2) / 10000
  result = result / cell_area_ha

  result = mask(result, vect(footprint))
  writeRaster(result, output_path, overwrite = TRUE)
}

composite_ids = unique(footprints$composite_id)

for (i in seq_len(length(composite_ids))) {
  cid = composite_ids[i]

  cat("  Heatmaps for", cid, "(", i, "/", length(composite_ids), ")\n")

  footprint = footprints |> filter(composite_id == cid)
  trees = all_trees |> filter(composite_id == cid)

  if (nrow(trees) == 0) {
    cat("    No trees, skipping\n")
    next
  }

  if (st_is_longlat(footprint)) {
    footprint = prepare_utm(footprint)
  }

  trees = st_transform(trees, st_crs(footprint))

  pine_path = file.path(HEATMAPS_PINE_DIR, paste0(cid, "_pine-proportion.tif"))
  generate_pine_heatmap(trees, footprint, pine_path)

  large_path = file.path(HEATMAPS_LARGE_TREES_DIR, paste0(cid, "_large-tree-density.tif"))
  generate_large_tree_heatmap(trees, footprint, large_path)

}

cat("Done generating heatmaps\n")
