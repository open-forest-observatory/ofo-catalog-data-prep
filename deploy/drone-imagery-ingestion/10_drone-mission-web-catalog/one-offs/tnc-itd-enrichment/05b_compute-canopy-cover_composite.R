# Purpose: For each composite mission with a CHM, compute overall canopy cover and mean canopy height
# (canopy defined as > 5 m), add as attributes to the composite mission GPKG, and produce a moving-window
# canopy cover raster (100 m diameter circular window) for each composite.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)
library(terra)
library(dplyr)
library(furrr)

plan(multisession, workers = availableCores()*2)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
COMPOSITE_POLYGONS_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/composite-drone-plot-summaries.gpkg")
CHM_DIR = file.path(DELIVERABLES_DIR, "composite-missions/canopy-height-models")
CANOPY_COVER_RASTER_DIR = file.path(DELIVERABLES_DIR, "composite-missions/canopy-cover-rasters")

CANOPY_HEIGHT_THRESHOLD = 5  # meters

dir.create(CANOPY_COVER_RASTER_DIR, recursive = TRUE, showWarnings = FALSE)


# --- Compute canopy stats for a single mission ---

compute_canopy_stats = function(chm_path) {
  chm = rast(chm_path)
  vals = values(chm, na.rm = TRUE)

  canopy_vals = vals[vals > CANOPY_HEIGHT_THRESHOLD]
  canopy_cover = length(canopy_vals) / length(vals)
  canopy_height = if (length(canopy_vals) > 0) mean(canopy_vals) else NA_real_

  list(canopy_cover = canopy_cover, canopy_height = canopy_height)
}


# --- Compute moving-window canopy cover raster ---

compute_canopy_cover_raster = function(chm_path, output_path) {
  chm = rast(chm_path)

  # Create binary canopy raster (1 = canopy, 0 = not canopy)
  canopy_binary = chm > CANOPY_HEIGHT_THRESHOLD

  # Resample to 1 m resolution for more efficient processing
  res_change_fact = 1 / res(chm)[1]
  if (res_change_fact > 1) {
    canopy_binary = aggregate(canopy_binary, fact = round(res_change_fact), fun = "mean", na.rm = TRUE)
  }

  # Circular moving window with 50 m radius (~100 m diameter)
  w = focalMat(canopy_binary, d = 50, type = "circle")
  canopy_cover_raster = focal(canopy_binary, w = w, fun = "mean", na.rm = TRUE, na.policy = "omit")

  writeRaster(canopy_cover_raster, output_path, overwrite = TRUE)
}


# --- Process all composite missions ---

composite_missions = st_read(COMPOSITE_POLYGONS_FILEPATH, quiet = TRUE)

chm_files = list.files(CHM_DIR, pattern = "_chm-mesh\\.tif$", full.names = TRUE)
chm_composite_ids = str_extract(basename(chm_files), "^[^_]+_[^_]+")

cat("Found", length(chm_files), "CHMs for", nrow(composite_missions), "composite missions\n")

process_mission = function(i) {
  library(terra)
  composite_id = chm_composite_ids[i]
  chm_path = chm_files[i]

  # Compute summary stats
  stats = tryCatch(
    compute_canopy_stats(chm_path),
    error = function(e) {
      warning("Failed to compute stats for ", composite_id, ": ", conditionMessage(e))
      NULL
    }
  )

  # Compute moving-window canopy cover raster
  output_path = file.path(CANOPY_COVER_RASTER_DIR, paste0(composite_id, "_canopy-cover.tif"))
  tryCatch(
    compute_canopy_cover_raster(chm_path, output_path),
    error = function(e) {
      warning("Failed to compute raster for ", composite_id, ": ", conditionMessage(e))
    }
  )

  if (is.null(stats)) {
    tibble(composite_id = composite_id, canopy_cover = NA_real_, canopy_height = NA_real_)
  } else {
    tibble(composite_id = composite_id, canopy_cover = stats$canopy_cover, canopy_height = stats$canopy_height)
  }
}

results = future_map(seq_along(chm_files), process_mission,
                     .progress = TRUE,
                     .options = furrr_options(scheduling = Inf)) |>
  bind_rows()

# Join results into composite missions
composite_missions = composite_missions |>
  left_join(results, by = "composite_id")

# Save updated GPKG with canopy stats
st_write(composite_missions, COMPOSITE_POLYGONS_FILEPATH, delete_dsn = TRUE)
cat("Updated composite mission GPKG with canopy_cover and canopy_height attributes\n")
