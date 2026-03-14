# Purpose: For each nadir mission with a CHM, compute overall canopy cover and mean canopy height
# (canopy defined as > 5 m), add as attributes to the nadir mission GPKG, and produce a moving-window
# canopy cover raster (100 m square window) for each mission.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)
library(terra)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
INDIVIDUAL_POLYGONS_FILEPATH = file.path(DELIVERABLES_DIR, "individual-missions/overall/individual-drone-plot-summaries.gpkg")
CHM_DIR = file.path(DELIVERABLES_DIR, "individual-missions/canopy-height-models")
CANOPY_COVER_RASTER_DIR = file.path(DELIVERABLES_DIR, "canopy-cover-rasters")

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

  # Determine window size in pixels for a 100 m square
  res_m = res(chm)[1]
  window_size = round(100 / res_m)
  # Ensure odd window size for symmetric focal window
  if (window_size %% 2 == 0) window_size = window_size + 1

  # Focal mean of binary raster = proportion canopy cover
  canopy_cover_raster = focal(canopy_binary, w = window_size, fun = "mean", na.rm = TRUE)

  writeRaster(canopy_cover_raster, output_path, overwrite = TRUE)
}


# --- Process all nadir missions ---

nadir_missions = st_read(INDIVIDUAL_POLYGONS_FILEPATH, quiet = TRUE)

chm_files = list.files(CHM_DIR, pattern = "_chm-mesh\\.tif$", full.names = TRUE)
chm_mission_ids = str_extract(basename(chm_files), "^[^_]+")

cat("Found", length(chm_files), "CHMs for", nrow(nadir_missions), "nadir missions\n")

# Initialize new columns
nadir_missions$canopy_cover = NA_real_
nadir_missions$canopy_height = NA_real_

for (i in seq_along(chm_files)) {
  mission_id = chm_mission_ids[i]
  chm_path = chm_files[i]
  cat("  Processing", mission_id, "(", i, "/", length(chm_files), ")\n")

  # Compute summary stats
  stats = tryCatch(
    compute_canopy_stats(chm_path),
    error = function(e) {
      cat("    WARNING: Failed to compute stats:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(stats)) {
    row_idx = which(nadir_missions$mission_id == mission_id)
    if (length(row_idx) == 1) {
      nadir_missions$canopy_cover[row_idx] = stats$canopy_cover
      nadir_missions$canopy_height[row_idx] = stats$canopy_height
    }
  }

  # Compute moving-window canopy cover raster
  output_path = file.path(CANOPY_COVER_RASTER_DIR, paste0(mission_id, "_canopy-cover.tif"))
  tryCatch(
    compute_canopy_cover_raster(chm_path, output_path),
    error = function(e) {
      cat("    WARNING: Failed to compute raster:", conditionMessage(e), "\n")
    }
  )
}

# Save updated GPKG with canopy stats
st_write(nadir_missions, INDIVIDUAL_POLYGONS_FILEPATH, delete_dsn = TRUE)
cat("Updated nadir mission GPKG with canopy_cover and canopy_height attributes\n")
