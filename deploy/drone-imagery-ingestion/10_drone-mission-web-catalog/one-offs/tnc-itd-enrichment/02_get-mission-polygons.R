# Purpose: Download composite and individual mission polygons from S3. For composites, intersect the
# 2 features (one per constituent mission) into a single polygon. For individual missions, filter to
# nadir only. Save each as a separate GPKG.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
TMP_DIR = file.path(DELIVERABLES_DIR, "tmp")
MISSIONS_S3_LISTING_FILEPATH = file.path(TMP_DIR, "s3-file-listing-missions.csv")
COMPOSITES_S3_LISTING_FILEPATH = file.path(TMP_DIR, "s3-file-listing-composites.csv")

COMPOSITE_POLYGONS_DIR = file.path(DELIVERABLES_DIR, "composite-missions/overall")
INDIVIDUAL_POLYGONS_DIR = file.path(DELIVERABLES_DIR, "individual-missions/overall")
COMPOSITE_POLYGONS_FILEPATH = file.path(COMPOSITE_POLYGONS_DIR, "composite-drone-plot-summaries.gpkg")
INDIVIDUAL_POLYGONS_FILEPATH = file.path(INDIVIDUAL_POLYGONS_DIR, "individual-drone-plot-summaries.gpkg")

dir.create(COMPOSITE_POLYGONS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INDIVIDUAL_POLYGONS_DIR, recursive = TRUE, showWarnings = FALSE)


# --- Helper functions ---

# Replace double slashes with single slashes in a URL, preserving the :// in the protocol
sanitize_url = function(url) {
  gsub("(?<!:)//", "/", url, perl = TRUE)
}

sf_from_url = function(url) {
  temp_file = tempfile(fileext = ".gpkg")
  download.file(url, temp_file, quiet = TRUE, method = "wget")
  sf = st_read(temp_file, quiet = TRUE)
  unlink(temp_file)
  return(sf)
}

# Force all non-geometry columns to character (for consistent bind_rows across missions with
# differing column types)
force_all_cols_to_character = function(df) {
  df |>
    mutate(across(-any_of(c("geometry", "geom")), as.character))
}


# --- Composite mission polygons ---

cat("Processing composite mission polygons...\n")

composites_listing = read_csv(COMPOSITES_S3_LISTING_FILEPATH, show_col_types = FALSE)

composite_polygon_files = composites_listing |>
  filter(str_detect(filepath, "metadata-mission/.*_mission-metadata\\.gpkg$")) |>
  filter(composite_id != "000000_000000")

cat("Found", nrow(composite_polygon_files), "composite polygon files\n")

# Download a composite polygon GPKG and intersect its 2 features into a single polygon
download_and_intersect_composite = function(composite_id, filepath) {
  url = sanitize_url(paste0(DATA_SERVER_COMPOSITES_BASE_URL, filepath))
  sf = sf_from_url(url)

  if (nrow(sf) > 1) {
    intersected = reduce(st_geometry(sf), st_intersection)
  } else {
    intersected = st_geometry(sf)[[1]]
  }

  st_sf(
    composite_id = composite_id,
    geometry = st_sfc(intersected, crs = st_crs(sf))
  )
}

composite_polygons = map2(
  composite_polygon_files$composite_id,
  composite_polygon_files$filepath,
  possibly(download_and_intersect_composite, otherwise = NULL)
)

composite_polygons = compact(composite_polygons)
composite_polygons = bind_rows(composite_polygons)

st_write(composite_polygons, COMPOSITE_POLYGONS_FILEPATH, delete_dsn = TRUE)
cat("Wrote", nrow(composite_polygons), "composite polygons to", COMPOSITE_POLYGONS_FILEPATH, "\n")


# --- Individual mission polygons (nadir only) ---

cat("\nProcessing individual mission polygons...\n")

missions_listing = read_csv(MISSIONS_S3_LISTING_FILEPATH, show_col_types = FALSE)

mission_polygon_files = missions_listing |>
  filter(str_detect(filepath, "metadata-mission/.*_mission-metadata\\.gpkg$"))

cat("Found", nrow(mission_polygon_files), "individual mission polygon files\n")

download_mission_polygon = function(filepath) {
  url = sanitize_url(paste0(DATA_SERVER_MISSIONS_BASE_URL, filepath))
  sf_from_url(url)
}

mission_polygons = map(
  mission_polygon_files$filepath,
  possibly(download_mission_polygon, otherwise = NULL)
)

mission_polygons = compact(mission_polygons)
mission_polygons = map(mission_polygons, force_all_cols_to_character)
mission_polygons = bind_rows(mission_polygons)

cat("Downloaded", nrow(mission_polygons), "individual mission polygons\n")

# Filter to nadir missions: at least one pitch attribute < 10, AND altitude > 90
mission_polygons = mission_polygons |>
  mutate(
    camera_pitch_nominal = as.numeric(camera_pitch_nominal),
    camera_pitch_derived = as.numeric(camera_pitch_derived),
    altitude_agl_nominal = as.numeric(altitude_agl_nominal)
  ) |>
  filter(
    (!is.na(camera_pitch_nominal) & camera_pitch_nominal < 10 |
     !is.na(camera_pitch_derived) & camera_pitch_derived < 10) &
    !is.na(altitude_agl_nominal) & altitude_agl_nominal > 90
  )

cat("Filtered to", nrow(mission_polygons), "nadir missions\n")

st_write(mission_polygons, INDIVIDUAL_POLYGONS_FILEPATH, delete_dsn = TRUE)
cat("Wrote nadir mission polygons to", INDIVIDUAL_POLYGONS_FILEPATH, "\n")
