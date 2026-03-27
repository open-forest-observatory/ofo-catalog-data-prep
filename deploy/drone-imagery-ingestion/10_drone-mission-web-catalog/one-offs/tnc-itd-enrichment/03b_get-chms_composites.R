# Purpose: Download CHM rasters for composite missions, using the composite mission polygon GPKG (from
# script 02) to determine which composites to include, and the S3 file listing (from script 01) to
# determine which CHMs are available.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
TMP_DIR = file.path(DELIVERABLES_DIR, "tmp")
COMPOSITES_S3_LISTING_FILEPATH = file.path(TMP_DIR, "s3-file-listing-composites.csv")
COMPOSITE_POLYGONS_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/composite-drone-plot-summaries.gpkg")
CHM_DIR = file.path(DELIVERABLES_DIR, "composite-missions/canopy-height-models")

dir.create(CHM_DIR, recursive = TRUE, showWarnings = FALSE)


# --- Helper functions ---

sanitize_url = function(url) {
  gsub("(?<!:)//", "/", url, perl = TRUE)
}


# --- Identify composite missions with CHMs available ---

composite_missions = st_read(COMPOSITE_POLYGONS_FILEPATH, quiet = TRUE)
composite_ids_in_gpkg = composite_missions$composite_id

composites_listing = read_csv(COMPOSITES_S3_LISTING_FILEPATH, show_col_types = FALSE)

# CHM file path pattern: {composite_id}/photogrammetry_01/full/{composite_id}_chm-mesh.tif
chm_files = composites_listing |>
  filter(str_detect(filepath, "_chm-mesh\\.tif$"),
         grepl("/full/", filepath),
         composite_id %in% composite_ids_in_gpkg)

cat("Found", nrow(chm_files), "CHMs available for", length(composite_ids_in_gpkg), "composite missions\n")


# --- Download CHMs ---

download_chm = function(filepath, composite_id) {
  url = sanitize_url(paste0(DATA_SERVER_COMPOSITES_BASE_URL, filepath))
  dest = file.path(CHM_DIR, paste0(composite_id, "_chm-mesh.tif"))

  if (file.exists(dest)) return(TRUE)

  download.file(url, dest, quiet = TRUE, method = "wget")
  TRUE
}

results = map2(
  chm_files$filepath,
  chm_files$composite_id,
  possibly(download_chm, otherwise = FALSE),
  .progress = TRUE
)

n_success = sum(unlist(results))
cat("Successfully downloaded", n_success, "of", nrow(chm_files), "CHMs to", CHM_DIR, "\n")
