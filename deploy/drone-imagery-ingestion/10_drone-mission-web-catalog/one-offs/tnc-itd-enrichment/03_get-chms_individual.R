# Purpose: Download CHM rasters for nadir missions only, using the nadir mission polygon GPKG (from
# script 02) to determine which missions are nadir, and the S3 file listing (from script 01) to
# determine which CHMs are available.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
TMP_DIR = file.path(DELIVERABLES_DIR, "tmp")
MISSIONS_S3_LISTING_FILEPATH = file.path(TMP_DIR, "s3-file-listing-missions.csv")
INDIVIDUAL_POLYGONS_FILEPATH = file.path(DELIVERABLES_DIR, "drone-plot-summaries/individual-missions/overall/individual-drone-plot-summaries.gpkg")
CHM_DIR = file.path(DELIVERABLES_DIR, "individual-missions/canopy-height-models")

dir.create(CHM_DIR, recursive = TRUE, showWarnings = FALSE)


# --- Helper functions ---

sanitize_url = function(url) {
  gsub("(?<!:)//", "/", url, perl = TRUE)
}


# --- Identify nadir missions with CHMs available ---

nadir_missions = st_read(INDIVIDUAL_POLYGONS_FILEPATH, quiet = TRUE)
nadir_mission_ids = nadir_missions$mission_id

missions_listing = read_csv(MISSIONS_S3_LISTING_FILEPATH, show_col_types = FALSE)

# CHM file path pattern: {mission_id}/photogrammetry_03/full/{mission_id}_chm-mesh.tif
chm_files = missions_listing |>
  filter(str_detect(filepath, "_chm-mesh\\.tif$")) |>
  filter(mission_id %in% nadir_mission_ids)

cat("Found", nrow(chm_files), "CHMs available for", length(nadir_mission_ids), "nadir missions\n")


# --- Download CHMs ---

download_chm = function(filepath, mission_id) {
  url = sanitize_url(paste0(DATA_SERVER_MISSIONS_BASE_URL, filepath))
  dest = file.path(CHM_DIR, paste0(mission_id, "_chm-mesh.tif"))

  if (file.exists(dest)) return(TRUE)

  download.file(url, dest, quiet = TRUE, method = "wget")
  TRUE
}

results = map2(
  chm_files$filepath,
  chm_files$mission_id,
  possibly(download_chm, otherwise = FALSE)
)

n_success = sum(unlist(results))
cat("Successfully downloaded", n_success, "of", nrow(chm_files), "CHMs to", CHM_DIR, "\n")
