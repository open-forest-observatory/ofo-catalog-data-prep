# Purpose: Download all individual nadir mission ITD (detected tree) data from S3 and compile into a
# single GPKG.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
TMP_DIR = file.path(DELIVERABLES_DIR, "tmp")
MISSIONS_S3_LISTING_FILEPATH = file.path(TMP_DIR, "s3-file-listing-missions.csv")
INDIVIDUAL_POLYGONS_FILEPATH = file.path(DELIVERABLES_DIR,
                                        "individual-missions/overall/individual-drone-plot-summaries.gpkg")
ALL_DETECTED_TREES_FILEPATH = file.path(DELIVERABLES_DIR,
                                        "individual-missions/overall/all-detected-trees.gpkg")

# Overriding the default constants for bandaid fix to be able to specify two different
# photogrammetry config IDs (composites vs individuals)
INDIVIDUAL_PHOTOGRAMMETRY_CONFIG_ID = "03"


dir.create(dirname(ALL_DETECTED_TREES_FILEPATH), recursive = TRUE, showWarnings = FALSE)


# --- Helper functions ---

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

force_all_cols_to_character = function(df) {
  df |>
    mutate(across(-any_of(c("geometry", "geom")), as.character))
}

# From a filtered file listing, find the most recent ITD folder (same logic as web catalog code)
get_most_recent_itd_folder = function(filepaths) {
  filepath_parts = str_split(filepaths, fixed("/"))
  # Not all filepaths have 3+ parts; safely extract part 3 where it exists
  part_3 = map_chr(filepath_parts, \(x) if (length(x) >= 3) x[3] else NA_character_)
  itd_folders = part_3[str_which(part_3, "^(itd_|detected-trees_)")] |> unique() |> sort(decreasing = TRUE)
  if (length(itd_folders) == 0) return(NA_character_)
  itd_folders[1]
}


# --- Find ITD treetop files for each mission ---

missions_listing = read_csv(MISSIONS_S3_LISTING_FILEPATH, show_col_types = FALSE)

nadir_missions = st_read(INDIVIDUAL_POLYGONS_FILEPATH, quiet = TRUE)
nadir_mission_ids = nadir_missions$mission_id

mission_ids = missions_listing |>
  pull(mission_id) |>
  unique()

mission_ids = intersect(mission_ids, nadir_mission_ids)

# For each mission, find the most recent treetops file
find_ttops_filepath = function(mission_id) {
  mission_files = missions_listing |>
    filter(mission_id == !!mission_id)

  itd_folder = get_most_recent_itd_folder(mission_files$filepath)
  if (is.na(itd_folder)) return(NULL)

  # Prefer v2 (classified) over v1
  ttops_v2 = file.path(mission_id, paste0("photogrammetry_", INDIVIDUAL_PHOTOGRAMMETRY_CONFIG_ID), itd_folder, paste0(mission_id, "_detected-tree-tops_classified.gpkg"))
  ttops_v1 = file.path(mission_id, paste0("photogrammetry_", INDIVIDUAL_PHOTOGRAMMETRY_CONFIG_ID), itd_folder, paste0(mission_id, "_treetops.gpkg"))

  if (ttops_v2 %in% mission_files$filepath) return(ttops_v2)
  if (ttops_v1 %in% mission_files$filepath) return(ttops_v1)
  return(NULL)
}

ttops_filepaths = map(mission_ids, find_ttops_filepath)
names(ttops_filepaths) = mission_ids
ttops_filepaths = compact(ttops_filepaths)

cat("Found ITD data for", length(ttops_filepaths), "of", length(mission_ids), "nadir missions\n")


# --- Download and compile ---

download_ttops = function(mission_id, filepath) {
  url = sanitize_url(paste0(DATA_SERVER_MISSIONS_BASE_URL, filepath))
  sf = sf_from_url(url)
  sf$mission_id = mission_id
  sf
}

all_trees = map2(
  names(ttops_filepaths),
  unlist(ttops_filepaths),
  possibly(download_ttops, otherwise = NULL),
  .progress = TRUE
)

all_trees = compact(all_trees)
all_trees = map(all_trees, force_all_cols_to_character)
all_trees = map(all_trees, st_transform, crs = 4326)
all_trees = bind_rows(all_trees)

cat("Downloaded", nrow(all_trees), "detected trees from", length(ttops_filepaths), "missions\n")

st_write(all_trees, ALL_DETECTED_TREES_FILEPATH, delete_dsn = TRUE)
cat("Wrote all detected trees to", ALL_DETECTED_TREES_FILEPATH, "\n")
