# Purpose: Download all composite ITD (detected tree) data from S3 and compile into a single GPKG.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
TMP_DIR = file.path(DELIVERABLES_DIR, "tmp")
COMPOSITES_S3_LISTING_FILEPATH = file.path(TMP_DIR, "s3-file-listing-composites.csv")
ALL_DETECTED_TREES_FILEPATH = file.path(DELIVERABLES_DIR, "drone-plot-summaries/composite-missions/overall/all-detected-trees.gpkg")

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


# --- Find ITD treetop files for each composite ---

composites_listing = read_csv(COMPOSITES_S3_LISTING_FILEPATH, show_col_types = FALSE)

composite_ids = composites_listing |>
  filter(composite_id != "000000_000000") |>
  pull(composite_id) |>
  unique()

# For each composite, find the most recent treetops file
find_ttops_filepath = function(composite_id) {
  composite_files = composites_listing |>
    filter(composite_id == !!composite_id)

  itd_folder = get_most_recent_itd_folder(composite_files$filepath)
  if (is.na(itd_folder)) return(NULL)

  # Prefer v2 (classified) over v1
  ttops_v2 = file.path(composite_id, paste0("photogrammetry_", PHOTOGRAMMETRY_CONFIG_ID), itd_folder, paste0(composite_id, "_detected-tree-tops_classified.gpkg"))
  ttops_v1 = file.path(composite_id, paste0("photogrammetry_", PHOTOGRAMMETRY_CONFIG_ID), itd_folder, paste0(composite_id, "_treetops.gpkg"))

  if (ttops_v2 %in% composite_files$filepath) return(ttops_v2)
  if (ttops_v1 %in% composite_files$filepath) return(ttops_v1)
  return(NULL)
}

ttops_filepaths = map(composite_ids, find_ttops_filepath)
names(ttops_filepaths) = composite_ids
ttops_filepaths = compact(ttops_filepaths)

cat("Found ITD data for", length(ttops_filepaths), "of", length(composite_ids), "composites\n")


# --- Download and compile ---

download_ttops = function(composite_id, filepath) {
  url = sanitize_url(paste0(DATA_SERVER_COMPOSITES_BASE_URL, filepath))
  sf = sf_from_url(url)
  sf$composite_id = composite_id
  sf
}

all_trees = map2(
  names(ttops_filepaths),
  unlist(ttops_filepaths),
  possibly(download_ttops, otherwise = NULL)
)

all_trees = compact(all_trees)
all_trees = map(all_trees, force_all_cols_to_character)
all_trees = bind_rows(all_trees)

cat("Downloaded", nrow(all_trees), "detected trees from", length(ttops_filepaths), "composites\n")

st_write(all_trees, ALL_DETECTED_TREES_FILEPATH, delete_dsn = TRUE)
cat("Wrote all detected trees to", ALL_DETECTED_TREES_FILEPATH, "\n")
