# Purpose: Combine all per-project within-mission duplicate log files into a single
# catalog-wide log. This merged log is used during curation to resolve image IDs
# that reference removed duplicates to their remaining counterparts.

library(tidyverse)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")

## Workflow

# List all within-mission duplicate files
within_mission_files = list.files(
  DUPLICATE_IMAGES_PATH,
  pattern = "_duplicate-images_within-mission\\.csv$",
  full.names = TRUE
)

if (length(within_mission_files) == 0) {
  stop("No within-mission duplicate files found in: ", DUPLICATE_IMAGES_PATH)
}

cat("Found", length(within_mission_files), "within-mission duplicate files\n")

# Read and combine all files, extracting project name from filename
read_duplicate_file = function(filepath) {
  # Extract project name from filename (e.g., "2019-focal_duplicate-images_within-mission.csv" -> "2019-focal")
  filename = basename(filepath)
  project_name = str_replace(filename, "_duplicate-images_within-mission\\.csv$", "")

  df = read_csv(filepath, show_col_types = FALSE)

  if (nrow(df) == 0) {
    return(NULL)
  }

  df |>
    mutate(project_name = project_name)
}

# Read all files and combine
cat("Reading and combining duplicate logs...\n")
combined_duplicates = map(within_mission_files, read_duplicate_file) |>
  bind_rows()

if (nrow(combined_duplicates) == 0) {
  cat("No duplicate images found across all projects. Creating empty output file.\n")
  # Create empty tibble with expected columns
  combined_duplicates = tibble(
    duplicate_group_id = character(),
    project_name = character(),
    image_id = character(),
    mission_id = character(),
    sub_mission_id = character(),
    lat = double(),
    lon = double(),
    file_size_gb = double(),
    image_path_contrib = character()
  )
} else {
  cat("Total duplicate images across all projects:", nrow(combined_duplicates), "\n")

  # Create catalog-wide unique duplicate_group_id by combining mission_id with
  # the per-project duplicate_group_id (which is unique within a mission)
  combined_duplicates = combined_duplicates |>
    mutate(duplicate_group_id = paste(mission_id, duplicate_group_id, sep = "_")) |>
    arrange(duplicate_group_id, image_path_contrib) |>
    select(
      duplicate_group_id,
      project_name,
      image_id,
      mission_id,
      sub_mission_id,
      lat,
      lon,
      file_size_gb,
      image_path_contrib,
      image_path_ofo
    )

  cat("Number of unique duplicate groups:", n_distinct(combined_duplicates$duplicate_group_id), "\n")
}

# Ensure output directory exists
create_dir(dirname(CATALOG_DUPLICATE_IMAGES_LOG_PATH))

# Write output
write_csv(combined_duplicates, CATALOG_DUPLICATE_IMAGES_LOG_PATH)
cat("Written to:", CATALOG_DUPLICATE_IMAGES_LOG_PATH, "\n")