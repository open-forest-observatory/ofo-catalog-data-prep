# Purpose: Find and record all images across all missions in a project that are duplicates,
# defined as having identical geographic coordinates (lat, lon) and the same file size (file_size_gb).
# These are the same image duplicated with a different file name.
# Produces two outputs:
#   1. Within-mission duplicates: duplicates that occur within the same mission
#   2. Cross-mission duplicates: duplicates across all missions regardless of mission ID

library(tidyverse)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")

## Workflow

# Determine which missions to process
missions_to_process = read_csv(MISSIONS_TO_PROCESS_RAW_IMAGERY_METADATA_LIST_PATH) |>
  pull(mission_id)

# Read and combine all parsed EXIF metadata for missions in the project
read_mission_metadata = function(mission_id_foc) {
  filepath = file.path(PARSED_EXIF_METADATA_PATH, paste0(mission_id_foc, ".csv"))
  if (file.exists(filepath)) {
    return(read_csv(filepath, show_col_types = FALSE))
  } else {
    warning(paste("Metadata file not found for mission:", mission_id_foc))
    return(NULL)
  }
}

# Read all mission metadata
all_metadata = map(missions_to_process, read_mission_metadata) |>
  bind_rows()

cat("Total images across all missions:", nrow(all_metadata), "\n")

# Filter to images with valid coordinates and file size
valid_metadata = all_metadata |>
  filter(!is.na(lat) & !is.na(lon) & !is.na(file_size_gb))

create_dir(DUPLICATE_IMAGES_PATH)

## 1. Within-mission duplicates
# Find duplicates within each mission (same mission_id, lat, lon, file_size_gb)
duplicates_within_mission = valid_metadata |>
  group_by(mission_id, lat, lon, file_size_gb) |>
  filter(n() > 1) |>
  mutate(duplicate_group_id = cur_group_id()) |>
  ungroup() |>
  arrange(mission_id, lat, lon, file_size_gb, image_id) |>
  select(duplicate_group_id, everything())

cat("\n-- Within-mission duplicates --\n")
cat("Number of duplicate images:", nrow(duplicates_within_mission), "\n")
cat("Number of unique duplicate groups:",
    n_distinct(duplicates_within_mission$duplicate_group_id), "\n")

output_filename_within = paste0(
  PROJECT_NAME_TO_PROCESS_RAW_IMAGERY_METADATA, "_duplicate-images_within-mission.csv"
)
output_filepath_within = file.path(DUPLICATE_IMAGES_PATH, output_filename_within)
write_csv(duplicates_within_mission, output_filepath_within)
cat("Written to:", output_filepath_within, "\n")

## 2. Cross-mission duplicates
# Find duplicates across all missions (same lat, lon, file_size_gb regardless of mission)
duplicates_cross_mission = valid_metadata |>
  group_by(lat, lon, file_size_gb) |>
  filter(n() > 1) |>
  mutate(duplicate_group_id = cur_group_id()) |>
  ungroup() |>
  arrange(lat, lon, file_size_gb, mission_id, image_id) |>
  select(duplicate_group_id, everything())

cat("\n-- Cross-mission duplicates --\n")
cat("Number of duplicate images:", nrow(duplicates_cross_mission), "\n")
cat("Number of unique duplicate groups:",
    n_distinct(duplicates_cross_mission$duplicate_group_id), "\n")

output_filename_cross = paste0(
  PROJECT_NAME_TO_PROCESS_RAW_IMAGERY_METADATA, "_duplicate-images_cross-mission.csv"
)
output_filepath_cross = file.path(DUPLICATE_IMAGES_PATH, output_filename_cross)
write_csv(duplicates_cross_mission, output_filepath_cross)
cat("Written to:", output_filepath_cross, "\n")
