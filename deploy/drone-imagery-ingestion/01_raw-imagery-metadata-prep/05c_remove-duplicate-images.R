# Purpose: Remove duplicate images from parsed EXIF metadata CSVs.
# For each duplicate group (within-mission), keep the image with the alphabetically first
# image_path_contrib. If ties, keep the first row. Remove the rest.
# Writes deduplicated CSVs to a new folder, preserving the originals.

library(tidyverse)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")

## Workflow

# Determine which missions to process
missions_to_process = read_csv(MISSIONS_TO_PROCESS_RAW_IMAGERY_METADATA_LIST_PATH) |>
  pull(mission_id)

# Read the within-mission duplicates list
duplicates_filepath = file.path(
  DUPLICATE_IMAGES_PATH,
  paste0(PROJECT_NAME_TO_PROCESS_RAW_IMAGERY_METADATA, "_duplicate-images_within-mission.csv")
)

if (!file.exists(duplicates_filepath)) {
  stop("Duplicates file not found. Run 05b_find-duplicate-images.R first.")
}

duplicates = read_csv(duplicates_filepath, show_col_types = FALSE)

# For each duplicate group, identify which images to remove
# Keep the one with alphabetically first image_path_contrib (first row if ties)
if (nrow(duplicates) > 0) {
  images_to_remove = duplicates |>
    group_by(duplicate_group_id) |>
    arrange(image_path_contrib) |>
    slice(-1) |>
    ungroup()

  cat("Total duplicate images:", nrow(duplicates), "\n")
  cat("Images to remove:", nrow(images_to_remove), "\n")
} else {
  images_to_remove = duplicates[0, ]
  cat("No duplicates found.\n")
}

# Create output directory
create_dir(PARSED_EXIF_DEDUPED_PATH)

# Track removals for reporting
removal_report <- list()

# Process each mission
for (mission_id_foc in missions_to_process) {
  input_filepath = file.path(PARSED_EXIF_METADATA_PATH, paste0(mission_id_foc, ".csv"))
  output_filepath = file.path(PARSED_EXIF_DEDUPED_PATH, paste0(mission_id_foc, ".csv"))

  if (!file.exists(input_filepath)) {
    warning(paste("Metadata file not found for mission:", mission_id_foc))
    next
  }

  # Read mission metadata
  mission_metadata = read_csv(input_filepath, show_col_types = FALSE)
  n_before = nrow(mission_metadata)

  # Get image_ids to remove for this mission
  ids_to_remove = images_to_remove |>
    filter(mission_id == mission_id_foc) |>
    pull(image_id)

  # Remove duplicate images
  mission_metadata_cleaned = mission_metadata |>
    filter(!(image_id %in% ids_to_remove))

  n_after = nrow(mission_metadata_cleaned)
  n_removed = n_before - n_after

  if (n_removed > 0) {
    removal_report[[mission_id_foc]] = n_removed
  }

  # Write to new location
  write_csv(mission_metadata_cleaned, output_filepath)
}

# Report missions with duplicates removed
if (length(removal_report) > 0) {
  cat("\n-- Missions with duplicates removed --\n")
  for (mission_id_foc in names(removal_report)) {
    n_removed <- removal_report[[mission_id_foc]]
    cat(mission_id_foc, ":", n_removed, "images removed\n")
  }
  cat("\nTotal:", length(removal_report), "missions affected,",
      sum(unlist(removal_report)), "images removed\n")
} else {
  cat("\nNo duplicates removed from any mission.\n")
}

cat("\nDeduplicated files written to:", PARSED_EXIF_DEDUPED_PATH, "\n")
