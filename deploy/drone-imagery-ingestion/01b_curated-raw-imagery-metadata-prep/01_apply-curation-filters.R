# deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/01_apply-curation-filters.R
# Purpose: Apply curation filters to remove extraneous images and update image IDs
# for missions that were curated using former metadata.
#
# Workflow:
# 1. Load curation notes and combine duplicate rows
# 2. For formerly-curated missions, create spatial crosswalk and update image ID references
# 3. Parse extraneous_images into lists of image IDs to exclude
# 4. Filter pre-curation image metadata to remove extraneous images
# 5. Output filtered image metadata per mission

library(tidyverse)
library(sf)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/curation-utils.R")
source("src/utils.R")

# ============================================================================
# Load curation notes and handle duplicates
# ============================================================================

cat("Loading and combining curation notes...\n")
curation_notes = load_and_combine_curation_notes(CURATION_NOTES_FILEPATH)

cat(sprintf("  Loaded %d sub-mission curation records\n", nrow(curation_notes)))

# ============================================================================
# Load formerly curated mission list and metadata
# ============================================================================

cat("Loading formerly curated mission data...\n")

# Check if formerly curated files exist
if (!file.exists(FORMERLY_CURATED_MISSION_LIST)) {
  cat("  No formerly curated mission list found - skipping crosswalk creation\n")
  formerly_curated_missions = character(0)
  formerly_curated_image_metadata = NULL
} else {
  formerly_curated_missions = read_csv(FORMERLY_CURATED_MISSION_LIST, col_types = cols(.default = "c")) |>
    pull(mission_id) |>
    unique()

  cat(sprintf("  Found %d formerly curated missions\n", length(formerly_curated_missions)))

  if (!file.exists(FORMERLY_CURATED_IMAGE_METADATA)) {
    warning("Formerly curated mission list exists but image metadata file not found")
    formerly_curated_image_metadata = NULL
  } else {
    formerly_curated_image_metadata = st_read(FORMERLY_CURATED_IMAGE_METADATA, quiet = TRUE)
  }
}

# ============================================================================
# Get list of all missions in current pre-curation metadata
# ============================================================================

current_image_files = list.files(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH, pattern = "\\.gpkg$", full.names = TRUE)
current_mission_ids = str_extract(basename(current_image_files), "^\\d{6}")

cat(sprintf("Found %d missions with pre-curation metadata\n", length(current_mission_ids)))

# ============================================================================
# Create image ID crosswalks for formerly curated missions
# ============================================================================

crosswalks = list()

if (length(formerly_curated_missions) > 0 && !is.null(formerly_curated_image_metadata)) {
  cat("Creating image ID crosswalks for formerly curated missions...\n")

  # Helper function to load current images for a mission
  load_current_images = function(mission_id) {
    filepath = file.path(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH, paste0(mission_id, "_image-metadata.gpkg"))
    if (file.exists(filepath)) {
      return(st_read(filepath, quiet = TRUE))
    }
    return(NULL)
  }

  # Only load current images for formerly curated missions that exist in current data
  missions_needing_crosswalk = intersect(formerly_curated_missions, current_mission_ids)

  if (length(missions_needing_crosswalk) > 0) {
    # Load all current images for missions needing crosswalk
    formerly_curated_current = missions_needing_crosswalk |>
      map(load_current_images) |>
      compact() |>
      bind_rows()

    # Create crosswalks for each formerly curated mission
    for (mission_id in missions_needing_crosswalk) {
      cat(sprintf("  Creating crosswalk for mission %s...\n", mission_id))
      crosswalks[[mission_id]] = create_image_id_crosswalk(
        formerly_curated_image_metadata,
        formerly_curated_current,
        mission_id
      )
    }
  }
}

# ============================================================================
# Update curation notes for formerly curated missions
# ============================================================================

cat("Updating curation notes with current image IDs...\n")

update_curation_row = function(row, crosswalks, formerly_curated_missions) {
  mission_id = row$mission_id

  if (!(mission_id %in% formerly_curated_missions)) {
    return(row)
  }

  crosswalk = crosswalks[[mission_id]]
  if (is.null(crosswalk) || nrow(crosswalk) == 0) {
    warning(paste("No crosswalk available for formerly curated mission", mission_id))
    return(row)
  }

  # Update all columns that might contain image ID references
  cols_to_update = c("extraneous_images", "collection_time_anomaly", "altitude_anomaly",
                     "spatial_anomaly", "camera_pitch_anomaly", "excess_images_anomaly",
                     "missing_images_anomaly", "other_anomalies", "anomaly_notes")

  for (col in cols_to_update) {
    if (col %in% names(row) && !is.na(row[[col]])) {
      row[[col]] = update_image_id_references(row[[col]], crosswalk)
    }
  }

  return(row)
}

# Update each row of curation notes
curation_notes_updated = curation_notes
for (i in seq_len(nrow(curation_notes))) {
  row_list = as.list(curation_notes[i, ])
  updated = update_curation_row(row_list, crosswalks, formerly_curated_missions)
  for (col in names(updated)) {
    curation_notes_updated[i, col] = updated[[col]]
  }
}

# ============================================================================
# Parse extraneous images and create exclusion lists per mission
# ============================================================================

cat("Parsing extraneous images...\n")

exclusion_lists = list()

for (i in seq_len(nrow(curation_notes_updated))) {
  row = curation_notes_updated[i, ]
  mission_id = row$mission_id

  if (!is.na(row$extraneous_images) && row$extraneous_images != "") {
    tryCatch({
      excluded_ids = parse_extraneous_images(row$extraneous_images, mission_id)

      if (length(excluded_ids) > 0) {
        if (!(mission_id %in% names(exclusion_lists))) {
          exclusion_lists[[mission_id]] = character(0)
        }
        exclusion_lists[[mission_id]] = unique(c(exclusion_lists[[mission_id]], excluded_ids))

        cat(sprintf("  Mission %s: %d images to exclude\n", mission_id, length(excluded_ids)))
      }
    }, error = function(e) {
      stop(paste("Error parsing extraneous_images for mission", mission_id, ":", e$message))
    })
  }
}

cat(sprintf("  Total: %d missions with images to exclude\n", length(exclusion_lists)))

# ============================================================================
# Filter image metadata and save to post-curation paths
# ============================================================================

cat("Filtering image metadata and saving post-curation versions...\n")

# Create output directory
create_dir(POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH)

# Process all current missions (applying curation filters where they exist)
all_missions_to_process = current_mission_ids

filter_and_save_mission = function(mission_id_foc) {
  # Load pre-curation image metadata
  input_filepath = file.path(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
                             paste0(mission_id_foc, "_image-metadata.gpkg"))

  if (!file.exists(input_filepath)) {
    warning(paste("No pre-curation metadata found for mission", mission_id_foc))
    return(FALSE)
  }

  image_metadata = st_read(input_filepath, quiet = TRUE)
  n_original = nrow(image_metadata)

  # Apply exclusion filter if this mission has exclusions
  if (mission_id_foc %in% names(exclusion_lists)) {
    excluded_ids = exclusion_lists[[mission_id_foc]]
    image_metadata = image_metadata |>
      filter(!(image_id %in% excluded_ids))
  }

  n_retained = nrow(image_metadata)
  n_excluded = n_original - n_retained

  if (n_excluded > 0) {
    cat(sprintf("  Mission %s: %d -> %d images (excluded %d)\n",
                mission_id_foc, n_original, n_retained, n_excluded))
  }

  # Save filtered image metadata
  output_filepath = file.path(POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
                              paste0(mission_id_foc, "_image-metadata.gpkg"))
  st_write(image_metadata, output_filepath, delete_dsn = TRUE, quiet = TRUE)

  return(TRUE)
}

# Process all missions in parallel
future::plan(multisession)
results = future_map_lgl(all_missions_to_process, filter_and_save_mission, .progress = TRUE)

cat(sprintf("\nFiltered %d missions successfully\n", sum(results)))

# ============================================================================
# Save updated curation notes for use by subsequent scripts
# ============================================================================

curation_notes_output_path = file.path(POST_CURATION_INTERMEDIATE_PATH,
                                        "curation-notes-processed.csv")
create_dir(dirname(curation_notes_output_path))
write_csv(curation_notes_updated, curation_notes_output_path)

cat(sprintf("Saved processed curation notes to: %s\n", curation_notes_output_path))

cat("\n**** Post-curation filtering complete ****\n")
