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
# Load catalog-wide duplicate log
# ============================================================================

cat("Loading catalog-wide duplicate log...\n")

if (!file.exists(CATALOG_DUPLICATE_IMAGES_LOG_PATH)) {
  stop("Catalog duplicate log not found. Run 01a_finalize-raw-imagery-metadata-prep/merge-duplicate-logs.R first: ",
       CATALOG_DUPLICATE_IMAGES_LOG_PATH)
}

duplicate_log = read_csv(CATALOG_DUPLICATE_IMAGES_LOG_PATH, col_types = cols(.default = "c"))
cat(sprintf("  Loaded %d duplicate image records\n", nrow(duplicate_log)))

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
      images = st_read(filepath, quiet = TRUE)
      # Extract lon/lat from geometry
      coords = st_coordinates(images)
      images$lon = coords[, 1]
      images$lat = coords[, 2]
      return(images)
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

    # Handle missions with minor ID swaps (<=4 differences):
    # Between the original curation and saving of the original curation image metadata,
    # some missions experienced accidental corruption where up to 4 images got swapped.
    # Since the swapping has been undone, the current IDs are correct for missions where
    # most IDs already correspond. For these missions, we treat the current IDs as
    # authoritative (matching the state at the time of curation, before the swap).
    # Missions with more extensive ID changes cannot be corrected this way.
    cat("Checking for missions with minor ID swaps (<=4 differences)...\n")
    for (mission_id in names(crosswalks)) {
      crosswalk = crosswalks[[mission_id]]
      if (!is.null(crosswalk) && nrow(crosswalk) > 0) {
        n_differences = sum(crosswalk$former_image_id != crosswalk$current_image_id)
        if (n_differences > 0 && n_differences <= 4) {
          cat(sprintf("  Mission %s: %d minor ID differences detected, treating current IDs as authoritative\n",
                      mission_id, n_differences))
          crosswalks[[mission_id]] = crosswalk |>
            mutate(former_image_id = current_image_id)
        }
      }
    }
  }
}

# ============================================================================
# Update curation notes for formerly curated missions
# ============================================================================

cat("Updating curation notes with current image IDs...\n")

update_curation_row = function(row, crosswalks, formerly_curated_missions, duplicate_log) {
  mission_id = row$mission_id
  is_formerly_curated = mission_id %in% formerly_curated_missions

  # Get crosswalk if this is a formerly curated mission
  crosswalk = NULL
  if (is_formerly_curated) {
    crosswalk = crosswalks[[mission_id]]
    if (is.null(crosswalk) || nrow(crosswalk) == 0) {
      warning(paste("No crosswalk available for formerly curated mission", mission_id))
      # Continue anyway - duplicate resolution can still be applied
    }
  }

  # Anomaly columns to update (extraneous_images handled separately in Step 9)
  anomaly_cols = c("collection_time_anomaly", "altitude_anomaly",
                   "spatial_anomaly", "camera_pitch_anomaly", "excess_images_anomaly",
                   "missing_images_anomaly", "other_anomalies", "anomaly_notes")

  for (col in anomaly_cols) {
    if (col %in% names(row) && !is.na(row[[col]])) {
      row[[col]] = update_anomaly_image_references(row[[col]], crosswalk, duplicate_log)
    }
  }

  return(row)
}

# Update each row of curation notes
curation_notes_updated = curation_notes
for (i in seq_len(nrow(curation_notes))) {
  row_list = as.list(curation_notes[i, ])
  updated = update_curation_row(row_list, crosswalks, formerly_curated_missions, duplicate_log)
  for (col in names(updated)) {
    curation_notes_updated[i, col] = updated[[col]]
  }
}

# ============================================================================
# Parse extraneous images and create exclusion lists per mission
# (with crosswalk application and duplicate resolution)
# ============================================================================

cat("Parsing extraneous images, applying crosswalks, and resolving duplicates...\n")

exclusion_lists = list()

# Helper to crosswalk and resolve an image ID
crosswalk_and_resolve = function(image_id, crosswalk, duplicate_log, current_image_ids) {
  working_id = image_id

  # Step 1: Apply crosswalk if available (formerly curated missions)
  if (!is.null(crosswalk) && nrow(crosswalk) > 0) {
    crosswalk_result = crosswalk |>
      filter(former_image_id == image_id) |>
      pull(current_image_id)

    if (length(crosswalk_result) > 0 && crosswalk_result[1] != "none" && crosswalk_result[1] != "multiple") {
      working_id = crosswalk_result[1]
    } else if (length(crosswalk_result) > 0 && crosswalk_result[1] == "none") {
      # Crosswalk says image doesn't exist in current data.
      # No duplicate lookup needed: if duplicates existed at this location,
      # one would remain after deduplication and the crosswalk would have found it.
      return(NA_character_)
    }
    # If "multiple" or not found in crosswalk, continue with original ID
  }

  # Step 2: If the working ID exists in current data, use it directly
  if (working_id %in% current_image_ids) {
    return(working_id)
  }

  # Step 3: Otherwise, try to find its remaining duplicate
  remaining = get_remaining_duplicate(working_id, duplicate_log)

  if (!is.na(remaining) && remaining %in% current_image_ids) {
    return(remaining)
  }

  # Neither exists - return NA (will be filtered out)
  return(NA_character_)
}

for (i in seq_len(nrow(curation_notes_updated))) {
  row = curation_notes_updated[i, ]
  mission_id = row$mission_id

  if (!is.na(row$extraneous_images) && row$extraneous_images != "") {
    tryCatch({
      # Parse extraneous images (handles range expansion)
      excluded_ids = parse_extraneous_images(row$extraneous_images, mission_id)

      if (length(excluded_ids) > 0) {
        # Load current image IDs for this mission to check existence
        current_filepath = file.path(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
                                     paste0(mission_id, "_image-metadata.gpkg"))
        if (file.exists(current_filepath)) {
          current_images = st_read(current_filepath, quiet = TRUE)
          current_image_ids = current_images$image_id

          # Get crosswalk if this is a formerly curated mission
          crosswalk = NULL
          is_formerly_curated = mission_id %in% formerly_curated_missions
          if (is_formerly_curated) {
            crosswalk = crosswalks[[mission_id]]
          }

          # Apply crosswalk and resolve duplicates for each ID
          resolved_ids_raw = map_chr(excluded_ids, ~crosswalk_and_resolve(.x, crosswalk, duplicate_log, current_image_ids))

          # Count statistics before filtering
          n_missing = sum(is.na(resolved_ids_raw))
          n_changed = sum(!is.na(resolved_ids_raw) & excluded_ids != resolved_ids_raw)

          # Remove NAs (images that don't exist) and deduplicate
          resolved_ids = unique(resolved_ids_raw[!is.na(resolved_ids_raw)])

          if (n_missing > 0) {
            cat(sprintf("  Mission %s: %d image(s) not found in current data (skipped)\n",
                        mission_id, n_missing))
          }

          if (n_changed > 0) {
            cat(sprintf("  Mission %s: resolved %d image(s) via crosswalk/duplicates\n",
                        mission_id, n_changed))
          }

          if (length(resolved_ids) > 0) {
            if (!(mission_id %in% names(exclusion_lists))) {
              exclusion_lists[[mission_id]] = character(0)
            }
            exclusion_lists[[mission_id]] = unique(c(exclusion_lists[[mission_id]], resolved_ids))

            cat(sprintf("  Mission %s: %d images to exclude\n", mission_id, length(resolved_ids)))
          }
        }
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
future::plan(multisession(workers = future::availableCores()*3))
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
