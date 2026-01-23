# Implementation Plan: Incorporate Curation into Drone Data Ingestion Pipeline

## Overview

This plan adds post-curation functionality to the drone data ingestion pipeline. The curation step filters out extraneous images identified by human curators and incorporates anomaly notes into the mission metadata.

**Key workflow:**
```
Pre-curation metadata (existing)
    ↓
[NEW] Apply curation filters (remove extraneous images, update IDs for formerly-curated missions)
    ↓
[NEW] Re-summarize metadata at mission/sub-mission level
    ↓
[NEW] Merge curation anomaly notes into metadata
    ↓
[MODIFIED] Copy/symlink images to sorted folders (post-curation)
    ↓
[MODIFIED] Fix EXIF, create thumbnails/zip, upload to S3
```

---

## Phase 1: Update Constants File

### Step 1.1: Add missing and new path constants to `00_set-constants.R`

Add the following constants:

```r
# Missing constant (used by scripts but not defined)
SORTED_IMAGERY_PATH = file.path(RAW_IMAGERY_INGESTION_PATH, "2_sorted")

# Post-curation intermediate paths (mirror pre-curation structure)
POST_CURATION_INTERMEDIATE_PATH = file.path(RAW_IMAGERY_INGESTION_PATH, "metadata/5_post-curation-intermediate")

POST_CURATION_IMAGE_EXIF_W_SORTING_PLAN_PATH = file.path(POST_CURATION_INTERMEDIATE_PATH, "1_exif-w-sorting-plan")
POST_CURATION_DERIVED_METADATA_PER_MISSION_PATH = file.path(POST_CURATION_INTERMEDIATE_PATH, "2_derived-metadata-per-mission")
POST_CURATION_DERIVED_METADATA_PER_SUB_MISSION_PATH = file.path(POST_CURATION_INTERMEDIATE_PATH, "3_derived-metadata-per-sub-mission")

# Post-curation final paths
POST_CURATION_FINAL_PATH = file.path(RAW_IMAGERY_INGESTION_PATH, "metadata/6_post-curation-final")

POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH = file.path(POST_CURATION_FINAL_PATH, "1_parsed-exif-per-image")
POST_CURATION_FULL_METADATA_PER_MISSION_PATH = file.path(POST_CURATION_FINAL_PATH, "2_full-metadata-per-mission")
POST_CURATION_FULL_METADATA_PER_SUB_MISSION_PATH = file.path(POST_CURATION_FINAL_PATH, "3_full-metadata-per-sub-mission")

# Combined post-curation metadata files
POST_CURATION_FULL_METADATA_PER_MISSION_COMBINED_FILEPATH = file.path(POST_CURATION_FINAL_PATH, "ofo-all-missions-metadata-curated.gpkg")
POST_CURATION_FULL_METADATA_PER_IMAGE_COMBINED_FILEPATH = file.path(POST_CURATION_FINAL_PATH, "ofo-all-images-metadata-curated.gpkg")
```

**File:** `deploy/drone-imagery-ingestion/00_set-constants.R`

**Lines to add after line 70** (after the existing curation constants)

---

## Phase 2: Create Shared Utility Functions

### Step 2.1: Create `src/curation-utils.R` with extraneous image parsing functions

This module handles parsing the various extraneous_images notation formats.

```r
# src/curation-utils.R
# Utility functions for parsing and applying curation notes

library(tidyverse)
library(sf)

#' Parse extraneous_images string into a vector of image IDs
#'
#' Handles formats:
#' - Single: "000030-01_001183"
#' - Comma-separated: "000036-01_000437, 000036-01_000436"
#' - Ranges: "000034-01_000912 to 000034-01_000933"
#' - Reverse ranges: "000032-01_003332 to 000032-01_003255"
#' - Mixed: all of the above combined
#'
#' @param extraneous_string Character string from curation notes
#' @param expected_mission_id The mission ID these images should belong to (for validation)
#' @return Character vector of image IDs to exclude
parse_extraneous_images = function(extraneous_string, expected_mission_id = NULL) {
  if (is.na(extraneous_string) || extraneous_string == "" || is.null(extraneous_string)) {
    return(character(0))
  }


  # Trim whitespace
  extraneous_string = trimws(extraneous_string)

  # Split by comma (with optional surrounding whitespace)
  parts = str_split(extraneous_string, ",\\s*")[[1]]
  parts = trimws(parts)

  all_image_ids = character(0)

  for (part in parts) {
    if (grepl(" to ", part, fixed = TRUE)) {
      # This is a range
      range_parts = str_split(part, "\\s+to\\s+")[[1]]
      if (length(range_parts) != 2) {
        warning(paste("Invalid range format:", part))
        next
      }

      start_id = trimws(range_parts[1])
      end_id = trimws(range_parts[2])

      # Extract mission-submission prefix and numeric suffix
      start_match = str_match(start_id, "^(\\d{6}-\\d{2})_(\\d+)$")
      end_match = str_match(end_id, "^(\\d{6}-\\d{2})_(\\d+)$")

      if (is.na(start_match[1]) || is.na(end_match[1])) {
        warning(paste("Invalid image ID format in range:", part))
        next
      }

      start_prefix = start_match[2]
      end_prefix = end_match[2]

      if (start_prefix != end_prefix) {
        stop(paste("Range spans different sub-missions:", part,
                   "- this appears to be a typo in the curation notes."))
      }

      # Validate mission ID if provided
      if (!is.null(expected_mission_id)) {
        range_mission = substr(start_prefix, 1, 6)
        if (range_mission != expected_mission_id) {
          stop(paste("Cross-mission reference detected:", part,
                     "references mission", range_mission,
                     "but curation row is for mission", expected_mission_id,
                     "- this appears to be a typo in the curation notes."))
        }
      }

      start_num = as.integer(start_match[3])
      end_num = as.integer(end_match[3])

      # Handle reverse ranges
      if (start_num > end_num) {
        num_seq = seq(end_num, start_num)
      } else {
        num_seq = seq(start_num, end_num)
      }

      # Generate all IDs in range
      range_ids = paste0(start_prefix, "_", str_pad(num_seq, 6, pad = "0"))
      all_image_ids = c(all_image_ids, range_ids)

    } else {
      # Single image ID
      single_id = trimws(part)

      # Validate format
      if (!grepl("^\\d{6}-\\d{2}_\\d+$", single_id)) {
        warning(paste("Invalid image ID format:", single_id))
        next
      }

      # Validate mission ID if provided
      if (!is.null(expected_mission_id)) {
        id_mission = substr(single_id, 1, 6)
        if (id_mission != expected_mission_id) {
          stop(paste("Cross-mission reference detected:", single_id,
                     "references mission", id_mission,
                     "but curation row is for mission", expected_mission_id,
                     "- this appears to be a typo in the curation notes."))
        }
      }

      all_image_ids = c(all_image_ids, single_id)
    }
  }

  return(unique(all_image_ids))
}


#' Load and combine curation notes, handling duplicate rows
#'
#' @param curation_filepath Path to curation notes CSV
#' @return Tibble with one row per sub-mission, combined extraneous_images and anomaly notes
load_and_combine_curation_notes = function(curation_filepath) {
  curation_raw = read_csv(curation_filepath, col_types = cols(.default = "c"))

  # Define anomaly columns to aggregate
  anomaly_cols = c("anomaly_severity", "collection_time_anomaly", "altitude_anomaly",
                   "spatial_anomaly", "camera_pitch_anomaly", "excess_images_anomaly",
                   "missing_images_anomaly", "other_anomalies", "anomaly_notes")

  # Group by mission_id and sub_mission_id, combine values
  curation_combined = curation_raw |>
    group_by(mission_id, sub_mission_id) |>
    summarise(
      # Combine extraneous_images with comma separator
      extraneous_images = paste(na.omit(extraneous_images[extraneous_images != ""]), collapse = ", "),
      # For anomaly columns, take first non-NA/non-empty or combine with comma
      across(all_of(anomaly_cols), ~ paste(na.omit(.x[.x != ""]), collapse = ", ")),
      .groups = "drop"
    ) |>
    # Replace empty strings with NA for cleaner handling
    mutate(across(everything(), ~ if_else(.x == "", NA_character_, .x)))

  return(curation_combined)
}


#' Map image IDs from formerly curated metadata to current metadata using spatial join
#'
#' @param formerly_curated_images SF object with formerly curated image points
#' @param current_images SF object with current image points
#' @param mission_id Mission ID to filter on
#' @param distance_threshold Maximum distance (meters) for spatial matching
#' @return Tibble with columns: former_image_id, current_image_id
create_image_id_crosswalk = function(formerly_curated_images, current_images,
                                      mission_id, distance_threshold = 1) {
  # Filter to focal mission
  former = formerly_curated_images |> filter(mission_id == !!mission_id)
  current = current_images |> filter(mission_id == !!mission_id)

  if (nrow(former) == 0 || nrow(current) == 0) {
    warning(paste("No images found for mission", mission_id, "in one or both datasets"))
    return(tibble(former_image_id = character(), current_image_id = character()))
  }

  # Transform to UTM for accurate distance calculation
  former_utm = former |> st_transform(st_crs(paste0("+proj=utm +zone=",
                                                     round((st_coordinates(former)[1,1] + 180) / 6) + 1)))
  current_utm = current |> st_transform(st_crs(former_utm))

  # Find nearest current image for each former image
  nearest_idx = st_nearest_feature(former_utm, current_utm)
  distances = st_distance(former_utm, current_utm[nearest_idx,], by_element = TRUE)

  crosswalk = tibble(
    former_image_id = former$image_id,
    current_image_id = current$image_id[nearest_idx],
    distance_m = as.numeric(distances)
  )

  # Warn about matches beyond threshold
  far_matches = crosswalk |> filter(distance_m > distance_threshold)
  if (nrow(far_matches) > 0) {
    warning(paste("Found", nrow(far_matches), "image matches in mission", mission_id,
                  "with distance >", distance_threshold, "meters. Max distance:",
                  round(max(far_matches$distance_m), 2), "m"))
  }

  return(crosswalk)
}


#' Update image ID references in a string using crosswalk
#'
#' @param text String potentially containing image ID references
#' @param crosswalk Tibble with former_image_id and current_image_id columns
#' @return String with updated image IDs
update_image_id_references = function(text, crosswalk) {
  if (is.na(text) || text == "") {
    return(text)
  }

  # Find all image ID patterns in the text
  pattern = "\\d{6}-\\d{2}_\\d{6}"
  matches = str_extract_all(text, pattern)[[1]]

  if (length(matches) == 0) {
    return(text)
  }

  # Replace each match with its current equivalent
  result = text
  for (old_id in matches) {
    new_id = crosswalk |> filter(former_image_id == old_id) |> pull(current_image_id)
    if (length(new_id) == 1) {
      result = str_replace(result, fixed(old_id), new_id)
    }
  }

  return(result)
}
```

**File:** `src/curation-utils.R`

---

## Phase 3: Create Post-Curation Metadata Scripts

Create new folder: `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/`

### Step 3.1: Create `01_apply-curation-filters.R`

This script:
1. Loads curation notes and combines duplicate rows
2. For formerly-curated missions, creates spatial crosswalk and updates image ID references
3. Parses extraneous_images into lists of image IDs to exclude
4. Filters pre-curation image metadata to remove extraneous images
5. Outputs filtered image metadata per mission

```r
# deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/01_apply-curation-filters.R
# Purpose: Apply curation filters to remove extraneous images and update image IDs
# for missions that were curated using former metadata.

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
formerly_curated_missions = read_csv(FORMERLY_CURATED_MISSION_LIST, col_types = cols(.default = "c")) |>
  pull(mission_id) |>
  unique()

cat(sprintf("  Found %d formerly curated missions\n", length(formerly_curated_missions)))

formerly_curated_image_metadata = st_read(FORMERLY_CURATED_IMAGE_METADATA, quiet = TRUE)

# ============================================================================
# Create image ID crosswalks for formerly curated missions
# ============================================================================

cat("Creating image ID crosswalks for formerly curated missions...\n")

# Get list of all missions in current pre-curation metadata
current_image_files = list.files(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH, pattern = "\\.gpkg$", full.names = TRUE)
current_mission_ids = str_extract(basename(current_image_files), "^\\d{6}")

# Load all current image metadata into one sf object for crosswalk creation
load_current_images = function(mission_id) {
  filepath = file.path(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH, paste0(mission_id, "_image-metadata.gpkg"))
  if (file.exists(filepath)) {
    return(st_read(filepath, quiet = TRUE))
  }
  return(NULL)
}

# Only load current images for formerly curated missions
formerly_curated_current = formerly_curated_missions |>
  intersect(current_mission_ids) |>
  map(load_current_images) |>
  compact() |>
  bind_rows()

# Create crosswalks for each formerly curated mission
crosswalks = list()
for (mission_id in intersect(formerly_curated_missions, current_mission_ids)) {
  cat(sprintf("  Creating crosswalk for mission %s...\n", mission_id))
  crosswalks[[mission_id]] = create_image_id_crosswalk(
    formerly_curated_image_metadata,
    formerly_curated_current,
    mission_id
  )
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

# Actually, let's do this more simply:
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
  sub_mission_id = paste0(mission_id, "-", str_pad(row$sub_mission_id, 2, pad = "0"))

  if (!is.na(row$extraneous_images) && row$extraneous_images != "") {
    tryCatch({
      excluded_ids = parse_extraneous_images(row$extraneous_images, mission_id)

      if (!(mission_id %in% names(exclusion_lists))) {
        exclusion_lists[[mission_id]] = character(0)
      }
      exclusion_lists[[mission_id]] = unique(c(exclusion_lists[[mission_id]], excluded_ids))

      cat(sprintf("  Mission %s: %d images to exclude\n", mission_id, length(excluded_ids)))
    }, error = function(e) {
      stop(paste("Error parsing extraneous_images for mission", mission_id, ":", e$message))
    })
  }
}

# ============================================================================
# Filter image metadata and save to post-curation paths
# ============================================================================

cat("Filtering image metadata and saving post-curation versions...\n")

# Create output directories
create_dir(POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH)
create_dir(POST_CURATION_IMAGE_EXIF_W_SORTING_PLAN_PATH)

# Get all missions to process (those with curation notes OR in formerly curated list)
missions_with_curation = unique(curation_notes_updated$mission_id)
all_missions_to_process = unique(c(missions_with_curation, current_mission_ids))

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

  # Also create filtered sorting plan (for script 08)
  sorting_plan_input = file.path(IMAGE_EXIF_W_SORTING_PLAN_PATH, paste0(mission_id_foc, ".csv"))
  if (file.exists(sorting_plan_input)) {
    sorting_plan = read_csv(sorting_plan_input, col_types = cols(.default = "c"))

    # Filter to retained images only
    # Need to match on image filename
    retained_image_ids = image_metadata$image_id
    sorting_plan_filtered = sorting_plan |>
      filter(paste0(mission_id, "-", sub_mission_id, "_",
                    str_pad(row_number(), 6, pad = "0")) %in% retained_image_ids |
             image_filename_out %in% paste0(retained_image_ids, ".JPG") |
             image_filename_out %in% paste0(retained_image_ids, ".jpg"))

    # Actually, let's match differently - extract image_id from filename
    sorting_plan = sorting_plan |>
      mutate(image_id_derived = tools::file_path_sans_ext(image_filename_out))

    sorting_plan_filtered = sorting_plan |>
      filter(image_id_derived %in% retained_image_ids)

    sorting_plan_output = file.path(POST_CURATION_IMAGE_EXIF_W_SORTING_PLAN_PATH,
                                    paste0(mission_id_foc, ".csv"))
    write_csv(sorting_plan_filtered |> select(-image_id_derived), sorting_plan_output)
  }

  return(TRUE)
}

# Process all missions
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
```

**File:** `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/01_apply-curation-filters.R`

---

### Step 3.2: Create `02_summarize-curated-metadata-per-mission.R`

This script re-runs the summarization logic (equivalent to pre-curation script 06) on the filtered images.

```r
# deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/02_summarize-curated-metadata-per-mission.R
# Purpose: Re-summarize EXIF metadata at mission and sub-mission level after curation filtering.
# This is equivalent to 01_raw-imagery-metadata-prep/06_summarize-exif-metadata-per-mission.R
# but operates on post-curation filtered images.

library(tidyverse)
library(sf)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/metadata-extraction_imagery_dataset.R")
source("src/utils.R")

# Reuse the polygon computation function from script 06
compute_polygons_and_images_retained = function(image_metadata, column_to_split_on, image_merge_distance) {
  split_metadata = split(image_metadata, image_metadata[[column_to_split_on]])
  dataset_ids = names(split_metadata)

  polygons_and_inds = purrr::map2(
    split_metadata,
    dataset_ids,
    extract_mission_polygon,
    image_merge_distance = image_merge_distance,
    identify_images_in_polygon = TRUE
  )

  intersection_image_ids = purrr::map(polygons_and_inds, "intersection_image_ids")
  image_counts = purrr::map(intersection_image_ids, length)
  too_few_images_bool = image_counts < 10

  if (any(too_few_images_bool)) {
    too_few_images_names = names(which(too_few_images_bool))
    warning(paste0("The following missions/sub-missions had fewer than 10 images retained: ",
                   paste(too_few_images_names, collapse = ", ")))
  }

  polygons_and_inds = polygons_and_inds[!too_few_images_bool]
  polygons_sfc = purrr::map(polygons_and_inds, "polygon")
  polygons_sf = purrr::map(polygons_sfc, st_as_sf)
  intersection_image_ids = purrr::map(polygons_and_inds, "intersection_image_ids")
  intersection_image_ids = unlist(intersection_image_ids, use.names = FALSE)

  return(list(polygons = polygons_sf, retained_image_IDs = intersection_image_ids))
}

# ============================================================================
# Workflow
# ============================================================================

# Get missions to process from post-curation filtered metadata
post_curation_files = list.files(POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
                                  pattern = "\\.gpkg$", full.names = FALSE)
missions_to_process = str_extract(post_curation_files, "^\\d{6}")

cat(sprintf("Processing %d missions...\n", length(missions_to_process)))

# Create output directories
create_dir(POST_CURATION_DERIVED_METADATA_PER_MISSION_PATH)
create_dir(POST_CURATION_DERIVED_METADATA_PER_SUB_MISSION_PATH)

summarize_curated_exif = function(mission_id_foc) {
  # Read post-curation filtered image metadata
  metadata_filepath = file.path(POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
                                paste0(mission_id_foc, "_image-metadata.gpkg"))
  image_metadata = st_read(metadata_filepath, quiet = TRUE)

  # Convert to data frame with lon/lat columns for compatibility
  coords = st_coordinates(image_metadata)
  image_metadata_df = image_metadata |>
    st_drop_geometry() |>
    mutate(lon = coords[, 1], lat = coords[, 2])

  # Compute polygons at mission level
  mission_res = tryCatch({
    compute_polygons_and_images_retained(
      image_metadata = image_metadata_df,
      column_to_split_on = "mission_id",
      image_merge_distance = IMAGE_MERGE_DISTANCE
    )
  }, error = function(e) {
    warning(paste("Failed to compute mission polygon for", mission_id_foc, ":", e$message))
    return(NULL)
  })

  if (is.null(mission_res)) return(FALSE)

  # Compute polygons at sub-mission level
  sub_mission_res = tryCatch({
    compute_polygons_and_images_retained(
      image_metadata = image_metadata_df,
      column_to_split_on = "sub_mission_id",
      image_merge_distance = IMAGE_MERGE_DISTANCE
    )
  }, error = function(e) {
    warning(paste("Failed to compute sub-mission polygons for", mission_id_foc, ":", e$message))
    return(NULL)
  })

  if (is.null(sub_mission_res)) return(FALSE)

  # Filter to images retained in both mission and sub-mission polygons
  images_retained_in_both = intersect(
    mission_res$retained_image_IDs,
    sub_mission_res$retained_image_IDs
  )

  image_metadata_df = image_metadata_df |> filter(image_id %in% images_retained_in_both)

  if (nrow(image_metadata_df) < 10) {
    warning(paste("Fewer than 10 images retained for mission", mission_id_foc, "after polygon filtering"))
    return(FALSE)
  }

  # Re-compute polygons with filtered images
  mission_res = compute_polygons_and_images_retained(
    image_metadata = image_metadata_df,
    column_to_split_on = "mission_id",
    image_merge_distance = IMAGE_MERGE_DISTANCE
  )

  sub_mission_res = compute_polygons_and_images_retained(
    image_metadata = image_metadata_df,
    column_to_split_on = "sub_mission_id",
    image_merge_distance = IMAGE_MERGE_DISTANCE
  )

  # Extract mission-level summary statistics
  summary_mission = extract_imagery_dataset_metadata(
    metadata = image_metadata_df,
    mission_polygon = mission_res$polygons[[1]],
    dataset_id = mission_id_foc
  )

  # Attribute polygon with metadata
  mission_poly = mission_res$polygons[[1]]
  mission_poly_attributed = bind_cols(mission_poly, summary_mission)
  mission_poly_attributed$mission_id = mission_id_foc

  # Write mission-level metadata
  metadata_per_mission_filepath = file.path(POST_CURATION_DERIVED_METADATA_PER_MISSION_PATH,
                                             paste0(mission_id_foc, ".gpkg"))
  st_write(mission_poly_attributed, metadata_per_mission_filepath, delete_dsn = TRUE, quiet = TRUE)

  # Process each sub-mission
  sub_mission_ids = unique(image_metadata_df$sub_mission_id)

  for (sub_mission_id_foc in sub_mission_ids) {
    image_metadata_sub = image_metadata_df |> filter(sub_mission_id == sub_mission_id_foc)

    if (!(sub_mission_id_foc %in% names(sub_mission_res$polygons))) {
      warning(paste("No polygon found for sub-mission", sub_mission_id_foc))
      next
    }

    polygon_sub = sub_mission_res$polygons[[sub_mission_id_foc]]

    summary_sub = extract_imagery_dataset_metadata(
      metadata = image_metadata_sub,
      mission_polygon = polygon_sub,
      dataset_id = sub_mission_id_foc
    )

    sub_poly_attributed = bind_cols(polygon_sub, summary_sub)
    sub_poly_attributed$mission_id = mission_id_foc
    sub_poly_attributed$sub_mission_id = sub_mission_id_foc

    metadata_per_sub_filepath = file.path(POST_CURATION_DERIVED_METADATA_PER_SUB_MISSION_PATH,
                                           paste0(sub_mission_id_foc, ".gpkg"))
    st_write(sub_poly_attributed, metadata_per_sub_filepath, delete_dsn = TRUE, quiet = TRUE)
  }

  return(TRUE)
}

# Run for each mission
future::plan(multisession)
results = future_walk(missions_to_process, summarize_curated_exif, .progress = TRUE)

cat("Summarization complete.\n")
```

**File:** `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/02_summarize-curated-metadata-per-mission.R`

---

### Step 3.3: Create `03_merge-metadata-and-curation-notes.R`

This script:
1. Merges derived EXIF metadata with Baserow contributed metadata (like script 07)
2. Adds curation anomaly columns to mission and sub-mission metadata

```r
# deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/03_merge-metadata-and-curation-notes.R
# Purpose: Merge post-curation derived metadata with Baserow metadata and curation notes.
# Adds anomaly columns from curation notes to mission and sub-mission metadata.

library(tidyverse)
library(sf)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")

# ============================================================================
# Load processed curation notes
# ============================================================================

curation_notes_path = file.path(POST_CURATION_INTERMEDIATE_PATH, "curation-notes-processed.csv")
curation_notes = read_csv(curation_notes_path, col_types = cols(.default = "c"))

# Prepare curation notes for joining
# At sub-mission level: use directly
# At mission level: aggregate across sub-missions

prepare_curation_for_sub_mission = function(curation_notes) {
  curation_notes |>
    mutate(
      sub_mission_id_full = paste0(mission_id, "-", str_pad(sub_mission_id, 2, pad = "0"))
    ) |>
    select(
      sub_mission_id = sub_mission_id_full,
      anomaly_severity,
      anomalies_collection_time = collection_time_anomaly,
      anomalies_altitude = altitude_anomaly,
      anomalies_spatial = spatial_anomaly,
      anomalies_camera_pitch = camera_pitch_anomaly,
      anomalies_excess_images = excess_images_anomaly,
      anomalies_missing_images = missing_images_anomaly,
      anomalies_other = other_anomalies,
      anomaly_notes
    )
}

prepare_curation_for_mission = function(curation_notes) {
  # Aggregate across sub-missions within each mission
  curation_notes |>
    group_by(mission_id) |>
    summarise(
      anomaly_severity = paste(na.omit(anomaly_severity[anomaly_severity != ""]), collapse = ", "),
      anomalies_collection_time = paste(na.omit(collection_time_anomaly[collection_time_anomaly != ""]), collapse = ", "),
      anomalies_altitude = paste(na.omit(altitude_anomaly[altitude_anomaly != ""]), collapse = ", "),
      anomalies_spatial = paste(na.omit(spatial_anomaly[spatial_anomaly != ""]), collapse = ", "),
      anomalies_camera_pitch = paste(na.omit(camera_pitch_anomaly[camera_pitch_anomaly != ""]), collapse = ", "),
      anomalies_excess_images = paste(na.omit(excess_images_anomaly[excess_images_anomaly != ""]), collapse = ", "),
      anomalies_missing_images = paste(na.omit(missing_images_anomaly[missing_images_anomaly != ""]), collapse = ", "),
      anomalies_other = paste(na.omit(other_anomalies[other_anomalies != ""]), collapse = ", "),
      anomaly_notes = paste(na.omit(anomaly_notes[anomaly_notes != ""]), collapse = ", "),
      .groups = "drop"
    ) |>
    # Replace empty strings with NA
    mutate(across(everything(), ~ if_else(.x == "", NA_character_, .x)))
}

curation_sub_mission = prepare_curation_for_sub_mission(curation_notes)
curation_mission = prepare_curation_for_mission(curation_notes)

# ============================================================================
# Get missions to process
# ============================================================================

post_curation_files = list.files(POST_CURATION_DERIVED_METADATA_PER_MISSION_PATH,
                                  pattern = "\\.gpkg$", full.names = FALSE)
missions_to_process = str_extract(post_curation_files, "^\\d{6}")

cat(sprintf("Merging metadata for %d missions...\n", length(missions_to_process)))

# Create output directories
create_dir(POST_CURATION_FULL_METADATA_PER_MISSION_PATH)
create_dir(POST_CURATION_FULL_METADATA_PER_SUB_MISSION_PATH)

# ============================================================================
# Merge function
# ============================================================================

merge_metadata_for_mission = function(mission_foc) {
  # Load derived (post-curation) metadata
  derived_mission_filepath = file.path(POST_CURATION_DERIVED_METADATA_PER_MISSION_PATH,
                                        paste0(mission_foc, ".gpkg"))
  derived_mission = st_read(derived_mission_filepath, quiet = TRUE)

  # Load contributed (Baserow) metadata
  baserow_mission_filepath = file.path(EXTRACTED_METADATA_PER_MISSION_PATH,
                                        paste0(mission_foc, ".csv"))

  if (!file.exists(baserow_mission_filepath)) {
    warning(paste("No Baserow metadata found for mission", mission_foc))
    return(FALSE)
  }

  baserow_mission = read_csv(baserow_mission_filepath, col_types = cols(.default = "c"))
  baserow_mission = baserow_mission |> rename(sub_mission_ids = sub_mission_id)

  # Merge derived and contributed
  full_metadata_mission = bind_cols(
    derived_mission |> select(-mission_id),
    baserow_mission
  ) |>
    select(!ends_with("_derived"), everything())

  # Add curation anomaly columns
  curation_row = curation_mission |> filter(mission_id == mission_foc)
  if (nrow(curation_row) == 1) {
    for (col in setdiff(names(curation_row), "mission_id")) {
      full_metadata_mission[[col]] = curation_row[[col]]
    }
  } else {
    # No curation notes for this mission - add NA columns
    for (col in setdiff(names(curation_mission), "mission_id")) {
      full_metadata_mission[[col]] = NA_character_
    }
  }

  # Write mission-level metadata
  output_filepath = file.path(POST_CURATION_FULL_METADATA_PER_MISSION_PATH,
                               paste0(mission_foc, "_mission-metadata.gpkg"))
  st_write(full_metadata_mission, output_filepath, delete_dsn = TRUE, quiet = TRUE)

  # Process sub-missions
  sub_mission_files = list.files(
    POST_CURATION_DERIVED_METADATA_PER_SUB_MISSION_PATH,
    pattern = paste0(mission_foc, "-\\d{2}\\.gpkg$"),
    full.names = TRUE
  )

  for (sub_file in sub_mission_files) {
    sub_mission_id_foc = str_extract(basename(sub_file), "\\d{6}-\\d{2}")

    derived_sub = st_read(sub_file, quiet = TRUE)

    baserow_sub_filepath = file.path(EXTRACTED_METADATA_PER_SUB_MISSION_PATH,
                                      paste0(sub_mission_id_foc, ".csv"))

    if (!file.exists(baserow_sub_filepath)) {
      warning(paste("No Baserow metadata found for sub-mission", sub_mission_id_foc))
      next
    }

    baserow_sub = read_csv(baserow_sub_filepath, col_types = cols(.default = "c"))

    full_metadata_sub = bind_cols(
      derived_sub |> select(-sub_mission_id, -mission_id),
      baserow_sub
    ) |>
      select(!ends_with("_derived"), everything())

    # Add curation anomaly columns
    curation_sub_row = curation_sub_mission |> filter(sub_mission_id == sub_mission_id_foc)
    if (nrow(curation_sub_row) == 1) {
      for (col in setdiff(names(curation_sub_row), "sub_mission_id")) {
        full_metadata_sub[[col]] = curation_sub_row[[col]]
      }
    } else {
      for (col in setdiff(names(curation_sub_mission), "sub_mission_id")) {
        full_metadata_sub[[col]] = NA_character_
      }
    }

    output_sub_filepath = file.path(POST_CURATION_FULL_METADATA_PER_SUB_MISSION_PATH,
                                     paste0(sub_mission_id_foc, "_sub-mission-metadata.gpkg"))
    st_write(full_metadata_sub, output_sub_filepath, delete_dsn = TRUE, quiet = TRUE)
  }

  return(TRUE)
}

# Run for all missions
future::plan(multisession)
results = future_map_lgl(missions_to_process, merge_metadata_for_mission, .progress = TRUE)

cat(sprintf("\nMerged metadata for %d missions successfully\n", sum(results)))
```

**File:** `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/03_merge-metadata-and-curation-notes.R`

---

### Step 3.4: Create control script `control_curated-metadata_01-to-03.R`

```r
# deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/control_curated-metadata_01-to-03.R
# Purpose: Run the post-curation metadata preparation pipeline

repo_root = "/ofo-share/repos/derek/ofo-catalog-data-prep"

cat("\n\n**** Starting post-curation script 01: Apply curation filters ****\n\n")
source(file.path(repo_root, "deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/01_apply-curation-filters.R"))

cat("\n\n**** Starting post-curation script 02: Summarize curated metadata ****\n\n")
source(file.path(repo_root, "deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/02_summarize-curated-metadata-per-mission.R"))

cat("\n\n**** Starting post-curation script 03: Merge metadata and curation notes ****\n\n")
source(file.path(repo_root, "deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/03_merge-metadata-and-curation-notes.R"))

cat("\n\n**** Post-curation metadata preparation complete! ****\n\n")
```

**File:** `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/control_curated-metadata_01-to-03.R`

---

## Phase 4: Modify Raw Imagery File Prep Scripts

### Step 4.1: Update `08_copy-images-to-standardized-folders.R`

Modify to:
1. Use symlinks instead of hardlinks by default
2. Support both pre-curation and post-curation paths via parameter

```r
# deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/08_copy-images-to-standardized-folders.R
# Purpose: Copy images to standardized folder structure using symlinks.
# For post-curation: uses POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH

library(tidyverse)
library(furrr)
library(sf)

#' Copy mission images to standardized folders
#'
#' @param mission_id_foc Mission ID to process
#' @param use_post_curation If TRUE, use post-curation metadata paths
copy_mission_images = function(mission_id_foc, use_post_curation = TRUE) {

  cat("\n **** Copying images for mission", mission_id_foc, "**** \n")

  # Select appropriate metadata path
  if (use_post_curation) {
    image_metadata_path = POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH
  } else {
    image_metadata_path = PARSED_EXIF_FOR_RETAINED_IMAGES_PATH
  }

  image_metadata_file = file.path(image_metadata_path, paste0(mission_id_foc, "_image-metadata.gpkg"))

  if (!file.exists(image_metadata_file)) {
    warning(paste("No image metadata found for mission", mission_id_foc))
    return(FALSE)
  }

  image_metadata = st_read(image_metadata_file, quiet = TRUE)

  # Determine absolute input and output paths
  image_metadata$image_path_contrib_abs = file.path(
    CONTRIBUTED_IMAGERY_PATH,
    image_metadata$image_path_contrib
  )

  image_metadata$image_path_ofo_abs = file.path(
    SORTED_IMAGERY_PATH,
    image_metadata$image_path_ofo
  )

  # Create output folders
  folders_out_abs = unique(dirname(image_metadata$image_path_ofo_abs))
  sapply(folders_out_abs, dir.create, recursive = TRUE, showWarnings = FALSE)

  # Create symlinks (not hardlinks)
  # Remove existing files/symlinks first
  existing = file.exists(image_metadata$image_path_ofo_abs)
  if (any(existing)) {
    file.remove(image_metadata$image_path_ofo_abs[existing])
  }

  # Create symlinks
  file.symlink(image_metadata$image_path_contrib_abs, image_metadata$image_path_ofo_abs)

  return(TRUE)
}
```

**File:** `deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/08_copy-images-to-standardized-folders.R`

---

### Step 4.2: Update `09_fix-exif.R`

Modify to replace symlinks with actual files when modifications are needed.

```r
# deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/09_fix-exif.R
# Purpose: Fix EXIF metadata. If file is a symlink and needs modification,
# replace symlink with actual file copy first.

library(furrr)
library(tidyverse)
library(exifr)

run_cmd_chunks = function(cmd, filepaths, chunk_size = 500) {
  chunks = split(filepaths, ceiling(seq_along(filepaths) / chunk_size))
  for (chunk in chunks) {
    file_string = paste0(shQuote(chunk), collapse = " ")
    system(paste0(cmd, " ", file_string))
  }
}

#' Replace symlink with actual file copy
#'
#' @param filepath Path to check and potentially replace
replace_symlink_with_copy = function(filepath) {
  if (Sys.readlink(filepath) != "") {
    # It's a symlink - get the target
    target = Sys.readlink(filepath)
    # If target is relative, make it absolute
    if (!startsWith(target, "/")) {
      target = normalizePath(file.path(dirname(filepath), target))
    }
    # Remove symlink and copy actual file
    file.remove(filepath)
    file.copy(target, filepath)
  }
}

fix_exif = function(mission_id_foc, use_post_curation = TRUE) {

  cat("\n **** Fixing EXIF for mission", mission_id_foc, "**** \n")

  folder = file.path(SORTED_IMAGERY_PATH, mission_id_foc)

  # Get all images in the folder
  image_filepaths = list.files(folder, full.names = TRUE, recursive = TRUE,
                                pattern = "(.jpg$)|(.jpeg$)|(.JPG$)|(.JPEG$)")

  if (length(image_filepaths) == 0) {
    warning(paste("No images found for mission", mission_id_foc))
    return(FALSE)
  }

  image_filenames = basename(image_filepaths)
  images = data.frame(filename = image_filenames, filepath = image_filepaths)

  # Select appropriate sorting plan path
  if (use_post_curation) {
    sorting_plan_path = POST_CURATION_IMAGE_EXIF_W_SORTING_PLAN_PATH
  } else {
    sorting_plan_path = IMAGE_EXIF_W_SORTING_PLAN_PATH
  }

  sorting_plan_file = file.path(sorting_plan_path, paste0(mission_id_foc, ".csv"))

  if (!file.exists(sorting_plan_file)) {
    warning(paste("No sorting plan found for mission", mission_id_foc))
    return(FALSE)
  }

  exif = read.csv(sorting_plan_file) |>
    select(any_of(c("image_filename_out", "Orientation", "GPSTimeStamp")))

  exif = left_join(images, exif, by = c("filename" = "image_filename_out"))

  # Determine which images need fixing
  exif$fix_orientation = ifelse(is.null(exif$Orientation), FALSE,
                                 !is.na(exif$Orientation) & exif$Orientation != 1)

  if (is.null(exif$GPSTimeStamp)) {
    exif$fix_gpstimestamp = FALSE
  } else {
    exif$fix_gpstimestamp = ifelse(is.na(exif$GPSTimeStamp), FALSE,
                                    grepl("[0-9]+:[0-9]+:[0-9]+\\.[0-9]+", exif$GPSTimeStamp))
  }

  exif = exif |>
    mutate(fix_both = (fix_orientation & fix_gpstimestamp)) |>
    mutate(fix_orientation = ifelse(fix_both, FALSE, fix_orientation),
           fix_gpstimestamp = ifelse(fix_both, FALSE, fix_gpstimestamp))

  # Get files that need any fixing
  files_to_fix = exif |> filter(fix_orientation | fix_gpstimestamp | fix_both)

  if (nrow(files_to_fix) > 0) {
    cat("Replacing symlinks with file copies for", nrow(files_to_fix), "files that need EXIF fixes...\n")
    walk(files_to_fix$filepath, replace_symlink_with_copy)
  }

  # Fix orientation only
  exif_to_fix_orientation = exif |> filter(fix_orientation == TRUE)
  if (nrow(exif_to_fix_orientation) > 0) {
    cat("Fixing orientation flag for", nrow(exif_to_fix_orientation), "images\n")
    command = "exiftool -n -r -fast4 -overwrite_original -Orientation=1"
    run_cmd_chunks(command, exif_to_fix_orientation$filepath)
  }

  # Fix GPSTimeStamp only
  exif_to_fix_gpstimestamp = exif |> filter(fix_gpstimestamp == TRUE)
  if (nrow(exif_to_fix_gpstimestamp) > 0) {
    cat("Fixing GPSTimeStamp for", nrow(exif_to_fix_gpstimestamp), "images\n")
    command = "exiftool -n -r -fast4 -overwrite_original -GPSTimeStamp="
    run_cmd_chunks(command, exif_to_fix_gpstimestamp$filepath)
  }

  # Fix both
  exif_to_fix_both = exif |> filter(fix_both == TRUE)
  if (nrow(exif_to_fix_both) > 0) {
    cat("Fixing both orientation and GPSTimeStamp for", nrow(exif_to_fix_both), "images\n")
    command = "exiftool -n -r -fast4 -overwrite_original -Orientation=1 -GPSTimeStamp="
    run_cmd_chunks(command, exif_to_fix_both$filepath)
  }

  if (nrow(exif_to_fix_orientation) == 0 && nrow(exif_to_fix_gpstimestamp) == 0 &&
      nrow(exif_to_fix_both) == 0) {
    cat("No EXIF fixes needed for", mission_id_foc, "\n")
  }

  # Check for incomplete exiftool operations
  exiftool_temp_files = list.files(folder, full.names = TRUE, recursive = TRUE,
                                    pattern = "_exiftool_tmp$")
  if (length(exiftool_temp_files) > 0) {
    warning("Exiftool temp files found in", mission_id_foc, "- incomplete operation")
    return(FALSE)
  }

  return(TRUE)
}
```

**File:** `deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/09_fix-exif.R`

---

### Step 4.3: Update `10_raw-imagery-thumbnails-and-zip.R`

Update to use post-curation metadata paths.

Add parameter `use_post_curation = TRUE` and update the metadata file paths:

```r
# In make_raw_imagery_thumbnails_and_zip function, update these lines:

# Select appropriate metadata path
if (use_post_curation) {
  points_filepath = file.path(POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
                              paste0(mission_id_foc, "_image-metadata.gpkg"))
  footprint_filepath = file.path(POST_CURATION_FULL_METADATA_PER_MISSION_PATH,
                                  paste0(mission_id_foc, "_mission-metadata.gpkg"))
} else {
  points_filepath = file.path(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
                              paste0(mission_id_foc, "_image-metadata.gpkg"))
  footprint_filepath = file.path(FULL_METADATA_PER_MISSION_PATH,
                                  paste0(mission_id_foc, "_mission-metadata.gpkg"))
}
```

**File:** `deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/10_raw-imagery-thumbnails-and-zip.R`

---

### Step 4.4: Update `11_copy-raw-imagery-to-upload-staging-dir.R`

Update to use post-curation metadata paths.

**File:** `deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/11_copy-raw-imagery-to-upload-staging-dir.R`

---

### Step 4.5: Update `08-13_prep-raw-imagery-files_per-mission.R`

Add `use_post_curation` parameter that propagates to all sub-functions:

```r
prep_raw_imagery_files_per_mission = function(mission_id_foc, use_post_curation = TRUE) {

  copy_mission_images(mission_id_foc, use_post_curation = use_post_curation)
  exif_success = fix_exif(mission_id_foc, use_post_curation = use_post_curation)

  if (!exif_success) {
    warning("EXIF fix failed for mission", mission_id_foc)
    return(FALSE)
  }

  make_raw_imagery_thumbnails_and_zip(mission_id_foc, use_post_curation = use_post_curation)
  copy_raw_imagery_to_upload_staging_dir(mission_id_foc, use_post_curation = use_post_curation)
  upload_raw_imagery_to_object_store(mission_id_foc)
  delete_prepped_raw_imagery(mission_id_foc)

  return(TRUE)
}
```

**File:** `deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/08-13_prep-raw-imagery-files_per-mission.R`

---

## Phase 5: Update Web Catalog Scripts

### Step 5.1: Update `10_create-drone-data-catalog-webpages.R`

Ensure it can use post-curation metadata. The script already reads from `MISSION_METADATA_FILEPATH` and `IMAGE_METADATA_FILEPATH`, so the approach is:

1. After running post-curation pipeline, combine post-curation metadata files into the combined files
2. OR update the constants to point to post-curation paths

**Recommendation:** Create a script to combine post-curation metadata files:

```r
# deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/04_combine-metadata-files.R
# Purpose: Combine individual mission metadata files into combined files for web catalog

library(tidyverse)
library(sf)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")

# Combine mission metadata
mission_files = list.files(POST_CURATION_FULL_METADATA_PER_MISSION_PATH,
                           pattern = "\\.gpkg$", full.names = TRUE)

cat(sprintf("Combining %d mission metadata files...\n", length(mission_files)))

mission_metadata_list = map(mission_files, ~ st_read(.x, quiet = TRUE))
mission_metadata_combined = bind_rows(mission_metadata_list)

st_write(mission_metadata_combined, POST_CURATION_FULL_METADATA_PER_MISSION_COMBINED_FILEPATH,
         delete_dsn = TRUE)
cat(sprintf("Written to: %s\n", POST_CURATION_FULL_METADATA_PER_MISSION_COMBINED_FILEPATH))

# Combine image metadata
image_files = list.files(POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
                         pattern = "\\.gpkg$", full.names = TRUE)

cat(sprintf("Combining %d image metadata files...\n", length(image_files)))

image_metadata_list = map(image_files, ~ st_read(.x, quiet = TRUE))
image_metadata_combined = bind_rows(image_metadata_list)

st_write(image_metadata_combined, POST_CURATION_FULL_METADATA_PER_IMAGE_COMBINED_FILEPATH,
         delete_dsn = TRUE)
cat(sprintf("Written to: %s\n", POST_CURATION_FULL_METADATA_PER_IMAGE_COMBINED_FILEPATH))

cat("Metadata combination complete.\n")
```

**File:** `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/04_combine-metadata-files.R`

---

## Phase 6: Testing and Validation

### Step 6.1: Create validation script

```r
# deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/validate_curation_pipeline.R
# Purpose: Validate that post-curation pipeline produced expected outputs

library(tidyverse)
library(sf)

source("deploy/drone-imagery-ingestion/00_set-constants.R")

cat("=== Validating Post-Curation Pipeline Outputs ===\n\n")

# Check 1: Filtered image metadata exists
cat("1. Checking filtered image metadata...\n")
filtered_files = list.files(POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH, pattern = "\\.gpkg$")
cat(sprintf("   Found %d filtered image metadata files\n", length(filtered_files)))

# Check 2: Derived metadata exists
cat("2. Checking derived metadata...\n")
derived_mission_files = list.files(POST_CURATION_DERIVED_METADATA_PER_MISSION_PATH, pattern = "\\.gpkg$")
derived_sub_files = list.files(POST_CURATION_DERIVED_METADATA_PER_SUB_MISSION_PATH, pattern = "\\.gpkg$")
cat(sprintf("   Found %d mission and %d sub-mission derived metadata files\n",
            length(derived_mission_files), length(derived_sub_files)))

# Check 3: Full metadata with curation notes exists
cat("3. Checking full metadata with curation notes...\n")
full_mission_files = list.files(POST_CURATION_FULL_METADATA_PER_MISSION_PATH, pattern = "\\.gpkg$")
full_sub_files = list.files(POST_CURATION_FULL_METADATA_PER_SUB_MISSION_PATH, pattern = "\\.gpkg$")
cat(sprintf("   Found %d mission and %d sub-mission full metadata files\n",
            length(full_mission_files), length(full_sub_files)))

# Check 4: Curation columns present
cat("4. Checking curation columns in metadata...\n")
sample_file = list.files(POST_CURATION_FULL_METADATA_PER_MISSION_PATH, pattern = "\\.gpkg$", full.names = TRUE)[1]
if (length(sample_file) > 0) {
  sample_metadata = st_read(sample_file, quiet = TRUE)
  expected_cols = c("anomaly_severity", "anomalies_altitude", "anomalies_spatial", "anomaly_notes")
  present_cols = expected_cols[expected_cols %in% names(sample_metadata)]
  cat(sprintf("   Curation columns present: %s\n", paste(present_cols, collapse = ", ")))
}

# Check 5: Image counts reduced appropriately
cat("5. Checking image count reduction...\n")
curation_notes = read_csv(CURATION_NOTES_FILEPATH, col_types = cols(.default = "c"))
missions_with_extraneous = curation_notes |>
  filter(!is.na(extraneous_images) & extraneous_images != "") |>
  pull(mission_id) |>
  unique()

for (mission_id in head(missions_with_extraneous, 3)) {
  pre_file = file.path(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH, paste0(mission_id, "_image-metadata.gpkg"))
  post_file = file.path(POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH, paste0(mission_id, "_image-metadata.gpkg"))

  if (file.exists(pre_file) && file.exists(post_file)) {
    pre_count = nrow(st_read(pre_file, quiet = TRUE))
    post_count = nrow(st_read(post_file, quiet = TRUE))
    cat(sprintf("   Mission %s: %d -> %d images (-%d)\n",
                mission_id, pre_count, post_count, pre_count - post_count))
  }
}

cat("\n=== Validation Complete ===\n")
```

**File:** `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/validate_curation_pipeline.R`

---

## Phase 7: Folder Renumbering (Post-Implementation)

After the implementation is tested and working:

### Step 7.1: Rename folders

```
01_raw-imagery-metadata-prep      → 01_pre-curation-metadata-prep
01b_curated-raw-imagery-metadata-prep → 02_post-curation-metadata-prep
02_raw-imagery-file-prep          → 03_imagery-file-prep
03_photogrammetry                 → 04_photogrammetry
04_itd                            → 05_itd
10_drone-mission-web-catalog      → 10_drone-mission-web-catalog (unchanged)
```

### Step 7.2: Update all `source()` calls

Search and replace in all R files:
- `01_raw-imagery-metadata-prep` → `01_pre-curation-metadata-prep`
- `01b_curated-raw-imagery-metadata-prep` → `02_post-curation-metadata-prep`
- `02_raw-imagery-file-prep` → `03_imagery-file-prep`
- etc.

### Step 7.3: Update constants file paths if necessary

---

## Summary: Files to Create/Modify

### New Files:
1. `src/curation-utils.R` - Utility functions for parsing curation notes
2. `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/01_apply-curation-filters.R`
3. `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/02_summarize-curated-metadata-per-mission.R`
4. `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/03_merge-metadata-and-curation-notes.R`
5. `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/04_combine-metadata-files.R`
6. `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/control_curated-metadata_01-to-03.R`
7. `deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/validate_curation_pipeline.R`

### Modified Files:
1. `deploy/drone-imagery-ingestion/00_set-constants.R` - Add new path constants
2. `deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/08_copy-images-to-standardized-folders.R` - Use symlinks, add post-curation support
3. `deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/09_fix-exif.R` - Handle symlink replacement
4. `deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/10_raw-imagery-thumbnails-and-zip.R` - Add post-curation support
5. `deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/11_copy-raw-imagery-to-upload-staging-dir.R` - Add post-curation support
6. `deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/08-13_prep-raw-imagery-files_per-mission.R` - Add post-curation parameter

---

## Execution Order

1. **Phase 1:** Update constants (Step 1.1)
2. **Phase 2:** Create utility functions (Step 2.1)
3. **Phase 3:** Create post-curation metadata scripts (Steps 3.1-3.4)
4. **Phase 4:** Modify file prep scripts (Steps 4.1-4.5)
5. **Phase 5:** Update/create web catalog support (Step 5.1)
6. **Phase 6:** Test and validate (Step 6.1)
7. **Phase 7:** Folder renumbering (after validation passes)

---

## Open Questions Resolved

| Question | Resolution |
|----------|------------|
| Image ID matching for formerly curated | Spatial join handles cases where IDs match but locations differ, or IDs added/dropped |
| Duplicate rows in curation notes | Combine additively |
| Cross-mission references | Error with descriptive message (typo) |
| SORTED_IMAGERY_PATH | Added to constants: `/z_OLD_ofo-share/.../2_sorted` |
| Symlinks vs hardlinks | Use symlinks; replace with file copy when modification needed |
| Script 06 still needed | Yes, reused in post-curation summarization |
| Folder renumbering | After implementation tested and working |
| Post-curation final paths | Yes, create `6_post-curation-final` structure |
