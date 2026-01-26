# src/curation-utils.R
# Utility functions for parsing and applying curation notes.
# Used by post-curation scripts to filter extraneous images and update image IDs.

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
#' @param expected_mission_id The mission ID these images should belong to (for validation).
#'   If provided, will error on cross-mission references (likely typos in curation notes).
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
    # Skip empty parts
    if (part == "") next

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
      # Format: NNNNNN-NN_NNNNNN (e.g., 000034-01_000912)
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

      # Handle reverse ranges (e.g., 3332 to 3255)
      if (start_num > end_num) {
        num_seq = seq(end_num, start_num)
      } else {
        num_seq = seq(start_num, end_num)
      }

      # Generate all IDs in range with consistent 6-digit padding
      range_ids = paste0(start_prefix, "_", str_pad(num_seq, 6, pad = "0"))
      all_image_ids = c(all_image_ids, range_ids)

    } else {
      # Single image ID
      single_id = trimws(part)

      # Validate format - must match NNNNNN-NN_NNNNNN pattern
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
#' Curation notes may have multiple rows for the same sub-mission (e.g., when
#' different types of anomalies are recorded separately). This function combines
#' them into one row per sub-mission.
#'
#' @param curation_filepath Path to curation notes CSV
#' @return Tibble with one row per sub-mission, combined extraneous_images and anomaly notes
load_and_combine_curation_notes = function(curation_filepath) {
  curation_raw = read_csv(curation_filepath, col_types = cols(.default = "c"))

  # Define anomaly columns to aggregate
  anomaly_cols = c("anomaly_severity", "collection_time_anomaly", "altitude_anomaly",
                   "spatial_anomaly", "camera_pitch_anomaly", "excess_images_anomaly",
                   "missing_images_anomaly", "other_anomalies", "anomaly_notes")

  # Filter to only columns that exist in the data
  anomaly_cols = intersect(anomaly_cols, names(curation_raw))

  # Group by mission_id and sub_mission_id, combine values
  curation_combined = curation_raw |>
    group_by(mission_id, sub_mission_id) |>
    summarise(
      # Combine extraneous_images with comma separator
      extraneous_images = paste(na.omit(extraneous_images[extraneous_images != ""]), collapse = ", "),
      # For anomaly columns, combine non-empty values with comma
      across(all_of(anomaly_cols), ~ paste(na.omit(.x[.x != ""]), collapse = ", ")),
      .groups = "drop"
    ) |>
    # Replace empty strings with NA for cleaner handling
    mutate(across(everything(), ~ if_else(.x == "", NA_character_, .x)))

  return(curation_combined)
}


#' Map image IDs from formerly curated metadata to current metadata using coordinate and attribute matching
#'
#' For missions that were curated using an older version of the image metadata,
#' this function creates a crosswalk between the former image IDs and the current
#' image IDs by matching images based on their geographic coordinates and additional
#' attributes (file size, camera settings, etc.) when available.
#'
#' @param formerly_curated_images SF object with formerly curated image points.
#'   Must have columns: mission_id, image_id, and geometry. Optional columns used
#'   for matching if present and non-null: file_size_gb, accuracy_x, accuracy_y,
#'   camera_pitch, camera_roll, camera_yaw, exposure, aperture, iso.
#' @param current_images SF object with current image points.
#'   Must have columns: mission_id, image_id, and geometry. Same optional columns.
#' @param mission_id Mission ID to filter on
#' @param coord_precision Number of decimal places for coordinate matching (default 8).
#'   Uses truncation to avoid rounding ambiguity.
#' @return Tibble with columns: former_image_id, current_image_id.
#'   current_image_id will be:
#'   - The actual ID if exactly one match at that location (regardless of sub-mission)
#'   - "none" if no matches at that location
#'   - "multiple" if more than one match at that location
create_image_id_crosswalk = function(formerly_curated_images, current_images,
                                     mission_id, coord_precision = 8) {
  # Filter to focal mission
  former = formerly_curated_images |> filter(mission_id == !!mission_id)
  current = current_images |> filter(mission_id == !!mission_id)

  if (nrow(former) == 0 || nrow(current) == 0) {
    warning(paste("No images found for mission", mission_id, "in one or both datasets"))
    return(tibble(former_image_id = character(), current_image_id = character()))
  }

  # Helper to extract sub-mission ID from image ID (format: NNNNNN-NN_NNNNNN)
  extract_sub_mission = function(image_id) substr(image_id, 1, 9)

  # Extract coordinates (always used)
  former_coords = st_coordinates(former)
  current_coords = st_coordinates(current)

  # Use truncation (floor of scaled value) to avoid rounding ambiguity
  coord_scale = 10^coord_precision

  # Start with coordinate-based key components
  former_key_parts = list(
    floor(former_coords[, 1] * coord_scale),
    floor(former_coords[, 2] * coord_scale)
  )
  current_key_parts = list(
    floor(current_coords[, 1] * coord_scale),
    floor(current_coords[, 2] * coord_scale)
  )

  # Define additional matching attributes and their precisions
  match_attrs = list(
    file_size_gb = 9,
    accuracy_x = 4,
    accuracy_y = 4,
    camera_pitch = 2,
    camera_roll = 2,
    camera_yaw = 2,
    exposure = 6,
    aperture = 2,
    iso = 0
  )

  # Add each attribute to the key if it exists and has no NULLs in both datasets
  for (attr in names(match_attrs)) {
    precision = match_attrs[[attr]]
    scale = 10^precision

    # Check if column exists in both and has no NAs
    former_has = attr %in% names(former) && !any(is.na(former[[attr]]))
    current_has = attr %in% names(current) && !any(is.na(current[[attr]]))

    if (former_has && current_has) {
      former_key_parts = c(former_key_parts, list(floor(former[[attr]] * scale)))
      current_key_parts = c(current_key_parts, list(floor(current[[attr]] * scale)))
    }
  }

  # Build final keys by pasting all components
  former_key = do.call(paste, c(former_key_parts, sep = "_"))
  current_key = do.call(paste, c(current_key_parts, sep = "_"))

  # Build lookup table: key -> list of current image IDs at that location
  current_lookup = split(current$image_id, current_key)

  # Match each former image
  matches = current_lookup[former_key]
  match_counts = lengths(matches)

  # Initialize result
  crosswalk = tibble(
    former_image_id = former$image_id,
    current_image_id = character(length(former_key))
  )

  # Handle zero matches
  no_match = match_counts == 0
  crosswalk$current_image_id[no_match] = "none"

  # Handle multiple matches
  multi_match = match_counts > 1
  crosswalk$current_image_id[multi_match] = "multiple"

  # Handle single matches - always return actual image ID
  single_match = match_counts == 1
  if (any(single_match)) {
    single_idx = which(single_match)
    single_current_ids = map_chr(matches[single_idx], ~ .x[1])
    crosswalk$current_image_id[single_idx] = single_current_ids
  }

  return(crosswalk)
}


#' Update image ID references in a string using crosswalk
#'
#' Finds all image ID patterns in a text string and replaces them with their
#' corresponding current image IDs from the crosswalk.
#'
#' @param text String potentially containing image ID references
#' @param crosswalk Tibble with former_image_id and current_image_id columns
#' @return String with updated image IDs
update_image_id_references = function(text, crosswalk) {
  if (is.na(text) || text == "") {
    return(text)
  }

  # Find all image ID patterns in the text (format: NNNNNN-NN_NNNNNN)
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
    } else {
      warning(paste("Former image ID", old_id, "could not be mapped to a current image ID"))
    }
  }

  return(result)
}


#' Find the remaining duplicate for an image that was removed during deduplication
#'
#' Given an image ID that may have been removed as a duplicate, look up its
#' duplicate group and return the image ID that was kept (alphabetically first
#' by image_path_contrib).
#'
#' @param image_id The image ID to look up
#' @param duplicate_log Tibble containing the catalog-wide duplicate log with columns:
#'   duplicate_group_id, image_id, image_path_contrib
#' @return Character: the kept image ID if found in a duplicate group,
#'   NA if not a duplicate, or the original image_id if it was the one kept
get_remaining_duplicate = function(image_id, duplicate_log) {
  # Find the duplicate group containing this image
  image_record = duplicate_log |>
    filter(image_id == !!image_id)

  if (nrow(image_record) == 0) {
    return(NA_character_)
  }

  group = image_record$duplicate_group_id[1]

  # Get all images in this group and find the one that was kept
  # (alphabetically first by image_path_contrib)
  group_images = duplicate_log |>
    filter(duplicate_group_id == group) |>
    arrange(image_path_contrib)

  # The first one was kept
  kept_id = group_images$image_id[1]
  return(kept_id)
}


#' Update image ID references in anomaly columns
#'
#' Updates image IDs in anomaly columns using the crosswalk (if provided) and
#' duplicate resolution. This function is for anomaly columns only (not extraneous_images).
#' Anomaly columns contain individual image IDs, not ranges.
#'
#' For formerly curated missions: applies crosswalk AND duplicate resolution.
#' For currently curated missions: applies only duplicate resolution (crosswalk = NULL).
#'
#' @param text String containing image ID references (individual IDs, not ranges)
#' @param crosswalk Tibble with former_image_id and current_image_id columns,
#'   or NULL for currently curated missions (no crosswalking needed)
#' @param duplicate_log Tibble with catalog-wide duplicate information
#' @return String with updated image IDs. IDs that no longer exist are marked as "(image deleted)".
update_anomaly_image_references = function(text, crosswalk, duplicate_log) {
  if (is.na(text) || text == "") {
    return(text)
  }

  # Find all image ID patterns in the text (format: NNNNNN-NN_NNNNNN)
  pattern = "\\d{6}-\\d{2}_\\d{6}"
  matches = str_extract_all(text, pattern)[[1]]

  if (length(matches) == 0) {
    return(text)
  }

  has_crosswalk = !is.null(crosswalk) && nrow(crosswalk) > 0
  result = text

  for (old_id in matches) {
    new_id = old_id  # Default: keep original

    # Step 1: Apply crosswalk if available (formerly curated missions)
    if (has_crosswalk) {
      crosswalk_result = crosswalk |>
        filter(former_image_id == old_id) |>
        pull(current_image_id)

      if (length(crosswalk_result) > 0) {
        new_id = crosswalk_result[1]
      }
    }

    # Step 2: Handle crosswalk results and check for duplicates
    if (new_id == "none") {
      # Crosswalk says image doesn't exist in current data.
      # No duplicate lookup needed: if duplicates existed at this location,
      # one would remain after deduplication and the crosswalk would have found it.
      result = str_replace(result, fixed(old_id), "(image deleted)")
    } else if (new_id == "multiple") {
      # Multiple matches - ambiguous, keep original with warning
      warning(paste("Multiple matches found for image", old_id, "- keeping original"))
    } else {
      # new_id is a valid ID (either crosswalked or original)
      # Check if this ID was a removed duplicate and resolve to remaining
      remaining = get_remaining_duplicate(new_id, duplicate_log)

      if (!is.na(remaining) && remaining != new_id) {
        # The ID we have was a removed duplicate; use the remaining one
        result = str_replace(result, fixed(old_id), remaining)
      } else if (new_id != old_id) {
        # Crosswalked to a different ID
        result = str_replace(result, fixed(old_id), new_id)
      }
      # If new_id == old_id and not a removed duplicate, no change needed
    }
  }

  return(result)
}
