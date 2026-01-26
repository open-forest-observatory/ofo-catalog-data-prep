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
#'   - The actual ID if exactly one match at that location and same sub-mission
#'   - "none" if no matches at that location
#'   - "multiple" if more than one match at that location
#'   - "different sub-mission" if exactly one match but in a different sub-mission
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

  # Handle single matches - check sub-mission
  single_match = match_counts == 1
  if (any(single_match)) {
    single_idx = which(single_match)
    single_current_ids = map_chr(matches[single_idx], ~ .x[1])

    single_former_subs = extract_sub_mission(crosswalk$former_image_id[single_idx])
    single_current_subs = extract_sub_mission(single_current_ids)

    same_sub = single_former_subs == single_current_subs

    crosswalk$current_image_id[single_idx] = if_else(
      same_sub,
      single_current_ids,
      "different sub-mission"
    )
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
