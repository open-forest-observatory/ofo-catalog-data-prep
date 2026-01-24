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


#' Map image IDs from formerly curated metadata to current metadata using spatial join
#'
#' For missions that were curated using an older version of the image metadata,
#' this function creates a crosswalk between the former image IDs and the current
#' image IDs by matching images based on their geographic coordinates.
#'
#' @param formerly_curated_images SF object with formerly curated image points.
#'   Must have columns: mission_id, image_id, and geometry.
#' @param current_images SF object with current image points.
#'   Must have columns: mission_id, image_id, and geometry.
#' @param mission_id Mission ID to filter on
#' @param distance_threshold Maximum distance (meters) for spatial matching.
#'   Matches beyond this distance will generate a warning.
#' @return Tibble with columns: former_image_id, current_image_id, distance_m
create_image_id_crosswalk = function(formerly_curated_images, current_images,
                                      mission_id, distance_threshold = 1) {
  # Filter to focal mission
  former = formerly_curated_images |> filter(mission_id == !!mission_id)
  current = current_images |> filter(mission_id == !!mission_id)

  if (nrow(former) == 0 || nrow(current) == 0) {
    warning(paste("No images found for mission", mission_id, "in one or both datasets"))
    return(tibble(former_image_id = character(), current_image_id = character(), distance_m = numeric()))
  }

  # Get centroid coordinates to determine appropriate UTM zone
  centroid_coords = st_coordinates(st_centroid(st_union(former)))

  # Calculate UTM zone from longitude
  utm_zone = floor((centroid_coords[1, 1] + 180) / 6) + 1

  # Determine if northern or southern hemisphere
  if (centroid_coords[1, 2] >= 0) {
    utm_epsg = 32600 + utm_zone  # Northern hemisphere
  } else {
    utm_epsg = 32700 + utm_zone  # Southern hemisphere
  }

  # Transform to UTM for accurate distance calculation
  former_utm = former |> st_transform(utm_epsg)
  current_utm = current |> st_transform(utm_epsg)

  # Find nearest current image for each former image
  nearest_idx = st_nearest_feature(former_utm, current_utm)
  distances = st_distance(former_utm, current_utm[nearest_idx, ], by_element = TRUE)

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
