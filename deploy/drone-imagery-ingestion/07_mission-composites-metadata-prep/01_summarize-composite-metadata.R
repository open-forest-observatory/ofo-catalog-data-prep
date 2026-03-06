# Purpose: Enrich composite mission metadata by:
# 1. Adding non-image-dependent attributes from each constituent mission's individual metadata on S3
# 2. Computing image-dependent derived attributes by summarizing composite image data per-mission
# 3. Computing photogrammetry-derived altitude metrics from altitude_agl
#
# The enriched metadata is re-uploaded to S3, overwriting the original.

library(tidyverse)
library(sf)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/summarization-utils.R")  # brings in metadata-extraction_imagery_dataset.R and utils.R

# ==============================================================================
# Constants
# ==============================================================================

GEOMETRY_COLUMNS = c("geom", "geometry", "x")

# ==============================================================================
# Helper Functions
# ==============================================================================

#' List composite mission folders on S3
#'
#' @param remote Name of the rclone remote
#' @param remote_dir Remote directory path for composites
#' @return Character vector of composite IDs (e.g., "000001_000002")
list_composites_on_s3 = function(remote, remote_dir) {
  cmd = paste0("rclone lsf ", remote, ":", remote_dir, " --dirs-only --config /dev/null")
  output = system(cmd, intern = TRUE)
  # Remove trailing slashes
  composite_ids = str_remove(output, "/$")
  return(composite_ids)
}


#' Download a single file from S3 via rclone
#'
#' @param remote_path Full remote path (remote:path)
#' @param local_path Local destination path
#' @return TRUE on success, FALSE on failure
download_s3_file = function(remote_path, local_path) {
  create_dir(dirname(local_path))
  cmd = paste0("rclone copyto ", remote_path, " ", local_path, " --config /dev/null")
  result = system(cmd)
  return(result == 0)
}


#' Fix missing polygon rows in composite mission metadata
#'
#' If the composite mission metadata has only 1 row (only 'hn'), downloads the
#' missing mission's individual metadata to get its polygon, and creates a 2-row
#' sf object.
#'
#' @param composite_mission_meta SF object with composite mission metadata
#' @param composite_id The composite ID (e.g., "000001_000002")
#' @param composite_image_meta Data frame with composite image metadata
#' @return 2-row sf object
fix_missing_polygon = function(composite_mission_meta, composite_id, composite_image_meta) {
  if (nrow(composite_mission_meta) == 2) {
    return(composite_mission_meta)
  }

  # Parse composite_id to get both mission IDs
  mission_ids = str_split(composite_id, "_")[[1]]

  # Identify which mission ID is missing
  existing_mission_id = composite_mission_meta$mission_id
  missing_mission_id = setdiff(mission_ids, existing_mission_id)

  if (length(missing_mission_id) == 0) {
    warning("Could not identify missing mission ID for composite ", composite_id)
    return(composite_mission_meta)
  }

  # Download the missing mission's individual metadata from S3
  remote_meta_path = paste0(
    RCLONE_REMOTE, ":", REMOTE_MISSIONS_DIR,
    missing_mission_id, "/metadata-mission/", missing_mission_id, "_mission-metadata.gpkg"
  )
  local_meta_path = file.path(COMPOSITE_METADATA_TEMPDIR, "individual-missions",
                               paste0(missing_mission_id, "_mission-metadata.gpkg"))

  success = download_s3_file(remote_meta_path, local_meta_path)
  if (!success) {
    warning("Failed to download individual metadata for mission ", missing_mission_id)
    return(composite_mission_meta)
  }

  individual_meta = st_read(local_meta_path, quiet = TRUE)

  # Determine mission_type from image metadata for the missing mission
  missing_images = composite_image_meta |> filter(mission_id == missing_mission_id)
  if (nrow(missing_images) > 0 && "mission_type" %in% names(missing_images)) {
    missing_type = unique(missing_images$mission_type)[1]
  } else {
    # Infer as complement: if existing is "hn", the other could be "ln" or vice versa
    existing_type = composite_mission_meta$mission_type[1]
    missing_type = ifelse(existing_type == "hn", "ln", "hn")
    warning("Inferred mission_type '", missing_type, "' for missing mission ", missing_mission_id)
  }

  # Get the active geometry column name to ensure it matches when binding rows
  geom_col = attr(composite_mission_meta, "sf_column")

  # Ensure CRS matches before extracting geometry
  individual_meta = st_transform(individual_meta, st_crs(composite_mission_meta))

  # Create new row using the same geometry column name as composite_mission_meta
  new_row_data = data.frame(
    composite_id = composite_id,
    mission_type = missing_type,
    mission_id = missing_mission_id,
    date = composite_mission_meta$date[1],
    area_m2 = as.numeric(st_area(individual_meta[1, ]))
  )
  new_row_data[[geom_col]] = st_geometry(individual_meta)[1]
  new_row = st_sf(new_row_data, sf_column_name = geom_col)

  # Combine
  result = bind_rows(composite_mission_meta, new_row)
  return(result)
}


#' Extract contributed (non-derived) columns from individual mission metadata
#'
#' Follows the pattern from 03_merge-metadata-and-curation-notes.R.
#' Extracts all columns NOT ending in _derived, PLUS keeps camera_pitch_derived
#' and smart_oblique_derived. Drops geometry.
#'
#' @param individual_mission_meta SF object with full individual mission metadata
#' @return Tibble with contributed columns (plus camera_pitch_derived and smart_oblique_derived)
extract_contributed_columns = function(individual_mission_meta) {
  all_cols = names(individual_mission_meta)

  # Identify derived columns (those ending in _derived)
  derived_cols = all_cols[str_detect(all_cols, "_derived$")]

  # Keep camera_pitch_derived and smart_oblique_derived (mission-setup properties)
  keep_derived = c("camera_pitch_derived", "smart_oblique_derived")
  derived_to_drop = setdiff(derived_cols, keep_derived)

  # Identify geometry columns
  geom_cols = intersect(all_cols, GEOMETRY_COLUMNS)

  # Contributed columns = all columns minus derived-to-drop and geometry
  cols_to_keep = setdiff(all_cols, c(derived_to_drop, geom_cols))

  individual_mission_meta |>
    st_drop_geometry() |>
    select(all_of(cols_to_keep))
}


#' Force all non-geometry columns to character (for safe binding)
force_all_cols_to_character = function(df) {
  df |>
    mutate(across(-any_of(c("geometry", "geom")), as.character))
}


# ==============================================================================
# Main Processing Function
# ==============================================================================

#' Process a single composite mission
#'
#' @param composite_id_foc The composite ID to process (e.g., "000001_000002")
#' @return TRUE on success, FALSE on failure
process_composite = function(composite_id_foc) {
  tryCatch({
    cat(sprintf("Processing composite %s...\n", composite_id_foc))

    # --- 1. Download from S3 ---
    composite_dir = file.path(COMPOSITE_METADATA_TEMPDIR, composite_id_foc)
    create_dir(composite_dir)

    remote_mission_meta_path = paste0(
      RCLONE_REMOTE, ":", REMOTE_COMPOSITES_DIR,
      composite_id_foc, "/metadata-mission/", composite_id_foc, "_mission-metadata.gpkg"
    )
    remote_image_meta_path = paste0(
      RCLONE_REMOTE, ":", REMOTE_COMPOSITES_DIR,
      composite_id_foc, "/metadata-images/", composite_id_foc, "_image-metadata.gpkg"
    )

    local_mission_meta_path = file.path(composite_dir, paste0(composite_id_foc, "_mission-metadata.gpkg"))
    local_image_meta_path = file.path(composite_dir, paste0(composite_id_foc, "_image-metadata.gpkg"))

    if (!download_s3_file(remote_mission_meta_path, local_mission_meta_path)) {
      warning("Failed to download mission metadata for composite ", composite_id_foc)
      return(FALSE)
    }
    if (!download_s3_file(remote_image_meta_path, local_image_meta_path)) {
      warning("Failed to download image metadata for composite ", composite_id_foc)
      return(FALSE)
    }

    # --- 2. Read files ---
    composite_mission_meta = st_read(local_mission_meta_path, quiet = TRUE)
    composite_image_meta = st_read(local_image_meta_path, quiet = TRUE)

    # --- 3. Fix missing polygon ---
    composite_mission_meta = fix_missing_polygon(composite_mission_meta, composite_id_foc, composite_image_meta)

    if (nrow(composite_mission_meta) != 2) {
      warning("Expected 2 rows in composite mission metadata for ", composite_id_foc,
              " but got ", nrow(composite_mission_meta))
      return(FALSE)
    }

    # --- 4. Prepare images: extract lon/lat from geometry, drop geometry ---
    coords = st_coordinates(composite_image_meta)
    composite_images_df = composite_image_meta |>
      st_drop_geometry() |>
      mutate(lon = coords[, 1], lat = coords[, 2])

    # --- 5. Process each of the 2 mission rows ---
    enriched_rows = list()

    for (i in 1:2) {
      mission_row = composite_mission_meta[i, ]
      mission_id_foc = mission_row$mission_id
      mission_polygon = st_geometry(mission_row)

      cat(sprintf("  Processing mission %s within composite %s...\n", mission_id_foc, composite_id_foc))

      # Download individual mission metadata from S3
      remote_indiv_meta_path = paste0(
        RCLONE_REMOTE, ":", REMOTE_MISSIONS_DIR,
        mission_id_foc, "/metadata-mission/", mission_id_foc, "_mission-metadata.gpkg"
      )
      local_indiv_meta_path = file.path(COMPOSITE_METADATA_TEMPDIR, "individual-missions",
                                         paste0(mission_id_foc, "_mission-metadata.gpkg"))

      if (!download_s3_file(remote_indiv_meta_path, local_indiv_meta_path)) {
        warning("Failed to download individual metadata for mission ", mission_id_foc)
        return(FALSE)
      }

      individual_meta = st_read(local_indiv_meta_path, quiet = TRUE)

      # Extract contributed columns (non-derived + camera_pitch_derived + smart_oblique_derived)
      contributed = extract_contributed_columns(individual_meta)

      # Filter composite images to this mission
      filtered_images = composite_images_df |> filter(mission_id == mission_id_foc)

      if (nrow(filtered_images) < 10) {
        warning("Fewer than 10 images for mission ", mission_id_foc,
                " in composite ", composite_id_foc, "; skipping")
        return(FALSE)
      }

      # Call extract_imagery_dataset_metadata() for standard derived attributes
      # and photogrammetry altitude metrics (via include_photogrammetry_altitude parameter)
      derived = extract_imagery_dataset_metadata(
        metadata = filtered_images |> as.data.frame(),
        mission_polygon = st_as_sf(mission_polygon),
        dataset_id = mission_id_foc,
        use_new_fidelity_col = TRUE,
        include_photogrammetry_altitude = TRUE
      )

      if (is.null(derived)) {
        warning("extract_imagery_dataset_metadata returned NULL for mission ", mission_id_foc)
        return(FALSE)
      }

      # Override camera_pitch_derived and smart_oblique_derived with values from
      # individual mission metadata (these are mission-setup properties)
      if ("camera_pitch_derived" %in% names(contributed)) {
        derived$camera_pitch_derived = contributed$camera_pitch_derived
      }
      if ("smart_oblique_derived" %in% names(contributed)) {
        derived$smart_oblique_derived = contributed$smart_oblique_derived
      }

      # Remove camera_pitch_derived and smart_oblique_derived from contributed
      # since they are now in derived (avoid duplication)
      contributed = contributed |>
        select(-any_of(c("camera_pitch_derived", "smart_oblique_derived")))

      # Also remove mission_id from contributed if present (it's already in composite meta)
      contributed = contributed |>
        select(-any_of(c("mission_id")))

      # Assemble enriched row: composite columns + contributed + derived (includes photogrammetry) + geometry
      composite_cols = mission_row |>
        st_drop_geometry() |>
        select(composite_id, mission_type, mission_id)

      # Force all to character before binding to prevent type mismatches
      composite_cols_chr = force_all_cols_to_character(composite_cols)
      contributed_chr = force_all_cols_to_character(contributed)
      derived_chr = force_all_cols_to_character(derived)

      enriched = bind_cols(composite_cols_chr, contributed_chr, derived_chr)

      # Re-attach geometry from the composite mission polygon
      enriched_sf = st_sf(enriched, geometry = st_geometry(mission_row))

      enriched_rows[[i]] = enriched_sf
    }

    # --- 6. Combine the 2 enriched rows ---
    enriched_combined = bind_rows(enriched_rows[[1]], enriched_rows[[2]])

    # --- 7. Write enriched gpkg locally, then upload to S3 ---
    output_path = file.path(composite_dir, paste0(composite_id_foc, "_mission-metadata.gpkg"))
    st_write(enriched_combined, output_path, delete_dsn = TRUE, quiet = TRUE)

    # Upload to S3, overwriting the original
    remote_upload_dir = paste0(
      RCLONE_REMOTE, ":", REMOTE_COMPOSITES_DIR,
      composite_id_foc, "/metadata-mission/"
    )
    upload_cmd = paste0("rclone copy ", output_path, " ", remote_upload_dir, " --config /dev/null")
    upload_result = system(upload_cmd)

    if (upload_result != 0) {
      warning("Failed to upload enriched metadata for composite ", composite_id_foc)
      return(FALSE)
    }

    # --- 8. Clean up temp files ---
    unlink(composite_dir, recursive = TRUE)

    cat(sprintf("  Successfully processed composite %s\n", composite_id_foc))
    return(TRUE)

  }, error = function(e) {
    warning("Error processing composite ", composite_id_foc, ": ", e$message)
    return(FALSE)
  })
}


# ==============================================================================
# Main Script Body
# ==============================================================================

composite_ids = list_composites_on_s3(RCLONE_REMOTE, REMOTE_COMPOSITES_DIR)
cat(sprintf("Found %d composites on S3\n", length(composite_ids)))

create_dir(COMPOSITE_METADATA_TEMPDIR)

# Set up for parallelizing across all composites
future::plan(future::multisession, workers = future::availableCores()*3)

results = furrr::future_map_lgl(composite_ids, process_composite, .progress = TRUE)

cat(sprintf("\nProcessed %d of %d composites successfully\n", sum(results), length(results)))
