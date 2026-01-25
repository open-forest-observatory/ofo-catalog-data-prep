# src/summarization-utils.R
# Shared functions for summarizing image metadata at mission/sub-mission level.
# Used by both pre-curation (script 06) and post-curation (script 02) workflows.

library(tidyverse)
library(sf)

source("src/metadata-extraction_imagery_dataset.R")
source("src/utils.R")


#' Compute mission/sub-mission polygons and identify retained images
#'
#' Splits image metadata by the specified column, computes polygons for each group,
#' and identifies which images fall within the computed polygons.
#'
#' @param image_metadata Data frame with image metadata including lon, lat columns
#' @param column_to_split_on Column name to group by ("mission_id" or "sub_mission_id")
#' @param image_merge_distance Distance threshold for polygon computation
#' @return List with `polygons` (named list of sf objects) and `retained_image_IDs` (character vector)
compute_polygons_and_images_retained = function(image_metadata, column_to_split_on, image_merge_distance) {
  # Split the metadata by the values in the requested column. Note that for the mission level, this
  # approach is overkill because it will only be passed one mission of metadata at a time. It is a
  # vestige of a former way it was being applied. But it is still useful for the sub-mission level.
  split_metadata = split(image_metadata, image_metadata[[column_to_split_on]])
  # Those unique values are used as the dataset_id for logging purposes
  dataset_ids = names(split_metadata)

  # Extract the polygons for each chunk (mission or sub-mission) of metadata
  polygons_and_inds = purrr::map2(
    split_metadata,
    dataset_ids,
    extract_mission_polygon,
    image_merge_distance = image_merge_distance,
    identify_images_in_polygon = TRUE
  )

  # get count of retained images per dataset
  intersection_image_ids = purrr::map(polygons_and_inds, "intersection_image_ids")
  image_counts = purrr::map(intersection_image_ids, length)

  # get the ones that had < 10 images retained
  too_few_images_bool = image_counts < 10

  if (any(too_few_images_bool)) {
    too_few_images_names = names(which(too_few_images_bool))
    too_few_images_names_string = paste(too_few_images_names, collapse = ", ")
    warning(
      paste0("The following missions/sub-missions had fewer than 10 images retained in the computed polygons and were dropped entirely: ", too_few_images_names_string)
    )
  }

  # Drop the polygons (and image IDs) that had too few images retained
  polygons_and_inds = polygons_and_inds[!too_few_images_bool]

  polygons_sfc = purrr::map(polygons_and_inds, "polygon")
  polygons_sf = purrr::map(polygons_sfc, st_as_sf)
  intersection_image_ids = purrr::map(polygons_and_inds, "intersection_image_ids")

  intersection_image_ids = unlist(intersection_image_ids, use.names = FALSE)

  return(list(polygons = polygons_sf, retained_image_IDs = intersection_image_ids))
}


#' Read image metadata from file, auto-detecting format
#'
#' @param filepath Path to CSV or gpkg file
#' @return Data frame with image metadata (includes lon/lat columns for CSV, or sf object for gpkg)
read_image_metadata = function(filepath) {
  ext = tools::file_ext(filepath)

  if (ext == "csv") {
    df = read_csv(filepath, show_col_types = FALSE)
  } else if (ext == "gpkg") {
    sf_data = st_read(filepath, quiet = TRUE)
    # Extract coordinates and add as columns for compatibility
    coords = st_coordinates(sf_data)
    df = sf_data |>
      st_drop_geometry() |>
      mutate(lon = coords[, 1], lat = coords[, 2])
  } else {
    stop(paste("Unsupported file format:", ext))
  }

  # If column preprocessed_exif_gpstimestamp exists, coerce to character (from hms which was the
  # default of read_csv)
  if ("preprocessed_exif_gpstimestamp" %in% colnames(df)) {
    df = df |> mutate(preprocessed_exif_gpstimestamp = as.character(preprocessed_exif_gpstimestamp))
  }
}


#' Summarize EXIF metadata for a single mission
#'
#' This function performs the core summarization workflow:
#' 1. Reads image metadata from input path (CSV or gpkg)
#' 2. Computes polygons at mission and sub-mission level
#' 3. Filters to retained images (intersection of both polygon computations)
#' 4. Optionally writes filtered image metadata to output path
#' 5. Re-computes polygons with filtered images
#' 6. Extracts summary statistics using extract_imagery_dataset_metadata()
#' 7. Writes attributed polygons to output paths
#'
#' @param mission_id_foc Mission ID to process
#' @param input_metadata_path Directory or file path containing input image metadata.
#'   If directory: looks for {mission_id}.csv or {mission_id}_image-metadata.gpkg
#'   If file: uses directly (must be CSV or gpkg)
#' @param output_derived_mission_path Directory for output mission-level derived metadata gpkg
#' @param output_derived_sub_mission_path Directory for output sub-mission-level derived metadata gpkg
#' @param output_retained_images_path Optional: directory for output filtered image metadata gpkg.
#'   If NULL, filtered images are not written (useful for post-curation where images are already filtered)
#' @param image_merge_distance Distance threshold for polygon computation
#' @return TRUE on success, FALSE on failure (with warning)
summarize_mission_exif = function(mission_id_foc,
                                   input_metadata_path,
                                   output_derived_mission_path,
                                   output_derived_sub_mission_path,
                                   output_retained_images_path = NULL,
                                   image_merge_distance) {

  # Determine input file path
  if (dir.exists(input_metadata_path)) {
    # It's a directory - look for the file
    csv_path = file.path(input_metadata_path, paste0(mission_id_foc, ".csv"))
    gpkg_path = file.path(input_metadata_path, paste0(mission_id_foc, "_image-metadata.gpkg"))

    if (file.exists(csv_path)) {
      input_filepath = csv_path
    } else if (file.exists(gpkg_path)) {
      input_filepath = gpkg_path
    } else {
      warning(paste("No image metadata file found for mission", mission_id_foc, "in", input_metadata_path))
      return(FALSE)
    }
  } else if (file.exists(input_metadata_path)) {
    # It's a file path
    input_filepath = input_metadata_path
  } else {
    warning(paste("Input path does not exist:", input_metadata_path))
    return(FALSE)
  }

  # Read image metadata
  image_metadata = tryCatch({
    read_image_metadata(input_filepath)
  }, error = function(e) {
    warning(paste("Failed to read image metadata for mission", mission_id_foc, ":", e$message))
    return(NULL)
  })

  if (is.null(image_metadata)) return(FALSE)

  # Compute polygons at mission level
  mission_res = tryCatch({
    compute_polygons_and_images_retained(
      image_metadata = image_metadata,
      column_to_split_on = "mission_id",
      image_merge_distance = image_merge_distance
    )
  }, error = function(e) {
    warning(paste("Failed to compute mission polygon for", mission_id_foc, ":", e$message))
    return(NULL)
  })

  if (is.null(mission_res)) return(FALSE)

  # Compute polygons at sub-mission level
  sub_mission_res = tryCatch({
    compute_polygons_and_images_retained(
      image_metadata = image_metadata,
      column_to_split_on = "sub_mission_id",
      image_merge_distance = image_merge_distance
    )
  }, error = function(e) {
    warning(paste("Failed to compute sub-mission polygons for", mission_id_foc, ":", e$message))
    return(NULL)
  })

  if (is.null(sub_mission_res)) return(FALSE)

  # Compute the images that were retained in both the mission polygons and the sub-mission polygons
  images_retained_in_both = intersect(
    mission_res$retained_image_IDs, sub_mission_res$retained_image_IDs
  )

  # Filter the image metadata to only include data for those images
  image_metadata = image_metadata |> filter(image_id %in% images_retained_in_both)

  if (nrow(image_metadata) < 10) {
    warning(paste("Fewer than 10 images retained for mission", mission_id_foc, "after polygon filtering"))
    return(FALSE)
  }

  # Write filtered image metadata if output path is provided
  if (!is.null(output_retained_images_path)) {
    create_dir(output_retained_images_path)
    metadata_perimage_subset_filepath = file.path(
      output_retained_images_path,
      paste0(mission_id_foc, "_image-metadata.gpkg")
    )
    image_metadata_sf = image_metadata |>
      st_as_sf(coords = c("lon", "lat"), crs = 4326)
    st_write(image_metadata_sf, metadata_perimage_subset_filepath, delete_dsn = TRUE, quiet = TRUE)
  }

  # Re-compute the polygons now that extraneous images have been filtered out
  mission_res = compute_polygons_and_images_retained(
    image_metadata = image_metadata,
    column_to_split_on = "mission_id",
    image_merge_distance = image_merge_distance
  )
  sub_mission_res = compute_polygons_and_images_retained(
    image_metadata = image_metadata,
    column_to_split_on = "sub_mission_id",
    image_merge_distance = image_merge_distance
  )

  # Extract the summary statistics for the mission
  summary_mission = extract_imagery_dataset_metadata(
    metadata = image_metadata,
    # only one set of polygons is expected, since only passed one mission of data, so take the first index
    mission_polygon = mission_res$polygons[[1]],
    dataset_id = mission_id_foc
  )

  # Attribute the polygon with the metadata
  mission_poly = mission_res$polygons[[1]]
  mission_poly_attributed = bind_cols(mission_poly, summary_mission)

  # Give it a mission ID column
  mission_poly_attributed$mission_id = mission_id_foc

  # Write mission-level metadata
  create_dir(output_derived_mission_path)
  metadata_per_mission_filepath = file.path(output_derived_mission_path, paste0(mission_id_foc, ".gpkg"))
  st_write(
    mission_poly_attributed,
    metadata_per_mission_filepath,
    delete_dsn = TRUE,
    quiet = TRUE
  )

  # Now repeat for each sub-mission within the mission
  sub_mission_ids = unique(image_metadata$sub_mission_id)
  create_dir(output_derived_sub_mission_path)

  for (sub_mission_id_foc in sub_mission_ids) {
    # Pull the metadata for the images from this sub-mission
    image_metadata_sub_mission = image_metadata |> filter(sub_mission_id == sub_mission_id_foc)

    # Pull the polygon for this sub-mission
    polygon_sub_mission_foc = sub_mission_res$polygons[[sub_mission_id_foc]]

    if (is.null(polygon_sub_mission_foc)) {
      warning(paste("No polygon found for sub-mission", sub_mission_id_foc))
      next
    }

    # Extract the summary statistics for this sub-mission
    summary_sub_mission = extract_imagery_dataset_metadata(
      metadata = image_metadata_sub_mission,
      mission_polygon = polygon_sub_mission_foc,
      dataset_id = sub_mission_id_foc
    )

    # Attribute the polygon with the metadata
    sub_mission_poly_attributed = bind_cols(polygon_sub_mission_foc, summary_sub_mission)

    # Give it a mission ID and sub-mission ID column
    sub_mission_poly_attributed$mission_id = mission_id_foc
    sub_mission_poly_attributed$sub_mission_id = sub_mission_id_foc

    metadata_per_sub_mission_filepath = file.path(
      output_derived_sub_mission_path,
      paste0(sub_mission_id_foc, ".gpkg")
    )

    # Write
    st_write(
      sub_mission_poly_attributed,
      metadata_per_sub_mission_filepath,
      delete_dsn = TRUE,
      quiet = TRUE
    )
  }

  return(TRUE)
}
