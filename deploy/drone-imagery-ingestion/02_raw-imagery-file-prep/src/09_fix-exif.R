# Purpose: Fix EXIF metadata. If file is a symlink and needs modification,
# replace symlink with actual file copy first.
#
# Uses preprocessed_exif_Orientation and preprocessed_exif_GPSTimeStamp columns
# from image metadata gpkg (added by script 05) instead of sorting plan CSV.

library(furrr)
library(tidyverse)
library(sf)

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

  # Select appropriate metadata path
  if (use_post_curation) {
    metadata_path = POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH
  } else {
    metadata_path = PARSED_EXIF_FOR_RETAINED_IMAGES_PATH
  }

  metadata_file = file.path(metadata_path, paste0(mission_id_foc, "_image-metadata.gpkg"))

  if (!file.exists(metadata_file)) {
    warning(paste("No image metadata found for mission", mission_id_foc))
    return(FALSE)
  }

  # Read image metadata gpkg which contains preprocessed EXIF columns
  image_metadata = st_read(metadata_file, quiet = TRUE) |> st_drop_geometry()

  # Construct actual file paths from the sorted imagery location
  image_metadata = image_metadata |>
    mutate(filepath = file.path(SORTED_IMAGERY_PATH, image_path_ofo))

  # Check that files exist
  existing_files = file.exists(image_metadata$filepath)
  if (!all(existing_files)) {
    warning(paste("Some images not found for mission", mission_id_foc,
                  "- missing", sum(!existing_files), "files"))
    image_metadata = image_metadata |> filter(existing_files)
  }

  if (nrow(image_metadata) == 0) {
    warning(paste("No images found for mission", mission_id_foc))
    return(FALSE)
  }

  # Determine which images need fixing using preprocessed EXIF columns
  # Orientation needs fixing if != 1
  image_metadata$fix_orientation = !is.na(image_metadata$preprocessed_exif_Orientation) &
                                    image_metadata$preprocessed_exif_Orientation != 1

  # GPSTimeStamp needs fixing if it has decimal seconds (causes Metashape errors)
  image_metadata$fix_gpstimestamp = !is.na(image_metadata$preprocessed_exif_GPSTimeStamp) &
                                     grepl("[0-9]+:[0-9]+:[0-9]+\\.[0-9]+",
                                           image_metadata$preprocessed_exif_GPSTimeStamp)

  image_metadata = image_metadata |>
    mutate(fix_both = (fix_orientation & fix_gpstimestamp)) |>
    mutate(fix_orientation = ifelse(fix_both, FALSE, fix_orientation),
           fix_gpstimestamp = ifelse(fix_both, FALSE, fix_gpstimestamp))

  # Get files that need any fixing
  files_to_fix = image_metadata |> filter(fix_orientation | fix_gpstimestamp | fix_both)

  if (nrow(files_to_fix) > 0) {
    cat("Replacing symlinks with file copies for", nrow(files_to_fix), "files that need EXIF fixes...\n")
    walk(files_to_fix$filepath, replace_symlink_with_copy)
  }

  # Fix orientation only
  exif_to_fix_orientation = image_metadata |> filter(fix_orientation == TRUE)
  if (nrow(exif_to_fix_orientation) > 0) {
    cat("Fixing orientation flag for", nrow(exif_to_fix_orientation), "images\n")
    command = "exiftool -n -r -fast4 -overwrite_original -Orientation=1"
    run_cmd_chunks(command, exif_to_fix_orientation$filepath)
  }

  # Fix GPSTimeStamp only
  exif_to_fix_gpstimestamp = image_metadata |> filter(fix_gpstimestamp == TRUE)
  if (nrow(exif_to_fix_gpstimestamp) > 0) {
    cat("Fixing GPSTimeStamp for", nrow(exif_to_fix_gpstimestamp), "images\n")
    command = "exiftool -n -r -fast4 -overwrite_original -GPSTimeStamp="
    run_cmd_chunks(command, exif_to_fix_gpstimestamp$filepath)
  }

  # Fix both
  exif_to_fix_both = image_metadata |> filter(fix_both == TRUE)
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
  folder = file.path(SORTED_IMAGERY_PATH, mission_id_foc)
  exiftool_temp_files = list.files(folder, full.names = TRUE, recursive = TRUE,
                                    pattern = "_exiftool_tmp$")
  if (length(exiftool_temp_files) > 0) {
    warning("Exiftool temp files found in", mission_id_foc, "- incomplete operation")
    return(FALSE)
  }

  return(TRUE)
}
