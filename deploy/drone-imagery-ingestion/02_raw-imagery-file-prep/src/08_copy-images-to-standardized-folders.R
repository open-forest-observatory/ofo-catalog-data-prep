# Purpose: Copy images to standardized folder structure using symlinks.
# For post-curation: uses POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH

library(tidyverse)
library(furrr)
library(sf)

source("src/utils.R")

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
  walk(folders_out_abs, create_dir)

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
