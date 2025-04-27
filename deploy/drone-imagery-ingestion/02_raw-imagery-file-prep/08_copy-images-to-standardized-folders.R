# Purpose: Read the processed image-level metadata data which contains the plan for sorting the
# images (origin and destination filepaths). Copy the images to the new folder structure based on
# the plan.

library(tidyverse)
library(furrr)
library(sf)

copy_mission_images = function(mission_id_foc) {
  image_metadata_file = file.path(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH, paste0(mission_id_foc, "_image-metadata.gpkg"))
  image_metadata = st_read(image_metadata_file)

  # Determine the absolute input and output paths
  image_metadata$image_path_contrib_abs = file.path(
    CONTRIBUTED_IMAGERY_PATH,
    image_metadata$image_path_contrib
  )

  image_metadata$image_path_ofo_abs = file.path(
    SORTED_IMAGERY_PATH,
    image_metadata$image_path_ofo
  )

  # Create the output folder(s)
  folders_out_abs = unique(dirname(image_metadata$image_path_ofo_abs))
  sapply(folders_out_abs, dir.create, recursive = TRUE, showWarnings = FALSE)

  # Copy files as hardlinks
  walk2(image_metadata$image_path_contrib_abs, image_metadata$image_path_ofo_abs, file.link)
}
