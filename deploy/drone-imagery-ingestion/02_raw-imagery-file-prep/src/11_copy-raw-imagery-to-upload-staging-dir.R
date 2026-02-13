# Purpose: Take the imagery and associated metadata products (image points and mission polygons) for
# each mission and hardlink them all into a single file tree (all under one mission folder) that is
# ready to be uploaded to CyVerse.

library(furrr)
library(tidyverse)

copy_raw_imagery_to_upload_staging_dir = function(mission_id_foc, use_post_curation = TRUE) {
  # For tracking down the mission ID(s) that produces warnings when this function is called inside a
  # map() function, you can include this line and see which mission ID warning was printed just
  # before the real warning. warning(paste("Parsing EXIF for mission ID:", mission_id_foc))
  cat("\n **** Copying raw imagery to upload staging dir for mission", mission_id_foc, "**** \n")

  # Select appropriate metadata paths
  if (use_post_curation) {
    mission_metadata_path = POST_CURATION_FULL_METADATA_PER_MISSION_PATH
    image_metadata_path = POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH
  } else {
    mission_metadata_path = FULL_METADATA_PER_MISSION_PATH
    image_metadata_path = PARSED_EXIF_FOR_RETAINED_IMAGES_PATH
  }

  # MISSION FOOTPRINTS

  # Symlink all files from the publishable mission footprints dir to the unified publishable tree
  # First make the directories
  infile = file.path(mission_metadata_path, paste0(mission_id_foc, "_mission-metadata.gpkg"))
  outfile = file.path(UPLOAD_STAGING_DIR_PATH, mission_id_foc, "metadata-mission", paste0(mission_id_foc, "_mission-metadata.gpkg"))
  outdir = dirname(outfile)
  dir.create(outdir, recursive = TRUE)
  
  file.symlink(infile, outfile)


  # MISSION POINTS

  # Symlink all files from the publishable mission points dir to the unified publishable tree
  # First make the directories
  infile = file.path(image_metadata_path, paste0(mission_id_foc, "_image-metadata.gpkg"))
  outfile = file.path(UPLOAD_STAGING_DIR_PATH, mission_id_foc, "metadata-images", paste0(mission_id_foc, "_image-metadata.gpkg"))
  outdir = dirname(outfile)
  dir.create(outdir, recursive = TRUE)

  file.symlink(infile, outfile)


  # IMAGES

  # Symlink all files from the publishable images dir to the unified publishable tree
  # First make the directories
  infiles = list.files(file.path(IMAGERY_ZIP_AND_EXAMPLES_PATH, mission_id_foc), full.names = FALSE, recursive = TRUE)
  infiles_full = file.path(IMAGERY_ZIP_AND_EXAMPLES_PATH, mission_id_foc, infiles)
  outfiles = file.path(UPLOAD_STAGING_DIR_PATH, mission_id_foc, infiles)
  outdirs = unique(dirname(outfiles))

  walk(outdirs, dir.create, recursive = TRUE)

  file.symlink(infiles_full, outfiles)

  return(TRUE)

}
