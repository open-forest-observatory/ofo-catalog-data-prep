# Purpose: Delete the folders of prepped raw imagery files (e.g. zipped imagery folder, example
# images, metadata) that have now been uploaded to CyVerse. TODO: Could delete even farther back the
# pipeline: the mission images in the 2_sorted folder. But (unless there was an EXIF issue to fix)
# these files are all hardlinks of the contributed imagery files, so deleting them would not save
# much space.

delete_prepped_raw_imagery = function(mission_id_foc) {

  # Raw imagery zip and example images
  dir = file.path(IMAGERY_ZIP_AND_EXAMPLES_PATH, mission_id_foc)
  unlink(dir, recursive = TRUE)

  # Upload staging directory
  dir = file.path(UPLOAD_STAGING_DIR_PATH, mission_id_foc)
  unlink(dir, recursive = TRUE)
}
