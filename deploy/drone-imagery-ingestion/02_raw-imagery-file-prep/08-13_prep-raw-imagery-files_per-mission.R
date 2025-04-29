# Purpose: Run the imagery prep steps 08 through 14 for a specified mission ID. This begins with (as
# input) the contributor-provided raw imagery folder and parsed metadata at the mission level
# (including which images belong with which mission and sub-mission and their source and destination
# file paths) and will result in the final curated set of raw imagery (zipped) and associated metada for
# the mission uploaded to the cyverse data store.

library(tidyverse)

source("deploy/drone-imagery-ingestion/00_set-constants.R")

# Load the defs of the fucntions used for steps 08 through 14
source("deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/08_copy-images-to-standardized-folders.R")
source("deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/09_fix-exif.R")
source("deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/10_raw-imagery-thumbnails-and-zip.R")
source("deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/11_copy-raw-imagery-to-upload-staging-dir.R")
source("deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/12_upload-raw-imagery-to-cyverse.R")
source("deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/13_delete-prepped-raw-imagery.R")


prep_raw_imagery_files_per_mission = function(mission_id_foc) {

  # # Check if it's already completed (based on existence of the zip file on cyverse) and skip if so
  # # But note that currently this won't work when parallelized because irods blocks the simultaneous
  # # queries. Such cases result in us assuming the mission is not completed and running the pipeline
  # # for it. 
  # already_run = check_if_zip_file_exists(mission_id_foc)
  # if (already_run) {
  #   cat("Mission", mission_id_foc, "already completed. Skipping.\n")
  #   return(FALSE)
  # }

  copy_mission_images(mission_id_foc)
  fix_exif(mission_id_foc)
  make_raw_imagery_thumbnails_and_zip(mission_id_foc)
  copy_raw_imagery_to_upload_staging_dir(mission_id_foc)
  upload_raw_imagery_to_cyverse(mission_id_foc)
  delete_prepped_raw_imagery(mission_id_foc)

  return(TRUE)
}