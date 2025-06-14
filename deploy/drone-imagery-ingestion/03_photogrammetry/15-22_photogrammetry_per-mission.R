# Purpose: Run the photogrammetry (and post-processing) pipeline for a given mission ID and upload
# the results to CyVerse

# If COPC conversion is desired, this script requires the existince of a conda environment named `untwine` with the `untwine` package
# installed. This can be created with the system command: `conda create -n untwine -c conda-forge untwine
# -y`

# IMPORTANT NOTE: You must have already configured an rclone remote on this machine for the object
# store. This is already done automatically on ofo dev images. The config should look like this:
# https://github.com/open-forest-observatory/ofo-ansible/blob/main/roles/ofo/files/rclone.conf In
# addition, you must have your S3 credentials set in the environment variables RCLONE_S3_ACCESS_KEY_ID and
# RCLONE_S3_SECRET_ACCESS_KEY. You can do this by running the following command (change the values to
# yours): 'export RCLONE_S3_ACCESS_KEY_ID=your_access_key_id; export
# RCLONE_S3_SECRET_ACCESS_KEY=your_secret_access_key' and then reboot.

library(tidyverse)
library(sf)
library(lidR)
library(terra)
library(purrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")

source("deploy/drone-imagery-ingestion/03_photogrammetry/src/17_download-unzip-images.R")
source("deploy/drone-imagery-ingestion/03_photogrammetry/src/18_prep-metashape-configs.R")
source("deploy/drone-imagery-ingestion/03_photogrammetry/src/19_run-metashape.R")
source("deploy/drone-imagery-ingestion/03_photogrammetry/src/19b_upload-metashape-products.R")
source("deploy/drone-imagery-ingestion/03_photogrammetry/src/20_postprocess-photogrammetry-products.R")
source("deploy/drone-imagery-ingestion/03_photogrammetry/src/21_upload-postprocessed-photogrammetry.R")


photogrammetry_per_mission = function(mission_id, config_id) {

  download_unzip_images(mission_id)
  config_filename = prep_metashape_config(mission_id, config_id)
  # config_filename has the form <config_id>_<mission_id> so it's redundant with the mission_id and
  # config_id arguments passed in to this function, but this should be more robust in case the
  # config filename convention changes.
  run_metashape(config_filename)
  upload_photogrammetry_outputs_to_object_store(mission_id, config_id)
  postprocess_photogrammetry(mission_id, config_id)
  upload_postprocessed_photogrammetry_to_object_store(mission_id, config_id)

  return(TRUE)

}
