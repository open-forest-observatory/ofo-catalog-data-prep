# Purpose: Run the photogrammetry (and post-processing) pipeline for a given mission ID and upload
# the results to CyVerse

# This script requires the existince of a conda environment named `untwine` with the `untwine` package
# installed. This can be created with the system command: `conda create -n untwine -c conda-forge untwine
# -y`

# IMPORTANT NOTE: You must have already authenticated irods with your CyVerse account on this
# machine using iinit. Note that by default, OFO dev instances come pre-authenticated for an
# anonymous (read-only) user, so you will need to run the following lines (change the username to
# yours) to authenticate as yourself and then type in your password when prompted
## echo '{"irods_host": "data.cyverse.org", "irods_port": 1247, "irods_user_name": "djyoung", "irods_zone_name": "iplant"}' > /home/exouser/.irods/irods_environment.json; iinit

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
source("deploy/drone-imagery-ingestion/03_photogrammetry/src/20_postprocess-photogrammetry-products.R")
source("deploy/drone-imagery-ingestion/03_photogrammetry/src/21_upload-postprocessed-photogrammetry.R")


photogrammetry_per_mission = function(mission_id, config_id) {

  download_unzip_images(mission_id)
  config_filename = prep_metashape_config(mission_id, config_id)
  # config_filename has the form <mission_id>_<config_id> so it's redundant with the mission_id and
  # config_id arguments passed in to this function, but this should be more robust in case the
  # config filename convention changes.
  run_metashape(config_filename)
  postprocess_photogrammetry(mission_id, config_id)
  upload_postprocessed_photogrammetry_to_data_store(mission_id, config_id)

  return(TRUE)

}
