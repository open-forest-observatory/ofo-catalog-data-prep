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


photogrammetry_per_mission = function(mission_id) {

  download_unzip_images(mission_id)
  prep_metashape_config(mission_id)
  run_metashape(mission_id)
  # TODO: Eventually, run_metashape should return the run_id, but for now we are inferring it as the
  # most recent (i.e. highest) run_id in the directory for the mission. This is being determine in
  # the next step, postprocess_photogrammetry, which then passes the run_id on to the next step.
  run_id = postprocess_photogrammetry(mission_id)
  upload_postprocessed_photogrammetry_to_cyverse(mission_id, run_id)

  return(TRUE)

}
