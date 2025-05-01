# Purpose: Run the photogrammetry (and post-processing) pipeline for a given mission ID and upload
# the results to CyVerse

# This script requires the existince of a conda environment named `untwine` with the `untwine` package
# installed. This can be created with the system command: `conda create -n untwine -c conda-forge untwine
# -y`

library(tidyverse)
library(sf)
library(lidR)
library(terra)
library(purrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")


photogrammetry_per_mission = function(mission_id) {

  source("deploy/drone-imagery-ingestion/03_photogrammetry/src/17_download-unzip-images.R")
  source("deploy/drone-imagery-ingestion/03_photogrammetry/src/18_prep-metashape-configs.R")
  source("deploy/drone-imagery-ingestion/03_photogrammetry/src/19_run-metashape.R")
  source("deploy/drone-imagery-ingestion/03_photogrammetry/src/20_postprocess-photogrammetry-products.R")

  download_unzip_images(mission_id)
  prep_metashape_config(mission_id)
  run_metashape(mission_id)
  postprocess_photogrammetry(mission_id)

  return()

}
