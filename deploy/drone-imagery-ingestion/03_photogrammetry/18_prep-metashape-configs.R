# Purpose: Prepare the metashape config and shell files for processing all image sets

library(dplyr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/photogrammetry-prep.R")
source("src/utils.R")

mission_id_foc = "000574"

# Folder relative to the path mounted to the docker container, which should be PHOTOGRAMMETRY_DIR
imagery_folder = file.path(PHOTOGRAMMETRY_DIR, INPUT_IMAGES_SUBDIR, mission_id_foc)

# Processing

# Add subdirs of the mission folder (sub-missions) to the config as separate folders so
# that they can be calibrated differently (e.g. for different cameras)
subdirs = list.dirs(imagery_folder, recursive = FALSE)

# Remove the PHOTOGRAMMETRY_DIR from the base of the path because we are going to mount that to the
# docker container and we just want to pass the path below that
photogrammetry_dir_nofinalslash = str_remove(PHOTOGRAMMETRY_DIR, "/$")
subdirs = str_remove(subdirs, fixed(paste0(photogrammetry_dir_nofinalslash, "/")))
photo_paths_in_container = paste0("/data/", subdirs)

config_overrides = data.frame(
  photo_path = I(list(photo_paths_in_container)),
  config_filename = mission_id_foc,
  output_path = file.path("/data", METASHAPE_OUTPUT_SUBDIR),
  project_path = file.path("/data", METASHAPE_PROJECT_SUBDIR)
)

# Make the folder that will hold the derived configs
create_dir(file.path(PHOTOGRAMMETRY_DIR, DERIVED_METASHAPE_CONFIG_SUBDIR))

make_derived_yaml(
  cfg_base_filepath = file.path(PHOTOGRAMMETRY_DIR, BASE_METASHAPE_CONFIG_FILE_SUBPATH),
  replacements = config_overrides,
  derived_yaml_dir = file.path(PHOTOGRAMMETRY_DIR, DERIVED_METASHAPE_CONFIG_SUBDIR)
)

# TODO next step: Docker run for automate-metashape:
# docker run -v /ofo-share/catalog-data-prep/02_photogrammetry:/data -e AGISOFT_FLS=$AGISOFT_FLS ghcr.io/open-forest-observatory/automate-metashape --config_file /data/04_derived-metashape-configs/01/000574.yml