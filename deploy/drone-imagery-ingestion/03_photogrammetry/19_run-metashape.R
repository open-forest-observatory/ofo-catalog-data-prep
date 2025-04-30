# TODO next step: Docker run for automate-metashape:
# docker run -v /ofo-share/catalog-data-prep/02_photogrammetry:/data -e AGISOFT_FLS=$AGISOFT_FLS ghcr.io/open-forest-observatory/automate-metashape --config_file /data/04_derived-metashape-configs/01/000574.yml

# Purpose: Prepare the metashape config and shell files for processing all image sets

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/photogrammetry-prep.R")
source("src/utils.R")

run_metashape = function(mission_id_foc) {

  config_filpath_in_container = paste0("/data/", DERIVED_METASHAPE_CONFIG_SUBDIR, "/", mission_id_foc, ".yml")

  # Template for docker run command: docker run -v /ofo-share/catalog-data-prep/02_photogrammetry:/data -e AGISOFT_FLS=$AGISOFT_FLS ghcr.io/open-forest-observatory/automate-metashape --config_file /data/04_derived-metashape-configs/01/000574.yml
  docker_command = paste0(
    "docker run -v ", PHOTOGRAMMETRY_DIR, ":/data -e AGISOFT_FLS=$AGISOFT_FLS ",
    "ghcr.io/open-forest-observatory/automate-metashape --config_file ", config_filpath_in_container
  )

  # Run it
  docker_result = system(docker_command)

  if (docker_result != 0) {
    warning("Automate-metashape run failed with error code: ", docker_result)
    return(FALSE)
  }

  return(TRUE)
}

mission_id_foc = "000574"
run_metashape(mission_id_foc)