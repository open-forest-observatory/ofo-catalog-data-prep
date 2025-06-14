# TODO next step: Docker run for automate-metashape:
# docker run -v /ofo-share/catalog-data-prep/02_photogrammetry:/data -e AGISOFT_FLS=$AGISOFT_FLS ghcr.io/open-forest-observatory/automate-metashape --config_file /data/04_derived-metashape-configs/01/000574.yml

# Purpose: Prepare the metashape config and shell files for processing all image sets

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/photogrammetry-prep.R")
source("src/utils.R")

run_metashape = function(config_filename) {

  config_filpath_in_container = paste0("/data/", DERIVED_METASHAPE_CONFIG_SUBDIR, "/", config_filename, ".yml")

  # Template for docker run command: docker run -v /ofo-share/catalog-data-prep/02_photogrammetry:/data -e AGISOFT_FLS=$AGISOFT_FLS ghcr.io/open-forest-observatory/automate-metashape --config_file /data/04_derived-metashape-configs/01/000574.yml
  docker_command = paste0(
    "docker run -v ", PHOTOGRAMMETRY_DIR, ":/data -e AGISOFT_FLS=$AGISOFT_FLS --gpus all --pull always ",
    "ghcr.io/open-forest-observatory/automate-metashape --config_file ", config_filpath_in_container
  )

  # Run it
  docker_result = system(docker_command)

  n_tries = 1

  while (docker_result != 0 & n_tries < 3) {
    warning("Automate-metashape run for mission ID", mission_id_foc, "failed with error code: ", docker_result, ". Retrying (attempt ", n_tries, "of 3) in 30 seconds.")
    
    Sys.sleep(30)

    docker_result = system(docker_command)
    n_tries = n_tries + 1
  }

  # If it still didn't work, return false with a warning
  if (docker_result != 0) {
    warning("Automate-metashape run for mission ID", mission_id_foc, "failed with error code: ", docker_result, ". Giving up after 3 attempts.")
    return(FALSE)
  }


  # Delete the metashape project file (and folder) for this mission and config
  project_files = list.files(
    path = file.path(PHOTOGRAMMETRY_DIR, METASHAPE_PROJECT_SUBDIR),
    pattern = paste0(config_id_foc, "_", mission_id_foc, "\\.*"),
    full.names = TRUE
  )

  # Surround project-files in single quotes
  project_files = paste0("'", project_files, "'")

  # Try the deletion as a system command
  command = paste0("sudo rm -rf ", paste(project_files, collapse = " "))
  result = system(command)

  # If the containing folder (photogrammetry-outputs) is now empty, delete that too
  mission_dir = file.path(PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_SUBDIR)
  files_in_mission_dir = list.files(mission_dir)
  if (length(files_in_mission_dir) == 0) {
    unlink(mission_dir, recursive = TRUE)
  }

  # Delete the downloaded raw imagery files for this metashape run. TODO: If we will ever run more
  # than one metashape config of the same imagery mission at once, need to assign the downloaded
  # imagery unique folder names so we only delete the imagery from the config that completed.

  # Delete the unzipped folder of images.
  folder_to_delete = file.path(PHOTOGRAMMETRY_DIR, INPUT_IMAGES_SUBDIR, mission_id_foc)
  unlink(folder_to_delete, recursive = TRUE)

  return(TRUE)
}
