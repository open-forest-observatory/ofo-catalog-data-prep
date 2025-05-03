# Purpose: Run the imagery prep steps 08 through 14 for a set of specified mission IDs.

# IMPORTANT NOTE: You must have already authenticated irods with your CyVerse account on this
# machine using iinit. Note that by default, OFO dev instances come pre-authenticated for an
# anonymous (read-only) user, so you will need to run the following lines (change the username to
# yours) to authenticate as yourself and then type in your password when prompted
## echo '{"irods_host": "data.cyverse.org", "irods_port": 1247, "irods_user_name": "djyoung", "irods_zone_name": "iplant"}' > /home/exouser/.irods/irods_environment.json; iinit

library(tidyverse)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")

# Load the function that calls steps 15-22 for a specified mission ID
source("deploy/drone-imagery-ingestion/03_photogrammetry/15-22_photogrammetry_per-mission.R")

# # Run this once per processing run: Determine the missions to process. This depends on having set the project name in the file
# # "imagery-project-to-process-raw-imagery-metadata.txt". It puts the list of missions to process in the file
# # "missions-to-process.csv". Note that the script that is sourced here does not define a function,
# # it actually runds the code complete the task described in this comment.
# source("deploy/drone-imagery-ingestion/03_photogrammetry/16_determine-missions-to-process.R")

# If there is a command line arg, use that as the chunk of missions to process and load the list of
# mission IDs in that chunk
args = commandArgs(trailingOnly = TRUE)
chunk = args[1]
cat("Processing chunk ", chunk, "\n")

if (!is.na(chunk)) {
  # Read the chunk file
  chunk_filepath = file.path(MISSIONS_TO_PROCESS_PHOTOGRAMMETRY_PATH, "chunks", paste0("missions-to-process_", chunk, ".csv"))
  missions_to_process = read_csv(chunk_filepath) |> pull(mission_id)
} else {
  # If no chunk is specified, use all the missions to process
  missions_to_process = read_csv(MISSIONS_TO_PROCESS_PHOTOGRAMMETRY_LIST_PATH) |> pull(mission_id)
}

# Use (and update) the following if you want to override the missions to process determined above, such as if
# the final check revealed some missions were not completed
# missions_to_process = c("000574", "000565", "000570", missions_to_process)
# missions_to_process = missions_to_process[!duplicated(missions_to_process)]

purrr::walk(
  missions_to_process,
  photogrammetry_per_mission,
  .progress = TRUE
)

# Of the missions_to_process, see which (if any) have not been uploaded (based on which ones still
# have local files -- since a successful upload is followed by local deletion) and return a vector
# of those mission IDs for re-upload.
get_unuploaded_missions = function(missions_to_process) {
  local_postprocessed_folders = list.files(
    path = file.path(
      PHOTOGRAMMETRY_DIR,
      PHOTOGRAMMETRY_POSTPROCESSED_SUBDIR
    ),
    pattern = "^processed_",
    full.names = FALSE,
    recursive = TRUE,
    include.dirs = TRUE
  )

  # Extract the mission_id and run_id from the returned folder names, which have the format
  # {mission_id}/processed_{run_id}

  folder_parts = str_split(local_postprocessed_folders, "/")
  mission_ids = map_chr(folder_parts, 1)
  run_ids = map_chr(folder_parts, 2) |> str_remove("processed_")
  local_mission_and_run = data.frame(
    mission_id = mission_ids,
    run_id = run_ids
  )

  # Which mission IDs are still present locally?
  unuploaded_missions = local_mission_and_run |>
    filter((mission_id %in% missions_to_process))

  return(unuploaded_missions)
}


# Get the unuploaded missions
unuploaded_missions = get_unuploaded_missions(missions_to_process)

# Re-attempt upload of the unuploaded missions
if (nrow(unuploaded_missions) > 0) {
  cat("\n\n **** Re-attempting upload of unuploaded missions **** \n")
  for (i in 1:nrow(unuploaded_missions)) {
    mission_id = unuploaded_missions$mission_id[i]
    run_id = unuploaded_missions$run_id[i]
    upload_postprocessed_photogrammetry_to_cyverse(mission_id, run_id)
  }
} else {
  cat("\n\n **** All missions have been uploaded successfully. **** \n")
}

# Are there still any unuploaded missions? If so, print a warning and `dput` the data frame of
# unuploaded missions.

if (nrow(unuploaded_missions) > 0) {
  cat(
    "\n\n **** WARNING: The following missions were still not uploaded successfully after a second attempt:",
    paste(unuploaded_missions$mission_id, collapse = ", "),
    " **** \n")
  dput(unuploaded_missions)
}
