# Purpose: Get the mission IDs of the missions to process, using one of several options, and write a
# CSV of the mission IDs to be referenced by the next scripts in the pipeline.

# TODO: This assumes that all missions that were processed for raw imagery metadata have also been
# processed for raw imagery files, which is not necessarily true so we should instead more directly
# find the processed raw files (by querying cyverse?) and run photogrammetry for those.

library(tidyverse)
library(purrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")

get_missions_per_project = function(project_name) {

  # Read the crosswalk file that contains the mission IDs for the project
  crosswalk_file = file.path(CONTRIBUTED_TO_SORTED_MISSION_ID_CROSSWALK_PATH, paste0(project_name, ".csv"))
  crosswalk = read_csv(crosswalk_file)

  # Extract the mission IDs from the crosswalk
  missions_to_process = crosswalk$mission_id |> unique()

  # Format as data frame
  missions_to_process_df = data.frame(mission_id = missions_to_process)

  return(missions_to_process_df)

}

missions_to_process_df = map_df(PROJECT_NAMES_TO_PROCESS_PHOTOGRAMMETRY, get_missions_per_project)


# Write
write_csv(missions_to_process_df, MISSIONS_TO_PROCESS_PHOTOGRAMMETRY_LIST_PATH)
