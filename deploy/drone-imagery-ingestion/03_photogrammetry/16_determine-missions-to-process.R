# Purpose: Get the mission IDs of the missions to process, using one of several options, and write a
# CSV of the mission IDs to be referenced by the next scripts in the pipeline.

# TODO: This assumes that all missions that were processed for raw imagery metadata have also been
# processed for raw imagery files, which is not necessarily true so we should instead more directly
# find the processed raw files (by querying cyverse?) and run photogrammetry for those.

library(tidyverse)
library(purrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")

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

## As an alternative, break into chunks and write to a separate file per chunk
create_dir(file.path(MISSIONS_TO_PROCESS_PHOTOGRAMMETRY_PATH, "chunks"))

chunks = chunk_up_collate(missions_to_process_df$mission_id, N_CHUNKS_PHOTOGRAMMETRY)

for (i in seq_along(chunks)) {
  chunk = chunks[[i]]
  chunk_df = data.frame(mission_id = chunk)
  chunk_filename = paste0("missions-to-process_", i, ".csv")
  chunk_filepath = file.path(MISSIONS_TO_PROCESS_PHOTOGRAMMETRY_PATH, "chunks", chunk_filename)
  write_csv(chunk_df, chunk_filepath)
}
