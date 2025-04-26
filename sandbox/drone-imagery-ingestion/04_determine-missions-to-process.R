# Purpose: Get the mission IDs of the missions to process, using one of several options, and write a
# CSV of the mission IDs to be referenced by the next scripts in the pipeline.

library(tidyverse)

source("sandbox/drone-imagery-ingestion/00_set-constants.R")

# Read the crosswalk file that contains the mission IDs for the project
crosswalk_file = file.path(CONTRIBUTED_TO_SORTED_MISSION_ID_CROSSWALK_PATH, paste0(IMAGERY_PROJECT_NAME, ".csv"))
crosswalk = read_csv(crosswalk_file)

# Extract the mission IDs from the crosswalk
missions_to_process = crosswalk$mission_id |> unique()

# If a subset of missions is specified, filter the missions to process
if (!is.null(IMAGERY_PROJECT_SUBSET_MISSIONS) && length(IMAGERY_PROJECT_SUBSET_MISSIONS) > 0) {
  missions_to_process = missions_to_process[missions_to_process %in% IMAGERY_PROJECT_SUBSET_MISSIONS]
}

# Format as data frame
missions_to_process_df = data.frame(mission_id = missions_to_process)

# Write
write_csv(missions_to_process_df, MISSIONS_TO_PROCESS_LIST_PATH)
