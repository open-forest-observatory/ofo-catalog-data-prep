# Purpose: At the mission level and the sub-mission level (separately), merge the human-provided Baserow metadata
# and the summarized (mission or sub-mission level) EXIF metadata extracted from the drone imagery.

library(tidyverse)
library(sf)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/merge-utils.R")


## Workflow

# Determine which missions to process
missions_to_process = read_csv(MISSIONS_TO_PROCESS_RAW_IMAGERY_METADATA_LIST_PATH) |>
  pull(mission_id)

# Create the output folders
create_dir(FULL_METADATA_PER_MISSION_PATH)
create_dir(FULL_METADATA_PER_SUB_MISSION_PATH)


# Parallelize across all selected missions using the shared merge function
future::plan(multisession)
future_walk(
  missions_to_process,
  ~ merge_derived_and_contributed_metadata(
    mission_foc = .x,
    derived_mission_path = DERIVED_METADATA_PER_MISSION_PATH,
    derived_sub_mission_path = DERIVED_METADATA_PER_SUB_MISSION_PATH,
    contributed_mission_path = EXTRACTED_METADATA_PER_MISSION_PATH,
    contributed_sub_mission_path = EXTRACTED_METADATA_PER_SUB_MISSION_PATH,
    output_mission_path = FULL_METADATA_PER_MISSION_PATH,
    output_sub_mission_path = FULL_METADATA_PER_SUB_MISSION_PATH
  ),
  .progress = TRUE
)
