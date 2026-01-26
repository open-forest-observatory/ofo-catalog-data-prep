# Purpose: Read a CSV file containing the deduplicated EXIF data for all images in a mission,
# *at the mission level* (i.e., combining both dates of a two-date mission; both orientations of a
# two-part grid mission), and extract/process into human-readable metadata at the mission level and
# sub-mission level using the metadata summarization functions of the ofo package.

library(tidyverse)
library(sf)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/summarization-utils.R")


## Workflow

# Determine which missions to process
missions_to_process = read_csv(MISSIONS_TO_PROCESS_RAW_IMAGERY_METADATA_LIST_PATH) |>
  pull(mission_id)

# Create the output folders
create_dir(DERIVED_METADATA_PER_MISSION_PATH)
create_dir(DERIVED_METADATA_PER_SUB_MISSION_PATH)
create_dir(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH)


# Run for each mission_id using the shared summarization function
future::plan(multisession)
future_walk(
  missions_to_process,
  ~ summarize_mission_exif(
    mission_id_foc = .x,
    input_metadata_path = PARSED_EXIF_DEDUPED_PATH,
    output_derived_mission_path = DERIVED_METADATA_PER_MISSION_PATH,
    output_derived_sub_mission_path = DERIVED_METADATA_PER_SUB_MISSION_PATH,
    output_retained_images_path = PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
    image_merge_distance = IMAGE_MERGE_DISTANCE
  ),
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)
