# Purpose: Generate a list of mission IDs that belong to the NRS project.
# This list is used by get_merge_distance_for_mission() to determine which missions
# need the larger IMAGE_MERGE_DISTANCE_NRS_OVERRIDE value.
#
# Run this script once to generate the list, then the summarization scripts can
# use the list for fast lookups without loading full mission metadata.

library(tidyverse)
library(sf)

source("deploy/drone-imagery-ingestion/00_set-constants.R")

# Read the combined mission metadata
if (!file.exists(FULL_METADATA_PER_MISSION_COMBINED_FILEPATH)) {
  stop(paste(
    "Combined mission metadata file not found:",
    FULL_METADATA_PER_MISSION_COMBINED_FILEPATH,
    "\nRun the metadata pipeline first to generate this file."
  ))
}

mission_metadata = st_read(FULL_METADATA_PER_MISSION_COMBINED_FILEPATH, quiet = TRUE) |>
  st_drop_geometry()

# Check that project_id column exists
if (!"project_id" %in% names(mission_metadata)) {
  stop("project_id column not found in mission metadata")
}

# Filter for NRS project and extract mission IDs
nrs_missions = mission_metadata |>
  filter(project_id == NRS_PROJECT_ID) |>
  select(mission_id) |>
  distinct()

cat(sprintf("Found %d missions from NRS project (project_id: %s)\n",
            nrow(nrs_missions), NRS_PROJECT_ID))

# Create output directory if needed
output_dir = dirname(NRS_MISSIONS_LIST_FILEPATH)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created directory:", output_dir, "\n")
}

# Write the list
write_csv(nrs_missions, NRS_MISSIONS_LIST_FILEPATH)
cat("Wrote NRS mission list to:", NRS_MISSIONS_LIST_FILEPATH, "\n")
