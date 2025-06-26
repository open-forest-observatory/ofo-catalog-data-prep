# Purpose: Get the mission IDs of the missions to process, using one of several sources of mission
# lists, and exclude missions for which the expected final products all exist already, and write a
# CSV of the mission IDs to process, to be referenced by the next scripts in the pipeline.

# TODO: This assumes that all missions that were processed for raw imagery metadata have also been
# processed for raw imagery files, which is not necessarily true so we should instead more directly
# find the processed raw files (by querying jetstream?) and run photogrammetry for those.

# TODO: Check if the final products already exist and 

library(tidyverse)
library(purrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")
source("src/photogrammetry-prep.R")

EXPECTED_FILE_SUFFIXES_FULL = c("_cameras.xml", "_chm-mesh.tif", "_chm-ptcloud.tif",
  "_dsm-mesh.tif", "_dsm-ptcloud.tif", "_dtm-ptcloud.tif",
  "_log.txt", "_model-georeferenced.ply", "_model-local.ply",
  "_ortho-dsm-mesh.tif", "_ortho-dsm-ptcloud.tif",
  "_points-copc.laz", "_report.pdf")

EXPECTED_FILE_SUFFIXES_THUMBNAIL = c(
  "_chm-mesh.png", "_chm-ptcloud.png", "_dsm-mesh.png",
  "_dsm-ptcloud.png", "_dtm-ptcloud.png", "_ortho-dsm-mesh.png",
  "_ortho-dsm-ptcloud.png"
)


# Option 1: By project name
get_missions_per_project = function(project_name) {

  # Read the crosswalk file that contains the mission IDs for the project
  crosswalk_file = file.path(CONTRIBUTED_TO_SORTED_MISSION_ID_CROSSWALK_PATH, paste0(project_name, ".csv"))
  crosswalk = read_csv(crosswalk_file)

  # Extract the mission IDs from the crosswalk
  missions_to_process = crosswalk$mission_id |> unique()

  return(missions_to_process)

}

# missions_to_process = map(PROJECT_NAMES_TO_PROCESS_PHOTOGRAMMETRY, get_missions_per_project) |> unlist()


# Option 2: By object store listing
get_objectstore_missions_w_raw_imagery = function() {

  # Query the object store for a file listing
  remote_dir = paste0(RCLONE_REMOTE, ":", REMOTE_MISSIONS_DIR)
  command = paste("rclone lsf", remote_dir, "-R --files-only", sep = " ")
  listing = system(command, intern = TRUE)

  # Filter for missions with raw imagery
  raw_images_filepaths = listing[grepl("^[0-9]{6}/images/[0-9]{6}_images.zip$", listing)]
  missions_w_raw_images = str_split(raw_images_filepaths, "/") |> map_chr(1) |> unique()

  return(missions_w_raw_images)

}

# missions_to_process = get_objectstore_missions_w_raw_imagery()


# Option 3: By object store raw imagery listing, but only if all expected products for these missions do not
# already exist
get_objectstore_missions_not_yet_processed = function() {
  expected_missions = get_objectstore_missions_w_raw_imagery()

  # FULL FILES
  expected_processed_full_path = paste0(
    expected_missions, "/processed_", PHOTOGRAMMETRY_CONFIG_ID, "/full/"
  )

  expected_files = paste0(
    rep(expected_missions, each = length(EXPECTED_FILE_SUFFIXES_FULL)),
    EXPECTED_FILE_SUFFIXES_FULL
  )

  expected_processed_full_filepaths = paste0(
    rep(expected_processed_full_path, each = length(EXPECTED_FILE_SUFFIXES_FULL)),
    expected_files
  )

  # THUMBNAILS
  expected_processed_thumbnail_path = paste0(
    expected_missions, "/processed_", PHOTOGRAMMETRY_CONFIG_ID, "/thumbnails/"
  )

  expected_thumbnail_files = paste0(
    rep(expected_missions, each = length(EXPECTED_FILE_SUFFIXES_THUMBNAIL)),
    EXPECTED_FILE_SUFFIXES_THUMBNAIL
  )

  expected_processed_thumbnail_filepaths = paste0(
    rep(expected_processed_thumbnail_path, each = length(EXPECTED_FILE_SUFFIXES_THUMBNAIL)),
    expected_thumbnail_files
  )

  # COMBINED

  expected_processed_filepaths = c(
    expected_processed_full_filepaths,
    expected_processed_thumbnail_filepaths
  )

  # Get existing processed files
  remote_dir = paste0(RCLONE_REMOTE, ":", REMOTE_MISSIONS_DIR)
  command = paste("rclone lsf", remote_dir, "-R --files-only", sep = " ")
  existing_files = system(command, intern = TRUE)

  missing_files = setdiff(expected_processed_filepaths, existing_files) |> sort()

  missions_to_reprocess = str_split(missing_files, "/") |> map(1) |> unlist() |> unique() |> sort()

  # Remove "0000" since that is for debugging/prototypeing
  missions_to_reprocess = missions_to_reprocess[missions_to_reprocess != "000000"]

  return(missions_to_reprocess)
}

missions_to_process = get_objectstore_missions_not_yet_processed()





# Convert to data frame
missions_to_process_df = data.frame(mission_id = missions_to_process)


# Write
write_csv(missions_to_process_df, MISSIONS_TO_PROCESS_PHOTOGRAMMETRY_LIST_PATH)

## As an alternative/addition, break into chunks and write to a separate file per chunk
create_dir(file.path(MISSIONS_TO_PROCESS_PHOTOGRAMMETRY_PATH, "chunks"))

chunks = chunk_up_collate(missions_to_process_df$mission_id, N_CHUNKS_PHOTOGRAMMETRY)

for (i in seq_along(chunks)) {
  chunk = chunks[[i]]
  chunk_df = data.frame(mission_id = chunk)
  chunk_filename = paste0("missions-to-process_", i, ".csv")
  chunk_filepath = file.path(MISSIONS_TO_PROCESS_PHOTOGRAMMETRY_PATH, "chunks", chunk_filename)
  write_csv(chunk_df, chunk_filepath)
}
