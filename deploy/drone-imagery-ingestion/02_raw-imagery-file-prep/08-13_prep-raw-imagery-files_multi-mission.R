# Purpose: Run the imagery prep steps 08 through 14 for a set of specified mission IDs.

# IMPORTANT NOTE: You must have already configure an rclone remote on this machine for the object
# store. See the upload task source script for details.

# NOTE: Requires 'exiftool' be installed as an ubuntu package

library(tidyverse)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")

# Load the function that calls steps 08 through 14 for a specified mission ID
source("deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/08-13_prep-raw-imagery-files_per-mission.R")

# Determine the missions to process. This depends on having set the project name in the file
# "/ofo-share/catalog-data-prep/00_missions-to-process/02_raw-imagery-file-prep/projects-to-process.txt".
# It puts the list of missions to process in the file "missions-to-process.csv". Note that the
# script that is sourced here does not define a function, it actually runds the code complete the
# task described in this comment.
source("deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/src/08_0pre_determine-missions-to-process.R")
missions_to_process = read_csv(MISSIONS_TO_PROCESS_RAW_IMAGERY_FILES_LIST_PATH) |> pull(mission_id)

# Use (and update) this if you want to override the missions to process determined above, such as if
# the final check revealed some missions were completed all the way through to cyverse upload
# missions_to_process = c("000162", "000186")

future::plan(multicore, workers = future::availableCores()*2)

furrr::future_map(
  missions_to_process,
  prep_raw_imagery_files_per_mission,
  .progress = TRUE
)


# ######## Check uploaded files to ensure all are present and return a vector of missing mission IDs
# # for re-processing
# ########

# # Example images
# coll_name = shQuote(paste0(CYVERSE_MISSIONS_DIR, "______/images/examples/%"))
# data_name = shQuote("example\\__.JPG")
# call = paste0("iquest --no-page '%s/%s' \"select COLL_NAME, DATA_NAME where COLL_NAME like ", coll_name, " and DATA_NAME like ", data_name, "\"")
# example_images = system(call, intern = TRUE)

# # Image zip
# coll_name = shQuote(paste0(CYVERSE_MISSIONS_DIR, "______/images"))
# data_name = shQuote("______\\_images.zip")
# call = paste0("iquest --no-page '%s/%s' \"select COLL_NAME, DATA_NAME where COLL_NAME like ", coll_name, " and DATA_NAME like ", data_name, "\"")
# zip_images = system(call, intern = TRUE)

# # Metadata images
# coll_name = shQuote(paste0(CYVERSE_MISSIONS_DIR, "______/metadata-images"))
# data_name = shQuote("______\\_image-metadata.gpkg")
# call = paste0("iquest --no-page '%s/%s' \"select COLL_NAME, DATA_NAME where COLL_NAME like ", coll_name, " and DATA_NAME like ", data_name, "\"")
# metadata_images = system(call, intern = TRUE)

# # Metadata missions
# coll_name = shQuote(paste0(CYVERSE_MISSIONS_DIR, "______/metadata-mission"))
# data_name = shQuote("______\\_mission-metadata.gpkg")
# call = paste0("iquest --no-page '%s/%s' \"select COLL_NAME, DATA_NAME where COLL_NAME like ", coll_name, " and DATA_NAME like ", data_name, "\"")
# metadata_missions = system(call, intern = TRUE)

# # For each of these, get the mission ID as the 6 digit number in the file path
# mission_id_example_images = str_extract(example_images, "/[0-9]{6}/") |> str_remove_all("/")
# mission_id_zip_images = str_extract(zip_images, "/[0-9]{6}/") |> str_remove_all("/")
# mission_id_metadata_images = str_extract(metadata_images, "/[0-9]{6}/") |> str_remove_all("/")
# mission_id_metadata_missions = str_extract(metadata_missions, "/[0-9]{6}/") |> str_remove_all("/")

# df_of_mission_id_count = function(mission_ids, count_col_name) {
#   table = table(mission_ids) |> as.data.frame(stringsAsFactors = FALSE)
#   colnames(table) = c("mission_id", count_col_name)
#   return(table)
# }

# # There should be 1 file for each mission ID, except 8 for the example images. So, summarize each
# # into a table
# mission_id_table_example_images = df_of_mission_id_count(mission_id_example_images, count_col_name = "example_images")
# mission_id_table_zip_images = df_of_mission_id_count(mission_id_zip_images, count_col_name = "zip_images")
# mission_id_table_metadata_images = df_of_mission_id_count(mission_id_metadata_images, count_col_name = "metadata_images")
# mission_id_table_metadata_missions = df_of_mission_id_count(mission_id_metadata_missions, count_col_name = "metadata_missions")

# # Comine the tables into one with a series of left_joins
# mission_id_counts = purrr::reduce(
#   list(
#     mission_id_table_example_images,
#     mission_id_table_zip_images,
#     mission_id_table_metadata_images,
#     mission_id_table_metadata_missions
#   ),
#   dplyr::left_join,
#   by = 'mission_id'
# )

# # Check the counts
# partially_missing_mission_ids = mission_id_counts |>
#   filter(example_images != 8 | zip_images != 1 | metadata_images != 1 | metadata_missions != 1) |>
#   pull(mission_id)

# # Make sure every mission is represented
# fully_missing_mission_ids = setdiff(missions_to_process, mission_id_counts$mission_id)

# missing_mission_ids = c(partially_missing_mission_ids, fully_missing_mission_ids)

# if (length(missing_mission_ids) > 0) {
#   cat("Some missions failed to upload some or all files. Here's a dput of the missing mission IDs so they can be re-run.\n")
#   dput(missing_mission_ids)
# }