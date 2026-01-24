# deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/02_summarize-curated-metadata-per-mission.R
# Purpose: Re-summarize EXIF metadata at mission and sub-mission level after curation filtering.
#
# This is equivalent to 01_raw-imagery-metadata-prep/06_summarize-exif-metadata-per-mission.R
# but operates on post-curation filtered images from script 01.
#
# Workflow:
# 1. Read post-curation filtered image metadata (already excludes extraneous images)
# 2. Compute polygons at mission and sub-mission level
# 3. Filter to images retained in polygon intersection
# 4. Re-compute polygons with final retained images
# 5. Extract summary statistics and write attributed polygons

library(tidyverse)
library(sf)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/summarization-utils.R")

# ============================================================================
# Determine which missions to process
# ============================================================================

post_curation_files = list.files(
  POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
  pattern = "\\.gpkg$",
  full.names = FALSE
)

missions_to_process = str_extract(post_curation_files, "^\\d{6}")

cat(sprintf("Processing %d missions for post-curation summarization...\n", length(missions_to_process)))

# ============================================================================
# Create output directories
# ============================================================================

create_dir(POST_CURATION_DERIVED_METADATA_PER_MISSION_PATH)
create_dir(POST_CURATION_DERIVED_METADATA_PER_SUB_MISSION_PATH)

# ============================================================================
# Run summarization for each mission using shared function
# ============================================================================

# Note: Script 01 filters images based on curation notes. This script further
# filters based on polygon retention (images outside computed polygons are excluded).
# We OVERWRITE the files in POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH
# with the doubly-filtered images (curation + polygon filtering).

future::plan(multisession)
results = future_walk(
  missions_to_process,
  ~ summarize_mission_exif(
    mission_id_foc = .x,
    input_metadata_path = POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
    output_derived_mission_path = POST_CURATION_DERIVED_METADATA_PER_MISSION_PATH,
    output_derived_sub_mission_path = POST_CURATION_DERIVED_METADATA_PER_SUB_MISSION_PATH,
    output_retained_images_path = POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,  # Overwrite with doubly-filtered images
    image_merge_distance = IMAGE_MERGE_DISTANCE
  ),
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)

cat("\n**** Post-curation summarization complete ****\n")
