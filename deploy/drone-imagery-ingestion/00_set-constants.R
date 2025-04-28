# Purpose: set global constants for the imagery data processing pipeline. They mostly consist of
# paths to files/dirs for input/output data.

library(tidyverse)

# What is the padding width for the dataset ID in the folder names of the imagery folders to be ingested? (New format) This is used to
# force the Baserow dataset ID column to conform to the image folder names, so this should reflect
# the padding used when naming the image folders in the 1_manually-cleaned folder.
FOLDER_DATASET_ID_PADDING = 6


IMAGE_MERGE_DISTANCE = 50
# If processing NRS datasets, use a larger merge distance because some datasets used very low overlap
IMAGE_MERGE_DISTANCE = 100

# For example image selection & image set zipping
N_EXAMPLE_IMAGES = 4
THUMBNAIL_SIZE = "800"
SKIP_EXISTING = FALSE # Skip processing for missions that already have all outputs


# Handle difference in how the current directory is set between debugging and command line call
if (file.exists("deploy/drone-imagery-ingestion/imagery_project_name.txt")) {
  IMAGERY_PROJECT_NAME_FILE = "deploy/drone-imagery-ingestion/imagery_project_name.txt"
} else {
  IMAGERY_PROJECT_NAME_FILE = "imagery_project_name.txt"
}
IMAGERY_PROJECT_NAME = read_lines(IMAGERY_PROJECT_NAME_FILE)


CONTRIBUTED_IMAGERY_PATH = "/ofo-share/drone-imagery-organization/1_manually-cleaned"
RAW_EXIF_PATH = "/ofo-share/drone-imagery-organization/metadata/1_reconciling-contributions/1_raw-exif/"

CONTRIBUTED_METADATA_PATH = "/ofo-share/drone-imagery-organization/ancillary/baserow-snapshots"

CONTRIBUTED_TO_SORTED_MISSION_ID_CROSSWALK_PATH = "/ofo-share/drone-imagery-organization/metadata/1_reconciling-contributions/3_contributed-to-sorted-id-crosswalk/"

IMAGE_EXIF_W_SORTING_PLAN_PATH = "/ofo-share/drone-imagery-organization/metadata/1_reconciling-contributions/2_exif-w-sorting-plan/"

# Should be "conributed"
EXTRACTED_METADATA_PER_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/2_intermediate/1_contributed-metadata-per-mission/"
EXTRACTED_METADATA_PER_SUB_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/2_intermediate/2_contributed-metadata-per-sub-mission/"

# # Optionally can specify a subset of missions to process (for the mission-level processing)
IMAGERY_PROJECT_SUBSET_MISSIONS = NULL
# IMAGERY_PROJECT_SUBSET_MISSIONS = c(000643:000900) |> str_pad(6, pad = "0", side = "left")

MISSIONS_TO_PROCESS_LIST_PATH = file.path("sandbox", "drone-imagery-ingestion", "missions-to-process.csv")

PARSED_EXIF_METADATA_PATH = "/ofo-share/drone-imagery-organization/metadata/2_intermediate/4_parsed-exif"

PARSED_EXIF_FOR_RETAINED_IMAGES_PATH = "/ofo-share/drone-imagery-organization/metadata/3_final/3_parsed-exif-per-image"
DERIVED_METADATA_PER_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/2_intermediate/6_derived-metadata-per-mission"
DERIVED_METADATA_PER_SUB_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/2_intermediate/7_derived-metadata-per-sub-mission"

FULL_METADATA_PER_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/3_final/1_full-metadata-per-mission/"
FULL_METADATA_PER_SUB_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/3_final/2_full-metadata-per-sub-mission/"

SORTED_IMAGERY_PATH = "/ofo-share/drone-imagery-organization/2_sorted"

IMAGERY_ZIP_AND_EXAMPLES_PATH = "/ofo-share/drone-imagery-organization/4_raw-imagery-zip-and-examples"

IN_PROCESS_PATH = "/ofo-share/tmp/raw-imagery-publish-prep-progress-tracking/"

UPLOAD_STAGING_DIR_PATH = "/ofo-share/drone-imagery-organization/9_temp-upload-staging"

# Remote patih in the CyVerse data store where all mission data is stored. Must end with a trailing
# slash.
CYVERSE_MISSIONS_DIR = paste0("/iplant/home/shared/ofo/public/missions_02/")
