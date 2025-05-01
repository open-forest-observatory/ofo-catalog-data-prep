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
if (file.exists("deploy/drone-imagery-ingestion/01_raw-imagery-metadata-prep/project-to-process.txt")) {
  PROJECT_TO_PROCESS_RAW_IMAGERY_METADATA_FILEPATH = "deploy/drone-imagery-ingestion/01_raw-imagery-metadata-prep/project-to-process.txt"
} else {
  PROJECT_TO_PROCESS_RAW_IMAGERY_METADATA_FILEPATH = "project-to-process.txt"
}
PROJECT_NAME_TO_PROCESS_RAW_IMAGERY_METADATA = read_lines(PROJECT_TO_PROCESS_RAW_IMAGERY_METADATA_FILEPATH)

MISSIONS_TO_PROCESS_RAW_IMAGERY_METADATA_LIST_PATH = file.path("deploy", "drone-imagery-ingestion", "01_raw-imagery-metadata-prep", "missions-to-process.csv")


CONTRIBUTED_IMAGERY_PATH = "/ofo-share/drone-imagery-organization/1_manually-cleaned"
RAW_EXIF_PATH = "/ofo-share/drone-imagery-organization/metadata/1_reconciling-contributions/1_raw-exif/"

CONTRIBUTED_METADATA_PATH = "/ofo-share/drone-imagery-organization/ancillary/baserow-snapshots"

CONTRIBUTED_TO_SORTED_MISSION_ID_CROSSWALK_PATH = "/ofo-share/drone-imagery-organization/metadata/1_reconciling-contributions/3_contributed-to-sorted-id-crosswalk/"

IMAGE_EXIF_W_SORTING_PLAN_PATH = "/ofo-share/drone-imagery-organization/metadata/1_reconciling-contributions/2_exif-w-sorting-plan/"

# Should be "conributed"
EXTRACTED_METADATA_PER_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/2_intermediate/1_contributed-metadata-per-mission/"
EXTRACTED_METADATA_PER_SUB_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/2_intermediate/2_contributed-metadata-per-sub-mission/"

PARSED_EXIF_METADATA_PATH = "/ofo-share/drone-imagery-organization/metadata/2_intermediate/4_parsed-exif"

PARSED_EXIF_FOR_RETAINED_IMAGES_PATH = "/ofo-share/drone-imagery-organization/metadata/3_final/3_parsed-exif-per-image"
DERIVED_METADATA_PER_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/2_intermediate/6_derived-metadata-per-mission"
DERIVED_METADATA_PER_SUB_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/2_intermediate/7_derived-metadata-per-sub-mission"

FULL_METADATA_PER_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/3_final/1_full-metadata-per-mission/"
FULL_METADATA_PER_SUB_MISSION_PATH = "/ofo-share/drone-imagery-organization/metadata/3_final/2_full-metadata-per-sub-mission/"




# Handle difference in how the current directory is set between debugging and command line call
if (file.exists("deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/projects-to-process.txt")) {
  PROJECTS_TO_PROCESS_RAW_IMAGERY_FILES_FILEPATH = "deploy/drone-imagery-ingestion/02_raw-imagery-file-prep/projects-to-process.txt"
} else {
  PROJECTS_TO_PROCESS_RAW_IMAGERY_FILES_FILEPATH = "projects-to-process.txt"
}
PROJECT_NAMES_TO_PROCESS_RAW_IMAGERY_FILES = read_lines(PROJECTS_TO_PROCESS_RAW_IMAGERY_FILES_FILEPATH)
# Remove any project names starting in "#" (commented out)
PROJECT_NAMES_TO_PROCESS_RAW_IMAGERY_FILES = PROJECT_NAMES_TO_PROCESS_RAW_IMAGERY_FILES[!grepl("^#", PROJECT_NAMES_TO_PROCESS_RAW_IMAGERY_FILES)]

MISSIONS_TO_PROCESS_RAW_IMAGERY_FILES_LIST_PATH = file.path("deploy", "drone-imagery-ingestion", "02_raw-imagery-file-prep", "missions-to-process.csv")






SORTED_IMAGERY_PATH = "/ofo-share/drone-imagery-organization/2_sorted"

IMAGERY_ZIP_AND_EXAMPLES_PATH = "/ofo-share/drone-imagery-organization/4_raw-imagery-zip-and-examples"

IN_PROCESS_PATH = "/ofo-share/tmp/raw-imagery-publish-prep-progress-tracking/"

UPLOAD_STAGING_DIR_PATH = "/ofo-share/drone-imagery-organization/9_temp-upload-staging"

# Remote patih in the CyVerse data store where all mission data is stored. Must end with a trailing
# slash.
CYVERSE_MISSIONS_DIR = paste0("/iplant/home/shared/ofo/public/missions_05/")

UPLOAD_ERROR_LOG = "/ofo-share/drone-imagery-organization/temp/cyverse-upload-log.txt"


# Handle difference in how the current directory is set between debugging and command line call
if (file.exists("deploy/drone-imagery-ingestion/03_photogrammetry/projects-to-process.txt")) {
  PROJECTS_TO_PROCESS_PHOTOGRAMMETRY_FILEPATH = "deploy/drone-imagery-ingestion/03_photogrammetry/projects-to-process.txt"
} else {
  PROJECTS_TO_PROCESS_PHOTOGRAMMETRY_FILEPATH = "projects-to-process.txt"
}
PROJECT_NAMES_TO_PROCESS_PHOTOGRAMMETRY = read_lines(PROJECTS_TO_PROCESS_PHOTOGRAMMETRY_FILEPATH)
# Remove any project names starting in "#" (commented out)
PROJECT_NAMES_TO_PROCESS_PHOTOGRAMMETRY = PROJECT_NAMES_TO_PROCESS_PHOTOGRAMMETRY[!grepl("^#", PROJECT_NAMES_TO_PROCESS_PHOTOGRAMMETRY)]

MISSIONS_TO_PROCESS_PHOTOGRAMMETRY_LIST_PATH = file.path("deploy", "drone-imagery-ingestion", "03_photogrammetry", "missions-to-process.csv")



# OFO_R_REPO_PATH = "/ofo-share/repos-derek/ofo-r"

AUTOMATE_METASHAPE_PATH = "/ofo-share/repos-derek/automate-metashape"

PHOTOGRAMMETRY_DIR = "/ofo-share/catalog-data-prep/02_photogrammetry"

# All of these are relative to PHOTOGRAMMETRY_DIR:
BASE_METASHAPE_CONFIG_FILE_SUBPATH = "01_base-metashape-config/base_01.yml"
DOWNLOADED_IMAGERY_ZIP_SUBDIR = "02_downloaded-imagery-zip"
INPUT_IMAGES_SUBDIR = "03_input-images"
DERIVED_METASHAPE_CONFIG_SUBDIR = "04_derived-metashape-configs/01"
METASHAPE_PROJECT_SUBDIR = "05_photogrammetry-projects"
METASHAPE_OUTPUT_SUBDIR = "06_photogrammetry-outputs"
PHOTOGRAMMETRY_POSTPROCESSED_SUBDIR = "07_photogrammetry-outputs-postprocessed"

# For photogrammetry post-processing

# What fraction of the system RAM can TERRA use. The terra default is 0.6, you cannot go above 0.9
# without a warning.
TERRA_MEMFRAC = 0.9
OUTPUT_MAX_DIM = 800