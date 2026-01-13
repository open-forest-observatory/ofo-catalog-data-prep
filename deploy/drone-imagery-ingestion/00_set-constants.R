# Purpose: set global constants for the imagery data processing pipeline. They mostly consist of
# paths to files/dirs for input/output data.

# nolint start

library(tidyverse)
library(unixtools) # To set tempdir during session

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


PROJECT_TO_PROCESS_RAW_IMAGERY_METADATA_FILEPATH = "/ofo-share/project-data/catalog-data-prep/00_missions-to-process/01_raw-imagery-meatadata-prep/project-to-process.txt"
PROJECT_NAME_TO_PROCESS_RAW_IMAGERY_METADATA = tryCatch(
  read_lines(PROJECT_TO_PROCESS_RAW_IMAGERY_METADATA_FILEPATH),
  error = function(e) {
    warning("Could not read file: ", PROJECT_TO_PROCESS_RAW_IMAGERY_METADATA_FILEPATH, "\n", e$message)
    character(0)
  }
)
MISSIONS_TO_PROCESS_RAW_IMAGERY_METADATA_LIST_PATH = "/ofo-share/project-data/catalog-data-prep/00_missions-to-process/01_raw-imagery-meatadata-prep/missions-to-process.csv"

CONTRIBUTED_IMAGERY_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/1_manually-cleaned"
RAW_EXIF_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/1_reconciling-contributions/1_raw-exif/"

CONTRIBUTED_METADATA_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/ancillary/baserow-snapshots"

CONTRIBUTED_TO_SORTED_MISSION_ID_CROSSWALK_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/1_reconciling-contributions/3_contributed-to-sorted-id-crosswalk/"

IMAGE_EXIF_W_SORTING_PLAN_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/1_reconciling-contributions/2_exif-w-sorting-plan/"

# Should be "conributed"
EXTRACTED_METADATA_PER_MISSION_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/2_intermediate/1_contributed-metadata-per-mission/"
EXTRACTED_METADATA_PER_SUB_MISSION_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/2_intermediate/2_contributed-metadata-per-sub-mission/"

PARSED_EXIF_METADATA_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/2_intermediate/4_parsed-exif"

PARSED_EXIF_FOR_RETAINED_IMAGES_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/3_final/3_parsed-exif-per-image"
DERIVED_METADATA_PER_MISSION_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/2_intermediate/6_derived-metadata-per-mission"
DERIVED_METADATA_PER_SUB_MISSION_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/2_intermediate/7_derived-metadata-per-sub-mission"

FULL_METADATA_PER_MISSION_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/3_final/1_full-metadata-per-mission/"
FULL_METADATA_PER_SUB_MISSION_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/3_final/2_full-metadata-per-sub-mission/"

# The following file should not be assumed to be there, since we are not necessarily keeping all
# mission data combined locally, and it should not be assumed to contain every mission, but this is where the file will be if it has been created
FULL_METADATA_PER_MISSION_COMBINED_FILEPATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/3_final/ofo-all-missions-metadata.gpkg"
FULL_METADATA_PER_IMAGE_COMBINED_FILEPATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/metadata/3_final/ofo-all-images-metadata.gpkg"

PROJECT_TO_PROCESS_RAW_IMAGERY_FILES_FILEPATH = "/ofo-share/project-data/catalog-data-prep/00_missions-to-process/02_raw-imagery-file-prep/projects-to-process.txt"
PROJECT_NAMES_TO_PROCESS_RAW_IMAGERY_FILES = tryCatch(
  read_lines(PROJECT_TO_PROCESS_RAW_IMAGERY_FILES_FILEPATH),
  error = function(e) {
    warning("Could not read file: ", PROJECT_TO_PROCESS_RAW_IMAGERY_FILES_FILEPATH, "\n", e$message)
    character(0)
  }
)
# Exclude "commented out" project names
PROJECT_NAMES_TO_PROCESS_RAW_IMAGERY_FILES = PROJECT_NAMES_TO_PROCESS_RAW_IMAGERY_FILES[!grepl("^#", PROJECT_NAMES_TO_PROCESS_RAW_IMAGERY_FILES)]

MISSIONS_TO_PROCESS_RAW_IMAGERY_FILES_LIST_PATH = "/ofo-share/project-data/catalog-data-prep/00_missions-to-process/02_raw-imagery-file-prep/missions-to-process.csv"



SORTED_IMAGERY_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/2_sorted"

IMAGERY_ZIP_AND_EXAMPLES_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/4_raw-imagery-zip-and-examples"

IN_PROCESS_PATH = "/ofo-share/tmp/raw-imagery-publish-prep-progress-tracking/"

UPLOAD_STAGING_DIR_PATH = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/9_temp-upload-staging"

# Name of configured rclone remote for object storage
RCLONE_REMOTE = "js2s3"
REMOTE_MISSIONS_DIR = "/ofo-public/drone/missions_02/"
REMOTE_PHOTOGRAMMETRY_DIR = "/ofo-internal/photogrammetry-outputs/"

UPLOAD_ERROR_LOG = "/ofo-share/project-data/catalog-data-prep/01_raw-imagery-ingestion/temp/cyverse-upload-log.txt"


PROJECTS_TO_PROCESS_PHOTOGRAMMETRY_FILEPATH = "/ofo-share/project-data/catalog-data-prep/00_missions-to-process/03_photogrammetry/projects-to-process.txt"
PROJECT_NAMES_TO_PROCESS_PHOTOGRAMMETRY = tryCatch(
  read_lines(PROJECTS_TO_PROCESS_PHOTOGRAMMETRY_FILEPATH),
  error = function(e) {
    warning("Could not read file: ", PROJECTS_TO_PROCESS_PHOTOGRAMMETRY_FILEPATH, "\n", e$message)
    character(0)
  }
)
# Exclude "commented out" project names
PROJECT_NAMES_TO_PROCESS_PHOTOGRAMMETRY = PROJECT_NAMES_TO_PROCESS_PHOTOGRAMMETRY[!grepl("^#", PROJECT_NAMES_TO_PROCESS_PHOTOGRAMMETRY)]

MISSIONS_TO_PROCESS_PHOTOGRAMMETRY_LIST_PATH = "/ofo-share/project-data/catalog-data-prep/00_missions-to-process/03_photogrammetry/missions-to-process.csv"
MISSIONS_TO_PROCESS_PHOTOGRAMMETRY_PATH = "/ofo-share/project-data/catalog-data-prep/00_missions-to-process/03_photogrammetry"




# OFO_R_REPO_PATH = "/ofo-share/repos-derek/ofo-r"

AUTOMATE_METASHAPE_PATH = "/ofo-share/repos-derek/automate-metashape"

PHOTOGRAMMETRY_DIR = "/ofo-share/project-data/catalog-data-prep/02_photogrammetry"

# All of these are relative to PHOTOGRAMMETRY_DIR:
BASE_METASHAPE_CONFIG_SUBPATH = "01_base-metashape-config"
DOWNLOADED_IMAGERY_ZIP_SUBDIR = "02_downloaded-imagery-zip"
INPUT_IMAGES_SUBDIR = "03_input-images"
DERIVED_METASHAPE_CONFIG_SUBDIR = "04_derived-metashape-configs"
METASHAPE_PROJECT_SUBDIR = "05_photogrammetry-projects"
METASHAPE_OUTPUT_SUBDIR = "06_photogrammetry-outputs"
# After the Metashape processing is done, the outputs are uploaded, deleted from the above
# folder, then downloaded to the folder below
METASHAPE_OUTPUT_DOWNLOADED_SUBDIR = "06_photogrammetry-outputs-downloaded"
PHOTOGRAMMETRY_POSTPROCESSED_SUBDIR = "07_photogrammetry-outputs-postprocessed"

# The photogrammetry base config file in BASE_METASHAPE_CONFIG_SUBPATH should have the following ID
# as its filename (with a .yml extension)
PHOTOGRAMMETRY_CONFIG_ID = "02"

# The number of chunks to break the photogrammetry processing into (one chunk for each instance)
N_CHUNKS_PHOTOGRAMMETRY = 24



# For photogrammetry post-processing

# For COPC conversion, the current approach of using Untwine is super slow for large point clouds
# due to the extensive rapid file caching and the files being stored on the NFS server. Therefore we
# can use Metashape for the conversion and disable here.
CONVERT_TO_COPC = FALSE

# What fraction of the system RAM can TERRA use. The terra default is 0.6, you cannot go above 0.9
# without a warning.
TERRA_MEMFRAC = 0.9
OUTPUT_MAX_DIM = 800

TEMPDIR = "/ofo-share/tmp/data-catalog-prep/"

## Some setup that should apply to any scripts we run

terra::terraOptions(memfrac = TERRA_MEMFRAC)

# Set the tempdir
unixtools::set.tempdir(TEMPDIR)
cat("Tempdir is:", tempdir(), "\n")


# For ITD

CHMS_FOR_ITD_LIST_PATH = "/ofo-share/project-data/catalog-data-prep/00_missions-to-process/04_itd/chms-for-itd.csv"

ITD_PARAMETERIZATION_ID = 1
CHM_RES = 0.25
CHM_SMOOTH_WIDTH = 7
LMF_A = 0
LMF_B = 0.11
LMF_C = 0
LMF_DIAM_MIN = 0.5
LMF_DIAM_MAX = 100

# The folder beneath the mission processed folder where the ITD outputs will be stored (in the
# Object Store))
ITD_FOLDER = paste0("itd_", str_pad(ITD_PARAMETERIZATION_ID, width = 4, pad = "0", side = "left"))




### Drone imagery web catalog

BASE_OFO_URL = "https://openforestobservatory.org/"
# BASE_OFO_URL = "http://localhost:1313/"
WEBSITE_REPO_PATH = "/ofo-share/repos-derek/ofo-website-3"

# Path to the plot details template page within this repo
MISSION_DETAILS_TEMPLATE_FILEPATH = fs::path("deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/templates/drone-mission-details.md")

# Path to plot details dir relative to the 'content' dir in the website repo. No leading slash but
# trailing slash
MISSION_DETAILS_PAGE_DIR = "data/drone/mission-details/"

# Path to static website files dirs within website repo relative to 'static'. Leading slash but no
# trailing slash
DATATABLE_HEADER_FILES_DIR = "/datatable-header-files_drone-mission-catalog"
LEAFLET_HEADER_FILES_DIR = "/leaflet-header-files_drone-mission-catalog"
MISSION_CATALOG_DATATABLE_DIR = "/drone-mission-catalog-datatable/"
MISSION_CATALOG_DATATABLE_FILENAME = "drone-mission-catalog-datatable.html"
MISSION_CATALOG_MAP_DIR = "/drone-mission-catalog-map/"
MISSION_CATALOG_MAP_FILENAME = "drone-mission-catalog-map.html"
MISSION_DETAILS_DATATABLE_DIR = "/drone-mission-details-datatables"
MISSION_DETAILS_MAP_DIR = "/drone-mission-details-maps"
ITD_MAP_DIR = "/itd-maps"

WEBSITE_STATIC_PATH = file.path(WEBSITE_REPO_PATH, "static", "")
WEBSITE_CONTENT_PATH = file.path(WEBSITE_REPO_PATH, "content", "")

DATA_SERVER_BASE_URL = "https://js2.jetstream-cloud.org:8001/swift/v1/"
DATA_SERVER_MISSIONS_BASE_URL = paste0(DATA_SERVER_BASE_URL, REMOTE_MISSIONS_DIR)

# Metadata pulled from S3 to use in catalog generation
MISSION_METADATA_FILEPATH = "/ofo-share/project-data/catalog-data-prep/05_drone-imagery-web-catalog/01_metadata/mission-metadata.gpkg"
IMAGE_METADATA_FILEPATH = "/ofo-share/project-data/catalog-data-prep/05_drone-imagery-web-catalog/01_metadata/image-metadata.gpkg"
S3_LISTING_FILEPATH = "/ofo-share/project-data/catalog-data-prep/05_drone-imagery-web-catalog/01_metadata/s3-listing.csv"

# Lists of mission and plot IDs to withhold from broad-scale ML training (save for testing)
MISSIONS_TO_WITHHOLD_FILEPATH = "/ofo-share/project-data/catalog-data-prep/stratification-data/strat-output/withheld_drone_mission_ids_v1.csv"

# nolint end