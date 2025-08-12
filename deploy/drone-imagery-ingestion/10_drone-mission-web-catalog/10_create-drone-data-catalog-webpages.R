# Purpose: Create the Hugo markdown pages for all drone mission datasets including
# visualizations/links to the actual data on CyVerse, as well as a dataset index
# page. NOTE that if we ever want to delete the whole directory of mission deatail pages and start
# fresh, we have to populate it with a _index.md file so that Hugo will still recognize the directory
# as a section of the website.

library(dplyr)
library(leaflet)
library(stringr)
library(sf)
library(htmlwidgets)
library(DT)
library(jinjar)
library(readr)
library(furrr)


source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")
source("src/web-catalog-creation_shared-functions.R")
source("src/web-catalog-creation_drone-imagery-catalog.R")



# Load the mission polygons and points metadata
mission_polygons_w_metadata = st_read(MISSION_METADATA_FILEPATH)
mission_points = st_read(IMAGE_METADATA_FILEPATH)

# Save header library files required by embedded HTML datatables and leaflet maps
save_dt_header_files(WEBSITE_STATIC_PATH, DATATABLE_HEADER_FILES_DIR)
save_leaflet_header_files(WEBSITE_STATIC_PATH, LEAFLET_HEADER_FILES_DIR)

# # ** Optionally keep only nadir missions (or keep only missions with photogrammetry products, etc)
# mission_polygons_w_metadata = mission_polygons_w_metadata |>
#   filter(abs(camera_pitch_derived) < 10)

# Compile relevant and human-readable values from mission attributes as additional columns of the
# mission polygons object
mission_polygons_w_summary_data = compile_mission_summary_data(
  mission_level_metadata = mission_polygons_w_metadata,
  base_ofo_url = BASE_OFO_URL,
  mission_details_dir = MISSION_DETAILS_PAGE_DIR
)

# Make a HTML data table of plot catalog
dt = make_mission_catalog_datatable(
  mission_summary = mission_polygons_w_summary_data,
  website_static_path = WEBSITE_STATIC_PATH,
  datatable_header_files_dir = DATATABLE_HEADER_FILES_DIR,
  mission_catalog_datatable_dir = MISSION_CATALOG_DATATABLE_DIR,
  mission_catalog_datatable_filename = MISSION_CATALOG_DATATABLE_FILENAME
)

# Make leaflet map of field data catalog
m = make_mission_catalog_map(
  mission_summary = mission_polygons_w_summary_data,
  website_static_path = WEBSITE_STATIC_PATH,
  leaflet_header_files_dir = LEAFLET_HEADER_FILES_DIR,
  mission_catalog_map_dir = MISSION_CATALOG_MAP_DIR,
  mission_catalog_map_filename = MISSION_CATALOG_MAP_FILENAME
)


############################# Detail pages

s3_file_listing = read_csv(S3_LISTING_FILEPATH)

mission_summary = mission_polygons_w_summary_data |> dplyr::arrange(mission_id)
mission_ids = mission_summary$mission_id

# Make mission details pages. Was unable to parallelize this -- possibly because of a conflict with
# creating HTML widgets in parallel (saw an error about a widget staging directory that was
# missing). It seems that it might work to set a temp dir separately for each process using the
# mission ID which should prevent conflicts. Like so:
  # temp_dir <- file.path(tempdir(), paste0("proc_", process_id))
  # dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # # Force R to use this temp dir
  # Sys.setenv(TMPDIR = temp_dir)


walk(
  mission_ids,
  make_mission_details_page,
  all_mission_ids = mission_ids,
  mission_summaries = mission_summary,
  mission_points = mission_points,
  s3_file_listing = s3_file_listing,
  website_static_path = WEBSITE_STATIC_PATH,
  website_content_path = WEBSITE_CONTENT_PATH,
  leaflet_header_files_dir = LEAFLET_HEADER_FILES_DIR,
  datatable_header_files_dir = DATATABLE_HEADER_FILES_DIR,
  mission_details_datatable_dir = MISSION_DETAILS_DATATABLE_DIR,
  mission_details_map_dir = MISSION_DETAILS_MAP_DIR,
  itd_map_dir = ITD_MAP_DIR,
  mission_details_template_filepath = MISSION_DETAILS_TEMPLATE_FILEPATH,
  mission_details_page_dir = MISSION_DETAILS_PAGE_DIR
)


# # Example of how to call the function to make a single mission page
# make_mission_details_page(
#   x = "001441",
#   y = mission_points_list[[1]], # This is the list of image points for that mission
#   all_mission_ids = mission_ids,
#   mission_summary = mission_summary,
#   s3_file_listing = s3_file_listing,
#   website_static_path = WEBSITE_STATIC_PATH,
#   website_content_path = WEBSITE_CONTENT_PATH,
#   leaflet_header_files_dir = LEAFLET_HEADER_FILES_DIR,
#   datatable_header_files_dir = DATATABLE_HEADER_FILES_DIR,
#   mission_details_datatable_dir = MISSION_DETAILS_DATATABLE_DIR,
#   mission_details_map_dir = MISSION_DETAILS_MAP_DIR,
#   itd_map_dir = ITD_MAP_DIR,
#   mission_details_template_filepath = MISSION_DETAILS_TEMPLATE_FILEPATH,
#   mission_details_page_dir = MISSION_DETAILS_PAGE_DIR
# )
