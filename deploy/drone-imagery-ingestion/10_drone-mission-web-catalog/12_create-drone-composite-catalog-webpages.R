# Purpose: Create the Hugo markdown pages for all composite drone mission datasets including
# visualizations/links to the actual data on S3, as well as a dataset index page.
# Composite missions are pairs of overlapping drone flights at different altitudes that are
# processed together for unified data products.

library(dplyr)
library(leaflet)
library(stringr)
library(sf)
library(htmlwidgets)
library(DT)
library(jinjar)
library(readr)
library(purrr)


source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")
source("src/web-catalog-creation_shared-functions.R")
source("src/web-catalog-creation_drone-imagery-catalog.R")



# Load the composite mission polygons and points metadata
composite_polygons_w_metadata = st_read(COMPOSITE_MISSION_METADATA_FILEPATH)
composite_points = st_read(COMPOSITE_IMAGE_METADATA_FILEPATH)

# Save header library files required by embedded HTML datatables and leaflet maps
save_dt_header_files(WEBSITE_STATIC_PATH, COMPOSITE_DATATABLE_HEADER_FILES_DIR)
save_leaflet_header_files(WEBSITE_STATIC_PATH, COMPOSITE_LEAFLET_HEADER_FILES_DIR)

# Compile relevant and human-readable values from composite mission attributes
composite_polygons_w_summary_data = compile_composite_summary_data(
  composite_mission_metadata = composite_polygons_w_metadata,
  base_ofo_url = BASE_OFO_URL,
  composite_details_dir = COMPOSITE_MISSION_DETAILS_PAGE_DIR,
  individual_mission_details_dir = MISSION_DETAILS_PAGE_DIR
)

# Make a HTML data table of composite catalog
dt = make_composite_catalog_datatable(
  composite_summary = composite_polygons_w_summary_data,
  website_static_path = WEBSITE_STATIC_PATH,
  datatable_header_files_dir = COMPOSITE_DATATABLE_HEADER_FILES_DIR,
  composite_catalog_datatable_dir = COMPOSITE_MISSION_CATALOG_DATATABLE_DIR,
  composite_catalog_datatable_filename = COMPOSITE_MISSION_CATALOG_DATATABLE_FILENAME
)

# Make leaflet map of composite catalog
m = make_composite_catalog_map(
  composite_summary = composite_polygons_w_summary_data,
  website_static_path = WEBSITE_STATIC_PATH,
  leaflet_header_files_dir = COMPOSITE_LEAFLET_HEADER_FILES_DIR,
  composite_catalog_map_dir = COMPOSITE_MISSION_CATALOG_MAP_DIR,
  composite_catalog_map_filename = COMPOSITE_MISSION_CATALOG_MAP_FILENAME
)

# Copy the listing page template to the website content directory
listing_template_source = COMPOSITE_MISSION_CATALOG_TEMPLATE_FILEPATH
listing_template_dest = file.path(WEBSITE_CONTENT_PATH, "data/drone-composites/_index.html")
dir.create(dirname(listing_template_dest), recursive = TRUE, showWarnings = FALSE)
file.copy(listing_template_source, listing_template_dest, overwrite = TRUE)


############################# Detail pages

s3_file_listing = read_csv(COMPOSITE_S3_LISTING_FILEPATH)

composite_summary = composite_polygons_w_summary_data |> dplyr::arrange(composite_id)
composite_ids = unique(composite_summary$composite_id)

# Make composite details pages using sequential walk() to avoid HTML widget staging directory
# conflicts that can occur with parallel processing
walk(
  composite_ids,
  make_composite_details_page,
  all_composite_ids = composite_ids,
  composite_summaries = composite_summary,
  composite_points = composite_points,
  s3_file_listing = s3_file_listing,
  website_static_path = WEBSITE_STATIC_PATH,
  website_content_path = WEBSITE_CONTENT_PATH,
  leaflet_header_files_dir = COMPOSITE_LEAFLET_HEADER_FILES_DIR,
  datatable_header_files_dir = COMPOSITE_DATATABLE_HEADER_FILES_DIR,
  composite_details_datatable_dir = COMPOSITE_MISSION_DETAILS_DATATABLE_DIR,
  composite_details_map_dir = COMPOSITE_MISSION_DETAILS_MAP_DIR,
  composite_itd_map_dir = COMPOSITE_ITD_MAP_DIR,
  composite_details_template_filepath = COMPOSITE_MISSION_DETAILS_TEMPLATE_FILEPATH,
  composite_details_page_dir = COMPOSITE_MISSION_DETAILS_PAGE_DIR,
  .progress = TRUE
)
