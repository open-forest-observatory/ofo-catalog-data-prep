# Purpose: Create internal curation website pages for drone mission datasets.
# This script generates simplified pages showing only mission-level metadata and image location maps
# for data curators to inspect datasets.
#
# Key features:
# - Generates pages from image- and mission-level metadata files (.gpkg format)
# - Supports secondary image metadata file for specific missions listed in override list
# - Outputs simplified pages with mission metadata table and leaflet map only
# - Reuses core functions from the public catalog creation code

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

# ============================================================================
# Load metadata
# ============================================================================

cat("Loading mission metadata...\n")
mission_polygons_w_metadata = st_read(MISSION_METADATA_FILEPATH, quiet = TRUE)

cat("Loading primary image metadata...\n")
primary_image_points = st_read(IMAGE_METADATA_FILEPATH, quiet = TRUE)

# Load post-curation metadata if available (for side-by-side comparison on curation pages)
post_curation_mission_metadata = NULL
post_curation_image_points = NULL
missions_with_post_curation_data = character(0)

if (file.exists(POST_CURATION_MISSION_METADATA_FILEPATH)) {
  cat("Loading post-curation mission metadata...\n")
  post_curation_mission_metadata = st_read(POST_CURATION_MISSION_METADATA_FILEPATH, quiet = TRUE)
  missions_with_post_curation_data = unique(post_curation_mission_metadata$mission_id)
  cat(sprintf("  Found post-curation data for %d missions\n", length(missions_with_post_curation_data)))
} else {
  cat("No post-curation mission metadata found at:", POST_CURATION_MISSION_METADATA_FILEPATH, "\n")
}

if (file.exists(POST_CURATION_IMAGE_METADATA_FILEPATH)) {
  cat("Loading post-curation image metadata...\n")
  post_curation_image_points = st_read(POST_CURATION_IMAGE_METADATA_FILEPATH, quiet = TRUE)
} else {
  cat("No post-curation image metadata found at:", POST_CURATION_IMAGE_METADATA_FILEPATH, "\n")
}

# Load secondary image metadata if specified (not empty string and file exists)
secondary_image_points = NULL
if (nzchar(SECONDARY_IMAGE_METADATA_FILEPATH) && file.exists(SECONDARY_IMAGE_METADATA_FILEPATH)) {
  cat("Loading secondary image metadata...\n")
  secondary_image_points = st_read(SECONDARY_IMAGE_METADATA_FILEPATH, quiet = TRUE)
}

# Load override list if specified (not empty string and file exists)
override_missions = character(0)
if (nzchar(IMAGERY_METADATA_MISSION_OVERRIDE_LIST_FILEPATH) &&
    file.exists(IMAGERY_METADATA_MISSION_OVERRIDE_LIST_FILEPATH)) {
  cat("Loading override list...\n")
  override_df = read_csv(IMAGERY_METADATA_MISSION_OVERRIDE_LIST_FILEPATH, show_col_types = FALSE)
  if ("mission_id" %in% names(override_df)) {
    override_missions = unique(override_df$mission_id)
    cat(sprintf("  Found %d missions in override list\n", length(override_missions)))
  } else {
    warning("Override list does not contain 'mission_id' column. Ignoring override list.")
  }
} else {
  cat("No override list specified. All missions will use primary image metadata.\n")
}

# ============================================================================
# Prepare combined image metadata
# ============================================================================

cat("Preparing image metadata...\n")

# Function to determine which image metadata to use for each mission
get_image_metadata_for_missions = function(mission_ids,
                                           primary_points,
                                           secondary_points = NULL,
                                           override_list = character(0)) {

  if (is.null(secondary_points) || length(override_list) == 0) {
    # No secondary metadata or override list, use primary for all
    return(primary_points)
  }

  # Filter missions that should use secondary metadata
  missions_to_override = mission_ids[mission_ids %in% override_list]
  missions_to_keep_primary = mission_ids[!(mission_ids %in% override_list)]

  cat(sprintf("  Using secondary metadata for %d missions\n", length(missions_to_override)))
  cat(sprintf("  Using primary metadata for %d missions\n", length(missions_to_keep_primary)))

  # Extract relevant points from each source
  primary_subset = primary_points |> filter(mission_id %in% missions_to_keep_primary)
  secondary_subset = secondary_points |> filter(mission_id %in% missions_to_override)


  # Convert all secondary columns to character to enable merge (was already done for primary apparently)
  secondary_subset = secondary_subset |> mutate(across(image_id:mission_id, as.character))

  # Combine
  combined_points = bind_rows(primary_subset, secondary_subset)

  return(combined_points)
}

# Get all mission IDs
all_mission_ids = unique(mission_polygons_w_metadata$mission_id)

# Get combined image metadata (primary + secondary where appropriate)
mission_points = get_image_metadata_for_missions(
  mission_ids = all_mission_ids,
  primary_points = primary_image_points,
  secondary_points = secondary_image_points,
  override_list = override_missions
)

cat(sprintf("Total image points loaded: %d\n", nrow(mission_points)))

# ============================================================================
# Save header library files
# ============================================================================

cat("Saving header library files...\n")
save_dt_header_files(WEBSITE_STATIC_PATH, CURATION_DATATABLE_HEADER_FILES_DIR)
save_leaflet_header_files(WEBSITE_STATIC_PATH, CURATION_LEAFLET_HEADER_FILES_DIR)

# ============================================================================
# Compile mission summary data
# ============================================================================

cat("Compiling mission summary data...\n")
mission_polygons_w_summary_data = compile_mission_summary_data(
  mission_level_metadata = mission_polygons_w_metadata,
  base_ofo_url = BASE_OFO_URL,
  mission_details_dir = CURATION_MISSION_DETAILS_PAGE_DIR
)

mission_summary = mission_polygons_w_summary_data |> dplyr::arrange(mission_id)
mission_ids = mission_summary$mission_id

cat(sprintf("Processing %d missions\n", length(mission_ids)))

# ============================================================================
# Create mission catalog overview (map and datatable)
# ============================================================================

cat("Creating mission catalog map and datatable...\n")

# Make HTML datatable of mission catalog
dt = make_mission_catalog_datatable(
  mission_summary = mission_summary,
  website_static_path = WEBSITE_STATIC_PATH,
  datatable_header_files_dir = CURATION_DATATABLE_HEADER_FILES_DIR,
  mission_catalog_datatable_dir = CURATION_MISSION_CATALOG_DATATABLE_DIR,
  mission_catalog_datatable_filename = CURATION_MISSION_CATALOG_DATATABLE_FILENAME
)

# Make leaflet map of mission catalog
m = make_mission_catalog_map(
  mission_summary = mission_summary,
  website_static_path = WEBSITE_STATIC_PATH,
  leaflet_header_files_dir = CURATION_LEAFLET_HEADER_FILES_DIR,
  mission_catalog_map_dir = CURATION_MISSION_CATALOG_MAP_DIR,
  mission_catalog_map_filename = CURATION_MISSION_CATALOG_MAP_FILENAME
)

# Copy catalog template to website content directory
catalog_page_path = file.path(WEBSITE_CONTENT_PATH, "data/drone-curation/_index.html")
dir.create(dirname(catalog_page_path), showWarnings = FALSE, recursive = TRUE)
file.copy(
  CURATION_MISSION_CATALOG_TEMPLATE_FILEPATH,
  catalog_page_path,
  overwrite = TRUE
)

cat(sprintf("  Catalog page saved to: %s\n", catalog_page_path))
cat(sprintf("  Catalog map saved to: %s\n",
            file.path(WEBSITE_STATIC_PATH, CURATION_MISSION_CATALOG_MAP_DIR,
                      CURATION_MISSION_CATALOG_MAP_FILENAME)))
cat(sprintf("  Catalog datatable saved to: %s\n\n",
            file.path(WEBSITE_STATIC_PATH, CURATION_MISSION_CATALOG_DATATABLE_DIR,
                      CURATION_MISSION_CATALOG_DATATABLE_FILENAME)))

# ============================================================================
# Create mission detail pages
# ============================================================================

# Function to create a single curation page for a mission
make_mission_curation_page = function(mission_id_foc,
                                      all_mission_ids,
                                      mission_summary,
                                      mission_points,
                                      post_curation_mission_metadata = NULL,
                                      post_curation_image_points = NULL,
                                      missions_with_post_curation_data = character(0)) {

  cat(sprintf("Making curation page for mission %s\n", mission_id_foc))

  # Extract the mission-level metadata
  mission_summary_foc = mission_summary |> filter(mission_id == mission_id_foc)

  # Get the mission points for this mission
  mission_points_foc = mission_points |> filter(mission_id == mission_id_foc)

  # Make pre-curation details map
  mission_details_map_path = make_mission_details_map(
    mission_summary_foc = mission_summary_foc,
    mission_points_foc = mission_points_foc,
    mission_polygons_for_mission_details_map = mission_summary,
    mission_centroids = st_centroid(mission_summary),
    website_static_path = WEBSITE_STATIC_PATH,
    leaflet_header_files_dir = CURATION_LEAFLET_HEADER_FILES_DIR,
    mission_details_map_dir = CURATION_MISSION_DETAILS_MAP_DIR
  )

  # Make pre-curation details datatable
  mission_details_datatable_path = make_mission_details_datatable(
    mission_summary_foc = mission_summary_foc,
    website_static_path = WEBSITE_STATIC_PATH,
    datatable_header_files_dir = CURATION_DATATABLE_HEADER_FILES_DIR,
    mission_details_datatable_dir = CURATION_MISSION_DETAILS_DATATABLE_DIR
  )

  # Check if this mission has post-curation data and generate post-curation widgets if so
  has_post_curation_data = mission_id_foc %in% missions_with_post_curation_data
  post_curation_map_path = NA
  post_curation_datatable_path = NA

  if (has_post_curation_data && !is.null(post_curation_image_points) && !is.null(post_curation_mission_metadata)) {
    cat(sprintf("  Generating post-curation widgets for mission %s\n", mission_id_foc))

    # Get post-curation mission metadata for this mission
    post_curation_mission_summary_foc = post_curation_mission_metadata |>
      filter(mission_id == mission_id_foc)

    # Get post-curation image points for this mission
    post_curation_points_foc = post_curation_image_points |>
      filter(mission_id == mission_id_foc)

    # Make post-curation details map
    post_curation_map_path = make_mission_details_map(
      mission_summary_foc = post_curation_mission_summary_foc,
      mission_points_foc = post_curation_points_foc,
      mission_polygons_for_mission_details_map = mission_summary,
      mission_centroids = st_centroid(mission_summary),
      website_static_path = WEBSITE_STATIC_PATH,
      leaflet_header_files_dir = CURATION_LEAFLET_HEADER_FILES_DIR,
      mission_details_map_dir = CURATION_POST_MISSION_DETAILS_MAP_DIR
    )

    # Make post-curation details datatable
    post_curation_datatable_path = make_mission_details_datatable(
      mission_summary_foc = post_curation_mission_summary_foc,
      website_static_path = WEBSITE_STATIC_PATH,
      datatable_header_files_dir = CURATION_DATATABLE_HEADER_FILES_DIR,
      mission_details_datatable_dir = CURATION_POST_MISSION_DETAILS_DATATABLE_DIR
    )
  }

  # Compute previous and next dataset for navigation
  current_index = which(all_mission_ids == mission_id_foc)
  if (current_index == 1) {
    previous_mission_id = all_mission_ids[length(all_mission_ids)]
  } else {
    previous_mission_id = all_mission_ids[current_index - 1]
  }
  if (current_index == length(all_mission_ids)) {
    next_mission_id = all_mission_ids[1]
  } else {
    next_mission_id = all_mission_ids[current_index + 1]
  }

  next_dataset_page_path = paste0(
    "/", CURATION_MISSION_DETAILS_PAGE_DIR, "/", next_mission_id
  )
  previous_dataset_page_path = paste0(
    "/", CURATION_MISSION_DETAILS_PAGE_DIR, "/", previous_mission_id
  )

  # Render curation page from template
  render_mission_details_page(
    template_filepath = CURATION_MISSION_DETAILS_TEMPLATE_FILEPATH,
    mission_summary_foc = mission_summary_foc,
    s3_file_listing = data.frame(),  # Empty, not needed for curation view
    mission_details_map_path = mission_details_map_path,
    itd_map_path = NA,  # Not included in curation view
    mission_details_datatable_path = mission_details_datatable_path,
    next_dataset_page_path = next_dataset_page_path,
    previous_dataset_page_path = previous_dataset_page_path,
    website_repo_content_path = WEBSITE_CONTENT_PATH,
    mission_details_page_dir = CURATION_MISSION_DETAILS_PAGE_DIR,
    display_data = FALSE,  # No S3 data products in curation view
    # Post-curation parameters for side-by-side display
    has_post_curation_data = has_post_curation_data,
    post_curation_map_html_path = post_curation_map_path,
    post_curation_datatable_html_path = post_curation_datatable_path
  )

  gc()
}

cat("Creating mission detail pages...\n\n")

# Ensure mission details page directory exists
mission_details_content_dir = file.path(WEBSITE_CONTENT_PATH, CURATION_MISSION_DETAILS_PAGE_DIR)
if (!dir.exists(mission_details_content_dir)) {
  dir.create(mission_details_content_dir, showWarnings = FALSE, recursive = TRUE)
}

# Set up parallel processing
plan(multisession, workers = parallel::detectCores()*2)

# Process missions in parallel with progress reporting
future_walk(
  mission_ids,
  make_mission_curation_page,
  all_mission_ids = mission_ids,
  mission_summary = mission_summary,
  mission_points = mission_points,
  post_curation_mission_metadata = post_curation_mission_metadata,
  post_curation_image_points = post_curation_image_points,
  missions_with_post_curation_data = missions_with_post_curation_data,
  .progress = TRUE
)

# Reset to sequential processing
plan(sequential)

cat("\n=== Curation site generation complete! ===\n")
cat(sprintf("Generated mission catalog page: %s\n",
            file.path(WEBSITE_CONTENT_PATH, "data/drone-curation/_index.html")))
cat(sprintf("Generated %d mission detail pages\n", length(mission_ids)))
cat(sprintf("Mission detail pages saved to: %s\n",
            file.path(WEBSITE_CONTENT_PATH, CURATION_MISSION_DETAILS_PAGE_DIR)))
cat(sprintf("Mission detail maps saved to: %s\n",
            file.path(WEBSITE_STATIC_PATH, CURATION_MISSION_DETAILS_MAP_DIR)))
cat(sprintf("Mission detail datatables saved to: %s\n",
            file.path(WEBSITE_STATIC_PATH, CURATION_MISSION_DETAILS_DATATABLE_DIR)))
