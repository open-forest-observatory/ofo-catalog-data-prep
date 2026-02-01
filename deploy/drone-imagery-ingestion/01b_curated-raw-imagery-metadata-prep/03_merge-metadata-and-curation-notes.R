# Purpose: Merge post-curation derived metadata with contributed metadata and curation notes.
#
# Key difference from pre-curation: Uses pre-curation full metadata gpkgs as source for
# contributed metadata (not raw baserow CSVs). This ensures consistency with what was
# used for pre-curation and avoids reconciliation issues.
#
# Workflow:
# 1. Load processed curation notes from script 01
# 2. For each mission:
#    - Load post-curation derived metadata (from script 02)
#    - Extract contributed columns from pre-curation full metadata gpkgs
#    - Merge derived and contributed metadata
#    - Add curation anomaly columns
#    - Write to post-curation full metadata paths

library(tidyverse)
library(sf)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")

# ============================================================================
# Define columns that are "contributed" vs "derived" in pre-curation metadata
# ============================================================================

# Columns from the pre-curation gpkg that are derived (computed from EXIF)
# These will be replaced by newly-computed post-curation derived values
DERIVED_COLUMN_SUFFIXES = c("_derived")
GEOMETRY_COLUMNS = c("geom", "geometry", "x")  # sf geometry column names

#' Extract contributed-only columns from pre-curation full metadata
#'
#' @param full_metadata_sf SF object with full pre-curation metadata
#' @return Tibble with only contributed columns (no derived or geometry)
extract_contributed_columns = function(full_metadata_sf) {
  # Get column names
  all_cols = names(full_metadata_sf)

  # Identify derived columns (those ending in _derived)
  derived_cols = all_cols[str_detect(all_cols, "_derived$")]

  # Identify geometry columns
  geom_cols = intersect(all_cols, GEOMETRY_COLUMNS)

  # Contributed columns = all columns minus derived and geometry
  contributed_cols = setdiff(all_cols, c(derived_cols, geom_cols))

  # Return as tibble (drop geometry)
  full_metadata_sf |>
    st_drop_geometry() |>
    select(all_of(contributed_cols))
}

# ============================================================================
# Load processed curation notes
# ============================================================================

curation_notes_path = file.path(POST_CURATION_INTERMEDIATE_PATH, "curation-notes-processed.csv")
curation_notes = read_csv(curation_notes_path, col_types = cols(.default = "c"))

cat(sprintf("Loaded %d rows of processed curation notes\n", nrow(curation_notes)))

# ============================================================================
# Prepare curation notes for joining at mission and sub-mission levels
# ============================================================================

#' Prepare curation notes for sub-mission level joining
#'
#' @param curation_notes Tibble with curation notes
#' @return Tibble with sub_mission_id and renamed anomaly columns
prepare_curation_for_sub_mission = function(curation_notes) {
  curation_notes |>
    mutate(
      sub_mission_id_full = paste0(mission_id, "-", str_pad(sub_mission_id, 2, pad = "0"))
    ) |>
    select(
      sub_mission_id = sub_mission_id_full,
      anomaly_severity,
      anomalies_collection_time = collection_time_anomaly,
      anomalies_altitude = altitude_anomaly,
      anomalies_spatial = spatial_anomaly,
      anomalies_camera_pitch = camera_pitch_anomaly,
      anomalies_excess_images = excess_images_anomaly,
      anomalies_missing_images = missing_images_anomaly,
      anomalies_other = other_anomalies,
      anomaly_notes
    )
}

#' Prepare curation notes for mission level joining
#'
#' Aggregates sub-mission level notes to mission level by combining text
#'
#' @param curation_notes Tibble with curation notes
#' @return Tibble with mission_id and aggregated anomaly columns
prepare_curation_for_mission = function(curation_notes) {
  # Aggregate across sub-missions within each mission
  curation_notes |>
    group_by(mission_id) |>
    summarise(
      anomaly_severity = paste(na.omit(anomaly_severity[anomaly_severity != ""]), collapse = "; "),
      anomalies_collection_time = paste(na.omit(collection_time_anomaly[collection_time_anomaly != ""]), collapse = "; "),
      anomalies_altitude = paste(na.omit(altitude_anomaly[altitude_anomaly != ""]), collapse = "; "),
      anomalies_spatial = paste(na.omit(spatial_anomaly[spatial_anomaly != ""]), collapse = "; "),
      anomalies_camera_pitch = paste(na.omit(camera_pitch_anomaly[camera_pitch_anomaly != ""]), collapse = "; "),
      anomalies_excess_images = paste(na.omit(excess_images_anomaly[excess_images_anomaly != ""]), collapse = "; "),
      anomalies_missing_images = paste(na.omit(missing_images_anomaly[missing_images_anomaly != ""]), collapse = "; "),
      anomalies_other = paste(na.omit(other_anomalies[other_anomalies != ""]), collapse = "; "),
      anomaly_notes = paste(na.omit(anomaly_notes[anomaly_notes != ""]), collapse = "; "),
      .groups = "drop"
    ) |>
    # Replace empty strings with NA
    mutate(across(everything(), ~ if_else(.x == "", NA_character_, .x)))
}

curation_sub_mission = prepare_curation_for_sub_mission(curation_notes)
curation_mission = prepare_curation_for_mission(curation_notes)

cat(sprintf("Prepared curation notes: %d mission-level, %d sub-mission-level\n",
            nrow(curation_mission), nrow(curation_sub_mission)))

# ============================================================================
# Get missions to process
# ============================================================================

post_curation_files = list.files(POST_CURATION_DERIVED_METADATA_PER_MISSION_PATH,
                                  pattern = "\\.gpkg$", full.names = FALSE)
missions_to_process = str_extract(post_curation_files, "^\\d{6}")

cat(sprintf("Merging metadata for %d missions...\n", length(missions_to_process)))

# Create output directories
create_dir(POST_CURATION_FULL_METADATA_PER_MISSION_PATH)
create_dir(POST_CURATION_FULL_METADATA_PER_SUB_MISSION_PATH)

# ============================================================================
# Merge function for a single mission
# ============================================================================

merge_metadata_for_mission = function(mission_foc) {
  # --- Mission-level merge ---

  # Load derived (post-curation) metadata
  derived_mission_filepath = file.path(POST_CURATION_DERIVED_METADATA_PER_MISSION_PATH,
                                        paste0(mission_foc, ".gpkg"))

  if (!file.exists(derived_mission_filepath)) {
    warning(paste("No post-curation derived metadata found for mission", mission_foc))
    return(FALSE)
  }

  derived_mission = st_read(derived_mission_filepath, quiet = TRUE)

  # Load contributed metadata from PRE-CURATION full metadata gpkg (not raw baserow CSV)
  pre_curation_mission_filepath = file.path(FULL_METADATA_PER_MISSION_PATH,
                                             paste0(mission_foc, "_mission-metadata.gpkg"))

  if (!file.exists(pre_curation_mission_filepath)) {
    warning(paste("No pre-curation full metadata found for mission", mission_foc))
    return(FALSE)
  }

  pre_curation_full = st_read(pre_curation_mission_filepath, quiet = TRUE)
  contributed_mission = extract_contributed_columns(pre_curation_full)

  # Merge derived (post-curation) and contributed
  # Keep geometry from derived, remove redundant mission_id
  full_metadata_mission = bind_cols(
    derived_mission |> select(-mission_id),
    contributed_mission
  ) |>
    # Put derived columns (ending in _derived) at the end
    select(!ends_with("_derived"), everything())

  # Add curation anomaly columns
  curation_row = curation_mission |> filter(mission_id == mission_foc)
  if (nrow(curation_row) == 1) {
    for (col in setdiff(names(curation_row), "mission_id")) {
      full_metadata_mission[[col]] = curation_row[[col]]
    }
  } else {
    # No curation notes for this mission - add NA columns
    for (col in setdiff(names(curation_mission), "mission_id")) {
      full_metadata_mission[[col]] = NA_character_
    }
  }

  # Write mission-level metadata
  output_filepath = file.path(POST_CURATION_FULL_METADATA_PER_MISSION_PATH,
                               paste0(mission_foc, "_mission-metadata.gpkg"))
  st_write(full_metadata_mission, output_filepath, delete_dsn = TRUE, quiet = TRUE)

  # --- Sub-mission-level merge ---

  sub_mission_files = list.files(
    POST_CURATION_DERIVED_METADATA_PER_SUB_MISSION_PATH,
    pattern = paste0(mission_foc, "-\\d{2}\\.gpkg$"),
    full.names = TRUE
  )

  for (sub_file in sub_mission_files) {
    sub_mission_id_foc = str_extract(basename(sub_file), "\\d{6}-\\d{2}")

    derived_sub = st_read(sub_file, quiet = TRUE)

    # Load contributed metadata from pre-curation full metadata gpkg
    pre_curation_sub_filepath = file.path(FULL_METADATA_PER_SUB_MISSION_PATH,
                                           paste0(sub_mission_id_foc, "_sub-mission-metadata.gpkg"))

    if (!file.exists(pre_curation_sub_filepath)) {
      warning(paste("No pre-curation full metadata found for sub-mission", sub_mission_id_foc))
      next
    }

    pre_curation_sub_full = st_read(pre_curation_sub_filepath, quiet = TRUE)
    contributed_sub = extract_contributed_columns(pre_curation_sub_full)

    # Merge derived and contributed
    full_metadata_sub = bind_cols(
      derived_sub |> select(-sub_mission_id, -mission_id),
      contributed_sub
    ) |>
      select(!ends_with("_derived"), everything())

    # Add curation anomaly columns
    curation_sub_row = curation_sub_mission |> filter(sub_mission_id == sub_mission_id_foc)
    if (nrow(curation_sub_row) == 1) {
      for (col in setdiff(names(curation_sub_row), "sub_mission_id")) {
        full_metadata_sub[[col]] = curation_sub_row[[col]]
      }
    } else {
      for (col in setdiff(names(curation_sub_mission), "sub_mission_id")) {
        full_metadata_sub[[col]] = NA_character_
      }
    }

    output_sub_filepath = file.path(POST_CURATION_FULL_METADATA_PER_SUB_MISSION_PATH,
                                     paste0(sub_mission_id_foc, "_sub-mission-metadata.gpkg"))
    st_write(full_metadata_sub, output_sub_filepath, delete_dsn = TRUE, quiet = TRUE)
  }

  return(TRUE)
}

# ============================================================================
# Run for all missions
# ============================================================================

future::plan(multisession(workers = future::availableCores() * 3))
results = future_map_lgl(missions_to_process, merge_metadata_for_mission, .progress = TRUE)

cat(sprintf("\nMerged metadata for %d missions successfully\n", sum(results)))

cat("\n**** Post-curation merge complete ****\n")
