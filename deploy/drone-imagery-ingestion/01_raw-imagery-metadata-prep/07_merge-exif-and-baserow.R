# Purpose: At the mission level and the sub-mission level (separately), merge the human-provided Baserow metadata
# and the summarized (mission or sub-mission level) EXIF metadata extracted from the drone imagery.

library(tidyverse)
library(sf)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")


# ============================================================================
# Merge function
# ============================================================================

#' Merge derived EXIF metadata with contributed (Baserow) metadata
#'
#' Combines derived metadata (computed from EXIF) with contributed metadata (from Baserow)
#' at both mission and sub-mission levels. Optionally adds curation anomaly columns.
#'
#' @param mission_foc Mission ID to process
#' @param derived_mission_path Directory containing derived mission-level gpkg files
#' @param derived_sub_mission_path Directory containing derived sub-mission-level gpkg files
#' @param contributed_mission_path Directory containing contributed mission-level CSV files
#' @param contributed_sub_mission_path Directory containing contributed sub-mission-level CSV files
#' @param output_mission_path Directory for output full mission-level gpkg files
#' @param output_sub_mission_path Directory for output full sub-mission-level gpkg files
#' @param curation_notes_mission Optional: tibble with mission-level curation anomaly columns.
#'   Must have 'mission_id' column for joining.
#' @param curation_notes_sub_mission Optional: tibble with sub-mission-level curation anomaly columns.
#'   Must have 'sub_mission_id' column for joining.
#' @return TRUE on success, FALSE on failure (with warning)
merge_derived_and_contributed_metadata = function(
    mission_foc,
    derived_mission_path,
    derived_sub_mission_path,
    contributed_mission_path,
    contributed_sub_mission_path,
    output_mission_path,
    output_sub_mission_path,
    curation_notes_mission = NULL,
    curation_notes_sub_mission = NULL
) {

  # --- Mission-level merge ---

  # Derived filepaths
  baserow_mission_filepath = file.path(contributed_mission_path, paste0(mission_foc, ".csv"))
  exif_metadata_mission_filepath = file.path(derived_mission_path, paste0(mission_foc, ".gpkg"))

  # Check files exist

  if (!file.exists(baserow_mission_filepath)) {
    warning(paste("Contributed mission metadata not found:", baserow_mission_filepath))
    return(FALSE)
  }
  if (!file.exists(exif_metadata_mission_filepath)) {
    warning(paste("Derived mission metadata not found:", exif_metadata_mission_filepath))
    return(FALSE)
  }

  # Load the mission-level metadata
  baserow_mission = read_csv(baserow_mission_filepath, show_col_types = FALSE)
  exif_metadata_mission = st_read(exif_metadata_mission_filepath, quiet = TRUE)

  # dataset_id is confusing so it is just dropped in general
  # In the context of a mission, the sub_mission_id has different meanings in baserow and the exif.
  # In baserow, it means the list of sub-missions included in the mission that a given image is part of.
  # In the exif it is just the sub-mission that a given image is a part of.
  if ("sub_mission_id" %in% names(baserow_mission)) {
    baserow_mission = baserow_mission |>
      # Rename the sub_mission_id to be more descriptive (its a comma-separated list of sub-mission IDs)
      rename(sub_mission_ids = sub_mission_id)
  }

  # Bind together the derived and contributed mission-level metadata
  full_metadata_mission = bind_cols(
    exif_metadata_mission |>
      select(-any_of(c("mission_id", "dataset_id"))), # remove redundant/confusing cols
    baserow_mission
  ) |>
    # Put all the derived columns (which end in _derived) at the end
    select(!ends_with("_derived"), everything())

  # Add curation anomaly columns if provided
  if (!is.null(curation_notes_mission) && nrow(curation_notes_mission) > 0) {
    curation_row = curation_notes_mission |> filter(mission_id == mission_foc)
    if (nrow(curation_row) == 1) {
      # Add each curation column (except mission_id which is for joining)
      for (col in setdiff(names(curation_row), "mission_id")) {
        full_metadata_mission[[col]] = curation_row[[col]]
      }
    } else {
      # No curation notes for this mission - add NA columns
      for (col in setdiff(names(curation_notes_mission), "mission_id")) {
        full_metadata_mission[[col]] = NA_character_
      }
    }
  }

  # Write mission-level output
  create_dir(output_mission_path)
  full_metadata_mission_filepath = file.path(
    output_mission_path,
    paste0(mission_foc, "_mission-metadata.gpkg")
  )
  st_write(full_metadata_mission, full_metadata_mission_filepath, delete_dsn = TRUE, quiet = TRUE)

  # --- Sub-mission-level merge ---

  # Get the sub-missions that make up the mission
  sub_mission_files = list.files(
    path = contributed_sub_mission_path,
    pattern = paste0(mission_foc, "-[0-9]{2}\\.csv$"),
    full.names = TRUE
  )
  sub_mission_ids = sub_mission_files |>
    basename() |>
    str_extract(paste0(mission_foc, "-[0-9]{2}"))

  create_dir(output_sub_mission_path)

  for (sub_mission_id_foc in sub_mission_ids) {
    # Check derived sub-mission file exists
    exif_sub_mission_filepath = file.path(
      derived_sub_mission_path,
      paste0(sub_mission_id_foc, ".gpkg")
    )
    if (!file.exists(exif_sub_mission_filepath)) {
      warning(paste("Derived sub-mission metadata not found:", exif_sub_mission_filepath))
      next
    }

    # Load the sub-mission metadata
    baserow_sub_mission = read_csv(
      file.path(contributed_sub_mission_path, paste0(sub_mission_id_foc, ".csv")),
      show_col_types = FALSE
    )
    exif_metadata_sub_mission = st_read(exif_sub_mission_filepath, quiet = TRUE)

    # Bind together the derived and contributed sub-mission-level metadata
    full_metadata_sub_mission = bind_cols(
      exif_metadata_sub_mission |>
        select(-any_of(c("sub_mission_id", "mission_id", "dataset_id"))),
      baserow_sub_mission
    ) |>
      # Put all the derived columns (which end in _derived) at the end
      select(!ends_with("_derived"), everything())

    # Add curation anomaly columns if provided
    if (!is.null(curation_notes_sub_mission) && nrow(curation_notes_sub_mission) > 0) {
      curation_sub_row = curation_notes_sub_mission |> filter(sub_mission_id == sub_mission_id_foc)
      if (nrow(curation_sub_row) == 1) {
        for (col in setdiff(names(curation_sub_row), "sub_mission_id")) {
          full_metadata_sub_mission[[col]] = curation_sub_row[[col]]
        }
      } else {
        for (col in setdiff(names(curation_notes_sub_mission), "sub_mission_id")) {
          full_metadata_sub_mission[[col]] = NA_character_
        }
      }
    }

    # Write sub-mission-level output
    full_metadata_sub_mission_filepath = file.path(
      output_sub_mission_path,
      paste0(sub_mission_id_foc, "_sub-mission-metadata.gpkg")
    )
    st_write(full_metadata_sub_mission, full_metadata_sub_mission_filepath, delete_dsn = TRUE, quiet = TRUE)
  }

  return(TRUE)
}


# ============================================================================
# Workflow
# ============================================================================

# Determine which missions to process
missions_to_process = read_csv(MISSIONS_TO_PROCESS_RAW_IMAGERY_METADATA_LIST_PATH) |>
  pull(mission_id)

# Create the output folders
create_dir(FULL_METADATA_PER_MISSION_PATH)
create_dir(FULL_METADATA_PER_SUB_MISSION_PATH)


# Parallelize across all selected missions
future::plan(multisession)
future_walk(
  missions_to_process,
  ~ merge_derived_and_contributed_metadata(
    mission_foc = .x,
    derived_mission_path = DERIVED_METADATA_PER_MISSION_PATH,
    derived_sub_mission_path = DERIVED_METADATA_PER_SUB_MISSION_PATH,
    contributed_mission_path = EXTRACTED_METADATA_PER_MISSION_PATH,
    contributed_sub_mission_path = EXTRACTED_METADATA_PER_SUB_MISSION_PATH,
    output_mission_path = FULL_METADATA_PER_MISSION_PATH,
    output_sub_mission_path = FULL_METADATA_PER_SUB_MISSION_PATH
  ),
  .progress = TRUE
)
