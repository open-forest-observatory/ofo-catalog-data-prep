# Purpose: For each sub-mission, pull the relevant Baserow metadata. For each mission, merge the
# pulled sub-mission data into a mission-level record by concatenating the multiple different values
# where they differ. Save out the mission-level and sub-mission-level data. If the dataset
# association records indicate that a mission is one part of a two-part grid, change the
# mission-level mission type from "normal" to "grid".

library(tidyverse)

source("deploy/drone-imagery-ingestion/00_set-constants.R")

# Derived constants
crosswalk_filepath = file.path(CONTRIBUTED_TO_SORTED_MISSION_ID_CROSSWALK_PATH, paste0(PROJECT_NAME_TO_PROCESS_RAW_IMAGERY_METADATA, ".csv"))

# Pull in the baserow (human-entered) metadata
baserow_datasets = read_csv(file.path(CONTRIBUTED_METADATA_PATH, "export - datasets-imagery.csv"))
baserow_projects = read_csv(file.path(CONTRIBUTED_METADATA_PATH, "export - acquisition-projects.csv"))
baserow_dataset_associations = read.csv(file.path(CONTRIBUTED_METADATA_PATH, "export - dataset-associations - Grid.csv"))

# Fix the formatting of the Baserow dataset_id, and remove the internal baserow 'id' column
baserow_datasets = baserow_datasets |>
  mutate(dataset_id = str_pad(dataset_id, 6, pad = "0", side = "left")) |>
  select(-id)

# Merge some project-level data into the dataset-level data:
# - contributor_names
# - contact_info
# - license
# - objectives
# - forest_focused
baserow_projects = baserow_projects |>
  select(project_id, contributor_names, contact_info, license, objectives, forest_focused)
baserow_datasets = baserow_datasets |>
  left_join(baserow_projects, by = c("project_id" = "project_id"))

# Apply any dataset-level contributor overrides
baserow_datasets = baserow_datasets |>
  mutate(contributor_names = ifelse(!is.na(contributor_names_override), contributor_names_override, contributor_names),
         contact_info = ifelse(!is.na(contact_info_override), contact_info_override, contact_info)) |>
  select(-contributor_names_override, -contact_info_override)


# Load the crosswalk linking sub-mission folder names to baserow rows. If there is more than one sub-mission for
# a mission, they may or may not have different entries in Baserow, and their extracted
# exif may or may not have differences.
crosswalk = read_csv(crosswalk_filepath)

# We need to create a set of attributes at the mission level and the sub-mission level. At the
# mission level, we can generate the *derived* attributes from the exif data, but the *manually
# provided* attributes from Baserow like base station location and drone model will need to be
# pulled from both matching baserow rows and concatenated

concat_unique = function(x) {
  x = x[!is.na(x)]
  unq = unique(x)
  if (length(unq) > 1) {
    return(paste(unq, collapse = ", "))
  } else if (length(unq) == 0) {
    return(NA)
  } else {
    return(unq)
  }
}

# Start with the sub-mission-level Baserow attributes
sub_mission_baserow = left_join(crosswalk, baserow_datasets, by = c("dataset_id_baserow" = "dataset_id"))

# For sub-missions that came from unsplittable _and_ folders (where multiple Baserow records
# correspond to the same images but couldn't be separated), merge in the additional Baserow
# records' attributes
if ("addl_dataset_ids_baserow" %in% colnames(sub_mission_baserow)) {

  # Identify rows with additional dataset IDs
  rows_with_addl = which(!is.na(sub_mission_baserow$addl_dataset_ids_baserow))

  if (length(rows_with_addl) > 0) {

    # Get all columns that came from baserow (everything except crosswalk columns)
    crosswalk_cols = c("dataset_id_baserow", "sub_mission_id", "mission_id", "n_images",
                       "project_name", "addl_dataset_ids_baserow", "addl_baserow_differ_by",
                       "why_not_separable")
    baserow_cols = setdiff(colnames(sub_mission_baserow), crosswalk_cols)

    for (i in rows_with_addl) {

      addl_ids = sub_mission_baserow$addl_dataset_ids_baserow[i] |>
        str_split(",") |>
        unlist() |>
        str_trim() |>
        str_pad(6, pad = "0", side = "left")

      # Look up the additional Baserow records
      addl_baserow = baserow_datasets |>
        filter(dataset_id %in% addl_ids)

      if (nrow(addl_baserow) > 0) {
        # For each baserow column, concatenate unique values from primary and additional records
        for (col in baserow_cols) {
          primary_val = sub_mission_baserow[[col]][i]
          addl_vals = addl_baserow[[col]]
          all_vals = c(primary_val, addl_vals)
          sub_mission_baserow[[col]][i] = concat_unique(all_vals)
        }
      }
    }
  }
}

sub_mission_baserow = sub_mission_baserow |>
  select(-dataset_id_baserow, -any_of(c("addl_dataset_ids_baserow", "addl_baserow_differ_by", "why_not_separable")))

# Next, the mission-level Baserow attributes, including concatenating the sub-mission-level
# attributes if they differ

mission_baserow = sub_mission_baserow |>
  # At the mission level, the dataset_id is the mission_id
  mutate(across(everything(), as.character)) |>
  group_by(mission_id) |>
  summarize(across(everything(), concat_unique)) |>
  ungroup()

# Determine if the combined mission comprises a grid, by looking up all the dataset_ids (these are
# mission IDs, not sub-mission IDs) from the dataset_associations table that have an association
# type of "multi-orientation"

grid_missions = baserow_dataset_associations |>
  filter(assoc_type == "multi-orientation") |>
  pull(dataset_ids) |>
  str_split(pattern = ",") |>
  unlist() |>
  unique() |>
  str_pad(6, pad = "0", side = "left")

# We have the Baserow dataset IDs for the missions. Look up the sub-mission IDs (folder names) for
# them.
grid_mission_ids = crosswalk |>
  filter(dataset_id_baserow %in% grid_missions) |>
  pull(mission_id) |>
  unique()

# Identify which mission-level baserow records match grid missions and set the flight pattern to
# "grid"
mission_baserow[mission_baserow$mission_id %in% grid_mission_ids, "flight_pattern"] = "grid"

# Save out the metadata, one file per mission or sub-mission
dir.create(EXTRACTED_METADATA_PER_MISSION_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(EXTRACTED_METADATA_PER_SUB_MISSION_PATH, recursive = TRUE, showWarnings = FALSE)

mission_ids = unique(mission_baserow$mission_id)

for (mission_id_foc in mission_ids) {

  mission_baserow_foc = mission_baserow |>
    filter(mission_id == mission_id_foc)

  file_out = file.path(EXTRACTED_METADATA_PER_MISSION_PATH, paste0(mission_id_foc, ".csv"))

  write_csv(mission_baserow_foc, file_out)

}

sub_mission_ids = unique(sub_mission_baserow$sub_mission_id)

for (sub_mission_id_foc in sub_mission_ids) {

  sub_mission_baserow_foc = sub_mission_baserow |>
    filter(sub_mission_id == sub_mission_id_foc)

  file_out = file.path(EXTRACTED_METADATA_PER_SUB_MISSION_PATH, paste0(sub_mission_id_foc, ".csv"))

  write_csv(sub_mission_baserow_foc, file_out)

}
