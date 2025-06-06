# Purpose: For each sub-mission, pull the relevant Baserow metadata. For each mission, merge the
# pulled sub-mission data into a mission-level record by concatenating the multiple different values
# where they differ. Save out the mission-level and sub-mission-level data. If the dataset
# association records indicate that a mission is one part of a two-part grid, change the
# mission-level mission type from "normal" to "grid".

# TODO: Currently, we are ignoring cases where there are two baserow records for a mission, but the
# mission could not be split out into sub-missions that we were confident corresponded to the two
# baserow records. The occurrence of this is recorded in the crosswalks in 1c_exif-for-sorting. When
# this occurred, all missions were crosswalked to the first Baserow record. We could improve this
# script's workflow to take the non-directly-assignable records into account and report both entries
# for any attributes that differ between the two records for an unsplittable dataset.

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
baserow_projects = baserow_projects |>
  select(project_id, contributor_names, contact_info, license, objectives)
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

# Start with the sub-mission-level Baserow attributes
sub_mission_baserow = left_join(crosswalk, baserow_datasets, by = c("dataset_id_baserow" = "dataset_id")) |>
  select(-dataset_id_baserow)

# Next, the mission-level Baserow attributes, including concatenating the sub-mission-level
# attributes if they differ

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
