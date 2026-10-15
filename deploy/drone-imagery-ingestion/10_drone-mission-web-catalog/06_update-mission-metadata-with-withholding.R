# Purpose: Add a column to the mission metadata that specifies whether it is to be withheld from
# broad-scale ML training

library(sf)
library(tidyverse)

source("deploy/drone-imagery-ingestion/00_set-constants.R")

mission_metadata = st_read(POST_CURATION_MISSION_METADATA_FILEPATH) # Changed this from MISSION_METADATA_FILEPATH

missions_to_withhold = read_csv(MISSIONS_TO_WITHHOLD_FILEPATH)$mission_id

# Columns carried over from the Baserow export that are empty, internal to Baserow, or superseded
# (withheld_from_training is superseded by withhold_from_training). Temporary fix: drop them here
# at the final export rather than upstream.
leftover_baserow_cols = c(
  "calibration_panel_reflectance", "calibration_panel_id", "x_offset", "y_offset", "curation",
  "gcps.ground", "gcps-ground", "dataset.associations...datasets.imagery",
  "dataset-associations - datasets-imagery", "Created.on", "Created on", "date_modified",
  "withheld_from_training",
  # Vestigial per-sub-mission image count from an earlier pipeline version
  "n_images"
)

# All columns arrive as text (05_get-all-mission-metadata.R casts them so that multi-valued
# contributed columns combine across missions). The _derived columns are never multi-valued, so
# restore their native types. Times and datetimes stay as text: they have no time zone and some
# missions recorded local time while others recorded UTC, so storing them as datetimes would imply
# a time zone that isn't known.
restore_derived_col_types = function(mission_metadata) {
  derived_cols = grep("_derived$", names(mission_metadata), value = TRUE)
  converted = mission_metadata |>
    st_drop_geometry() |>
    select(all_of(derived_cols)) |>
    type_convert(col_types = cols(.default = col_guess()))
  keep_as_text = names(converted)[sapply(converted, function(x) inherits(x, c("POSIXct", "hms", "difftime")))]
  derived_cols = setdiff(derived_cols, keep_as_text)
  mission_metadata[derived_cols] = converted[derived_cols]
  mission_metadata
}

mission_metadata = restore_derived_col_types(mission_metadata)

mission_metadata = mission_metadata |>
  mutate(
    withhold_from_training = mission_id %in% missions_to_withhold
  ) |>
  select(-any_of(leftover_baserow_cols)) |>
  select(everything()) # to get geom back at the end

st_write(mission_metadata, POST_CURATION_MISSION_METADATA_FILEPATH, delete_dsn = TRUE) # Changed this from MISSION_METADATA_FILEPATH
