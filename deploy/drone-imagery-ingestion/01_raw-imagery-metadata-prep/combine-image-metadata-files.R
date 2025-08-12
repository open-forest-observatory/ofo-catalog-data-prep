# Purpose: Take all the image-level metadata files in PARSED_EXIF_FOR_RETAINED_IMAGES_PATH and
# combine them into a single file that makes it easier to filter and intersect image metadata NOTE:
# This assumes that all the image metadata files that are potentially relevant are in the directory
# PARSED_EXIF_FOR_RETAINED_IMAGES_PATH. This shouldn't be assumed, and therefore if you're not sure,
# you may wish to download them to this folder from the object store.

library(tidyverse)
library(sf)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")

## Functions

# Force all cols except geometry to character
force_all_cols_to_character = function(df) {
  df |>
    mutate(across(-any_of(c("geometry", "geom")), as.character))
}

metadata_files = list.files(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH, pattern = "_image-metadata.gpkg$", full.names = TRUE)

plan(multisession)
metadata_sf_list = future_map(metadata_files, st_read, .progress = TRUE)

# Need to set all cols to character, because for some missions, a given column is inferred to be
# character, whereas for other, it is inferred to be numeric. This happens for example when one
# mission has the image count as "121, 324, 134" because it is has multiple sub-missions, but
# another has image count of "121" because it is a single mission.
metadata_sf_list_char = future_map(metadata_sf_list, force_all_cols_to_character)


metadata_sf = bind_rows(metadata_sf_list_char)

# The following code was to convert data types back to their native types, but this is not necessary
# (it can be done in the downstream analysis if necessary) and it's problematic here because it
# converts some columns incorrectly to a "difftime" format that is not compatible with gpkg.
# # Remove geometry, re-infer column data types, and re-add geometry
# metadata_geom = st_geometry(metadata_sf)
# metadata_nogeom = st_drop_geometry(metadata_sf)
# metadata_nogeom = type_convert(metadata_nogeom)
# metadata_sf2 = st_sf(metadata_nogeom, metadata_geom)


# Write
st_write(metadata_sf, FULL_METADATA_PER_IMAGE_COMBINED_FILEPATH, delete_dsn = TRUE)
