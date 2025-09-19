# Purpose: Compare the image metadata across two different versions of the metadata creation to
# determine which missions change, as these will be the ones to re-process in downstram steps.

library(tidyverse)
library(sf)

METADATA_OLD_FILEPATH = "/ofo-share/catalog-data-prep/01_raw-imagery-ingestion/metadata_backup20250918/3_final/ofo-all-images-metadata_backup20250918.gpkg"
METADATA_NEW_FILEPATH = "/ofo-share/catalog-data-prep/01_raw-imagery-ingestion/metadata/3_final/ofo-all-images-metadata.gpkg"

old = st_read(METADATA_OLD_FILEPATH) |> st_drop_geometry()
new = st_read(METADATA_NEW_FILEPATH) |> st_drop_geometry()

missions_old = unique(old$mission_id)
missions_new = unique(new$mission_id)

setdiff(missions_old, missions_new)
setdiff(missions_new, missions_old)

missions_changed = c()
missions_unchanged = c()

for (mission_foc in missions_new) {
  print(mission_foc)

  old_foc = old |>
    filter(mission_id == mission_foc) |>
    arrange(image_path_contrib) |>
    select(original_file_name, image_path_contrib, image_path_ofo)

  new_foc = new |>
    filter(mission_id == mission_foc) |>
    arrange(image_path_contrib) |>
    select(original_file_name, image_path_contrib, image_path_ofo)

  changed = !identical(old_foc, new_foc)
  if (changed) {
    missions_changed = c(missions_changed, mission_foc)
  } else {
    missions_unchanged = c(missions_unchanged, mission_foc)}
}
