# Purpose: Get all composite mission polygon and image point metadata from the object store and compile it in
# a single gpkg each (one for polygon and one for points). Follows the same pattern as 05_get-all-mission-metadata.R.

library(sf)
library(tidyverse)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")


# Force all cols except geometry to character
force_all_cols_to_character = function(df) {
  df |>
    mutate(across(-any_of(c("geometry", "geom")), as.character))
}



# ---- Processing

# Query the object store for a file listing
remote_dir = paste0(RCLONE_REMOTE, ":", REMOTE_COMPOSITES_DIR)
command = paste("rclone lsf", remote_dir, "-R --files-only", sep = " ")
listing = system(command, intern = TRUE)
listing_df = tibble(filepath = listing)

filepath_parts = str_split(listing_df$filepath, "/")
listing_df$composite_id = map_chr(filepath_parts, 1)


write_csv(listing_df, COMPOSITE_S3_LISTING_FILEPATH)



# But don't use the file listing now, just download all the polygons using filtering
tempdir = file.path(TEMPDIR, "composite-polygons")
unlink(tempdir, recursive = TRUE, force = TRUE)
dir.create(tempdir, recursive = TRUE, showWarnings = FALSE)

command = paste("rclone copy", remote_dir, tempdir, "--include '**/*_mission-metadata.gpkg'", "--progress --transfers 32 --checkers 32 --stats 1s --retries 5 --retries-sleep=15s --s3-upload-cutoff 100Mi --s3-chunk-size 100Mi --s3-upload-concurrency 16 --multi-thread-streams 2", sep = " ")
result = system(command)
# Get paths to all the downloaded files
composite_polygon_files = list.files(tempdir, pattern = "_mission-metadata.gpkg$", full.names = TRUE, recursive = TRUE)

# Exclude composite 000000_000000 which is a dummy composite used for testing
composite_polygon_files = composite_polygon_files[!grepl("000000_000000_mission-metadata.gpkg", composite_polygon_files)]

# Open and combine all the downloaded composite polygons in single SF object
composite_polygons = future_map(composite_polygon_files, st_read)

# Need to set all cols to character, because for some missions, a given column is inferred to be
# character, whereas for other, it is inferred to be numeric. This happens for example when one
# mission has the image count as "121, 324, 134" because it is has multiple sub-missions, but
# another has image count of "121" because it is a single mission.
composite_polygons = future_map(composite_polygons, force_all_cols_to_character)
composite_polygons = bind_rows(composite_polygons)

st_write(composite_polygons, COMPOSITE_MISSION_METADATA_FILEPATH, delete_dsn = TRUE)

unlink(tempdir, recursive = TRUE, force = TRUE)



# And now for composite points
tempdir = file.path(TEMPDIR, "composite-points")
unlink(tempdir, recursive = TRUE, force = TRUE)
dir.create(tempdir, recursive = TRUE, showWarnings = FALSE)

command = paste("rclone copy", remote_dir, tempdir, "--include '**/*_image-metadata.gpkg'", "--progress --transfers 32 --checkers 32 --stats 1s --retries 5 --retries-sleep=15s --s3-upload-cutoff 100Mi --s3-chunk-size 100Mi --s3-upload-concurrency 16 --multi-thread-streams 2", sep = " ")
result = system(command)
# Get paths to all the downloaded files
composite_point_files = list.files(tempdir, pattern = "_image-metadata.gpkg$", full.names = TRUE, recursive = TRUE)

# Exclude composite 000000_000000 which is a dummy composite used for testing
composite_point_files = composite_point_files[!grepl("000000_000000_image-metadata.gpkg", composite_point_files)]

# Open and combine all the downloaded composite points in single SF object
composite_points = future_map(composite_point_files, st_read)

# Need to set all cols to character, because for some missions, a given column is inferred to be
# character, whereas for other, it is inferred to be numeric. This happens for example when one
# mission has the image count as "121, 324, 134" because it is has multiple sub-missions, but
# another has image count of "121" because it is a single mission.
composite_points = future_map(composite_points, force_all_cols_to_character)
composite_points = bind_rows(composite_points)

# TODO: TEMPORARY: Override image_id to be the same as the image filename (sans extension)
composite_points$image_id = tools::file_path_sans_ext(basename(composite_points$image_path_ofo))

st_write(composite_points, COMPOSITE_IMAGE_METADATA_FILEPATH, delete_dsn = TRUE)

unlink(tempdir, recursive = TRUE, force = TRUE)
