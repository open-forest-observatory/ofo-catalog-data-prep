# Purpose: Get all mission polygon and image point metadata from the object store and compile it in
# a single gpkg each (onr for point and one for poly). TODO: Eventually we will replace this with database queries.

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
remote_dir = paste0(RCLONE_REMOTE, ":", REMOTE_MISSIONS_DIR)
command = paste("rclone lsf", remote_dir, "-R --files-only", sep = " ")
listing = system(command, intern = TRUE)
listing_df = tibble(filepath = listing)

filepath_parts = str_split(listing_df$filepath, "/")
listing_df$mission_id = map_chr(filepath_parts, 1)


write_csv(listing_df, S3_LISTING_FILEPATH)



# But don't use the file listing now, just download all the polygons using filtering
tempdir = file.path(TEMPDIR, "mission-polygons")
unlink(tempdir, recursive = TRUE, force = TRUE)
dir.create(tempdir, recursive = TRUE, showWarnings = FALSE)

command = paste("rclone copy", remote_dir, tempdir, "--include '**/*_mission-metadata.gpkg'", "--progress --transfers 32 --checkers 32 --stats 1s --retries 5 --retries-sleep=15s --s3-upload-cutoff 100Mi --s3-chunk-size 100Mi --s3-upload-concurrency 16 --multi-thread-streams 2", sep = " ")
result = system(command)
# Get paths to all the downloaded files
mission_polygon_files = list.files(tempdir, pattern = "_mission-metadata.gpkg$", full.names = TRUE, recursive = TRUE)

# Exclude mission 000000 which is a dummy mission used for testing
mission_polygon_files = mission_polygon_files[!grepl("000000_mission-metadata.gpkg", mission_polygon_files)]

# Open and combine all the downloaded mission polygons in single SF object
mission_polygons = future_map(mission_polygon_files, st_read)

# Need to set all cols to character, because for some missions, a given column is inferred to be
# character, whereas for other, it is inferred to be numeric. This happens for example when one
# mission has the image count as "121, 324, 134" because it is has multiple sub-missions, but
# another has image count of "121" because it is a single mission.
mission_polygons = future_map(mission_polygons, force_all_cols_to_character)
mission_polygons = bind_rows(mission_polygons)

st_write(mission_polygons, MISSION_METADATA_FILEPATH, delete_dsn = TRUE)

unlink(tempdir, recursive = TRUE, force = TRUE)



# And now for mission points
tempdir = file.path(TEMPDIR, "mission-points")
unlink(tempdir, recursive = TRUE, force = TRUE)
dir.create(tempdir, recursive = TRUE, showWarnings = FALSE)

command = paste("rclone copy", remote_dir, tempdir, "--include '**/*_image-metadata.gpkg'", "--progress --transfers 32 --checkers 32 --stats 1s --retries 5 --retries-sleep=15s --s3-upload-cutoff 100Mi --s3-chunk-size 100Mi --s3-upload-concurrency 16 --multi-thread-streams 2", sep = " ")
result = system(command)
# Get paths to all the downloaded files
mission_point_files = list.files(tempdir, pattern = "_image-metadata.gpkg$", full.names = TRUE, recursive = TRUE)

# Exclude mission 000000 which is a dummy mission used for testing
mission_point_files = mission_point_files[!grepl("000000_image-metadata.gpkg", mission_point_files)]

# Open and combine all the downloaded mission points in single SF object
mission_points = future_map(mission_point_files, st_read)

# Need to set all cols to character, because for some missions, a given column is inferred to be
# character, whereas for other, it is inferred to be numeric. This happens for example when one
# mission has the image count as "121, 324, 134" because it is has multiple sub-missions, but
# another has image count of "121" because it is a single mission.
mission_points = future_map(mission_points, force_all_cols_to_character)
mission_points = bind_rows(mission_points)

# TODO: TEMPORARY: Override image_id to be the same as the image filename (sans extension)
mission_points$image_id = tools::file_path_sans_ext(basename(mission_points$image_path_ofo))

st_write(mission_points, IMAGE_METADATA_FILEPATH, delete_dsn = TRUE)

unlink(tempdir, recursive = TRUE, force = TRUE)
