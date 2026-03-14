# Purpose: Pull S3 file listings for both missions and composites directories and save them as CSVs
# for use by subsequent scripts in this pipeline.

source("deploy/drone-imagery-ingestion/00_set-constants.R")

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
TMP_DIR = file.path(DELIVERABLES_DIR, "tmp")
MISSIONS_S3_LISTING_FILEPATH = file.path(TMP_DIR, "s3-file-listing-missions.csv")
COMPOSITES_S3_LISTING_FILEPATH = file.path(TMP_DIR, "s3-file-listing-composites.csv")

dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)


# --- Missions listing ---

cat("Getting missions file listing from S3...\n")
remote_dir = paste0(RCLONE_REMOTE, ":", REMOTE_MISSIONS_DIR)
command = paste("rclone lsf", remote_dir, "-R --files-only")
listing = system(command, intern = TRUE)
listing_df = tibble(filepath = listing)

# Exclude loose files at the root level (not inside a mission subdirectory)
listing_df = listing_df |> filter(str_detect(filepath, "/"))

filepath_parts = str_split(listing_df$filepath, "/")
listing_df$mission_id = map_chr(filepath_parts, 1)

write_csv(listing_df, MISSIONS_S3_LISTING_FILEPATH)
cat("Wrote", nrow(listing_df), "mission files to", MISSIONS_S3_LISTING_FILEPATH, "\n")


# --- Composites listing ---

cat("Getting composites file listing from S3...\n")
remote_dir = paste0(RCLONE_REMOTE, ":", REMOTE_COMPOSITES_DIR)
command = paste("rclone lsf", remote_dir, "-R --files-only")
listing = system(command, intern = TRUE)
listing_df = tibble(filepath = listing)

# Exclude loose files at the root level (not inside a composite subdirectory)
listing_df = listing_df |> filter(str_detect(filepath, "/"))

filepath_parts = str_split(listing_df$filepath, "/")
listing_df$composite_id = map_chr(filepath_parts, 1)

write_csv(listing_df, COMPOSITES_S3_LISTING_FILEPATH)
cat("Wrote", nrow(listing_df), "composite files to", COMPOSITES_S3_LISTING_FILEPATH, "\n")
