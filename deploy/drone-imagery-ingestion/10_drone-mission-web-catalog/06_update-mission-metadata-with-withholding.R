# Purpose: Add a column to the mission metadata that specifies whether it is to be withheld from
# broad-scale ML training

library(sf)
library(tidyverse)

source("deploy/drone-imagery-ingestion/00_set-constants.R")

mission_metadata = st_read(MISSION_METADATA_FILEPATH)

missions_to_withhold = read_csv(MISSIONS_TO_WITHHOLD_FILEPATH)$mission_id

mission_metadata = mission_metadata |>
  mutate(
    withhold_from_training = mission_id %in% missions_to_withhold
  ) |>
  select(everything()) # to get geom back at the end

st_write(mission_metadata, MISSION_METADATA_FILEPATH, delete_dsn = TRUE)
