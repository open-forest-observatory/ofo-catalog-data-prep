# Purpose: Compute crown area (sq m) for each detected tree crown and merge it into the detected
# trees with DBH GPKG.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
ALL_DETECTED_CROWNS_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-tree-crowns.gpkg")
ALL_DETECTED_TREES_WITH_DBH_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-trees-with-dbh.gpkg")


# --- Load data ---

all_crowns = st_read(ALL_DETECTED_CROWNS_FILEPATH, quiet = TRUE)
all_trees = st_read(ALL_DETECTED_TREES_WITH_DBH_FILEPATH, quiet = TRUE)


# --- Pre-process crowns: remove duplicates ---

# Remove duplicated crown IDs per composite_id (take the first), matching the approach in script 06
# for treetop points. Crowns link to treetops via treetop_unique_ID.
dups = all_crowns |>
  sf::st_drop_geometry() |>
  dplyr::group_by(composite_id, treetop_unique_ID) |>
  dplyr::filter(dplyr::n() > 1) |>
  dplyr::summarise(n_removed = dplyr::n() - 1, .groups = "drop")

if (nrow(dups) > 0) {
  msgs = dplyr::mutate(dups, msg = paste0("  composite_id=", composite_id, ", treetop_unique_ID=", treetop_unique_ID, ": ", n_removed, " duplicate(s) removed"))
  warning(paste(c("Duplicated crown treetop_unique_IDs found per composite (keeping first):", msgs$msg), collapse = "\n"))
  all_crowns = all_crowns |>
    dplyr::group_by(composite_id, treetop_unique_ID) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
}


# --- Compute crown area ---

# Project each crown to its local UTM zone to compute area in sq m
all_crowns$crown_area_sqm = as.numeric(st_area(all_crowns))


# --- Merge crown area into trees ---

# Build a lookup table of crown area keyed by composite_id + treetop_unique_ID
crown_areas = all_crowns |>
  sf::st_drop_geometry() |>
  dplyr::select(composite_id, treetop_unique_ID, crown_area_sqm)

# The trees file uses unique_ID which corresponds to the crowns' treetop_unique_ID
all_trees = all_trees |>
  dplyr::left_join(crown_areas, by = c("composite_id" = "composite_id", "unique_ID" = "treetop_unique_ID"))

n_matched = sum(!is.na(all_trees$crown_area_sqm))
cat("Matched crown area for", n_matched, "of", nrow(all_trees), "trees\n")


# --- Save ---

st_write(all_trees, ALL_DETECTED_TREES_WITH_DBH_FILEPATH, delete_dsn = TRUE)
cat("Updated detected trees GPKG with crown area\n")
