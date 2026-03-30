# Purpose: Estimate DBH for each detected tree using cloud2trees and add as a column to the
# detected trees GPKG.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)
library(cloud2trees)
library(furrr)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
ALL_DETECTED_TREES_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-trees.gpkg")
ALL_DETECTED_TREES_WITH_DBH_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-trees-with-dbh.gpkg")


# --- Load data ---

all_trees = st_read(ALL_DETECTED_TREES_FILEPATH, quiet = TRUE)


# --- Pre-process

# Remove duplicated tree IDs per composite_id (take the first), but providing a warning that says,
# for each duplicated tree ID per plot, its ID and how many duplicates were removed.
dups = all_trees |>
  sf::st_drop_geometry() |>
  dplyr::group_by(composite_id, unique_ID) |>
  dplyr::filter(dplyr::n() > 1) |>
  dplyr::summarise(n_removed = dplyr::n() - 1, .groups = "drop")

if (nrow(dups) > 0) {
  msgs = dplyr::mutate(dups, msg = paste0("  composite_id=", composite_id, ", unique_ID=", unique_ID, ": ", n_removed, " duplicate(s) removed"))
  warning(paste(c("Duplicated tree IDs found per plot (keeping first):", msgs$msg), collapse = "\n"))
  all_trees = all_trees |>
    dplyr::group_by(composite_id, unique_ID) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
}




# Split trees by plot

# Process each plot separately because the function derives the allometry for the local area of the
# plot

all_trees_by_plot = split(all_trees, all_trees$composite_id)


# --- Estimate DBH using cloud2trees ---

add_dbh = function(all_trees_foc) {
  # Re-establish sf geometry attribute after parallel serialization
  all_trees_foc = sf::st_as_sf(all_trees_foc)
  cat("Estimating DBH for plot with", nrow(all_trees_foc), "trees\n")

  all_trees_foc_conformed = all_trees_foc |>
    dplyr::select(treeID = unique_ID, tree_height_m = height, species_prediction, live_dead_prediction) |>
    sf::st_as_sf()

  res = cloud2trees::trees_dbh(all_trees_foc_conformed)
  all_trees_foc$dbh = res$fia_est_dbh_cm

  all_trees_foc
}

plan(multisession, workers = availableCores()/4)
res = future_map(all_trees_by_plot, add_dbh, .progress = TRUE, .options = furrr_options(seed = TRUE, scheduling = Inf))
# res = map(all_trees_by_plot, add_dbh, .progress = TRUE)

# Combine results
all_trees_with_dbh = bind_rows(res)


# --- Save updated detected trees with DBH ---

st_write(all_trees_with_dbh, ALL_DETECTED_TREES_WITH_DBH_FILEPATH, delete_dsn = TRUE)
cat("Updated detected trees GPKG with DBH\n")

