# Purpose: Estimate DBH for each detected tree using cloud2trees and add as a column to the
# detected trees GPKG.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)
library(cloud2trees)
library(furrr)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
ALL_DETECTED_TREES_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-trees.gpkg")


# --- Load data ---

all_trees = st_read(ALL_DETECTED_TREES_FILEPATH, quiet = TRUE)


# --- Pre-process: split trees by plot ---

# Process each plot separately because the function derives the allometry for the local area of the
# plot

all_trees_by_plot = split(all_trees, all_trees$composite_id)


# --- Estimate DBH using cloud2trees ---

add_dbh = function(all_trees_foc) {
  cat("Estimating DBH for plot with", nrow(all_trees_foc), "trees\n")

  all_trees_foc_conformed = all_trees_foc |>
    select(treeID = unique_ID, tree_height_m = height, species_prediction, live_dead_prediction)

  res = cloud2trees::trees_dbh(all_trees_foc_conformed)
  all_trees_foc$dbh = res$fia_est_dbh_cm

  all_trees_foc
}

plan(multisession, workers = availableCores()/4)
res = future_map(all_trees_by_plot, add_dbh, .progress = TRUE, .options = furrr_options(seed = TRUE, scheduling = Inf))




# --- Save updated detected trees with DBH ---

st_write(all_trees, ALL_DETECTED_TREES_FILEPATH, delete_dsn = TRUE)
cat("Updated detected trees GPKG with DBH\n")
