# Purpose: Summarize the ground reference data in various ways to inform modeling decisions. Not part of a pipeline.

library(sf)
library(tidyverse)

GROUND_REF_PLOTS_OUTPUT_GPKG = file.path(REPO_DATA_LOCAL_PATH, "ofo_ground-reference_plots.gpkg")
GROUND_REF_TREES_OUT_GPKG = file.path(REPO_DATA_LOCAL_PATH, "ofo_ground-reference_trees.gpkg")