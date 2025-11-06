# Purpose: Take the strat input data frame created in previous step and run the stratificaiton
# functions, once for each of the 5 drone_pairing_tiers

library(tidyverse)

source("deploy/ground-ref-data/determine-dataset-withholding/src/greedy-forward-simple.R")


d = read_csv("/ofo-share/catalog-data-prep/stratification-data/for_strat.csv")
# Note that in the most recent run of the above, a random (unpaired) drone footprint (mission
# 000346) got lumped into Group 1 (the Emerald Point flights). Did not determine why but left it
#         included, not much impact.

d = d |>
  rename(group_id = drone_group_id, 
         pairing_tier = drone_pairing_tier,
         area_ha = plot_area_ha,
         n_trees = n_trees_live)


table(d$pairing_tier)


# Remove a fraction of some ground plots that are highly-clustered and very small (NEON and Eshom)
# so they don't skew the distribution we're aiming for, or use up too many of the plots we want to
# withhold. If the ones we keep get selected, the ones around them will get selected anyway because
# they are in the same drone groups
d_unfiltered = d

table(d$project_name)
overabundant_indexes = which(d$project_name %in% c("NEON2023", "Eshom2020"))
# select every fourth index
every_fourth = seq(1, length(overabundant_indexes), by = 4)
indexes_to_remove = overabundant_indexes[-every_fourth]
d = d[-indexes_to_remove, ]
table(d$project_name)


# Start out by stratifying only plots that have drone data
d_w_drone = d |>
  filter(pairing_tier %in% c("aligned_paired_drone_footprint", "non_paired_drone_footprint", 
                              "paired_drone_footprint"))




# # Run order: aligned_paired_drone_footprint, paired_drone_footprint, non_paired_drone_footprint,
# # no_drone_footprint, drone_no_plots (the latter is drone footprints with no plots in them, the
# # others are one record for each plot)





# Run all of the 3 plot-drone pairings together as an alternative approach

# Use the reference distribution from all drone-paired plots
res_joint <- select_withheld_groups(plots_df = d_w_drone, preselected_groups = c())

print_selection_report(res_joint)
create_factorial_plots(res_joint)

groups_withheld1 = unique(res_joint$withheld_plots$group_id)

save.image(file = "/ofo-share/scratch-derek/strat_workspace.RData")

# TODOs:
# Check what % of drone footprints this is
# Reduce % plots constraint
# Think about adding a groups % constraint

# Run for no_drone_footprint (plots without drone data)
# Run for drone_no_plots (drone footprints without plots)





## Check what percent of drone footprints this is
library(sf)

# Load the drone footprints with group IDs and determine which overlap any plots (so we use only
# those as the denominator, as those are the only ones we've given a chance to be selected)
all_footprints = st_read("/ofo-share/catalog-data-prep/stratification-data/all_drone_footprints_with_group_ids.gpkg")

footprints_with_plots = all_footprints |>
  filter(group_id %in% d_w_drone$group_id)

footprints_selected = footprints_with_plots |>
  filter(group_id %in% res_joint$withheld_plots$group_id)

selected_missions = unique(c(footprints_selected$mission_id_hn, footprints_selected$mission_id_lo))
all_missions = unique(c(footprints_with_plots$mission_id_hn, footprints_with_plots$mission_id_lo))
percent_missions_selected = length(selected_missions)/length(all_missions) * 100
percent_missions_selected

# Determine what percent of groups this is
percent_groups_selected = length(unique(res_joint$withheld_plots$group_id)) / length(unique(d_w_drone$group_id)) * 100
percent_groups_selected

# Make shapefiles of the selected and unselected footprints and plots. Just one layer for plots and
# one for footprints, with an attribute indicating selected vs unselected

all_footprints = all_footprints |>
  mutate(selected_for_withholding = group_id %in% res_joint$withheld_plots$group_id)


## Plots

# Load all field plots and determine if withheld based on whether they are physically overlapped by
# a withheld drone footprint (in which case, withhold them too)

ground_plots = st_read("/ofo-share/catalog-data-prep/stratification-data/downloaded-from-gdrive/ofo_ground-reference_plots.gpkg")

footprints_selected = all_footprints |>
  filter(selected_for_withholding == TRUE) |>
  st_union()

ground_plots = st_transform(ground_plots, st_crs(footprints_selected))

ground_plots_intersect = st_intersects(ground_plots, footprints_selected, sparse = FALSE)[, 1]
ground_plots$selected_for_withholding = ground_plots_intersect

st_write(ground_plots, "/ofo-share/catalog-data-prep/stratification-data/strat-output/ground_plots_with_selection_only-drone-ground-paired.gpkg", delete_dsn = TRUE)
st_write(all_footprints, "/ofo-share/catalog-data-prep/stratification-data/strat-output/drone_footprints_with_selection_only-drone-ground-paired.gpkg", delete_dsn = TRUE)


# Now stratify the plots_no_drone tier, passing it the list of groups already withheld above to
# force them to be withheld (want all overlaps withheld)

# Continuous variables to stratify
CONTINUOUS_VARS <- c("ppt", "trees_per_ha", "mean_ba_live", "area_ha")

# Categorical variables to stratify
CATEGORICAL_VARS <- c("ecoregion", "sp_comp_group", "project_name", "pairing_tier")

# Binning parameters
TARGET_PLOTS_PER_BIN <- 3 # default was 10
MAX_BINS <- 3 # was 5
N_BINS_FACTORIAL <- 3

# Selection parameters
MIN_PCT <- 15
MAX_PCT <- 22
MIN_GROUPS_PCT <- 15
MAX_GROUPS_PCT <- 22
FACTORIAL_WEIGHT <- 0.5

# Factorial combinations to track
FACTORIAL_COMBINATIONS <- list(
  list(var1 = "mean_ba_live", var2 = "sp_comp_group"),
  list(var1 = "trees_per_ha", var2 = "sp_comp_group"),
  list(var1 = "ppt", var2 = "sp_comp_group"),
  list(var1 = "mean_ba_live", var2 = "ecoregion"),
  list(var1 = "trees_per_ha", var2 = "ecoregion"),
  list(var1 = "ppt", var2 = "ecoregion"),
  list(var1 = "trees_per_ha", var2 = "mean_ba_live"),
  list(var1 = "pairing_tier", var2 = "sp_comp_group")
)

d_plot_no_drone = d |>
  filter(pairing_tier %in% c("no_drone_footprint"))

# Can only pass it the withheld groups that are present in the data we're asking it to stratify
groups_withheld_pre = groups_withheld1[groups_withheld1 %in% d_plot_no_drone$group_id]
# OK actually no groups overlap here, so it's an empty vector but that's OK

res_plot_no_drone <- select_withheld_groups(plots_df = d_plot_no_drone, 
                                              preselected_groups = groups_withheld_pre)
print_selection_report(res_plot_no_drone)

groups_withheld2 = unique(res_plot_no_drone$withheld_plots$group_id)
plots_withheld2 = unique(res_plot_no_drone$withheld_plots$plot_id)

groups_withheld_cum = unique(c(groups_withheld1, groups_withheld2))

# Now stratify the drone_no_plots tier, passing it the list of groups already withheld above (since
# ther will be some drone_no_plots footprints that are members of groups already being withheld
# because they overlap drone plots that do have ground plots and were selected for withholding)

# Continuous variables to stratify
CONTINUOUS_VARS <- c("ppt")

# Categorical variables to stratify
CATEGORICAL_VARS <- c("ecoregion", "project_name")

# Binning parameters
TARGET_PLOTS_PER_BIN <- 10
MAX_BINS <- 5
N_BINS_FACTORIAL <- 3

# Selection parameters
MIN_PCT <- 15
MAX_PCT <- 22
MIN_GROUPS_PCT <- 15
MAX_GROUPS_PCT <- 22
FACTORIAL_WEIGHT <- 0.5

# Factorial combinations to track
FACTORIAL_COMBINATIONS <- list(
  list(var1 = "ppt", var2 = "ecoregion"),
  list(var1 = "ppt", var2 = "project_name"),
  list(var1 = "ecoregion", var2 = "project_name")
)

d_drone_no_plot = d |>
  filter(pairing_tier %in% c("drone_no_plots")) |>
  # have to give each "plot" (drone footprint) an equal number of (dummy) trees for the algorithm to
  # work since it tries to match a % trees selection criterion
  mutate(n_trees = 10)

# Can only pass it the withheld groups that are present in the data we're asking it to stratify
groups_withheld_pre = groups_withheld_cum[groups_withheld_cum %in% d_drone_no_plot$group_id]
# Apparently there are also no drone missions without plots that are in any withheld groups.
# The way the to-stratify data frame was created, we excluded no-plot drone missions that were in
# groups containing other drone missions with plots, so this makes sense.

res_drone_no_plots <- select_withheld_groups(plots_df = d_drone_no_plot, 
                                              preselected_groups = groups_withheld_pre)

print_selection_report(res_drone_no_plots)

groups_withheld3 = unique(res_drone_no_plots$withheld_plots$group_id)
groups_withheld_final = unique(c(groups_withheld1, groups_withheld2, groups_withheld3)) |> sort()


# Get the plot_ids and mission_ids of the withheld missions and plots overall, and also get
# percentages selected

# Load the drone footprints with group IDs and determine which overlap any plots (so we use only
# those as the denominator, as those are the only ones we've given a chance to be selected)
all_footprints = st_read("/ofo-share/catalog-data-prep/stratification-data/all_drone_footprints_with_group_ids.gpkg")

footprints_selected = all_footprints |>
  filter(group_id %in% groups_withheld_final)

selected_missions = unique(c(footprints_selected$mission_id_hn, footprints_selected$mission_id_lo))
all_missions = unique(c(all_footprints$mission_id_hn, all_footprints$mission_id_lo))
percent_missions_selected = length(selected_missions)/length(all_missions) * 100
percent_missions_selected

# Determine what percent of groups this is
percent_groups_selected = length(unique(groups_withheld_final)) / length(unique(d$group_id)) * 100
percent_groups_selected

# Make shapefiles of the selected and unselected footprints and plots. Just one layer for plots and
# one for footprints, with an attribute indicating selected vs unselected

all_footprints = all_footprints |>
  mutate(selected_for_withholding = group_id %in% groups_withheld_final)


## Plots

# Load all field plots and determine if withheld based on whether they are physically overlapped by
# a withheld drone footprint (in which case, withhold them too)

ground_plots = st_read("/ofo-share/catalog-data-prep/stratification-data/downloaded-from-gdrive/ofo_ground-reference_plots.gpkg")

footprints_selected = all_footprints |>
  filter(selected_for_withholding == TRUE) |>
  st_union()

ground_plots = st_transform(ground_plots, st_crs(footprints_selected))

ground_plots_intersect = st_intersects(ground_plots, footprints_selected, sparse = FALSE)[, 1]
ground_plots$selected_for_withholding = ground_plots_intersect

# Also select the plots that we selected in Step 2, which is for plots that don't have drone
# footprints, so we can't use overlap with selected drone footprints to select them
ground_plots$selected_for_withholding = ifelse(ground_plots$plot_id %in% plots_withheld2, TRUE, ground_plots$selected_for_withholding)

st_write(ground_plots, "/ofo-share/catalog-data-prep/stratification-data/strat-output/ground_plots_with_selection.gpkg", delete_dsn = TRUE)
st_write(all_footprints, "/ofo-share/catalog-data-prep/stratification-data/strat-output/drone_footprints_with_selection.gpkg", delete_dsn = TRUE)

# Get the list of the plot and drone mission IDs
withheld_plot_ids = ground_plots |>
  filter(selected_for_withholding == TRUE) |>
  st_set_geometry(NULL) |>
  pull(plot_id) |>
  sort()

footprints_selected = all_footprints |>
  filter(selected_for_withholding == TRUE)
withheld_mission_ids = unique(c(footprints_selected$mission_id_hn, footprints_selected$mission_id_lo)) |>
  sort()

# Save to CSVs (two single-column files)
write_csv(tibble(plot_id = withheld_plot_ids),
          "/ofo-share/catalog-data-prep/stratification-data/strat-output/withheld_ground_plot_ids.csv")
write_csv(tibble(mission_id = withheld_mission_ids),
          "/ofo-share/catalog-data-prep/stratification-data/strat-output/withheld_drone_mission_ids.csv")
