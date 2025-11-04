# Purpose: Take the strat input data frame created in previous step and run the stratificaiton
# functions, once for each of the 5 drone_pairing_tiers

library(tidyverse)

source("deploy/ground-ref-data/determine-dataset-withholding/src/greedy-forward-simple.R")


d = read_csv("/ofo-share/catalog-data-prep/stratification-data/for_strat.csv")

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
res_joint <- select_withheld_groups(plots_df = d_w_drone)

print_selection_report(res_joint)
create_factorial_plots(res_joint)






print(res_joint$diagnostics$factorial_plots)

groups_withheld1 = unique(res_joint$withheld_plots$group_id)

save.image(file = "/ofo-share/scratch-derek/strat_workspace.RData")

# TODOs:
# Check what % of drone footprints this is
# Reduce % plots constraint
# Think about adding a groups % constraint

# Run for no_drone_footprint (plots without drone data)
# Run for drone_no_plots (drone footprints without plots)





# ## Check what percent of drone footprints this is

# # Load the drone footprints with group IDs and determine which overlap any plots (so we use only
# # those as the denominator, as those are the only ones we've given a chance to be selected)
# all_footprints = st_read("/ofo-share/catalog-data-prep/stratification-data/all_drone_footprints_with_group_ids.gpkg")

# footprints_with_plots = all_footprints |>
#   filter(group_id %in% d_w_drone$group_id)

# footprints_selected = footprints_with_plots |>
#   filter(group_id %in% res_joint$withheld_plots$group_id)

# percent_footprints_with_plots_selected = nrow(footprints_selected) / nrow(footprints_with_plots) * 100
# percent_footprints_with_plots_selected

# # Determine what percent of groups this is
# percent_groups_selected = length(res_joint$withheld_plots$group_id) / length(unique(d_w_drone$group_id)) * 100
# percent_groups_selected
