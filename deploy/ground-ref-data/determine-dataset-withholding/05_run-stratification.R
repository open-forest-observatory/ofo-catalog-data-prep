# Purpose: Take the strat input data frame created in previous step and run the stratificaiton
# functions, once for each of the 5 drone_pairing_tiers

library(tidyverse)

source("deploy/ground-ref-data/determine-dataset-withholding/src/proportional_sampling_simplified.R")


d = read_csv("/ofo-share/catalog-data-prep/stratification-data/for_strat.csv")

d = d |>
  rename(group_id = drone_group_id, 
         pairing_tier = drone_pairing_tier,
         area_ha = plot_area_ha,
         n_trees = n_trees_live)


table(d$pairing_tier)

# Run order: aligned_paired_drone_footprint, paired_drone_footprint, non_paired_drone_footprint,
# no_drone_footprint, drone_no_plots (the latter is drone footprints with no plots in them, the
# others are one record for each plot)

# Set cols to names expected by stratification function

# Start with aligned_paired_drone_footprint

plots_df = d |>
  filter(pairing_tier == "aligned_paired_drone_footprint")

res = select_withheld_groups(plots_df)

print_selection_report(res)

print(res$diagnostics$factorial_plots)


# Get the plots in the withheld groups

groups_withheld = res$withheld_plots$group_id

plots_withheld = plots_df |>
  filter(group_id %in% groups_withheld)

table(plots_withheld$project_name)
table(plots_df$project_name)






# Figuring out why we got very divergent mean and SD for mean_ba_live for this selection:

withheld_plots = plots_df |>
  filter(group_id %in% res$withheld_plots$group_id)

hist(withheld_plots$mean_ba_live)
hist(plots_df$mean_ba_live)
