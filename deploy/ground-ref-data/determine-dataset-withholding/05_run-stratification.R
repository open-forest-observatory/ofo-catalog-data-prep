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

res1 = select_withheld_groups(plots_df, required_groups = c()) # The STEF 2018 standalone one is group 69

print_selection_report(res1)

print(res1$diagnostics$factorial_plots)

groups_withheld1 = res1$withheld_plots$group_id

# Prep the next run: "paired_drone_footprint"

plots_df = d |>
  filter(pairing_tier == "paired_drone_footprint")

res2 = select_withheld_groups(plots_df, required_groups = groups_withheld1)

print_selection_report(res2)
print(res2$diagnostics$factorial_plots)

# Combine the withheld groups from both runs
groups_withheld2 = unique(c(groups_withheld1, res2$withheld_plots$group_id))

# Next run: "non_paired_drone_footprint"
plots_df = d |>
  filter(pairing_tier == "non_paired_drone_footprint")

res3 = select_withheld_groups(plots_df, required_groups = groups_withheld2)

print_selection_report(res3)

print(res3$diagnostics$factorial_plots)


# Combine the withheld groups from all three runs
groups_withheld3 = unique(c(groups_withheld2, res3$withheld_plots$group_id))

# Next run: "no_drone_footprint"
plots_df = d |>
  filter(pairing_tier == "no_drone_footprint")

res4 = select_withheld_groups(plots_df, required_groups = groups_withheld3)

print_selection_report(res4)
print(res4$diagnostics$factorial_plots)

# Combine the withheld groups from all four runs
groups_withheld4 = unique(c(groups_withheld3, res4$withheld_plots$group_id))





plots_withheld = plots_df |>
  filter(group_id %in% groups_withheld)

table(plots_withheld$project_name)
table(plots_df$project_name)






# Figuring out why we got very divergent mean and SD for mean_ba_live for this selection:

withheld_plots = plots_df |>
  filter(group_id %in% res$withheld_plots$group_id)

hist(withheld_plots$mean_ba_live)
hist(plots_df$mean_ba_live)
