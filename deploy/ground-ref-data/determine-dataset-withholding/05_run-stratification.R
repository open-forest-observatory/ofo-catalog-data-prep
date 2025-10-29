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


# Prepare reference from ALL data

# UNFINISHED:
# # Average the merics across all plots in a group so that groups with tons of plots (e.g. NEON) don't
# # sway the average, and the average is more reflective of the drone mision coverage than the plot
# # abundance

# d_groupmeans = d |>
#   group_by(group_id) |>
#   summarise(
#     mean_ba_live = mean(mean_ba_live, na.rm = TRUE),
#     ppt = mean(ppt, na.rm = TRUE),
#     trees_per_ha = mean(trees_per_ha, na.rm = TRUE),
#     ba_sqm_per_ha = mean(ba_sqm_per_ha, na.rm = TRUE),
#     ecoregion = first(ecoregion),
#     n_plots = n(),
#     first_plot_id = first(plot_id)

#   ) |>
#   ungroup()

# # Create a reference distribution data frame with just the plot-level attributes that are truly only
# # relevant to ground plot data
# d_plots = d |>
#   select(plot_id, pairing_tier, project_name, area_ha, n_trees,  )

# The reference distribution will be based only on plots that have drone data
d_w_drone = d |>
  filter(pairing_tier %in% c("aligned_paired_drone_footprint", "non_paired_drone_footprint", 
                              "paired_drone_footprint"))
ref <- prepare_reference_distributions(d_w_drone)




# Run order: aligned_paired_drone_footprint, paired_drone_footprint, non_paired_drone_footprint,
# no_drone_footprint, drone_no_plots (the latter is drone footprints with no plots in them, the
# others are one record for each plot)

# Set cols to names expected by stratification function

# Start with aligned_paired_drone_footprint

plots_df = d |>
  filter(pairing_tier == "aligned_paired_drone_footprint")

# Original using only within-tier reference
# res1 = select_withheld_groups(plots_df, required_groups = c()) # The STEF 2018 standalone one is group 69

# Use the reference distribution from all drone-paired plots
res1 <- select_withheld_groups(
  plots_df = plots_df #,
  # reference_distribution = ref$reference_distribution,
  # reference_factorial = ref$reference_factorial,
  # reference_plots_df = ref$reference_plots_df
)

print_selection_report(res1)

print(res1$diagnostics$factorial_plots)

groups_withheld1 = unique(res1$withheld_plots$group_id)

tier1_selected_plots <- ref$reference_plots_df |>
  filter(group_id %in% res1$withheld_group_ids)





# Prep the next run: "paired_drone_footprint"

plots_df = d |>
  filter(pairing_tier == "paired_drone_footprint")


res2 <- select_withheld_groups(
  plots_df = plots_df,
  # reference_distribution = ref$reference_distribution,
  # reference_factorial = ref$reference_factorial,
  # reference_plots_df = ref$reference_plots_df,
  # previous_selected_plots = tier1_selected_plots,
  required_groups = groups_withheld1
)

print_selection_report(res2)
print(res2$diagnostics$factorial_plots)

# Combine the withheld groups from both runs
groups_withheld2 = unique(c(groups_withheld1, res2$withheld_plots$group_id))

tier2_selected_plots <- ref$reference_plots_df |>
  filter(group_id %in% c(res1$withheld_group_ids, res2$withheld_group_ids))

# Next run: "non_paired_drone_footprint"
plots_df = d |>
  filter(pairing_tier == "non_paired_drone_footprint")

res3 <- select_withheld_groups(
  plots_df = plots_df,
  # reference_distribution = ref$reference_distribution,
  # reference_factorial = ref$reference_factorial,
  # reference_plots_df = ref$reference_plots_df,
  # previous_selected_plots = tier2_selected_plots,
  required_groups = groups_withheld2
)


print_selection_report(res3)

print(res3$diagnostics$factorial_plots)


# Combine the withheld groups from all three runs
groups_withheld3 = unique(c(groups_withheld2, res3$withheld_plots$group_id))


# Check how well these selected goroups across 3 tiers match the ref distrib

d_relevant = d |>
  filter(pairing_tier %in% c("aligned_paired_drone_footprint", "paired_drone_footprint",
                              "non_paired_drone_footprint"))

combined_eval = evaluate_combined_selection(d_relevant, groups_withheld3, ref)

print_combined_evaluation_report(combined_eval)

# Print combined factorial plots
print(combined_eval$diagnostics$factorial_plots)














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


# When we get to stratifying drone_no_plots and plots_no_drone, use only that tier as a reference