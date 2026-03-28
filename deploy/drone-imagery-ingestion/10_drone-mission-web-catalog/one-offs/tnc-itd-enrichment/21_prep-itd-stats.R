# Read in the already-generated ITD stats CSV, filter to focal Yuba plots, and summarize/visualize

library(tidyverse)
library(sf)


stats = read_csv("/ofo-share/repos/david/tree-detection-parameterization/data/scores_full_grid_correct_heights.csv")
yuba_missions_list = read_csv("/ofo-share/project-data/tnc-yuba-deliverables/yuba-drone-missions.csv")
ground_plots = st_read("/ofo-share/project-data/tnc-yuba-deliverables/ground-reference/ofo_ground-reference_plots.gpkg")
ground_trees = st_read("/ofo-share/project-data/tnc-yuba-deliverables/ground-reference/ofo_ground-reference_trees.gpkg")



yuba_missions = yuba_missions_list |>
  filter(group == "yuba") |>
  pull(mission_id)

# Split the stats 'dataset' column into ground_plot_id, nadir_mission_id, and ground_mission_id by
# underscores

stats_parsed = stats |>
  separate(dataset, into = c("ground_plot_id", "nadir_mission_id", "ground_mission_id"), sep = "_") |>
  filter(nadir_mission_id %in% yuba_missions | ground_mission_id %in% yuba_missions)

# Filter the param combo to the best one, which was selected for the ITD run:
# "sigma_0.3333333333333333__b_0.025_c_0.0625"

stats_filtered = stats_parsed |>
  filter(param_combo == "sigma_0.3333333333333333__b_0.025_c_0.0625") |>
  select(ground_plot_id, nadir_mission_id, ground_mission_id, F1, precision, recall, n_core_field, n_core_drone)


# Next, pull in the ground plot trees and bounds, filter to trees > 10 m ht, compute density and
# mean height, then visualize F scores vs density and mean height

ground_trees = ground_trees |>
  filter(height > 10, live_dead == "L") |>
  st_drop_geometry()

# Compute plot-level metrics

ground_trees_summ = ground_trees |>
  group_by(plot_id) |>
  summarize(
    tree_count = n(),
    basal_area_sqm = sum(pi * (dbh / 100 / 2)^2),
    mean_height = mean(height)
  )

# Compute area for ground plot bounds

ground_plots = ground_plots |>
  st_drop_geometry() |>
  select(plot_id, plot_area_ha, num_ohvis_trees_excluded, min_ht)

# Join the tree metrics with the ground plot metrics and compute density and basal area per hectare
ground_plots_summ = ground_plots |>
  left_join(ground_trees_summ, by = "plot_id") |>
  mutate(
    tph = tree_count / plot_area_ha,
    bah = basal_area_sqm / plot_area_ha
  )

# Join the summary metrics onto the stats
stats_w_ground = stats_filtered |>
  left_join(ground_plots_summ, by = c("ground_plot_id" = "plot_id"))

# Filter to yuba missions
stats_w_ground = stats_w_ground |>
  filter(nadir_mission_id %in% yuba_missions) |>
  filter(min_ht <= 10)

# Visualize F1 vs tree density and mean height

p = ggplot(stats_w_ground, aes(x = tph, y = F1, color = mean_height)) +
  geom_point() +
  labs(x = "Trees per hectare", y = "F1 score") +
  scale_color_viridis_c(name = "Mean tree height (m)") +
  theme_bw(16)

png("/ofo-share/project-data/tnc-yuba-deliverables/figures/itd_f1_vs_density.png", width = 800, height = 600)
print(p)
dev.off()

# RESUME BY LOOKING IN QGIS FOR BAD ALIGNMENT FOR THE BAD F SCORE PLOTS, POSSIBLY FIRST SAVING OUT
# THE FILERED TREE DATA


# Next, shift ground plot bounds with shifts specified by shifts_per_dataset.json, get the drone
# trees within the shifted bounds, compute the drone-based density and BA, and compare to the ground
# plot values. Plot on a scatterplot with points colored by mean tree height? and with a 1:1 line