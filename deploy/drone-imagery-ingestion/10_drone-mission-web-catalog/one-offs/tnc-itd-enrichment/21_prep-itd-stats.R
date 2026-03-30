# Read in the already-generated ITD stats CSV, filter to focal Yuba plots, and summarize/visualize
# TODO: Clarify the difference between the two stats data frames saved

library(tidyverse)
library(sf)


stats = read_csv("/ofo-share/repos/david/tree-detection-parameterization/data/scores_full_grid_correct_heights.csv")
yuba_missions_list = read_csv("/ofo-share/project-data/tnc-yuba-deliverables/yuba-drone-missions.csv")
ground_plots = st_read("/ofo-share/project-data/tnc-yuba-deliverables/ground-reference/ofo_ground-reference_plots.gpkg")
ground_trees = st_read("/ofo-share/project-data/tnc-yuba-deliverables/ground-reference/ofo_ground-reference_trees.gpkg")

# Ground plot 0108 has a tree that is impossibly large DBH (894.08), must be typo, this is 352
# inches, probalby should be 35.2 inches, which is 89.408, so make this correction
ground_trees[ground_trees$plot_id == "0108" & ground_trees$dbh == 894.08, "dbh"] = 89.408

# Plot 42 tree with DBH 82.296 also appears erroneous because its height is very short for that
# diameter and there are no defects recorded, so remove it
ground_trees = ground_trees |> filter(!(plot_id == "0042" & dbh == 82.296))


yuba_missions = yuba_missions_list |>
  filter(group == "yuba") |>
  pull(mission_id)

# Split the stats 'dataset' column into ground_plot_id, nadir_mission_id, and ground_mission_id by
# underscores

stats_parsed = stats |>
  separate(dataset, into = c("ground_plot_id", "nadir_mission_id", "oblique_mission_id"), sep = "_") |>
  filter(nadir_mission_id %in% yuba_missions | oblique_mission_id %in% yuba_missions)

# Filter the param combo to the best one, which was selected for the ITD run:
# "sigma_0.3333333333333333__b_0.025_c_0.0625"

stats_filtered = stats_parsed |>
  filter(param_combo == "sigma_0.3333333333333333__b_0.025_c_0.0625") |>
  select(ground_plot_id, nadir_mission_id, oblique_mission_id, F1, precision, recall, n_core_field, n_core_drone)


# Next, pull in the ground plot trees and bounds, filter to trees > 10 m ht, compute density and
# mean height, then visualize F scores vs density and mean height

ground_trees_foc = ground_trees |>
  filter(height > 10, live_dead == "L")

# Compute plot-level metrics

ground_trees_summ = ground_trees_foc |>
  st_drop_geometry() |>
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

# Join the stats onto the summary metrics
stats_w_ground = ground_plots_summ |>
  left_join(stats_filtered, by = c("plot_id" = "ground_plot_id"))

# Filter to yuba missions
stats_w_ground = stats_w_ground |>
  filter(nadir_mission_id %in% yuba_missions) |>
  filter(min_ht <= 10)

# Filter to missions with quality shifts (i.e. we were able to identify the correct shift
# unambiguously). This means that we are only looking at missions that were evaluated in the species
# project, even though there are more TNC missions than that.
shift_quality = read_csv("/ofo-share/project-data/tnc-yuba-deliverables/ground-reference/shift_quality.csv")
quality_shifts = shift_quality |>
  filter(Quality >= 3) |>
  pull(Dataset) |>
  str_remove(".tif")

stats_w_ground = stats_w_ground |>
  mutate(dataset_id = paste(plot_id, nadir_mission_id, oblique_mission_id, sep = "_")) |>
  filter(dataset_id %in% quality_shifts)


# For inspection of the focal trees of each plot to understand poor alignment etc
st_write(ground_trees_foc, "/ofo-share/project-data/tnc-yuba-deliverables/ground-reference/ofo_ground-reference_trees_foc.gpkg", delete_dsn = TRUE)

# Remove 3 obviously poorly aligned plots (based on manual inspection, this is all pots with F1 < 0.51)
stats_w_ground_foc = stats_w_ground |>
  filter(F1 >= 0.51)

# Include only one row per plot_id, keeping the highest F1 score (best-aligned mission pair)
stats_w_ground_foc = stats_w_ground_foc |>
  group_by(plot_id) |>
  slice_max(F1, n = 1) |>
  ungroup()

# Visualize F1 vs tree density and mean height
p = ggplot(stats_w_ground_foc, aes(x = tph, y = F1, color = mean_height)) +
  geom_point() +
  labs(x = "Trees per hectare", y = "F1 score") +
  scale_color_viridis_c(name = "Mean tree height (m)") +
  theme_bw(20)

png("/ofo-share/project-data/tnc-yuba-deliverables/figures/itd_f1_vs_density.png", width = 800, height = 600)
print(p)
dev.off()

# Save out the full table
write_csv(stats_w_ground_foc, "/ofo-share/project-data/tnc-yuba-deliverables/tables/itd_stats_with_ground_metrics.csv")

# Get the mean F score, precision, and recall
mean_scores = tibble(
  mean_F1 = mean(stats_w_ground_foc$F1, na.rm = TRUE),
  mean_precision = mean(stats_w_ground_foc$precision, na.rm = TRUE),
  mean_recall = mean(stats_w_ground_foc$recall, na.rm = TRUE)
)
mean_scores

# Next, shift ground plot bounds with shifts specified by shifts_per_dataset.json, get the drone
# trees within the shifted bounds, compute the drone-based density and BA, and compare to the ground
# plot values. Plot on a scatterplot with points colored by mean tree height and with a 1:1 line

library(jsonlite)

# Read shifts and drone trees
shifts = fromJSON("/ofo-share/project-data/tnc-yuba-deliverables/ground-reference/shift_per_dataset.json")
drone_trees = st_read("/ofo-share/project-data/tnc-yuba-deliverables/composite-missions/overall/all-detected-trees_analysis-ready.gpkg", quiet = TRUE)

# Re-read ground plots with geometry for shifting
ground_plots_geom = st_read("/ofo-share/project-data/tnc-yuba-deliverables/ground-reference/ofo_ground-reference_plots.gpkg", quiet = TRUE) |>
  select(plot_id, geom) |>
  st_transform(32610)  # match drone trees CRS

# Filter drone trees to live and height > 10 m (matching ground tree criteria)
drone_trees_foc = drone_trees |>
  filter(live_dead_prediction == "Live", height > 10)

# For each dataset, shift plot bounds, clip drone trees, compute metrics
drone_metrics = stats_w_ground_foc |>
  rowwise() |>
  mutate(
    composite_id = paste0(nadir_mission_id, "_", oblique_mission_id)
  ) |>
  ungroup()

compute_drone_metrics = function(dataset_id, plot_id, nadir_mission_id, oblique_mission_id, ...) {
  # Get shift for this dataset
  shift_xy = shifts[[dataset_id]]
  if (is.null(shift_xy)) {
    return(tibble(dataset_id = dataset_id, drone_tree_count = NA, drone_bah = NA, drone_tph = NA))
  }
  shift_x = shift_xy[1, 1]
  shift_y = shift_xy[1, 2]

  # Get ground plot geometry and shift it
  plot_geom = ground_plots_geom |> filter(plot_id == !!plot_id)
  if (nrow(plot_geom) == 0) {
    return(tibble(dataset_id = dataset_id, drone_tree_count = NA, drone_bah = NA, drone_tph = NA))
  }
  shifted_geom = plot_geom |> st_geometry() + c(shift_x, shift_y)
  shifted_geom = st_sfc(shifted_geom, crs = 32610)
  shifted_plot = st_sf(geometry = shifted_geom)

  # Save out the shifted plot to /ofo-share/project-data/tnc-yuba-deliverables/tmp/shifted_plots/shifted_plot_{dataset_id}.gpkg for inspection
  st_write(shifted_plot, file.path("/ofo-share/project-data/tnc-yuba-deliverables/tmp/shifted_plots", paste0("shifted_plot_", dataset_id, ".gpkg")), delete_dsn = TRUE)

  # Get drone trees for this composite within the shifted bounds
  comp_id = paste0(nadir_mission_id, "_", oblique_mission_id)
  drone_foc = drone_trees_foc |> filter(composite_id == comp_id)
  if (nrow(drone_foc) == 0) {
    return(tibble(dataset_id = dataset_id, drone_tree_count = 0, drone_bah = 0, drone_tph = 0))
  }

  # Clip drone trees to shifted plot bounds
  drone_in_plot = st_intersection(drone_foc, shifted_plot)

  # Compute metrics
  n_drone = nrow(drone_in_plot)
  plot_area_ha = as.numeric(st_area(shifted_plot)) / 10000
  drone_ba_sqm = sum(pi * (drone_in_plot$dbh / 100 / 2)^2, na.rm = TRUE)

  tibble(
    dataset_id = dataset_id,
    drone_tree_count = n_drone,
    drone_tph = n_drone / plot_area_ha,
    drone_bah = drone_ba_sqm / plot_area_ha
  )
}

results = pmap_dfr(drone_metrics, compute_drone_metrics)

# Join drone metrics onto ground metrics
stats_comparison = drone_metrics |>
  left_join(results, by = "dataset_id") |>
  filter(!is.na(drone_tree_count))

# Remove pairings with zero drone trees (removed due to poor shift, or not flown)
stats_comparison = stats_comparison |>
  filter(drone_tree_count > 0)

# Exclude actual poor alignment (bad shift) plot
stats_comparison = stats_comparison |>
  filter(dataset_id != "0091_001415_001414")

# Include only one row for each plot_id (there are multiple if there are multiple drone mission
# pairs overlapping it). Keep the hightest F1 score of them, since that should correspond to the
# best-aligned mission pair with best photogrammetry
stats_comparison = stats_comparison |>
  group_by(plot_id) |>
  slice_max(F1, n = 1) |>
  ungroup()

# # Scatterplot: ground vs drone tree density, colored by mean height, with 1:1 line
# p_density = ggplot(stats_comparison, aes(x = tph, y = drone_tph, color = mean_height)) +
#   geom_point(size = 3) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
#   labs(x = "Ground trees per hectare", y = "Drone trees per hectare", title = "Tree density: ground vs drone") +
#   scale_color_viridis_c(name = "Mean tree\nheight (m)") +
#   theme_bw(16)

# png("/ofo-share/project-data/tnc-yuba-deliverables/figures/ground_vs_drone_density.png", width = 800, height = 700)
# print(p_density)
# dev.off()

# Compute R2
r2_ba = cor(stats_comparison$bah, stats_comparison$drone_bah)^2
r2_ba

# Scatterplot: ground vs drone basal area per hectare, colored by mean height, with 1:1 line
p_ba = ggplot(stats_comparison, aes(x = bah, y = drone_bah)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  annotate("text", x = min(stats_comparison$bah), y = max(stats_comparison$drone_bah),
           label = sprintf("R² = %.2f", r2_ba), hjust = 0, vjust = 1, size = 8) +
  labs(x = "Ground BA (sq m / ha)", y = "Drone BA (sq m / ha)") +
  # scale_color_viridis_c(name = "F1 score") +
  theme_bw(20)

png("/ofo-share/project-data/tnc-yuba-deliverables/figures/ground_vs_drone_ba.png", width = 800, height = 700)
print(p_ba)
dev.off()

# Save comparison table
write_csv(stats_comparison, "/ofo-share/project-data/tnc-yuba-deliverables/tables/itd_ground_vs_drone_comparison.csv")

# Get the mean F score, precision, and recall
mean_scores = tibble(
  mean_F1 = mean(stats_comparison$F1, na.rm = TRUE),
  mean_precision = mean(stats_comparison$precision, na.rm = TRUE),
  mean_recall = mean(stats_comparison$recall, na.rm = TRUE)
)
mean_scores

