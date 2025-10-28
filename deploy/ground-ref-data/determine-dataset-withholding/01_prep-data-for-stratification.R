# Purpose: Select a stratified random subset of paired datasets (drone-high, drone-low, ground) to
# exclude from model training. Stratify across tree density, basal area/height, contributing
# project, survey year, (forest type?). About 20% of plots and trees.

library(tidyverse)
library(sf)
library(terra)


# Function to get only the mission_id column from a footprint layer
get_mission_id = function(footprint_layer) {
  footprint_layer |>
    select("mission_id", "project_name")
}

# Read predictor layers to stratify across
ppt = rast("/ofo-share/catalog-data-prep/stratification-data/data-layers-for-strat/prism_ppt_us_30s_2020_avg_30y.tif")
ecoregion = st_read("/ofo-share/catalog-data-prep/stratification-data/data-layers-for-strat/epa-ecoregions-l3.gpkg")

# Read the ground plots and trees
ground_plots = st_read("/ofo-share/catalog-data-prep/stratification-data/downloaded-from-gdrive/ofo_ground-reference_plots.gpkg")
ground_trees = st_read("/ofo-share/catalog-data-prep/stratification-data/downloaded-from-gdrive/ofo_ground-reference_trees.gpkg") |>
  st_drop_geometry()

# For WADNR plots we'll assume they're all live for purposes of species grouping (not specified in
# the data). Plots 186-196.
ground_trees = ground_trees |>
  mutate(live_dead = ifelse(as.numeric(plot_id) %in% c(186:196), "L", live_dead))

# Remove 3 trees that are impossibly large and must be erroneous
ground_trees = ground_trees |>
  filter(is.na(height) | is.na(dbh) | (!(height < 50 & dbh > 300) & !(height > 80)))

# Read the pairings definitions
pairings = read_csv("/ofo-share/catalog-data-prep/stratification-data/ground_plot_drone_mission_matches.csv")

# Read the selected HN-LO pairings from Amritha's model input file (this is a subset that includes
# only the pairings where a plot was successfully aligned to both the HN and LO drone data). We are
# only using this to classify plots as being in a paired drone footprint or not.
aligned_pairings = st_read("/ofo-share/scratch-amritha/tree-species-scratch/model_input_plots_with_split.gpkg") |>
  st_drop_geometry() |>
  select(plot_id, mission_id_hn, mission_id_lo) |>
  distinct()

# Get the unique drone mission pairings (dropping those that are replicated in the pairings because they cover
# multiple plots)
pairings_unique = pairings |>
  select(mission_id_hn, mission_id_lo) |>
  distinct()


# Manually exclude 000934 and 000935 from being selected for test (force them in train) because they connect nearly all STEF missions into one giant group
# and we want to allow a part of them to be selected as test. Excluding them will disconnect the
# chain of drone missions and allow part of STEF2018 to be selected (which we will manually force).
# pairings_unique = pairings_unique |>
#   filter(!(mission_id_hn %in% c("000934", "000935") | mission_id_lo %in% c("000934", "000935")))


# Read the drone mission footrpints and merge to single layer with one feature per mission
footprint_paths = list.files("/ofo-share/catalog-data-prep/01_raw-imagery-ingestion/metadata/3_final/1_full-metadata-per-mission", full.names = TRUE)
footprints = lapply(footprint_paths, st_read)
footrpints = lapply(footprints, get_mission_id)
footprints = do.call(rbind, footrpints)
footprints = st_transform(footprints, 5070)

# # Manually exclude 000934 and 000935 with same logic as above
# footprints = footprints |>
#   filter(!(mission_id %in% c("000934", "000935")))

# Temporary write to inspect
st_write(footprints, "/ofo-share/catalog-data-prep/stratification-data/all_drone_footprints.gpkg", delete_dsn = TRUE)

## For each record, get the max footprint that incorporates the plots and the drone missions and buffers by 50 m

footprint_list = list()

for (i in 1:nrow(pairings_unique)) {

  record = pairings_unique[i,]
  footprint_hn = footprints |> filter(mission_id == record$mission_id_hn)
  footprint_lo = footprints |> filter(mission_id == record$mission_id_lo)
  footprint_combined = st_union(footprint_hn, footprint_lo)

  footprint_list[[i]] = st_geometry(footprint_combined)

}

footprint_sfc = do.call(c, footprint_list)
if (length(footprint_sfc) != nrow(pairings_unique)) {
  stop("Number of footprints does not match number of pairings")
}

st_geometry(pairings_unique) = footprint_sfc

pairings_unique$hn_lo_pairing = TRUE

pairings_unique = pairings_unique |> st_cast("MULTIPOLYGON")

# Append to this the remaining footprints that are not in the pairings, so we can cluster
# them into groups based on spatial overlap
other_footprints = footprints |>
  filter(!(mission_id %in% c(pairings_unique$mission_id_hn, pairings_unique$mission_id_lo)))
other_footprints$hn_lo_pairing = FALSE

other_footprints$mission_id_hn = other_footprints$mission_id
other_footprints$mission_id_lo = other_footprints$mission_id

other_footprints = other_footprints |>
  select(mission_id_hn, mission_id_lo, hn_lo_pairing)

st_geometry(other_footprints) = "geometry"

all_footprints = rbind(pairings_unique, other_footprints)

## Define unique groups from the records here such that if there is any spatial overlap between two records, they are in the same group. Applies to chains of overlap.

# Create a spatial intersection matrix, using buffered plots to ensure separate groups are truly isolated
intersects_matrix = st_intersects(all_footprints |> st_buffer(10), sparse = FALSE)

# Function to find connected components using depth-first search
find_connected_components = function(adj_matrix) {
  n = nrow(adj_matrix)
  visited = rep(FALSE, n)
  groups = rep(NA, n)
  current_group = 1
  
  # Depth-first search helper function
  dfs = function(node) {
    visited[node] <<- TRUE
    groups[node] <<- current_group
    
    # Visit all unvisited neighbors
    neighbors = which(adj_matrix[node, ])
    for (neighbor in neighbors) {
      if (!visited[neighbor]) {
        dfs(neighbor)
      }
    }
  }
  
  # Process each node
  for (i in 1:n) {
    if (!visited[i]) {
      dfs(i)
      current_group = current_group + 1
    }
  }
  
  return(groups)
}

# Find connected components (groups)
all_footprints$group_id = find_connected_components(intersects_matrix)

# # Determine if any polygon in the group is an HN-LO pairing
# group_hn_lo = all_footprints |>
#   st_drop_geometry() |>
#   group_by(group_id) |>
#   summarize(has_hn_lo_pairing = any(hn_lo_pairing))

# # Join back to all_footprints
# # ... incomplete: only need to complete if we want plots in non-paired drone polys that touch paired polys to be assiged the hn-lo pairing

# Write to temp to inspect
st_write(pairings_unique, "/ofo-share/catalog-data-prep/stratification-data/pairings_grouped.gpkg", delete_dsn = TRUE)















## Get plot attributes assigned to support stratification

# For trees without a DBH, estimate it from the height using allometric equation
allometric_dbh_func = function(height) {
  exp((log(height - 1.3) + 0.3136489123372108) / 0.84623571)
}

ground_trees = ground_trees |>
  mutate(dbh_est = ifelse(is.na(dbh), allometric_dbh_func(height), dbh))

# Keep only the trees that have a DBH and it is > 20 cm, except for SSMT plots include all trees
# because they are mostly very small planted trees in a grid
ground_trees_filtered = ground_trees |>
  filter((!is.na(dbh_est) & dbh_est > 20) | as.numeric(plot_id) %in% c(277:285))

# Summarize the tree data

ground_trees_summ = ground_trees_filtered |>
  group_by(plot_id) |>
  summarize(
    n_trees_live = sum(live_dead == "L"),
    pct_dead = sum(live_dead == "D") / n() * 100,
    ba_live_sqm = sum(ifelse(live_dead == "L", (dbh_est / 200)^2 * pi, 0)),
    mean_ba_live = ba_live_sqm / n_trees_live,
  ) |>
  ungroup()

## Compute species comp groups

# Compute the top species and filter to those
table(ground_trees_filtered$species_code) |> sort(decreasing = TRUE)

# Make PIPJ be PIPO and AB be ABCO
ground_trees_filtered = ground_trees_filtered |>
  mutate(species_code = ifelse(species_code == "PIPJ", "PIPO", species_code),
         species_code = ifelse(species_code == "AB", "ABCO", species_code))

# Keep all species that have counts > 20, unless their names start with "UNK"
top_species = ground_trees_filtered |>
  filter(!str_starts(species_code, "UNK")) |>
  group_by(species_code) |>
  summarize(n = n()) |>
  filter(n > 20) |>
  pull(species_code)

trees_for_comp = ground_trees_filtered |>
  filter(species_code %in% top_species, live_dead == "L")


# Assuming your data frame is called 'trees' with columns: plot_id, species
# Create a species composition matrix (plots x species)
species_matrix = trees_for_comp |>
  count(plot_id, species_code) |>
  pivot_wider(names_from = species_code, 
              values_from = n, 
              values_fill = 0)

# Extract plot IDs and convert to matrix
plot_ids = species_matrix$plot_id
species_matrix = species_matrix |>
  select(-plot_id) |>
  as.matrix()

# Optional: Convert to relative abundance (proportions)
species_matrix = species_matrix / rowSums(species_matrix)

# Perform k-means clustering (k=8 in hopes of classifying OHDS as its own group)
set.seed(123)  # For reproducibility
k = 8
clusters = kmeans(species_matrix, centers = k, nstart = 25)

# Create result data frame
plot_clusters = data.frame(
  plot_id = plot_ids,
  sp_comp_group = clusters$cluster
)

# Join back to original data
ground_trees_summ = ground_trees_summ |>
  left_join(plot_clusters, by = "plot_id")

# # If cluster is NA, set to 99
# ground_trees_summ = ground_trees_summ |>
#   mutate(sp_comp_group = ifelse(is.na(sp_comp_group), 99, sp_comp_group))

# Join summarized tree data to plot data
ground_plots_summ = ground_plots |>
  left_join(ground_trees_summ, by = "plot_id")

# Compute tree density
ground_plots_summ = ground_plots_summ |>
  mutate(trees_per_ha = n_trees_live / (plot_area_ha), 
         ba_sqm_per_ha = ba_live_sqm / (plot_area_ha))

# Only considering plots with an area over 0.03 ha (to exclude the RMR mini-plots)
ground_plots_summ = ground_plots_summ |>
  filter(plot_area_ha >= 0.03)

# Also remove the very outlier SBGEO plot and rancho marino plots which are all under one drone polygon
ground_plots_summ = ground_plots_summ |>
  filter(!(plot_id %in% c("0081", str_pad(242:269, width = 4, pad = "0", side = "left"))))


# Write to temp to inspect
st_write(ground_plots_summ, "/ofo-share/catalog-data-prep/stratification-data/ground_plots_summ.gpkg", delete_dsn = TRUE)

# Also save as centroids so it's easier to see in GIS when zoomed out and also to extract ppt and
# ecoregion values
ground_plots_summ_centroids = ground_plots_summ |>
  st_centroid()
st_write(ground_plots_summ_centroids, "/ofo-share/catalog-data-prep/stratification-data/ground_plots_summ_centroids.gpkg", delete_dsn = TRUE)


## Pull in PRISM ppt and EPA ecoregion

ppt_extracted = extract(ppt, ground_plots_summ_centroids)
ecoregion_extracted = st_intersection(ecoregion, st_transform(ground_plots_summ_centroids, st_crs(ecoregion)))$US_L3NAME

ground_plots_summ$ecoregion = ecoregion_extracted
ground_plots_summ$ppt = ppt_extracted[,2]



# Determine for each plot whether it is: not under a drone polygon, under a non-paired drone
# polygon, or under a paired drone polygon. Keep in mind they could be under multiple polygons and
# in such a case the paired drone polygon takes precedence. Loop through each plot. Also extract the
# drone group_id.
ground_plots_summ$drone_pairing_tier = "no_drone_footprint"
ground_plots_summ$drone_group_id = NA

for (i in 1:nrow(ground_plots_summ)) {
  plot = ground_plots_summ[i,]
  # Find all footprints that intersect with this plot
  intersecting_footprints = st_intersection(all_footprints, st_transform(plot, st_crs(all_footprints)))
  if (nrow(intersecting_footprints) > 0) {

    # Assign the group_id (if multiple, assign the first one)
    ground_plots_summ$drone_group_id[i] = intersecting_footprints$group_id[1]

    # If any of the intersecting footprints are paired, assign "paired_drone_footprint"
    if (any(intersecting_footprints$hn_lo_pairing)) {
      ground_plots_summ$drone_pairing_tier[i] = "paired_drone_footprint"
    } else {
      ground_plots_summ$drone_pairing_tier[i] = "non_paired_drone_footprint"
    }
  }
}

# Top precedence: was the plot successfully aligned to a HN-LO drone pairing? Get the aligned plot
# IDs from Amritha's model input file.
ground_plots_summ[ground_plots_summ$plot_id %in% aligned_pairings$plot_id, "drone_pairing_tier"] = "aligned_paired_drone_footprint"


# For the plots that do not overlap drone footprints, assign a new group_id that is unique to each
# plot and doesn't overlap any other plot
max_existing_group_id = max(all_footprints$group_id)
for (i in 1:nrow(ground_plots_summ)) {
  plot = ground_plots_summ[i,]
  if (plot$drone_pairing_tier == "no_drone_footprint") {
    max_existing_group_id = max_existing_group_id + 1
    ground_plots_summ$drone_group_id[i] = max_existing_group_id
  }
}


# Find the drone missions that are in groups that do not contain any plots, and convert them to
# "plots" for purposes of stratification, so we can stratify the drone-only footprints as well

drone_groups_with_plots = unique(ground_plots_summ$drone_group_id)
drone_footprints_without_plots = all_footprints |>
  filter(!(group_id %in% drone_groups_with_plots))

drone_no_plot = drone_footprints_without_plots |>
  st_centroid() |>
  rename(drone_group_id = group_id) |>
  mutate(plot_id = mission_id_hn) # Using mission_id_hn as a stand-in for plot_id here

# Extract ppt, ecoregion, project and drone_pairing_tier = "drone_no_plots"

ppt_extracted_drone = extract(ppt, drone_no_plot)
ecoregion_extracted_drone = st_intersection(ecoregion, st_transform(drone_no_plot, st_crs(ecoregion)))$US_L3NAME
drone_no_plot$ppt = ppt_extracted_drone[,2]
drone_no_plot$ecoregion = ecoregion_extracted_drone
drone_no_plot$drone_pairing_tier = "drone_no_plots"

project_lookup = footprints |>
  select(mission_id, project_name) |>
  st_drop_geometry() |>
  distinct()

drone_no_plot = drone_no_plot |>
  left_join(project_lookup, by = c("mission_id_hn" = "mission_id"))



st_geometry(drone_no_plot) = "geom"

# Temporary write to inspect
st_write(drone_no_plot, "/ofo-share/catalog-data-prep/stratification-data/drone_footprints_without_plots.gpkg", delete_dsn = TRUE)


# Merge together with the ground plots dataframe for stratification
for_strat = bind_rows(ground_plots_summ |> st_transform(st_crs(drone_no_plot)), drone_no_plot)


# Keep only the columns needed for stratification
for_strat = for_strat |>
  select(
    plot_id, drone_group_id, drone_pairing_tier, project_name, survey_date, plot_area_ha, n_trees_live,
    ecoregion, ppt, trees_per_ha, ba_sqm_per_ha, mean_ba_live, sp_comp_group
  ) |>
  mutate(survey_year = str_sub(survey_date, 1, 4))

# Temporary write to inspect
st_write(for_strat, "/ofo-share/catalog-data-prep/stratification-data/for_strat.gpkg", delete_dsn = TRUE)

# Remove geometry for stratification processing and save
for_strat_df = for_strat |>
  st_drop_geometry()

write_csv(for_strat_df, "/ofo-share/catalog-data-prep/stratification-data/for_strat.csv")
