# Purpose: Read the ground reference data from Google Sheets and plot boundaries from local files,
# and write the data to gpkg files for publishing. Initially written and being used to create field
# refence data files to share with OFO collaborators.

library(dplyr)
library(googlesheets4)
library(stringr)
library(sf)

source(file.path("src", "web-catalog-creation_ground-ref-data.R"))

PLOT_BOUNDARIES_PATH = "~/repo-data-local/ofo-catalog-data-prep/field-plot-boundaries"
GOOGLE_SHEET_ID = "1GjDseDCR1BX_EIkJAni7rk2zvK6nHmZz1nOFBd1d6k4"

# ---- Processing

# Load and prep field ref data

tabular_data = read_and_standardize_tabular_field_ref_data(GOOGLE_SHEET_ID)
bounds = read_and_merge_plot_boundaries(plot_boundaries_dir = PLOT_BOUNDARIES_PATH, base_ofo_url = "", plot_details_dir = "")

# Remove some fields that were added to the bounds table for web catalog creation that are not
# needed for the publishable data files
bounds = bounds |>
  select(-c(bounds_plot_id_link))

check_field_ref_data(tabular_data, bounds)

trees = prep_trees(trees = tabular_data$trees, species_codes = tabular_data$species_codes)
# Any warning of duplicate coedes probably due to JUOC/JUGR 64

plots = prep_plots(tabular_data$plots)

# # For some reason, subplot_shape, project_id, contributor_plot_id, plot_area, survey_date_approx is being read as a list column (just one entry per plot though), so
# # unlist it, first setting NULL to NA. Seems to be a bug in googlesheets4 that we should find a
# # better solution to.
# plots = plots |>
#   mutate(subplot_shape = lapply(subplot_shape, function(x) ifelse(is.null(x), NA, x)))
# plots = plots |>
#   mutate(subplot_shape = unlist(subplot_shape))
# plots = plots |>
#   mutate(project_id = lapply(project_id, function(x) ifelse(is.null(x), NA, x)))
# plots = plots |>
#   mutate(project_id = unlist(project_id))
# plots = plots |>
#   mutate(contributor_plot_id = lapply(contributor_plot_id, function(x) ifelse(is.null(x), NA, x)))
# plots = plots |>
#   mutate(contributor_plot_id = unlist(contributor_plot_id))
# plots = plots |>
#   mutate(plot_area = lapply(plot_area, function(x) ifelse(is.null(x), NA, x)))
# plots = plots |>
#   mutate(plot_area = unlist(plot_area) |> as.numeric())
# plots = plots |>
#   mutate(survey_date_approx = lapply(survey_date_approx, function(x) ifelse(is.null(x), NA, x)))
# plots = plots |>
#   mutate(survey_date_approx = unlist(survey_date_approx))


  
# For some reason, contributor_tree_id is being read as a list column (just one entry per tree
# though), so unlist it, first setting NULL to NA
trees = trees |>
  mutate(contributor_tree_id = lapply(contributor_tree_id, function(x) ifelse(is.null(x), NA, x)))
trees = trees |>
  mutate(contributor_tree_id = unlist(contributor_tree_id))
  
  
# Remove plots that have a "display message" which is used to indicate something is anomalous about them
plots = plots |>
  filter(is.na(display_message))

trees_for_plot_summary = prep_trees_for_plot_summary(trees)

plot_level_tree_summary = summarize_trees_by_plot(trees_for_plot_summary)

plot_summary = compile_plot_summary_table(plots = plots,
                                          projects = tabular_data$projects,
                                          plot_level_tree_summary = plot_level_tree_summary,
                                          bounds = bounds,
                                          base_ofo_url = "",
                                          plot_details_dir = "")

# Remove some fields that were added to the plot summary table for web catalog creation that are not
# needed for the publishable data files
plot_summary = plot_summary |>
  select(-c(plot_id_link))
  
# Remove embargoed plots
plot_summary = plot_summary |>
  filter(!embargoed)

# Imprecise test for whether the data is still in the process of entry and should be skipped. TODO:
# Consider turning into a function
plot_summary = plot_summary |>
  dplyr::filter(!is.na(plot_id) & !is.na(survey_date))

# Merge the plot summary to the bounds data to get the geometry
bounds = bounds |> select(-area_ha_sf)
plot_bounds_w_summary = right_join(bounds, plot_summary, by = c("plot_id" = "plot_id"))

# Remove summary data without any gemoetry (no plot bounds data -- presumably partially entered
# plots still in progress)
plot_bounds_w_summary = plot_bounds_w_summary[!st_is_empty(plot_bounds_w_summary), ]

# Keep the columns that are relevant
plot_bounds_w_summary = plot_bounds_w_summary |>
  select(
    plot_id,
    project_id,
    project_name = name_short,
    investigator_names,
    investigator_contacts,
    license_short,
    license,
    hyperplot_id,
    survey_date_approx,
    survey_date,
    plot_shape,
    plot_area_ha,
    subplots,
    subplot_shape,
    subplot_area,
    top_species,
    height_measured,
    includes_snags,
    includes_damage,
    damage_codes_inspected,
    min_dbh,
    min_dbh_live,
    max_dbh_of_primary_trees,
    min_ht,
    min_dbh_ohvis,
    min_ht_ohvis,
    plot_lon,
    plot_lat,
    num_ohvis_trees_excluded,
    contributor_plot_id
  )




# Write the plot-level data to a gpkg
st_write(plot_bounds_w_summary, "~/temp/ofo_ground-reference_plots.gpkg", delete_dsn = TRUE)


### Trees

# Which trees to we not have plot data for (either because no polygon, or no plot-level tabular data)?
trees_no_plot_data = trees |>
  filter(!plot_id %in% plot_bounds_w_summary$plot_id)
table(trees_no_plot_data$plot_id)
# It is just embargoed data and ones with no polygons (in process of being entered)
# 53-56: OSU
# 57: FOCAL incomplete ABCO plot
# 82, 83, 87: embargoed Lamping plots
# 9999: dummy plot

# Filter to the trees that are in the field plots that we have bounds & summary data for
trees = trees |>
  filter(plot_id %in% plot_bounds_w_summary$plot_id) |>
  filter(plot_id != 9999)

# Keep just the columns that are relevant
trees = trees |>
  select(
    plot_id,
    subplot_id,
    tree_lat,
    tree_lon,
    height,
    height_allometric,
    height_above_plot_center,
    dbh,
    species_code = sp_code,
    growth_form,
    live_dead,
    crown_position,
    ohvis,
    crown_ratio,
    crown_ratio_compacted,
    height_to_crown,
    height_to_needle,
    scorch_height,
    percent_prefire_crown_green,
    precent_postfire_crown_green,
    crown_width_1,
    crown_width_2,
    crown_width_allometric,
    decay_class,
    damage_1,
    damage_2,
    damage_3,
    damage_4,
    damage_5,
    contributor_tree_id,
    notes
  )


# Make spatial
trees_sf = sf::st_as_sf(trees, coords = c("tree_lon", "tree_lat"), crs = 4326)

# Write to gpkg
st_write(trees_sf, "~/temp/ofo_ground-reference_trees.gpkg", delete_dsn = TRUE)
