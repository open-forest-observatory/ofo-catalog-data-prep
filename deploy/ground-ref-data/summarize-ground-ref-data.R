# Purpose: Summarize the ground reference data in various ways to inform modeling decisions. Not part of a pipeline.

library(sf)
library(tidyverse)

REPO_DATA_LOCAL_PATH = "~/repo-data-local/ofo-catalog-data-prep"
GROUND_REF_PLOTS_OUTPUT_GPKG = file.path(REPO_DATA_LOCAL_PATH, "ofo_ground-reference_plots.gpkg")
GROUND_REF_TREES_OUT_GPKG = file.path(REPO_DATA_LOCAL_PATH, "ofo_ground-reference_trees.gpkg")

GROUND_REF_SUMMARY_OUT_CSV = file.path(REPO_DATA_LOCAL_PATH, "ofo_ground-reference_plot_hardwood-prop-summary.csv")

plots = st_read(GROUND_REF_PLOTS_OUTPUT_GPKG)
trees = st_read(GROUND_REF_TREES_OUT_GPKG)

## Classify each tree as hardwood true/false
table(trees$species_code) |> sort(decreasing = TRUE)

hardwood_codes = c("NODE3", "QUWI2", "QUKE", "QUCH2", "ALRU2", "ARME", "PREM", "UMCA", "FRACAL", "SASC", "POTR5", "QUAG", "QUDO", "SALIX", "ARCVISM", "ACMA3", "ARVI4", "RHPU", "ACCI", "CONU4", "CELE", "CEAPAL", "QUDE", "QUEXMO", "RHOOCC", "SESE3", "RHCR", "ARCGLA", "POBAT", "AECA", "LUNSUB", "FRACAC5", "QUEV", "ACRU", "AMAL2", "ARCPUN", "CEALEU", "FRPU7", "QUGA4", "ACNE2", "BEOC", "CERMONG", "FRLA", "RHIALI")

trees = trees |>
  mutate(is_hardwood = ifelse(species_code %in% hardwood_codes, TRUE, FALSE))
table(trees$is_hardwood)

# Ignoring tree species UNKSNAG, compute the percent of trees per plot that are hardwoods, by basal
# area, total height, and tree count

trees_no_snags = trees |>
  filter(species_code != "UNKSNAG")
hardwood_summary = trees_no_snags |>
  mutate(ba = (dbh / 2)^2 * pi) |>
  group_by(plot_id) |>
  summarize(
    total_ba = sum(ba),
    hardwood_ba = sum(ba[is_hardwood]),
    prop_hardwood_ba = hardwood_ba / total_ba,
    total_height = sum(height),
    hardwood_height = sum(height[is_hardwood]),
    prop_hardwood_height = hardwood_height / total_height,
    total_count = n(),
    hardwood_count = sum(is_hardwood),
    prop_hardwood_count = hardwood_count / total_count
  )

# Compute an overall prop_hardwood for each plot that prioritizes basal area, then height, then count
hardwood_summary = hardwood_summary |>
  mutate(
    overall_prop_hardwood = case_when(
      !is.na(prop_hardwood_ba) ~ prop_hardwood_ba,
      is.na(prop_hardwood_ba) & !is.na(prop_hardwood_height) ~ prop_hardwood_height,
      is.na(prop_hardwood_ba) & is.na(prop_hardwood_height) & !is.na(prop_hardwood_count) ~ prop_hardwood_count,
      TRUE ~ NA
    )
  )


# Make sure we have a summary for every plot
plots_with_summ = unique(hardwood_summary$plot_id)
plots_all = unique(plots$plot_id)

setdiff(plots_all, plots_with_summ)
setdiff(plots_with_summ, plots_all)
# Just missing tree data for one NEON plot (presumably no trees in it)

hw_summ_simp = hardwood_summary |>
  select(plot_id, prop_hardwood = overall_prop_hardwood) |>
  st_drop_geometry()

plots_w_hw_summ = plots |>
  left_join(hw_summ_simp, by = "plot_id") |>
  st_drop_geometry() |>
  select(plot_id, prop_hardwood) |> # removed: project_id, project_name,
  arrange(plot_id)

# Write out summary table
write_csv(plots_w_hw_summ, GROUND_REF_SUMMARY_OUT_CSV)
