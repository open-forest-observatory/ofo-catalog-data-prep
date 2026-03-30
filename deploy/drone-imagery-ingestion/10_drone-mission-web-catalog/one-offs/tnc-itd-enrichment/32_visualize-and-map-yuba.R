# Purpose: Map the drone mission footprints and ground plots for the Yuba focal area

library(sf)
library(tidyverse)
library(ceramic)
library(tidyterra)
library(rnaturalearth)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"

YUBA_BOUNDARY_FILEPATH   = file.path(DELIVERABLES_DIR, "raw-input/north_yuba_area.kml")
FOOTPRINTS_FILEPATH      = file.path(DELIVERABLES_DIR, "composite-missions/footprint-boundaries-delivery/composite-drone-plot-summaries_yuba.gpkg")
NADIR_MISSIONS_FILEPATH  = file.path(DELIVERABLES_DIR, "composite-missions/footprint-boundaries-delivery/individual-drone-plot-summaries_yuba.gpkg")
GROUND_PLOTS_FILEPATH    = file.path(DELIVERABLES_DIR, "ground-reference-delivery/ofo_ground-reference_plots_yuba.gpkg")

FIGURES_DIR = file.path(DELIVERABLES_DIR, "figures")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

MAP_FILEPATH = file.path(FIGURES_DIR, "yuba-drone-and-ground-plots-map.png")

MAP_BUFFER_M = 8000  # padding around data extent for map display


# --- Load data ---

yuba_boundary = st_read(YUBA_BOUNDARY_FILEPATH, quiet = TRUE) |>
  st_union() |>
  st_as_sf() |>
  st_transform(3857)

footprints = st_read(FOOTPRINTS_FILEPATH, quiet = TRUE) |>
  st_transform(3857)

nadir_missions = st_read(NADIR_MISSIONS_FILEPATH, quiet = TRUE) |>
  st_transform(3857)

ground_plots = st_read(GROUND_PLOTS_FILEPATH, quiet = TRUE) |>
  st_transform(3857)

# Represent all layers as centroids for the regional map
nadir_pts       = st_centroid(nadir_missions)
footprint_pts   = st_centroid(footprints)
ground_plot_pts = st_centroid(ground_plots)


# --- Compute map extent ---

# Base extent on the Yuba boundary only, buffered for padding
map_extent = st_bbox(yuba_boundary) |>
  st_as_sfc() |>
  st_as_sf() |>
  st_buffer(MAP_BUFFER_M) |>
  st_bbox() |>
  st_as_sfc() |>
  st_as_sf()


# --- Fetch basemap and state boundaries ---

basemap = ceramic::cc_location(loc = map_extent)

states = ne_states(country = "United States of America", returnclass = "sf") |>
  st_transform(3857) |>
  st_intersection(map_extent)

# counties = ne_download(scale = "large", type = "counties", category = "cultural", returnclass = "sf") |>
#   st_transform(3857) |>
#   st_intersection(map_extent)


# --- Map ---

p = ggplot() +
  geom_spatraster_rgb(data = basemap, alpha = 0.65) +
  geom_sf(data = states,   color = "grey20", fill = NA, linewidth = 0.7) +
  geom_sf(data = yuba_boundary,
          color = "black", fill = NA, linewidth = 1.1) +
  geom_sf(data = nadir_pts,
          aes(fill = "Nadir-only missions"),
          color = "black", size = 4.5, shape = 21) +
  geom_sf(data = footprint_pts,
          aes(fill = "Composite missions"),
          color = "black", size = 4.5, shape = 21) +
  geom_sf(data = ground_plot_pts,
          aes(fill = "Ground reference plots"),
          color = "black", size = 2.5, shape = 21) +
  scale_fill_manual(
    name = NULL,
    breaks = c("Nadir-only missions", "Composite missions", "Ground reference plots"),
    values = c(
      "Nadir-only missions"    = "#1f77b4",
      "Composite missions"     = "#ff9100",
      "Ground reference plots" = "#2ca02c"
    )
  ) +
  guides(fill = guide_legend(
    override.aes = list(
      shape = c(21, 21, 21),
      size  = c(4.5, 4.5, 2.5)
    )
  )) +
  coord_sf(crs = 4326) +
  theme_bw(14) +
  theme(
    panel.grid    = element_blank(),
    legend.position  = "bottom",
    legend.text   = element_text(size = 12),
    plot.title    = element_text(size = 14, face = "bold"),
    plot.caption  = element_text(size = 10, color = "grey40")
  )

png(filename = MAP_FILEPATH, width = 10, height = 8, units = "in", res = 300)
print(p)
dev.off()

cat("Saved map to", MAP_FILEPATH, "\n")


# Compute total flight area summed across the nadir missions
total_nadir_area_ha = sum(st_area(nadir_missions)) / 10000

# Compute average flight area
average_nadir_area_ha = mean(st_area(nadir_missions)) / 10000
cat("Total area covered by nadir missions:", round(total_nadir_area_ha, 2), "ha\n")
cat("Average area covered per nadir mission:", round(average_nadir_area_ha, 2), "ha\n")
