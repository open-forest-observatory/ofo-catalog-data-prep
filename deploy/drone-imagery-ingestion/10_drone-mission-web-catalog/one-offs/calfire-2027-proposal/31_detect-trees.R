# Purpose: For each focal-area polygon (one per CHM mission), crop the corresponding CHM to that
# focal area, detect treetops, and save the results to a file.

library(tidyverse)
library(sf)
library(lidR)
library(terra)
library(nngeo)
library(smoothr)
library(furrr)

## Set constants
# source("deploy/drone-imagery-ingestion/00_set-constants.R")

CHM_FOLDER = "/ofo-share/scratch/derek/calfire-2027-proposal-prep/chms"
ITD_FOLDER = "/ofo-share/scratch/derek/calfire-2027-proposal-prep/itd"
FOCAL_AREAS_FILE = "/ofo-share/scratch/derek/calfire-2027-proposal-prep/focal-areas/focal-areas-calfireprop.gpkg"

# Distance (m) to buffer the focal-area polygon outward before cropping the CHM, so tree detection
# has real CHM context beyond the focal boundary (avoiding edge artifacts). Detected treetops are
# clipped back to the unbuffered focal polygon.
FOCAL_BUFFER = 10

## Functions

# Function to create a variable radius window function for LMF
make_win_fun <- function(a, b, c, diam_min, diam_max) {
  win_fun <- function(x) {
    win <- a + b*x + c*x^2
    win[win < diam_min] = diam_min
    win[win > diam_max] = diam_max
    return(win)
  }
  return(win_fun)
}

# Function to resample and smooth a CHM
resample_and_smooth_chm = function(chm, res, smooth_width) {
  chm_resamp <- terra::project(chm, terra::crs(chm), res = res, method = "bilinear")
  chm_smooth <- terra::focal(chm_resamp, w = matrix(1, smooth_width, smooth_width), mean, na.rm = TRUE)
  return(chm_smooth)
}

# Function to predict trees from a prepped (resampled and smoothed) CHM
predict_trees_from_chm <- function(chm, lmf_a, lmf_b, lmf_c, lmf_diam_min, lmf_diam_max) {
  win_fun <- make_win_fun(lmf_a, lmf_b, lmf_c, lmf_diam_min, lmf_diam_max)
  ttops <- lidR::locate_trees(chm, algorithm = lmf(ws = win_fun, shape = "circular", hmin = 5), )
  return(ttops)
}



#### Workflow

ITD_PARAMETERIZATION_ID = 1
CHM_RES = 0.25
CHM_SMOOTH_WIDTH = 7
LMF_A = 0
LMF_B = 0.11
LMF_C = 0
LMF_DIAM_MIN = 0.5
LMF_DIAM_MAX = 100

# Detect treetops within a single focal-area polygon, using the CHM for that polygon's mission ID.
# `mission_id` is the integer mission ID; `focal_areas` is the full sf of focal-area polygons.
detect_ttops_and_crowns = function(mission_id, focal_areas) {

  # CHM filenames use a zero-padded 6-digit mission ID
  mission_id_padded = sprintf("%06d", mission_id)
  chm_file_foc = file.path(CHM_FOLDER, paste0(mission_id_padded, "_chm-mesh.tif"))

  if (!file.exists(chm_file_foc)) {
    warning("No CHM found for mission ", mission_id_padded, " at ", chm_file_foc, "; skipping.")
    return(invisible(NULL))
  }

  # The focal-area polygon for this mission
  focal = focal_areas[focal_areas$mission_id == mission_id, ]

  chm = terra::rast(chm_file_foc)

  # Reproject the focal polygon to the CHM's CRS, then crop and mask the CHM to the focal area
  # (buffered outward so detection has context beyond the focal boundary)
  focal = st_transform(focal, terra::crs(chm))
  focal_buff = st_buffer(focal, FOCAL_BUFFER)
  chm = terra::crop(chm, terra::vect(focal_buff), mask = TRUE)

  # Prep the CHM
  chm_smooth = resample_and_smooth_chm(chm, CHM_RES, CHM_SMOOTH_WIDTH)

  # Detect trees
  ttops = predict_trees_from_chm(chm_smooth, LMF_A, LMF_B, LMF_C, LMF_DIAM_MIN, LMF_DIAM_MAX)

  # Extract tree height from the non-smoothed CHM
  ttops$Z = extract(chm, ttops)[,2]

  # NOTE: Crown delineation is currently excluded due to long compute time, but the code is left here for reference and potential future use.
  # # Delineate crowns: silva
  crowns_silva = lidR::silva2016(chm_smooth, ttops, max_cr_factor = 0.24, exclusion = 0.2)()
  crowns_silva <- as.polygons(crowns_silva)
  crowns_silva <- st_as_sf(crowns_silva)
  # crowns_silva <- st_simplify(crowns_silva, preserveTopology = TRUE, dTolerance = 0.1)
  crowns_silva <- st_cast(crowns_silva, "MULTIPOLYGON")
  crowns_silva <- st_cast(crowns_silva, "POLYGON")
  crowns_silva <- st_remove_holes(crowns_silva)
  crowns_silva <- st_make_valid(crowns_silva)
  crowns_silva <- smooth(crowns_silva, method = "ksmooth", smoothness = 3)
  crowns_silva <- st_simplify(crowns_silva, preserveTopology = TRUE, dTolerance = 0.1)

  # Delineate crowns: watershed
  crowns_watershed = lidR::watershed(chm_smooth, th_tree = 2, tol = 0, ext = 1)()
  crowns_watershed <- as.polygons(crowns_watershed)
  crowns_watershed <- st_as_sf(crowns_watershed)
  crowns_watershed <- st_cast(crowns_watershed, "MULTIPOLYGON")
  crowns_watershed <- st_cast(crowns_watershed, "POLYGON")
  crowns_watershed <- st_remove_holes(crowns_watershed)
  crowns_watershed <- st_make_valid(crowns_watershed)
  crowns_watershed <- smooth(crowns_watershed, method = "ksmooth", smoothness = 3)
  crowns_watershed <- st_simplify(crowns_watershed, preserveTopology = TRUE, dTolerance = 0.1)

  # Clip the treetops and crowns to the (unbuffered) focal-area polygon
  ttops = st_intersection(ttops, focal)
  crowns_silva = st_intersection(crowns_silva, focal)
  crowns_watershed = st_intersection(crowns_watershed, focal)

  # Write predicted treetops and crowns to the ITD folder
  if (!dir.exists(ITD_FOLDER)) dir.create(ITD_FOLDER, recursive = TRUE)
  ttops_outfile = file.path(ITD_FOLDER, paste0(mission_id_padded, "_treetops.gpkg"))
  crowns_silva_outfile = file.path(ITD_FOLDER, paste0(mission_id_padded, "_crowns-silva.gpkg"))
  crowns_watershed_outfile = file.path(ITD_FOLDER, paste0(mission_id_padded, "_crowns-watershed.gpkg"))

  st_write(ttops, ttops_outfile, delete_dsn = TRUE)
  st_write(crowns_silva, crowns_silva_outfile, delete_dsn = TRUE)
  st_write(crowns_watershed, crowns_watershed_outfile, delete_dsn = TRUE)

  gc()
}


# Read the focal-area polygons (one per CHM mission, identified by the mission_id field)
focal_areas = st_read(FOCAL_AREAS_FILE, quiet = TRUE)

# Process each focal area / mission
future_walk(
  focal_areas$mission_id,
  detect_ttops_and_crowns,
  focal_areas = focal_areas,
  .progress = TRUE,
  .options = furrr_options(scheduling = Inf)
)
