# Purpose: For all CHMs within a specified mission ID range, detect treetops and save results to a
# file.

library(tidyverse)
library(sf)
library(lidR)
library(terra)
library(nngeo)
library(smoothr)
library(furrr)

## Set constants
source("deploy/drone-imagery-ingestion/00_set-constants.R")

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

detect_ttops_and_crowns = function(chm_file_foc) {

  chm_identifier = tools::file_path_sans_ext(basename(chm_file_foc))
  # Get the mission ID from the filename (the first 6 characters)
  mission_id = str_sub(basename(chm_file_foc), 1, 6)
  # Get the photogrammetry processing run ID from the filename
  photogrammetry_run_id = str_match(chm_file_foc, "/processed_(\\d{2})")[,2]
  # Construct the output directory path within the object store bucket: where to upload the results
  output_directory = paste0(REMOTE_MISSIONS_DIR, "/", mission_id, "/processed_", photogrammetry_run_id, "/", ITD_FOLDER)

  # Download the CHM from JS2 Object Store and load
  remote_file = paste0(RCLONE_REMOTE, ":", REMOTE_MISSIONS_DIR, chm_file_foc)
  temp_chm_file = tempfile(chm_identifier, fileext = ".tif")
  command = paste("rclone copyto", remote_file, temp_chm_file, "--progress --transfers 32 --checkers 32 --stats 1s --retries 5 --retries-sleep=15s --s3-upload-cutoff 100Mi --s3-chunk-size 100Mi --s3-upload-concurrency 16 --multi-thread-streams 2", sep = " ")
  system(command)
  chm = terra::rast(temp_chm_file)


  # Prep the CHM
  chm_smooth = resample_and_smooth_chm(chm, CHM_RES, CHM_SMOOTH_WIDTH)

  # Detect trees
  ttops = predict_trees_from_chm(chm_smooth, LMF_A, LMF_B, LMF_C, LMF_DIAM_MIN, LMF_DIAM_MAX)

  # Extract tree height from the non-smoothed CHM
  ttops$Z = extract(chm, ttops)[,2]

  # For removing edge trees:
  # Get the bounds of the CHM
  chm_non_na = chm
  chm_non_na[!is.na(chm)] = 1
  chm_poly = as.polygons(chm_non_na, values = FALSE)

  # Buffer in by 10 m
  chm_poly_buff = buffer(chm_poly, width = -10)

  # Delineate crowns: silva
  crowns_silva = lidR::silva2016(chm_smooth, ttops, max_cr_factor = 0.24, exclusion = 0.1)()
  crowns_silva <- as.polygons(rast(crowns_silva))
  crowns_silva <- st_as_sf(crowns_silva)
  crowns_silva <- st_cast(crowns_silva, "MULTIPOLYGON")
  crowns_silva <- st_cast(crowns_silva, "POLYGON")
  crowns_silva <- st_remove_holes(crowns_silva)
  crowns_silva <- st_make_valid(crowns_silva)
  crowns_silva <- smooth(crowns_silva, method = "ksmooth", smoothness = 3)
  crowns_silva <- st_simplify(crowns_silva, preserveTopology = TRUE, dTolerance = 0.1)

  crowns_watershed = lidR::watershed(chm_smooth, th_tree = 2, tol = 0, ext = 1)()
  crowns_watershed <- as.polygons(rast(crowns_watershed))
  crowns_watershed <- st_as_sf(crowns_watershed)
  crowns_watershed <- st_cast(crowns_watershed, "MULTIPOLYGON")
  crowns_watershed <- st_cast(crowns_watershed, "POLYGON")
  crowns_watershed <- st_remove_holes(crowns_watershed)
  crowns_watershed <- st_make_valid(crowns_watershed)
  crowns_watershed <- smooth(crowns_watershed, method = "ksmooth", smoothness = 3)
  crowns_watershed <- st_simplify(crowns_watershed, preserveTopology = TRUE, dTolerance = 0.1)

  # Crop the ttops
  ttops = st_intersection(ttops, chm_poly_buff |> st_as_sf())

  # Assign crowns the treetop height and remove those that have no treetops in them
  crowns_silva = st_join(crowns_silva, ttops)
  crowns_silva = crowns_silva[, -1]
  crowns_silva = crowns_silva[!is.na(crowns_silva$Z),]

  crowns_watershed = st_join(crowns_watershed, ttops)
  crowns_watershed = crowns_watershed[, -1]
  crowns_watershed = crowns_watershed[!is.na(crowns_watershed$Z),]

  # Write predicted treetops and crowns to a temp folder for uploading
  random_number = sample(1:1000000, 1) |> str_pad(width = 6, pad = "0", side = "left")
  temp_folder = file.path(TEMPDIR, paste0("itdtemp_", random_number))
  dir.create(temp_folder, recursive = TRUE, showWarnings = FALSE)

  ttops_tempfile = file.path(temp_folder, paste0(mission_id, "_treetops.gpkg"))
  crowns_watershed_tempfile = file.path(temp_folder, paste0(mission_id, "_crowns-watershed.gpkg"))
  crowns_silva_tempfile = file.path(temp_folder, paste0(mission_id, "_crowns-silva.gpkg"))

  st_write(ttops, ttops_tempfile, delete_dsn = TRUE)
  st_write(crowns_watershed, crowns_watershed_tempfile, delete_dsn = TRUE)
  st_write(crowns_silva, crowns_silva_tempfile, delete_dsn = TRUE)

  # Upload the results to the Object Store
  # Construct data transfer command line call
  remote_dir = paste0(RCLONE_REMOTE, ":", output_directory)
  command = paste("rclone copy", temp_folder, remote_dir, "--progress --transfers 32 --checkers 32 --stats 1s --retries 5 --retries-sleep=15s --s3-upload-cutoff 100Mi --s3-chunk-size 100Mi --s3-upload-concurrency 16", sep = " ")
  result = system(command)


  gc()
  file.remove(ttops_tempfile, crowns_watershed_tempfile, crowns_silva_tempfile, temp_chm_file)
  unlink(temp_folder, recursive = TRUE)
}


chm_paths = read_csv(CHMS_FOR_ITD_LIST_PATH)$chm_path

# Keep only the ones made from the mesh
chm_paths = chm_paths[grepl("_chm-mesh.tif$", chm_paths)]

plan(multicore(workers = availableCores() / 2))

# Process each CHM
future_walk(chm_paths, detect_ttops_and_crowns, .progress = TRUE, .options = furrr_options(scheduling = Inf))
