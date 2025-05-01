# Purpose: Take the photogrammetry products and compute the "deliverable" versions of the outputs
# (e.g. CHM, cloud-optimized, thumbnails) by postprocessing. Perform at the mission level.

# This script requires the existince of a conda environment named `untwine` with the `untwine` package
# installed. This can be created with the system command: `conda create -n untwine -c conda-forge untwine
# -y`

library(tidyverse)
library(sf)
library(lidR)
library(terra)
library(purrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/utils.R")

## Workflow

terra::terraOptions(memfrac = TERRA_MEMFRAC)



# Function to crop raster to mission polygon and write as COG
crop_raster_save_cog = function(raster_filepath_foc, mission_polygon, output_path) {

  # Read, crop, and write the raster as COG
  raster = terra::rast(raster_filepath_foc)
  mission_polygon_matchedcrs = st_transform(mission_polygon, st_crs(raster))
  raster_cropped = terra::crop(raster, mission_polygon_matchedcrs, mask = TRUE, extend = TRUE)
  output_file_path = file.path(
    output_path, "full",
    paste0(basename(raster_filepath_foc))
  )

  terra::writeRaster(
    raster_cropped,
    output_file_path,
    overwrite = TRUE,
    filetype = "COG",
    gdal = "BIGTIFF=IF_SAFER",
    todisk = TRUE
  )
}


# Function to make a CHM from two raster file paths
make_chm = function(dsm_filepath_foc, dtm_filepath_foc) {
  # Read the rasters
  dsm = terra::rast(dsm_filepath_foc)
  dtm = terra::rast(dtm_filepath_foc)

  # Make sure the rasters are in the same CRS
  dtm = terra::project(dtm, dsm)

  # Compute the CHM
  chm = dsm - dtm

  return(chm)

}


postprocess_photogrammetry = function(mission_id_foc) {

  # Change ownership of all these files to the current user
  command = paste0(
    "sudo chown ", Sys.getenv("USER"), " ",
    file.path(PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_SUBDIR, paste0(mission_id_foc, "_*"))
  )
  system(command)

  # Get the list of photogrammetry outputs from the most recent processing run of this mission. TODO:
  # note that if the naming convention for the processing run ID changes, taking the `max` may not
  # necessarily get the most recent run.

  output_filenames_allruns = list.files(
    path = file.path(PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_SUBDIR),
    pattern = paste0("^", mission_id_foc, "_")
  )

  output_files_allruns = data.frame(
    output_filenames = output_filenames_allruns
  )

  # Extract the processing run ID from all products. The filename convention is
  # <mission_id>_<processing_run_id>_<output_type>.<ext>
  output_files_allruns = output_files_allruns |>
    dplyr::mutate(run_id = str_split(output_filenames, "_") |> map(2) |> unlist()) |>
    # Extract the extension
    dplyr::mutate(extension = str_split(output_filenames, "\\.") |> map(2) |> unlist()) |>
    dplyr::mutate(type = str_split(output_filenames, "_") |> map(3) |> unlist()) |>
    dplyr::mutate(type = str_split(type, "\\.") |> map(1) |> unlist())

  # The highest value is the most recent processing run
  run_id_foc = max(output_files_allruns$run_id)

  # Get the list of files in the most recent run. These are the ones to post-process.
  output_files_foc = output_files_allruns |>
    dplyr::filter(run_id == run_id_foc)

  # Create output folder
  output_path = file.path(
    PHOTOGRAMMETRY_DIR, PHOTOGRAMMETRY_POSTPROCESSED_SUBDIR,
    mission_id_foc, paste0("processed_", run_id_foc)
  )

  # Create the output file dirs
  create_dir(output_path)
  create_dir(file.path(output_path, "full"))
  create_dir(file.path(output_path, "thumbnails"))


  ## Get the mission polygon (for cropping) by downloading from cyverse

  cyverse_url = paste0(
    "https://data.cyverse.org/dav-anon/",
    CYVERSE_MISSIONS_DIR, mission_id_foc,
    "/metadata-mission/",
    mission_id_foc, "_mission-metadata.gpkg"
  )

  tempfile = tempfile(fileext = ".gpkg")

  download.file(
    url = cyverse_url,
    destfile = tempfile,
    method = "curl",
    extra = "-L" # Follow redirects
  )

  # Confirm that file download worked
  if (file.exists(tempfile)) {
    size = file.info(tempfile)$size

    if (size < 5000) {
      warning("Downloaded mission polygon file is too small (is it an error message?) for mission: ", mission_id_foc)
      return(FALSE)
    }
  } else {
      warning("Mission polygon download failed (no server response or no remote file) for mission: ", mission_id_foc)
      return(FALSE)
  }

  mission_polygon = st_read(tempfile)


  ## Crop DSMs, DTM, ortho to mission polygon and write as COG

  rast_file_foc = output_files_foc |>
    dplyr::filter(extension %in% c("tif", "tiff")) |>
    dplyr::pull(output_filenames)

  rast_filepath_foc = file.path(
    PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_SUBDIR,
    rast_file_foc
  )

  # Apply the raster crop & write function to all the rasters
  purrr::walk(
    rast_filepath_foc,
    crop_raster_save_cog,
    output_path = output_path,
    mission_polygon = mission_polygon
  )

  ## Make CHMs

  # Determine what would be the filepaths of the potential DSM and DTM files (if they exist):
  dem_filepath_foc = output_files_foc |>
    dplyr::filter(extension %in% c("tif", "tiff")) |>
    dplyr::filter(type %in% c("dsm-ptcloud", "dsm-mesh", "dtm-ptcloud")) |>
    # Add the paths of the cropped, COG versions
    dplyr::mutate(full_filepath = file.path(
      output_path, "full", output_filenames
    ))

  # If there is both a dsm-mesh and a dtm-ptcloud, make a chm-mesh
  if("dsm-mesh" %in% dem_filepath_foc$type && 
      "dtm-ptcloud" %in% dem_filepath_foc$type) {

    dsm_filepath_foc = dem_filepath_foc |>
      dplyr::filter(type == "dsm-mesh") |>
      dplyr::pull(full_filepath)

    dtm_filepath_foc = dem_filepath_foc |>
      dplyr::filter(type == "dtm-ptcloud") |>
      dplyr::pull(full_filepath)

    chm = make_chm(dsm_filepath_foc, dtm_filepath_foc)

    # Write the CHM as a COG
    output_file_path = file.path(
      output_path, "full",
      paste0(mission_id_foc, "_", run_id_foc, "_chm-mesh.tif")
    )

    terra::writeRaster(
      chm,
      output_file_path,
      overwrite = TRUE,
      filetype = "COG",
      gdal = "BIGTIFF=IF_SAFER",
      todisk = TRUE
    )
  }

  # If there is both a dsm-ptcloud and a dtm-ptcloud, make a chm-ptcloud
  if("dsm-ptcloud" %in% dem_filepath_foc$type && 
      "dtm-ptcloud" %in% dem_filepath_foc$type) {

    dsm_filepath_foc = dem_filepath_foc |>
      dplyr::filter(type == "dsm-ptcloud") |>
      dplyr::pull(full_filepath)

    dtm_filepath_foc = dem_filepath_foc |>
      dplyr::filter(type == "dtm-ptcloud") |>
      dplyr::pull(full_filepath)

    chm = make_chm(dsm_filepath_foc, dtm_filepath_foc)

    # Write the CHM as a COG
    output_file_path = file.path(
      output_path, "full",
      paste0(mission_id_foc, "_", run_id_foc, "_chm-ptcloud.tif")
    )

    terra::writeRaster(
      chm,
      output_file_path,
      overwrite = TRUE,
      filetype = "COG",
      gdal = "BIGTIFF=IF_SAFER",
      todisk = TRUE
    )
  }


  ## Crop and save point cloud as COPC (if it exists)

  if ("points" %in% output_files_foc$type) {

    # Get the file path of the point cloud
    point_cloud_filename_foc = output_files_foc |>
      dplyr::filter(type == "points") |>
      dplyr::pull(output_filenames)

    # Just in case there is anomalously more than one, we will take the first
    point_cloud_filename_foc = point_cloud_filename_foc[1]

    # Get the full path of the point cloud
    input_filepath = file.path(
      PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_SUBDIR,
      point_cloud_filename_foc
    )

    # Get the output file path
    output_filepath = file.path(
      output_path, "full",
      paste0(mission_id_foc, "_", run_id_foc, "_points.laz")
    )

    # Read the pointcloud
    time = system.time(point_cloud <- lidR::readLAS(input_filepath))
    print(paste0("Point cloud loading took ", time[[3]]))
    point_cloud_crs = sf::st_crs(point_cloud)
    mission_polygon_in_point_cloud_crs = sf::st_transform(mission_polygon, point_cloud_crs)

    time = system.time(cropped_point_cloud <- lidR::clip_roi(point_cloud, mission_polygon_in_point_cloud_crs))
    print(paste0("Cropping point cloud took ", time[[3]]))
    cropped_file = tempfile(fileext = ".laz")
    time = system.time(lidR::writeLAS(cropped_point_cloud, cropped_file))
    print(paste0("Saving point cloud took ", time[[3]]))

    untwine_command = paste0("untwine -i ", cropped_file, " -o ", output_filepath)
    conda_command = paste0("/home/exouser/miniconda3/condabin/conda run -n untwine ", untwine_command)
    time = system.time(system(conda_command))
    print(paste0("Converting to COPC took ", time[[3]]))

    file.remove(cropped_file)

    if (!file.exists(output_filepath)) {
      warning(paste0("Failed to create ", output_file_path))
    }
  }

  # Copy any other files that are not raster or point clouds straight to the output folder. TODO:
  # Spatially clip the mesh to the mission polygon before copying.
  other_files = output_files_foc |>
    dplyr::filter(!(extension %in% c("tif", "laz"))) |>
    dplyr::mutate(output_filepath_full = file.path(
      output_path, "full", output_filenames
    )) |>
    dplyr::mutate(input_filepath_full = file.path(
      PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_SUBDIR,
      output_filenames
    ))

  file.link(
    other_files$input_filepath_full,
    other_files$output_filepath_full
  )


  ## Make thumbnails

  # List all tifs within the output folder
  tif_files = list.files(file.path(output_path, "full"), "*.tif", full.names = FALSE)

  # For each full-resolution tif file
  for (tif_file in tif_files) {
    # Full path to the tif file
    tif_file_path = file.path(output_path, "full", tif_file)
    # Create the output file in the thumbnails folder with the same name but png extension
    output_file = str_replace(file.path(output_path, "thumbnails", tif_file), "tif$", "png")

    # if ((skip_existing && file.exists(output_file))) {
    #   print(paste0("Skipping creation of existing thumbnail: ", output_file))
    #   next()
    # }

    print(paste0("Creating thumbnail: ", output_file))

    # Read the raster
    raster = terra::rast(tif_file_path)
    # Compute the number of rows and columns
    n_row = terra::nrow(raster)
    n_col = terra::ncol(raster)
    # Compute the maximum dimension and determine the scale factor to make it match the specified
    # maximum size
    max_dim = max(n_row, n_col)
    scale_factor = OUTPUT_MAX_DIM / max_dim
    new_n_row = floor(n_row * scale_factor)
    new_n_col = floor(n_col * scale_factor)

    # Specify a PNG file as the output device
    # Make sure the background is transperant
    png(output_file, width = new_n_col, height = new_n_row, bg = "transparent")

    # Determine whether this is scalar or RGB data
    n_lyr = terra::nlyr(raster)
    if (n_lyr == 1) {
      # The mar argument ensures that there is not an excessive white border around the image
      plot(raster, axes = FALSE, legend = FALSE, mar = c(0, 0, 0, 0)) # , bg = "transparent")
    } else if (n_lyr %in% c(3, 4)) {
      # Make sure the background is transperant
      terra::plotRGB(raster, bgalpha = 0)
    } else {
      stop(paste0("Input data had an unexpected number of layers: ", str(n_lyr)))
    }

    # Close the PNG
    dev.off()
  }


  # Delete the files from this mission * processing_run from the photogrammetry-outputs folder
  command = paste0(
    "sudo rm -rf ", file.path(PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_SUBDIR, paste0(mission_id_foc, "_", run_id_foc, "_*"))
  )
  system(command)

  # TODO: Also delete Metashape project?

}
