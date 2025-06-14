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

# Function to crop raster to mission polygon and write as COG
crop_raster_save_cog = function(raster_filepath_foc, output_filename, mission_polygon, output_path) {

  # Read, crop, and write the raster as COG
  raster = terra::rast(raster_filepath_foc)
  mission_polygon_matchedcrs = st_transform(mission_polygon, st_crs(raster))
  raster_cropped = terra::crop(raster, mission_polygon_matchedcrs, mask = TRUE)
  output_file_path = file.path(
    output_path, "full",
    output_filename
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


postprocess_photogrammetry = function(mission_id_foc, config_id_foc) {


  # Download the files to process

  # Construct data download command line call
  local_dir = file.path(PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_DOWNLOADED_SUBDIR)
  remote_dir = paste0(RCLONE_REMOTE, ":", REMOTE_PHOTOGRAMMETRY_DIR, paste0("config_", config_id_foc))
  filter_phrase = paste0("--include ", config_id_foc, "_", mission_id_foc, "_*")
  command = paste("rclone copy", remote_dir, local_dir, filter_phrase, "--max-depth 1 --progress --transfers 32 --checkers 32 --stats 1s --retries 5 --retries-sleep=15s --s3-upload-cutoff 100Mi --s3-chunk-size 100Mi --s3-upload-concurrency 16", sep = " ")
  result = system(command)

  # Check if we got an error, and if so, retry up to 10 times, waiting 1 minute between each
  tries = 0
  while (result != 0 && tries < 10) {
    tries = tries + 1
    Sys.sleep(60) # Wait 1 minute before retrying
    cat("\n **** Uploading photogrammetry outputs to object store for mission", mission_id_foc, "failed. Retrying... (attempt", tries, ") **** \n")
    result = system(command)
  }

  # If we still got an upload error, print a warning and save to log and return false (don't do the
  # next step which is deleting the local directory)
  if (result != 0) {
    toprint = (paste(Sys.time(), "- Error uploading photogrammetry outputs to object store for mission", mission_id_foc, "(all 10 tries failed).\n"))
    warning(toprint)
    write(toprint, file = UPLOAD_ERROR_LOG, append = TRUE)
    return(FALSE)
  }

  # Get the list of photogrammetry outputs from the specified processing run of this mission
  # (config_id_foc) (the files we just downloaded)

  photogrammetry_output_filenames = list.files(
    path = local_dir,
    pattern = paste0("^", config_id_foc, "_", mission_id_foc, "_"),
  )

  photogrammetry_output_files = data.frame(
    photogrammetry_output_filename = photogrammetry_output_filenames
  )

  # Make a data frame with the product type, mission ID, processing config ID, and extension for
  # each file, and also the postprocessed filename.
  photogrammetry_output_files = photogrammetry_output_files |>
    dplyr::mutate(extension = str_split(photogrammetry_output_filename, "\\.") |> map(2) |> unlist()) |>
    dplyr::mutate(type = str_split(photogrammetry_output_filename, "_") |> map(3) |> unlist()) |>
    dplyr::mutate(type = str_split(type, "\\.") |> map(1) |> unlist()) |>
    dplyr::mutate(postprocessed_filename = paste0(
      mission_id_foc, "_", type, ".", extension
    ))

  # Create output folder
  postprocessed_path = file.path(
    PHOTOGRAMMETRY_DIR, PHOTOGRAMMETRY_POSTPROCESSED_SUBDIR,
    mission_id_foc, paste0("processed_", config_id_foc)
  )

  # Create the output file dirs
  create_dir(postprocessed_path)
  create_dir(file.path(postprocessed_path, "full"))
  create_dir(file.path(postprocessed_path, "thumbnails"))

  ## Get the mission polygon (for cropping) by downloading from Object Store
  polygon_temp_file = tempfile(fileext = ".gpkg")
  remote_file = paste0(RCLONE_REMOTE, ":", REMOTE_MISSIONS_DIR, mission_id_foc, "/metadata-mission/", mission_id_foc, "_mission-metadata.gpkg")
  command = paste("rclone copyto", remote_file, polygon_temp_file, "--progress --transfers 32 --checkers 32 --stats 1s --retries 5 --retries-sleep=15s --s3-upload-cutoff 100Mi --s3-chunk-size 100Mi --s3-upload-concurrency 16 --multi-thread-streams 2", sep = " ")
  system(command)

  # Confirm that file download worked
  if (file.exists(polygon_temp_file)) {
    size = file.info(polygon_temp_file)$size

    if (size < 5000) {
      warning("Downloaded mission polygon file is too small (is it an error message?) for mission: ", mission_id_foc)
      return(FALSE)
    }
  } else {
      warning("Mission polygon download failed (no server response or no remote file) for mission: ", mission_id_foc)
      return(FALSE)
  }

  mission_polygon = st_read(polygon_temp_file)


  ## Crop DSMs, DTM, ortho to mission polygon and write as COG

  photogrammetry_outputs_rast = photogrammetry_output_files |>
    dplyr::filter(extension %in% c("tif", "tiff"))

  photogrammetry_output_rast_filepaths = file.path(
    PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_DOWNLOADED_SUBDIR,
    photogrammetry_outputs_rast$photogrammetry_output_filename
  )

  photogrammetry_postprocessed_filenames = photogrammetry_outputs_rast$postprocessed_filename

  # Apply the raster crop & write function to all the rasters
  purrr::walk2(
    photogrammetry_output_rast_filepaths,
    photogrammetry_postprocessed_filenames,
    crop_raster_save_cog,
    output_path = postprocessed_path,
    mission_polygon = mission_polygon
  )

  ## Make CHMs

  # Determine what would be the filepaths of the potential DSM and DTM files (if they exist):
  dem_filepaths = photogrammetry_output_files |>
    dplyr::filter(extension %in% c("tif", "tiff")) |>
    dplyr::filter(type %in% c("dsm-ptcloud", "dsm-mesh", "dtm-ptcloud")) |>
    # Add the paths of the cropped, COG versions
    dplyr::mutate(postprocessed_filepath = file.path(
      postprocessed_path, "full", postprocessed_filename
    ))

  # If there is both a dsm-mesh and a dtm-ptcloud, make a chm-mesh
  if("dsm-mesh" %in% dem_filepaths$type &&
       "dtm-ptcloud" %in% dem_filepaths$type) {

    dsm_filepaths = dem_filepaths |>
      dplyr::filter(type == "dsm-mesh") |>
      dplyr::pull(postprocessed_filepath)

    dtm_filepaths = dem_filepaths |>
      dplyr::filter(type == "dtm-ptcloud") |>
      dplyr::pull(postprocessed_filepath)

    chm = make_chm(dsm_filepaths, dtm_filepaths)

    # Write the CHM as a COG
    postprocessed_file_path = file.path(
      postprocessed_path, "full",
      paste0(mission_id_foc, "_chm-mesh.tif")
    )

    terra::writeRaster(
      chm,
      postprocessed_file_path,
      overwrite = TRUE,
      filetype = "COG",
      gdal = "BIGTIFF=IF_SAFER",
      todisk = TRUE
    )
  }

  # If there is both a dsm-ptcloud and a dtm-ptcloud, make a chm-ptcloud
  if ("dsm-ptcloud" %in% dem_filepaths$type && 
      "dtm-ptcloud" %in% dem_filepaths$type) {

    dsm_filepaths = dem_filepaths |>
      dplyr::filter(type == "dsm-ptcloud") |>
      dplyr::pull(postprocessed_filepath)

    dtm_filepaths = dem_filepaths |>
      dplyr::filter(type == "dtm-ptcloud") |>
      dplyr::pull(postprocessed_filepath)

    chm = make_chm(dsm_filepaths, dtm_filepaths)

    # Write the CHM as a COG
    output_file_path = file.path(
      postprocessed_path, "full",
      paste0(mission_id_foc, "_chm-ptcloud.tif")
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


  # Crop and save point cloud as COPC (if it exists) (if specified) (NOTE that currently, we are
  # having Metashape save a COPC so this is skipped) TODO: Consider naming the output to
  # *_points-copc.laz instead of *_points.laz, like we are having Metashape do. But note that if we
  # do this, then the code below that copies all other files will have to be modified to not copy
  # the *_points-copc.laz file  (why? I think because then it would overwrite the just-converted
  # one, but I'm not sure whether the photogrammetry outputs folder would ever have a copc and
  # non-copc pointcloud).

  if (CONVERT_TO_COPC & ("points" %in% photogrammetry_output_files$type)) {

    # Get the file path of the point cloud
    point_cloud_filename = photogrammetry_output_files |>
      dplyr::filter(type == "points") |>
      dplyr::pull(photogrammetry_output_filename)

    # Just in case there is anomalously more than one, we will take the first
    point_cloud_filename = point_cloud_filename[1]

    # Get the full path of the point cloud
    input_filepath = file.path(
      PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_DOWNLOADED_SUBDIR,
      point_cloud_filename
    )

    # Get the output file path
    output_filepath = file.path(
      postprocessed_path, "full",
      paste0(mission_id_foc, "_points.laz")
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

  # Copy any other files that are not raster or point clouds straight to the output folder. Except
  # that COPC-format points (filename ending *_points-copc.laz) should get copied because they are
  # skipped by the COPC conversion above. TODO: Spatially clip the COPC and mesh to the mission polygon
  # before copying.
  other_files = photogrammetry_output_files |>
    dplyr::filter((!(extension %in% c("tif", "laz"))) | (type == "points-copc")) |>
    dplyr::mutate(postprocessed_filepath_full = file.path(
      postprocessed_path, "full", postprocessed_filename
    )) |>
    dplyr::mutate(input_filepath_full = file.path(
      PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_DOWNLOADED_SUBDIR,
      photogrammetry_output_filename
    ))

  file.link(
    other_files$input_filepath_full,
    other_files$postprocessed_filepath_full
  )


  ## Make thumbnails

  # List all tifs within the output folder
  tif_files = list.files(file.path(postprocessed_path, "full"), "*.tif", full.names = FALSE)

  # For each full-resolution tif file
  for (tif_file in tif_files) {
    # Full path to the tif file
    tif_file_path = file.path(postprocessed_path, "full", tif_file)
    # Create the output file in the thumbnails folder with the same name but png extension
    thumbnail_filepath = str_replace(file.path(postprocessed_path, "thumbnails", tif_file), "tif$", "png")

    # if ((skip_existing && file.exists(output_file))) {
    #   print(paste0("Skipping creation of existing thumbnail: ", output_file))
    #   next()
    # }

    print(paste0("Creating thumbnail: ", thumbnail_filepath))

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
    png(thumbnail_filepath, width = new_n_col, height = new_n_row, bg = "transparent")

    # Determine whether this is single-channel or RGB data
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

  ## Delete the intermediate files no longer needed
  
  # TODO: Consider breaking out the deletion steps below into a separate function(s); some (the ones
  # that delete the photogrammetry inputs) can come even sooner than here.

  # Delete the files from this mission * processing_run from the photogrammetry-outputs-downloaded folder. Need
  # to run this as a system command with sudo in case the files were owned by root due to being
  # created by Docker.
  command = paste0(
    "sudo rm -rf ",
    file.path(
      PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_DOWNLOADED_SUBDIR,
      paste0(config_id_foc, "_", mission_id_foc, "_*")
    )
  )
  system(command)

  # Delete the mission polygon tempfile
  file.remove(polygon_temp_file)

  return(TRUE)

}
